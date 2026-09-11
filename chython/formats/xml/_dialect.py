# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""The engine that reads and writes an XML dialect declared as a table.

CML and MRV differ only in spelling, so there is one walker and the vocabulary is data: a
:class:`Dialect` is element names plus a :class:`Field` per attribute.  An attribute with no row, and a
child element no handler claims, are logged ``unsupported: `` here rather than per dialect.  The neutral
intermediate is :class:`~chython.formats.ctfile._ctab.Ctab`, so stereo, hydrogen counts and S-groups are
the same code as in the MDL readers.
"""

from collections.abc import Callable
from re import DOTALL, compile as re_compile
from typing import NamedTuple
from xml.etree.ElementTree import Element, SubElement, indent as ETIndent, tostring as ETToString

from ._errors import MalformedXml, UnsupportedXml, XmlError
from ._tree import local, namespace_of, parse_xml, text_of
from ..ctfile import Ctab, CtabAtom, CtabBond
from ..ctfile._ctab import LABEL_ELEMENT
from ..ctfile._hydrogens import StatedChannels
from ..ctfile._sgroup import UNSUPPORTED
from ...core import LogRecord, LOST, REPAIRED
from ...core._core import element_symbols, R_INDEX_MAX


__all__ = ['Field', 'Tags', 'Dialect', 'Record', 'NotModelled', 'register', 'dialect', 'dialects',
           'sniff', 'apply_fields', 'read_document', 'read_molecule', 'write_molecule',
           'molecule_nodes', 'parse_records', 'read_molecules', 'build_records',
           'parse_xml_document', 'read_xml',
           'decode_int', 'decode_float', 'encode_nonzero', 'encode_stated', 'encode_from_table',
           'encode_coordinate', 'encode_element', 'xml_text', 'synthetic_node', 'outermost',
           'resolve_element',
           'choose_coordinates', 'emit_coordinates', 'read_properties', 'emit_properties', 'cdata',
           'serialize']


class NotModelled(ValueError):
    """Raised by a :attr:`Field.decode` for a value the format has and this library does not.

    Logged with the ``unsupported: `` prefix.  The other three answers a decode may give: a value;
    ``ValueError`` for a *malformed* value, logged as a repair and costing the attribute; an
    :class:`~._errors.XmlError`, which refuses the whole record and is what a pseudo-atom raises.
    CML's ``order="partial12"`` is the motivating case -- a legal delocalised order with no
    representation in a molecule.
    """


class Field(NamedTuple):
    """One compiled table row: a file attribute, the slot it fills, and a codec both ways.

    `slot` names a :class:`~chython.formats.ctfile._ctab.CtabAtom` or ``CtabBond`` attribute, or -- when
    the quantity has no home on the intermediate, like CML's *total* ``hydrogenCount`` -- a key in the
    per-item spill, which both directions use under the same name.  `decode` raises ``ValueError`` on
    garbage; `encode` returns ``None`` to omit the attribute, which is how a default stays unwritten.
    Either may be ``None``, making the row write-only or read-only.  `write_slot` is where the *writer*
    reads a value the reader parked raw for :attr:`Dialect.finish` to interpret -- CML's spin
    multiplicity -- since one attribute may appear only once in the index.
    """
    attribute: str
    slot: str
    decode: Callable | None = None
    encode: Callable | None = None
    write_slot: str | None = None


class Tags(NamedTuple):
    """The element names of one dialect.  Local names only -- see :func:`~._tree.local` for why."""
    molecule: str = 'molecule'
    atom_array: str = 'atomArray'
    atom: str = 'atom'
    bond_array: str = 'bondArray'
    bond: str = 'bond'
    #: Local names this dialect will accept as a document root.  Used by :func:`sniff` only.
    roots: frozenset = frozenset(('cml', 'molecule'))


class Dialect(NamedTuple):
    """One compiled dialect: a vocabulary, plus the handlers a vocabulary cannot express.

    Every handler defaults to ``None``, so a dialect declares only what it has:

    * ``atom_child(node, record, position, log) -> bool`` -- a child of an ``<atom>``; ``True`` means
      claimed.  CML's ``<atomParity>``, CML 1's ``<string builtin="elementType">``.
    * ``bond_child(node, record, position, log) -> bool`` -- a child of a ``<bond>``; ``<bondStereo>``.
    * ``molecule_child(node, record, log) -> bool`` -- a ``<molecule>`` child that is neither array.
    * ``atom_array_hook(node, record, log) -> bool`` / ``bond_array_hook`` -- the array element itself
      carrying CML's whitespace-separated column form; ``True`` means the engine must not also walk the
      children.
    * ``finish(record, log)`` -- once, after the walk and before the build; where a field depending on
      the finished graph resolves (a total hydrogen count, a parity over atom ids).
    * ``emit_molecule(record, element, log)`` / ``emit_atom(record, position, element, log)`` /
      ``emit_bond(record, position, element, log)`` -- the write-side mirror of the child handlers.
    * ``document(root, log)`` -- once per document, before the molecules; the acceptance rule *above* a
      molecule, since :func:`molecule_nodes` returns a flat list and only the dialect knows what that
      loses (CML's ``<reaction>`` states its molecules' roles).

    `namespaces` is every URI :func:`sniff` matches -- real files declare several historical ones for
    the same vocabulary -- and `ns` is the single one this dialect writes.
    """
    name: str
    #: How this dialect spells a stated hydrogen count and a stated total valence, for the advice in
    #: :func:`~..ctfile._hydrogens.calc_implicit`'s unknown-count lines.  Has no default on purpose: the
    #: derivation is shared with the CTAB versions, so a silent dialect would inherit MDL's spellings and
    #: advise an XML reader to add a data S-group its document cannot hold.  `valence` is ``None`` for a
    #: dialect with no total-valence channel, as CML is.
    channels: StatedChannels
    namespaces: frozenset = frozenset()
    ns: str = ''
    tags: Tags = Tags()
    #: The attribute holding an atom's identity, which the engine resolves bond references against.
    atom_id: str = 'id'
    #: Attribute names that may hold a bond's two endpoint ids, most preferred first.
    bond_refs: tuple = ('atomRefs2',)
    #: Two attributes holding one endpoint id each, tried when none of `bond_refs` is present.  CML
    #: spells its endpoints this way in the array form and, legally, on a `<bond>` element too.
    bond_ref_pair: tuple = ()
    molecule_fields: tuple = ()
    atom_fields: tuple = ()
    bond_fields: tuple = ()
    #: Attributes consumed elsewhere or carrying no chemistry, so *not* worth an `unsupported:` line.
    #: An attribute belongs here when a reader that honoured it would produce the same molecule.
    molecule_ignored: frozenset = frozenset()
    atom_ignored: frozenset = frozenset()
    bond_ignored: frozenset = frozenset()
    #: Child element names of `<molecule>` that carry no chemistry, likewise.
    molecule_children_ignored: frozenset = frozenset()
    atom_child: Callable | None = None
    bond_child: Callable | None = None
    molecule_child: Callable | None = None
    atom_array_hook: Callable | None = None
    bond_array_hook: Callable | None = None
    finish: Callable | None = None
    document: Callable | None = None
    emit_molecule: Callable | None = None
    emit_atom: Callable | None = None
    emit_bond: Callable | None = None

    def index(self, which):
        """``{attribute: Field}`` for ``'molecule'``, ``'atom'`` or ``'bond'``.

        Built per call -- a ``NamedTuple`` has nowhere to cache -- so every caller must hoist the result
        out of its loop.  The write path takes the row sequence directly instead.
        """
        rows = {'molecule': self.molecule_fields, 'atom': self.atom_fields,
                'bond': self.bond_fields}[which]
        return {row.attribute: row for row in rows}


class Record:
    """One molecule as this layer holds it: a ``Ctab``, the file's atom ids, and a scratch dict.

    One type for both directions -- the reader fills it from a tree for :meth:`Ctab.build`, the writer
    from a molecule for :func:`write_molecule` -- so the two cannot diverge on what an order or a wedge
    means.  `ids` is positional, ``ids[i]`` naming ``ctab.atoms[i]``, and `index_of` is its inverse.
    `extras`, `atom_extras` and `bond_extras` are the dialect's own space: ``CtabAtom`` is
    ``__slots__``-ed, so a value that cannot go onto it lands in the parallel dict.
    """
    __slots__ = ('ctab', 'ids', 'index_of', 'extras', 'atom_extras', 'bond_extras')

    def __init__(self, ctab=None):
        self.ctab = Ctab() if ctab is None else ctab
        self.ids = []
        self.index_of = {}
        self.extras = {}
        self.atom_extras = []
        self.bond_extras = []

    def add_atom(self, atom, ident=None, log=None):
        """Append `atom`, register `ident`, return its 0-based position.

        `log` takes the line a duplicate `ident` earns, and is optional because a caller synthesising
        ids cannot collide.
        """
        position = len(self.ctab.atoms)
        self.ctab.atoms.append(atom)
        if ident is None:
            ident = f'a{position + 1}'
        self.ids.append(ident)
        self.atom_extras.append({})
        # First declaration wins: letting the later atom claim the name would silently move every bond
        # written before it.  Logged unprefixed -- an XML `id` is unique by definition, so the file is
        # broken -- and the atom is kept, because a bond may have resolved to the wrong one of the two
        # and nothing downstream would look odd.
        if ident in self.index_of:
            if log is not None:
                log.append(LogRecord('xml:duplicate-atom-id', (),
                                     f'atom {ident}: id declared twice; the first declaration keeps the name, so '
                                     f'a reference to it names that atom and not this one', LOST))
        else:
            self.index_of[ident] = position
        return position

    def add_bond(self, bond):
        """Append `bond`, return its 0-based position."""
        position = len(self.ctab.bonds)
        self.ctab.bonds.append(bond)
        self.bond_extras.append({})
        return position

    def __len__(self):
        return len(self.ctab.atoms)

    def __repr__(self):
        return (f'Record({len(self.ctab.atoms)} atoms, {len(self.ctab.bonds)} bonds, '
                f'extras={sorted(self.extras)})')


# the dialect registry

#: name -> Dialect.  Populated by :func:`_load` on first use and never at import.
_DIALECT_CACHE = {}

#: Whether :func:`_load` has run.  A flag and not ``if _DIALECT_CACHE``: `register` calls `_load` so an
#: outside dialect cannot register ahead of the shipped ones, and emptiness as the guard would make
#: those two recurse forever.
_LOADED = False


def _load():
    """Import and register every dialect shipped here.  Idempotent, and lazy by design."""
    global _LOADED
    if _LOADED:
        return
    _LOADED = True  # set first: `register` calls back into here and must find the door shut
    from ._cml import CML
    register(CML)
    from ._mrv import MRV
    register(MRV)


def register(dial):
    """Register `dial` under its own name, replacing any dialect already there.

    Public so a dialect outside this package -- a house format, a vendor variant of CML -- is a table
    someone registers rather than a fork of this module.
    """
    _load()
    _DIALECT_CACHE[dial.name] = dial
    return dial


def dialect(name):
    """The registered dialect called `name`."""
    _load()
    try:
        return _DIALECT_CACHE[name]
    except KeyError:
        raise ValueError(f'no XML dialect named {name!r}; have {sorted(_DIALECT_CACHE)}') from None


def dialects():
    """Every registered dialect name, sorted."""
    _load()
    return tuple(sorted(_DIALECT_CACHE))


def sniff(root, log=None):
    """The dialect that fits document `root`, by namespace first and by root element name second.

    A document whose namespace no dialect claims is read by the dialect claiming its root name, with a
    line naming the unrecognised namespace -- vendor variants spell atoms and bonds nearly as CML does,
    so refusing would lose a readable file.  A document declaring no namespace names no vocabulary, so
    it gets no line and is still chosen by root name.
    """
    _load()
    seen = {namespace_of(node.tag) for node in root.iter()}
    seen.discard('')
    for name in sorted(_DIALECT_CACHE):
        dial = _DIALECT_CACHE[name]
        if seen & dial.namespaces:
            return dial
    tag = local(root.tag)
    for name in sorted(_DIALECT_CACHE):
        dial = _DIALECT_CACHE[name]
        if tag in dial.tags.roots:
            if log is not None and seen:
                # Prefixed, so a caller feeding in a corpus of vendor files can learn mechanically that
                # they were read as another dialect.  Only when a namespace was declared and went
                # unclaimed: a document declaring none is evidence of nothing, and prefixing a
                # losslessly read file answers "did we lose something?" with an untrue yes.
                log.append(LogRecord('xml:unclaimed-namespace', (),
                                     f'{UNSUPPORTED}record: no dialect claims namespace(s) '
                                     f'{", ".join(sorted(seen))}; <{tag}> read as {dial.name}', LOST))
            return dial
    known = ', '.join(sorted(_DIALECT_CACHE))
    raise MalformedXml(f'<{tag}> is not a document root any registered dialect ({known}) reads')


# reading

def apply_fields(node, table, target, log, where, ignored=(), structural=(), spill=None):
    """Apply every attribute of `node` that `table` has a row for; log the rest ``unsupported: ``.

    The acceptance rule for every dialect, enforced once here.  `ignored` is the third case and is
    narrow: an attribute belongs there only when a reader that honoured it would build the same
    molecule.  A row's ``decode`` may answer four ways:

    * a value -- applied to `target`, or to `spill` when `target` has no such slot;
    * :class:`NotModelled` -- an ``unsupported: `` line;
    * ``ValueError`` or ``TypeError`` -- a plain line; a malformed field costs the attribute, not the
      atom;
    * :class:`~._errors.XmlError` -- re-raised.  An ``elementType`` naming a pseudo-atom cannot degrade;
      defaulting it to carbon would be a silent wrong answer.
    """
    for name, value in node.attrib.items():
        key = local(name)
        row = table.get(key)
        if row is None or row.decode is None:
            if key in ignored or key in structural:
                continue
            log.append(LogRecord('xml:attribute-not-modelled', (),
                                 f'{UNSUPPORTED}{where}: attribute {key}={value!r:.40} is not modelled', LOST))
            continue
        try:
            decoded = row.decode(value)
        except NotModelled as e:
            log.append(LogRecord('xml:attribute-not-modelled', (),
                                 f'{UNSUPPORTED}{where}: attribute {key}={value!r:.40}: {e}', LOST))
            continue
        except XmlError:
            raise
        except (ValueError, TypeError) as e:
            log.append(LogRecord('xml:attribute-not-read', (),
                                 f'{where}: attribute {key}={value!r:.40} not read: {e}', LOST))
            continue
        try:
            setattr(target, row.slot, decoded)
        except AttributeError:
            # The row's value has no home on the intermediate -- CML's third coordinate, its spin
            # multiplicity -- so it lands in the spill and `finish` picks it up.
            if spill is None:
                raise
            spill[row.slot] = decoded


# the shared half of a dialect -- codecs and handlers that belong to no dialect
#
# The test for putting something here is that its body reads no per-dialect name; `resolve_element`
# passes it by taking the pseudo-atom table as an argument, since that table *is* vocabulary.  A
# dialect's bond-order table, radical names and `<bondStereo>` spellings fail it and stay where they are.

#: Element symbols by atomic number, for :func:`encode_element`.  One copy for the package.
_SYMBOLS = element_symbols()


def decode_int(text):
    """An integer attribute.  ``1.0`` is accepted: writers with one float formatter for everything."""
    text = text.strip()
    try:
        return int(text)
    except ValueError:
        return int(float(text))  # raises ValueError itself on real garbage, which the engine logs


def decode_float(text):
    """A float attribute."""
    return float(text.strip())


def encode_nonzero(value):
    """`value` as text, or ``None`` for a zero.

    For a field whose zero is the format's way of saying nothing: a charge, an isotope, a map number.
    ``mrvMap="0"`` is MRV for unmapped, so neither it nor ``formalCharge="0"`` is written.  Contrast
    :func:`encode_stated`, where a zero is a statement.
    """
    return str(value) if value else None


def encode_stated(value):
    """`value` as text, or ``None`` for ``None``.

    For a field whose zero is a *stated* zero -- a hydrogen count, a total valence -- so ``0`` is written
    and only "nothing determined this" is left off.  ``H_UNKNOWN`` reaches here as ``None``, and an
    absent attribute is every dialect's way of saying "work it out".
    """
    return None if value is None else str(value)


def encode_from_table(table):
    """A codec looking `value` up in `table`, answering ``None`` for a value the table has no key for.

    ``None`` leaves the attribute off, which a dialect that cannot spell an order must do: the dative
    order 8 has no CML or MRV spelling, and ``order="1"`` would turn a coordination contact into a
    covalent bond downstream.  Reporting it is the writer's job -- a codec has no log.
    """
    return table.get


def encode_coordinate(value):
    """One coordinate as text, to 4 decimal places -- MDL's own precision, enough for a drawing.

    Whether there is a layout at all is asked once per record by :func:`emit_coordinates`, the only
    caller; a per-atom test would write ``x2`` for one atom of a flat molecule and not its neighbour.
    """
    return f'{value or 0.0:.4f}'


def encode_element(element):
    """An atomic number as an ``elementType``, the name every dialect spells it with.

    ``element`` is an ``int`` on the write path and a ``str`` on the read path -- a reader parks the
    file's raw text in ``atom_extras`` and its ``finish`` resolves it -- so both are accepted.
    """
    if isinstance(element, str):
        return element
    if element == 0:
        # `R` and not `R#`: `R` is what `molconvert mrv` writes for a marker, and CML's own vocabulary
        # lists it too.  The INDEX is not in this value -- MRV spells it `rgroupRef`, which is a row of
        # its own, and CML has no spelling for one.
        return 'R'
    return _SYMBOLS[element]


def xml_text(text, what, log=None):
    """`text` with every byte no codec accepts replaced by U+FFFD, reporting how many.

    XML 1.0 ADMITS NO LONE SURROGATE -- not escaped, not as a numeric reference, not at all -- so a name
    line that `MoleculeContainer.title` hands over with a byte in it cannot be written.  The loss is the
    format's, so it is taken here and not in the container; every other writer in this tree puts the byte
    back.  A title whose text genuinely contains U+FFFD is untouched and logs nothing, the strict decode
    below being what separates the two.
    """
    raw = text.encode('utf8', 'surrogateescape')
    try:
        raw.decode('utf8')
    except UnicodeDecodeError:
        out = raw.decode('utf8', 'replace')
        if log is not None:
            n = out.count('�')
            log.append(LogRecord('xml:text-not-utf8', (),
                                 f'{UNSUPPORTED}{what}: {n} byte(s) are not valid UTF-8; the record is written '
                                 f'with {n} replacement character(s)', LOST))
        return out
    return text


def synthetic_node(tag, attributes):
    """An ``Element`` standing for one row of an array, so the real table can decode it.

    A synthetic node rather than a second decode path: an array form is rewritten into the element
    form's attribute names and handed to :func:`apply_fields`, so a charge is parsed by one function
    however the file spelled it, and an attribute with no row is reported rather than dropped here.
    """
    node = Element(tag)
    for key, value in attributes.items():
        node.set(key, value)
    return node


def outermost(root, want):
    """Every element named `want` under `root`, outermost only, in document order."""
    found = []
    stack = [root]
    while stack:
        node = stack.pop(0)
        if local(node.tag) == want:
            found.append(node)
            continue
        stack = list(node) + stack
    return found


def resolve_element(text, where, log, pseudo):
    """An ``elementType`` as ``(symbol, isotope, label, r_index)``, or a refusal.

    The same four-value answer as the CTfile resolvers, and the same recoveries -- ``D``/``T``, an
    upper-cased symbol -- because a document converted out of a molfile inherits the molfile's damage.
    `pseudo` is the dialect's own ``{symbol: what it is}`` for symbols a molecule cannot hold, and is an
    argument because it is vocabulary: CML inherits MDL's ``A``/``Q``/``X``, MRV has none.

    `label` is the text where it names no element, with `symbol` then :data:`LABEL_ELEMENT`; `r_index`
    is the group number of an ``R<n>``, which an ``rgroupRef`` attribute may also carry.
    """
    symbol = text.strip()
    if not symbol:
        raise MalformedXml(f'{where}: elementType is empty')
    if symbol == 'D':
        log.append(LogRecord('xml:element-folded', (), f'{where}: D read as hydrogen isotope 2', REPAIRED))
        return 'H', 2, None, 0
    if symbol == 'T':
        log.append(LogRecord('xml:element-folded', (), f'{where}: T read as hydrogen isotope 3', REPAIRED))
        return 'H', 3, None, 0
    # The R family before the pseudo table, and before the element symbols: `R` is Marvin's own spelling
    # for a marker -- `elementType="R" rgroupRef="1"` is what `molconvert mrv` writes for an `M  RGP`
    # atom -- and `*` is the same attachment point unindexed.  Tested EXACTLY: `symbol[0] == 'R'` alone
    # would capture Rb, Re, Rh, Rn, Ra, Rf, Rg and Ru.
    upper = symbol.upper()
    if upper in ('R', 'R#', '*') or (upper[0] == 'R' and upper[1:].isdigit()):
        # Case-folded like a symbol, and for the same reason: a document converted out of a molfile
        # inherits the molfile's upper-casing, and `r1` is nobody's element.
        if symbol != upper:
            log.append(LogRecord('xml:element-folded', (),
                                 f'{where}: elementType {symbol!r} read as {upper!r}', REPAIRED))
        index = int(upper[1:]) if upper[1:].isdigit() else 0
        if index > R_INDEX_MAX:
            log.append(LogRecord('xml:r-index-too-wide', (),
                                 f'{where}: R index {index} is past R_INDEX_MAX ({R_INDEX_MAX}), so the '
                                 f'marker is left unindexed', LOST))
            index = 0
        return 'R', 0, None, index
    # The pseudo table is consulted before the element symbols because a dialect's table claims symbols
    # a molecule cannot hold, and `A`, `Q`, `X` and `M` genuinely are query primitives.
    if symbol in pseudo:
        raise UnsupportedXml(f'{where}: elementType {symbol!r} is {pseudo[symbol]}, which a molecule '
                             f'cannot represent. Read this file with a query reader')
    if symbol in _SYMBOLS:
        return symbol, 0, None, 0
    folded = symbol.capitalize()
    if folded in pseudo:
        raise UnsupportedXml(f'{where}: elementType {symbol!r} is {pseudo[folded]}, which a molecule '
                             f'cannot represent. Read this file with a query reader')
    if folded in _SYMBOLS:
        log.append(LogRecord('xml:element-folded', (), f'{where}: elementType {symbol!r} read as {folded!r}', REPAIRED))
        return folded, 0, None, 0
    # Free text where an element belongs -- `Pol`, `OMe`, a registry identifier.  The same answer the
    # CTfile readers give it: the marker carrying the text, since the record is worth more than the one
    # field nobody can resolve.
    log.append(LogRecord('xml:element-type-as-label', (),
                         f'{where}: elementType {symbol!r} names no element, so it is read as the display '
                         f'label it is: the atom is kept as the marker {LABEL_ELEMENT!r} carrying '
                         f'{symbol!r} as its alias. Whatever the label abbreviates is not in the '
                         f'structure', LOST))
    return LABEL_ELEMENT, 0, symbol, 0


def choose_coordinates(record, log):
    """Choose between ``x2/y2`` and ``x3/y3/z3``, once for the record.

    2D wins where both are present: the wedges and the double-bond geometry are all measured against the
    drawing, and mixing the sets would produce a geometry in no file and then read stereo out of it.  The
    discarded conformer earns an ``unsupported: `` line.  The two tests are asymmetric because
    ``x2``/``y2`` have rows landing on ``CtabAtom.x``/``y`` while the third set spills, so 2D is a
    question about values and 3D one about keys; a dialect growing an ``x3`` row must revisit this.
    """
    atoms = record.ctab.atoms
    extras = record.atom_extras
    flat = any(atom.x or atom.y for atom in atoms)
    solid = any('x3' in s or 'y3' in s or 'z3' in s for s in extras)
    if not solid:
        record.ctab.dimensionality = '2D' if flat else ''
        return
    if flat:
        log.append(LogRecord('xml:coordinates-conflict', (),
                             f'{UNSUPPORTED}coordinates: both 2D and 3D coordinates are present; the 2D drawing '
                             f'is used, because every stereo statement in the record is measured against it, and '
                             f'the 3D conformer is dropped', LOST))
        record.ctab.dimensionality = '2D'
        return
    for atom, spill in zip(atoms, extras):
        atom.x = spill.get('x3', 0.0)
        atom.y = spill.get('y3', 0.0)
        atom.z = spill.get('z3', 0.0)
    record.ctab.dimensionality = '3D'


def emit_coordinates(record, atom, node):
    """One atom's coordinates, in the set ``Ctab.dimensionality`` says the record has.

    ``''`` gets no attributes: ``x2="0.0000" y2="0.0000"`` on every atom is a drawing nobody made, and a
    reader would read stereo out of it.  ``'3D'`` gets ``x3``/``y3``/``z3``, a conformer written as
    ``x2``/``y2`` being a projection relabelled as a drawing with its third coordinate gone.
    """
    dimensionality = record.ctab.dimensionality
    if dimensionality == '3D':
        node.set('x3', encode_coordinate(atom.x))
        node.set('y3', encode_coordinate(atom.y))
        node.set('z3', encode_coordinate(atom.z))
    elif dimensionality:
        node.set('x2', encode_coordinate(atom.x))
        node.set('y2', encode_coordinate(atom.y))


def read_properties(node, record, log, *, rule):
    """``<propertyList>`` or a lone ``<property>`` into ``record.ctab.meta``.  Returns ``True``.

    Shared by both XML dialects: MRV and CML spell a data field the same way, and Marvin's own
    ``sdf``-to-``mrv`` conversion writes exactly this element inside ``<molecule>``.

    Named by ``title`` first and ``dictRef`` second: ``dictRef`` is a reference into a dictionary and
    carries that dictionary's prefix, which is not part of the field name.  Measured on one SD field
    through ``molconvert`` 25.1.3 -- ``mrv`` writes ``dictRef="ID" title="ID"``, ``cml`` writes
    ``dictRef="marvin:ID" title="ID"`` -- so ``title`` is the spelling that agrees across both, and
    ``dictRef`` names the field only in a document that states no title.  The value is the
    ``<scalar>``'s text kept as a
    ``str``: ``dataType`` says ``xsd:double`` on a field whose value is ``'>100'`` often enough that
    coercing would lose data.  ``Ctab.build`` copies the mapping onto ``mol.meta``, so a reader sees the
    properties and not only a parser.
    """
    # Built locally and merged at the end, so a `<propertyList>` from which nothing survived leaves no
    # empty mapping behind for a caller to test the truth of rather than the presence of.
    store = {}
    entries = [node] if local(node.tag) == 'property' else [c for c in node]
    for entry in entries:
        tag = local(entry.tag)
        if tag != 'property':
            log.append(LogRecord(f'{rule}:property-list-child-unknown', (),
                                 f'{UNSUPPORTED}record: <{tag}> in <propertyList> is not modelled', LOST))
            continue
        name = entry.get('title') or entry.get('dictRef')
        scalars = [c for c in entry if local(c.tag) == 'scalar']
        others = [local(c.tag) for c in entry if local(c.tag) != 'scalar']
        if others:
            # `<array>` and `<matrix>` hold a vector or a table per property: a real construct nothing
            # downstream can consume, so it is named rather than flattened.
            log.append(LogRecord(f'{rule}:property-non-scalar', (),
                                 f'{UNSUPPORTED}record: property {name!r:.30} carries '
                                 f'{", ".join(sorted(set(others)))}, which is not modelled', LOST))
        if name is None:
            log.append(LogRecord(f'{rule}:property-no-key', (),
                                 f'record: a <property> with neither dictRef nor title, dropped', LOST))
            continue
        if not scalars:
            continue
        if len(scalars) > 1:
            log.append(LogRecord(f'{rule}:property-multiple-scalars', (),
                                 f'{UNSUPPORTED}record: property {name!r:.30} has {len(scalars)} scalars and '
                                 f'one value is stored; several per property is not modelled', LOST))
        store[name] = text_of(scalars[0])
    if store:
        # Merged rather than assigned: a document may write several `<propertyList>` elements.
        record.ctab.meta.update(store)
    return True


def emit_properties(record, node, ns, *, names=('dictRef',)):
    """``<propertyList>`` for whatever ``record.ctab.meta`` holds, or nothing for an empty one.

    `names` is the attribute or attributes the field name goes out under: CML writes ``dictRef``, MRV
    both, which is what Marvin 25.1.3 writes and what its reader looks for.  Values go out marked for
    :func:`serialize` to write as CDATA -- see :func:`cdata` for the measurement behind that.
    """
    properties = record.ctab.meta
    if not properties:
        return
    plist = SubElement(node, _qualify('propertyList', ns))
    for name, value in properties.items():
        prop = SubElement(plist, _qualify('property', ns))
        for attribute in names:
            prop.set(attribute, str(name))
        scalar = SubElement(prop, _qualify('scalar', ns))
        scalar.text = cdata('' if value is None else str(value))


#: The pair :func:`cdata` wraps a value in and :func:`serialize` turns into a CDATA section.  U+0001 and
#: U+0002 are illegal in XML 1.0 -- as text, escaped, or as a numeric reference -- so no document can
#: carry either and no value can be mistaken for the marker; :func:`cdata` drops them from a value that
#: somehow holds one.
_CDATA_IN, _CDATA_OUT = '\x01', '\x02'
_CDATA_SPAN = re_compile(f'{_CDATA_IN}(.*?){_CDATA_OUT}', DOTALL)

#: What ``ElementTree`` escapes in element text, and how to put each back.  ``&amp;`` last: the others
#: reintroduce no ampersand, and reversing it first would decode ``&amp;lt;`` into ``<``.
_UNESCAPE = (('&lt;', '<'), ('&gt;', '>'), ('&#13;', '\r'), ('&amp;', '&'))


def cdata(text):
    """`text`, marked so that :func:`serialize` writes it as a CDATA section.

    A round trip through Marvin 25.1.3 measures the difference: a ``<scalar>`` holding ``R at C2`` as
    element text is read back as ``R``, the same value inside a CDATA section is read back whole, and a
    numeric character reference in that position is read back as its own source text.  A CDATA section
    and escaped text are one document to a conforming parser, so the section costs nothing to write.
    """
    return f'{_CDATA_IN}{text.replace(_CDATA_IN, "").replace(_CDATA_OUT, "")}{_CDATA_OUT}'


def _unwrap_cdata(match):
    body = match.group(1)
    for escaped, raw in _UNESCAPE:
        body = body.replace(escaped, raw)
    # `]]>` cannot appear inside a section, so a value holding one is written as two sections; that is
    # the sequence's only spelling and every parser rejoins them into the one value.
    return f'<![CDATA[{body.replace("]]>", "]]]]><![CDATA[>")}]]>'


def serialize(root, indent='  '):
    """`root` as an XML document: the declaration, the tree, and every marked value as a CDATA section.

    `indent` is the pretty-printing step; ``None`` writes one line.  Indented by default because these
    files are read by people at least as often as by programs, and an indented document diffs.
    """
    if indent is not None:
        ETIndent(root, space=indent)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            + _CDATA_SPAN.sub(_unwrap_cdata, ETToString(root, encoding='unicode')))


def molecule_nodes(root, dial):
    """Every ``<molecule>`` in `root`, outermost only, in document order.

    A molecule nested inside another is CML's way of writing an assembly, so reading both would double
    every atom.  The nesting is reported by the molecule walker, which has no handler for the child.
    """
    return outermost(root, dial.tags.molecule)


def read_document(root, dial, log):
    """Every molecule under `root` as a list of :class:`Record`.

    :attr:`Dialect.document` runs first, so a line about what the *document* structure loses precedes
    the lines about the molecules in it -- the order a caller reads the log in.
    """
    if dial.document is not None:
        dial.document(root, log)
    return [read_molecule(node, dial, log) for node in molecule_nodes(root, dial)]


def parse_records(source, dial, *, log=None, **kwargs):
    """Every ``<molecule>`` in `source` as a list of :class:`Record`, read as `dial` spells it.

    The body of every named ``parse_`` reader, which takes its dialect as an argument rather than
    sniffing: :func:`~._cml.parse_cml` knows it is the CML reader.  `kwargs` reach
    :func:`~._tree.parse_xml`, so no keyword above can bypass the entity and depth policy.
    """
    out = [] if log is None else log
    root = parse_xml(source, log=out, **kwargs)
    return read_document(root, dial, out)


def build_records(records, log, ignore_stereo=False):
    """Build `records`, extending `log` with what each build had to say.

    A record that is merely wrong is built and logged; one that cannot be built raises rather than
    shortening the list, a short list being indistinguishable from a file with fewer molecules in it.
    """
    molecules = []
    for record in records:
        # `Ctab.build` opens its log with a copy of `ctab.log`, which :func:`read_molecule` has already
        # given the caller; only what the build itself added is new here.  Sliced by POSITION and not
        # filtered by string: two molecules of one document legitimately log the same sentence.
        already = len(record.ctab.log)
        mol, _, record_log = record.ctab.build(ignore_stereo=ignore_stereo)
        log.extend(record_log[already:])
        molecules.append(mol)
    return molecules


def read_molecules(source, dial, *, log=None, ignore_stereo=False, **kwargs):
    """Every molecule in `source`, built, read as `dial` spells it.

    The body of every named ``read_`` reader, and the pair of :func:`parse_records`.
    """
    out = [] if log is None else log
    return build_records(parse_records(source, dial, log=out, **kwargs), out, ignore_stereo)


def parse_xml_document(source, *, log=None, **kwargs):
    """Every molecule in `source` as a list of :class:`Record`, with the dialect chosen by the file.

    The dialect-agnostic entry point, and :func:`sniff`'s only caller.  `kwargs` reach
    :func:`~._tree.parse_xml`, so ``engine``, ``max_depth`` and ``allow_dtd`` behave as on the named
    readers.  Prefer a named reader when the format is known; this is for a file whose vocabulary the
    caller has not been told, and it says in the log when it fell back to a more general dialect.
    """
    out = [] if log is None else log
    root = parse_xml(source, log=out, **kwargs)
    return read_document(root, sniff(root, out), out)


def read_xml(source, *, log=None, ignore_stereo=False, **kwargs):
    """Every molecule in `source`, built, with the dialect chosen by the file.

    :func:`parse_xml_document` plus the build.  Not :func:`read_molecules`, because the dialect is a
    property of the parsed root rather than an argument.
    """
    out = [] if log is None else log
    return build_records(parse_xml_document(source, log=out, **kwargs), out, ignore_stereo)


def read_molecule(node, dial, log):
    """One ``<molecule>`` element as a :class:`Record`, table-driven throughout.

    THE WALK WRITES TO THE RECORD'S OWN ``Ctab.log``, not to `log` directly, and `log` gets a copy at
    the end.  A document holds many ``<molecule>`` elements and one flat list cannot say which of them
    a line is about; the record's own list is what :meth:`~chython.formats.ctfile._ctab.Ctab.build`
    folds onto ``mol.log``, so molecule 3's lines end up on molecule 3.
    """
    record = Record()
    # Handed over here, the one function every read path goes through, rather than in each dialect's
    # `finish`: `Ctab.build` composes its unknown-count advice from it.
    record.ctab.channels = dial.channels
    own = record.ctab.log
    tags = dial.tags
    # `record` rather than `molecule` as the location word, here and in the two array walkers below: the
    # log-prefix convention has a closed set of leading tokens and `molecule`, `atomArray` and
    # `bondArray` are not in it.  The element name still appears inside the message.
    apply_fields(node, dial.index('molecule'), record.ctab, own, 'record',
                 dial.molecule_ignored, (dial.atom_id,), spill=record.extras)

    for child in node:
        tag = local(child.tag)
        if tag == tags.atom_array:
            _read_atom_array(child, record, dial, own)
        elif tag == tags.bond_array:
            _read_bond_array(child, record, dial, own)
        elif dial.molecule_child is not None and dial.molecule_child(child, record, own):
            continue
        elif tag in dial.molecule_children_ignored:
            continue
        else:
            own.append(LogRecord('xml:element-not-modelled', (),
                                 f'{UNSUPPORTED}record: <{tag}> in <{tags.molecule}> is not modelled', LOST))

    if dial.finish is not None:
        dial.finish(record, own)
    log.extend(own)
    return record


def _read_atom_array(node, record, dial, log):
    tags = dial.tags
    table = dial.index('atom')
    if dial.atom_array_hook is not None and dial.atom_array_hook(node, record, log):
        return
    for child in node:
        tag = local(child.tag)
        if tag != tags.atom:
            log.append(LogRecord('xml:element-not-modelled', (),
                                 f'{UNSUPPORTED}atom: <{tag}> in <{tags.atom_array}> is not modelled', LOST))
            continue
        atom = CtabAtom()
        position = record.add_atom(atom, child.get(dial.atom_id), log)
        atom.file_index = position + 1
        where = f'atom {record.ids[position]}'
        apply_fields(child, table, atom, log, where, dial.atom_ignored, (dial.atom_id,),
                     spill=record.atom_extras[position])
        for grandchild in child:
            if dial.atom_child is None or not dial.atom_child(grandchild, record, position, log):
                log.append(LogRecord('xml:element-not-modelled', (),
                                     f'{UNSUPPORTED}{where}: <{local(grandchild.tag)}> is not modelled', LOST))


def _read_bond_array(node, record, dial, log):
    tags = dial.tags
    table = dial.index('bond')
    if dial.bond_array_hook is not None and dial.bond_array_hook(node, record, log):
        return
    for child in node:
        tag = local(child.tag)
        if tag != tags.bond:
            log.append(LogRecord('xml:element-not-modelled', (),
                                 f'{UNSUPPORTED}bond: <{tag}> in <{tags.bond_array}> is not modelled', LOST))
            continue
        refs = next((child.get(name) for name in dial.bond_refs if child.get(name)), None)
        ident = child.get(dial.atom_id) or f'b{len(record.ctab.bonds) + 1}'
        if refs is None:
            # CML's two-attribute endpoint spelling (`atomRef1`/`atomRef2`), legal on a `<bond>` element
            # as well as in the array form -- so it lives here and not in a dialect's array hook, or the
            # element form would drop a bond whose endpoints it holds.
            pair = [child.get(name) for name in dial.bond_ref_pair]
            if len(pair) != 2 or not all(pair):
                log.append(LogRecord('xml:bond-no-refs', (), f'bond {ident}: no {dial.bond_refs[0]}, dropped', LOST))
                continue
            names = pair
        else:
            names = refs.split()
        if len(names) != 2:
            log.append(LogRecord('xml:bond-bad-refs', (),
                                 f'bond {ident}: {dial.bond_refs[0]}={refs!r:.40} does not name two atoms, '
                                 f'dropped', LOST))
            continue
        try:
            a, b = (record.index_of[name] for name in names)
        except KeyError as e:
            log.append(LogRecord(
                'xml:bond-unknown-atom', (), f'bond {ident}: references unknown atom {e.args[0]!r}, dropped', LOST))
            continue
        bond = CtabBond(a, b)
        position = record.add_bond(bond)
        where = f'bond {ident}'
        apply_fields(child, table, bond, log, where, dial.bond_ignored,
                     (dial.atom_id,) + tuple(dial.bond_refs) + tuple(dial.bond_ref_pair),
                     spill=record.bond_extras[position])
        for grandchild in child:
            if dial.bond_child is None or not dial.bond_child(grandchild, record, position, log):
                log.append(LogRecord('xml:element-not-modelled', (),
                                     f'{UNSUPPORTED}{where}: <{local(grandchild.tag)}> is not modelled', LOST))


# writing

def _qualify(tag, ns):
    return f'{{{ns}}}{tag}' if ns else tag


def write_molecule(record, dial, log, parent=None):
    """`record` as a ``<molecule>`` element in `dial`'s vocabulary.  Returns the element.

    Driven by the same table as the reader, through :attr:`Field.encode`, so a bond order cannot be read
    as ``A`` and written as ``4``.  An empty ``<bondArray/>`` is not written for a single-atom molecule:
    real files omit it, and an empty array is a statement a reader has to decide about.
    """
    ns = dial.ns
    ctab = record.ctab
    element = (Element(_qualify(dial.tags.molecule, ns)) if parent is None
               else SubElement(parent, _qualify(dial.tags.molecule, ns)))
    _encode_into(element, dial.molecule_fields, ctab)
    if dial.emit_molecule is not None:
        dial.emit_molecule(record, element, log)

    atom_array = SubElement(element, _qualify(dial.tags.atom_array, ns))
    for position, atom in enumerate(ctab.atoms):
        node = SubElement(atom_array, _qualify(dial.tags.atom, ns))
        node.set(dial.atom_id, record.ids[position])
        _encode_into(node, dial.atom_fields, atom, spill=record.atom_extras[position])
        if dial.emit_atom is not None:
            dial.emit_atom(record, position, node, log)

    if ctab.bonds:
        bond_array = SubElement(element, _qualify(dial.tags.bond_array, ns))
        for position, bond in enumerate(ctab.bonds):
            node = SubElement(bond_array, _qualify(dial.tags.bond, ns))
            node.set(dial.atom_id, f'b{position + 1}')
            node.set(dial.bond_refs[0], f'{record.ids[bond.a]} {record.ids[bond.b]}')
            _encode_into(node, dial.bond_fields, bond, spill=record.bond_extras[position])
            if dial.emit_bond is not None:
                dial.emit_bond(record, position, node, log)
    return element


def _encode_into(node, rows, source, spill=None):
    """Set every attribute `rows` encodes to something other than ``None``, in declaration order.

    Declaration order, not the dict's, so two runs over one molecule produce byte-identical output.
    `rows` is a sequence and never looked up by attribute, so there is no index for a caller to build per
    atom.  :attr:`Field.write_slot` wins over :attr:`Field.slot` here and only here -- the writer wants
    the value :attr:`Dialect.finish` interpreted, not the raw one the reader parked.  `spill` is the
    per-item dict of the read path, and is the fallback when the slot names nothing on `source`; without
    it a quantity the intermediate cannot hold would have to borrow another row's slot.
    """
    for row in rows:
        if row.encode is None:
            continue
        slot = row.write_slot or row.slot
        value = getattr(source, slot, None)
        if value is None and spill is not None:
            value = spill.get(slot)
        text = row.encode(value)
        if text is not None:
            node.set(row.attribute, text)
