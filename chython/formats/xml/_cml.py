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
"""CML -- Chemical Markup Language -- as a table over :mod:`._dialect`, read and written.

Three spellings are read, all routed through the one table: the CML 2/3 element form (``<atom
id="a1" elementType="C" x2="0"/>``, the only form written), the array form (``<atomArray atomID="a1 a2"
elementType="C O"/>``, one whitespace-separated column per attribute) and CML 1 (``<atom><string
builtin="elementType">C</string></atom>``, with ``<stringArray builtin=...>`` as its array form).
Stereo is stated onto ``Ctab``, so the core interprets it exactly as it does an MDL record's.
"""

from xml.etree.ElementTree import Element, SubElement

from ._dialect import (Dialect, Field, NotModelled, Record, Tags, apply_fields, choose_coordinates,
                       decode_float, decode_int, emit_coordinates, emit_properties, encode_element,
                       encode_from_table, encode_nonzero, encode_stated, outermost, parse_records,
                       read_molecules, read_properties, resolve_element, serialize, synthetic_node,
                       write_molecule, xml_text)
from ._tree import local, text_of
from ..ctfile import Ctab, CtabAtom, CtabBond
from ..ctfile._ctab import WEDGE_FROM_V2000, WEDGE_TO_V2000
from ..ctfile._hydrogens import StatedChannels
from ..ctfile._sgroup import UNSUPPORTED, resolve_output
from ...core.wedge import SU_TETRA, _permutation_is_odd, cis_trans_for_write, wedges_for_write
from ...core import LOST, LogRecord, MoleculeContainer, REPAIRED, WEDGE_DOWN, WEDGE_NONE, WEDGE_UP


__all__ = ['CML', 'CML_NS', 'parse_cml', 'read_cml', 'record_from_molecule', 'write_cml',
           'write_cml_element']


#: The namespace this writer emits.  CML 3's, which is also what CML 2 documents declare in practice.
CML_NS = 'http://www.xml-cml.org/schema'

#: Every namespace URI seen on a real CML document: the CML 2 core and the CML 1 DTD-era URIs are both
#: still in deposited data, and matching only the current one sends those files down
#: :func:`~._dialect.sniff`'s fallback with a log line for no reason.
CML_NAMESPACES = frozenset((
    CML_NS,
    'http://www.xml-cml.org/schema/cml2/core',
    'http://www.xml-cml.org/schema/cml3/core',
    'http://www.xml-cml.org/dict/cml',
    'http://www.xml-cml.org',
))

#: CML bond orders a molecule can hold, both spellings of each.  ``A`` is aromatic and is stored as
#: chython's order 4 as stated, not kekulised -- the CTfile reader's posture to bond type 4.
_ORDERS = {'1': 1, 'S': 1, '2': 2, 'D': 2, '3': 3, 'T': 3, 'A': 4}

#: Orders CML has and a molecule does not.  A delocalised fractional order has no field, and ``unknown``
#: is the file declining to say, which is not single and must not be read as one.
_PARTIAL_ORDERS = frozenset(('partial01', 'partial12', 'partial23', 'unknown', 'other'))

#: Order back out.  A separate table because `_ORDERS` is many-to-one -- reading accepts ``S`` and ``1``
#: -- and inverting a many-to-one dict silently picks whichever key came last.
_ORDERS_OUT = {1: '1', 2: '2', 3: '3', 4: 'A'}

#: Element symbols CML uses for a SET of elements: MDL's ``A``/``Q``/``X`` query types, which a CML file
#: converted from a molfile inherits.  A molecule can hold none of them, so each earns a refusal naming a
#: query reader rather than a carbon.  ``R``/``R1``.., ``Du``/``Dummy``/``*`` are NOT here: each names one
#: atom, so each is the marker -- `resolve_element` takes the R family and every other unresolvable
#: symbol reads as a label on the marker.
_PSEUDO = {'A': 'the any-atom query type', 'Q': 'the any-heteroatom query type',
           'X': 'the halogen query type', 'AH': 'the any-atom-or-hydrogen query type',
           'QH': 'the any-heteroatom-or-hydrogen query type', 'M': 'the any-metal query type'}


# codecs -- one per table row, and only the ones that read CML's own vocabulary; an integer, a float, a
# coordinate and an element symbol are the same in every dialect and come from `_dialect`.

def _order_in(text):
    """CML ``order`` as a chython bond order, or :class:`~._dialect.NotModelled`."""
    key = text.strip()
    if key in _ORDERS:
        return _ORDERS[key]
    if key.upper() in _ORDERS:
        return _ORDERS[key.upper()]
    if key.lower() in _PARTIAL_ORDERS:
        raise NotModelled(f'bond order {key!r} has no representation in a molecule; read as single')
    raise ValueError(f'{key!r} is not a CML bond order')


#: A bond order as CML text, or ``None`` to leave the attribute off.  Order 8, the dative bond, gets
#: ``None``: CML has no coordination order, and :func:`record_from_molecule` says so in the log.
_order_out = encode_from_table(_ORDERS_OUT)


def _spin_in(text):
    """CML ``spinMultiplicity`` as-stated, for :func:`_finish` to turn into a radical bit."""
    return decode_int(text)


def _spin_out(radical):
    """A radical bit as a spin multiplicity.  Doublet, because one bit is one unpaired electron."""
    return '2' if radical else None


# the table

_MOLECULE_FIELDS = (
    Field('title', 'title', lambda t: t, lambda t: t or None),
)

_ATOM_FIELDS = (
    # `element_type` is no `CtabAtom` slot, so the engine spills it and `_finish` resolves it: resolution
    # must tell "stated nothing" from "stated carbon", and `CtabAtom.element` defaults to carbon.
    # `write_slot` sends the writer to the real slot.
    Field('elementType', 'element_type', lambda t: t.strip(), encode_element, write_slot='element'),
    Field('formalCharge', 'charge', decode_int, encode_nonzero),
    # CML 2 spells the mass number `isotope` and CML 3 spells it `isotopeNumber`; the schema admits both
    # and readers differ over which they take, so both go out on an atom that has one.  MEASURED, one
    # document per spelling: `isotope` is read by Indigo 1.45 and Marvin 25.1.3, `isotopeNumber` by CDK
    # 2.12 and Indigo, both by all three.
    Field('isotopeNumber', 'isotope', decode_int, encode_nonzero),
    Field('isotope', 'isotope', decode_int, encode_nonzero),
    Field('spinMultiplicity', 'spin', _spin_in, _spin_out, write_slot='radical'),
    # CML's `hydrogenCount` is the TOTAL, explicit neighbours included -- unlike MRV's, and unlike
    # `CtabAtom.stated_h`, which is the implicit count `calc_implicit` takes.  `_finish` subtracts.  The
    # slot holds the total in both directions, in `atom_extras` since `CtabAtom` has no field for it;
    # writing it from `stated_h` instead costs one hydrogen per round trip for every drawn hydrogen.
    Field('hydrogenCount', 'hydrogen_total', decode_int, encode_stated),
    # Coordinates are read by the table and written by `_emit_atom`: *which* set to write is a question
    # about the record (see `emit_coordinates`) and a row is handed one atom's one value.  A row encoding
    # `x2` unconditionally writes `x2="0.0000"` on a record with no layout, a drawing nobody made.
    Field('x2', 'x', decode_float, None),
    Field('y2', 'y', decode_float, None),
    Field('x3', 'x3', decode_float, None),
    Field('y3', 'y3', decode_float, None),
    Field('z3', 'z3', decode_float, None),
)

_BOND_FIELDS = (
    Field('order', 'order', _order_in, _order_out),
)

#: Attributes carrying no chemistry, so no `unsupported:` line.  Each entry claims that a reader
#: honouring it would build the same molecule.
_MOLECULE_IGNORED = frozenset((
    'convention',        # names the dictionary a `<bondStereo>` or `<scalar>` is read against
    'ref', 'role',       # document-structure pointers
    'formula',           # derivable from the atoms, and not authoritative when it disagrees
))
_ATOM_IGNORED = frozenset(('ref',))
#: `atomRef1`/`atomRef2` are deliberately not here: they are the endpoints, which is structural rather
#: than ignorable.  `bond_ref_pair` reads them in both forms and `structural` is what silences them.
_BOND_IGNORED = frozenset(('ref',))

#: Children of `<molecule>` with no chemistry in them.
_MOLECULE_CHILDREN_IGNORED = frozenset(('metadataList', 'metadata'))

#: The array form's column names mapped onto the element form's attribute names.  Only the identity
#: column differs, which is what lets one table serve both forms.
_ARRAY_ALIASES = {'atomID': 'id', 'bondID': 'id', 'atomId': 'id', 'bondId': 'id'}


def _columns(node, log, where, member):
    """``{name: [values]}`` for the array form, from attributes and CML 1 ``<*Array>`` children.

    Returns ``None`` when `node` is not in the array form, which must stay distinguishable from "array
    form with no columns".  Ragged columns are truncated to the shortest and reported; dropping the array
    instead loses every atom in the file over one bad column.

    The form is decided by the children and never by the attributes, `member` naming the per-item element:
    an ``<atomArray>`` holding ``<atom>`` children is the element form whatever attributes it carries.
    Attributes cannot decide it -- ``id`` and ``dictRef`` are legal CML *global* attributes on every
    element.  In a node carrying both the children win, and the dropped columns are named.
    """
    out = {}
    unnamed = []
    for name, value in node.attrib.items():
        plain = local(name)
        key = _ARRAY_ALIASES.get(plain, plain)
        # An unaliased `id` is the element's own name, not a column: the array form spells its identity
        # column `atomID`/`bondID`, and a bare `id` on an `<atomArray>` is the global attribute every CML
        # element may carry.  Read as a column it invents an atom out of `<atomArray id="aa1"/>`.
        if plain == 'id' or key in ('title', 'convention', 'ref'):
            continue
        out[key] = value.split()
    members = sum(1 for child in node if local(child.tag) == member)
    if members:
        if out:
            log.append(LogRecord('cml:array-columns-beside-elements', (),
                                 f'{where}: {len(out)} array column(s) ({", ".join(sorted(out))}) beside '
                                 f'{members} <{member}> element(s); the elements are read and the columns '
                                 f'dropped', LOST))
        return None
    for child in node:
        tag = local(child.tag)
        if not tag.endswith('Array'):
            continue
        builtin = child.get('builtin')
        if builtin is None:
            # Not logged yet: this function may still decline the node, and then the engine walks the
            # same children and reports them itself.  Logging here as well gives one construct two lines.
            unnamed.append(tag)
            continue
        out[_ARRAY_ALIASES.get(builtin, builtin)] = text_of(child).split()
    if not out:
        return None
    for tag in unnamed:
        log.append(LogRecord('cml:array-child-no-builtin', (),
                             f'{UNSUPPORTED}{where}: <{tag}> without a builtin attribute is not modelled',
                             LOST))
    width = min(len(v) for v in out.values())
    ragged = {k: len(v) for k, v in out.items() if len(v) != width}
    if ragged:
        log.append(LogRecord('cml:array-columns-ragged', (),
                             f'{where}: array columns have different lengths ({ragged}); truncated to {width}',
                             REPAIRED))
    return {k: v[:width] for k, v in out.items()}, width


def _atom_array(node, record, log):
    """The array form of ``<atomArray>``.  ``True`` when it claimed the node."""
    found = _columns(node, log, 'atom', CML.tags.atom)
    if found is None:
        return False
    columns, width = found
    ids = columns.pop('id', None)
    table = CML.index('atom')
    for i in range(width):
        atom = CtabAtom()
        position = record.add_atom(atom, ids[i] if ids else None, log)
        atom.file_index = position + 1
        where = f'atom {record.ids[position]}'
        apply_fields(synthetic_node('atom', {k: v[i] for k, v in columns.items()}), table, atom, log,
                     where, _ATOM_IGNORED, ('id',), spill=record.atom_extras[position])
    return True


def _bond_array(node, record, log):
    """The array form of ``<bondArray>``.  ``True`` when it claimed the node.

    The endpoints are two columns here -- ``atomRef1`` and ``atomRef2`` -- rather than the element form's
    single ``atomRefs2`` holding a pair, the one place the two forms differ in more than a name.
    """
    found = _columns(node, log, 'bond', CML.tags.bond)
    if found is None:
        return False
    columns, width = found
    ids = columns.pop('id', None)
    first = columns.get('atomRef1')
    second = columns.get('atomRef2')
    refs = columns.get('atomRefs2')
    table = CML.index('bond')
    for i in range(width):
        ident = ids[i] if ids else f'b{i + 1}'
        if first is not None and second is not None:
            names = [first[i], second[i]]
        elif refs is not None:
            names = refs[i].split()
        else:
            log.append(LogRecord('cml:bond-no-endpoints', (),
                                 f'bond {ident}: the bond array names no endpoints, dropped', LOST))
            continue
        if len(names) != 2:
            log.append(LogRecord('cml:bond-bad-endpoint-count', (),
                                 f'bond {ident}: {names} does not name two atoms, dropped', LOST))
            continue
        try:
            a, b = (record.index_of[name] for name in names)
        except KeyError as e:
            log.append(LogRecord('cml:bond-unknown-atom', (),
                                 f'bond {ident}: references unknown atom {e.args[0]!r}, dropped', LOST))
            continue
        bond = CtabBond(a, b)
        position = record.add_bond(bond)
        apply_fields(synthetic_node('bond', {k: v[i] for k, v in columns.items()}), table, bond, log,
                     f'bond {ident}', _BOND_IGNORED,
                     ('id',) + CML.bond_refs + CML.bond_ref_pair,
                     spill=record.bond_extras[position])
    return True


def _atom_child(node, record, position, log):
    """A child element of ``<atom>``: CML 1's ``<string builtin=...>``, and ``<atomParity>``."""
    tag = local(node.tag)
    if tag == 'atomParity':
        return _atom_parity(node, record, position, log)
    builtin = node.get('builtin')
    if builtin is None or tag not in ('string', 'float', 'integer'):
        return False
    atom = record.ctab.atoms[position]
    key = _ARRAY_ALIASES.get(builtin, builtin)
    if key == 'id':
        return True  # consumed by the engine from the attribute; a builtin id restates it
    apply_fields(synthetic_node('atom', {key: text_of(node)}), CML.index('atom'), atom, log,
                 f'atom {record.ids[position]}', _ATOM_IGNORED, ('id',),
                 spill=record.atom_extras[position])
    return True


def _atom_parity(node, record, position, log):
    """Park an ``<atomParity>`` for :func:`_resolve_parities`.  Always ``True``.

    Parked and not resolved: an ``<atomParity>`` is a child of ``<atom>``, read while the atom array is
    still being walked, so ``atomRefs4`` on an early atom names atoms that do not exist yet -- which on a
    real file is nearly every descriptor, a stereocentre normally being drawn before its substituents.
    """
    ident = record.ids[position]
    # Which spelling was read is parked with the references: `_resolve_parities` reports a wrong-length
    # quadruple, and naming `atomRefs4` for a file that wrote `atomRefs` sends the reader looking for an
    # attribute their document does not contain.
    refs, attribute = (node.get('atomRefs4') or '').split(), 'atomRefs4'
    if not refs and node.get('atomRefs'):
        refs, attribute = node.get('atomRefs').split(), 'atomRefs'
    value = text_of(node)
    try:
        sign = int(float(value))
    except ValueError:
        log.append(LogRecord('cml:atom-parity-bad-value', (),
                             f'stereo: atom {ident}: atomParity value {value!r:.20} is not a number, ignored',
                             LOST))
        return True
    if sign:  # `0` is CML for "no configuration", which is what the atom already says
        record.atom_extras[position]['parity'] = (sign, refs, attribute)
    return True


def _resolve_parities(record, log):
    """Every parked ``<atomParity>`` as ``CtabAtom.parity``, once every atom id is known.

    CML measures its sign against the order ``atomRefs4`` lists; ``CtabAtom.parity`` against ascending
    atom-block position, so the sign flips iff the permutation between the two is odd.  Naming the centre
    itself is CML's spelling for an implicit hydrogen or lone pair, and it ranks last at ``len(atoms) + 1``
    -- the rank ``stated_parity`` gives its own undrawn direction, so the two cancel on a round trip.
    Exactly four references are required: a parity depends on *where* a missing direction sits.  An
    unresolvable frame costs the descriptor, not the atom.
    """
    high = len(record.ctab.atoms) + 1
    for position, spill in enumerate(record.atom_extras):
        if 'parity' not in spill:
            continue
        sign, refs, attribute = spill['parity']
        ident = record.ids[position]
        if len(refs) != 4:
            log.append(LogRecord('cml:atom-parity-wrong-ref-count', (),
                                 f'stereo: atom {ident}: {attribute} names {len(refs)} atoms, not four; ignored, '
                                 f'since a tetrahedral parity is defined only against four references',
                                 LOST))
            continue
        keys = []
        for name in refs:
            if name == ident:
                keys.append(high)
                continue
            try:
                keys.append(record.index_of[name])
            except KeyError:
                log.append(LogRecord('cml:atom-parity-unknown-atom', (),
                                     f'stereo: atom {ident}: atomParity references unknown atom {name!r}, '
                                     f'ignored', LOST))
                keys = None
                break
        if keys is None:
            continue
        if keys.count(high) > 1:
            log.append(LogRecord('cml:atom-parity-phantom-twice', (),
                                 f'stereo: atom {ident}: atomParity names the phantom direction twice, ignored',
                                 LOST))
            continue
        record.ctab.atoms[position].parity = _parity_field(sign, keys)


def _parity_field(sign, keys):
    """CML's ``<atomParity>`` sign, stated in the frame `keys`, as a ``CtabAtom.parity``.

    Calibrated, not read off the specification: ``<atomParity atomRefs4="a b c d">1</atomParity>`` is MDL
    atom-parity field 1 measured in the ``atomRefs4`` frame, with the centre's own id -- or a missing
    fourth reference -- ranking last.  `keys` is that frame as atom-block positions, and the field flips
    iff `keys` is an odd permutation.  The formula is its own inverse, so
    :func:`~chython.core.wedge.stated_parity` translates back over the same `keys`.
    """
    base = 1 if sign > 0 else 2
    return (3 - base) if _permutation_is_odd(keys) else base


def _parity_refs(record, position):
    """The four references and the sign an ``<atomParity>`` states for the atom at `position`.

    ``None`` when the atom states no configuration, or when its bonds cannot name four directions.  Reads
    ``CtabAtom.parity`` and nothing else, so a record from a file and one from a molecule are written from
    the same field.  The frame is ascending atom-block position -- the frame the field is already measured
    in -- so the sign *is* the field.  A three-bonded centre names itself last; anything else is a frame a
    quadruple cannot describe and the caller reports it.
    """
    parity = record.ctab.atoms[position].parity
    if parity not in (1, 2):
        return None
    refs = sorted({b.b if b.a == position else b.a
                   for b in record.ctab.bonds if position in (b.a, b.b)})
    if len(refs) == 3:
        refs.append(position)
    elif len(refs) != 4:
        return None
    return refs, 1 if parity == 1 else -1


def _bond_child(node, record, position, log):
    """``<bondStereo>`` as ``CtabBond.wedge``, plus CML 1's ``<string builtin=...>`` on a bond."""
    tag = local(node.tag)
    bond = record.ctab.bonds[position]
    ident = f'b{position + 1}'
    if tag != 'bondStereo':
        # CML 1 spells every scalar as a typed child with a `builtin` name, a bond's `order` being the one
        # that matters; routed through the same table, so `A` and `partial12` mean on a child what they
        # mean on an attribute.
        builtin = node.get('builtin')
        if builtin is None or tag not in ('string', 'float', 'integer'):
            return False
        apply_fields(synthetic_node('bond', {builtin: text_of(node)}), CML.index('bond'), bond, log,
                     f'bond {ident}', _BOND_IGNORED, ('id', 'atomRefs2'),
                     spill=record.bond_extras[position])
        return True
    convention = node.get('convention', '')
    if convention and convention.upper() != 'MDL':
        # A `convention` names the dictionary its content is defined in, and MDL's is the only one read
        # here: `W` is a wedge in CML's own vocabulary and could mean anything in somebody else's.
        # Checked before the MDL branch, so a foreign `conventionValue` is not decoded as a CTfile code
        # either -- that number belongs to that dictionary too.
        log.append(LogRecord('cml:bond-stereo-unknown-convention', (),
                             f'{UNSUPPORTED}bond {ident}: bondStereo convention {convention!r:.30} names a '
                             f'dictionary this reader has not read; nothing applied', LOST))
        return True
    if convention or node.get('conventionValue') is not None:
        # ChemAxon and every molfile-derived converter write the CTfile bond-stereo number straight
        # through, dressed as a dictionary reference.  Decoded with the CTfile reader's own table, so what
        # 1, 4 and 6 mean is stated once in this tree.
        raw = node.get('conventionValue', '')
        try:
            code = int(raw)
        except ValueError:
            # One line on purpose: `test_unmodelled_constructs_are_prefixed` reads the first string
            # literal of an `append(...)`, so a marker word carried into a second f-string fragment is a
            # message the convention check never sees.
            log.append(LogRecord('cml:bond-stereo-convention-value-bad', (),
                                 f'bond {ident}: bondStereo conventionValue {raw!r:.20} is not a number, ignored',
                                 LOST))
            return True
        if code in WEDGE_FROM_V2000:
            bond.wedge = WEDGE_FROM_V2000[code]
        elif code == 3:
            # MDL 3 on a double bond is "cis or trans, unknown which": a real construct with no field in a
            # molecule, an unset configuration and an explicitly unknown one being one arena value.
            log.append(LogRecord('cml:bond-stereo-cis-or-trans-unknown', (),
                                 f'{UNSUPPORTED}bond {ident}: "cis or trans, unknown which" is not modelled',
                                 LOST))
        else:
            log.append(LogRecord('cml:bond-stereo-convention-value-bad', (),
                                 f'bond {ident}: bondStereo conventionValue {code} is not 0, 1, 3, 4 or 6, ignored',
                                 LOST))
        return True

    text = text_of(node).upper()
    if text in ('W', 'WEDGE'):
        bond.wedge = WEDGE_UP
    elif text in ('H', 'HATCH', 'HASH'):
        bond.wedge = WEDGE_DOWN
    elif text in ('C', 'T'):
        # CML's own double-bond descriptor.  Recorded verbatim and resolved in `_finish` by
        # `_stated_cis_trans`, the frame being written as atom *names* and nothing guaranteeing that
        # `<atomArray>` precedes `<bondArray>`.  Ranking the letter against the drawing is the core's job,
        # in `assign_parities`, for every format at once.
        record.bond_extras[position]['cis_trans'] = (text, (node.get('atomRefs4') or '').split())
    elif not text:
        log.append(LogRecord('cml:bond-stereo-empty', (),
                             f'bond {ident}: empty bondStereo, ignored', LOST))
    else:
        log.append(LogRecord('cml:bond-stereo-unknown-letter', (),
                             f'{UNSUPPORTED}bond {ident}: bondStereo {text!r:.20} is not modelled', LOST))
    return True


def _molecule_child(node, record, log):
    """``<name>`` as the record's title, ``<propertyList>`` or a bare ``<scalar>`` as its data fields."""
    tag = local(node.tag)
    if tag in ('propertyList', 'property'):
        return read_properties(node, record, log, rule='cml')
    if tag == 'scalar':
        # A molecule-level `<scalar>` is the other spelling of a data field: CDK 2.12 writes one per
        # molecule property, `<scalar dictRef="cdk:molecularProperty" title="BATCH_ID">`.  Here `title`
        # is the field name and `dictRef` names the kind of property, so it is not read as the name.
        name = node.get('title')
        if not name:
            log.append(LogRecord('cml:molecule-scalar-no-title', (),
                                 f'{UNSUPPORTED}record: a molecule-level <scalar> states no title to '
                                 f'name the field by', LOST))
            return True
        record.ctab.meta[name] = text_of(node)
        return True
    if tag != 'name':
        return False
    text = text_of(node)
    if not text:
        return True
    if record.ctab.title and record.ctab.title != text:
        # Two names for one molecule: the first wins, a `title` attribute being read before any child, and
        # the second is named rather than dropped silently.  Prefixed, not a repair line -- CML allows a
        # molecule several names under different `convention` attributes, so the file is correct and the
        # one-title container is the limitation.
        log.append(LogRecord('cml:molecule-extra-name', (),
                             f'{UNSUPPORTED}record: <name>{text!r:.30} not stored; the molecule is already '
                             f'titled {record.ctab.title!r:.30}', LOST))
    else:
        record.ctab.title = text
    return True


def _document(root, log):
    """One ``unsupported: `` line per ``<reaction>``, with the molecule count of each role it names.

    A reaction's content *is* the roles, and the molecule walker reads straight through one, so the
    molecules are still returned flat and the counts are what lets a caller recover the record's shape
    from the log.  Roles are read off the file's own element names rather than mapped onto a fixed
    vocabulary: CML has more of them than a reaction has sides (``<substanceList>``, ``<spectator>``).
    Outermost reactions only, or a reaction nested in a scheme counts its molecules twice.
    """
    for node in outermost(root, 'reaction'):
        counts = {}
        _roles(node, 'unplaced', counts)
        ident = node.get('id') or '(unnamed)'
        if counts:
            parts = [f'{n} {role}' for role, n in counts.items()]
            what = f'{" and ".join([", ".join(parts[:-1]), parts[-1]] if len(parts) > 1 else parts)} ' \
                   f'molecule(s) read as a flat list'
        else:
            what = 'it holds no molecules'
        log.append(LogRecord('cml:reaction-roles', (),
                             f'{UNSUPPORTED}record: <reaction> {ident} roles are not modelled; {what}', LOST))


def _roles(node, role, counts):
    """Count the ``<molecule>`` descendants of `node` by the role element that encloses each.

    `role` is the nearest enclosing element's name with a trailing ``List`` stripped, so
    ``<productList><product>`` and the bare ``<productList><molecule>`` real files write count as one
    role.  The descent stops at a ``<molecule>``, so an assembly's nested parts count once, as they do to
    the walker that reads them.
    """
    for child in node:
        tag = local(child.tag)
        if tag == CML.tags.molecule:
            counts[role] = counts.get(role, 0) + 1
        else:
            _roles(child, tag[:-4] if tag.endswith('List') else tag, counts)


# finish -- everything that needs the whole molecule

def _finish(record, log):
    """Resolve everything that could not be decided one attribute at a time.

    1. the element, whose recoveries and refusals need the atom's name for the message;
    2. the radical bit, from a spin multiplicity whose meaning above a doublet has to be reported;
    3. the coordinates, 2D and 3D being two attribute sets and the choice a property of the record;
    4. the implicit hydrogen count, CML stating the total and the explicit hydrogens being neighbours;
    5. the atom parities, since ``atomRefs4`` may name atoms the walker has not reached yet;
    6. a cis/trans ``<bondStereo>``, whose frame is atom names and nothing orders the two arrays.
    """
    ctab = record.ctab
    atoms = ctab.atoms
    extras = record.atom_extras
    # Aggregated: a document stating no `elementType` anywhere states none on every atom, so a line each
    # is the whole atom list repeated.  One line, a count and the first atom it happened on.
    carbons, first_carbon = 0, None

    for position, (atom, spill) in enumerate(zip(atoms, extras)):
        where = f'atom {record.ids[position]}'
        if 'element_type' in spill:
            atom.element, isotope, label, r_index = resolve_element(spill['element_type'], where, log, _PSEUDO)
            if isotope and not atom.isotope:
                atom.isotope = isotope
            atom.r_index = r_index
            if label is not None:
                # The alias is where a display label goes, as in every other reader; CML has no alias
                # attribute of its own, so nothing can already be there to outrank it.
                ctab.aliases[position] = label
        else:
            # `<atom>` with no `elementType` is a file that did not say, and carbon is what every CML
            # reader assumes; said out loud rather than chosen silently.
            carbons += 1
            if first_carbon is None:
                first_carbon = where

        if 'spin' in spill:
            spin = spill['spin']
            atom.radical = spin > 1
            if spin > 2:
                log.append(LogRecord('cml:spin-high-multiplicity', (),
                                     f'{UNSUPPORTED}{where}: spin multiplicity {spin} is not modelled; read '
                                     f'as one radical centre', LOST))
            elif spin < 1:
                log.append(LogRecord('cml:spin-out-of-range', (),
                                     f'{where}: spin multiplicity {spin} is out of range, read as a singlet',
                                     REPAIRED))

    if carbons:
        log.append(LogRecord('cml:atom-no-element-type', (),
                             f'atom: {carbons} atom(s) with no elementType, read as carbon (first {first_carbon})',
                             REPAIRED))
    choose_coordinates(record, log)

    # The hydrogen total, after the elements resolve and so after "is this neighbour a hydrogen" has an
    # answer.  Counted over the bonds, the file's own total being the number being corrected.
    explicit = [0] * len(atoms)
    for bond in ctab.bonds:
        for end, other in ((bond.a, bond.b), (bond.b, bond.a)):
            if 0 <= other < len(atoms) and atoms[other].element == 'H' and 0 <= end < len(atoms):
                explicit[end] += 1
    for position, (atom, spill) in enumerate(zip(atoms, extras)):
        if 'hydrogen_total' not in spill:
            continue
        total = spill['hydrogen_total']
        implicit = total - explicit[position]
        if implicit < 0:
            log.append(LogRecord('cml:hydrogen-count-below-drawn', (),
                                 f'atom {record.ids[position]}: hydrogenCount {total} is below the '
                                 f'{explicit[position]} hydrogen(s) drawn on it; recomputed', REPAIRED))
            # The refused total leaves the spill, which is what the writer emits: keeping it re-publishes
            # a count this reader has just said it does not believe.
            del spill['hydrogen_total']
            continue
        atom.stated_h = implicit

    _resolve_parities(record, log)
    _stated_cis_trans(record, log)


def _stated_cis_trans(record, log):
    """Resolve every ``<bondStereo>C``/``T`` onto its ``CtabBond.configuration``.  Ranks nothing.

    This dialect owns the spelling only: a letter and a frame written as atom *names*, which only this
    record can turn into positions.  Whether the letter or the drawing wins is
    :func:`chython.core.wedge.stated_cis_trans`'s business, called from ``Ctab.build`` for every format.
    Deferred to ``_finish`` rather than done in ``_bond_child`` because the frame names atoms and nothing
    guarantees ``<atomArray>`` precedes ``<bondArray>``.  A missing ``atomRefs4`` is silent and normal --
    Marvin writes the letter bare, the frame being the drawing the record also carries.
    """
    bonds = record.ctab.bonds
    for position, spill in enumerate(record.bond_extras):
        if 'cis_trans' not in spill or position >= len(bonds):
            continue
        letter, refs = spill['cis_trans']
        ident = f'b{position + 1}'
        if refs and len(refs) != 4:
            log.append(LogRecord('cml:bond-stereo-wrong-ref-count', (),
                                 f'stereo: bond {ident}: cis/trans bondStereo names {len(refs)} atoms, not '
                                 f'four; not read', LOST))
            continue
        try:
            quad = tuple(record.index_of[name] for name in refs)
        except KeyError as e:
            log.append(LogRecord('cml:bond-stereo-unknown-atom', (),
                                 f'stereo: bond {ident}: cis/trans bondStereo references unknown atom '
                                 f'{e.args[0]!r}, not read', LOST))
            continue
        bonds[position].configuration = (letter, quad)


# writing -- the child elements a table cannot express

def _emit_molecule(record, node, log):
    """``<propertyList>`` for whatever ``record.ctab.meta`` holds.

    The record's title goes out as the ``title`` attribute of ``<molecule>``, by the table, so no
    ``<name>`` is written: both spellings read back the same, and writing both is what this reader's
    ``<name>`` handler reports as a contradiction.

    A property states its name in both ``dictRef`` and ``title``, which is what ``molconvert cml``
    writes for an SD field: a reader keying on either one gets the name.
    """
    emit_properties(record, node, CML_NS if node.tag.startswith('{') else '',
                    names=('dictRef', 'title'))


def _emit_atom(record, position, node, log):
    """The coordinates, and an ``<atomParity>`` for an atom that states a configuration."""
    emit_coordinates(record, record.ctab.atoms[position], node)
    if record.ctab.atoms[position].r_index:
        # `elementType="R"` is written by the table; the INDEX has no CML spelling, and `molconvert cml`
        # drops it the same way.  MRV keeps it, in `rgroupRef`.
        log.append(LogRecord('cml:r-index-not-written', (record.ids[position],),
                             f'{UNSUPPORTED}atom {record.ids[position]}: R index '
                             f'{record.ctab.atoms[position].r_index} is not written; CML spells an '
                             f'R-group label and not its number', LOST))
    frame = _parity_refs(record, position)
    if frame is None:
        if record.ctab.atoms[position].parity in (1, 2):
            log.append(LogRecord('cml:atom-parity-unwritable', (),
                                 f'stereo: atom {record.ids[position]}: configuration not written as an '
                                 f'atomParity; its bonds cannot name four directions', LOST))
        return
    refs, sign = frame
    child = SubElement(node, f'{{{CML_NS}}}atomParity' if node.tag.startswith('{') else 'atomParity')
    child.set('atomRefs4', ' '.join(record.ids[i] for i in refs))
    child.text = '1' if sign > 0 else '-1'


def _emit_bond(record, position, node, log):
    """``<bondStereo>`` for a wedged bond or a configured double bond, in the spellings this reader reads.

    ``W``/``H`` rather than ``convention="MDL"``: the plain letters are CML's own vocabulary and the MDL
    dictionary reference is a vendor accommodation accepted on the way in.  An ``either`` wedge goes out as
    the MDL number, CML having no letter for it.  One element name, two constructs, never both on one bond
    -- a wedge is on a single bond and a configuration on a double one.  ``atomRefs4`` is always written,
    unlike Marvin, because a document this writer produces may carry no drawing to be the frame.
    """
    bond = record.ctab.bonds[position]
    tag = f'{{{CML_NS}}}bondStereo' if node.tag.startswith('{') else 'bondStereo'
    if bond.configuration is not None:
        letter, refs = bond.configuration
        child = SubElement(node, tag)
        child.set('atomRefs4', ' '.join(record.ids[i] for i in refs))
        child.text = letter
        return
    wedge = bond.wedge
    if wedge == WEDGE_NONE:
        return
    child = SubElement(node, tag)
    if wedge == WEDGE_UP:
        child.text = 'W'
    elif wedge == WEDGE_DOWN:
        child.text = 'H'
    else:
        child.set('convention', 'MDL')
        child.set('conventionValue', str(WEDGE_TO_V2000[wedge]))


CML = Dialect(
    name='cml',
    # `valence=None` on purpose: CML states a hydrogen count -- as a *total*, resolved in `_finish` -- and
    # has no total-valence field at all, so `calc_implicit` drops that clause of its advice rather than
    # naming a field the schema does not have.
    channels=StatedChannels(count='a `hydrogenCount` attribute', valence=None),
    namespaces=CML_NAMESPACES,
    ns=CML_NS,
    # `reaction` is a root this dialect reads: a reaction document declaring the CML namespace reaches
    # here anyway, since `sniff` matches a namespace before a root name, so without it the namespace-less
    # spelling of the same document would be refused instead.
    tags=Tags(roots=frozenset(('cml', 'molecule', 'list', 'moleculeList', 'reaction'))),
    atom_id='id',
    # `atomRefs2` is the element form's spelling and `atomRefs` appears in CML 1; both name a pair.
    bond_refs=('atomRefs2', 'atomRefs'),
    # CML also spells the two endpoints as one attribute each.  Declared for the whole dialect, not just
    # the array hook, so the element form cannot drop a bond the array form reads.
    bond_ref_pair=('atomRef1', 'atomRef2'),
    molecule_fields=_MOLECULE_FIELDS,
    atom_fields=_ATOM_FIELDS,
    bond_fields=_BOND_FIELDS,
    molecule_ignored=_MOLECULE_IGNORED,
    atom_ignored=_ATOM_IGNORED,
    bond_ignored=_BOND_IGNORED,
    molecule_children_ignored=_MOLECULE_CHILDREN_IGNORED,
    atom_child=_atom_child,
    bond_child=_bond_child,
    molecule_child=_molecule_child,
    atom_array_hook=_atom_array,
    bond_array_hook=_bond_array,
    finish=_finish,
    document=_document,
    emit_molecule=_emit_molecule,
    emit_atom=_emit_atom,
    emit_bond=_emit_bond,
)


#: The same dialect with no namespace on its tags, which is what the writer uses.  The namespace comes off
#: because CML declares itself with a *default* namespace -- unqualified attributes included -- and
#: ``ElementTree.tostring(default_namespace=...)`` refuses exactly that, so :func:`write_cml_element` writes
#: the plain ``xmlns`` attribute instead and leaves the global ``register_namespace`` alone.
_WRITE = CML._replace(ns='')


# reading

def parse_cml(source, *, log=None, **kwargs):
    """Every ``<molecule>`` in `source` as a list of :class:`~._dialect.Record`.

    `source` is a ``str``, ``bytes``, a path or an open file.  `kwargs` reach
    :func:`~._tree.parse_xml`, where ``engine``, ``max_depth`` and ``allow_dtd`` live; nothing here can
    bypass the entity and depth policy.  Returns records rather than molecules so a caller can see the
    file's own atom ids and title before deciding to build.  :func:`read_cml` is the one-step version.
    """
    return parse_records(source, CML, log=log, **kwargs)


def read_cml(source, *, log=None, ignore_stereo=False, **kwargs):
    """Every molecule in `source`, built.  Returns a list of ``MoleculeContainer``.

    A record that cannot be built at all raises; one that is merely wrong is built and logged.  So in a
    multi-molecule document one unstorable atom refuses the whole call rather than returning a short list,
    which is indistinguishable from a file with fewer molecules in it.
    """
    return read_molecules(source, CML, log=log, ignore_stereo=ignore_stereo, **kwargs)


# writing

def record_from_molecule(mol, *, title=None, log=None):
    """`mol` as a :class:`~._dialect.Record` ready for :func:`~._dialect.write_molecule`.

    Every value lands on the field the *reader* fills -- the layout on ``Ctab.dimensionality``, the
    configuration on ``CtabAtom.parity``, the total on the ``hydrogenCount`` row's own slot -- never in a
    private scratch key, or a record straight from a file loses them on write.  Three decisions need the
    molecule and so are not a codec's: the wedges come from
    :func:`~chython.core.wedge.wedges_for_write`, the same chooser the MDL emitters use; the parities are
    written as well, ``<atomParity>`` being the only channel a molecule with no drawing has; and a dative
    bond is written with no ``order`` attribute and reported, CML having no coordination order.
    """
    out = [] if log is None else log
    record = Record(Ctab())
    ctab = record.ctab
    # Through the CTfile resolver because it hands back the molecule's S-groups as well as its title,
    # which is how this writer knows there are any to report.
    ctab.title, store = resolve_output(mol, title, None, log=out)
    ctab.meta.update(mol.meta)
    # XML cannot carry a byte that is not text, so the loss is taken by the format that cannot carry it.
    ctab.title = xml_text(ctab.title, 'title', out)
    if store is not None and (store.records or store.aliases):
        out.append(LogRecord('cml:sgroup-not-written', (),
                             f'{UNSUPPORTED}sgroup: {len(store.records) + len(store.aliases)} S-group(s) are '
                             f'not written; CML\'s S-group vocabulary is not modelled', LOST))

    sids = list(mol.atom_numbers)
    # No `*_of` accessor for the marker's index, so it comes off the atom views once, here.  CML has no
    # spelling for one -- `_emit_atom` reports the loss -- and the value still travels, so a record built
    # rather than written keeps it.
    r_indices = {a.n: a.r_index for a in mol.atoms() if a.is_r and a.r_index}
    position = {}
    has_xy = mol.has_coordinates
    # Set once and read by `emit_coordinates`.  The arena holds x and y only, so a molecule is never the
    # 3D case; a record read from a file carrying `x3`/`y3`/`z3` is, hence the ask rather than an assumption.
    ctab.dimensionality = '2D' if has_xy else ''
    # Aggregated: a molecule built in code and never given hydrogen counts has an unknown one on every
    # atom, so one line with a count and the first atom, not the whole molecule twice over.
    unknown, first_unknown = 0, None
    for sid in sids:
        atom = CtabAtom()
        atom.element = mol.element_of(sid)
        atom.charge = mol.charge_of(sid)
        atom.isotope = mol.isotope_of(sid)
        atom.radical = mol.radical_of(sid)
        atom.map_number = mol.map_number_of(sid)
        atom.r_index = r_indices.get(sid, 0)
        if has_xy:
            atom.x, atom.y = mol.xy_of(sid)
        # and otherwise the `CtabAtom` default, which nothing writes: `ctab.dimensionality` above is the
        # record's one statement about whether there is a layout.
        position[sid] = record.add_atom(atom)
        atom.file_index = position[sid] + 1
        # `is None`, never `== H_UNKNOWN`: the sentinel is what the arena stores, but what the accessor
        # *answers* for an atom holding it is `None`.  `_hydrogens.valence_for_write` tests the same way.
        implicit = mol.implicit_h_of(sid)
        if implicit is None:
            # No `hydrogenCount` at all, so no `hydrogen_total` parked: an absent count is CML for
            # "derive it", and any number here would be invented.
            unknown += 1
            if first_unknown is None:
                first_unknown = sid
        else:
            # The TOTAL under the row's own slot, which is what the attribute means, while `stated_h`
            # gets the implicit count -- the shape a reader would have left, so `Ctab.build` agrees.
            record.atom_extras[position[sid]]['hydrogen_total'] = mol.total_h_of(sid)
            atom.stated_h = implicit

    if unknown:
        out.append(LogRecord('cml:hydrogen-count-unknown', (first_unknown,),
                             f'atom: {unknown} atom(s) with hydrogen count unknown, no hydrogenCount written '
                             f'(first atom {first_unknown})'))
    # One line for the whole record, not one per atom.  The numbers are still carried on the intermediate
    # -- `Ctab.build` reads them -- and it is only the CML *document* that has nowhere to put them.
    mapped = sum(1 for a in ctab.atoms if a.map_number)
    if mapped:
        out.append(LogRecord('cml:map-numbers-not-written', (),
                             f'{UNSUPPORTED}record: {mapped} atom map number(s) are not written; CML has no '
                             f'atom-map attribute', LOST))

    # The double-bond descriptors first, because `wedges_for_write` has to be told which anchors this
    # writer is about to state some other way.  Framed: `<bondStereo>` takes an `atomRefs4`, so every
    # configured cis/trans bond is writable here and none of them is a loss.
    descriptors = cis_trans_for_write(mol, framed=True)
    wedges, _ = wedges_for_write(mol, out, cis_trans_stated={a for a, _, _ in descriptors})
    wedge_of = {(narrow, wide): code for narrow, wide, code in wedges}
    ctab_bond = {}
    for bond in mol.bonds():
        a, b = bond.n, bond.m
        code = wedge_of.get((a, b))
        if code is None and wedge_of.get((b, a)) is not None:
            a, b = b, a  # the wedge's narrow end is the bond's first atom, so it is written that way
            code = wedge_of[(a, b)]
        if bond.order == 8:
            out.append(LogRecord('cml:dative-bond-no-order', (a, b),
                                 f'{UNSUPPORTED}bond {position[a] + 1}-{position[b] + 1}: CML has no '
                                 f'coordination bond order; the bond is written with no order', LOST))
        made = CtabBond(position[a], position[b], bond.order, code or WEDGE_NONE)
        ctab_bond[frozenset((a, b))] = made
        record.add_bond(made)

    # Onto `CtabBond.configuration`, the field a file's own `<bondStereo>` lands on, so `_emit_bond` has
    # one thing to read and a record that came from a file writes its descriptors back.  The letter is
    # derived from the stored parity, never replayed from what a reader put here -- see `CtabBond`.
    for anchor, letter, frame in descriptors:
        made = ctab_bond.get(frozenset(frame[1:3]))
        if made is not None:
            made.configuration = (letter, tuple(position[i] for i in frame))

    high = len(sids) + 1   # the phantom direction's rank, as `_resolve_parities` and `core.wedge` give it
    for unit in mol.stereo_units():
        anchor = unit['anchor']
        parity = mol.parity_of(anchor)
        if unit['kind'] != SU_TETRA or not parity:
            continue
        # Onto `CtabAtom.parity`, the same field a file's own `<atomParity>` lands on, so `_emit_atom`
        # has one thing to read and a record that came from a file writes its descriptors too.
        keys = [high if r is None else position[r] for r in unit['refs']]
        if keys.count(high) > 1:
            # Two directions with no atom of their own: a quadruple cannot name them apart, so no
            # descriptor is written -- the same case `stated_parity` declines to read on the way in.
            out.append(LogRecord('cml:atom-parity-unwritable', (anchor,),
                                 f'stereo: atom {anchor}: configuration not written as an atomParity; two of '
                                 f'its four directions have no atom to name', LOST))
            continue
        # Sign `+1` for core parity 1, per `_parity_field`'s calibration, in the core's own `refs` order;
        # `_emit_atom` states the field in the frame it is measured in, so no permutation inverts a sign.
        ctab.atoms[position[anchor]].parity = _parity_field(1 if parity == 1 else -1, keys)
    return record


def write_cml_element(molecules, *, log=None, title=None):
    """`molecules` as a ``<cml>`` ``Element``.  Accepts one molecule or an iterable of them.

    A :class:`~._dialect.Record` may be passed in place of a molecule, which is how a caller writes
    back what :func:`parse_cml` read -- properties and file atom ids included -- rather than only what
    survives on the container.
    """
    out = [] if log is None else log
    if isinstance(molecules, (MoleculeContainer, Record)):
        molecules = [molecules]
    root = Element('cml')
    # The namespace as the attribute it is on the wire; see `_WRITE` for why not `default_namespace`.
    root.set('xmlns', CML_NS)
    for i, mol in enumerate(molecules, 1):
        record = mol if isinstance(mol, Record) else record_from_molecule(
            mol, title=title, log=out)
        element = write_molecule(record, _WRITE, out, parent=root)
        element.set('id', f'm{i}')
    return root


def write_cml(molecules, *, log=None, title=None, indent='  '):
    """`molecules` as a CML document.  Returns ``str``.

    `indent` is the pretty-printing step; ``None`` writes one line.
    """
    root = write_cml_element(molecules, log=log, title=title)
    return serialize(root, indent)
