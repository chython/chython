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
"""ChemAxon's MRV, the second dialect on the engine in :mod:`._dialect`: a table and a few handlers.

Discriminated by NAMESPACE (``http://www.chemaxon.com``) and never by root name, MRV's root also being
called ``cml``.  Vocabulary from ChemAxon's documentation and ``mrvSchema_18_11_0.xsd``, value sets
cross-checked against RDKit's BSD Marvin parser.  Records land in ``ctfile._ctab.Ctab``."""

from xml.etree.ElementTree import Element, SubElement

from ._dialect import (Dialect, Field, Record, Tags, apply_fields, choose_coordinates, decode_float,
                       decode_int, emit_coordinates, emit_properties, encode_element,
                       encode_from_table, encode_nonzero, encode_stated, outermost, parse_records,
                       read_molecules, read_properties, resolve_element, serialize, synthetic_node,
                       write_molecule, xml_text)
from ._tree import local, text_of
from ..ctfile import Ctab, CtabAtom, CtabBond
from ..ctfile._ctab import WEDGE_FROM_V2000, WEDGE_TO_V2000
from ..ctfile._hydrogens import H_MAX, StatedChannels
from ..ctfile._sgroup import NO_INDEX, SGroup, UNSUPPORTED, resolve_output
from ...core.wedge import cis_trans_for_write, wedges_for_write
from ...core import (LogRecord, LOST, MoleculeContainer, REPAIRED, STEREO_ABS, STEREO_AND, STEREO_OR,
                     WEDGE_DOWN, WEDGE_NONE, WEDGE_UP)


__all__ = ['MRV', 'MRV_NS', 'parse_mrv', 'read_mrv', 'record_from_molecule', 'write_mrv',
           'write_mrv_element']


#: The namespace every Marvin document declares, and the one this writer writes.
MRV_NS = 'http://www.chemaxon.com'

#: What :func:`~._dialect.sniff` matches a document against.  One URI: unlike CML, which accumulated a
#: dozen historical schema locations, ChemAxon has published exactly this one.
MRV_NAMESPACES = frozenset((MRV_NS,))

#: MRV's four bond orders.  Case-folded on read: the format states them upper-case and ``a`` appears.
#: ``A`` is aromatic and is stored as order 4 as stated, not kekulised.
_ORDERS = {'1': 1, '2': 2, '3': 3, 'A': 4}

#: Order back out.  A separate table rather than an inversion: an inverted many-to-one dict picks
#: whichever key came last, silently.
_ORDERS_OUT = {1: '1', 2: '2', 3: '3', 4: 'A'}

#: ``convention`` on a ``<bond>``: MRV's spelling of a coordination bond, and of a hydrogen bond.
#: A dative bond carries this attribute and **no** ``order``, which is why the two cannot share a row.
_COORD, _HYDROGEN = 'cxn:coord', 'cxn:hydrogen'

#: MRV's radical names, as unpaired-electron counts.  Eight spellings for four counts -- ``divalent1``
#: and ``divalent3`` are the singlet and triplet carbene, which differ in a way a molecule's one
#: radical bit cannot hold -- so the count is what is read and the multiplicity is reported.
_RADICALS = {'monovalent': 1, 'divalent': 2, 'divalent1': 2, 'divalent3': 2, 'trivalent': 3,
             'trivalent2': 3, 'trivalent4': 3, '4': 4}

#: One radical bit back out.  ``monovalent`` is one unpaired electron, which is what the bit means.
_RADICAL_OUT = 'monovalent'

#: ``mrvStereoGroup``'s three spellings, as the prefix each uses and the stereo class it names.  The
#: group *number* is inside the token -- ``or1``, ``and2`` -- unlike CTfile, where it is a separate
#: collection index.  Tokens are case-folded, real exports writing ``or1`` lower case.  Mind the
#: direction, as easy to invert here as in CTfile's collections: an ``and`` group is a racemate (both
#: enantiomers present) and an ``or`` group one enantiomer of unknown identity.  See
#: :data:`~chython.formats.ctfile._ctab.STEREO_FROM_COLLECTION`.
_STEREO_GROUPS = {'abs': STEREO_ABS, 'and': STEREO_AND, 'or': STEREO_OR}

#: A stereo class back out, as the token's prefix.  A separate table rather than an inversion, as above.
#: ``STEREO_UNSPECIFIED`` is absent on purpose: it is the class of a centre nobody has classified, which
#: MRV spells by writing no group at all.
_STEREO_GROUPS_OUT = {STEREO_ABS: 'abs', STEREO_AND: 'and', STEREO_OR: 'or'}

#: ``elementType`` values a molecule cannot hold.  Empty, and that is the finding: MRV's non-element
#: vocabulary is ``R`` with an ``rgroupRef`` and ``*``, each of which names ONE atom, so each is the
#: marker.  A set-valued query type -- CML inherits MDL's ``A``, ``Q``, ``X`` -- has no MRV spelling.
_PSEUDO = {}

#: The column form's null, in every column that has one.  ``-`` is what a Marvin writer puts where an
#: atom states nothing and its neighbours state something; the ``radical`` column spells the same thing
#: ``0``, which its own codec handles.
_NULL = '-'


# codecs -- only those reading MRV's own vocabulary; the dialect-agnostic ones come from `_dialect`

#: A bond order as MRV text, or ``None`` to leave the attribute off.  Order 8 gets ``None`` here and
#: its ``convention="cxn:coord"`` from the `convention` row, which is how MRV itself writes a
#: coordination bond: the convention attribute and no order at all.
_order_out = encode_from_table(_ORDERS_OUT)


def _order_in(text):
    """MRV ``order`` as a chython bond order."""
    key = text.strip().upper()
    if key in _ORDERS:
        return _ORDERS[key]
    raise ValueError(f'{text.strip()!r} is not one of MRV\'s bond orders 1, 2, 3, A')


def _convention_in(text):
    """A ``<bond>``'s ``convention``, parked as-stated for :func:`_finish`.

    Parked and not applied: ``convention`` outranks ``order`` and the two attributes arrive in whatever
    sequence the file wrote them, so a row writing straight onto ``CtabBond.order`` would give a dative
    bond order 1 whenever the writer put ``order`` second.
    """
    key = text.strip().lower()
    if key in (_COORD, _HYDROGEN):
        return key
    raise ValueError(f'{text.strip()!r} is not one of MRV\'s bond conventions {_COORD}, {_HYDROGEN}')


def _convention_out(order):
    """``cxn:coord`` for a coordination bond, and nothing for any other order."""
    return _COORD if order == 8 else None


def _radical_in(text):
    """An MRV ``radical`` name as an unpaired-electron count, or ``None`` for a column-form null.

    ``0`` is the null the ``radical`` *column* uses and ``-`` the null every other column uses, so both
    answer ``None``.  The count rather than the name, eight names mapping onto four counts.
    """
    key = text.strip()
    if not key or key in ('0', _NULL):
        return None
    try:
        return _RADICALS[key.lower()]
    except KeyError:
        raise ValueError(f'{key!r} is not one of MRV\'s radical names '
                         f'{", ".join(sorted(_RADICALS))}') from None


def _radical_out(radical):
    """A radical bit as an MRV radical name.  One bit is one unpaired electron, so ``monovalent``."""
    return _RADICAL_OUT if radical else None


def _h_in(text):
    """``hydrogenCount`` as the implicit count it is.

    Range-checked here rather than in `finish`: a count outside what an atom can hold is a malformed
    attribute, which the engine logs while keeping the atom.
    """
    count = decode_int(text)
    if count < 0:
        raise ValueError(f'hydrogenCount {count} is negative')
    if count > H_MAX:
        raise ValueError(f'hydrogenCount {count} is above the {H_MAX} an atom can hold')
    return count


def _valence_in(text):
    """``mrvValence`` as the stated total valence.  Negative is refused; RDKit refuses it too."""
    value = decode_int(text)
    if value < 0:
        raise ValueError(f'mrvValence {value} is negative')
    return value


def _alias_in(text):
    """``mrvAlias`` as the display label it is, or ``None`` for the ``0`` placeholder.

    Otherwise verbatim, as the V2000 ``A  <n>`` path is: an alias is free text, and an XML attribute value
    has no fixed-width padding to strip.  An **empty** alias is stored rather than dropped, matching V2000,
    where the header is the statement that an alias exists and the text only its content.

    ``0`` is "no alias here", the placeholder :func:`_stereo_group_in` reads in ``mrvStereoGroup="0 and1
    0"``.  The column form needs one, every cell of a column being filled; MEASURED at Marvin 25.1.3, the
    element form reads it the same way, ``<atom mrvAlias="0"/>`` coming back out of ``molconvert`` with no
    ``mrvAlias``.  So the placeholder is the attribute's and not the column's, and one codec serves both.

    The rest of ChemAxon's escape convention (``"zero"`` for the string ``0``, ``"."`` for empty) is **not**
    applied: the same measurement returns both verbatim, so reading them as escapes would rewrite a label
    the drawing needs.
    """
    return None if text == '0' else text


def _alias_out(text):
    """An alias back out, or ``None`` -- an atom with no alias -- to leave the attribute off.

    ``encode_stated`` rather than a truthiness test: see :func:`_alias_in` for why an empty alias is a
    statement.  The one text with no spelling here is ``0`` itself, which :func:`_emit_atom` reports.
    """
    return encode_stated(text)


def _stereo_group_in(text):
    """``mrvStereoGroup`` as the ``(kind, group)`` pair :attr:`Ctab.groups` is keyed to, or ``None``.

    ``None`` for the format's three ways of saying "in no group": ``0`` (what a real export writes in
    every cell but the one with a group), ``-`` and an empty value.  The number lives inside the token, so
    the split is prefix-plus-digits and a bare ``abs`` decodes to the core's group 0, which is what an
    absolute centre already is.  Out-of-domain numbers are checked by :meth:`Ctab.build`, which has the
    atom to name in the message.
    """
    key = text.strip().lower()
    if not key or key in ('0', _NULL):
        return None
    prefix = key.rstrip('0123456789')
    if prefix not in _STEREO_GROUPS:
        raise ValueError(f'{text.strip()!r} is not one of MRV\'s stereo groups '
                         f'{", ".join(sorted(_STEREO_GROUPS))}<n>')
    digits = key[len(prefix):]
    return _STEREO_GROUPS[prefix], int(digits) if digits else 0


def _stereo_group_out(pair):
    """A ``(kind, group)`` pair as an MRV stereo-group token, or ``None`` to leave the attribute off.

    An unclassified centre gets no token, :data:`_STEREO_GROUPS_OUT` having no row for it.  The number is
    written only when there is one: an absolute group carries the core's group 0 and MRV spells it bare
    ``abs``.
    """
    if pair is None:
        return None
    kind, group = pair
    prefix = _STEREO_GROUPS_OUT.get(kind)
    if prefix is None:
        return None
    return f'{prefix}{group}' if group else prefix


# the table

_MOLECULE_FIELDS = (
    Field('title', 'title', lambda t: t, lambda t: t or None),
)

_ATOM_FIELDS = (
    # `element_type` is not a `CtabAtom` slot, so the engine spills it and `_finish` resolves it: the
    # resolution has to tell "stated nothing" from "stated carbon", and a row writing onto
    # `CtabAtom.element` -- whose default is carbon -- cannot.
    Field('elementType', 'element_type', lambda t: t.strip(), encode_element, write_slot='element'),
    # The marker's index, which `elementType` does not carry: `elementType="R" rgroupRef="1"` is what
    # `molconvert mrv` writes for an `M  RGP` atom.  `encode_nonzero`, an unindexed marker having no
    # group to reference -- and `_finish` must not let the resolver's own zero overwrite what this row set.
    Field('rgroupRef', 'r_index', decode_int, encode_nonzero),
    Field('formalCharge', 'charge', decode_int, encode_nonzero),
    Field('isotope', 'isotope', decode_int, encode_nonzero),
    # Parked and resolved, like CML's spin multiplicity: a count above one unpaired electron has to be
    # reported, and a codec has no log.
    Field('radical', 'radical_name', _radical_in, _radical_out, write_slot='radical'),
    # `encode_nonzero` and not `encode_stated`: `mrvMap="0"` is MRV for unmapped, so a zero is silence.
    Field('mrvMap', 'map_number', decode_int, encode_nonzero),
    # `encode_stated` and not `encode_nonzero`: `mrvValence="0"` is a *stated* zero valence.
    Field('mrvValence', 'valence', _valence_in, encode_stated),
    # Straight onto `stated_h`: MRV's `hydrogenCount` is the IMPLICIT count -- the hydrogens *not* drawn
    # -- so it is already what `Ctab.build` passes to the hydrogen derivation.  CML's is the total and
    # needs its explicit neighbours subtracted, which is the one place the two dialects differ in meaning
    # rather than in spelling.
    Field('hydrogenCount', 'stated_h', _h_in, encode_stated),
    # Not a `CtabAtom` attribute, so the engine spills it and `_finish` moves it to `Ctab.aliases`, the
    # position-keyed dict the V2000 `A  <n>` path fills.  The spill is left in place rather than popped,
    # which is what lets a `Record` out of `parse_mrv` be written straight back.  Consequence: `_cells`
    # hands a `0` cell in this column to the codec instead of skipping it, which is why `_alias_in`
    # answers `None` for one.
    Field('mrvAlias', 'alias', _alias_in, _alias_out),
    # Parked and resolved like the alias: the home is `Ctab.groups`, a dict keyed by atom position, and a
    # `Field` can only set an attribute on the one atom it is handed.  Both forms are real and one row
    # serves both: Marvin Sketch states it as an `<atomArray mrvStereoGroup="or1 0 0 0 0 0 0">` column,
    # MarvinJS as an `<atom>` attribute.  The row's existence makes `_cells` hand a `0` cell in that column
    # to the codec instead of skipping it, which is why `_stereo_group_in` answers `None` for one.
    Field('mrvStereoGroup', 'stereo_group', _stereo_group_in, _stereo_group_out),
    # Coordinates are read by the table and written by `_emit_atom`: *which* set to write is a question
    # about the record, and a row is handed one atom's one value, so a row encoding `x2` unconditionally
    # writes a drawing nobody made onto a record with no layout.
    Field('x2', 'x', decode_float, None),
    Field('y2', 'y', decode_float, None),
    Field('x3', 'x3', decode_float, None),
    Field('y3', 'y3', decode_float, None),
    Field('z3', 'z3', decode_float, None),
)

_BOND_FIELDS = (
    Field('order', 'order', _order_in, _order_out),
    # One attribute, two directions, and the slots differ: the reader parks the file's text for
    # `_finish` to resolve against `order`, and the writer reads the finished order and spells order 8
    # as MRV spells it.  `queryType` has no row at all -- a query bond order is not a bond order -- and
    # the engine names it in the log.
    Field('convention', 'convention', _convention_in, _convention_out, write_slot='order'),
)

#: Attributes that carry no chemistry, so no `unsupported:` line.  Each is a claim that a reader
#: honouring it would build the same molecule.  `molID` is the document's own name for the molecule and
#: is dropped, `write_mrv_element` inventing `m<i>` from its own enumeration; `id` on an *atom* is
#: preserved, a `<bond atomRefs2=...>` pointing at it.
_MOLECULE_IGNORED = frozenset(('molID',))
#: The Marvin GUI's selection state.  A drawing-program flag: the same atoms either way.
#:
#: Two sets holding one name, not one set read twice: each entry is a claim about the construct it sits
#: on, and the atom and bond vocabularies are independent, so sharing the object would make an atom-only
#: attribute silently silent on bonds as well.
_ATOM_IGNORED = frozenset(('isSelected',))
_BOND_IGNORED = frozenset(('isSelected',))

#: The array form's identity column, under the element form's name for it.  Every other column is
#: spelled exactly as the per-atom attribute is, which is what lets one table serve both forms.
_ARRAY_ALIASES = {'atomID': 'id'}

#: Attributes an ``<atomArray>`` or ``<bondArray>`` may carry that are **not columns**:
#:
#: * ``id`` -- the array element's own XML identity.  The array form spells its identity *column*
#:   ``atomID``, so reading a bare ``id`` as a one-entry column invents an atom out of an empty array;
#: * ``title`` -- a human label on the array; a reader honouring it builds the same atoms;
#: * ``convention`` -- **reported**, not silent.  It names the dictionary the array's content is defined
#:   in, and a dictionary this reader has not read can redefine what every column means.  Same reasoning
#:   `_bond_child` applies to a ``convention`` on a ``<bondStereo>``.
_ARRAY_SILENT = frozenset(('id', 'title'))
_ARRAY_REPORTED = frozenset(('convention',))


# the column form -- an <atomArray> holding one attribute per column

def _columns(node, log, where, member):
    """``{name: [values]}`` for the array form, or ``None`` when `node` is in the element form.

    The form is decided by the children, not the attributes: an ``<atomArray>`` holding ``<atom>`` children
    is the element form whatever attributes it carries, since the vocabulary a global XML attribute may come
    from is open and a skip list would eventually lose every atom of an element-form array.  A node carrying
    both contradicts itself; the children win and the columns are named.  Called for **both** forms, which
    is why non-column attributes are handled here -- in the element form this returns ``None``.
    """
    out = {}
    for name, value in node.attrib.items():
        plain = local(name)
        key = _ARRAY_ALIASES.get(plain, plain)
        if plain in _ARRAY_SILENT:
            continue
        if plain in _ARRAY_REPORTED:
            log.append(LogRecord('mrv:array-foreign-dict', (),
                                 f'{UNSUPPORTED}{where}: <{local(node.tag)}> {plain}={value!r:.40} names a '
                                 f'dictionary this reader has not read; the columns are read under MRV\'s own '
                                 f'meanings', LOST))
            continue
        out[key] = value.split()
    members = sum(1 for child in node if local(child.tag) == member)
    if members:
        if out:
            log.append(LogRecord('mrv:array-column-conflict', (),
                                 f'{where}: {len(out)} array column(s) ({", ".join(sorted(out))}) beside '
                                 f'{members} <{member}> element(s); the elements are read and the columns '
                                 f'dropped', LOST))
        return None
    if not out:
        return None
    width = min(len(v) for v in out.values())
    ragged = {k: len(v) for k, v in out.items() if len(v) != width}
    if ragged:
        # Truncated rather than dropped: one short column must not lose every atom in the file.
        log.append(LogRecord(
            'mrv:array-ragged-columns', (),
            f'{where}: array columns have different lengths ({ragged}); truncated to {width}', REPAIRED))
    return {k: v[:width] for k, v in out.items()}, width


def _cells(columns, i, table):
    """One row of the array as ``{attribute: text}``, with the nulls left out.

    Two nulls, and the second is why `table` -- the dialect's index for this member -- is an argument.
    ``-`` is the null every column has.  A cell of ``0`` in a column with **no row** carries nothing:
    ``lonePair``, ``sgroupRef`` and ``rgroupRef`` are written for every atom with ``0`` meaning "not this
    atom's business" (``sgroupRef="0 0 0 0 sg1"``).  A column that *does* have a row keeps its zeroes,
    ``hydrogenCount="0"`` being a stated absence, so a column gaining a row must handle its own null in
    the codec -- as :func:`_radical_in` and :func:`_stereo_group_in` both do.
    """
    out = {}
    for key, values in columns.items():
        value = values[i]
        if value == _NULL or (value == '0' and key not in table):
            continue
        out[key] = value
    return out


def _atom_array(node, record, log):
    """The column form of ``<atomArray>``.  ``True`` when it claimed the node.

    Not optional for this dialect: real Marvin files write the column form, and the engine walks
    ``<atom>`` children only, so without this hook such a file comes back with no atoms and no bonds and a
    log holding nothing but dropped-bond lines.
    """
    found = _columns(node, log, 'atom', MRV.tags.atom)
    if found is None:
        return False
    columns, width = found
    ids = columns.pop('id', None)
    table = MRV.index('atom')  # once for the array, not once per row: `index` rebuilds the dict
    for i in range(width):
        atom = CtabAtom()
        position = record.add_atom(atom, ids[i] if ids else None, log)
        atom.file_index = position + 1
        where = f'atom {record.ids[position]}'
        apply_fields(synthetic_node('atom', _cells(columns, i, table)), table, atom, log, where,
                     _ATOM_IGNORED, ('id',), spill=record.atom_extras[position])
    return True


def _bond_array(node, record, log):
    """The column form of ``<bondArray>``, which is named rather than read.  ``True`` when claimed.

    No source describes the column vocabulary for bonds, so such an array is reported once and its bonds
    lost loudly rather than guessed from the atom array's shape.  It also collects the bond ids, which is
    why it runs for the element form it declines: an S-group's ``bondList`` names bonds by ``id`` while the
    engine registers only *atom* ids, so a ``{bond id: (atom id, atom id)}`` table is built here and
    resolved in :func:`_sgroups`.  Only ids the file **states** are collected -- a positional fallback could
    resolve a ``bondList`` to the wrong bond rather than to nothing.
    """
    names = record.extras.setdefault('bond_names', {})
    for child in node:
        if local(child.tag) != MRV.tags.bond:
            continue
        ident = child.get(MRV.atom_id)
        refs = (child.get(MRV.bond_refs[0]) or '').split()
        if ident is not None and len(refs) == 2:
            names[ident] = tuple(refs)
    found = _columns(node, log, 'bond', MRV.tags.bond)
    if found is None:
        return False
    columns, width = found
    log.append(LogRecord('mrv:bond-array-column-form', (),
                         f'{UNSUPPORTED}bond: <bondArray> states {width} bond(s) as {len(columns)} column(s) '
                         f'({", ".join(sorted(columns))}); the column form of a bond array is not modelled', LOST))
    return True


# <bondStereo>, in the three spellings MRV writes it in

def _bond_child(node, record, position, log):
    """``<bondStereo>`` as ``CtabBond.wedge``.  ``True`` when it claimed the node.

    Three spellings, each a closed set in both sources: the bare letters ``W``/``H`` as the element's
    text, ``dictRef="cml:W"``/``"cml:H"``, and ``convention="MDL" conventionValue="1|3|4|6"``.  ``C`` and
    ``T`` are the fourth thing the text may be and are a *double-bond* configuration, landing on
    ``CtabBond.configuration``.
    """
    if local(node.tag) != 'bondStereo':
        return False
    bond = record.ctab.bonds[position]
    ident = f'b{position + 1}'
    convention = node.get('convention', '')
    if convention and convention.upper() != 'MDL':
        # A `convention` names the dictionary its content is defined in, and MDL's is the only one this
        # reader has.  Checked before the MDL branch, so a `conventionValue` under a foreign dictionary is
        # not decoded as a CTfile code either -- that number is that dictionary's too.
        log.append(LogRecord('mrv:bond-stereo-foreign-dict', (),
                             f'{UNSUPPORTED}bond {ident}: bondStereo convention {convention!r:.30} names a '
                             f'dictionary this reader has not read; nothing applied', LOST))
        return True
    if convention or node.get('conventionValue') is not None:
        # The CTfile bond-stereo number written straight through, dressed as a dictionary reference.
        # Decoded with the CTfile reader's own table, so 1, 4 and 6 have one meaning in this tree.
        raw = node.get('conventionValue', '')
        try:
            code = int(raw)
        except ValueError:
            log.append(LogRecord(
                'mrv:bond-stereo-bad-value', (),
                f'bond {ident}: bondStereo conventionValue {raw!r:.20} is not a number, dropped', LOST))
            return True
        if code in WEDGE_FROM_V2000:
            bond.wedge = WEDGE_FROM_V2000[code]
        elif code == 3:
            # MDL 3 on a double bond is "cis or trans, unknown which": a real construct with no field in a
            # molecule, an unset configuration and an explicitly unknown one being the same arena value.
            log.append(LogRecord('mrv:bond-stereo-cis-trans-unknown', (),
                                 f'{UNSUPPORTED}bond {ident}: "cis or trans, unknown which" is not modelled', LOST))
        else:
            # `WEDGE_FROM_V2000` has 0 as well as 1, 4 and 6 -- 0 is "no wedge", a stated absence -- so
            # the numbers this message names are the four it reads plus the 3 the branch above answers.
            log.append(LogRecord('mrv:bond-stereo-bad-code', (),
                                 f'bond {ident}: bondStereo conventionValue {code} is not 0, 1, 3, 4 or 6, '
                                 f'dropped', LOST))
        return True

    dictref = node.get('dictRef', '')
    text = (dictref.split(':')[-1] if dictref else text_of(node)).upper()
    if text == 'W':
        bond.wedge = WEDGE_UP
    elif text == 'H':
        bond.wedge = WEDGE_DOWN
    elif text in ('C', 'T'):
        # Stored with an EMPTY frame, which is the one substantive difference from CML: no source
        # describes an `atomRefs4` on an MRV `<bondStereo>`, so the letter names nothing it is measured
        # over.  `stated_cis_trans` reads a bare letter only where a terminal cannot carry a second
        # substituent, which is exactly where "which pair is C" has one answer.
        bond.configuration = (text, ())
        return True
    elif not text:
        log.append(LogRecord('mrv:bond-stereo-empty', (), f'bond {ident}: empty bondStereo, dropped', LOST))
    else:
        # A bare line and not `unsupported: `: MRV's `<bondStereo>` text is a closed set -- `W`, `H`, `C`,
        # `T` plus the `convention`/`dictRef` spellings above -- so a fifth letter is a value outside the
        # dialect's vocabulary rather than a construct this reader declines to model.
        log.append(LogRecord('mrv:bond-stereo-bad-letter', (),
                             f'bond {ident}: bondStereo {text!r:.20} is not one of W, H, C or T, dropped', LOST))
    return True


# the document level -- everything between the root and a <molecule>

#: Document elements that are pure containers: they hold molecules and state nothing else about them,
#: so the flat list the walker returns loses nothing.  `MDocument` is the wrapper every MRV file has and
#: `MChemicalStruct` is the structure half of it.
_CONTAINERS = frozenset(('MDocument', 'MChemicalStruct', 'molecule'))

#: Role elements of an ``<reaction>``, which :func:`_roles` reports rather than these naming themselves.
_ROLE_CONTAINERS = frozenset(('reaction', 'reactantList', 'agentList', 'productList',
                              'reactant', 'agent', 'product'))

#: Attributes the *container* elements may carry that state nothing about the molecules inside them.
#: Both are what the real exports in ``test/`` carry: ``version="ChemAxon file format v18.11.0, generated
#: by v19.7.0"`` names the writer that produced the document -- provenance this library declines to claim
#: on the way out, see :func:`write_mrv_element` -- and ``schemaLocation`` the schema it validates
#: against.  ``local`` sees ``schemaLocation`` bare, its ``xsi:`` prefix being a namespace.
_ROOT_IGNORED = frozenset(('version', 'schemaLocation'))


def _document(root, log):
    """Report what the document states about its molecules that a flat list of them cannot hold.

    Three things.  A ``<reaction>``'s content *is* the roles, and the walker reads straight through one,
    so the counts are what lets a caller holding only the log recover the record's shape.  Marvin's
    document furniture -- text boxes, arrows, reaction signs, electron containers, polylines -- is a real
    construct with nothing here to hold it.  And the container elements' own attributes, ``<MDocument>``
    and ``<MChemicalStruct>`` needing no handler to be read *through*.  All three are aggregated per
    document rather than named per element.
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
        log.append(LogRecord('mrv:reaction-roles-not-modelled', (),
                             f'{UNSUPPORTED}record: <reaction> {ident} roles are not modelled; {what}', LOST))

    furniture = {}
    _furniture(root, furniture)
    if furniture:
        named = ', '.join(f'{n} <{tag}>' for tag, n in sorted(furniture.items()))
        log.append(LogRecord('mrv:document-furniture', (),
                             f'{UNSUPPORTED}record: {named} outside any <molecule>; MRV\'s document furniture is '
                             f'not modelled', LOST))

    carried = {}
    for node in _containers(root):
        tag = local(node.tag)
        for name in node.attrib:
            plain = local(name)
            if plain not in _ROOT_IGNORED:
                carried.setdefault(tag, set()).add(plain)
    for tag, names in sorted(carried.items()):
        log.append(LogRecord('mrv:display-settings-not-modelled', (),
                             f'{UNSUPPORTED}record: <{tag}> carries {", ".join(sorted(names))}; MRV\'s document '
                             f'display settings are not modelled', LOST))


def _containers(root):
    """`root` and every pure-container element under it, stopping at a ``<molecule>``.

    The three elements reading *through* loses nothing on -- the root, ``<MDocument>``,
    ``<MChemicalStruct>`` -- and no further.  A ``<reaction>`` and its role elements are excluded, the
    construct as a whole being reported already; a ``<molecule>``'s attributes reach
    :func:`~._dialect.apply_fields`, which enforces the same rule, so a fragment rooted at one yields
    nothing.
    """
    found = []
    stack = [root]
    while stack:
        node = stack.pop()
        if local(node.tag) == MRV.tags.molecule:
            continue
        found.append(node)
        stack.extend(child for child in node if local(child.tag) in _CONTAINERS)
    return found


def _roles(node, role, counts):
    """Count the ``<molecule>`` descendants of `node` by the role element that encloses each.

    `role` is the nearest enclosing element's name with a trailing ``List`` stripped, so that
    ``<productList><product><molecule>`` and ``<productList><molecule>`` -- both of which real files
    write -- count as one role rather than two.  The descent stops at a ``<molecule>``, as the walker does.
    """
    for child in node:
        tag = local(child.tag)
        if tag == MRV.tags.molecule:
            counts[role] = counts.get(role, 0) + 1
        else:
            _roles(child, tag[:-4] if tag.endswith('List') else tag, counts)


def _furniture(node, counts):
    """Count the outermost document elements that are neither containers nor part of a reaction."""
    for child in node:
        tag = local(child.tag)
        if tag in _CONTAINERS:
            if tag != 'molecule':
                _furniture(child, counts)
        elif tag in _ROLE_CONTAINERS:
            _furniture(child, counts)
        else:
            counts[tag] = counts.get(tag, 0) + 1


# finish -- everything that needs the whole molecule

def _finish(record, log):
    """Resolve everything that could not be decided one attribute at a time.

    The element (its messages need the atom's name), the radical bit (a count above one has to be
    reported), the alias and the stereo group (both live in dicts on the ``Ctab``, which a
    :class:`Field` cannot reach), the coordinates (2D and 3D are two attribute sets), a bond's
    ``convention`` (it outranks ``order``, so it cannot be applied while attributes are still arriving)
    and the S-groups, which name atoms and bonds and so need the whole molecule.  **The S-groups run
    last, after the conventions**: a ``cxn:hydrogen`` bond is dropped there, and an S-group naming one
    must not come out holding a pair that is no longer a bond.
    """
    ctab = record.ctab
    atoms = ctab.atoms
    # Aggregated, not one line per atom: an `<atomArray>` with no `elementType` anywhere is one defect,
    # and a line per atom would bury every other finding in the record.
    carbons, first_carbon = 0, ''

    for position, (atom, spill) in enumerate(zip(atoms, record.atom_extras)):
        where = f'atom {record.ids[position]}'
        if 'element_type' in spill:
            atom.element, isotope, label, r_index = resolve_element(spill['element_type'], where, log, _PSEUDO)
            if isotope and not atom.isotope:
                atom.isotope = isotope
            if r_index:  # `or`, not a plain assignment: `rgroupRef` is the index's own spelling here
                atom.r_index = r_index
            if label is not None:
                # `setdefault`, so an `mrvAlias` -- which says "label" outright -- outranks a label
                # recovered from the element column, exactly as V2000's `A  <n>` line does.
                ctab.aliases.setdefault(position, label)
        else:
            # MRV requires an `elementType`, so an `<atom>` without one is a broken file; carbon is the
            # only thing to do with it, and it is said out loud rather than defaulted silently.
            carbons += 1
            if not first_carbon:
                first_carbon = where

        electrons = spill.get('radical_name')
        if electrons:
            atom.radical = True
            if electrons > 1:
                log.append(LogRecord('mrv:radical-not-modelled', (),
                                     f'{UNSUPPORTED}{where}: {electrons} unpaired electrons are not modelled; '
                                     f'read as one radical centre', LOST))

        # `is not None` and not `in`: the `0` placeholder decodes to `None`, and an empty alias -- which is
        # a statement -- must still land.  The dict is the one the V2000 `A  <n>` path fills, keyed by atom
        # position; `Ctab.build` translates those positions into stable ids and calls `set_aliases`.
        if spill.get('alias') is not None:
            ctab.aliases[position] = spill['alias']

        group = spill.get('stereo_group')
        if group is not None:
            # The same dict MDL's two spellings fill -- V3000's COLLECTION block and V2000's Sgroup
            # encoding -- so `Ctab.build` is the single caller of `set_stereo_group`, which refuses a
            # number outside the core's domain with the atom named.
            ctab.groups[position] = group

    if carbons:
        log.append(LogRecord('mrv:no-element-type', (),
                             f'atom: {carbons} atom(s) with no elementType, read as carbon '
                             f'(first {first_carbon})', REPAIRED))
    choose_coordinates(record, log)
    _conventions(record, log)
    _sgroups(record, log)


def _conventions(record, log):
    """Apply every parked bond ``convention``, which outranks the bond's ``order``.

    A ``cxn:hydrogen`` bond is **dropped** -- the one place this dialect removes something it read.  A
    hydrogen bond is not a covalent bond of any order, and kept as the single bond its missing ``order``
    would default to it fuses two molecules the file drew as two.
    """
    ctab = record.ctab
    keep, keep_extras, dropped = [], [], 0
    for position, (bond, spill) in enumerate(zip(ctab.bonds, record.bond_extras)):
        convention = spill.get('convention')
        if convention == _COORD:
            bond.order = 8
        elif convention == _HYDROGEN:
            dropped += 1
            continue
        keep.append(bond)
        keep_extras.append(spill)
    if dropped:
        log.append(LogRecord('mrv:hydrogen-bond-dropped', (),
                             f'{UNSUPPORTED}bond: {dropped} hydrogen bond(s) dropped; a molecule has no bond '
                             f'order for one, and reading it as single would join two molecules the file drew '
                             f'apart', LOST))
        ctab.bonds[:] = keep
        record.bond_extras[:] = keep_extras


# S-groups -- a <molecule> inside a <molecule>, which is MRV's whole spelling for one.  ChemAxon's page
# lists `role` and glosses it "S-group type like SRU" with no enumeration, so the role names below are
# the ones RDKit's Marvin reader dispatches on (BSD-3, read for the vocabulary only), each paired with
# the CTfile type RDKit's own molfile writer emits for it.

#: ``role`` -> the CTfile Sgroup type it is the MRV spelling of.  The model is
#: :class:`~chython.formats.ctfile._sgroup.SGroup`, the same store the V2000 and V3000 parsers fill, so an
#: S-group read here writes back out as an SDF Sgroup with no second translation.
#:
#: Two roles are absent: ``MulticenterSgroup`` has no CTfile type at all, naming a point that is the mean
#: of several atoms rather than a set of them, and contracted ``MolTemplateSgroup``-style groups are
#: refused by shape rather than by name -- see :func:`_sgroups`.
_SGROUP_ROLES = {'SruSgroup': 'SRU', 'CopolymerSgroup': 'COP', 'ModificationSgroup': 'MOD',
                 'MultipleSgroup': 'MUL', 'DataSgroup': 'DAT', 'GenericSgroup': 'GEN',
                 'MonomerSgroup': 'MON', 'SuperatomSgroup': 'SUP'}

#: A type back out, as the ``role`` MRV spells it.  A separate table rather than an inversion, and here it
#: is also the list of types MRV can express, so a molecule carrying an Sgroup type absent from it earns a
#: line rather than silently losing a record.
_SGROUP_ROLES_OUT = {'SRU': 'SruSgroup', 'COP': 'CopolymerSgroup', 'MOD': 'ModificationSgroup',
                     'MUL': 'MultipleSgroup', 'DAT': 'DataSgroup', 'GEN': 'GenericSgroup',
                     'MON': 'MonomerSgroup', 'SUP': 'SuperatomSgroup'}

#: Per type, the nested ``<molecule>``'s own attributes and the CTfile keyword each one is.  Keyed by
#: type rather than shared, because ``title`` means two different keywords: on a polymer bracket or a
#: superatom it is the bracket's ``LABEL``, and on a multiple group it is ``MULT``, the repeat *count*.
#:
#: A type with no entry for an attribute does not model it there, so the attribute is reported -- which is
#: how ``fieldName`` on an ``SruSgroup`` gets named, CTfile having no ``FIELDNAME`` on a repeat unit.
_SGROUP_KEYWORDS = {
    'SRU': {'title': 'LABEL', 'connect': 'CONNECT'},
    'COP': {'title': 'LABEL', 'connect': 'CONNECT'},
    'MOD': {'title': 'LABEL', 'connect': 'CONNECT'},
    'SUP': {'title': 'LABEL'},
    'MON': {'title': 'LABEL'},
    'MUL': {'title': 'MULT'},
    'GEN': {},
    # `fieldName` and `fieldData` are NOT here: they have dedicated slots on the record (`name` and
    # `data`) rather than a keyword bag entry, and `_sgroups` fills those directly.  `queryType` and
    # `queryOp` do not, so they ride in `fields` under their CTfile keywords.
    'DAT': {'queryType': 'QUERYTYPE', 'queryOp': 'QUERYOP'},
}

#: The nested element's structural attributes, consumed by :func:`_sgroups` itself rather than by a
#: keyword.  ``molID`` and ``id`` are the document's names for the group -- dropped for the reason
#: :data:`_MOLECULE_IGNORED` drops the outer ``molID``, since the writer regenerates them.
_SGROUP_STRUCTURAL = frozenset(('role', 'id', 'molID', 'atomRefs', 'bondList'))


def _molecule_child(node, record, log):
    """``<propertyList>`` as the record's data fields; a nested ``<molecule>`` -- MRV's S-group -- parked
    for :func:`_finish`.

    Parked and not resolved: ``atomRefs`` and ``bondList`` name atoms and bonds by the file's own ids and
    nothing says the S-group element follows the arrays declaring them, so a hook resolving as it arrived
    would read a forward reference as a dangling one.
    """
    tag = local(node.tag)
    if tag in ('propertyList', 'property'):
        # The same reader CML uses: `molconvert mrv` on an SDF writes the SD fields here, so this is the
        # channel a data field arrives on when a Marvin document is the source.
        return read_properties(node, record, log, rule='mrv')
    if tag != MRV.tags.molecule:
        return False
    record.extras.setdefault('sgroups', []).append(node)
    return True


def _sgroups(record, log):
    """Every parked nested ``<molecule>`` as an :class:`~chython.formats.ctfile._sgroup.SGroup`.

    The references land in the alphabet a ``Ctab`` uses -- 0-based atom positions, bonds as ``(a, b)``
    position pairs rather than bond numbers.  ``bondList`` names bonds by *id*, so it resolves through the
    id table :func:`_bond_array` builds and then through the atom ids, never through a bond's position:
    this dialect has one bond it deletes.

    The record is left **unnumbered** (``NO_INDEX``): MRV names an S-group with a string (``sg1``) and
    states no Sgroup number, and both CTfile writers fall back to the record's position for one with none.
    """
    nodes = record.extras.get('sgroups')
    if not nodes:
        return
    ctab = record.ctab
    names = record.extras.get('bond_names', {})
    # The bonds that still exist, as a membership test for a `bondList` entry.  Built after `_conventions`
    # has run, which is why this function is last: a `cxn:hydrogen` bond is gone by now.
    present = {frozenset((b.a, b.b)) for b in ctab.bonds}

    for i, node in enumerate(nodes, 1):
        where = f'sgroup {node.get("id") or i}'
        role = (node.get('role') or '').strip()
        if not role:
            log.append(LogRecord('mrv:sgroup-no-role', (),
                                 f'{where}: <molecule> in <molecule> states no role, so nothing says what kind of '
                                 f'S-group it is; dropped', LOST))
            continue
        # By shape and not by role, which makes it one test instead of a list: a contracted abbreviation
        # carries its own <atomArray>, so its atoms are not in this molecule and its `atomRefs` would
        # resolve against the wrong alphabet.
        own = sorted({local(c.tag) for c in node} & {MRV.tags.atom_array, MRV.tags.bond_array})
        if own:
            log.append(LogRecord('mrv:sgroup-contracted', (),
                                 f'{UNSUPPORTED}{where}: a {role} carrying its own <{">, <".join(own)}> states a '
                                 f'contracted group whose atoms are not in this molecule; not modelled', LOST))
            continue
        stype = _SGROUP_ROLES.get(role)
        if stype is None:
            log.append(LogRecord('mrv:sgroup-role-not-modelled', (),
                                 f'{UNSUPPORTED}{where}: role {role!r:.40} is not modelled; the roles this '
                                 f'reader places are {", ".join(sorted(_SGROUP_ROLES))}', LOST))
            continue

        sg = SGroup(stype)
        lost = 0
        for name in (node.get('atomRefs') or '').split():
            position = record.index_of.get(name)
            if position is None:
                lost += 1
            else:
                sg.atoms.append(position)
        for name in (node.get('bondList') or '').split():
            pair = names.get(name)
            positions = None if pair is None else tuple(record.index_of.get(x) for x in pair)
            if positions is None or None in positions or frozenset(positions) not in present:
                lost += 1
            else:
                sg.bonds.append(positions)
        if lost:
            log.append(LogRecord('mrv:sgroup-bad-refs', (),
                                 f'{where}: {lost} reference(s) name an atom or bond this molecule does not '
                                 f'have, dropped', LOST))

        keywords = _SGROUP_KEYWORDS[stype]
        for attribute, keyword in keywords.items():
            value = node.get(attribute)
            if value is not None:
                sg.fields.setdefault(keyword, []).append(value)
        consumed = _SGROUP_STRUCTURAL | set(keywords)
        if stype == 'DAT':
            consumed |= {'fieldName', 'fieldData', 'x', 'y'}
            sg.name = node.get('fieldName', '')
            # `x`/`y` are the drawn label's anchor, which is what CTfile states in FIELDDISP's first two
            # columns, so they land there and no styling tail is invented for them.
            try:
                sg.disp = (float(node.get('x')), float(node.get('y')), '')
            except (TypeError, ValueError):
                if node.get('x') is not None or node.get('y') is not None:
                    log.append(LogRecord('mrv:sgroup-anchor-not-a-number', (),
                                         f'{where}: x/y {node.get("x")!r:.20}/{node.get("y")!r:.20} do '
                                         f'not parse as a label anchor; no FIELDDISP taken from them',
                                         LOST))
            data = node.get('fieldData')
            if data is not None:
                # Encoded because `SGroup.data` is `list[bytes]`: an SDF data field need not be UTF-8, so
                # the store holds bytes.  An XML attribute arrived decoded, so the encoding is known here.
                sg.data.append(data.encode('utf8'))

        unmodelled = sorted({local(k) for k in node.attrib} - consumed)
        if unmodelled:
            log.append(LogRecord('mrv:sgroup-attr-not-modelled', (),
                                 f'{UNSUPPORTED}{where}: {", ".join(unmodelled)} on a {role} '
                                 f'{"is" if len(unmodelled) == 1 else "are"} not modelled', LOST))
        ctab.sgroups.append(sg)


def _emit_sgroups(record, element, count, log):
    """Write `record`'s S-groups as nested ``<molecule>`` elements of `element`.  Returns `count`.

    `count` is the document's running molecule number, in and out: MRV requires a ``molID`` on every
    ``<molecule>`` and a nested one is a ``<molecule>``, so nested groups take numbers from the same
    sequence.  Called after :func:`~._dialect.write_molecule` rather than from ``emit_molecule`` for
    element order -- Marvin writes S-groups after the arrays, the hook runs before them.  Bond ids are read
    back off the handed ``<bondArray>``, the engine deciding what a bond is called, so a ``bondList`` cannot
    point at a name nothing wrote.
    """
    ctab = record.ctab
    if not ctab.sgroups:
        return count
    bond_id = {}
    for child in element:
        if local(child.tag) == MRV.tags.bond_array:
            for position, node in enumerate(child):
                bond = ctab.bonds[position]
                bond_id[frozenset((bond.a, bond.b))] = node.get(MRV.atom_id)
            break

    for i, sg in enumerate(ctab.sgroups, 1):
        where = f'sgroup {i} {sg.type}'
        role = _SGROUP_ROLES_OUT.get(sg.type)
        if role is None:
            log.append(LogRecord('mrv:sgroup-no-mrv-role', (),
                                 f'{UNSUPPORTED}{where}: MRV has no role for a {sg.type} group, not written', LOST))
            continue
        count += 1
        child = SubElement(element,
                           f'{{{MRV_NS}}}molecule' if element.tag.startswith('{') else 'molecule')
        child.set('molID', f'm{count}')
        child.set('id', f'sg{i}')
        child.set('role', role)
        child.set('atomRefs', ' '.join(record.ids[a] for a in sg.atoms))
        if sg.bonds:
            # `SGroup.bonds` holds atom references, so a pair whose bond an edit removed is storable; this
            # is where it stops being expressible, and it is reported rather than dropped silently.
            named = [bond_id.get(frozenset(pair)) for pair in sg.bonds]
            child.set('bondList', ' '.join(x for x in named if x is not None))
            missing = sum(1 for x in named if x is None)
            if missing:
                log.append(LogRecord('mrv:sgroup-bond-not-written', (),
                                     f'{where}: {missing} bond reference(s) name a pair this molecule has no bond '
                                     f'for, not written', LOST))
        if sg.type == 'DAT':
            if sg.name:
                child.set('fieldName', sg.name)
            # The label anchor, always written: Marvin 25.1.3 returns an empty `<MDocument>` for a
            # `DataSgroup` with no `x`/`y`, and its own writer states `x="0.0000" y="0.0000"` for a group
            # whose source file carried no FIELDDISP.  A zero pair anchors a label at the frame's origin
            # and states nothing about an atom, so it is not the invented drawing `x2`/`y2` would be.
            x, y = sg.disp[:2] if sg.disp is not None else (0.0, 0.0)
            child.set('x', f'{x:.4f}')
            child.set('y', f'{y:.4f}')
            if sg.data:
                child.set('fieldData', sg.field_data)

        # The keyword bag, through the same table the reader reads, so a keyword MRV can express is written
        # under the attribute it arrived as.  Popped from a copy: the record is the caller's.
        fields = {k: list(v) for k, v in sg.fields.items()}
        for attribute, keyword in _SGROUP_KEYWORDS[sg.type].items():
            values = fields.pop(keyword, None)
            if values:
                child.set(attribute, values[0])
                if len(values) > 1:
                    log.append(LogRecord('mrv:sgroup-keyword-repeated', (),
                                         f'{where}: {keyword} stated {len(values)} times; MRV has one '
                                         f'{attribute} attribute, so the first is written', LOST))
        unwritten = sorted(fields)
        if sg.subtype:
            unwritten.append('SUBTYPE')
        if sg.patoms:
            unwritten.append('PATOMS')
        if sg.cstates:
            unwritten.append('CSTATE')
        if sg.disp is not None and (sg.type != 'DAT' or sg.disp[2].strip()):
            # A DAT group's anchor goes out as `x`/`y` above; FIELDDISP's styling columns after it do not.
            unwritten.append('FIELDDISP styling' if sg.type == 'DAT' else 'FIELDDISP')
        if sg.parent != NO_INDEX:
            unwritten.append('PARENT')
        if unwritten:
            log.append(LogRecord('mrv:sgroup-field-not-written', (),
                                 f'{UNSUPPORTED}{where}: {", ".join(unwritten)} '
                                 f'{"has" if len(unwritten) == 1 else "have"} no MRV spelling, not written', LOST))
    return count


# writing -- the mirror of the handlers above

def _emit_molecule(record, node, log):
    """``<propertyList>`` for whatever ``record.ctab.meta`` holds, before the arrays as Marvin writes it.

    Both ``dictRef`` and ``title`` carry the field name, which is what Marvin 25.1.3 writes converting an
    SDF; a reader taking either one gets the name.
    """
    emit_properties(record, node, MRV_NS if node.tag.startswith('{') else '',
                    names=('dictRef', 'title'))


def _emit_atom(record, position, node, log):
    """The coordinates, plus the one alias text this dialect has no spelling for."""
    emit_coordinates(record, record.ctab.atoms[position], node)
    if node.get('mrvAlias') == '0':
        # `0` is the placeholder in both forms, so this label reads back as no label -- here and in Marvin.
        # Written anyway, the document then carrying the text for a reader that takes it literally.
        log.append(LogRecord('mrv:alias-text-is-the-placeholder', (record.ids[position],),
                             f'{UNSUPPORTED}atom {record.ids[position]}: an alias whose text is "0" '
                             f'has no MRV spelling, "0" being mrvAlias\'s placeholder for no alias; '
                             f'written as stated', LOST))


def _emit_bond(record, position, node, log):
    """``<bondStereo>`` for a wedged bond or a double-bond configuration, in the spelling this reader's
    own text branch reads.

    The bare letters, ``dictRef`` and the MDL dictionary reference being accommodations this reader accepts
    only on the way in; an ``either`` wedge has no letter and goes out as its MDL number.  One element name,
    two constructs, never both on a bond: ``W``/``H`` is which end of a single bond is nearer the viewer,
    ``C``/``T`` how a double bond is configured, and the configuration is checked first because a bond that
    has one carries no wedge.  No ``atomRefs4``, Marvin writing none -- so only the unambiguous descriptors
    reach here; CML, whose files carry the frame, writes it.
    """
    bond = record.ctab.bonds[position]
    if bond.configuration is not None:
        child = SubElement(node,
                           f'{{{MRV_NS}}}bondStereo' if node.tag.startswith('{') else 'bondStereo')
        child.text = bond.configuration[0]
        return
    wedge = bond.wedge
    if wedge == WEDGE_NONE:
        return
    child = SubElement(node, f'{{{MRV_NS}}}bondStereo' if node.tag.startswith('{') else 'bondStereo')
    if wedge == WEDGE_UP:
        child.text = 'W'
    elif wedge == WEDGE_DOWN:
        child.text = 'H'
    else:
        child.set('convention', 'MDL')
        child.set('conventionValue', str(WEDGE_TO_V2000[wedge]))


MRV = Dialect(
    name='mrv',
    # Both channels are attributes rather than MDL constructs: an `.mrv` states an implicit count as
    # `hydrogenCount` and a total valence as `mrvValence`, so advice about an `MRV_IMPLICIT_H` data
    # S-group -- a molfile's way of carrying the same quantity -- would be wrong here.
    channels=StatedChannels(count='a `hydrogenCount` attribute', valence='an `mrvValence` attribute'),
    namespaces=MRV_NAMESPACES,
    ns=MRV_NS,
    # `cml` is the root name here, so `roots` does not discriminate this dialect -- the namespace does.  It
    # is populated only so `sniff`'s second pass has an answer for a Marvin document declaring no namespace;
    # `molecule` is here because a fragment pulled out of one is still this vocabulary.
    tags=Tags(roots=frozenset(('cml', 'molecule'))),
    atom_id='id',
    bond_refs=('atomRefs2',),
    molecule_fields=_MOLECULE_FIELDS,
    atom_fields=_ATOM_FIELDS,
    bond_fields=_BOND_FIELDS,
    molecule_ignored=_MOLECULE_IGNORED,
    atom_ignored=_ATOM_IGNORED,
    bond_ignored=_BOND_IGNORED,
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


#: The same dialect with no namespace on its tags, which is what the writer uses.  Not a second table:
#: ``_replace`` copies the one above, so the two cannot drift.  The namespace comes off because MRV
#: declares itself with a *default* namespace, under which attribute names are unqualified, and
#: ``ElementTree.tostring(default_namespace=...)`` refuses a document with unqualified attributes -- which
#: is every MRV document.  So the declaration is written as the plain ``xmlns`` attribute it is on the
#: wire, leaving the process-wide ``register_namespace`` untouched.
_WRITE = MRV._replace(ns='')


# reading

def parse_mrv(source, *, log=None, **kwargs):
    """Every ``<molecule>`` in `source` as a list of :class:`~._dialect.Record`.

    `source` is a ``str``, ``bytes``, a path or an open file.  `kwargs` reach :func:`~._tree.parse_xml`,
    where the entity and depth policy lives; no keyword here can bypass it.

    Returns records rather than molecules so a caller can see the file's own atom ids and its title before
    deciding to build.  :func:`read_mrv` is the one-step version.
    """
    return parse_records(source, MRV, log=log, **kwargs)


def read_mrv(source, *, log=None, ignore_stereo=False, **kwargs):
    """Every molecule in `source`, built.  Returns a list of ``MoleculeContainer``.

    A record that cannot be built at all raises; a record that is merely wrong is built and logged.  For a
    multi-molecule document one unstorable atom therefore refuses the whole call rather than returning a
    short list, which would be indistinguishable from a file with fewer molecules in it.
    """
    return read_molecules(source, MRV, log=log, ignore_stereo=ignore_stereo, **kwargs)


# writing

def record_from_molecule(mol, *, title=None, log=None):
    """`mol` as a :class:`~._dialect.Record` ready for :func:`~._dialect.write_molecule`.

    Every value lands in a form a :class:`~._dialect.Field` encodes.  Three decisions need the molecule:
    wedges from :func:`~chython.core.wedge.wedges_for_write`, the chooser the MDL emitters use; the
    implicit hydrogen count, which is the one MRV states, written only when known (``hydrogenCount="0"``
    is a statement); and the double-bond configurations from
    :func:`~chython.core.wedge.cis_trans_for_write` with ``framed=False``, MRV's bare
    ``<bondStereo>C``/``T`` having nowhere to name a frame, so the ambiguous ones are reported as lost.
    """
    out = [] if log is None else log
    record = Record(Ctab())
    ctab = record.ctab
    # Through the CTfile resolver because it hands back the molecule's S-groups and aliases as well as
    # its title, and decodes their `bytes` the same way.
    ctab.title, store = resolve_output(mol, title, None, log=out)
    ctab.meta.update(mol.meta)
    # XML cannot carry a byte that is not text, so the loss is taken by the format that cannot carry it.
    ctab.title = xml_text(ctab.title, 'title', out)
    aliases = store.aliases if store is not None else {}

    sids = list(mol.atom_numbers)
    # The marker's index has no `*_of` accessor, so it comes off the atom views once, here.
    r_indices = {a.n: a.r_index for a in mol.atoms() if a.is_r and a.r_index}
    position = {}
    has_xy = mol.has_coordinates
    # The record's own answer to "have I a layout", set once and read by `emit_coordinates`.  The arena
    # holds x and y only, so a molecule is never the 3D case; a record read from an MRV file carrying
    # `x3`/`y3`/`z3` is, which is why the writer asks rather than assumes.
    ctab.dimensionality = '2D' if has_xy else ''
    # Aggregated, as `_finish` aggregates its own: a molecule built in code has an unknown count on every
    # atom, so one line, a count, and the first atom it happened on.
    unknown, first_unknown = 0, None
    # The canonical form of the groups, which the V3000 emitter asks for too: a group's *number* is
    # arbitrary, so two writes of one molecule must not disagree about which arbitrary number it got.
    groups = {}
    if mol.has_stereo_groups:
        for (kind, group), members in mol.canonical_stereo_groups().items():
            for sid in members:
                groups[sid] = (kind, group)
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
        position[sid] = record.add_atom(atom)
        atom.file_index = position[sid] + 1
        if sid in aliases:
            # Both places, matching the read path: the spill is what `_encode_into` writes from, while
            # `Ctab.aliases` is what `Ctab.build` puts back if the caller builds this record rather than
            # writing it.  Same for the stereo group below.
            record.atom_extras[position[sid]]['alias'] = aliases[sid]
            ctab.aliases[position[sid]] = aliases[sid]
        if sid in groups:
            record.atom_extras[position[sid]]['stereo_group'] = groups[sid]
            ctab.groups[position[sid]] = groups[sid]
        # `is None` and not `== H_UNKNOWN`: the sentinel is what the arena stores, and what the accessor
        # *answers* for an atom holding it is `None`.
        implicit = mol.implicit_h_of(sid)
        if implicit is None:
            unknown += 1
            if first_unknown is None:
                first_unknown = sid
        else:
            atom.stated_h = implicit

    if unknown:
        out.append(LogRecord('mrv:hydrogen-count-unknown', (),
                             f'atom: {unknown} atom(s) with hydrogen count unknown, no hydrogenCount written '
                             f'(first atom {first_unknown})', LOST))

    # Translated into the Ctab's alphabet, not written here: the mirror of `Ctab.build`, so a caller who
    # builds this record rather than writing it gets the S-groups back.  What each one can be *spelled* as
    # is `_emit_sgroups`' question.
    if store is not None and store.records:
        translated = store.translate(position)
        out.extend(translated.log)
        ctab.sgroups.extend(translated.records)
    # The double-bond descriptors first: `wedges_for_write` has to be told which anchors are stated some
    # other way.  `framed=False` because MRV's `<bondStereo>C` is a bare letter -- Marvin writes no
    # `atomRefs4` on one -- so a configuration on a bond whose terminal carries a second substituent is not
    # writable here, and `wedges_for_write` keeps its line for that anchor.
    descriptors = cis_trans_for_write(mol, framed=False)
    wedges, _ = wedges_for_write(mol, out, cis_trans_stated={a for a, _, _ in descriptors})
    wedge_of = {(narrow, wide): code for narrow, wide, code in wedges}
    ctab_bond = {}
    for bond in mol.bonds():
        a, b = bond.n, bond.m
        code = wedge_of.get((a, b))
        if code is None and wedge_of.get((b, a)) is not None:
            a, b = b, a  # the wedge's narrow end is the bond's first atom, so it is written that way
            code = wedge_of[(a, b)]
        made = CtabBond(position[a], position[b], bond.order, code or WEDGE_NONE)
        ctab_bond[frozenset((a, b))] = made
        record.add_bond(made)

    # Onto `CtabBond.configuration`, where this reader's own `<bondStereo>C` lands, so `_emit_bond` has one
    # thing to read.  The letter is derived from the stored parity, never replayed from the reader's -- see
    # the `CtabBond` docstring.  The frame is dropped: MRV has nowhere to put one.
    for _, letter, frame in descriptors:
        made = ctab_bond.get(frozenset(frame[1:3]))
        if made is not None:
            made.configuration = (letter, ())
    return record


def write_mrv_element(molecules, *, log=None, title=None):
    """`molecules` as a ``<cml>`` ``Element`` in MRV's vocabulary.  One molecule or an iterable.

    ``<cml><MDocument><MChemicalStruct>`` is the wrapper every Marvin document has; the two inner elements
    are pure containers, so they need no handler on the way in and two lines here.  No ``version``
    attribute is written: real files name their writer there (``version="ChemAxon file format v18.11.0,
    ..."``), and this is not that writer.  A :class:`~._dialect.Record` may be passed in place of a
    molecule, which writes back what :func:`parse_mrv` read, the file's own atom ids included.
    """
    out = [] if log is None else log
    if isinstance(molecules, (MoleculeContainer, Record)):
        molecules = [molecules]
    root = Element('cml')
    # The namespace as the attribute it is on the wire; see `_WRITE` for why not `default_namespace`.
    root.set('xmlns', MRV_NS)
    struct = SubElement(SubElement(root, 'MDocument'), 'MChemicalStruct')
    # A running count and not an enumerate: an S-group is a nested `<molecule>` and takes a `molID` from the
    # same sequence, as a real Marvin document does.
    count = 0
    for mol in molecules:
        record = mol if isinstance(mol, Record) else record_from_molecule(mol, title=title, log=out)
        element = write_molecule(record, _WRITE, out, parent=struct)
        # `molID` and not `id`, MRV's own spelling.  The value is this loop's count, not the file's -- see
        # `_MOLECULE_IGNORED`, which drops that on the way in.
        count += 1
        element.set('molID', f'm{count}')
        count = _emit_sgroups(record, element, count, out)
    return root


def write_mrv(molecules, *, log=None, title=None, indent='  '):
    """`molecules` as an MRV document.  Returns ``str``.

    `indent` is the pretty-printing step; ``None`` writes one line.  Indented by default because these
    files are read by people at least as often as by programs, and an indented document diffs.
    """
    root = write_mrv_element(molecules, log=log, title=title)
    return serialize(root, indent)
