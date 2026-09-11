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
"""The MRV name census as a ratchet: for every name ChemAxon documents, silence is forbidden.

For each name the reader must **model** it (the record differs from one built without it), **report** it
(a log line names it, or the construct it is aggregated under), or **skip it under a written claim** (an
``*_IGNORED`` / ``_ARRAY_SILENT`` frozenset in :mod:`.._mrv`).  Anything else is a silent drop."""

from re import escape, search
from typing import NamedTuple

from pytest import mark

from .._mrv import (_ARRAY_REPORTED, _ARRAY_SILENT, _ATOM_IGNORED, _BOND_IGNORED, _MOLECULE_IGNORED,
                    _ROOT_IGNORED, parse_mrv)


class Name(NamedTuple):
    """One documented name, and enough to build a document carrying it and one that does not.

    `kind` selects the builder.  `value` must *differ* from `control` or a model change is invisible
    and the row misclassifies as a silent drop; `control` of ``None`` means the control omits the name.
    `reports_as` is for a name the log names only through an enclosing construct -- a ``<scalar>`` is
    reported as part of its ``<propertyList>``.
    """
    name: str
    kind: str
    reason: str
    value: str = '1'
    control: str = None
    reports_as: str = ''
    spelling: str = ''


#: ChemAxon's documented name inventory, from
#: ``https://docs.chemaxon.com/latest/formats_marvin-documents-mrv.html``, spellings verbatim.
CENSUS = (
    # atom
    Name('elementType', 'atom', 'the element; the one atom attribute with no default', 'N', 'C'),
    Name('id', 'atom', 'the file\'s own atom identifier, which bonds refer to', 'a7', 'a3'),
    Name('formalCharge', 'atom', 'formal charge'),
    Name('hydrogenCount', 'atom', 'MRV\'s IMPLICIT count, unlike CML\'s total', '2'),
    Name('isotope', 'atom', 'mass number', '13'),
    Name('radical', 'atom', 'radical state; the page lists no values', 'monovalent'),
    Name('mrvMap', 'atom', 'atom-atom map number', '3'),
    Name('mrvValence', 'atom', 'a stated valence, which pins the implicit count', '2'),
    Name('x2', 'atom', '2D coordinate', '1.25', '2.0'),
    Name('y2', 'atom', '2D coordinate', '1.25', '0.0'),
    Name('x3', 'atom', '3D coordinate', '1.25'),
    Name('y3', 'atom', '3D coordinate', '1.25'),
    Name('z3', 'atom', '3D coordinate', '1.25'),
    # The page pairs these two and never says whether either is an attribute or a child element, so
    # they are probed as attributes on `<atom>`.  The property holds either way.
    Name('atomParity', 'atom', 'MDL parity; measured as reported, not modelled, and the page does '
                               'not say whether it is an attribute or a child element'),
    Name('atomRefs4', 'atom', 'the four neighbours a parity is stated against; measured as reported',
         'a1 a2 a3 a3'),
    # atom, silent by written claim
    Name('isSelected', 'atom', 'a GUI selection flag: drawing state, no molecule content'),
    # atom, reported onto a home that exists
    Name('mrvAlias', 'atom', 'a display label; modelled', 'OMe'),
    Name('mrvStereoGroup', 'atom', 'MDL enhanced stereo, as abs/and<N>/or<N>; modelled, onto '
                                   'Ctab.groups, in both the element and the column form', 'and1'),
    # atom, reported
    Name('mrvQueryProps', 'atom', 'query features; MRV does not read into a QueryContainer', 'X4'),
    Name('rgroupRef', 'atom', 'Markush/R-group reference; out of scope for the container model'),
    Name('attachmentPoint', 'atom', 'Markush attachment point'),
    Name('attachmentOrder', 'atom', 'Markush attachment ordering'),
    Name('ligandOrder', 'atom', 'Markush ligand ordering'),
    Name('residueId', 'atom', 'PDB residue annotation; the container stores no residue identity'),
    Name('residueType', 'atom', 'PDB residue annotation', 'ALA'),
    Name('residueAtomName', 'atom', 'PDB residue annotation', 'CA'),
    Name('mrvPseudo', 'atom', 'a pseudo-atom label; refused rather than defaulted to carbon', 'R'),
    Name('mrvExtraLabel', 'atom', 'an extra drawn label', 'x'),
    Name('mrvSetSeq', 'atom', 'atom set sequence number'),
    Name('mrvSetExtraLabelSeq', 'atom', 'extra-label set sequence number'),
    Name('chargeAngle', 'atom', 'where the charge is drawn: drawing state'),
    Name('lonePair', 'atom', 'a DRAWING count, not a chemical one; chython derives lone pairs from '
                             'the valence rules, so honouring it would let a drawing contradict them'),
    Name('mrvSpecIsotopeSymbolPreferred', 'atom', 'D/T rather than 2H/3H when drawn', 'true'),
    Name('mrvLinkNodeRep', 'atom', 'link node repetition range', '1-3'),
    Name('mrvLinkNodeOut', 'atom', 'link node outer bonds', '1 2'),
    Name('atomBicycloStereo', 'atom', 'bicyclic stereo annotation'),
    Name('sgroupRef', 'atom', 'S-group membership from the atom side; reported even though S-groups '
                              'are now read, because membership is stated by the nested molecule\'s '
                              'atomRefs and reading it twice gives the file two ways to say one '
                              'thing and this reader a way to disagree with itself', 'sg1'),
    Name('sgroupAttachmentPoint', 'atom', 'S-group attachment point; belongs to the contracted '
                                          'abbreviation this reader declines, whose atoms are not in '
                                          'the molecule that holds it'),
    Name('oneLetterName', 'atom', 'biopolymer residue letter; no container home', 'A'),
    Name('threeLetterName', 'atom', 'biopolymer residue code; no container home', 'ALA'),
    Name('reactionStereo', 'atom', 'reaction stereo annotation; the page states no owning element, '
                                   'so this row probes it on <atom>'),
    Name('correspondence', 'atom', 'a cross-reference; the page states no owning element, so this '
                                   'row probes it on <atom>', 'a1'),

    # bond
    Name('order', 'bond', 'bond order; the page says "single, double, triple, aromatic etc."',
         '2', '1'),
    Name('atomRefs2', 'bond', 'the two atoms; reversed rather than removed, a bond without it being '
                              'malformed rather than a bond with a default', 'a2 a1', 'a1 a2'),
    Name('convention', 'bond', 'cxn:coord is the dative bond and cxn:hydrogen the hydrogen bond',
         'cxn:coord'),
    Name('bondStereo', 'bond_child', 'the wedge; Marvin writes <bondStereo>C</bondStereo> with no '
                                     'atomRefs4 and the reader handles that', 'W'),
    # bond, silent by written claim
    Name('isSelected', 'bond', 'a GUI selection flag: drawing state, no molecule content'),
    # bond, reported
    Name('queryType', 'bond', 'a query bond; the page lists no values', 'SD'),
    Name('topology', 'bond', 'ring/chain topology, a query constraint'),
    Name('mrvBold', 'bond', 'drawn bold: drawing state', 'true'),
    Name('mrvHashed', 'bond', 'drawn hashed: drawing state', 'true'),
    Name('mrvReactingCenter', 'bond', 'MDL reacting-centre status'),

    # molecule
    Name('molID', 'molecule', 'the document\'s molecule identifier; dropped by written claim, and a '
                              'round trip renumbers', 'm9', 'm1'),
    Name('title', 'molecule', 'the molecule name', 'benzene'),
    Name('absStereo', 'molecule', 'the one enumeration the page gives (true/false), and still '
                                  'reported: two real exports carry a stereo group and no absStereo '
                                  'at all, so its absence is the ordinary case, and neither source '
                                  'says what its presence would mean', 'true'),
    Name('propertyList', 'mol_child', 'a molecule-level property list; modelled onto `mol.meta`, by the '
                                      'reader CML uses -- Marvin writes an SDF\'s data fields here'),
    # No `reports_as`: an entry with neither dictRef nor title is named in a line of its own, so this
    # row is the one property-list name the log spells out under its own word.
    Name('property', 'prop_child', 'one entry of a propertyList; modelled, and named in the log when '
                                   'it states no dictRef or title'),
    Name('scalar', 'prop_child', 'a scalar property value; modelled as the value of its property',
         reports_as='propertyList'),
    Name('array', 'prop_child', 'an array property value', reports_as='propertyList'),
    Name('delimiter', 'prop_child', 'the separator of an array property value',
         reports_as='propertyList'),

    # the nested <molecule> that is an S-group.  The eight roles this reader places come from two real
    # vendor exports and RDKit's Marvin reader, the page listing no values at all; a role neither
    # places is reported by name.
    Name('role', 'nested', 'the S-group type; placed onto the CTfile Sgroup type it spells, for the '
                           'eight roles two vendor exports and RDKit between them name', 'SruSgroup'),
    Name('atomRefs', 'nested', 'the S-group\'s atoms, as Ctab positions', 'a1 a2'),
    Name('bondList', 'nested', 'the S-group\'s bonds, resolved by bond id and stored as endpoint '
                               'pairs -- never as a bond position', 'b1'),
    Name('fieldName', 'nested', 'a data S-group\'s field name; SGroup.name, i.e. FIELDNAME', 'F'),
    Name('fieldData', 'nested', 'a data S-group\'s value; SGroup.data, i.e. FIELDDATA', 'v'),
    # nested, reported
    Name('fieldType', 'nested', 'a data S-group\'s field type; the page lists no values and the '
                                'reference implementation reads it nowhere, so equating it with '
                                'CTfile FIELDTYPE would be an inference and not a source', 'F'),
    Name('queryType', 'nested', 'a data S-group\'s query type; QUERYTYPE on a DataSgroup, so this '
                                'row is modelled -- reported only on a role that has no query',
         'Q'),
    Name('queryOp', 'nested', 'a data S-group\'s query operator; QUERYOP, as queryType', '='),

    # the array (column) form of an atomArray
    Name('id', 'array', 'the array element\'s own identifier; silent by written claim'),
    Name('title', 'array', 'the array element\'s own title; silent by written claim', 'x'),
    Name('convention', 'array', 'the array element\'s convention; reported rather than dropped, '
                                'because it can change what a column means', 'cml:x'),

    # the container elements' own attributes
    Name('version', 'root', 'the writer that produced the document: provenance, not molecule '
                            'content, and silent by written claim',
         'ChemAxon file format v18.11.0'),
    Name('schemaLocation', 'root', 'the schema the document validates against; silent by written '
                                   'claim, and `local` sees it bare', 'x.xsd', spelling='xsi:'),
    Name('multipageSet', 'container', 'ONE REPRESENTATIVE of MDocument\'s fifteen '
                                      'multipage*/*Set* display attributes, which the page gives '
                                      'only as globs and never enumerates'),

    # document furniture
    Name('MBracket', 'furniture', 'a drawn bracket'),
    Name('MEFlow', 'furniture', 'a drawn electron-flow arrow'),
    Name('MEllipse', 'furniture', 'a drawn ellipse'),
    Name('MElectron', 'furniture', 'a drawn electron'),
    Name('MElectronContainer', 'furniture', 'a drawn electron container'),
    Name('MPolyline', 'furniture', 'a drawn polyline'),
    Name('MRectangle', 'furniture', 'a drawn rectangle'),
    Name('MRoundedRectangle', 'furniture', 'a drawn rounded rectangle'),
    Name('MTextBox', 'furniture', 'a drawn text box'),
    Name('MNameTextBox', 'furniture', 'a drawn name text box'),
    Name('NoStructure', 'furniture', 'an explicit "no structure" placeholder'),
    Name('MMoleculeMovie', 'furniture', 'an animation'),
    Name('MPoint', 'furniture', 'a point; the page files it under the shapes that own it rather '
                                'than under MDocument, so this row probes it directly'),
    Name('MHead', 'furniture', 'document header'),
    Name('MarvinGUI', 'furniture', 'GUI state'),
    Name('mprop', 'furniture', 'a document property'),
    Name('Rgroup', 'furniture', 'a Markush R-group definition'),
    Name('RgroupBridge', 'furniture', 'a Markush R-group bridge'),
    Name('AttachmentPointArray', 'furniture', 'a Markush attachment point array'),

    # reaction
    Name('reaction', 'struct_child', 'the reaction wrapper; its molecules come back as a flat list '
                                     'with their roles reported'),
    Name('reactantList', 'struct_child', 'the reactant role', reports_as='reactant'),
    Name('agentList', 'struct_child', 'the agent role', reports_as='agent'),
    Name('productList', 'struct_child', 'the product role', reports_as='product'),
    # An arrow is document furniture.  Its `type` is glossed "e.g. EQUILIBRIUM" with no values listed
    # and the log names the element rather than the attribute, so `type` has no row of its own.
    Name('arrow', 'struct_child', 'the drawn reaction arrow; measured as reported furniture rather '
                                  'than modelled'),
)


_NS = 'http://www.chemaxon.com'
_XSI = 'http://www.w3.org/2001/XMLSchema-instance'

#: The third atom, which every ``atom``-kind row is probed on.  Unbonded because two rows vary ``id``
#: and ``elementType``, which on a bonded atom would test the response to a dangling reference.
_ATOM3 = (('id', 'a3'), ('elementType', 'C'), ('x2', '2.0'), ('y2', '0.0'))

#: The sanctioned silent skips: the frozenset, the kinds it governs, its name in :mod:`.._mrv`.  The
#: complete set of claims, kept so by the reverse ratchet below.  ``_ROOT_IGNORED`` governs two kinds
#: because a container element and the root are the same case -- read *through*, so nowhere to land.
_CLAIMS = ((_ATOM_IGNORED, ('atom',), '_ATOM_IGNORED'),
           (_BOND_IGNORED, ('bond',), '_BOND_IGNORED'),
           (_MOLECULE_IGNORED, ('molecule',), '_MOLECULE_IGNORED'),
           (_ARRAY_SILENT, ('array',), '_ARRAY_SILENT'),
           (_ROOT_IGNORED, ('root', 'container'), '_ROOT_IGNORED'))

_SANCTIONED = {kind: names for names, kinds, _ in _CLAIMS for kind in kinds}

MODELLED = 'modelled'
REPORTED = 'reported'
SILENT_BY_CLAIM = 'silent by written claim'
SILENT_DROP = 'SILENTLY DROPPED'


def _attrs(pairs):
    return ''.join(f' {k}="{v}"' for k, v in pairs)


def _build(*, root=(), mdoc=(), molecule=(), atomarray=(), atom3=None, bond=(), bond_children='',
           mol_children='', furniture='', struct_children=None, column=False):
    """A minimal MRV document with one knob per probe kind.

    Two atoms and one bond -- the least that lets a bond attribute be probed -- plus the optional third
    atom the ``atom`` rows use.  `struct_children` replaces the whole ``<molecule>``, for reaction rows.
    """
    ra = (('xmlns', _NS), ('xmlns:xsi', _XSI)) + tuple(root)
    if struct_children is None:
        if column:
            aa = ((('atomID', 'a1 a2'), ('elementType', 'C C'), ('x2', '0.0 1.0'),
                   ('y2', '0.0 0.0')) + tuple(atomarray))
            atoms = ''
        else:
            aa = tuple(atomarray)
            atoms = ('<atom id="a1" elementType="C" x2="0.0" y2="0.0"/>'
                     '<atom id="a2" elementType="C" x2="1.0" y2="0.0"/>')
            if atom3 is not None:
                atoms += f'<atom{_attrs(_override(_ATOM3, atom3))}/>'
        ba = (('id', 'b1'), ('atomRefs2', 'a1 a2'), ('order', '1'))
        body = (f'<molecule{_attrs(_override((("molID", "m1"),), molecule))}>'
                f'<atomArray{_attrs(aa)}>{atoms}</atomArray>'
                f'<bondArray><bond{_attrs(_override(ba, bond))}>{bond_children}</bond></bondArray>'
                f'{mol_children}</molecule>')
    else:
        body = struct_children
    return (f'<cml{_attrs(ra)}><MDocument{_attrs(mdoc)}>{furniture}'
            f'<MChemicalStruct>{body}</MChemicalStruct></MDocument></cml>')


def _override(base, extra):
    """`base` with `extra`'s keys replaced in place and its new keys appended.

    Order is kept so a probe and its control differ in one attribute value and not in attribute order.
    """
    extra = dict(extra)
    out = [(k, extra.pop(k)) for k, _ in base if k in extra]
    out += [(k, v) for k, v in base if k not in dict(out)]
    return tuple(sorted(out, key=lambda kv: [k for k, _ in base].index(kv[0]))) + tuple(extra.items())


_MINI_MOL = ('<molecule molID="m{n}"><atomArray>'
             '<atom id="a1" elementType="C" x2="0.0" y2="0.0"/></atomArray></molecule>')


def _pair(row):
    """The document carrying `row`'s name, and the otherwise identical one that does not."""
    spelling = row.spelling + row.name if row.spelling else row.name
    if row.kind in ('atom', 'bond', 'molecule', 'array', 'root', 'container'):
        knob = {'atom': 'atom3', 'bond': 'bond', 'molecule': 'molecule', 'array': 'atomarray',
                'root': 'root', 'container': 'mdoc'}[row.kind]
        probe = {knob: ((spelling, row.value),), 'column': row.kind == 'array'}
        if row.control is None:
            control = {knob: (), 'column': row.kind == 'array'}
        else:
            control = {knob: ((spelling, row.control),), 'column': row.kind == 'array'}
        if row.kind == 'atom':
            probe['atom3'] = probe['atom3'] or ()
            control.setdefault('atom3', ())
        return _build(**probe), _build(**control)
    if row.kind == 'bond_child':
        return (_build(bond_children=f'<{row.name}>{row.value}</{row.name}>'), _build())
    if row.kind == 'mol_child':
        # A `<propertyList/>` holding nothing states nothing, so the probe carries one entry: the empty
        # element would classify a modelled name as a silent drop.
        inner = ('<property title="F"><scalar>v</scalar></property>' if row.name == 'propertyList'
                 else '')
        return (_build(mol_children=f'<{row.name}>{inner}</{row.name}>'), _build())
    if row.kind == 'prop_child':
        return (_build(mol_children=f'<propertyList><{row.name}/></propertyList>'),
                _build(mol_children='<propertyList/>'))
    if row.kind == 'nested':
        # A role in both documents for every row but `role` itself: a nested `<molecule>` without one is
        # dropped and named, so a role-less probe would classify all of these as *reported* whether
        # S-groups are read or not.  `DataSgroup` models the most of these attributes, so it
        # distinguishes hardest -- an attribute still reported under it is one no role places.
        if row.name == 'role':
            return (_build(mol_children=f'<molecule{_attrs(((spelling, row.value),))}/>'),
                    _build(mol_children='<molecule/>'))
        role = (('role', 'DataSgroup'),)
        return (_build(mol_children=f'<molecule{_attrs(role + ((spelling, row.value),))}/>'),
                _build(mol_children=f'<molecule{_attrs(role)}/>'))
    if row.kind == 'furniture':
        return (_build(furniture=f'<{row.name}/>'), _build())
    if row.kind == 'struct_child':
        two = _MINI_MOL.format(n=1) + _MINI_MOL.format(n=2)
        if row.name == 'reaction':
            return _build(struct_children=f'<reaction>{two}</reaction>'), _build(struct_children=two)
        if row.name == 'arrow':
            return (_build(struct_children=f'<reaction>{two}<arrow type="DEFAULT"/></reaction>'),
                    _build(struct_children=f'<reaction>{two}</reaction>'))
        inner = f'<{row.name}>{two}</{row.name}>'
        return (_build(struct_children=f'<reaction>{inner}</reaction>'),
                _build(struct_children=f'<reaction>{two}</reaction>'))
    raise AssertionError(f'no builder for kind {row.kind!r}')


def _fingerprint(doc):
    """Everything the reader took from `doc`, and deliberately not the log.

    Over every ``Ctab`` slot (``sgroups`` and ``groups`` included, so a row reclassifies rather than
    fails when a feature lands), the record's atom ids and the spill dicts.  Including the log would
    call every reported name "modelled" and the module would assert nothing.
    """
    out = []
    for record in parse_mrv(doc, log=[]):
        ctab = record.ctab
        out.append((ctab.title, ctab.program, ctab.comment, ctab.dimensionality, ctab.chiral,
                    tuple(_slots(a) for a in ctab.atoms), tuple(_slots(b) for b in ctab.bonds),
                    repr(ctab.sgroups), repr(sorted(ctab.meta.items())),
                    repr(sorted(ctab.groups.items())),
                    repr(sorted(ctab.aliases.items())), tuple(record.ids),
                    tuple(repr(sorted(e.items())) for e in record.atom_extras),
                    tuple(repr(sorted(e.items())) for e in record.bond_extras),
                    repr(sorted(record.extras.items()))))
    return tuple(out)


def _slots(obj):
    return tuple(repr(getattr(obj, slot)) for slot in obj.__slots__)


def _named(lines, token):
    """`token` as a whole word in one of `lines` -- a substring would false-match names as short as
    ``id``, and the assertion is that the name appears, not that a sentence is worded a given way."""
    pattern = rf'(?<![\w:]){escape(token)}(?![\w])'
    return any(search(pattern, str(line)) for line in lines)


def classify(row):
    """What the reader does with `row`'s name: one of the four constants above.

    Reporting is tested first, since a name may be reported *and* change the record (``convention`` on
    an ``<atomArray>`` is) and the rule under test is "never neither", not "exactly one".  Two log
    scopes: a name reported under its own spelling must appear in a line the control did not produce,
    while one reported through an enclosing construct is looked for in the whole log -- the control
    carries that construct too, and demanding a new line would demand one per child.
    """
    probe, control = _pair(row)
    log, base = [], []
    parse_mrv(probe, log=log)
    parse_mrv(control, log=base)
    if row.reports_as:
        if _named(log, row.reports_as):
            return REPORTED
    elif _named([line for line in log if line not in base], row.name):
        return REPORTED
    if _fingerprint(probe) != _fingerprint(control):
        return MODELLED
    if row.name in _SANCTIONED.get(row.kind, ()):
        return SILENT_BY_CLAIM
    return SILENT_DROP


@mark.parametrize('row', CENSUS, ids=[f'{r.kind}:{r.name}' for r in CENSUS])
def test_a_documented_name_is_never_dropped_in_silence(row):
    """Model it, report it, or claim it -- silence is the only failure.

    No expected outcome is stated per name, so a reader that gets better does not fail its own census.
    """
    outcome = classify(row)
    assert outcome != SILENT_DROP, (
        f'MRV {row.kind} name {row.name!r} is read, produces no log line and is in no written claim '
        f'set -- a silent drop, which is the one outcome the input posture forbids.  Reason it is in '
        f'the census: {row.reason}.  Either give it a Field row, let the walker report it, or add it '
        f'to the frozenset in _mrv.py that claims honouring it would build the same molecule.')


@mark.parametrize('claims,kinds,label', _CLAIMS, ids=[label for _, _, label in _CLAIMS])
def test_every_sanctioned_silent_skip_is_in_the_inventory(claims, kinds, label):
    """A claim set cannot grow without the census witnessing the name.

    Otherwise silencing a name is a one-line change to a frozenset and the ratchet stops noticing it.
    """
    inventory = {row.name for row in CENSUS if row.kind in kinds}
    assert not claims - inventory, (
        f'_mrv.py\'s {label} claims these names may be skipped in silence, and the census above does '
        f'not name them: {sorted(claims - inventory)}.  Add a row with the reason a reader honouring '
        f'the name would build the same molecule.')


def test_the_reported_array_attributes_are_disjoint_from_the_silent_ones():
    """One name, one outcome, for the one kind that has both a silent set and a reported set."""
    assert not _ARRAY_SILENT & _ARRAY_REPORTED


def test_every_outcome_is_represented_so_the_ratchet_is_armed():
    """A classifier that answered one thing always would pass the census and mean nothing.

    Counts rather than per-name expectations, so a name moving between outcomes does not touch this.
    Measured 2026-09-05 over 98 names: 28 modelled, 63 reported, 7 silent by claim, 0 silent drops.
    The floors sit well under those on purpose -- they catch a stuck classifier, not today's reader.
    """
    outcomes = [classify(row) for row in CENSUS]
    assert outcomes.count(MODELLED) >= 15
    assert outcomes.count(REPORTED) >= 40
    assert outcomes.count(SILENT_BY_CLAIM) >= 5
    assert SILENT_DROP not in outcomes


def test_the_fingerprint_is_stable_across_two_reads_of_one_document():
    """Otherwise every row would look modelled and the census would assert nothing."""
    doc = _build(atom3=(('mrvAlias', 'OMe'),))
    assert _fingerprint(doc) == _fingerprint(doc)


def test_the_census_names_no_name_twice_for_one_kind():
    """`queryType` and `id` and `title` each appear twice, on different elements, and only that."""
    seen = [(row.kind, row.name) for row in CENSUS]
    assert len(seen) == len(set(seen))
