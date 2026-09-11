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
"""
chython 2 as a correctness oracle for these converters, run in a subprocess.

A property chython 2's converter keeps and this one loses is a defect and fails here; a loss both have
is a gap listed in `KNOWN_GAPS`.  Compared as multisets of atoms and bonds plus a configuration count,
never as a canonical SMILES string -- that oscillates on symmetric stereocentres.
"""
from pytest import mark, skip

from .conftest import requires_indigo, requires_openbabel, requires_rdkit
from ...core.test.oracle import ask


#: Spawning chython 2 is `chython.core.test.oracle`'s job: it owns the interpreter path, the `-I` that
#: keeps this repository off the child's `sys.path`, the 2.24 pin and the outside-the-checkout check.

#: Public compounds, one group per property a converter can silently drop.
CORPUS = [
    # constitution and aromaticity
    'CCO', 'CC(=O)O', 'CC#N', 'c1ccccc1', 'C1=CC=CC=C1', 'c1ccncc1', 'c1cc[nH]c1', 'CN1C=CN=C1',
    # charge
    '[Na+].[Cl-]', 'C[N+](C)(C)C.[Cl-]', '[O-][N+](=O)c1ccccc1', 'CC(=O)[O-]', 'NCC(=O)[O-]',
    # radical, stated with the CXSMILES `|^1:|` extension: all three parsers read it identically, while
    # `[CH3]` alone they disagree on, which would be a finding about the parsers and not this package.
    '[CH3] |^1:0|', 'CC[CH2] |^1:2|', 'c1ccccc1[O] |^1:6|',
    # isotope
    '[13CH4]', '[2H]O[2H]', '[13CH3]C(=O)O', '[15NH3]',
    # tetrahedral configuration
    'N[C@@H](C)C(=O)O', 'N[C@H](C)C(=O)O', 'F[C@](Cl)(Br)I', 'C[C@H](O)[C@@H](N)CC',
    'O[C@H]1CC[C@@H](N)CC1',
    # double-bond configuration
    'C/C=C/C', 'C/C=C\\C', 'F/C=C/F', 'F/C=C\\F', 'O=C(O)/C=C\\C(=O)O',
    # heteroatom oxidation states
    'O=S(=O)(O)O', 'FC(F)(F)S(=O)(=O)O', 'CS(C)=O',
]

#: Losses both generations have, as `(smiles, property)`; `test_known_gaps_are_still_gaps` fails when
#: one closes, because a stale entry misdescribes what the converters do.
KNOWN_GAPS = ()


# The child.  A string and not a file so nothing has to be declared as package data: it is source
# handed to another interpreter, not a resource.
CHILD = r'''
from chython import smiles, MoleculeContainer

# `_payload` (the corpus, already decoded) and `_emit` come from `oracle.PREAMBLE`, prepended for us.
records = {}


def fingerprint(mol):
    """Structural content of a chython 2 molecule: atoms, bonds, and how many configurations it holds."""
    atoms = sorted((a.atomic_symbol, a.charge, a.isotope or 0, int(a.is_radical),
                    -1 if a.implicit_hydrogens is None else a.implicit_hydrogens)
                   for _, a in mol.atoms())
    syms = {n: a.atomic_symbol for n, a in mol.atoms()}
    bonds = sorted((b.order,) + tuple(sorted((syms[n], syms[m]))) for n, m, b in mol.bonds())
    tetra = sum(a.stereo is not None for _, a in mol.atoms())
    cis_trans = sum(b.stereo is not None for *_, b in mol.bonds())
    return {'atoms': atoms, 'bonds': bonds, 'tetra': tetra, 'cis_trans': cis_trans}


def inchi(rd):
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog('rdApp.*')
    rd = Chem.Mol(rd)
    rd.RemoveAllConformers()
    return Chem.MolToInchi(rd)


for text in _payload:
    rec = {}
    try:
        # Parsed and nothing else: canonicalize() would thielize a Kekule input where the core does
        # not, so the two sides would differ on a normalization this comparison is not about.
        mol = smiles(text)
    except Exception as e:
        records[text] = {'parse_error': f'{type(e).__name__}: {e}'}
        continue
    rec['input'] = fingerprint(mol)

    try:
        rd = mol.to_rdkit()
        rec['export_inchi'] = inchi(rd)
        rec['roundtrip'] = fingerprint(MoleculeContainer.from_rdkit(rd))
    except Exception as e:
        rec['rdkit_error'] = f'{type(e).__name__}: {e}'

    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog('rdApp.*')
        rd0 = Chem.MolFromSmiles(text)
        rec['import'] = None if rd0 is None else fingerprint(MoleculeContainer.from_rdkit(rd0))
    except Exception as e:
        rec['import_error'] = f'{type(e).__name__}: {e}'

    # Indigo, export only: chython 2 has no from_indigo at all, so there is no round trip to compare.
    try:
        rec['export_indigo'] = mol.to_indigo().canonicalSmiles()
    except Exception as e:
        rec['indigo_error'] = f'{type(e).__name__}: {e}'

    # OpenBabel, export only, for the same reason.  The toolkit's own canonical SMILES is handed back
    # rather than an InChI: the parent computes both InChIs with one RDKit, so the comparison cannot
    # turn into a difference between two InChI builds.
    try:
        from openbabel.openbabel import OBConversion

        conv = OBConversion()
        conv.SetOutFormat('can')
        conv.AddOption('n')  # no molecule title in the output
        rec['export_openbabel'] = conv.WriteString(mol.to_openbabel()).strip()
    except Exception as e:
        rec['openbabel_error'] = f'{type(e).__name__}: {e}'

    # CDK, export only, same reason again.  The direction chython 3 reimplemented rather than
    # re-hosted, so two independent implementations of one job are checked against each other.
    try:
        from jpype import JClass

        container = mol.to_cdk()
        flavor = JClass('org.openscience.cdk.smiles.SmiFlavor').Absolute
        rec['export_cdk'] = str(JClass('org.openscience.cdk.smiles.SmilesGenerator')(flavor)
                                .create(container))
    except Exception as e:
        rec['cdk_error'] = f'{type(e).__name__}: {e}'

    records[text] = rec

_emit(records)
'''


def _oracle():
    """The oracle's answers for the whole corpus, or a skip.  One subprocess for the module.

    An absent oracle skips; a present-but-wrong one raises, since skipping on a version mismatch would
    report a green run in which none of these comparisons happened.
    """
    return ask(CHILD, CORPUS)


_CACHE = {}


def oracle():
    if 'records' not in _CACHE:
        _CACHE['records'] = _oracle()
    return _CACHE['records']


def v3_fingerprint(mol):
    """Structural content of a core molecule, in the same shape the child emits for chython 2."""
    from chython.core import SU_TETRA
    from chython.core._core import element_symbols

    syms = element_symbols()
    # Lists and not tuples: the oracle's half arrives through JSON, which has no tuple, and a shape
    # difference would read as a fidelity difference.
    atoms = sorted([syms[a.element], a.charge, a.isotope or 0, int(a.is_radical),
                    -1 if a.implicit_h is None else a.implicit_h]
                   for a in mol.atoms())
    sym_of = {a.n: syms[a.element] for a in mol.atoms()}
    bonds = sorted([b.order, *sorted((sym_of[b.n], sym_of[b.m]))] for b in mol.bonds())
    units = [u for u in mol.stereo_units() if u['parity'] != 0]
    return {'atoms': atoms, 'bonds': bonds,
            'tetra': sum(u['kind'] == SU_TETRA for u in units),
            'cis_trans': sum(u['kind'] != SU_TETRA for u in units)}


def v3_read(text):
    from chython.core import read_smiles

    return read_smiles(text)


def constitution_layers(inchi):
    """
    The InChI layers that describe what a molecule is made of, without the configuration layers.

    Formula (which carries the hydrogen count, so a lost radical shows up in it) plus `/c`, `/h`, `/q`,
    `/p` and `/i`.
    """
    parts = inchi.split('/')
    return (parts[1], {p[0]: p for p in parts[2:] if p[0] in 'chqpi'})


def same_input_or_skip(text, new_in, old):
    """
    Both generations start from the same molecule, or this record is skipped naming the difference.

    The two SMILES parsers are not the same program, so comparing a record they read differently would
    file a parser difference as a converter defect.  One is known, and is why `CORPUS` states its
    radicals as `|^1:|`: chython 2 infers a radical from `[CH3]` and the core takes the bracket
    literally.
    """
    for key in ('atoms', 'bonds'):
        if new_in[key] != old['input'][key]:
            skip(f'{text}: the two SMILES parsers disagree about {key}, so a converter comparison would '
                 f'not be about the converters.\n  chython 2: {old["input"][key]}\n'
                 f'  chython 3: {new_in[key]}')


@requires_rdkit
@mark.parametrize('text', CORPUS)
def test_rdkit_round_trip_keeps_everything_v2_kept(text):
    """
    `X -> RDKit -> X` loses no atom property in V3 that it did not lose in chython 2.

    Compared as multisets: the two generations number atoms differently and neither is a fidelity claim.
    """
    from .._rdkit import from_rdkit, to_rdkit

    old = oracle()[text]
    if 'roundtrip' not in old:
        skip(f'chython 2 could not round-trip this record: {old}')

    mol = v3_read(text)
    new_in = v3_fingerprint(mol)
    same_input_or_skip(text, new_in, old)
    new_out = v3_fingerprint(from_rdkit(to_rdkit(mol)))

    for key in ('atoms', 'bonds'):
        old_kept = old['roundtrip'][key] == old['input'][key]
        new_kept = new_out[key] == new_in[key]
        if old_kept and not new_kept:
            raise AssertionError(
                f'{text}: chython 2 kept its {key} through an RDKit round trip and chython 3 does not.\n'
                f'  before: {new_in[key]}\n  after:  {new_out[key]}'
            )
        if not old_kept and not new_kept:
            # A loss both generations have: it has to be written down as a known gap before this passes.
            assert (text, key) in KNOWN_GAPS, (
                f'{text}: both generations lose their {key} through an RDKit round trip and this is not '
                f'in KNOWN_GAPS.\n  before: {new_in[key]}\n  after:  {new_out[key]}'
            )


@requires_rdkit
@mark.parametrize('text', CORPUS)
def test_rdkit_round_trip_keeps_every_configuration_v2_kept(text):
    """
    A round trip through RDKit loses no configuration in V3 that chython 2 kept.

    Each generation is compared against itself -- its own count before its own round trip -- because
    V2's per-atom signs and the core's units are different models.
    """
    from .._rdkit import from_rdkit, to_rdkit

    old = oracle()[text]
    if 'roundtrip' not in old:
        skip(f'chython 2 could not round-trip this record: {old}')

    mol = v3_read(text)
    new_in = v3_fingerprint(mol)
    same_input_or_skip(text, new_in, old)
    new_out = v3_fingerprint(from_rdkit(to_rdkit(mol)))

    for key in ('tetra', 'cis_trans'):
        old_lost = old['input'][key] - old['roundtrip'][key]
        new_lost = new_in[key] - new_out[key]
        assert new_lost <= max(old_lost, 0), (
            f'{text}: chython 3 lost {new_lost} {key} configuration(s) through an RDKit round trip '
            f'where chython 2 lost {old_lost} (held {new_in[key]} before, {new_out[key]} after)'
        )


@requires_rdkit
@mark.parametrize('text', CORPUS)
def test_rdkit_export_says_the_same_thing_to_inchi(text):
    """
    RDKit, asked what it received, gives the same InChI from either generation's export.

    Computed by RDKit from RDKit's own molecule, so neither chython writer appears on both sides.
    """
    from rdkit import Chem, RDLogger

    from .._rdkit import to_rdkit

    old = oracle()[text]
    if 'export_inchi' not in old:
        skip(f'chython 2 could not export this record: {old}')

    RDLogger.DisableLog('rdApp.*')
    mol = v3_read(text)
    same_input_or_skip(text, v3_fingerprint(mol), old)
    rd = Chem.Mol(to_rdkit(mol))
    rd.RemoveAllConformers()
    assert Chem.MolToInchi(rd) == old['export_inchi'], (
        f'{text}: RDKit reads a different molecule from the two exports'
    )


@requires_rdkit
@mark.parametrize('text', CORPUS)
def test_rdkit_import_keeps_everything_v2_kept(text):
    """
    Reading an RDKit molecule keeps every property chython 2's reader kept.

    The RDKit molecule is built by RDKit from the same SMILES, so both readers get the same input.
    """
    from rdkit import Chem, RDLogger

    from .._rdkit import from_rdkit

    old = oracle()[text]
    if old.get('import') is None:
        skip(f'chython 2 could not import this record: {old}')

    RDLogger.DisableLog('rdApp.*')
    rd = Chem.MolFromSmiles(text)
    new = v3_fingerprint(from_rdkit(rd))

    assert new['atoms'] == old['import']['atoms'], (
        f'{text}: the two readers disagree about the atoms.\n'
        f'  chython 2: {old["import"]["atoms"]}\n  chython 3: {new["atoms"]}'
    )
    assert new['bonds'] == old['import']['bonds'], (
        f'{text}: the two readers disagree about the bonds.\n'
        f'  chython 2: {old["import"]["bonds"]}\n  chython 3: {new["bonds"]}'
    )


@requires_rdkit
@requires_indigo
@mark.parametrize('text', CORPUS)
def test_indigo_export_agrees_on_constitution_and_charge(text):
    """
    Indigo, asked what it received, describes the same skeleton from either generation's export.

    Export only, because chython 2 has no `from_indigo`; constitution only, because its Indigo exporter
    writes no bond configuration at all, so the `/b`, `/t`, `/m` and `/s` layers would compare
    capabilities rather than readings.
    """
    from rdkit import Chem, RDLogger

    from .._indigo import to_indigo

    old = oracle()[text]
    if 'export_indigo' not in old:
        # For the radicals in the corpus: chython 2's exporter calls `setRadical(2)` and Indigo answers
        # `Unknown radical type`.  V3 writes 102, Indigo's doublet, and
        # `test_indigo.py::test_radical_round_trip` is where that is asserted.
        skip(f'chython 2 cannot export this record to Indigo, so it is no oracle here: '
             f'{old.get("indigo_error", old)}')

    mol = v3_read(text)
    same_input_or_skip(text, v3_fingerprint(mol), old)

    RDLogger.DisableLog('rdApp.*')

    def skeleton(indigo_smiles):
        rd = Chem.MolFromSmiles(indigo_smiles)
        if rd is None:
            skip(f'RDKit will not read Indigo\'s own output {indigo_smiles!r}')
        return constitution_layers(Chem.MolToInchi(rd))

    assert skeleton(to_indigo(mol).canonicalSmiles()) == skeleton(old['export_indigo']), (
        f'{text}: Indigo describes a different skeleton from the two exports\n'
        f'  chython 2: {old["export_indigo"]}\n'
        f'  chython 3: {to_indigo(mol).canonicalSmiles()}'
    )


@requires_rdkit
@requires_openbabel
@mark.parametrize('text', CORPUS)
def test_openbabel_export_agrees_on_everything_including_configuration(text):
    """
    OpenBabel, asked what it received, describes the same molecule from either generation's export.

    The strongest differential here: chython 2's OpenBabel exporter does write bond and tetrahedral
    configuration, so the whole InChI is the key rather than the constitution layers alone.  Export
    only, because chython 2 has no `from_openbabel`.
    """
    from openbabel.openbabel import OBConversion
    from rdkit import Chem, RDLogger

    from .._openbabel import to_openbabel

    old = oracle()[text]
    if 'export_openbabel' not in old:
        skip(f'chython 2 cannot export this record to OpenBabel, so it is no oracle here: '
             f'{old.get("openbabel_error", old)}')

    mol = v3_read(text)
    same_input_or_skip(text, v3_fingerprint(mol), old)

    RDLogger.DisableLog('rdApp.*')
    conv = OBConversion()
    conv.SetOutFormat('can')
    conv.AddOption('n')
    new_smiles = conv.WriteString(to_openbabel(mol)).strip()

    def described(ob_smiles):
        rd = Chem.MolFromSmiles(ob_smiles)
        if rd is None:
            skip(f'RDKit will not read OpenBabel\'s own output {ob_smiles!r}')
        return Chem.MolToInchi(rd)

    assert described(new_smiles) == described(old['export_openbabel']), (
        f'{text}: OpenBabel describes a different molecule from the two exports\n'
        f'  chython 2: {old["export_openbabel"]}\n'
        f'  chython 3: {new_smiles}'
        f'\nRead the direction before calling it a regression: chython 3 keeping something chython 2 '
        f'dropped lands here too, and is an improvement to record rather than a defect to fix.'
    )


@requires_rdkit
@mark.parametrize('text', CORPUS)
def test_cdk_export_agrees_on_everything_including_configuration(cdk, text):
    """
    CDK, asked what it received, describes the same molecule from either generation's export.

    Two independent implementations of one job, with the whole InChI as the key because chython 2's CDK
    exporter does write configuration.  The child inherits the jar through `CDK_PATH`, which `-I` does
    not suppress.
    """
    from jpype import JClass
    from rdkit import Chem, RDLogger

    from .._cdk import to_cdk

    old = oracle()[text]
    if 'export_cdk' not in old:
        # Aromatic input reaches this skip from both generations equally: CDK has no aromatic bond
        # order, so the bond goes out `Order.UNSET` with the aromatic flag and CDK's SMILES writer will
        # not write a Kekule string from an unset order.  A caller who wants CDK SMILES kekulizes first.
        skip(f'chython 2 cannot export this record to CDK, so it is no oracle here: '
             f'{old.get("cdk_error", old)}')

    mol = v3_read(text)
    same_input_or_skip(text, v3_fingerprint(mol), old)

    RDLogger.DisableLog('rdApp.*')
    flavor = JClass('org.openscience.cdk.smiles.SmiFlavor').Absolute
    generator = JClass('org.openscience.cdk.smiles.SmilesGenerator')(flavor)
    new_smiles = str(generator.create(to_cdk(mol)))

    def described(cdk_smiles):
        rd = Chem.MolFromSmiles(cdk_smiles)
        if rd is None:
            skip(f'RDKit will not read CDK\'s own output {cdk_smiles!r}')
        return Chem.MolToInchi(rd)

    assert described(new_smiles) == described(old['export_cdk']), (
        f'{text}: CDK describes a different molecule from the two exports\n'
        f'  chython 2: {old["export_cdk"]}\n'
        f'  chython 3: {new_smiles}'
        f'\nRead the direction before calling it a regression: chython 2 drops allene stereo silently '
        f'here, so chython 3 keeping it would land in this assertion as well.'
    )


def test_known_gaps_are_still_gaps():
    """Every entry in `KNOWN_GAPS` still describes a real loss; a closed gap is a stale entry."""
    assert KNOWN_GAPS == (), (
        'KNOWN_GAPS is documented as verified by this test; add the verification before adding entries'
    )
