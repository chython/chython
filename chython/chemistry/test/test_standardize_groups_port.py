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
"""chython 2's per-rule standardization gate, ported: one drawn-wrong input, one repaired answer, per
rule.  V2's original is `git show 483ddaa:chython/algorithms/standardize/test/test_groups.py`; its `==`
compared SMILES strings, here it is the canonical form, so graphs are compared and not spellings.  The
three divergences are `REDERIVED` in DATA, `NEEDS_INFERRED_RADICAL` and `NOT_A_SMILES_ROUNDTRIP`.  The
`fix_tautomers` flag is gated here too: only a per-rule case list can say which rows it withholds.
"""
import pytest

from . import gen_standardize_rules as gen
from .. import standardize
from .._standardize import LogRecord
from .._tables import groups_rules, metals_rules
from ...core import read_smiles


def _repairs(log):
    """The rules `standardize()` fired.

    `mol.log` is one storage and the reader writes there too -- an input drawn wrong enough to reach a
    rule is often drawn wrong enough for the reader to have said so first -- so a record's stage is what
    says which call produced it.
    """
    return [record for record in log if record.stage != 'read']


# Verbatim from V2 apart from the two REDERIVED expectations, V2's own comments included, so that a
# diff against `git show 483ddaa:...` shows only those two.
DATA = [
    ('CP(C)(C)C', 'C[P+](C)(C)C'),
    ('CB1(C)[H]B(C)(C)[H]1', 'CB1(C)~[H]B(C)(C)~[H]1'),
    ('[O]N(C)[NH]', '[O-][N+](C)=N'), ('[O]N(C)[CH2]', '[O-][N+](C)=C'),
    ('[O]S(C)(C)[O]', 'O=S(C)(C)=O'), ('[O]S(C)(C)[S]', 'O=S(C)(C)=S'),
    ('BN(C)=C', 'B~N(C)=C'),
    ('B=N(C)(C)C', 'B~N(C)(C)C'),
    ('BS(C)C', 'B~S(C)C'), ('BO(C)C', 'B~O(C)C'),
    ('[B-]=[N+](C)C', 'BN(C)C'), ('C[B-]=[N+]C', 'CBNC'), ('[B-]=[N+]', 'BN'),
    ('[O-][B+3]([O-])([O-])[O-]', 'O[B-](O)(O)O'),
    ('[O-]B(O)(O)O', 'O[B-](O)(O)O'),
    ('OB(O)(O)O', 'O[B-](O)(O)O'),
    ('CN(C)(C)C', 'C[N+](C)(C)C'),
    ('C=N(=O)O', 'C[N+](=O)[O-]'),
    ('C=N(=O)C', 'C=[N+]([O-])C'), ('O=N(=O)C', 'O=[N+]([O-])C'), ('N=N(=O)C', 'N=[N+]([O-])C'),
    ('C=[N+]([O-])O', 'C[N+](=O)[O-]'),
    ('CN(=O)=N(=O)C', 'C[N+]([O-])=[N+]([O-])C'),
    ('CN(=O)=N(=N)C', 'C[N+]([O-])=[N+]([NH-])C'), ('CN(=O)=N(=NC)C', 'C[N+]([O-])=[N+]([N-]C)C'),
    # REDERIVED.  V2 wrote `C=[N+]([N-])C` and `N=[N+]([N-])C`; the anion keeps the one hydrogen it came
    # in with, and V2's reader could not tell `[N-]` (0 H in V3) from `[NH-]` (1 H).
    ('C=N(=N)C', 'C=[N+]([NH-])C'), ('C=N(=NC)C', 'C=[N+]([N-]C)C'),
    ('N=N(=N)C', 'N=[N+]([NH-])C'),
    ('[N-][N+](=O)C', 'N=[N+]([O-])C'), ('C[N-][N+](=O)C', 'CN=[N+]([O-])C'),
    ('[O-]N(=O)=O', '[O-][N+](=O)[O-]'),
    ('CN(:O):O', 'C[N+](=O)[O-]'),
    ('O=[N-]=O', '[O-]N=O'),
    ('O=N#N', 'O=[N+]=[N-]'), ('C=N#N', 'C=[N+]=[N-]'), ('N=N#N', 'N=[N+]=[N-]'),
    ('[O-][N+]#N', 'O=[N+]=[N-]'), ('C[CH-][N+]#N', 'CC=[N+]=[N-]'), ('[NH-][N+]#N', 'N=[N+]=[N-]'),
    ('C[N+]#N=[N-]', 'CN=[N+]=[N-]'),
    ('CN=N=N', 'CN=[N+]=[N-]'),
    ('CNN#N', 'CN=[N+]=[N-]'),
    ('[N-]#N=NC', '[N-]=[N+]=NC'),
    ('[N-]=N#N', '[N-]=[N+]=[N-]'),
    ('CC#N=N', 'CC=[N+]=[N-]'),
    ('CC#N=NC', 'CC#[N+][N-]C'), ('CC#N=O', 'CC#[N+][O-]'),
    ('NN#C', '[NH-][N+]#C'), ('ON#CC', '[O-][N+]#CC'), ('SN#CC', '[S-][N+]#CC'),
    ('CNN#C', 'C[N-][N+]#C'),
    ('N[N+]#[C-]', '[NH-][N+]#C'), ('O[N+]#[C-]', '[O-][N+]#C'), ('S[N+]#[C-]', '[S-][N+]#C'),
    ('CN[N+]#[C-]', 'C[N-][N+]#C'),
    ('CN#C', 'C[N+]#[C-]'),
    ('C[C-]=[N+]=N', 'CC=[N+]=[N-]'),
    ('CN(C)(C)=O', 'C[N+](C)(C)[O-]'), ('CN(C)(C)=NC', 'C[N+](C)(C)[N-]C'),
    ('C[N](C)=O |^1:1|', 'CN(C)[O] |^1:3|'),
    ('C=N(C)[O] |^1:3|', 'C=[N+](C)[O-]'),
    ('CN(C)#N', 'C[N+](C)=[N-]'),
    ('CN=[N+]', 'C[N+]#N'),
    ('[NH2+][C-]=O', 'N=C=O'), ('C[NH+][C-]=O', 'CN=C=O'),
    ('N#CO', 'N=C=O'),
    ('N#C[O-]', '[N-]=C=O'),
    ('CC(C)(C)[N+][O-]', 'CC(C)(C)N=O'),
    ('CN=O', 'C=NO'),
    ('NC=[O+]C', '[NH2+]=COC'), ('CNC=[O+]C', 'C[NH+]=COC'), ('CN(C)C=[O+]C', 'C[N+](C)=COC'),
    # amide rule: N=C-OH >> NH-C=O
    ('N=CO', 'NC=O'), ('N=CS', 'NC=S'),
    # ring amidation (6-membered): OH-C=N in ring >> O=C-NH in ring
    ('OC1=CC=CC=N1', 'O=C1NC=CC=C1'),
    ('OC1=CC=NC=N1', 'O=C1NC=NC=C1'),
    ('OC1=NC=CC=N1', 'O=C1N=CC=CN1'),
    ('OC1=C(O)N=CC=N1', 'O=C1NC=CNC1=O'),
    ('OC1=NC=CN=C1O', 'O=C1NC=CNC1=O'),
    ('OC1=CC=NC(=O)N1', 'O=C1NC=CC(=O)N1'),
    ('OC1=CC=NC(O)=N1', 'O=C1NC=CC(=O)N1'),
    ('OC1=NC(O)=NC=C1', 'O=C1NC=CC(=O)N1'),
    # 5-membered ring amidation (short flip)
    ('CN1C=CC(O)=N1', 'CN1NC(=O)C=C1'),
    # 5-membered ring amidation (long flip)
    ('CN1N=CC=C1O', 'CN1NC=CC1=O'),
    # hydroxypyridine to pyridone
    ('OC1=CC=NC=C1', 'O=C1C=CNC=C1'),
    # 6-membered ring N=C-CH adjacent to C=O (short flip)
    ('O=C1C=CCC=N1', 'O=C1NC=CC=C1'),
    # 6-membered ring N=C-C=C-CH adjacent to C=O (long flip)
    ('O=C1CC=CC=N1', 'O=C1NC=CC=C1'),
    # 6-membered ring N=C-CH adjacent to C=O (C=O between CH and ring end)
    ('O=C1CC=NC=C1', 'O=C1C=CNC=C1'),
    # 5-membered ring N=C-CH >> NH-C=C
    ('C1C=NC=N1', 'N1C=CN=C1'),
    # 6-membered ring exocyclic C=N to C-NH (cytosine-like)
    ('N=C1NC=CC(=O)N1', 'NC1=NC=CC(=O)N1'),
    # 5-membered ring N=C-CH with sp3 N,O closure
    ('CN1N=CCC1=O', 'CN1NC=CC1=O'),
    # acyclic enol (51)
    ('OC=C', 'O=CC'), ('OC(C)=C', 'O=C(C)C'),
    # P rules
    ('[O-][P+](C)(C)C', 'O=P(C)(C)C'),
    ('[CH2+][P-](C)(C)C', 'C=P(C)(C)C'),
    ('FP(F)(F)(F)(F)F', 'F[P-](F)(F)(F)(F)F'),
    ('C[P-](C)(=O)=O', 'CP(C)(=O)[O-]'),
    ('CP(C)(=O)[S-]', 'CP(C)(=S)[O-]'),
    ('CP(C)(=O)S', 'CP(C)(=S)O'),
    # S rules
    ('CS(=O)(=O)[S-]', 'CS(=O)(=S)[O-]'),
    ('CS(=O)(=O)S', 'CS(=O)(=S)O'),
    ('C[S+](C)[O-]', 'CS(C)=O'),
    ('O=[S+][O-]', 'O=S=O'), ('O=[S+](C)(C)[O-]', 'O=S(C)(C)=O'),
    ('O=[S+2]([O-])[O-]', 'O=S(=O)=O'),
    ('C[S+2](C)([O-])[O-]', 'CS(C)(=O)=O'),
    ('C[S-](C)[CH2+]', 'CS(C)=C'),
    ('O=[S-](C)=O', 'O=S(C)[O-]'),
    ('O=[S-](C)(=O)=O', 'O=S(C)(=O)[O-]'),
    ('O=[S-](=O)[S-]', 'S=S([O-])[O-]'),
    ('CS(C)(=O)O', 'CS(C)(=O)=O'),
    ('CS(C)(=O)N', 'CS(C)(=O)=N'),
    ('N=S(C)O', 'NS(C)=O'), ('N=S(C)(C)(C)O', 'NS(C)(C)(C)=O'),
    # C rules
    ('C#CO', 'C=C=O'),
    ('C#CNC', 'C=C=NC'),
    ('C=O |^1:0|', '[C-]#[O+]'),
    ('[O]O[O] |^1:0,2|', 'O=[O+][O-]'),
    ('[CH2+]N(C)C', 'C=[N+](C)C'),
    ('[CH2+]N(C)O', 'C=[N+](C)O'),
    ('[CH2+]=NC', 'C#[N+]C'),
    # Cl rules
    ('O[Cl+][O-]', 'OCl=O'),
    ('O[Cl+2]([O-])[O-]', 'OCl(=O)=O'),
    ('O[Cl+3]([O-])([O-])[O-]', 'OCl(=O)(=O)=O'),
    ('[Cl-]=O', 'Cl[O-]'),
    # S=N double rules
    ('OS(=N)(=N)O', 'O=S(N)(N)=O'), ('OS(=N)(=N)C', 'O=S(N)(=N)C'),
]

# The four inputs whose radical V2's reader invented from a short bracket valence and V3's does not, so
# `groups:01`/`groups:02` never match them.  That is a core reader question, not a table one, so these
# are xfailed rather than edited; when the reader ruling lands the strict xfails turn red.  Keyed by
# input, because the expectation is not what is in question.
NEEDS_INFERRED_RADICAL = frozenset((
    '[O]N(C)[NH]', '[O]N(C)[CH2]', '[O]S(C)(C)[O]', '[O]S(C)(C)[S]',
))

# The same four with the spin written out, which is what a V3 caller has to write.  These do fire, so
# `groups:01` and `groups:02` are gated even while the four above are xfailed.
SPELLED_RADICAL = [
    ('[O]N(C)[NH] |^1:0,3|', '[O-][N+](C)=N'),
    ('[O]N(C)[CH2] |^1:0,3|', '[O-][N+](C)=C'),
    ('[O]S(C)(C)[O] |^1:0,4|', 'O=S(C)(C)=O'),
    ('[O]S(C)(C)[S] |^1:0,4|', 'O=S(C)(C)=S'),
]

# Repairs correctly to a graph no SMILES comparison can express -- see
# `test_the_hypervalent_sulfur_product_is_right`.
NOT_A_SMILES_ROUNDTRIP = frozenset(('N=S(C)(C)(C)O',))

# The rules with no case here: `groups:44`, `groups:50`, the three rows written after the port
# (`gen_standardize_rules.ADDED`, gated by their own `examples` and by
# `test_standardize_overvalent_nitrogen.py` -- this file is V2's corpus and cannot grow a case for a rule
# V2 never had), and the metal table entire, since V2 shipped no metal-organic test.  Asserted so that a
# template edit shadowing a rule shows up as this set growing.  The metal ids are read from the table
# rather than written as a literal `range(22)` -- the collection was collapsed to 19 rows, and a stale
# range would fail here as three phantom rules.
UNREACHED = (frozenset(('groups:44', 'groups:50')) | frozenset(gen.ADDED)
             | frozenset(rule.id for rule in metals_rules()))


def _cases():
    for raw, expected in DATA:
        if raw in NEEDS_INFERRED_RADICAL:
            marks = [pytest.mark.xfail(strict=True, reason='V3 does not infer a radical from a '
                                       'short bracket valence; see SPELLED_RADICAL')]
        elif raw in NOT_A_SMILES_ROUNDTRIP:
            marks = [pytest.mark.xfail(strict=True, reason='product carries an H_UNKNOWN sulfur the '
                                       'reader would give 0 H; asserted atom by atom instead')]
        else:
            marks = []
        yield pytest.param(raw, expected, marks=marks, id=raw)


@pytest.mark.parametrize('raw,expected', list(_cases()))
def test_group(raw, expected):
    """One rule, one drawn-wrong molecule, one repaired answer.  V2's `test_group`, with V2's `==`
    replaced by the canonical form."""
    molecule = read_smiles(raw)
    log = molecule.log
    standardize(molecule)
    assert molecule == read_smiles(expected), (
        f'{raw} > {molecule.smiles} != {expected}  (fired: '
        f'{", ".join(record.rule for record in _repairs(log)) or "nothing"})')


@pytest.mark.parametrize('raw,expected', SPELLED_RADICAL, ids=[r for r, _ in SPELLED_RADICAL])
def test_group_with_the_radical_spelled_out(raw, expected):
    """The four radical cases as V3 obliges them to be written.  `groups:01` and `groups:02`."""
    molecule = read_smiles(raw)
    log = molecule.log
    standardize(molecule)
    assert molecule == read_smiles(expected), f'{raw} > {molecule.smiles} != {expected}'
    assert _repairs(log), f'{raw}: nothing fired, so the spelled radical did not reach the rule either'


def test_the_hypervalent_sulfur_product_is_right():
    """`N=S(C)(C)(C)O` >> `NS(C)(C)(C)=O`, asserted where a SMILES comparison cannot reach.

    The product's sulfur has no valence row for its environment, so `calc_implicit` answers
    `H_UNKNOWN` while the reader handed the same structure back answers 0 -- the molecule is therefore
    unequal to a re-read of its own SMILES.  That is a core disagreement between the two ways V3 gets
    a hydrogen count, not something `standardize()` can fix; it writes no hydrogen count at all.
    """
    molecule = read_smiles('N=S(C)(C)(C)O')
    assert standardize(molecule)

    heavy = {(atom.element, atom.charge, atom.is_radical) for atom in molecule.atoms()}
    assert heavy == {(6, 0, False), (7, 0, False), (8, 0, False), (16, 0, False)}
    hydrogens = sorted((atom.element, atom.implicit_h) for atom in molecule.atoms())
    assert hydrogens == [(6, 3), (6, 3), (6, 3), (7, 2), (8, 0), (16, None)], hydrogens

    # the graph is what V2 asked for; only the sulfur's hydrogen count is unanswerable
    expected = read_smiles('NS(C)(C)(C)=O')
    assert molecule.smiles == expected.smiles
    assert molecule != expected, ('the reader and `calc_implicit` now agree about a hypervalent '
                                  'sulfur; move this case back into DATA')


def test_the_case_count_is_v2s():
    """127 cases, as `483ddaa` has.  A silent loss of one is a rule losing its only test."""
    assert len(DATA) == 127
    assert len({raw for raw, _ in DATA}) == 127, 'an input appears twice, so one case is shadowed'


def _fired(raw, **kwargs):
    """`(molecule, [Rule])` -- the repaired molecule and the rules that repaired it, in order."""
    by_id = {rule.id: rule for rule in groups_rules() + metals_rules()}
    molecule = read_smiles(raw)
    log = molecule.log
    changed = standardize(molecule, **kwargs)
    return molecule, changed, [by_id[record.rule] for record in _repairs(log)]


def test_the_flag_withholds_exactly_the_flagged_rules():
    """`fix_tautomers=False` and no log record names a `tautomer` row, across all 127 inputs.

    The engine's half of the flag; which rows are withheld is the TSV's claim, tested below.
    """
    leaked = []
    for raw, _ in DATA:
        _, _, fired = _fired(raw, fix_tautomers=False)
        leaked += [rule.id for rule in fired if rule.tautomer]
    assert not leaked, f'flagged rules fired with the flag off: {sorted(set(leaked))}'


# The two inputs a general rule catches once its tautomer-picking neighbour is withheld, as
# `(specific, specific answer, general, general answer)`.
HANDED_OFF = {
    'C=N(=O)O': ('groups:11', 'C[N+](=O)[O-]', 'groups:12', 'C=[N+]([O-])O'),
    'CC#N=N': ('groups:27', 'CC=[N+]=[N-]', 'groups:28', 'CC#[N+][NH-]'),
}


@pytest.mark.parametrize('raw', sorted(raw for raw, _ in DATA
                                       if raw not in NEEDS_INFERRED_RADICAL))
def test_a_case_the_flag_withholds_is_repaired_or_left_alone_but_never_half_done(raw):
    """Every input, both ways round: no input is repaired by flagged and unflagged rules together.

    That partition (of the 123 inputs a rule reaches, 39 flagged, 84 unflagged, none mixed) is what
    stops the flag from leaving a molecule half-repaired.  With the flag off, 84 come out identical
    either way, 37 untouched, and the two in `HANDED_OFF` get a different but equally legal answer.
    """
    with_flag, _, fired = _fired(raw)
    without_flag, changed, _ = _fired(raw, fix_tautomers=False)

    if any(rule.tautomer for rule in fired):
        assert all(rule.tautomer for rule in fired), (
            f'{raw} is repaired by flagged and unflagged rules together '
            f'({", ".join(rule.id for rule in fired)}); the flag can now leave it half-repaired and '
            'this test no longer describes the collection')
        if raw not in HANDED_OFF:
            assert not changed, (f'{raw}: withholding the tautomer rules still changed it to '
                                 f'{without_flag.smiles}.  If an unflagged rule has legitimately '
                                 'taken over, add it to HANDED_OFF with both answers')
            assert without_flag == read_smiles(raw)
    else:
        assert without_flag == with_flag, (
            f'{raw}: the flag changed a repair no flagged rule performed '
            f'({with_flag.smiles} vs {without_flag.smiles})')


@pytest.mark.parametrize('raw', sorted(HANDED_OFF))
def test_a_general_rule_takes_over_when_the_tautomer_rule_is_withheld(raw):
    """The two molecules where both answers are correct: each is drawn with an illegal valence and a
    tautomer somebody may have meant, and only the second is the flag's business.

    `C=N(=O)O` is nitromethane drawn aci-nitro with a pentavalent N; `groups:11` returns nitromethane,
    withholding it leaves `groups:12` to charge-separate in place.  `CC#N=N` is the same shape:
    `groups:27` gives the diazo form, `groups:28` the nitrilimine.  Having an unflagged fallback is not
    a property of all 26 flagged rows -- `OS(=N)(=N)O` has none, and keeps its violation.
    """
    specific, specific_answer, general, general_answer = HANDED_OFF[raw]

    on, changed_on, fired_on = _fired(raw, fix_tautomers=True)
    assert changed_on
    assert [rule.id for rule in fired_on] == [specific]
    assert on == read_smiles(specific_answer), on.smiles

    off, changed_off, fired_off = _fired(raw, fix_tautomers=False)
    assert changed_off
    assert [rule.id for rule in fired_off] == [general]
    assert off == read_smiles(general_answer), off.smiles


def test_which_rules_the_gate_reaches():
    """92 of the 116 rules fire somewhere in here, and the 24 that do not are named in `UNREACHED`.

    22 of them are the metal table, which V2 never tested.  The other two are group rules: `groups:44`
    (pyrylium dearomatization) has no case, and `groups:50`'s own documented example is taken by
    `groups:47` before it can match, making it the collection's one rule shadowed on its stated input.
    Coverage of the 24 comes from the per-row `examples` cells run by
    `test_standardize_rules_examples.py`; this set is only a statement about the ported gate's reach.
    """
    fired = set()
    for raw, _ in DATA:
        molecule = read_smiles(raw)
        log = molecule.log
        standardize(molecule)
        fired |= {record.rule for record in _repairs(log)}
    for raw, _ in SPELLED_RADICAL:
        molecule = read_smiles(raw)
        log = molecule.log
        standardize(molecule)
        fired |= {record.rule for record in _repairs(log)}

    every = {rule.id for rule in groups_rules() + metals_rules()}
    assert fired <= every, f'a log named a rule no table declares: {sorted(fired - every)}'
    assert every - fired == UNREACHED, (
        'the set of rules this gate never reaches has changed.\n'
        f'  newly unreached: {sorted((every - fired) - UNREACHED)}\n'
        f'  newly reached:   {sorted(UNREACHED - (every - fired))}\n'
        'A rule that stops firing has been shadowed by an earlier one -- check the order before '
        'editing this set.')
