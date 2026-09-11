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
"""The corpus reads, composes and stays internally consistent.

The tables are a seed, so most rows that will ever be here are added by someone reading these gates: each
fails with the edit that fixes it, naming the row, the unresolved group or the offending number.  What is
NOT testable here is whether a row is right about chemistry -- a row that composes, lexes and fires may
still describe a reaction that does not happen, and that is chemist review.
"""
import pytest
from .._tables import (SLOT_STRIDE, _attach_ring_sizes as _ring, _offset_map_numbers,
                       _parse_ring_sizes, compose_smirks, functional_rules, read_table,
                       reaction_rules)
from ...core import IncorrectSmirks, read_smarts, read_smiles, read_smirks


# The file-shaped gates hold for all three.  Everything further down is composition, which
# `protective.tsv` does not use -- its rows are whole SMIRKS, so its gates are in `test_protective.py`.
TABLES = ('functional.tsv', 'protective.tsv', 'reactions.tsv')


def _rows():
    """Every composed row, flat.  `reaction_rules()` is keyed by name and a name names a row FAMILY, so
    a gate that is about rows rather than about families flattens first."""
    return [rule for family in reaction_rules().values() for rule in family]

# --- the tables as files ------------------------------------------------------------------------

@pytest.mark.parametrize('table', TABLES)
def test_every_table_reads_and_has_rows(table):
    rows = read_table(table)
    assert rows, f'{table} has a header and nothing else'


@pytest.mark.parametrize('table', TABLES)
def test_ids_are_unique_within_a_table(table):
    ids = [row['id'] for row in read_table(table)]
    assert len(set(ids)) == len(ids), f'{table} repeats an id: a log record could name two rows'


def test_the_merge_left_one_id_space():
    """One id space, numbered 1..n.  An id is positional, so a gap or a repeat is a dropped line."""
    ids = [int(row['id']) for row in read_table('reactions.tsv')]
    assert ids == list(range(1, len(ids) + 1)), 'reactions.tsv ids are not 1..n in order'
    assert sum(len(family) for family in reaction_rules().values()) == len(ids)


def test_a_functional_group_is_named_once():
    """The composer resolves a slot by name, so a duplicate name would silently shadow a pattern."""
    names = [row['name'] for row in read_table('functional.tsv')]
    assert len(set(names)) == len(names)


# --- composition --------------------------------------------------------------------------------

def test_the_slot_stride_is_what_the_tables_were_written_against():
    """A constant the table rows encode by hand, so changing it invalidates every product column."""
    assert SLOT_STRIDE == 100


def test_offsetting_moves_map_numbers_and_leaves_aromatic_bonds_alone():
    """`:` is two tokens and only the one inside a bracket is a map number.

    An aromatic bond with a ring-closure digit after it (`[C;a]:1:[C;D2]...`) is what a naive scan
    corrupts.
    """
    assert _offset_map_numbers('[C:1][N:2]', 100) == '[C:101][N:102]'
    assert _offset_map_numbers('[C;a:1]:[C;a:2]', 100) == '[C;a:101]:[C;a:102]'
    assert _offset_map_numbers('[C;a:1]:1:[C;D2]:[C;D2]:1', 200) == '[C;a:201]:1:[C;D2]:[C;D2]:1'
    assert _offset_map_numbers('[C:1]', 0) == '[C:1]', 'slot 0 is the identity'


def test_a_group_keeps_its_own_numbers_in_slot_zero():
    """Which is the property that lets a group gain an atom without moving anybody else's numbers."""
    smirks = compose_smirks(('carboxylic_acid', 'primary_amine'), '[A:1](=[A:2])-[A:101]-[A:102]')
    reactants = smirks.split('>>')[0]
    assert functional_rules()['carboxylic_acid'].smarts in reactants


def test_the_second_slot_is_offset_by_the_stride_and_not_by_the_first_group():
    """The offset is a constant of the SLOT, so `primary_amine`'s atoms are 101 and 102 whether the group
    before it has two atoms or twenty."""
    short = compose_smirks(('carboxylic_acid', 'primary_amine'), '[A:1]')
    long = compose_smirks(('vicinal_diol', 'primary_amine'), '[A:1]')
    assert '[N;D1;z1;x0:101]' in short and '[N;D1;z1;x0:101]' in long


def test_an_unknown_group_names_itself():
    with pytest.raises(ValueError) as exc:
        compose_smirks(('no_such_group',), '[A:1]')
    assert 'no_such_group' in str(exc.value)


def test_the_components_are_grouped_and_never_joined_unconstrained():
    """`(A).(B)` or `(A.B)`, and never a bare `.`.

    `.` says only "not bonded", so an unconstrained join also matches the two groups BONDED TOGETHER in
    one molecule -- a cyclization with no ring size stated, firing on a four-membered lactam as happily
    as on a forty-membered one.
    """
    inter = compose_smirks(('aryl_bromide', 'aryl_boronic_acid'), '[A:1]-[A:101]').split('>>')[0]
    intra = compose_smirks(('aryl_bromide', 'aryl_boronic_acid'), '[A:1]-[A:101]',
                           intramolecular=True).split('>>')[0]
    assert inter.startswith('(') and ').(' in inter, inter
    assert intra.startswith('(') and intra.endswith(')') and ').(' not in intra, intra
    # and one slot is not grouped at all, there being nothing to constrain
    assert not compose_smirks(('aryl_bromide',), '[A:1]').startswith('(')


def test_the_two_groupings_are_mutually_exclusive():
    """What stops a row that composed both templates from yielding one product twice.

    A claim about the core's component groups: if `(A).(B)` ever also matched one component, every
    `ring_sizes` row would double every intramolecular outcome.
    """
    acid_amine = ('carboxylic_acid', 'primary_amine')
    product = '[A:1](=[A:2])-[A:101]-[A:102]'
    inter = read_smirks(compose_smirks(acid_amine, product))
    intra = read_smirks(compose_smirks(acid_amine, _ring(product, 1, (5, 6, 7)),
                                       intramolecular=True))

    together = read_smiles('NCCCC(=O)O')            # 4-aminobutanoic acid: both groups, one component
    apart = read_smiles('CC(=O)O.CCN')              # the same two groups, two components
    assert not list(inter(together)) and len(list(intra(together))) == 1
    assert len(list(inter(apart))) == 1 and not list(intra(apart))


# --- the ring-size column -----------------------------------------------------------------------

def test_the_ring_sizes_column_composes_a_product_side_r():
    """`1:5,6,7` -> `[A;r5,r6,r7:1]`, and the disjunction lands on the atom the cell names."""
    assert _ring('[A:1](=[A:2])-[A:101]', 1, (5, 6, 7)) == '[A;r5,r6,r7:1](=[A:2])-[A:101]'
    assert _ring('[A:1](=[A:2])-[A:101]', 101, (5,)) == '[A:1](=[A:2])-[A;r5:101]'


def test_the_map_number_is_not_part_of_the_and_or_grammar():
    """WHY THE `r` DISJUNCTION MAY GO LAST IN THE BRACKET.

    A high AND after a `,` list binds to that list's LAST alternative only.  A map number is not a
    primitive the lexer ANDs -- it is an atom-level field -- so `[C;r5,r6:1]` is a 5-or-6-ring carbon
    numbered 1.  If that changes, every `ring_sizes` row silently loses its numbering on all but the
    last size, so it is pinned here rather than trusted.
    """
    for smarts in ('[C;r5,r6:1]', '[C;r5,r6;D2:1]', '[C;D2;r5,r6:1]'):
        query = read_smarts(smarts)
        assert query.map_numbers() == {1: 1}, smarts
        assert query.box_counts() == [2], f'{smarts} should be two boxes, one per ring size'


def test_a_ring_size_outside_the_primitives_range_is_refused():
    """`r` spans 3-14, so a 20 composes a template that reads and can never match -- refused at load,
    because a row that never fires is invisible."""
    with pytest.raises(ValueError) as exc:
        _parse_ring_sizes('1:20', 'reactions:99')
    assert '20-membered' in str(exc.value) and 'reactions:99' in str(exc.value)


def test_a_ring_atom_the_product_does_not_mention_is_refused():
    with pytest.raises(ValueError) as exc:
        _ring('[A:1]-[A:101]', 7, (5,))
    assert ':7' in str(exc.value)


def test_a_malformed_ring_sizes_cell_names_itself():
    with pytest.raises(ValueError) as exc:
        _parse_ring_sizes('5,6,7', 'reactions:99')
    assert '<atom>:<sizes>' in str(exc.value)
    assert _parse_ring_sizes('', 'reactions:99') == (0, ()), 'an empty cell is the common case'


def test_a_one_slot_row_may_not_name_a_ring():
    """The one thing a row's slot count still decides.

    A one-slot row is already one molecule, so there is no intermolecular reading for `ring_sizes` to
    distinguish it from and the second template would copy the first.  The loader refuses the column on
    such a row and says to write the `r` into `product` instead.
    """
    for rule in _rows():
        if len(rule.groups) == 1:
            assert not rule.ring_sizes, f'{rule.id} is one slot and names a ring'
            assert rule.intramolecular is None, f'{rule.id} composed an intramolecular template'


def test_every_intramolecular_template_is_a_multi_slot_row_with_a_ring():
    """The two fields agree, per row, in both directions -- no template without a ring, none without."""
    for rule in _rows():
        assert (rule.intramolecular is not None) == bool(rule.ring_sizes), rule.id
        if rule.ring_sizes:
            assert len(rule.groups) >= 2, rule.id
            assert rule.ring_atom, f'{rule.id} has sizes and no ring atom'


# --- every row, compiled ------------------------------------------------------------------------

def test_every_row_composes_and_lexes():
    """`reaction_rules` raises while composing, so this passing means every row produced a sealed
    template -- a bad row fails here and not on a molecule."""
    for rule in _rows():
        assert rule.template is not None
        assert rule.template.rule_id == rule.id


def test_a_rows_rule_id_is_table_qualified():
    """`reactions:13`, so a log record names the row a chemist can edit and not the composed SMIRKS
    string nobody wrote.  Qualified like `chemistry/`'s ids, so a log mixing the two stays readable."""
    for rule in _rows():
        table, _, number = rule.id.partition(':')
        assert table == 'reactions', rule.id
        assert number.isdigit(), rule.id


def test_both_arities_live_in_one_table():
    """One-slot and multi-slot rows are the same kind of thing and `react()` reaches both, so a table
    holding only one of them would be an arity filter in a filename."""
    slots = {len(rule.groups) for rule in _rows()}
    assert 1 in slots and 2 in slots, f'reactions.tsv holds only {slots}-slot rows'


def test_every_product_number_comes_from_a_slot():
    """A product number no reactant produces builds an atom that pairs with nothing.

    Legal in the notation, but always a typo in this corpus, where a created atom is written without a
    number.  So every `:N` on a product side must be `slot * 100 + m` for an atom `m` of that slot.
    """
    available = {}
    for name, group in functional_rules().items():
        available[name] = {int(n) for n in _map_numbers(group.smarts)}
    for rule in _rows():
        allowed = {slot * SLOT_STRIDE + n
                   for slot, group in enumerate(rule.groups) for n in available[group]}
        used = {int(n) for n in _map_numbers(rule.product)}
        assert used <= allowed, (
            f'{rule.id} ({rule.name}) references {sorted(used - allowed)}, which no slot produces. '
            f'Slot numbers available: {sorted(allowed)}.  Composed: {rule.template.smirks}')


#: Reaction rows that delete a NUMBERED reactant atom on purpose, and which numbers.  Keyed by name
#: because a row's chemistry is what justifies the deletion: `appel` swaps an alcohol's oxygen for a
#: halide, `nitro_to_amine` drops both nitro oxygens.  Every other row must keep every atom its groups
#: number -- see `test_a_numbered_reactant_atom_is_kept_unless_the_row_says_otherwise`.
_DELIBERATE_DELETIONS = {
    'appel': {1},                        # the alcohol oxygen, replaced by the halide
    'appel_chloride': {1},
    'amide_to_amine': {3},               # the carbonyl oxygen, reduced away
    'nitro_to_amine': {2, 3},            # both nitro oxygens
    'sulfoxide_to_thioether': {2},       # the sulfoxide oxygen
    'amidation': {3},                    # the carboxylic acid hydroxyl O, the leaving group
    'esterification': {3},
    'weinreb_amidation': {3},
    'hydrazide_formation': {3},
    'acid_chlorination': {3},            # the carboxylic acid hydroxyl O, replaced by Cl
    'acid_to_alcohol': {3},              # the carboxylic acid hydroxyl O, released as water
    'deoxygenative_coupling': {1},       # the alcohol oxygen, the whole point of the coupling
    'decarboxylative_coupling': {1, 2, 3},  # the carboxyl carbon and its oxygen, leaving as CO2 -- which
                                            # numbers those are is per group: an acid numbers them 1 and 2,
                                            # a redox-active ester 2 and 3
    'mitsunobu': {1},                    # the alcohol oxygen: the nucleophile supplies the one that stays,
                                         # and losing this one with inversion is the reaction
    'reductive_amination': {1},          # the carbonyl oxygen, leaving as water before the reduction
    'knoevenagel': {1},                  # the carbonyl oxygen, leaving as water
    'hwe': {1},                          # the carbonyl oxygen, leaving on the phosphorus
    'wittig': {1},
    'ugi_4cr': {1, 3},                   # the aldehyde oxygen and the acid hydroxyl, both as water
    # A ring former condenses out every oxygen the new ring does not need.  The acid keeps neither of its
    # two, the 1,4-diketone neither of its two, and the oxadiazole's ring O comes from the hydrazide.
    'imidazopyridine': {1, 3},
    'benzimidazole': {1, 2, 3},
    'paal_knorr': {1, 6},
    'oxadiazole': {2, 3},
}


def test_a_numbered_reactant_atom_is_kept_unless_the_row_says_otherwise():
    """The other direction of the rule above, and the one that breaks at a distance.

    `functional.tsv` leaves a leaving group UNNUMBERED, so a numbered reactant atom is one the group
    means to keep -- but a product side deletes by absence, so keeping it is the product side's job.
    Adding an atom to a group therefore changes what every row referencing that group deletes, silently:
    numbering `tertiary_amine`'s three substituents turned `nitrogen_oxidation` on trimethylamine from
    the N-oxide into `[NH3+][O-]`, and no test named that product.

    A row that does mean to drop a numbered atom says so in `_DELIBERATE_DELETIONS` with the chemistry
    beside it.
    """
    available = {name: {int(n) for n in _map_numbers(group.smarts)}
                 for name, group in functional_rules().items()}
    for rule in _rows():
        reactant = {slot * SLOT_STRIDE + n
                    for slot, group in enumerate(rule.groups) for n in available[group]}
        kept = {int(n) for n in _map_numbers(rule.product)}
        dropped = {n % SLOT_STRIDE for n in reactant - kept}
        assert dropped <= _DELIBERATE_DELETIONS.get(rule.name, set()), (
            f'{rule.id} ({rule.name}) silently deletes reactant atom(s) {sorted(dropped)} of '
            f'{rule.groups}: they are numbered, so the groups mean to keep them, and this product side '
            f'does not restate them.  Composed: {rule.template.smirks}')


def _map_numbers(smarts):
    """Every `:N` inside a bracket, as strings.  The scanner `_offset_map_numbers` uses, read-only."""
    out = []
    depth = 0
    i = 0
    while i < len(smarts):
        c = smarts[i]
        if c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
        elif c == ':' and depth > 0:
            j = i + 1
            while j < len(smarts) and smarts[j].isdigit():
                j += 1
            if j > i + 1:
                out.append(smarts[i + 1:j])
                i = j
                continue
        i += 1
    return out


# --- the key structure each accessor hands out --------------------------------------------------
#
# Three tables, three key structures, and that is the DATA rather than an inconsistency: a functional
# group's name names one row, a protecting group's name names one row, a reaction's name names a row
# FAMILY.  Collapsing them to one shape would lose the last fact.

def test_reaction_rules_groups_a_family_under_its_name():
    """`reaction=` selects by name, so the name is the key and its value is every spelling.

    294 rows under 72 names -- `amidation` is three, one per way the acid is activated -- so a
    name-to-row map would drop rows and a flat tuple makes every `reaction=` call a linear scan that
    rebuilds the known-name set to say what it did not find.
    """
    rules = reaction_rules()
    assert isinstance(rules, dict)
    assert all(isinstance(rows, tuple) and rows for rows in rules.values())
    assert all(rule.name == name for name, rows in rules.items() for rule in rows)
    assert sum(len(rows) for rows in rules.values()) == len(read_table('reactions.tsv'))
    assert len(rules['amidation']) > 1, 'the case a name-keyed row map could not hold'


def test_reaction_rules_keeps_table_order_inside_a_family_and_across_them():
    """Insertion order, so the dict IS the table read in order and the ids inside a family ascend."""
    rows = read_table('reactions.tsv')
    rules = reaction_rules()
    seen = []
    for row in rows:
        if row['name'] not in seen:
            seen.append(row['name'])
    assert list(rules) == seen
    for name, family in rules.items():
        ids = [int(rule.id.split(':')[1]) for rule in family]
        assert ids == sorted(ids), name


# --- laziness and caching -----------------------------------------------------------------------

def test_the_composed_templates_are_cached():
    """One `read_smirks` per row per process, pinned by identity -- the only thing that distinguishes a
    cache from a fast recompilation."""
    assert reaction_rules() is reaction_rules()
    assert reaction_rules()['amidation'][0].template is reaction_rules()['amidation'][0].template


def test_the_groups_are_cached_too():
    assert functional_rules() is functional_rules()


def test_the_table_accessor_and_the_pass_do_not_share_a_name():
    """`functional_rules()` reads the table; `functional_groups(molecule)` asks a molecule.  Two
    questions, and one name for both is a collision a caller resolves by import order."""
    from chython.reactions import functional_groups, functional_rules, roles

    table = functional_rules()
    assert isinstance(table, dict) and 'carboxylic_acid' in table
    assert table['carboxylic_acid'].id.startswith('functional:')

    hits = functional_groups(read_smiles('CC(=O)O'))
    assert hits['carboxylic_acid'] == 1

    assert 'aryl_halide' in roles()


def test_nothing_loads_at_import():
    """A process that only writes SMILES must not pay for the corpus.  Checked in a fresh subprocess,
    since a module-level compile would fill the caches before any accessor ran."""
    from subprocess import run
    from sys import executable
    script = ('import chython.reactions._tables as t;'
              'print(len(t._FUNCTIONAL_CACHE), t._RULES_CACHE, t._PROTECTIVE_CACHE)')
    result = run([executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ['0', 'None', 'None'], (
        'importing the loader compiled a table: ' + result.stdout)


# --- what a bad row does ------------------------------------------------------------------------

def test_a_product_side_that_cannot_be_read_is_refused_at_load_time():
    """And not on the first molecule.  `*` is the narrowest example: neither build nor check.

    Composing eagerly per table proves the whole table readable the first time any of it is used.
    """
    with pytest.raises(IncorrectSmirks):
        read_smirks(compose_smirks(('aryl_bromide',), '[A;*:1]'))
