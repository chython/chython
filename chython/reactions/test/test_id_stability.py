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
"""The four corpora's ids, pinned to the rows they name.

An id is what a consumer persists -- `functional_group_hits()` hands it out for exactly that -- so a
renumber corrupts stored data the way a moved bit index would.  `test_tables.py` checks id FORMAT
(`table:number`); this checks id IDENTITY.

ADDING A ROW IS FREE.  To record new rows, regenerate the frozen file:

    python -c "from chython.reactions.test.test_id_stability import regenerate; regenerate()"

and read the diff -- a new line is a new row, a CHANGED line is the failure this file exists to make
visible.
"""
from pathlib import Path
from chython.reactions import functional_rules, protective_rules, reaction_rules, roles


#: The four corpora, in the order `_frozen_ids.py` writes them: table file, frozen dict name.
TABLES = (('functional.tsv', 'FUNCTIONAL'), ('protective.tsv', 'PROTECTIVE'), ('roles.tsv', 'ROLES'),
          ('reactions.tsv', 'REACTIONS'))


def _live() -> dict[str, dict[str, str]]:
    """`{table: {id: label}}`, where a LABEL identifies the row and not merely its family.

    A name where a name is unique, a name plus its discriminator where it is not: `roles.tsv` names one
    role over several groups and `reactions.tsv` one reaction over several spellings, so a label of the
    name alone would let two rows of one family trade ids unseen.
    """
    return {'functional.tsv': {g.id: g.name for g in functional_rules().values()},
            'protective.tsv': {r.id: r.name for r in protective_rules().values()},
            'roles.tsv': {r.id: '%s/%s' % (r.name, r.group)
                          for rows in roles().values() for r in rows},
            'reactions.tsv': {r.id: '%s/%s' % (r.name, '+'.join(r.groups))
                              for rows in reaction_rules().values() for r in rows}}


def _pairs() -> list[tuple[str, dict[str, str], dict[str, str]]]:
    """`(table, live, frozen)` per corpus.

    The frozen module is imported HERE and not at file scope so `regenerate()` runs on a tree where
    `_frozen_ids.py` does not exist yet -- which is the state the docstring's command is for.
    """
    from . import _frozen_ids
    live = _live()
    return [(table, live[table], getattr(_frozen_ids, name)) for table, name in TABLES]


def test_no_frozen_id_names_a_different_row_now():
    for table, live, frozen in _pairs():
        moved = {i: (was, live[i]) for i, was in frozen.items() if i in live and live[i] != was}
        assert not moved, (f'{table}: {len(moved)} id(s) name a different row than they did: {moved}. An '
                           'id is persisted identity -- give the new row a new id and leave this one '
                           'where it is, or a stored column means something other than what it recorded.')


def test_no_frozen_id_disappeared():
    for table, live, frozen in _pairs():
        gone = {i: n for i, n in frozen.items() if i not in live}
        assert not gone, (f'{table}: {len(gone)} frozen id(s) are absent: {gone}. Deleting a row leaves '
                          'stored data pointing at nothing -- retire the row by emptying its pattern and '
                          'keeping its id, or record the retirement in _frozen_ids.py deliberately.')


def test_every_live_id_is_unique_within_its_table():
    """A dict cannot hold a collision, so the compiled count is compared against the row count."""
    counts = {'functional.tsv': len(functional_rules()), 'protective.tsv': len(protective_rules()),
              'roles.tsv': sum(len(rows) for rows in roles().values()),
              'reactions.tsv': sum(len(rows) for rows in reaction_rules().values())}
    for table, live, _ in _pairs():
        assert len(live) == counts[table], f'{table} compiled two rows onto one id'


def test_no_id_is_shared_across_the_four_tables():
    """What makes an id storable on its own, with no column saying which corpus it came from."""
    seen = {}
    for table, live, _ in _pairs():
        for i in live:
            assert i not in seen, f'{i} is in both {seen[i]} and {table}'
            seen[i] = table


# --- the generator the docstring names ------------------------------------------------------------

_HEADER = '''# -*- coding: utf-8 -*-
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
"""GENERATED.  Every rule id of the four corpora, mapped to the row it names.

Read by `test_id_stability.py` and by nothing else.  ADDING A ROW IS FREE: regenerate with

    python -c "from chython.reactions.test.test_id_stability import regenerate; regenerate()"

and read the diff.  A new line is a new row; a CHANGED line is a renumber, which moves what a stored id
means.  Never hand-edit a line to make the ratchet pass.
"""
'''


def _sort_key(rule_id: str) -> tuple[str, int, str]:
    """By the numeric part, so `functional:9` precedes `functional:10` and an append diffs as a block."""
    table, _, tail = rule_id.partition(':')
    return (table, int(tail), '') if tail.isdigit() else (table, 0, tail)


def regenerate():
    """Rewrite `_frozen_ids.py` from the live tables.  Called by hand; read the diff afterwards."""
    out = [_HEADER]
    for table, label in TABLES:
        live = _live()[table]
        out.append('#: `%s`: %d rows, id -> the row it names.\n%s = {\n%s}\n'
                   % (table, len(live), label,
                      ''.join('    %r: %r,\n' % (i, live[i]) for i in sorted(live, key=_sort_key))))
    (Path(__file__).parent / '_frozen_ids.py').write_text('\n'.join(out))
