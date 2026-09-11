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
"""No comment in `chython/formats` may cite a review document the reader cannot open.

"see finding 3" reads like a citation with nothing at the other end; name the behaviour instead.  A
number indexing something real -- an atom, a column, a spec page -- is not matched.  Known false
positive: `round` before a number.  Write "to 4 decimal places" with no bare `round` in front of it.
"""

from pathlib import Path
from re import IGNORECASE, finditer


#: Words that, immediately followed by a number, name a position in a review or planning document.
#: The optional plural is load-bearing: `rulings 2 and 3` reads right past a singular-only pattern.
_BOOKKEEPING = r'\b(?:ruling|finding|round|task|plan|wave)s?[ -]#?[0-9]'

#: Fewer Python files than this under `chython/formats` means the sweep found the wrong tree.  Well
#: below the count today (32), and far enough above zero to catch a scan of nothing.
_MINIMUM_SOURCES = 20


def _sources():
    """Every Python file under `chython/formats`, except this one, which has to quote the words."""
    root = Path(__file__).resolve().parent.parent  # chython/formats/
    here = Path(__file__).resolve()
    return sorted(p for p in root.rglob('*.py') if p.resolve() != here)


def test_no_source_references_a_review_document():
    """Nothing under `chython/formats`, tests included, cites a round, task, plan or finding number."""
    offenders = []
    for path in _sources():
        for n, line in enumerate(path.read_text(encoding='utf8').splitlines(), 1):
            for hit in finditer(_BOOKKEEPING, line, IGNORECASE):
                offenders.append(f'{path.name}:{n}: {hit.group()!r} in {line.strip()[:70]}')
    assert not offenders, ('references to review bookkeeping -- name the behaviour instead:\n'
                           + '\n'.join(offenders))


def test_the_pattern_fires():
    """A control: the pattern catches the phrasings that have appeared, and nothing else.

    Without it, a typo in the regex leaves the scan above green while matching nothing.
    """
    caught = ['the ruling-2 branch concatenates it onto the name',
              '# Finding 6: CTAB body lines do not arm continuation',
              'superseded in round 3',
              'see task 11',
              "the plan's wave-0 preamble",
              'Plan #4 covers this',
              # Plurals: the singular-only pattern read straight past this line.
              '# superseded by rulings 2 and 3; see also findings 4 and 5, and rounds 2 and 3',
              # An acknowledged false positive, kept visible: prose about decimal places, not a
              # review round.  The pattern cannot tell them apart; see the module docstring.
              'round 4 decimal places, after finding 2 non-numeric fields']
    for probe in caught:
        assert next(finditer(_BOOKKEEPING, probe, IGNORECASE), None), probe

    # Numbers that index something the reader can actually look at.
    allowed = ['to 4 decimal places, with no bare word in front of the number',
               'CTfile specification p.46 gives the rule',
               'atom 3 carries the wedge',
               'ruling F26 names an invariant, not a review round',
               'V3000 counts line, field 2',
               'a plan for the value: concatenate, do not replace']
    for probe in allowed:
        assert next(finditer(_BOOKKEEPING, probe, IGNORECASE), None) is None, probe


def test_the_scan_covers_the_package():
    """A floor under the file set: a count, plus one production and one test file known to be there.

    Move this file or restructure the package and `_sources` can come back nearly empty while the
    other two tests stay green.  The test-file check catches a glob that misses `test/`.
    """
    names = {p.name for p in _sources()}
    assert len(_sources()) >= _MINIMUM_SOURCES, sorted(names)
    assert '_v2000.py' in names, sorted(names)
    assert 'test_log_prefix.py' in names, sorted(names)
