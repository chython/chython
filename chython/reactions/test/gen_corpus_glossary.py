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
"""Render `functional.tsv` and `protective.tsv` into `docs/glossary.rst`.

    python -m chython.reactions.test.gen_corpus_glossary

Writes between the two markers in that page and leaves the prose above them alone.  No verbs and no
options: the TSVs are the authority, the page is the output, and there is nothing to go the other way.
`test_corpus_glossary.py` fails when the page and the corpora disagree.
"""
import re
import sys
from pathlib import Path

from .._tables import functional_rules, protective_rules


ROOT = Path(__file__).resolve().parents[3]
PAGE = ROOT / 'docs' / 'glossary.rst'
BEGIN = '.. BEGIN GENERATED GLOSSARY: python -m chython.reactions.test.gen_corpus_glossary'
END = '.. END GENERATED GLOSSARY'

# A description spells a group name or a pattern in single backticks, which is Markdown's code span and
# rst's title reference.  Rewritten to a literal rather than escaped, so the TSV cell stays the readable
# thing a chemist edits.
CODE_SPAN = re.compile(r'`([^`]+)`')

# What rst reads as markup and the corpora do not write today: emphasis, a trailing reference, and a
# substitution.  Refused rather than escaped -- the moment one appears, the author gets told how to
# spell it instead of finding the page rendered wrong.
MARKUP = ('*', '_', '|')


def rst_text(text):
    """A description cell as rst: code spans become literals, and anything else markup-like is refused."""
    if text.count('`') % 2:
        raise ValueError(f'{text!r} reads as rst markup: an unpaired backtick')
    outside = CODE_SPAN.sub('', text)
    for ch in MARKUP:
        if ch in outside:
            raise ValueError(f'{text!r} reads as rst markup: a bare {ch!r} outside a code span.  Put the '
                             f'name or pattern in single backticks -- the generator turns those into rst '
                             f'literals')
    return CODE_SPAN.sub(r'``\1``', text)


def _literals(names):
    """A tuple of names as one cell: `` ``a``, ``b`` ``."""
    return ', '.join(f'``{n}``' for n in names)


def _table(widths, header, rows):
    """One `list-table`.  Cells are emitted verbatim, so a caller passes rst, not raw TSV text."""
    lines = ['.. list-table::',
             '   :header-rows: 1',
             '   :widths: ' + ' '.join(str(w) for w in widths),
             '']
    for cells in [header, *rows]:
        for i, cell in enumerate(cells):
            lines.append(('   * - ' if not i else '     - ') + cell)
    lines.append('')
    return lines


def compile_glossary():
    """The generated half of `docs/glossary.rst`, markers included."""
    functional = functional_rules()
    protective = protective_rules()

    lines = [BEGIN,
             '.. Generated from chython/reactions/tables/functional.tsv and protective.tsv, which are',
             '   the authority.  Do not edit below this line -- run the command above.',
             '']

    lines += ['Functional groups',
              '-----------------',
              '',
              f'{len(functional)} functional groups, alphabetically.  The name is the key',
              ':meth:`chython.MoleculeContainer.functional_groups` returns and the key',
              ':func:`chython.functional_rules` is keyed on; the id is what a consumer stores.',
              '']
    lines += _table((22, 12, 30, 36), ('Name', 'Id', 'SMARTS', 'What it matches'),
                    [(f'``{name}``', f'``{rule.id}``', f'``{rule.smarts}``', rst_text(rule.description))
                     for name, rule in sorted(functional.items())])

    lines += ['Protecting groups',
              '-----------------',
              '',
              f'{len(protective)} protecting groups, alphabetically.  *Protects* names the functional',
              'group the row reveals, which is a row of the table above.',
              ':meth:`chython.MoleculeContainer.protective_groups` reports a match and',
              ':meth:`chython.MoleculeContainer.deprotect` applies the patch.',
              '',
              'Two rows can match one substructure -- a Boc is also a tert-butyl -- and the larger',
              'pattern is served first, so a name here is the most specific group that fits, not',
              'every group that could.',
              '']
    lines += _table((22, 11, 15, 26, 26),
                    ('Name', 'Id', 'Protects', 'SMARTS', 'What it removes'),
                    [(f'``{name}``', f'``{rule.id}``', _literals(rule.protects), f'``{rule.smarts}``',
                      rst_text(rule.description))
                     for name, rule in sorted(protective.items())])

    lines.append(END)
    return '\n'.join(lines)


def rewrite_page(block, path=PAGE):
    text = path.read_text(encoding='utf-8')
    start = text.index(BEGIN)
    stop = text.index(END) + len(END)
    if text[start:stop] == block:
        return False
    path.write_text(text[:start] + block + text[stop:])
    return True


def main(argv):
    if argv:                                     # there is one direction, so there are no verbs
        print(__doc__)
        return 2
    if rewrite_page(compile_glossary()):
        print(f'{PAGE.name} updated')
    else:
        print(f'{PAGE.name} is already the corpora')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
