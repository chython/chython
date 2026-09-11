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
"""The R across every boundary at once: the full format chain, and the three CTfile rules with teeth."""
from importlib.resources import files
from chython import smiles
from chython.core import MoleculeContainer
from chython.core._core import element_symbols     # `chython.core`'s `__all__` does not name it
from chython.formats import mol


def test_the_index_survives_the_whole_chain():
    # SMILES -> bytes -> V2000 -> V3000 -> SMILES, which is every writer that knows about an R.
    start = smiles('[R1]c1ccc([R2])cc1')
    through_bytes = MoleculeContainer.from_bytes(start.to_bytes())
    through_v2000 = mol(mol(through_bytes, version=2000))
    through_v3000 = mol(mol(through_v2000, version=3000))
    assert str(through_v3000) == str(start)
    assert sorted(a.r_index for a in through_v3000.atoms() if a.is_r) == [1, 2]


def test_rgp_wins_over_the_symbol_column_and_says_so():
    text = ('probe\n  chython\n\n'
            '  2  1  0  0  0  0  0  0  0  0999 V2000\n'
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
            '    0.0000    0.0000    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0\n'
            '  1  2  1  0  0  0  0\n'
            'M  RGP  1   2   5\n'
            'M  END\n')
    log = []
    m = mol(text, log=log)
    assert next(a for a in m.atoms() if a.is_r).r_index == 5
    assert any('RGP' in record.message for record in log)


def test_an_alias_of_r1_over_a_carbon_stays_a_carbon():
    # An `A  aaa` alias is display text.  Promoting it would rewrite chemistry off a label.
    text = ('probe\n  chython\n\n'
            '  2  1  0  0  0  0  0  0  0  0999 V2000\n'
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
            '  1  2  1  0  0  0  0\n'
            'A    2\n'
            'R1\n'
            'M  END\n')
    m = mol(text)
    assert not any(a.is_r for a in m.atoms())
    assert m.brutto['C'] == 2


def test_a_record_with_sap_reads_and_logs_its_attachment_point():
    # `M  SAP` names an S-group attachment point, which no container field holds, so it is logged as an
    # unmodelled property rather than refused -- and the `*` beside it is the marker, not a query type.
    text = ('probe\n  chython\n\n'
            '  2  1  0  0  0  0  0  0  0  0999 V2000\n'
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n'
            '    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n'
            '  1  2  1  0  0  0  0\n'
            'M  SAP  1   2   1   1\n'
            'M  END\n')
    log = []
    m = mol(text, log=log)
    assert sum(a.is_r for a in m.atoms()) == 1
    assert m.atom_count == 2
    assert any(record.rule == 'v2000:unmodelled-property' for record in log), log


def test_elements_tsv_gains_no_row():
    # R is not an element.  It enters at `element_symbols()`, below the generated block.
    rows = [line for line in files('chython.core').joinpath('elements.tsv').read_text(encoding='utf-8').split('\n')
            if line and not line.startswith('#')]
    assert len(rows) - 1 == 118           # a header and 118 elements
    assert len(element_symbols()) == 119  # plus R at index 0
