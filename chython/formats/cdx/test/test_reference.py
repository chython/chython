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
"""What we write, read by Revvity's open-source ChemDraw library and handed back in the other spelling.

| Chain | Must hold |
|---|---|
| our CDX bytes -> the library -> its CDXML -> chython | the record we wrote |
| our CDXML text -> the library -> its CDX -> chython | likewise |

The CDXML half of the first chain is itself checked against RDKit's reader in ``test_interop.py``, so a
pass here says the library accepts our binary as we meant it.  Two things are left out, each because the
library as built by the recipe below loses it on its own files too:

| Not compared | Library behaviour, on Revvity-written files as on ours |
|---|---|
| ``GenericNickname`` (R markers) | read back as empty |
| text, CDX -> CDXML chain | each text written as its run of suffixes: ``COOH`` becomes ``OOHOHH`` |

THE ORACLE IS OPTIONAL.  ``$CHYTHON_CDXCONVERT`` names the ``cdxconvert`` binary built from
``cdxconvert.cpp`` beside this module; unset skips, set to something unusable fails.  To build it,
with CMake and Boost (``brew install boost``) at hand::

    git clone https://github.com/Glysade/chemdraw && cd chemdraw
    # CMakeLists.txt as of 2026-06: line 24 reads `else(WORDS_BIGENDIAN)CM`; drop the `CM`.  The
    # ChemDraw target must also link `expatpp EXPAT::EXPAT` and include
    # `expatpp/expatpp-code-r6-trunk/src_pp`.
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
    clang++ -std=c++17 -DTARGET_API_LIB -D__linux -DHAVE_EXPAT_CONFIG_H -DPLATFORM_LITTLEENDIAN \\
        -Ichemdraw -Iexpatpp/expatpp-code-r6-trunk/src_pp -I"$(brew --prefix boost)/include" \\
        <chython>/chython/formats/cdx/test/cdxconvert.cpp -Lbuild -lChemDraw -Lbuild/expatpp -lexpatpp \\
        -lexpat -Wl,-rpath,"$PWD/build" -Wl,-rpath,"$PWD/build/expatpp" -o cdxconvert
    export CHYTHON_CDXCONVERT="$PWD/cdxconvert"
"""

from os import access, environ, X_OK
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from pytest import fail, fixture, mark, skip

from chython.core import smiles
from chython.core.reaction import ReactionContainer

from .. import cdx, cdxml
from ...ctfile import mol
from .test_cdx import CHFClBr_V3000, chfclbr, document, esterification


@fixture(scope='module')
def convert():
    binary = environ.get('CHYTHON_CDXCONVERT')
    if not binary:
        skip('cdxconvert oracle absent; set $CHYTHON_CDXCONVERT (build recipe in this module)')
    if not access(binary, X_OK):
        fail(f'$CHYTHON_CDXCONVERT={binary} is not an executable: a provisioned oracle that cannot run '
             f'answers nothing')

    def through(data, suffix_in, suffix_out):
        with TemporaryDirectory() as d:
            source, target = Path(d, 'in' + suffix_in), Path(d, 'out' + suffix_out)
            source.write_bytes(data if isinstance(data, bytes) else data.encode())
            done = run([binary, str(source), str(target)], capture_output=True, text=True, timeout=60)
            assert done.returncode == 0, f'the reference library refused our {suffix_in}: {done.stderr}'
            return target.read_bytes()
    return through


def _back(convert, record):
    """The record through both chains, as chython reads each result."""
    return (cdxml(convert(cdx(record), '.cdx', '.cdxml')),
            cdx(convert(cdxml(record), '.cdxml', '.cdx')))


MOLECULES = {
    'wedged centre': lambda: mol(CHFClBr_V3000),
    'hashed wedge': lambda: cdxml(chfclbr(display='Display="WedgedHashBegin"')),
    'OR group': lambda: cdxml(chfclbr().replace('<n id="4" p="0 0" />',
                                                '<n id="4" p="0 0" EnhancedStereoType="Or" '
                                                'EnhancedStereoGroupNum="2"/>')),
    'charge, isotope, radical, stated H': lambda: cdxml(document(
        [(4, 0, 0, 'Element="7" Charge="1" NumHydrogens="3"'), (5, 30, 0, 'Isotope="13"'),
         (6, 60, 0, 'Radical="Doublet"')], [(7, 4, 5, ''), (8, 5, 6, '')])),
    'aromatic ring': lambda: cdxml(document([(i, 30 * i, 0, '') for i in range(4, 10)],
                                            [(10 + i, 4 + i, 4 + (i + 1) % 6, 'Order="1.5"') for i in range(6)])),
    'dative to a metal': lambda: cdxml(document([(4, 0, 0, 'Element="78"'), (5, 30, 0, 'Element="7"')],
                                                [(6, 4, 5, 'Order="dative"')])),
    'two components': lambda: cdxml(document([(4, 0, 0, 'Element="11" Charge="1"'),
                                              (5, 60, 0, 'Element="17" Charge="-1"')])),
}


@mark.parametrize('name', sorted(MOLECULES))
def test_the_reference_library_reads_our_molecule(convert, name):
    record = MOLECULES[name]()
    for back in _back(convert, record):
        assert str(back) == str(record)


def test_the_reference_library_reads_our_reaction(convert):
    record = cdxml(esterification())
    via_cdx, via_cdxml = _back(convert, record)
    assert via_cdxml.meta == record.meta        # text: the CDXML -> CDX chain only, see the module
    for back in (via_cdx, via_cdxml):
        assert str(back) == str(record)
        assert back.reactants[0].atom(2).map_number == back.products[0].atom(2).map_number == 1
        assert back.reactants[1].atom(2).map_number == back.products[0].atom(4).map_number == 2


def test_the_reference_library_keeps_a_salt_participant_whole(convert):
    record = ReactionContainer([smiles('C[NH3+].[Cl-]')], [smiles('CN')], [])
    for back in _back(convert, record):
        assert len(back.reactants) == 1 and str(back) == str(record)
