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
"""Which package the façade's names actually came from.

`chython/__init__.py` does `from .core import *` and then `from .formats import *`.  On a name collision
the second import wins **silently** -- no error, no warning, and `__all__` is empty by design so there is
nothing to diff.  A name that exists in both packages is therefore served from `formats` and nobody finds
out.  Measured today: `mol` and `rxn` resolve to `chython.formats.ctfile._facade`, which is correct and
intended.  This file's job is to notice the day that changes, or the day a `core` name starts arriving
from `formats` because someone added a same-named function one layer up.
"""

from pytest import mark


#: `(attribute, defining module)` -- the module the façade's copy of the name must come from.
#: `inchi`/`inchikey` are `core`'s: an InChI is a string, not a file, so `formats` never sees them.
_OWNERS = [('inchi', 'chython.core._core'),
           ('inchikey', 'chython.core._core'),
           # the core's own bidirectional doors.  `smiles` is the one most likely to be shadowed --
           # every layer above deals in SMILES -- and `unpack` is `unpach` under chython 2's name, so
           # it must resolve to the same module and not to a wrapper somebody added beside it.
           ('smiles', 'chython.core._facade'),
           ('pach', 'chython.core._facade'),
           ('unpach', 'chython.core._facade'),
           ('unpack', 'chython.core._facade'),
           ('mol', 'chython.formats.ctfile._facade'),
           ('rxn', 'chython.formats.ctfile._facade'),
           ('SDFRead', 'chython.formats.ctfile._stream'),
           ('RDFRead', 'chython.formats.ctfile._rdf'),
           ('SGroup', 'chython.formats.ctfile._sgroup'),
           ('SGroupStore', 'chython.formats.ctfile._sgroup'),
           ('add_data_sgroup', 'chython.formats.ctfile._sgroup'),
           ('data_sgroups', 'chython.formats.ctfile._sgroup'),
           ('xyz', 'chython.formats.xyz'),
           ('XYZFrame', 'chython.formats.xyz'),
           ('XYZAtom', 'chython.formats.xyz'),
           ('mol2', 'chython.formats.mol2'),
           ('read_mol2', 'chython.formats.mol2'),
           ('mol2_mol', 'chython.formats.mol2'),
           ('Mol2ParseError', 'chython.formats.mol2'),
           ('build_molecule', 'chython.formats.pdb._builder'),
           ('pdb', 'chython.formats.pdb._legacy'),
           ('read_pdb', 'chython.formats.pdb._legacy'),
           ('mmcif', 'chython.formats.pdb._mmcif'),
           ('read_mmcif', 'chython.formats.pdb._mmcif'),
           ('PDBRecord', 'chython.formats.pdb._records'),
           ('mrv', 'chython.formats.xml._facade'),
           ('cml', 'chython.formats.xml._facade'),
           ('read_cml', 'chython.formats.xml._cml'),
           ('write_cml', 'chython.formats.xml._cml'),
           ('read_mrv', 'chython.formats.xml._mrv'),
           ('write_mrv', 'chython.formats.xml._mrv'),
           ('read_xml', 'chython.formats.xml._dialect'),
           ('XmlError', 'chython.formats.xml._errors'),
           ('MalformedXml', 'chython.formats.xml._errors'),
           ('UnsupportedXml', 'chython.formats.xml._errors'),
           ('ForbiddenXml', 'chython.formats.xml._errors'),
           ('saturate', 'chython.chemistry._saturate'),
           ('perceive_bonds', 'chython.chemistry._perceive'),
           # the corpus accessors.  Every other way to reach the reaction corpus goes through a
           # molecule and answers about that molecule; these four answer what the corpus HOLDS.
           ('functional_rules', 'chython.reactions._tables'),
           ('protective_rules', 'chython.reactions._tables'),
           ('reaction_rules', 'chython.reactions._tables'),
           ('roles', 'chython.reactions._tables')]


@mark.parametrize('name, owner', _OWNERS)
def test_the_facade_serves_each_name_from_the_package_that_owns_it(name, owner):
    import chython

    attribute = getattr(chython, name)
    assert attribute.__module__ == owner, f'chython.{name} now comes from {attribute.__module__}'


def test_core_and_formats_export_no_common_name():
    """The collision itself, not just its current outcome.

    A name in both packages is served from `formats` with no diagnostic, so the useful assertion is
    that the overlap is empty rather than that today's winner is the one we expected.
    """
    from chython import core, formats

    common = set(core.__all__) & set(formats.__all__)
    assert not common, f'{sorted(common)} exists in both core and formats; the façade serves formats'


def test_every_formats_export_reaches_the_facade():
    """A reader that landed in `formats` but never got wired into `chython`.

    `formats.__all__` is the package's own statement of its entry points, and `chython/__init__.py`
    re-exports it wholesale -- so the two can only disagree if the star-import stops being wholesale or
    something above shadows a name.  Both are silent.  Identity, not presence: a name served from
    somewhere else is the shadowing case and looks fine to `hasattr`.
    """
    import chython
    from chython import formats

    for name in formats.__all__:
        assert hasattr(chython, name), f'chython.formats exports {name} but the façade does not serve it'
        assert getattr(chython, name) is getattr(formats, name), \
            f'chython.{name} is not chython.formats.{name}; something above formats shadows it'


def test_the_pdb_subpackage_still_means_the_subpackage():
    """The collision `formats.__all__` is written to avoid, asserted rather than described.

    `pdb` names both the legacy reader's string form and the subpackage holding it, so importing the
    function into `chython.formats` would leave `chython.formats.pdb` meaning the function.  Six of that
    package's names are the STAR/CIF tokeniser and are reachable only through the subpackage, so the
    subpackage has to keep winning there.  The façade is the other half: `chython.pdb` is the function,
    because `chython` has no `pdb` submodule for it to collide with.
    """
    from types import ModuleType

    import chython
    from chython import formats
    from chython.formats.pdb import parse_star

    assert isinstance(formats.pdb, ModuleType), \
        'chython.formats.pdb is no longer the subpackage; the STAR tokeniser is now unreachable by path'
    assert formats.pdb.parse_star is parse_star
    assert 'pdb' not in formats.__all__, \
        "'pdb' in formats.__all__ puts the function on chython.formats and shadows the subpackage"

    assert callable(chython.pdb) and not isinstance(chython.pdb, ModuleType)
    for name in ('build_molecule', 'pdb', 'read_pdb', 'mmcif', 'read_mmcif', 'PDBRecord', 'PDBAtom',
                 'PDBBond'):
        assert getattr(chython, name) is getattr(formats.pdb, name), \
            f'chython.{name} is not chython.formats.pdb.{name}'


def test_the_facade_has_no_name_that_shadows_a_stdlib_top_level_module():
    """`chython.xml` must not exist, and the reason is the standard library rather than the subpackage.

    `chython.formats.xml` is reached by its own path, so exporting the name would put `chython.xml` beside
    `xml.etree` in every reader's head permanently -- the one collision this facade cannot annotate its way
    out of.  Asserted over the whole surface rather than for `xml` alone, because the next XML-family
    dialect (`json`? `csv`?) is where the rule gets forgotten.
    """
    from sys import stdlib_module_names

    import chython
    from chython import formats

    assert 'xml' not in formats.__all__, "'xml' in formats.__all__ puts chython.xml beside xml.etree"

    shadowed = {name for name in formats.__all__ if name in stdlib_module_names}
    assert not shadowed, f'{sorted(shadowed)} on the façade shadow standard-library module names'
    assert not hasattr(chython, 'xml')


def test_the_whole_corpus_is_reachable_without_a_molecule():
    """`mol.functional_groups()` answers about one molecule; these answer what the corpus holds.

    Every other door into the two corpora is a container method, so the full inventory was reachable
    only through `chython.reactions._tables`, a private module.  Keyed by the name a caller selects by,
    which is what `deprotect(protective=...)` and `react(reaction=...)` take.
    """
    import chython

    functional, protective = chython.functional_rules(), chython.protective_rules()
    assert len(functional) == 249 and len(protective) == 103
    assert all(isinstance(k, str) and v.name == k for k, v in functional.items())
    assert all(isinstance(k, str) and v.name == k for k, v in protective.items())
    assert functional['carboxylic_acid'].smarts and protective['amine_boc'].protects == ('amine',)

    # the other two are name-to-rows, because a reaction name names a row family.
    reactions, handles = chython.reaction_rules(), chython.roles()
    assert len(reactions) == 72 and sum(len(v) for v in reactions.values()) == 294
    assert len(handles) == 53 and sum(len(v) for v in handles.values()) == 87
    assert all(r.name == name for name, rows in reactions.items() for r in rows)


def test_each_corpus_accessor_returns_its_cached_object():
    """A caller may hold the inventory; it is not rebuilt per call and must not be mutated."""
    import chython

    for accessor in (chython.functional_rules, chython.protective_rules, chython.reaction_rules,
                     chython.roles):
        assert accessor() is accessor()
