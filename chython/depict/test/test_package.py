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
"""`chython.depict` must stand on its own.

Drawing code sits at the far end of the one-way dependency order: it may be imported by the containers and
must not import them back, nor reach up to the facade.  `import chython.depict` always executes
`chython/__init__.py` first, so the property is checked the only way it can be -- with a stub parent
standing in for the facade, which nothing reaching for a facade attribute could be satisfied by.
"""
import ast
from pathlib import Path
from subprocess import run
from sys import executable


# What `chython/__init__.py` pulls in, and what `depict/` is not allowed to need.
_facade_payload = ('chython.containers', 'chython.algorithms', 'chython.files', 'chython.formats',
                   'chython.reactor', 'chython.interop', 'chython.periodictable')

_probe = '''
import importlib, sys, types

stub = types.ModuleType('chython')
stub.__path__ = [%r]
stub.__package__ = 'chython'
sys.modules['chython'] = stub   # the facade's __init__ is never executed

importlib.import_module('chython.depict')
importlib.import_module('chython.depict.scene')
importlib.import_module('chython.depict.style')
importlib.import_module('chython.depict.metrics')
importlib.import_module('chython.depict.label')
importlib.import_module('chython.depict.bonds')
importlib.import_module('chython.depict.render.svg')

print(' '.join(sorted(m for m in sys.modules if m.startswith('chython.'))))
print(' '.join(sorted(vars(stub))))
'''


def _import_under_stub_parent():
    chython_dir = Path(__file__).resolve().parents[2]   # .../chython
    result = run([executable, '-c', _probe % str(chython_dir)], capture_output=True, text=True,
                 cwd=chython_dir.parent)
    assert result.returncode == 0, f'importing chython.depict without the facade failed:\n{result.stderr}'
    modules, touched = result.stdout.splitlines()[:2]
    return set(modules.split()), set(touched.split())


def test_importable_without_the_facade():
    """the drawing code loads with nothing but a stub for its own parent package"""
    modules, _ = _import_under_stub_parent()
    assert 'chython.depict' in modules


def test_does_not_import_the_facade_payload():
    modules, _ = _import_under_stub_parent()
    for name in _facade_payload:
        assert not any(m == name or m.startswith(f'{name}.') for m in modules), \
            f'chython.depict imported {name}; depiction must not depend on it'


def test_reads_nothing_off_the_facade():
    """no `from chython import ...` anywhere under `depict/`

    An upward import must resolve a name on the stub parent, where the only names are dunders and one
    attribute per submodule imported; anything else could only have come from the real facade.
    """
    modules, touched = _import_under_stub_parent()
    expected = {m.split('.')[1] for m in modules}   # submodule bindings the import system adds
    assert not {n for n in touched if not n.startswith('__')} - expected, \
        f'depict/ read {touched - expected} off the facade'


#: The packages `depict` may not name AT ALL.  Both sit BESIDE it, not below it: `chemistry`, `formats`,
#: `depict` and `interop` are peers above one core, so an import between any two of them is sideways and
#: `core` is the one home below every consumer.  THERE ARE NO EXCEPTIONS -- a ratchet with an exception is
#: a ratchet with a hole in the shape of the last bug; an entry here means something did not get
#: repointed.  `chython.interop` is deliberately absent, being the one peer `depict` is supposed to talk
#: to (a `clean2d` engine IS a third-party toolkit); `test_interop_is_named_only_inside_a_function` holds
#: the different property that applies there.
_SIDEWAYS = ('chython.chemistry', 'chython.formats')


def _depict_modules():
    """Every non-test module under `depict/`, DISCOVERED and not enumerated.

    `_probe` above names seven by hand, and a hand-written list cannot cover a module nobody added to it.
    """
    root = Path(__file__).resolve().parent.parent          # .../chython/depict
    return sorted(p for p in root.rglob('*.py') if 'test' not in p.parts)


def _named_packages(path, *, module_scope_only=False):
    """The `chython.*` packages a source file imports, absolute and relative alike.

    Both spellings, being one edge written two ways; a relative level is resolved against the file's own
    position.  Import STATEMENTS, not text, so a package name inside a string is not a hit.
    `module_scope_only` stops at the first function or class, separating a load-time dependency from one
    taken to answer a call.
    """
    root = Path(__file__).resolve().parents[3]             # the directory holding `chython/`
    parts = path.relative_to(root).with_suffix('').parts   # ('chython', 'depict', ...)
    package = parts[:-1] if path.name != '__init__.py' else parts
    found = set()

    def visit(node):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package[:len(package) - node.level + 1]
                found.add('.'.join((*base, node.module)) if node.module else '.'.join(base))
            elif node.module:
                found.add(node.module)
        elif module_scope_only and isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                                                     ast.ClassDef)):
            return                                        # a deferred import is not a load-time one
        for child in ast.iter_child_nodes(node):
            visit(child)

    # `encoding='utf-8'`: `depict/` is UTF-8 -- `field.py` documents its solver with `∇` -- and
    # `read_text` without it asks the locale, cp1252 on the Windows runner, which cannot decode that.
    visit(ast.parse(path.read_text(encoding='utf-8')))
    return {name for name in found if name.startswith('chython.')}


def test_depict_names_no_package_beside_it():
    """the layering, read off the source rather than off one import run

    `test_does_not_import_the_facade_payload` sees only what an import EXECUTES, so a dependency behind a
    deferred import is invisible to it until something calls that function.
    """
    offenders = {}
    for path in _depict_modules():
        for name in _named_packages(path):
            if any(name == p or name.startswith(f'{p}.') for p in _SIDEWAYS):
                offenders.setdefault(str(path.name), set()).add(name)
    assert not offenders, f'depict/ names a package beside it: {offenders}'


def test_interop_is_named_only_inside_a_function():
    """the one peer `depict` may talk to, and only when a caller asks for it

    `layout/molecule.py` puts every `to_*` import inside its own engine branch, for two reasons: a
    toolkit nobody named must not be imported to lay a molecule out, and `chython.depict` must load
    without `chython.interop`.  This fails if one of those five imports is hoisted for tidiness.
    """
    hoisted = {path.name: named for path in _depict_modules()
               if (named := {n for n in _named_packages(path, module_scope_only=True)
                             if n == 'chython.interop' or n.startswith('chython.interop.')})}
    assert not hoisted, f'interop named at module scope, so importing depict pulls a toolkit: {hoisted}'
    # and the premise: it IS named somewhere, or the test above is vacuous
    deferred = {n for path in _depict_modules() for n in _named_packages(path)
                if n.startswith('chython.interop')}
    assert deferred, 'no depict module names interop at all; this test has lost its subject'


def test_the_scan_covers_the_modules_it_claims_to():
    """a scan that silently found nothing would pass the test above

    The failure mode is an empty loop: a broken `rglob`, a renamed directory, a `parts` slice off by one.
    ASSERTED AS RELATIVE PATHS, NOT BARE FILENAMES, because `render/svg.py` and a deleted top-level
    `svg.py` answer to the same bare name, and a ratchet can stay green while losing its subject.
    """
    found = {p.relative_to(Path(__file__).resolve().parent.parent).as_posix()
             for p in _depict_modules()}
    assert {'render/svg.py', 'wedge.py', 'bonds.py', '__init__.py'} <= found, found
    assert 'svg.py' not in found, 'the V2 renderer is back at the top level of depict/'
    assert not any('test' in p.parts for p in _depict_modules()), 'test files are not production code'
    # and the resolver really does see relative imports, which is the spelling `depict/` uses
    wedge = next(p for p in _depict_modules() if p.name == 'wedge.py')
    assert 'chython.core.wedge' in _named_packages(wedge), _named_packages(wedge)


def test_settings_have_one_home():
    """`chython.clean2d_engine` is a view of `chython.depict._config`, not a second copy"""
    import chython
    from chython.depict import _config, get_clean2d_engine

    before = chython.clean2d_engine
    try:
        chython.clean2d_engine = 'rdkit'
        assert _config.clean2d_engine == 'rdkit', 'the facade kept its own copy of the setting'
        assert get_clean2d_engine() == 'rdkit'

        _config.set_clean2d_engine('smilesdrawer')
        assert chython.clean2d_engine == 'smilesdrawer', 'the facade did not see the change'
    finally:
        chython.clean2d_engine = before


def test_engine_is_validated_at_assignment():
    from pytest import raises

    import chython

    before = chython.clean2d_engine
    try:
        with raises(ValueError):
            chython.clean2d_engine = 'no-such-engine'
        assert chython.clean2d_engine == before, 'a rejected engine was stored anyway'
    finally:
        chython.clean2d_engine = before


def test_the_v2_renderer_is_gone():
    """a ratchet, in the same shape as `chython/test/test_v2_boundary.py`

    These modules were replaced, not wrapped; restoring one to read a constant out of it would put a
    second, untested renderer in the tree.
    """
    from importlib import import_module

    for name in ('grid', 'retro', 'svg', 'vector'):
        try:
            import_module(f'chython.depict.{name}')
        except ImportError:
            continue
        raise AssertionError(f'chython.depict.{name} is importable again; read it out of git instead')


def test_the_global_settings_dict_is_gone():
    """`depict_settings` mutated one process-wide dict; `DepictStyle` is immutable and passed in"""
    from chython.depict import _config

    assert not hasattr(_config, '_render_config')
    assert not hasattr(_config, 'depict_settings')
    assert len(_config.cpk) == 118  # kept: a palette is not a settings dict


def test_the_facade_no_longer_exports_the_dropped_names():
    import chython

    for name in ('GridDepict', 'RetroDepict', 'grid_depict', 'retro_depict', 'depict_settings'):
        assert not hasattr(chython, name), f'chython.{name} survived the depict rewrite'


def test_the_facade_still_exports_what_replaced_them():
    """`DepictStyle` and the two accessors replace `depict_settings`, which was a `chython.depict`
    export, so they are re-exported from the same address while `chython.depict.style` stays the one
    home."""
    import chython

    assert chython.get_clean2d_engine() is not None
    assert chython.set_clean2d_engine
    assert chython.Clean2DEngine
    from chython.depict import DepictStyle, get_depict_style, set_depict_style     # noqa: F401


def test_x3dom_is_live_and_reads_its_own_parameter_set():
    """The 3D side is two module functions registered onto the container, not a mixin: a `cdef class`
    cannot be extended from outside, so there is nothing for one to attach to.  Its parameters are its
    own frozen dict and never `DepictStyle`, which is 2D throughout -- see `chython/depict/x3dom.py`."""
    from chython.depict import molecule_depict3d, molecule_view3d, x3dom

    assert x3dom.molecule_depict3d is molecule_depict3d
    assert x3dom.molecule_view3d is molecule_view3d
    assert 'X3domMolecule' not in dir(x3dom)
    assert '_render_config' not in x3dom.__dict__
