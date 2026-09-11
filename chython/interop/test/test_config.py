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
Where the external-tool configuration lives, and that it cannot drift back to the facade.

`chython.interop.config` owns `class_paths` (the CDK/OPSIN JVM classpath) and `conformer_engine`, and
must own them rather than read them off `chython`, which would invert the dependency direction.
"""
from pathlib import Path
from subprocess import run
from sys import executable
from pytest import raises
from chython.interop import config


ROOT = Path(__file__).resolve().parent.parent


def test_defaults():
    """The two knobs exist with the documented defaults."""
    assert config.conformer_engine == 'rdkit'
    assert len(config.class_paths) == 2


def test_class_paths_reads_env(monkeypatch):
    """`CDK_PATH` and `OPSIN_PATH` still select the jars.

    Read at each `class_paths` access rather than once at import, so setting the variable after
    `import chython` works.
    """
    monkeypatch.setenv('CDK_PATH', '/tmp/some-cdk.jar')
    monkeypatch.setenv('OPSIN_PATH', '/tmp/some-opsin.jar')
    monkeypatch.setattr(config, 'class_paths', None)  # None = "ask the environment"
    assert config.class_paths == ['/tmp/some-cdk.jar', '/tmp/some-opsin.jar']


def test_class_paths_explicit_wins_over_env(monkeypatch):
    """An explicit assignment beats the environment; that is what makes it an override."""
    monkeypatch.setenv('CDK_PATH', '/tmp/from-env.jar')
    monkeypatch.setattr(config, 'class_paths', ['/tmp/explicit.jar'])
    assert config.class_paths == ['/tmp/explicit.jar']


def test_conformer_engine_rejects_unknown_value(monkeypatch):
    """A bad engine name fails where it is written, not later inside a converter."""
    with raises(ValueError, match='conformer_engine'):
        config.conformer_engine = 'no-such-engine'
    assert config.conformer_engine == 'rdkit'  # unchanged


def test_conformer_engine_accepts_known_values(monkeypatch):
    for engine in ('rdkit', 'cdpkit'):
        monkeypatch.setattr(config, 'conformer_engine', engine)
        assert config.conformer_engine == engine


def test_facade_alias_reads_and_writes_through():
    """`chython.conformer_engine` reads and writes the real knob, validation included.

    The write half is the point: a copied value would leave `chython.conformer_engine = 'cdpkit'`
    running and doing nothing.
    """
    import chython

    assert chython.conformer_engine == config.conformer_engine
    assert chython.class_paths == config.class_paths
    try:
        chython.conformer_engine = 'cdpkit'
        assert config.conformer_engine == 'cdpkit'
        with raises(ValueError, match='conformer_engine'):
            chython.conformer_engine = 'nonsense'
        assert config.conformer_engine == 'cdpkit'
    finally:
        chython.conformer_engine = 'rdkit'


def test_facade_alias_leaves_other_attributes_alone():
    """The alias must not swallow a genuine `AttributeError` on the facade."""
    import chython

    with raises(AttributeError, match='no attribute'):
        chython.definitely_not_a_chython_name


def test_interop_never_imports_the_facade():
    """No module under `chython/interop/` may import `chython` itself.

    A grep and not an import check, because a facade import deferred inside a function is invisible to
    any `sys.modules` snapshot.  Shipped modules only -- a test is a caller, not a layer.
    """
    offenders = []
    for path in sorted(p for p in ROOT.rglob('*.py') if 'test' not in p.parts):
        for i, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
            s = line.strip()
            if s.startswith('from chython import ') or s == 'import chython' or \
                    s.startswith('import chython '):
                offenders.append(f'{path.relative_to(ROOT)}:{i}: {s}')
    assert not offenders, 'interop must not import the facade:\n' + '\n'.join(offenders)


def test_importing_interop_does_not_load_a_toolkit():
    """Importing `chython.interop` must not start a JVM or load RDKit."""
    code = ('import sys; import chython.interop; '
            "print(','.join(sorted(m for m in sys.modules if '.' not in m "
            "or m.startswith('chython.interop'))))")
    out = run([executable, '-c', code], capture_output=True, text=True, cwd=ROOT.parent.parent)
    assert out.returncode == 0, out.stderr
    loaded = set(out.stdout.strip().split(','))
    assert not loaded & {'rdkit', 'jpype', 'openbabel', 'CDPL', 'indigo', 'openclatura'}
    # Interop's own leaves stay lazy too -- only the dispatch module and `config` are loaded.
    assert 'chython.interop._rdkit' not in loaded
    assert 'chython.interop._cdk' not in loaded
    assert 'chython.interop' in loaded


# The chython 2 packages this test measures against.  They stay named although deleted: restoring one
# out of git is how the regression comes back.
_V2_PACKAGES = ('chython.containers', 'chython.algorithms', 'chython.files', 'chython.reactor')


def test_interop_loads_without_chython_2_being_importable_at_all():
    """No module in this package may import chython 2.

    The parent package is replaced by a bare module carrying only `__path__`, so submodules resolve
    while `chython/__init__.py` never runs and `sys.modules` shows only what *this* package reached
    for.  Every submodule is imported, not just the dispatcher: a top-level count never reaches a leaf.
    """
    code = (
        'import sys, types, importlib\n'
        'pkg = types.ModuleType("chython")\n'
        f'pkg.__path__ = [{str(ROOT.parent)!r}]\n'
        'sys.modules["chython"] = pkg\n'
        'import chython.interop\n'
        'for m in ("_rdkit", "_cdk", "_openbabel", "_cdpkit", "_indigo", "_iupac", "_java", "_records",'
        '          "_stereo", "config", "conformers"):\n'
        '    importlib.import_module("chython.interop." + m)\n'
        'print(",".join(sorted(m for m in sys.modules if m.startswith("chython."))))\n'
    )
    out = run([executable, '-c', code], capture_output=True, text=True, cwd=ROOT.parent.parent)
    assert out.returncode == 0, out.stderr
    loaded = set(out.stdout.strip().split(','))
    offenders = sorted(m for m in loaded if m.startswith(_V2_PACKAGES))
    assert not offenders, 'interop reached into chython 2:\n' + '\n'.join(offenders)


def test_the_predicate_answers_no_without_the_facade():
    """`is_container` is True for a core container and False for everything foreign, facade or not.

    Run in a subprocess with `chython` stubbed out: choosing a direction must never need the facade and
    must never raise on an unfamiliar object.
    """
    code = (
        'import sys, types\n'
        'pkg = types.ModuleType("chython")\n'
        f'pkg.__path__ = [{str(ROOT.parent)!r}]\n'
        'sys.modules["chython"] = pkg\n'
        'from chython.interop import is_container\n'
        'from chython.core import MoleculeContainer\n'
        'assert not is_container("CCO") and not is_container(42) and not is_container(None)\n'
        'assert is_container(MoleculeContainer())\n'
        'print("ok")\n'
    )
    out = run([executable, '-c', code], capture_output=True, text=True, cwd=ROOT.parent.parent)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == 'ok'
