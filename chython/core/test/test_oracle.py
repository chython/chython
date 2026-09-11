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
"""The oracle's own guards, including the NEGATIVE CONTROL that makes the isolation guard checkable.

WHY THIS FILE EXISTS.  `oracle.py` is test infrastructure, and test infrastructure that is wrong does
not fail -- it passes.  A missing `-I` makes the oracle interpreter import `./chython/`, so every
differential in the repository compares the tree under test to itself, agrees on everything, and
reports success.  It survives review and it survives a green test run, because a worktree has no
`chython.egg-info` for the version pin to trip over.

So the guards get tested like anything else, and the important one is tested by BREAKING IT: drop
`-I` and the tree-under-test guard must fire.  A guard nobody has watched fail is a comment.

Every test here skips when the oracle is not provisioned -- it is an optional test dependency, and a
witness that has not been installed is not a regression.
"""
import ast
from pathlib import Path

from pytest import raises, skip

from . import oracle


def _provisioned():
    if oracle.interpreter() is None:
        skip('the chython 2 oracle is not provisioned; see chython.core.test.oracle.__doc__')


# --- the guards, working -------------------------------------------------------------------------

def test_the_oracle_is_the_pinned_version():
    """A version bump has to be LOUD, not silently absorbed.

    If this fails, re-run the differentials deliberately and then move `VERSION` -- do not assume a
    newer chython 2 is a better oracle, because some of what it changed may be the very defects the
    callers have written down as expected divergences.
    """
    _provisioned()
    version, _ = oracle.probe()
    assert version == oracle.VERSION


def test_the_oracle_is_not_the_tree_under_test():
    """The differential must not be the core against itself. Stated directly, not by inference."""
    _provisioned()
    _, source = oracle.probe()
    assert not source.startswith(str(oracle.ROOT)), source
    assert Path(source).is_file()


def test_verify_passes_on_a_provisioned_oracle():
    """Both guards together, through the entry point every caller actually uses."""
    _provisioned()
    oracle.verify()


# --- the negative control ------------------------------------------------------------------------

def test_WITHOUT_ISOLATION_THE_ORACLE_IMPORTS_THIS_TREE(monkeypatch):
    """THE POINT OF THE FILE.  Drop `-I`, and the guard that matters must fail.

    Run from the repository root -- which is where pytest runs, and the reason the defect was live --
    `python -c` puts the working directory at the head of `sys.path`, so the child imports `./chython/`
    instead of its own site-packages.  `cwd` is forced here rather than assumed, so the control holds
    when the suite is invoked from somewhere else and does not quietly become a no-op.

    If this test ever passes trivially -- that is, if the unisolated child stops importing this tree --
    then either `run` has grown a second layer of isolation (fine, say so here) or the repository has
    stopped being importable from its own root (not fine).  Do not delete it: without it, a missing
    `-I` is invisible.
    """
    _provisioned()
    monkeypatch.chdir(oracle.ROOT)
    version, source = oracle.probe(isolation=None)
    assert source.startswith(str(oracle.ROOT)), (
        'the unisolated oracle did NOT import this tree, so this control is no longer measuring '
        f'anything: it imported {source}')
    # and the version pin catches the same leak, but only where an egg-info exists -- so it is the
    # weaker of the two guards and guard 3 stands next to it rather than instead of it
    assert version != oracle.VERSION or not (oracle.ROOT / 'chython.egg-info').is_dir()


def test_verify_REFUSES_a_leaked_oracle_and_names_the_cause(monkeypatch):
    """`verify` on the leak, driven through the real function rather than a copy of its assertion.

    The probe result is substituted for what an unisolated child actually returns -- measured by the
    control above -- so this exercises `verify` itself.  The message must name `-I`, because the
    symptom of the defect is "every differential agrees", which reads like success.

    The first case is also what PINS THE ORDER of the two guards.  A leak trips both -- the child
    reads `./chython.egg-info` and calls itself 3.0 -- and with the version pin checked first the
    reader is told the symptom instead of the cause.  Written the wrong way round, this assertion is
    what caught it.
    """
    monkeypatch.setattr(oracle, '_PROBE', ('3.0', str(oracle.ROOT / 'chython' / '__init__.py')))
    with raises(AssertionError, match='TREE UNDER TEST'):
        oracle.verify()
    monkeypatch.setattr(oracle, '_PROBE', ('2.23', '/somewhere/site-packages/chython/__init__.py'))
    with raises(AssertionError, match='pinned to'):
        oracle.verify()


# --- absence is a skip, never a failure ----------------------------------------------------------

def test_an_absent_interpreter_skips_and_says_how_to_provision_one(monkeypatch):
    """The oracle is OPTIONAL.  A machine without it runs the whole suite and every direct assertion.

    Asserted through `require` rather than by trusting the decorators, and the message is checked for
    the two things somebody who hit it needs: the environment variable and the pinned version.
    """
    import pytest

    monkeypatch.setattr(oracle, '_RESOLVED', None)
    with raises(pytest.skip.Exception) as caught:
        oracle.require()
    message = str(caught.value)
    assert oracle.ENV_VAR in message
    assert oracle.VERSION in message


def test_the_environment_variable_wins_over_the_cache_path(monkeypatch, tmp_path):
    """So the oracle is not a fact about one developer's home directory."""
    fake = tmp_path / 'python'
    fake.write_text('')
    monkeypatch.setenv(oracle.ENV_VAR, str(fake))
    monkeypatch.setattr(oracle, '_RESOLVED', ...)
    assert oracle.interpreter() == fake


def test_a_path_that_is_not_a_file_resolves_to_absent(monkeypatch, tmp_path):
    """A directory, or a stale path, must read as "not provisioned" and not as a broken invocation."""
    monkeypatch.setenv(oracle.ENV_VAR, str(tmp_path))
    monkeypatch.setattr(oracle, '_RESOLVED', ...)
    assert oracle.interpreter() is None


# --- and no second copy --------------------------------------------------------------------------

#: the two spellings a file resolves the oracle with: the environment variable, and the default path
NAMES = (oracle.ENV_VAR, 'chython2-oracle')


def spawner_literals(path):
    """The oracle's name where it is CODE -- string literals and attribute access, never prose.

    A grep is wrong here for the reason the sibling import guard writes down at length: the name has
    legitimate non-code uses, and a DOCSTRING telling the reader which variable to set to provision the
    oracle trips one.  A guard that punishes accurate documentation gets exceptions bolted onto it until
    it means nothing, so this reads the AST instead: comments never enter one at all, module and
    function docstrings are skipped explicitly, and what is left is the shape a second spawner actually
    takes -- `environ.get('CHYTHON2_ORACLE')`, or the cache path written out as a literal.
    """
    tree = ast.parse(path.read_text(errors='replace'), filename=str(path))
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            first = node.body[0] if node.body else None
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) \
                    and isinstance(first.value.value, str):
                docstrings.add(id(first.value))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and id(node) not in docstrings:
            found.extend(name for name in NAMES if name in node.value)
    return found


def test_THIS_MODULE_IS_THE_ONLY_PLACE_THE_ORACLE_IS_SPAWNED():
    """ONE SPAWNER, and this is what stops it being re-scattered.

    Consolidation is not a property that stays true on its own -- the next differential will need
    chython 2, and writing three lines of `subprocess.run` is easier than finding this module.  So the
    interpreter path and the environment variable are asserted to appear, as code, in exactly one
    file, which is the cheapest signal that a fifth copy has appeared.

    If this fails, the fix is not to widen it: it is to route the new caller through `oracle.ask`,
    `oracle.ask_text` or `oracle.Session`, which is what it needed anyway.
    """
    # `oracle.py` defines them and this file tests the definition -- the two are the surface, not
    # copies of it.  Everything else must go through the module.
    allowed = {Path(oracle.__file__).resolve(), Path(__file__).resolve()}
    root = Path(oracle.__file__).resolve().parents[3]
    offenders = []
    for path in sorted(root.joinpath('chython').rglob('*.py')):
        if path.resolve() in allowed:
            continue
        if spawner_literals(path):
            offenders.append(str(path.relative_to(root)))
    assert offenders == [], (
        'these files name the oracle interpreter themselves instead of going through '
        f'chython.core.test.oracle: {offenders}')


def test_the_fifth_copy_scan_catches_code_and_ignores_prose(tmp_path):
    """The control.  Both halves, because either half alone is a scan that always agrees.

    Written over synthetic sources rather than over a file in the tree: a control pinned to a real
    file goes stale the day somebody fixes that file, which is the best possible reason for a test to
    start failing and the worst possible way to hear about it.
    """
    planted = tmp_path / 'planted.py'
    planted.write_text('from os import environ\n'
                       'from pathlib import Path\n'
                       'py = Path(environ.get("CHYTHON2_ORACLE", "~/.cache/chython2-oracle/bin/python"))\n')
    assert sorted(set(spawner_literals(planted))) == ['CHYTHON2_ORACLE', 'chython2-oracle']

    # prose: a module docstring, a function docstring and a comment, all naming it legitimately
    prose = tmp_path / 'prose.py'
    prose.write_text('"""Set CHYTHON2_ORACLE to point at chython 2."""\n'
                     'from chython.core.test.oracle import ask\n'
                     '# provisioned under ~/.cache/chython2-oracle by default\n'
                     'def f():\n'
                     '    """Skips unless CHYTHON2_ORACLE names an interpreter."""\n'
                     '    return ask("_emit(1)")\n')
    assert spawner_literals(prose) == []

    # and a class docstring, since the walk lists ClassDef separately from the two function kinds
    in_class = tmp_path / 'in_class.py'
    in_class.write_text('class C:\n    """CHYTHON2_ORACLE selects the interpreter."""\n')
    assert spawner_literals(in_class) == []
