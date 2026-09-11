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
"""One way to ask chython 2 a behavioural question: does a ported rule fire on the same molecules?

The oracle is a separately installed chython 2.24, spawned as a subprocess under `-I`.  A subprocess
is the only sanctioned way to consult it, because an in-process import would resolve to this worktree
and silently compare V3 to itself.  Absent oracle is a skip; point `CHYTHON2_ORACLE` at an
interpreter that has chython 2 to run these.
"""

from ...core.test.oracle import VERSION as ORACLE_VERSION, requires_oracle, run, verify


#: `ORACLE_VERSION` is re-exported rather than restated -- two copies of a pin are one pin and one lie.
__all__ = ['ORACLE_VERSION', 'ask', 'needs_oracle']

#: Decorate any test that shells out.  Absent oracle is a skip, never a failure.
needs_oracle = requires_oracle


def ask(script: str, stdin: str = '') -> list[str]:
    """Run `script` under the pinned chython 2 and hand back its stdout lines, identity lines removed.

    An adapter, not the oracle: a test enforces that `chython/core/test/oracle.py` is the only place
    the interpreter is spawned, so the path, version pin, `-I` flag and identity guards exist once.
    Keep scripts tiny and make them print data -- analysis done on that side is untestable here.
    """
    verify()                        # version pinned, and the child is not this checkout
    # `import sys` and nothing else: every caller writes its answer with `sys.stdout.write`, and the
    # identity lines are `verify`'s job, not the payload's.
    result = run('-c', 'import sys\n' + script, input=stdin)
    assert result.returncode == 0, f'the oracle failed:\n{result.stderr}'
    return result.stdout.splitlines()
