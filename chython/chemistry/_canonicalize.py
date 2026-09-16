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
"""`canonicalize()` -- the pre-pass that makes `canonical_bytes` a compound identity rather than a
drawing identity, so equal compounds hash equal and a corpus deduplicates by hash.  The stage order
is a correctness constraint, not taste: `kekule()`, `validate_stereo()`, `standardize()`,
`implicify_hydrogens()`, `neutralize()`, `standardize_kekule()`, `thiele()`, `standardize_isomers()`, and
`kekule()` again only under `keep_kekule=True`.  Each ordering pair is justified at its numbered step
below.

The order alone is not enough, because the last stage can unblock the third one: a `tautomer` row
wanting a free ring nitrogen cannot fire while the mobile hydrogen sits on it, and the placement is
what moves that hydrogen off.  So steps 3 to 8 run to a fixed point rather than once -- step 9.
"""
from ._hydrogens import implicify_hydrogens
from ._isomers import standardize_isomers
from ._kekule_form import standardize_kekule
from ._protomers import neutralize
from ._standardize import standardize
from ..core import LOST, LogRecord, MoleculeContainer, recording


__all__ = ['canonicalize']


_RULE_ROUNDS = 'canonicalize:rounds'

#: How many times steps 3 to 8 may be re-run before the pipeline gives up and says so.  Every shape
#: measured converges in two, and the loop cannot cycle in principle: the `tautomer` rows move a
#: hydrogen from oxygen or sulfur to nitrogen and never back, and the placement stage never makes an
#: oxygen or a sulfur a site.  The cap is here so a rule table that breaks either half is reported as
#: a loss rather than hanging.
_ROUNDS_MAX = 5


def canonicalize(molecule: MoleculeContainer, *, fix_tautomers: bool = True,
                 keep_kekule: bool = False) -> bool:
    """Bring `molecule` to the representation two drawings of one compound share.  Did it change?

    Run before deduplicating by `canonical_bytes`, `__hash__` or `__eq__`, which otherwise answer "same
    drawing" rather than "same compound".

    `molecule.log` gets one record per thing done or declined, tagged with the stage that wrote it, so
    `mol.log.by_stage('standardize')`, `mol.log.repaired()` and `mol.log.lost()` all answer afterwards.
    `fix_tautomers` is forwarded to `standardize()`, and switching it off gives up part of the guarantee
    (`Oc1ccccn1` and `O=c1cccc[nH]1` stop hashing equal); it deliberately does not reach
    `standardize_isomers()`, which picks between two valid annular forms rather than repairing one.
    `keep_kekule=True` costs a second kekulisation rather than a skipped `thiele()`, the placement stage
    needing the aromatic form.

    CHARGES ARE PAIRED OFF, not preserved atom by atom: glycine's zwitterion and its neutral drawing
    share a key, because step 5 runs `neutralize()`.  `neutralize()` leaves the NET charge untouched, so
    sodium acetate stays sodium acetate -- there is no proton in it to move -- while ammonium acetate
    becomes acetic acid and ammonia, both drawings of one salt.

    One stage does move the net charge, `standardize()`'s organometallic completion: a zinc or magnesium
    holding one carbon and no halide is charged rather than left as a neutral one-coordinate metal.  It
    is the only place in the pipeline where the total changes, and it says so in the log.
    """
    # A refusal at the answer boundary, which is the only place one belongs.  `smiles()` returns
    # whichever container its string describes, so a `>>` in a structure column arrives here as a
    # reaction; without this the first read of `canonical_bytes` fails with an `AttributeError` naming
    # a private attribute, which tells a caller nothing about what the pass takes.
    if not isinstance(molecule, MoleculeContainer):
        raise TypeError(f'canonicalize() takes a MoleculeContainer, got '
                        f'{type(molecule).__name__}; a ReactionContainer has its own canonicalize(), '
                        f'which runs this pass on each of its molecules')

    # the bool is measured, not accumulated: steps 1 and 6 are a round trip, so summing the stages'
    # own flags would report a change through both ends of a no-op on an already-canonical molecule.
    before = molecule.canonical_bytes

    # 1. Kekule first, so the group rules see definite bond orders.  An aromatic system with no Kekule
    #    form does not stop the pipeline: it is reported and the rest still runs.  `kekule()` also
    #    heals the hydrogen counts its own orders made derivable, so nothing here does that -- the
    #    pipeline is literally the manual steps.
    #    Nothing is recorded here about `result.unresolved`: the kekuliser writes a `LOST` record per
    #    system to `molecule.log` itself, naming it in `atoms`, so a second one here would summarise an
    #    event already reported.  `check_valence()` is still how those atoms are found.
    molecule.kekule()

    # 2. Drop the parities the constitution does not justify, on the localised form and before step 6
    #    hides a bond order behind an aromatic one.  A cis/trans sign on a double bond an alternating
    #    cycle can move states which Kekule form was drawn, not the compound's geometry, and it enters
    #    `canonical_bytes` whether or not the unit is stereogenic -- so the two bond-shift drawings of
    #    1,2-dimethylcyclooctatetraene keep two keys until it is gone.  The report is not recorded: the
    #    parities it names are the reader's own input and `validate_stereo()` returns them to a caller
    #    that wants them.
    molecule.validate_stereo()

    # 3. Repair the drawing.
    standardize(molecule, fix_tautomers=fix_tautomers)

    # 4. Explicit hydrogens are a `canonical_bytes` difference, so they have to go.
    implicify_hydrogens(molecule)

    # 5. Pair off the charges an acid/base row can pair off, so a zwitterion and its neutral drawing
    #    hash equal.  AFTER step 4, because `acids.tsv` reads implicit hydrogens: a cation drawn with
    #    hydrogen ATOMS is invisible to it until they have been folded in.  `keep_charge` stays at its
    #    default -- the net charge is part of the compound, so a canonical form may move a proton but
    #    never create or destroy one.  A quaternary ammonium keeps its counterion: it has no proton to
    #    give, so the pass finds no donor and declines.
    neutralize(molecule)

    # 6. Which Kekule form, decided before step 7 reads one.  `thiele()` refuses a candidate ring whose
    #    atom holds its double bond outside the ring, so a compound with several Kekule forms has one
    #    aromatic form per form until this stage picks between them; a ring too big for `thiele()` to
    #    consider gets a canonical alternation here instead, nothing later collapsing its two.  AFTER
    #    step 3, whose group rules are written against the orders the drawing carried, and after step 2,
    #    whose parities would otherwise pin the very bond this stage moves.
    standardize_kekule(molecule)

    # 7. Back to the aromatic form, which is the representation callers compare.  Ahead of step 8,
    #    because a mobile hydrogen is a property of the aromatic form: step 1's definite orders already
    #    say where the hydrogen is, leaving the placement stage nothing to choose.
    #    `result.refused` is not recorded on top of the pass's own records either, and for the same
    #    reason -- with the severity the aromatiser itself states, which is `REFUSED` and not a loss.
    molecule.thiele()

    # 8. Canonical placement of mobile hydrogens and charges -- what makes the two N-H forms of
    #    4-methylimidazole hash equal.  Not gated by `fix_tautomers`: that flag withholds local repair
    #    rules, and this picks which of two valid drawings to keep rather than repairing one.
    moved = standardize_isomers(molecule)

    # 9. Steps 3 to 8 again, while the placement keeps unblocking a repair.  `Oc1[nH]cnc2nncc1-2` is
    #    the shape: its mobile hydrogen sits on the one ring nitrogen the hydroxy-azine rows need free,
    #    so step 3 declines, and by the time step 8 has moved it the repair is behind us -- the drawing
    #    kept its hydroxy form and the same compound drawn the other way got the oxo form and a
    #    different key.  Only the placement is re-entered from, since it is the one stage that can put
    #    the molecule back into a shape an earlier stage would have acted on.
    #
    #    `kekule()` leads, and not for the reason step 1 does: the `tautomer` rows are written against
    #    definite bond orders, so on the aromatic form step 7 left behind they match nothing at all and
    #    re-running step 3 would be a guaranteed no-op.  Step 5 is re-entered only behind a repair,
    #    which is the only thing that can hand it a charged site it has not already seen.  Step 6 rides
    #    with the aromatisation and not with the repair: the kekulisation above is free to come back with
    #    a different Kekule form than the one it was handed, which is the form step 7 would then read.
    #    Step 2 is not re-entered: the parities it can justify are a property of the constitution, which
    #    no stage from here on changes.
    for _ in range(_ROUNDS_MAX):
        if not moved:
            break
        molecule.kekule()
        changed = standardize(molecule, fix_tautomers=fix_tautomers)
        if changed:
            implicify_hydrogens(molecule)
            neutralize(molecule)
        standardize_kekule(molecule)
        molecule.thiele()               # unconditional: step 6's form is what a caller compares, and
        if not changed:                 # the kekulisation above has to be undone either way
            break
        moved = standardize_isomers(molecule)
    else:
        with recording(molecule, stage='canonicalize') as log:
            log.append(LogRecord(_RULE_ROUNDS, (), f'repair and placement were still changing the '
                                 f'molecule after {_ROUNDS_MAX} rounds; the pipeline stopped there, so '
                                 f'this molecule is not a fixed point and two drawings of it may not '
                                 f'share a key', LOST))

    # 10. `keep_kekule` undoes step 7 rather than skipping it: skipping 7 would skip 8 with it, and the
    #     flag would then decide which tautomer the caller gets.  Step 6 is not re-run behind it and
    #     would have nothing to do: an aromatic ring is unwound into a ring holding every double bond it
    #     has room for, and which of its alternations comes back is the kekuliser's answer to give.
    if keep_kekule:
        molecule.kekule()

    return molecule.canonical_bytes != before
