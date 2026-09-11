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
# THE ONE IMPLICIT-HYDROGEN DERIVATION, AND WHY IT IS ITS OWN LAYER
#
# Any molecule parsed from any supported format has its implicit hydrogen counts set AT READ TIME --
# before `kekule()`, before `standardize()`, before `canonicalize()` -- by ONE algorithm, with the same
# answer for the same atom whatever format it arrived in.  Only a genuinely ambiguous atom is left as
# `H_UNKNOWN`.
#
# ONE ALGORITHM AND NOT ONE PER READER, because a second copy drifts and the aromatic half is where it
# drifts first.  `arom_classify_atom` is `cdef` in `_kekule.pxi`, so a derivation written in plain
# Python cannot reach the one function that settles the aromatic case and has nothing left but to
# answer None for every aromatic atom -- benzene carbon, a naphthalene fusion carbon, a pyridine N and
# a thiophene S alike, none of which the valence table's real limitation ("no row admits order 4")
# excuses.  The derivation therefore lives where the pieces already do -- `arom_classify_atom` in
# `_kekule.pxi`, `val_implicit_h` in `_valence.pxi` -- and every caller asks it.  This file is
# included after both.
#
# WHICH VALENCE MODEL, AND WHY SMILES IS NOT A CALLER OF THE WHOLE THING
#
# There are two valence models in this core and `_smiles_read.pxi`'s header forbids merging them.
# `smv_default_h` is the SMILES NOTATION model: eleven elements, charge ignored, environment ignored,
# and it answers "what is a reader REQUIRED to infer when the bracket is absent".  `val_implicit_h`
# is the CHEMISTRY collection: ~1036 rows over element, charge, radical, bond-order sum and, for the
# 0.6% that need it, the neighbour environment.
#
# The universal derivation here uses the CHEMISTRY collection, and it has to.  Build it on the
# notation model instead and every MDL record loses exactly what the baseline demands: a charged
# aromatic (pyridinium N+, an azolium) because the notation model ignores charge, and every element
# outside the organic subset -- Se, Si, Sn, As, Te -- because the notation model has no row for one.
#
# SMILES therefore calls `hyd_arom_takes` and not `hyd_derive_atom`, and that is not a leftover.  A
# bare SMILES atom can only BE one of the eleven, its count is fixed by OpenSMILES rather than by
# chemistry, and a bracket atom states the number outright -- so the notation model is not an
# approximation there, it is the definition.  What SMILES must not duplicate is the valence lookup's
# other half, the aromatic classification: that is the half that must not drift, so it is shared, and
# the two valence models stay apart.
#
# WHAT "AMBIGUOUS" MEANS, EXACTLY, AND WHY IT IS NOT A JUDGEMENT CALL
#
# One class of atom has no locally derivable count: the pnictogen with two aromatic neighbours and no
# hydrogen stated.  With one hydrogen it donates its lone pair and takes no ring double bond
# (pyrrole); without one it contributes a single pi electron and must take one (pyridine).  Both
# readings are valences; which one holds is decided by the RING, which is a global question this file
# deliberately does not ask.
#
# The test for it is not a pattern and not a list of elements -- it is the classifier's own answer.
# `arom_classify_atom` returns AROM_MAY for exactly the atoms whose class the ring decides, and MUST
# or MUST_NOT for every atom it can settle alone; "pyrrole or pyridine, the classic free choice" is
# the comment on the arm.  So the criterion is one comparison, it needs no second copy of the
# chemistry, and an arm that becomes undecidable later is covered the day it is written.  Today three
# arms answer MAY: neutral N/P/As with two aromatic neighbours, neutral P/As with three, and cationic
# N/P/As with two.  Nothing else in the table is ambiguous, which is why nothing else goes unknown.
#
# An ambiguous atom gets `H_UNKNOWN`, which is not a hole but the input to the mechanism that
# resolves it: `arom_prepare` reads an unknown nibble back as `AROM_H_UNSTATED`, so `kekule()`
# receives exactly the freedom it needs and logs the count it derives.  Storing a guessed 0 there
# instead is strictly worse -- it looks like a fact, and it takes the freedom away.
#
# `count_stated` IS THE ONE PLACE FORMATS LEGITIMATELY DIFFER, and it is an argument rather than a
# branch.  A SMILES string always states the count: the language's own rules make a bare lowercase
# `n` a nitrogen with no hydrogen and require `[nH]` to say otherwise, so no SMILES atom is ever
# ambiguous in the sense above and the reader passes True.  A CTfile states nothing of the kind -- a
# bond block with type 4 and no `MRV_IMPLICIT_H` has no channel for it -- so it passes False and the
# pyrrole-versus-pyridine nitrogen comes back unknown, to be settled by `kekule()`.  Same function,
# same table, same answer for every atom where the formats say the same thing.


cdef enum:
    # WHY THE CALLER GETS A CODE AND NOT A SENTENCE.  `chython/formats/` reports damage as a plain
    # `str` on a caller-supplied list, under a ratcheted prefix convention (`test_log_prefix.py`):
    # `unsupported: ` means the file is fine and we are the limitation, no prefix means the file was
    # broken or the gap is our own storage.  The same underlying fact sits on different sides of
    # that line depending on the format, so a message composed HERE would put half of them on the
    # wrong side.  The code says what happened; the caller says whose fault it is.
    HYD_DERIVED = 0                 # answered
    HYD_DERIVED_OTHER_READING = 1   # answered, but not by the reading the classifier picked
    HYD_NO_VALENCE_RULE = 2         # the collection describes nothing here -- our gap
    HYD_AMBIGUOUS_AROMATIC = 3      # pyrrole-versus-pyridine; needs a Kekule form, not a table
    HYD_REASON_MASK = 7

    # A FLAG BIT AND NOT A FIFTH CODE, because it answers a different question.  "This element in
    # this state has no aromatic form" describes the INPUT; the four codes above describe what the
    # derivation managed.  Both can be true at once -- an aromatic-bonded helium still gets read as
    # saturated and may well have a valence -- and the SMILES reader logs them as two sentences,
    # which a single-valued field cannot express.  So: `reason & HYD_REASON_MASK` for the outcome,
    # `reason & HYD_NO_AROMATIC_FORM` for the observation, and no caller has to rank them.
    HYD_NO_AROMATIC_FORM = 8


cdef bint hyd_arom_takes(uint32_t element, int charge, bint radical, uint32_t nbrs, bint exo,
                         bint count_stated, int stated_h, uint8_t *takes,
                         uint8_t *flag) noexcept nogil:
    """Does this aromatic atom take a ring double bond?  False when nothing local can say.

    THE SHARED HALF, and the only half worth sharing: `arom_classify_atom` is `cdef`, so this is what
    a Python caller could not reach and what every second copy of the derivation got wrong.  Both
    valence models call it and neither owns it.

    `flag` comes back with `HYD_NO_AROMATIC_FORM` set when this element in this state has no aromatic
    form at all.  That is not a failure -- the atom is still read as saturated and may well have a
    valence -- so it is a flag on the side rather than a False.

    `stated_h` is `AROM_H_UNSTATED` for a caller that is DERIVING the count and the count itself for
    one that already has it.  The distinction matters because the classifier's ambiguity is a
    consequence of not knowing: given a number it settles the pnictogen outright -- 0 is pyridine, 1
    is pyrrole NH -- so a caller with a count never sees AROM_MAY and never needs the two-reading
    fallback below.  `valence_report` is that caller, which is why an aromatic pyridine nitrogen is a
    checkable atom and not a shrug.

    `count_stated` is the WEAKER form of the same claim, for a caller whose notation fixes the count
    but which has not computed it yet: it withdraws the ambiguity gate without supplying a number.
    See the file header for the one format that may.
    """
    cdef uint8_t invalid = 0
    cdef uint8_t cls
    flag[0] = 0
    cls = arom_classify_atom(element, charge, radical, nbrs, exo, stated_h, &invalid)
    if invalid:
        flag[0] = HYD_NO_AROMATIC_FORM
    # THE AMBIGUITY GATE, AND THE CLASSIFIER STATES IT IN ONE WORD.  Asked with no hydrogen count,
    # `arom_classify_atom` answers AROM_MAY for exactly the atoms whose class the ring decides and a
    # local look cannot -- "pyrrole or pyridine, the classic free choice" is the comment on the arm
    # itself -- and MUST or MUST_NOT for every atom it can settle alone.  So the test is the answer,
    # not a re-interrogation of it: no pattern, no element list, and a new undecidable arm is covered
    # the day it is written.  Today the arms that answer MAY are neutral N/P/As with two neighbours,
    # neutral P/As with three, and cationic N/P/As with two -- pyrrole-versus-pyridine and its charge
    # and heavy-pnictogen analogues, which is precisely the set that has no local answer.
    elif not count_stated and cls == AROM_MAY:
        takes[0] = 1
        return False
    # `AROM_MAY` collapses to "takes one", which is right for the caller that gets here: only a
    # format whose notation fixes the count passes `count_stated`, and both such notations -- a bare
    # aromatic atom in SMILES, a `VAL=`-free aromatic atom the caller has vouched for -- mean the
    # no-hydrogen reading.  The kekuliser makes the same collapse for the same reason.
    takes[0] = 0 if cls == AROM_MUST_NOT else 1
    return True


cdef bint hyd_derive_atom(uint32_t element, int charge, bint radical, uint32_t osum, uint32_t arom,
                          uint32_t nbrs, bint exo, bint count_stated,
                          const uint16_t *env, uint32_t env_len,
                          uint32_t *h_out, uint32_t *sum_out, uint8_t *reason) noexcept nogil:
    """THE derivation, against the chemistry collection.  True with `h_out` written.  `reason` always.

    `osum` is the sum of non-aromatic bond orders, order 8 excluded, explicit hydrogens included.
    `arom` counts the aromatic bonds.  `nbrs` is `arom_classify_atom`'s neighbour count -- every bond
    of any order except 8 -- and is NOT the arena's degree.  `exo` is a double or triple bond outside
    the aromatic set.  `env` holds the VAL_ENV tokens of the NON-aromatic neighbours only, because an
    aromatic neighbour's order is not known until the ring is kekulised and the collection refuses
    order 4 outright; the aromatic bonds reach the collection through `osum` and nowhere else.  An
    atom whose row needs an environment it cannot be given comes back unanswered, not guessed.

    `sum_out` is the bond-order sum the chosen reading charged the atom -- `osum` when nothing is
    aromatic, `osum + arom + takes` when something is.  It exists so a caller can name the number in
    a message without re-deriving `takes`, which would mean a second copy of the decision.

    THE RADICAL IS A PARAMETER HERE AND A SUM TERM IN THE SMILES READER, and the difference is not
    an inconsistency.  `val_implicit_h` takes `radical` as its own argument and the collection's rows
    are written against it, so charging an extra unit of bond order for it as well would count the
    unpaired electron twice -- measured: a nitroxide oxygen, radical with one single bond, answers 0
    hydrogens at sum 1 and no rule at sum 2, so the double count turned `CN([O])C` into an atom whose
    count nobody could derive.  `smv_default_h` has no radical parameter at all, so over there adding
    it to the sum is the only way to express it, and `smi_cx` measured that the two agree.
    """
    cdef uint8_t takes = 0, flag = 0
    cdef int h

    # HYDROGEN CARRIES NO IMPLICIT HYDROGEN, whatever it is drawn bonded to, and that is a property of
    # the element rather than a lookup: one valence electron, so there is never a second one to hold a
    # partner with.  It is stated HERE because the collection is not a statement of it -- it holds rows
    # for hydrogen at a bond-order sum of one and nowhere else, so a bare `H` atom, an `H` on nothing
    # but a dative contact, and a diborane-style bridging `H` all fall off the table and come back
    # undeterminable.
    #
    # IT IS STATED HERE RATHER THAN IN A READER because it was in a reader, and that is exactly how the
    # formats came to disagree: `chython/formats/ctfile/_hydrogens.py` carried this short-circuit and
    # answered 0, while `chython/formats/mol2.py`, which delegates straight to the derivation, answered
    # `None` for the same lone `H` -- measured.  Two readers, one element fact, two answers.  A fact
    # about an element belongs where every reader shares it.
    if element == 1:
        sum_out[0] = osum
        h_out[0] = 0
        reason[0] = HYD_DERIVED
        return True

    # A MARKER CARRIES NO IMPLICIT HYDROGEN, and its count is derived rather than unknown: an
    # attachment point holds a fragment, and there is nothing about it left to determine.  The
    # question the env interning answers is its NEIGHBOUR's count, where it reads as carbon.
    if element == 0:
        sum_out[0] = osum
        h_out[0] = 0
        reason[0] = HYD_DERIVED
        return True

    if not arom:
        sum_out[0] = osum
        h = val_implicit_h(element, charge, radical, osum, env, env_len)
        if h != VAL_NO_RULE:
            h_out[0] = <uint32_t> h
            reason[0] = HYD_DERIVED
            return True
        reason[0] = HYD_NO_VALENCE_RULE
        return False

    if not hyd_arom_takes(element, charge, radical, nbrs, exo, count_stated, AROM_H_UNSTATED,
                          &takes, &flag):
        sum_out[0] = osum + arom + takes
        reason[0] = HYD_AMBIGUOUS_AROMATIC
        return False

    sum_out[0] = osum + arom + takes
    h = val_implicit_h(element, charge, radical, sum_out[0], env, env_len)
    if h != VAL_NO_RULE:
        h_out[0] = <uint32_t> h
        reason[0] = HYD_DERIVED | flag
        return True
    # The classification the ring implies has no valence and the other one does.  Reported, because
    # the stored count then disagrees with the class the kekuliser will pick.
    h = val_implicit_h(element, charge, radical, osum + arom + (1 - takes), env, env_len)
    if h != VAL_NO_RULE:
        h_out[0] = <uint32_t> h
        reason[0] = HYD_DERIVED_OTHER_READING | flag
        return True
    reason[0] = HYD_NO_VALENCE_RULE | flag
    return False


DEF HYD_ENV_MAX = 16    # the arena's degree ceiling is 14; two spare so the bound is never the bug


cdef void hyd_atom_context(atom_t *atoms, uint32_t *ptr, halfedge_t *edges, uint32_t i,
                           uint32_t *osum, uint32_t *arom, uint32_t *nbrs, bint *exo,
                           uint16_t *env, uint32_t *env_len) noexcept nogil:
    """Gather one atom's arguments off the CSR.  The aromatic set is every stored order 4.

    `arom_prepare` asks a narrower question -- is this bond a LIVE edge of the set the caller
    stated -- because a caller may kekulise a subset.  Nobody derives hydrogens against a subset:
    the count is a property of what the arena holds, so here the aromatic set is simply the
    aromatic bonds, and `exo` is a double or triple that is not one of them.

    `env` is written with the VAL_ENV token of each non-aromatic, non-dative neighbour and must have
    room for `HYD_ENV_MAX`.  The two excluded orders are excluded for different reasons and both are
    `_valence.pxi`'s: it REFUSES order 4 (no aromatic rows exist) and it refuses order 8 (whether a
    donated lone pair counts is the caller's policy, and this package's policy, stated in
    `chemistry/_implicit.py`, is that it does not).
    """
    cdef uint32_t j
    cdef uint8_t order, nbr_element
    osum[0] = 0
    arom[0] = 0
    nbrs[0] = 0
    exo[0] = False
    env_len[0] = 0
    for j in range(ptr[i], ptr[i + 1]):
        order = edges[j].order
        if order == 8:
            # A dative contact carries no electron pair, so it is in neither the bond-order sum nor
            # the classifier's neighbour count.  `_valence.pxi` and `arom_classify_atom` agree.
            continue
        nbrs[0] += 1
        if order == 4:
            arom[0] += 1
            continue
        osum[0] += order
        if order >= 2:
            exo[0] = True
        if env_len[0] < HYD_ENV_MAX:
            nbr_element = atoms[edges[j].to].element
            # R (element 0) reads as carbon for a neighbour's environment: element 0 appears in no
            # valence row, so a marker neighbour matches nothing and reports a false violation wherever
            # a state is enumerated by neighbour element with no `env=*` fallback AT THAT KEY.  64
            # elements have such a state; dimethyl sulfone's `S 0 0 6 0 -C -C =O =O` is one.
            if nbr_element == 0:
                nbr_element = 6
            env[env_len[0]] = <uint16_t> ((<uint32_t> order << VAL_ENV_SHIFT) | nbr_element)
            env_len[0] += 1


cdef int hyd_check_sum(uint32_t element, int charge, bint radical, uint32_t osum,
                       const uint16_t *env, uint32_t env_len, int want_h) noexcept nogil:
    """`val_check` for one bond-order sum, plus the case `val_check` cannot express: no count at all.

    `want_h == VAL_ANY_H` asks "does ANY row accept this element in this charge and radical state at
    this sum, whatever the hydrogen count" -- which is the only complete question available about an
    atom holding `H_UNKNOWN`.  `val_check` cannot be handed that: its parameter is a `uint32_t` count
    and the sentinel is 15, which is inside `H_NIBBLE_MAX`, so passing the raw nibble through would
    have the collection look for a row with FIFTEEN hydrogens and report a violation on every atom
    whose count is merely unrecorded.  The Python version this replaces passed 0 instead, which is
    the mirror-image error -- a claim that the atom has no hydrogens, made about an atom whose
    hydrogens nobody derived.

    Both wrong answers came from the same place: a count is not optional in the collection's
    vocabulary, so the caller has to widen the question rather than invent a value for it.  The
    skeleton is still checkable, and that is worth checking -- a hexavalent neutral carbon is a
    violation no hydrogen count could rescue.
    """
    if want_h != VAL_ANY_H:
        return val_check(element, charge, radical, osum, env, env_len, <uint32_t> want_h)
    if val_scan(element, charge, radical, osum, env, env_len, VAL_ANY_H) != VAL_NO_RULE:
        return VAL_VALID
    # the same three-way split `val_check` makes, and for the same reason: "no row matched" and "no
    # row exists" are different claims and only the first is about the molecule
    return VAL_VIOLATION if val_described(element, charge, radical) else VAL_UNKNOWN


cdef int hyd_check_atom(uint32_t element, int charge, bint radical, uint32_t osum, uint32_t arom,
                        uint32_t nbrs, bint exo, const uint16_t *env, uint32_t env_len,
                        int want_h) noexcept nogil:
    """VAL_VALID, VAL_VIOLATION or VAL_UNKNOWN -- the aromatic atom read the way the DERIVATION reads it.

    THE CHECK READS AN AROMATIC ATOM EXACTLY AS THE DERIVATION DOES.  `check_valence` and
    `calc_implicit` ask one question about one atom against one collection, so answering it in two
    places -- one consulting `arom_classify_atom`, one short-circuiting every atom with an aromatic
    bond to `'unknown'` -- makes benzene's carbon 1 hydrogen to one half and "the collection says
    nothing about this atom" to the other.  That no valence row admits order 4 is true and is not the
    whole story: order 4 does not need a row, it needs the classifier's decision about whether the atom
    takes a ring double bond, and then the ordinary rows answer.

    EITHER READING ACQUITS, and the asymmetry is deliberate.  A violation is a claim ABOUT THE
    MOLECULE, so it has to hold under every Kekule form the ring could still take;
    an unkekulised ring means the atom's true order sum is not yet a fact, and the classifier's
    answer is a strong prediction rather than one.  Reporting a violation that a legal Kekule
    assignment would dissolve is inventing bad input.  It also keeps the two halves from
    contradicting each other outright: `hyd_derive_atom` falls back to the second reading when the
    classifier's pick has no row, so a check that ignored that fallback would report a violation on a
    count THIS FILE wrote.

    An atom the classifier cannot settle is `VAL_UNKNOWN` -- not because the collection is silent
    about the element, but because we have not established which state to ask about.  That is the only
    aromatic atom `'unknown'` covers, and not every aromatic atom in the molecule.  IT IS NARROWER
    STILL BECAUSE THE COUNT IS AN ANSWER: a checker,
    unlike a deriver, usually knows the hydrogen count already, and handing it to
    `arom_classify_atom` settles the pnictogen the deriver had to leave open -- 0 is pyridine, 1 is
    pyrrole NH.  So the ambiguous case here is only the atom that is BOTH aromatic-ambiguous and
    missing its count, which after a `kekule()` is nothing at all.
    """
    cdef uint8_t takes = 0, flag = 0
    cdef int verdict, other
    # TRANSLATED RATHER THAN PASSED THROUGH.  Both sentinels happen to be -1 today -- `VAL_ANY_H` in
    # `_valence.pxi`, `AROM_H_UNSTATED` in `_kekule.pxi` -- and relying on that coincidence would make
    # a future change to either one silently classify every countless atom as pyridine.
    cdef int stated_h = AROM_H_UNSTATED if want_h == VAL_ANY_H else want_h

    if not arom:
        return hyd_check_sum(element, charge, radical, osum, env, env_len, want_h)
    if not hyd_arom_takes(element, charge, radical, nbrs, exo, False, stated_h, &takes, &flag):
        return VAL_UNKNOWN

    verdict = hyd_check_sum(element, charge, radical, osum + arom + takes, env, env_len, want_h)
    if verdict == VAL_VALID:
        return verdict
    other = hyd_check_sum(element, charge, radical, osum + arom + (1 - takes), env, env_len, want_h)
    if other == VAL_VALID:
        return VAL_VALID
    # neither reading has a row.  VIOLATION from either one means the collection DESCRIBES this
    # element here and rejected what it was shown, which outranks a silence.
    return VAL_VIOLATION if verdict == VAL_VIOLATION or other == VAL_VIOLATION else VAL_UNKNOWN


def valence_report(MoleculeContainer molecule):
    """`[(n, verdict)]` for every atom the collection does not call `'valid'`.

    The body behind `MoleculeContainer.check_valence()` and `chython.chemistry.check_valence`, which is
    registered rather than compiled in and is one line long.  It lives here, beside the derivation,
    because it is the SAME question -- what does the valence collection make of this atom in this
    environment -- and two copies of that answer drift.  One context walk (`hyd_atom_context`), one
    aromatic classification (`hyd_arom_takes`), two verbs.

    Two verdicts, and they are not the same claim.  `'violation'` means the collection describes this
    element in this charge and radical state and no row accepts what the molecule has -- a statement
    about the molecule.  `'unknown'` means we could not put a complete question: the collection
    describes nothing there, or the atom's aromatic class is the ring's to decide, or its hydrogen
    count was never derived.  Reporting those under one word is how a coverage hole gets mistaken for
    bad input.

    Returns a report.  It never raises and never edits -- refusing a record belongs at the answer
    boundary, not here.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef list out = []
    if not n_atoms:
        return out

    cdef list numbers = molecule._numbers
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, osum = 0, arom = 0, nbrs = 0, env_len = 0
    cdef bint exo = False
    cdef int want_h, verdict
    cdef uint16_t env[HYD_ENV_MAX]      # one scratch buffer for the whole loop, per RULES.md

    for i in range(n_atoms):
        a = atoms + i
        hyd_atom_context(atoms, ptr, edges, i, &osum, &arom, &nbrs, &exo, env, &env_len)
        want_h = VAL_ANY_H if at_implicit_h_unknown(a) else <int> at_implicit_h(a)
        verdict = hyd_check_atom(a.element, a.charge, at_radical(a), osum, arom, nbrs, exo,
                                 env, env_len, want_h)
        if verdict != VAL_VALID:
            out.append((numbers[i], VAL_VERDICT_NAMES[verdict]))
    return out


def derive_implicit_hydrogens(MoleculeContainer molecule, stated=None, bint fill_only=False):
    """Write every derivable implicit hydrogen count on `molecule`.  Returns what it could not.

    THE READ-TIME PASS.  A format reader calls this once, on a molecule it has finished building
    and before it hands it to anybody, and every atom whose count follows from what the file drew
    then holds that count.  It is not a repair and not a mutation of chemistry: nothing here changes
    an element, a charge, a bond or a hydrogen the file STATED.

    `stated` is an iterable of the stable ids whose count the record gave outright -- an MDL
    `MRV_IMPLICIT_H` data S-group, a SMILES bracket, MRV's `hydrogenCount`.  Those atoms are not
    touched and not reported: the file said the number and nothing here second-guesses it.

    `fill_only=True` restricts the pass to the atoms currently holding `H_UNKNOWN` and is the mode
    for AFTER a repair rather than at read time.  Its reason for existing is `kekule()`: the atoms
    this pass has to leave undecided are exactly the ones whose class the ring settles, and once the
    ring HAS been kekulised there is no aromatic bond left and the ordinary rows answer -- so the
    pyrrole-versus-pyridine nitrogen gets its count from the pass that resolved it.  A blanket
    recompute there is measured to LOSE counts: it overwrites the reader's correct 0 on diborane's
    bridging hydrogens and on ferrocene, where the arena holds a count no valence row can reproduce.
    Fill-only cannot, because it writes only where nothing is claimed.

    It also overlaps `stated` on purpose and does not replace it.  A count the record gave outright is
    already in the arena, so a reader that has stored its statements can simply pass `fill_only=True`;
    a reader mid-flight that has not yet needs to name them.

    Returns `{n: reason}` for every atom left undecided, each reason one of
    `HYD_NO_VALENCE_RULE`, `HYD_AMBIGUOUS_AROMATIC` or `HYD_NOT_AROMATIC_ELEMENT`.  Those atoms are
    written `H_UNKNOWN`, never a guessed zero -- a caller cannot tell an invented 0 from a real one
    and can tell `None`.  An atom answered by the second aromatic reading IS written and IS reported,
    with `HYD_DERIVED_OTHER_READING`, because the count then disagrees with the class `kekule()`
    will pick and a caller logging its input wants to know.

    TWO PASSES, and the split is not stylistic.  Reading the CSR needs a sealed arena; writing a
    count is a journal op.  So every decision is taken first, against one consistent structure, and
    the writes go out afterwards in a single edit scope -- which also means the seal recomputes
    whatever derives from a hydrogen count exactly once.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef dict undecided = {}
    if not n_atoms:
        return undecided

    cdef list numbers = molecule._numbers
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef dict index_of = molecule._index_of
    # a plain loop and not a set comprehension: a comprehension target lives in its own scope and
    # cannot be `cdef`-ed, so it comes out as `implicit declaration of`, which RULES.md §9.7 makes a
    # hard failure rather than a note
    cdef set skip = set()
    cdef uint32_t i, osum = 0, arom = 0, nbrs = 0, charged = 0, hn = 0, env_len = 0
    # `exo`, `reason` and the four above are all written through a pointer, and Cython cannot see
    # through one, so these initialisers are what keep the maybe-uninitialized warnings off
    cdef bint exo = False
    cdef uint8_t reason = HYD_DERIVED
    cdef uint16_t env[HYD_ENV_MAX]
    cdef list writes = []
    cdef object n, count
    if stated is not None:
        for n in stated:
            skip.add(<uint32_t> index_of[n])

    for i in range(n_atoms):
        if i in skip:
            continue
        a = atoms + i
        if fill_only and not at_implicit_h_unknown(a):
            continue
        # ONE SCRATCH BUFFER FOR THE WHOLE LOOP, per RULES.md: `env` is hoisted above and rewritten
        # from index 0 on every atom, so the pass allocates nothing per atom.
        hyd_atom_context(atoms, ptr, edges, i, &osum, &arom, &nbrs, &exo, env, &env_len)
        if hyd_derive_atom(a.element, a.charge, at_radical(a), osum, arom, nbrs, exo, False,
                           env, env_len, &hn, &charged, &reason):
            writes.append((numbers[i], hn))
            if reason != HYD_DERIVED:
                undecided[numbers[i]] = reason
        else:
            writes.append((numbers[i], H_UNKNOWN))
            undecided[numbers[i]] = reason

    with molecule.edit():
        for n, count in writes:
            molecule.set_hydrogens(n, count)
    return undecided


def derive_implicit_hydrogen(MoleculeContainer molecule, n, bint count_stated=False):
    """`(count, reason)` for one atom, WITHOUT writing it.  `count` is None when undecidable.

    The question-shaped half of `derive_implicit_hydrogens`, for a reader that ranks the derivation
    against other sources before it commits -- the CTfile reader puts a stated `MRV_IMPLICIT_H`
    above it and a stated total valence below it, and cannot use a function that writes.

    `count_stated=True` asserts that the caller's format states this atom's hydrogens, which
    withdraws the ambiguity gate; see the header.  A caller that does not know what that means wants
    the default.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t i = <uint32_t> molecule._index_of[n]
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a = atoms + i
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    # initialised for the same reason as in the sweep above: all of these are written by pointer
    cdef uint32_t osum = 0, arom = 0, nbrs = 0, charged = 0, hn = 0, env_len = 0
    cdef bint exo = False
    cdef uint8_t reason = HYD_DERIVED
    cdef uint16_t env[HYD_ENV_MAX]
    hyd_atom_context(atoms, ptr, edges, i, &osum, &arom, &nbrs, &exo, env, &env_len)
    if hyd_derive_atom(a.element, a.charge, at_radical(a), osum, arom, nbrs, exo, count_stated,
                       env, env_len, &hn, &charged, &reason):
        return hn, reason
    return None, reason


def _hyd_export():
    # Same route and same reason as `H_UNKNOWN` in `_molecule_container.pxi`: a `DEF`-or-enum name is
    # substituted textually wherever it appears as a name, assignment targets included, so
    # `HYD_DERIVED = HYD_DERIVED` compiles to `0 = 0`.  A string key is the one place the name
    # survives.  Exported because a format reader has to branch on these, and a caller restating
    # them as literals is a caller that drifts.
    globals()['HYD_DERIVED'] = HYD_DERIVED
    globals()['HYD_DERIVED_OTHER_READING'] = HYD_DERIVED_OTHER_READING
    globals()['HYD_NO_VALENCE_RULE'] = HYD_NO_VALENCE_RULE
    globals()['HYD_AMBIGUOUS_AROMATIC'] = HYD_AMBIGUOUS_AROMATIC
    globals()['HYD_REASON_MASK'] = HYD_REASON_MASK
    globals()['HYD_NO_AROMATIC_FORM'] = HYD_NO_AROMATIC_FORM


_hyd_export()
