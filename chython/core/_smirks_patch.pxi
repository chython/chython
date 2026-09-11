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
# The SMIRKS patcher: a ReactionTemplate applied to molecules, on the container's edit session.
#
# WHAT IT DOES, IN ONE PARAGRAPH
#
# Union the inputs into one working container.  For every embedding of the reactant side, copy the
# working container and edit the copy: delete the atoms that pair with nothing, build the atoms the
# product side creates, write the explicit fields of the atoms that pair, delete the bonds the product
# side does not restate, add or re-order the ones it does.  Then recompute the hydrogen counts the
# patch invalidated, report the parities the arena could not carry, drop the ones the reaction centre
# itself holds, split the result and hand back a ReactionContainer.  Dedupe on structure, guard every
# candidate, never mutate an input.
#
# WHAT IT DOES NOT DO
#
# It does not rebuild the molecule.  Under core decision D1 a parity is stored against CSR
# ascending-neighbour order, and the arena re-bases or drops every parity itself on each apply
# (`_harvest_parities` / `rebase_parity`) -- so there is no fresh container copied out in the
# replacement's connectivity order, no insertion-order property to preserve, and -- this is the
# design's N5 -- NO NEIGHBOUR-SET COMPARISON HERE AT ALL.  Guarding a parity carry with
# `sbonds[n].keys() == nbonds[n].keys()` says yes to a substitution replacing one neighbour with a
# differently-numbered atom, which is exactly the reaction at a stereocentre such a guard is reached
# for.  The arena's guard is positional instead, so this layer OBSERVES what the arena decided (N7),
# clears the group membership of anything it dropped (N8), and drops what the arena DID carry across a
# unit the patch wrote part of (`smk_rc_stereo`) -- a positional frame that still reads is not a
# configuration that survived the reaction.
#
# It does not derive stereogenicity either (N10): every question about which atoms can carry a
# configuration is the core's, and there is no candidate rule in this file.  The reaction-centre drop
# asks the core's own unit table which atoms a configuration is stated against; it asks nothing about
# which of them could hold one.
#
# HYDROGENS
#
# `h` is a CHECK primitive, so a product side cannot state a hydrogen count -- there is nothing to
# copy from the template.  The rule this layer picked, and the reason it is a rule and not a guess:
#
#   * an atom the patch WROTE (its element, charge, radical, or one of its bonds) and every surviving
#     neighbour of a deleted atom gets its count recomputed from the valence collection;
#   * an atom that merely sat inside the match keeps the count its input stated, exactly;
#   * where the collection has no answer -- no row for the element in that charge and radical state,
#     or an aromatic bond reaching the atom -- the count is stored as `H_UNKNOWN`.
#
# Never a guessed zero, and never an aromatic guess: `chython.chemistry.calc_implicit` writes
# H_UNKNOWN in the same two places for the same reason, and this file mirrors its policy in C because
# `core` cannot import `chemistry`.  The documented repair is the caller's: `kekule()` and then
# `calc_implicit`.
#
# THE THREE NUMBER SPACES (N6)
#
# A TEMPLATE atom is a stable id of one of the template's two sides.  A MOLECULE atom is a stable id
# of the working container.  A MAP NUMBER pairs the template's two sides, and a SECOND map number
# space is what this file writes onto the reaction it yields: contiguous from 1 over the atoms present
# on both sides, 0 for a leaving or an incoming one.  The two are never the same number -- the
# template's numbering is arbitrary and sparse, the output's is imposed.  Four dicts keep the spaces
# apart -- `react_to_prod`, `prod_to_new`, the matcher's `mapping` and `numbers` -- and no expression
# in this file uses one where another belongs.


# THE VALENCE-ENVIRONMENT WALK BELONGS TO `_hydrogens.pxi`, which declares its own bound as
# `HYD_ENV_MAX` and explains it there.  This file has no copy of the walk and no bound of its own:
# one derivation, one place its limits are declared.


# `ReactionContainer`, bound on first use.  `chython/core/__init__.py` imports `._core`
# FIRST and `.reaction` after it, because `reaction.py` imports from the extension -- so a module-level
# import here would run halfway through the package's own initialisation.  A log record is built by
# `mc_record`, which is in this translation unit, so nothing else needs binding.
cdef object _SMK_REACTION = None


cdef int smk_lazy_imports() except -1:
    global _SMK_REACTION
    # declared before they are imported, for the reason `_core.pyx` gives about `warn`: an import
    # statement binds a name Cython never saw declared, and `warn.undeclared` is on
    cdef object ReactionContainer
    if _SMK_REACTION is None:
        from chython.core.reaction import ReactionContainer
        _SMK_REACTION = ReactionContainer
    return 0


cdef inline str smk_rule_id(ReactionTemplate t):
    """What a log record calls this template: whatever `read_smirks` was told, or the string itself.

    A table-driven pass writes a table-qualified id (`'reactions:13'`) so the record names the ROW its
    author can go and edit; a template read from a string by hand has no such name, and the string is
    its identity, so it is what the record names.
    """
    return t.rule_id


cdef set smk_closure(MoleculeContainer work, set doomed, set remain):
    """Widen `doomed` by every fragment that only hung off it.

    A template deletes the atoms whose map number pairs with nothing.  Deleting one leaves whatever
    was bonded to it behind, and a leaving group is usually more than one atom: `[C:1][O:2][C:3]` with
    `:3` unmapped deletes the ether carbon and would otherwise emit that carbon's methyl as a
    free-floating CH3.  So from each neighbour of a deleted atom this walks outward WITHOUT CROSSING A
    DELETED ATOM, and if the fragment it reaches touches no surviving matched atom, the fragment goes
    too.  A template that means to KEEP an alkyl maps it -- absence is what says "delete", so the
    thing to write is `[C:1][O:2][C:3]>>[C:1][O;D1:2].[C:3][O;D1:4]`.

    THIS CANNOT REACH A COUNTER-ION, which is the property that matters for salts: the walk only ever
    starts at a NEIGHBOUR of a deleted atom, so a component with no bond into the reaction centre is
    never a candidate, however small it is.

    THE SEED IS TESTED TOO: a neighbour that is itself deleted starts no walk, so the "crosses no
    deleted atom" rule holds at the first step and not only inside the loop.
    """
    cdef set out = set(doomed)
    cdef set global_seen = set()
    cdef set seen
    cdef list stack
    cdef object x, n, current, other
    cdef bint reached
    for x in doomed:
        for n in work.neighbors_of(<uint32_t> x):
            if n in global_seen or n in remain or n in doomed:
                continue
            seen = {n}
            global_seen.add(n)
            stack = []
            for other in work.neighbors_of(<uint32_t> n):
                if other not in global_seen:
                    stack.append(other)
            reached = False
            while stack:
                current = stack.pop()
                if current in remain:
                    reached = True
                    break
                if current in doomed:
                    continue
                seen.add(current)
                global_seen.add(current)
                for other in work.neighbors_of(<uint32_t> current):
                    if other not in global_seen:
                        stack.append(other)
            if not reached:
                out |= seen
    return out


cdef int smk_hydrogens(MoleculeContainer m, uint32_t n) except -2:
    """The implicit hydrogen count the valence collection gives atom `n`, or `H_UNKNOWN`.

    ONE LINE OF POLICY AND NONE OF ITS OWN.  The walk over the CSR and the decision an aromatic atom
    needs both live in `_hydrogens.pxi`, which is included above this file, so the patcher asks
    instead of deciding.  Order 4 is not a state the collection has no answer for: it needs
    `arom_classify_atom`'s decision about whether the atom takes a ring double bond, and then the
    ordinary rows answer -- benzene's carbon is 1, so a template substituting an aryl ring gets a
    count.

    THE AMBIGUOUS ATOM STILL COMES BACK `H_UNKNOWN`, and so does a state with no row.  Those are the
    two honest failures and they are not zero: an atom whose hydrogens nobody can derive is not an
    atom with no hydrogens.  `kekule()` closes the first for a caller who runs it, which is what
    `chython/chemistry/test/test_reaction_hydrogen_repair.py` pins.
    """
    cdef uint32_t i = <uint32_t> m._index_of[n]
    cdef uint32_t *ptr = csr_ptr(m._structure)
    cdef halfedge_t *edges = csr_edges(m._structure)
    cdef atom_t *atoms = m._structure.atoms()
    cdef uint16_t env[HYD_ENV_MAX]
    # every one of these is written through a pointer, and Cython cannot see through one
    cdef uint32_t env_len = 0, osum = 0, arom = 0, nbrs = 0, hn = 0, charged = 0
    cdef bint exo = False
    cdef uint8_t reason = HYD_DERIVED
    hyd_atom_context(atoms, ptr, edges, i, &osum, &arom, &nbrs, &exo, env, &env_len)
    if hyd_derive_atom(atoms[i].element, atoms[i].charge, at_radical(&atoms[i]), osum, arom, nbrs,
                       exo, False, env, env_len, &hn, &charged, &reason):
        return <int> hn
    return H_UNKNOWN


# "the element is the R marker", in `smk_atom_fields`'s element slot, where 0 already means INHERIT.
# A distinct value because those two are distinct claims, and 118 elements sit between them.
DEF SMK_ELEMENT_R = -1


cdef tuple smk_atom_fields(tuple build, bint cx_radical):
    """`(element, charge, isotope, radical)` for one product atom.  EXPLICIT-ONLY.

    An unstated field is the DEFAULT and not the matched atom's value: charge zero, no mass number,
    not a radical.  That is what `read_smirks` documents and what its N1 lint reports a line about, so
    there is nothing to decide here -- only to apply.  Element 0 means INHERIT: the atom wrote `A`, or
    wrote no element at all and pairs with a reactant atom to take one.

    THE R MARKER IS ELEMENT 0 STATED OUTRIGHT -- `#0` -- so it cannot be spelled 0 here (RULES.md
    6.4).  It comes back as SMK_ELEMENT_R and `smk_element` translates it at the two write sites.

    `cx_radical` is the extension tail's `^N:` on this product atom, which is a radical statement in a
    second syntax and not a second field.
    """
    cdef int element = 0
    cdef int charge = 0
    cdef int isotope = 0
    cdef bint radical = cx_radical
    cdef tuple prim
    cdef int kind
    for prim in build:
        kind = <int> prim[0]
        if kind == PRIM_ELEMENT:
            element = <int> prim[1]
        elif kind == PRIM_R_MARKER:
            element = SMK_ELEMENT_R
        elif kind == PRIM_CHARGE:
            charge = <int> prim[1]
        elif kind == PRIM_ISOTOPE:
            isotope = <int> prim[1]
        elif kind == PRIM_RADICAL:
            radical = True
        # PRIM_ANY leaves the element at 0, which is the inheritance; PRIM_NO_ISOTOPE leaves the
        # isotope at 0, which is what it asks for; PRIM_STEREO is not a field of the atom's identity
        # and is applied afterwards, against a frame -- see `smk_stereo`
    return (element, charge, isotope, radical)


cdef inline int smk_element(tuple fields) noexcept:
    """The atomic number to write, with the R marker's sentinel translated back to element 0."""
    cdef int element = <int> fields[0]
    return 0 if element == SMK_ELEMENT_R else element


# ------------------------------------------------------------------------------------------------
# THE STEREO DIRECTIVES
# ------------------------------------------------------------------------------------------------
#
# `read_smirks` has already decided what the product side SAYS -- see the header above
# `smk_stereo_reading` for the four statements and for the refusals.  What is left here is:
#
#   `@~`   `3 - parity` at the unit's anchor.  Every kind chython models is two-state, so that is the
#          whole of the inversion, and it holds for a geometry as readily as for a parity.
#   `&<n>` a group, which CONFIGURES the unit and groups it, because a racemic unit is
#          configured-and-grouped and never unconfigured -- an unconfigured unit says "not known",
#          which is a different claim.  Any kind, since every kind is two-state: a mixture of a
#          cis/trans anchor's two states is an E/Z mixture.  Written with the groups, further down.
#   a drawn correlated set -> `smk_stereo_correlated`, the one place a frame is needed at all: a
#          relative configuration is a statement about several centres at once, so it has to be drawn
#          somewhere, and the product side's own arm order is that somewhere.
#
# `@=` NEEDS NO PASS.  Sparing the unit from `smk_rc_stereo`'s drop IS carrying the configuration
# through, because the arena re-based it already.  For the same reason `@~` consults no frame:
# `rebase_parity` did the work, and a flip of a value in a frame is a flip in every frame.
#
# NOTHING HERE INVENTS A CONFIGURATION, and `@~` is SILENT rather than loud where there is none to
# invert.  It is a conditional statement -- "whatever came in comes out the other way" -- so a
# substrate that came in unconfigured has nothing for it to be about, and that is what lets one corpus
# row carry it and still match an unconfigured substrate.  A template that states nothing at all
# reaches none of this: away from the reaction centre the atom keeps whatever the arena re-based (N10),
# and at the centre `smk_rc_stereo` drops it.

cdef set smk_stereo_invert(ReactionTemplate t, MoleculeContainer new, dict prod_to_new):
    """Apply every `@~`.  Returns the new-graph anchors whose parity was flipped.

    A stated atom is resolved to a unit the way `smk_rc_stereo` resolves one: either terminal of a bond
    kind names it, since which one anchors is a fact about slot order rather than about chemistry and a
    template cannot address it.  A unit holding no parity is left alone and not reported.
    """
    cdef set out = set()
    cdef set targets = set()
    cdef list numbers
    cdef list writes = []
    cdef Structure structure
    cdef stereo_unit_t *units
    cdef uint32_t k, count, partner
    cdef object sid, anchor, far
    cdef tuple write
    if not t.product_stereo_invert:
        return out
    for sid in t.product_stereo_invert:
        targets.add(prod_to_new[sid])
    numbers = new._numbers
    structure = new._structure
    new._require_clean()
    ensure_stereo_units(structure)
    # both taken after `ensure_stereo_units`, which reallocates the arena
    count = structure_stereo_unit_count(structure)
    units = structure_stereo_units(structure)
    for k in range(count):
        if not structure_parity_at(structure, units[k].anchor):
            continue
        anchor = numbers[units[k].anchor]
        partner = stereo_unit_partner(structure, &units[k])
        far = None if partner == SU_NO_REF else <object> numbers[partner]
        if anchor in targets or (far is not None and far in targets):
            out.add(anchor)
            # the OTHER value, of the two: the parity is three-state and 3 - parity swaps 1 and 2
            writes.append((anchor, 3 - structure_parity_at(structure, units[k].anchor)))
    if writes:
        # read out of the arena above and written here, because the parities the loop reads are the
        # ones the edit is about to move
        with new.edit():
            for write in writes:
                new.set_parity(<uint32_t> <int> write[0], <int> write[1])
    return out


# ------------------------------------------------------------------------------------------------
# N12: RELATIVE CONFIGURATION, WHICH IS A CORRELATED GROUP AND NOT A NEW TOKEN
# ------------------------------------------------------------------------------------------------
#
# An AND group already means "as drawn, or all members flipped", which IS a fixed relative
# configuration presented as a racemate.  A diastereoselective reaction is that same statement with
# one member's branch already settled by the substrate, so N12 needs no token and no kernel change --
# only the two cases below, and which one applies is read off the MOLECULE rather than the string:
#
#   * no member arrives configured  -> write the drawn set as drawn, keep the group.  The honest
#     product of an achiral substrate: one diastereomer, both enantiomers.
#   * some member arrives configured -> the set is RESOLVED against it.  The drawn signs are mirrored
#     as a whole if that is what agreeing with the carried centre takes, and no group is written: the
#     answer is a single diastereomer of known absolute configuration, because the substrate's was.
#
# Mirroring "as a whole" is the entire content of the requirement.  Flipping one member would change
# which diastereomer the template describes; flipping all of them keeps the relative configuration and
# changes only which enantiomer of it, which is exactly the freedom the substrate is being allowed to
# remove.  So one member's carried configuration decides a single bit for the whole group.
#
# The drawn value is read through `translate_stereo` rather than by re-deriving a permutation sign,
# for the reason the deleted frame path had: the arena already knows how to read a stored parity in a
# caller's order, and the map between the stored frame and the caller's is a bijection on {1,2}.

cdef tuple smk_stereo_correlated(ReactionTemplate t, MoleculeContainer new, dict prod_to_new,
                                 str rule, object log):
    """Apply every correlated group.  Returns `(new ids a sign landed on, product atoms to leave ungrouped)`.

    A member the patched molecule cannot configure is logged and skipped, as everywhere else in this
    file.  When that leaves the group with fewer than two members it is no longer a correlation, and
    what remains is written anyway: a single drawn configuration in a group is N3 case 1, a racemate,
    and dropping the statement entirely would lose the fact that the centre is a mixture at all.
    """
    cdef set stated = set()
    cdef set resolved = set()
    cdef dict corr = t.product_stereo_correlated
    cdef dict groups = t.product_stereo_groups
    cdef dict buckets = {}
    cdef dict unit
    cdef list members, order
    cdef object sid, ref, spec, key
    cdef tuple entry, drawn_spec
    cdef uint32_t nid
    cdef int drawn, here, want, stored, carried, flip
    if not corr:
        return stated, resolved
    for sid in sorted(corr):
        spec = groups[sid]
        members = <list> buckets.get(spec)
        if members is None:
            members = []
            buckets[spec] = members
        members.append(sid)

    for key in sorted(buckets):
        # every member's drawn value and what the arena currently says, read BEFORE anything is written
        members = []
        for sid in <list> buckets[key]:
            drawn_spec = <tuple> corr[sid]
            drawn = <int> drawn_spec[0]
            nid = <uint32_t> prod_to_new[sid]
            unit = <dict> new.unit_of(nid)
            if unit is None or <int> unit['kind'] != SU_TETRA:
                log.append(mc_record(rule, (<object> nid,),
                                    'the template draws atom %d as one member of a correlated '
                                    'stereo group and the patched molecule holds no tetrahedral '
                                    'centre there, so the member was skipped' % nid))
                continue
            order = []
            for ref in <tuple> drawn_spec[1]:
                order.append(None if ref is None else prod_to_new[ref])
            try:
                here = new.translate_stereo(nid, tuple(order))
            except (KeyError, ValueError) as exc:
                log.append(mc_record(rule, (<object> nid,),
                                    'the template draws atom %d against directions the patched '
                                    'molecule does not have, so the member was skipped: %s'
                                    % (nid, exc)))
                continue
            members.append((sid, <object> nid, drawn, tuple(order), here))
        if not members:
            continue

        # ONE BIT FOR THE WHOLE SET.  The lowest-numbered member that arrives configured decides it;
        # which member is asked cannot matter, because the drawn signs already fix the others relative
        # to it, and a set with no configured member is free to be drawn as written.
        carried = 0
        flip = 0
        for entry in members:
            if <int> entry[4]:
                carried = 1
                flip = 1 if <int> entry[4] != <int> entry[2] else 0
                break
        if carried:
            # the whole GROUP is resolved, not just the members that could be written: it is one
            # statement about a set, and leaving a skipped member grouped would claim a mixture the
            # rest of the set no longer is
            for sid in <list> buckets[key]:
                resolved.add(sid)

        for entry in members:
            nid = <uint32_t> <int> entry[1]
            drawn = <int> entry[2]
            here = <int> entry[4]
            want = (3 - drawn) if flip else drawn
            if here != want:
                stored = new.parity_of(nid)
                if here:
                    # the map between the frames is settled by this one read, and there are two values
                    with new.edit():
                        new.set_parity(nid, 3 - stored)
                else:
                    # nothing was stored, so the read said nothing about the permutation: write and ask
                    with new.edit():
                        new.set_parity(nid, 1)
                    if new.translate_stereo(nid, <tuple> entry[3]) != want:
                        with new.edit():
                            new.set_parity(nid, 2)
            stated.add(<object> nid)
    return stated, resolved


# ------------------------------------------------------------------------------------------------
# A DRAWN GEOMETRY: THE ONE ABSOLUTE CONFIGURATION A TEMPLATE STATES
# ------------------------------------------------------------------------------------------------
#
# `read_smirks` has already found the two terminals and the two marked substituents and reduced the
# `/` and `\` to one bit, `trans` (see the header above `smk_directions` for why they are read from the
# terminals and for the refusals).  What is left is the same write-and-ask this file uses for a drawn
# correlated set, over a frame of four instead of a permutation: in the order
# `(marked, other, marked, other)` the convention reads straight off -- even is trans -- so the frame is
# built by putting each terminal's MARKED substituent first within its own pair.
#
# The pair ORDER within `refs` is not a choice: ruling F26 puts the anchor's own directions first, so
# which terminal anchors decides which of the two marked atoms leads.  Exchanging the two pairs is even
# and would read the same, but taking the anchor's is what makes the frame a permutation of `refs` that
# `translate_stereo` accepts at all.

cdef set smk_stereo_geometry(ReactionTemplate t, MoleculeContainer new, dict prod_to_new,
                             str rule, object log):
    """Apply every geometry the product side drew.  Returns the new-graph anchors one landed on.

    A chain the patched molecule holds no cis/trans unit for -- an allene, a ring too small, a terminal
    the patch left with two identical directions -- is logged and skipped, as everywhere else here.
    """
    cdef set stated = set()
    cdef dict geom = t.product_stereo_geometry
    cdef dict unit
    cdef tuple key, spec, refs, marks, want
    cdef list order
    cdef object nid, other, near, far
    cdef uint32_t k, n1, n2
    cdef int desired, here, stored
    cdef bint ok
    cdef object exc                      # the except-as target; `warn.undeclared` counts it
    if not geom:
        return stated
    for key in sorted(geom):
        spec = <tuple> geom[key]
        n1 = <uint32_t> prod_to_new[key[0]]
        n2 = <uint32_t> prod_to_new[key[1]]
        # either terminal may anchor -- which one is a fact about slot order, not about chemistry
        nid = <object> n1
        unit = <dict> new.unit_of(n1)
        if unit is None or <int> unit['kind'] != SU_CIS_TRANS:
            nid = <object> n2
            unit = <dict> new.unit_of(n2)
        if unit is None or <int> unit['kind'] != SU_CIS_TRANS or <int> unit['n_refs'] != 4:
            log.append(mc_record(rule, (<object> n1, <object> n2),
                                'the template draws a geometry across the bond between atoms %d and '
                                '%d and the patched molecule holds no cis/trans unit there, so it was '
                                'skipped' % (n1, n2)))
            continue
        refs = <tuple> unit['refs']
        if nid == <object> n1:
            marks = (<object> prod_to_new[spec[0]], <object> prod_to_new[spec[1]])
        else:
            marks = (<object> prod_to_new[spec[1]], <object> prod_to_new[spec[0]])
        order = []
        ok = True
        for k in range(2):
            near = refs[2 * k]
            far = refs[2 * k + 1]
            if marks[k] == near:
                order.append(near)
                order.append(far)
            elif marks[k] == far:
                order.append(far)
                order.append(near)
            else:
                ok = False
                log.append(mc_record(rule, (<object> nid, marks[k]),
                                    'the template marks atom %s as a direction of the geometry '
                                    'anchored at atom %d and the patched molecule does not have it '
                                    'there, so the geometry was skipped' % (marks[k], <uint32_t> nid)))
                break
        if not ok:
            continue
        want = tuple(order)
        # `_smiles_read.pxi`'s convention, in the frame just built: even means positions 0 and 2 --
        # the two marked substituents -- are trans
        desired = 1 if <bint> spec[2] else 2
        try:
            here = new.translate_stereo(<uint32_t> nid, want)
        except (KeyError, ValueError) as exc:
            log.append(mc_record(rule, (<object> nid,),
                                'the template draws a geometry at atom %d against directions the '
                                'patched molecule does not order that way, so it was skipped: %s'
                                % (<uint32_t> nid, exc)))
            continue
        if here != desired:
            stored = new.parity_of(<uint32_t> nid)
            if here:
                # the map between the frames is settled by this one read, and there are two values
                with new.edit():
                    new.set_parity(<uint32_t> nid, 3 - stored)
            else:
                # nothing was stored, so the read said nothing about the permutation: write and ask
                with new.edit():
                    new.set_parity(<uint32_t> nid, 1)
                if new.translate_stereo(<uint32_t> nid, want) != desired:
                    with new.edit():
                        new.set_parity(<uint32_t> nid, 2)
        stated.add(nid)
    return stated


cdef dict smk_group_anchors(MoleculeContainer new, dict groups, str rule, object log):
    """`groups` re-keyed onto the anchor of whatever stereogenic unit each atom belongs to.

    A bond kind's FAR TERMINAL anchors nothing, so `unit_of` answers None there while the unit is
    perfectly addressable -- the fact `smk_rc_stereo` handles by testing the partner.  A group has to
    move onto the anchor because a parity and an enhanced-stereo byte are one statement about one unit
    and `_stereo_record_flips` reads both off the anchor; a byte left on the far terminal would be a
    mixture nothing honours.  An atom no stereogenic unit reaches is left where it is, for the
    `stereogenic` test below to drop with its own message.

    Both terminals grouped is one statement made twice: the lower id's group wins, and the collision is
    reported rather than resolved silently, because the two brackets may name different numbers.
    """
    cdef list numbers = new._numbers
    cdef Structure structure = new._structure
    cdef stereo_unit_t *units
    cdef uint32_t k, count, partner
    cdef dict anchors = {}
    cdef dict out = {}
    cdef object sid, anchor
    new._require_clean()
    ensure_stereo_units(structure)
    count = structure_stereo_unit_count(structure)
    units = structure_stereo_units(structure)
    for k in range(count):
        if not (units[k].spare & SU_STEREOGENIC):
            continue
        partner = stereo_unit_partner(structure, &units[k])
        if partner != SU_NO_REF:
            anchors[<object> numbers[partner]] = <object> numbers[units[k].anchor]
    for sid in sorted(groups):
        anchor = anchors.get(sid, sid)
        if anchor in out:
            log.append(mc_record(rule, (sid,),
                                'the template puts atom %d in an enhanced-stereo group and another '
                                'atom of the same stereo unit in one too; one unit takes one group, so '
                                'the group on atom %d was not written' % (sid, sid)))
            continue
        out[anchor] = groups[sid]
    return out


cdef dict smk_stereo_group_ids(ReactionTemplate t, MoleculeContainer new, dict prod_to_new,
                               set resolved):
    """`{new id: (kind, group)}` for the product side's `&<n>` / `o<n>` brackets and `|a:|` tail.

    `resolved` names the product atoms N12 settled against a configuration the substrate carried: they
    state a single diastereomer rather than a mixture, so they are skipped here and no id is allocated
    for a group that would have had no members.

    N3: a template's group NUMBERS are template-local and mean nothing outside it, so every distinct
    `(kind, number)` the template names is allocated an id that no group in the patched molecule is
    using.  Two atoms the template puts in one group land in one group here; a number the product
    happens to share with the input does not merge them.

    `abs` carries no number and needs no allocation.  Running out of ids -- all 63 in use, which no
    real molecule reaches -- drops the group rather than reusing an occupied one: a wrong group is a
    false claim about a mixture, and no group is only a missing one.
    """
    cdef dict out = {}
    cdef dict wanted = t.product_stereo_groups
    cdef object sid, key, taken
    cdef tuple spec
    cdef int kind, number
    cdef dict allocated = {}
    cdef set used = set()
    if not wanted:
        return out
    for key in new.stereo_groups():
        used.add(<int> (<tuple> key)[1])
    cdef int free = 1
    for sid in sorted(wanted):
        spec = <tuple> wanted[sid]
        kind = <int> spec[0]
        number = <int> spec[1]
        if sid not in prod_to_new or sid in resolved:
            continue
        if kind == SMI_SG_ABS:
            out[prod_to_new[sid]] = (kind, 0)
            continue
        key = (kind, number)
        taken = allocated.get(key)
        if taken is None:
            while free < 64 and free in used:
                free += 1
            if free > 63:
                continue
            taken = free
            allocated[key] = taken
            used.add(free)
        out[prod_to_new[sid]] = (kind, <int> taken)
    return out


# ------------------------------------------------------------------------------------------------
# THE CONFIGURATION AT THE REACTION CENTRE IS DROPPED
# ------------------------------------------------------------------------------------------------
#
# A patch that made or broke a bond at a stereo unit leaves that unit's configuration UNSTATED, for
# both kinds and every template.  The arena's re-base carries a parity whenever the frame still reads
# (N5's positional guard), and carrying it is the wrong answer where the reaction happened AT the
# centre: an SN2 substitution keeps the frame -- one neighbour replaced by another in the same slot --
# so the re-based sign would claim retention that the template never stated and that most courses do
# not have.  Unstated is the honest answer to a question the template did not answer.
#
# A template that knows better says so, in one of three ways, all of which take precedence here:
#
#   * `@=` on any atom of the unit  -> the configuration comes through unchanged.  What a cut states:
#     a fragment leaving takes no configuration with it.
#   * `@~` on any atom of the unit  -> it comes through as the unit's other state.
#   * `&<n>` / `o<n>` -> a mixture, which N3 writes as a configured-and-grouped centre.
#   * a sign, correlated with another product atom's -> a drawn relative configuration.
#
# WHAT COUNTS AS "AT THE CENTRE" is `changed`, the same set the hydrogen recompute runs on -- an atom
# whose element, charge, radical or bonds the patch wrote, or that lost a neighbour to a deletion.  A
# unit is at the centre when one of its OWN atoms is: the anchor for an atom kind, either terminal for
# a bond kind.  A substituent is not, and that is the line that matters -- acylating an amine two bonds
# from a stereocentre writes a neighbour of one of its directions and no part of the unit, so the
# configuration stands.  An atom that merely sat inside the match is not in `changed` at all.

cdef list smk_rc_stereo(MoleculeContainer new, set changed, set spared):
    """The anchors of every configured stereo unit the patch wrote a part of, ascending.

    `spared` names the atoms some directive already spoke for -- `@=`, `@~`, a sign, a group -- and spares
    the whole unit each is part of, since which terminal of a bond kind anchors it is a fact about slot
    order rather than about chemistry (`chiral_bonds` says why) and a template cannot address it.
    """
    cdef list numbers = new._numbers
    cdef Structure structure = new._structure
    cdef stereo_unit_t *units
    cdef uint32_t k, count, partner
    cdef list out = []
    cdef object anchor, far
    if not changed:
        return out
    new._require_clean()
    ensure_stereo_units(structure)
    # both taken after `ensure_stereo_units`, which reallocates the arena
    count = structure_stereo_unit_count(structure)
    units = structure_stereo_units(structure)
    for k in range(count):
        if not structure_parity_at(structure, units[k].anchor):
            continue
        anchor = numbers[units[k].anchor]
        partner = stereo_unit_partner(structure, &units[k])
        far = None if partner == SU_NO_REF else <object> numbers[partner]
        if anchor in spared or (far is not None and far in spared):
            continue
        if anchor in changed or (far is not None and far in changed):
            out.append(anchor)
    out.sort()
    return out


# ---------------------------------------------------------------------------
# N4: THE PRODUCT-SIDE POST-FILTER
# ---------------------------------------------------------------------------
#
# A product primitive either builds or checks, and the reader has already split them (`smk_place`).
# The build half is written above.  This is the other half: `D`, `h`, `H`, `x`, `z`, `r`, `R`, `M`
# the metal test and `@` the ring bond, asked of the PATCHED molecule.  That is what makes a
# cyclization state its ring size as a product-side `r` instead of as an argument -- the ring the
# primitive describes exists in the graph the primitive is read against, so there is no relational
# vocabulary and no new kernel.
#
# WHY THE BITS COME FROM `prim_apply` AND NOT FROM THE CONTAINER'S ACCESSORS.  `D`, `x` and `z`
# exclude a dative bond and `degree_of()` / `heteroatoms_of()` do not (both counts are correct;
# `_features.pxi` says why).  An `h` the input never recorded is H_UNKNOWN, which answers NEITHER
# `h2` nor `!h2`.  A ring size above 24 is bucketed.  Every one of those is a decision already made
# once, in the box compiler, so this filter compiles the check clauses into `wbox_t`s with the same
# `prim_apply` the matcher's boxes are built with and runs the kernel's own four ANDs against the
# molecule's feature words.  A product-side `D` therefore means exactly what a reactant-side `D`
# means, and there is no second definition of any primitive in this file.
#
# `box_fill_defaults` is deliberately NOT called.  It supplies a query atom's unstated defaults --
# neutral charge above all -- and a check box has no business stating them: the charge is the build
# half's to write, and a filled box would reject every product atom the template deliberately
# charged.  A box built only from check primitives forbids only what those primitives forbid.
#
# A rejection is not an error.  It is the ordinary negative outcome of a test, exactly as a
# reactant side that fails to match is, and the candidate is simply not yielded.  It does get ONE
# log line, unlike the matcher's silence, for a reason particular to arriving late: a template whose
# reactant side never matched produces no reaction and the author can see that from the pattern,
# while a template that matched and was then rejected here is otherwise indistinguishable from one
# that never matched at all.  One line per rejected candidate, naming the primitive in the author's
# own notation -- the same granularity as the exception path in `smk_enumerate`.


cdef inline bint smk_wbox_admits(wbox_t *box, uint64_t *f, uint64_t w0) noexcept nogil:
    """One check box against one patched atom or bond: the kernel's test, on a box with no arena.

    `box_admits` cannot be reused directly -- it reads a sealed `q_box_t` out of a query arena and
    there is no query here -- so this is the same four ANDs and the same multi-hot loop over the
    construction-time struct.  Kept adjacent to that comment on purpose: if the kernel's test ever
    grows a fifth term, both copies have to grow it.
    """
    cdef uint32_t k
    if w0 & box.neg[0]:
        return False
    if f[1] & box.neg[1] or f[2] & box.neg[2] or f[3] & box.neg[3]:
        return False
    for k in range(box.any_count):
        if not (f[box.any_word[k]] & box.any_mask[k]):
            return False
    return True


cdef int smk_admits(tuple clauses, uint64_t *f, uint64_t w0) except -1:
    """Every clause of one key's check, against one patched atom or bond.  1 passes, 0 fails.

    Returns the index of the failing clause plus one on failure -- 0 when every clause holds -- so
    the caller can spell the primitive that did it.  A clause is a `,` disjunction of alternatives
    and each alternative a `&` conjunction, so one box per alternative and the clause passes on the
    first box that admits.
    """
    cdef wbox_t box
    cdef Py_ssize_t i
    cdef object alt, prim
    cdef bint ok
    for i in range(len(clauses)):
        ok = False
        for alt in <tuple> clauses[i]:
            memset(&box, 0, sizeof(wbox_t))
            for prim in <tuple> alt:
                prim_apply(&box, <uint32_t> <int> (<tuple> prim)[0], <int32_t> <int> (<tuple> prim)[1],
                           <bint> (<tuple> prim)[2])
            if smk_wbox_admits(&box, f, w0):
                ok = True
                break
        if not ok:
            return <int> i + 1
    return 0


cdef tuple smk_post_filter(ReactionTemplate t, MoleculeContainer new, dict prod_to_new, set alive):
    """N4, applied.  `(atoms, message)` for the first check that fails, or None when they all hold.

    Runs on the FINISHED molecule and nowhere earlier: `r` and `R` need the ring perception of the
    graph the patch produced, and `h` and `H` need the counts the hydrogen recompute wrote.  Both
    are re-derived on every seal, so by the time this is called there is nothing left to wait for.
    """
    cdef Structure st = new._structure
    cdef uint64_t *feat = structure_features(st) + 4      # past `fill_features`' union row
    cdef uint64_t *edge_words = structure_edge_words(st)
    cdef uint32_t *ptr = csr_ptr(st)
    cdef halfedge_t *edges = csr_edges(st)
    cdef halfedge_t *e
    cdef uint64_t zero[4]
    cdef uint64_t *f
    cdef dict index_of = new._index_of
    cdef object key, sid, u, v
    cdef tuple clauses
    cdef uint32_t idx
    cdef int bad
    memset(&zero[0], 0, sizeof(zero))

    for sid in sorted(t.product_atom_check):
        clauses = <tuple> t.product_atom_check[sid]
        if not clauses:
            continue
        key = prod_to_new.get(sid)
        if key is None or key not in alive:
            continue
        idx = <uint32_t> <int> index_of[key]
        f = feat + 4 * idx
        # the atom's own aggregate word 0, which is the ROOT case of `atom_admits`' contract: no
        # bond is folded into a check box (only `M` and an element touch word 0, and only in the
        # element span), so the aggregate word's exact element bits settle everything word 0 says
        bad = smk_admits(clauses, f, f[0])
        if bad:
            return ((key,), 'the patched molecule fails the product-side check `%s` at atom %d '
                            '(%s), so this match produced nothing'
                    % (smk_check_spelling(<tuple> clauses[bad - 1]), <int> key,
                       smk_atom_name(<uint32_t> <int> sid, t.product_map_numbers)))

    for key in t.product_bonds:
        clauses = <tuple> t.product_bond_check[key]
        if not clauses:
            continue
        u = prod_to_new.get((<tuple> key)[0])
        v = prod_to_new.get((<tuple> key)[1])
        if u is None or v is None or u not in alive or v not in alive:
            continue
        e = csr_find_at(ptr, edges, <uint32_t> <int> index_of[u], <uint32_t> <int> index_of[v])
        if e is NULL:
            continue
        # the half-edge word, never the atom's aggregate: a bond primitive asks about ONE bond and
        # the aggregate ORs every incident bond's topology bits together
        bad = smk_admits(clauses, &zero[0], edge_words[e - edges])
        if bad:
            return ((u, v), 'the patched molecule fails the product-side check `%s` on the bond '
                            'between atoms %d and %d, so this match produced nothing'
                    % (smk_check_spelling(<tuple> clauses[bad - 1]), <int> u, <int> v))
    return None


cdef tuple smk_numbering(dict origin, set alive, set touched):
    """The imposed mapping and how many atoms it could not reach: `({work atom id: map number}, over)`.

    Contiguous from 1 over the atoms present on BOTH sides -- a work atom of a touched input that the
    patch left alive.  A leaving atom is absent from `alive` and a created one from `origin`, and both
    stay 0.  Ascending atom id, which groups the numbering by input for free: `union` allocates one
    ascending block per input and an id is never reused.

    CAPPED AT `MAP_NUMBER_MAX`, and the overflow is left at 0 rather than raised.  `atom_t.map_number`
    is 16 bits with a declared ceiling, so a substrate above it cannot be mapped 1-1 at all -- and
    refusing there would mean a peptide-sized input produces no reaction rather than an unmapped one.
    0 already means "paired with nothing known", so the partial answer is in the domain; the caller
    gets a log line and the reaction.
    """
    cdef dict numbers = {}
    cdef Py_ssize_t over = 0
    cdef object sid
    for sid in sorted(alive):
        if origin.get(sid) in touched:
            if len(numbers) < MAP_NUMBER_MAX:
                numbers[sid] = len(numbers) + 1
            else:
                over += 1
    return (numbers, over)


cdef MoleculeContainer smk_write_numbers(MoleculeContainer molecule, dict numbers):
    """A copy of `molecule` with `numbers` written and 0 everywhere else.

    One arena clone and nothing re-derived: `OP_SET_MAP_NUMBER` is on `_apply`'s non-invalidating list,
    so this drops no CIP label and harvests no parity.  Values are read before the scope opens, since a
    container with a pending journal refuses a read.
    """
    cdef MoleculeContainer out = molecule.copy()
    cdef list writes = []
    cdef object sid, pair
    for sid in out.atom_numbers:
        writes.append((sid, numbers.get(sid, 0)))
    with out.edit():
        for pair in writes:
            out.set_map_number(<uint32_t> (<tuple> pair)[0], <int> (<tuple> pair)[1])
    return out


cdef tuple smk_one(ReactionTemplate t, MoleculeContainer work, set work_bonds, dict origin,
                   dict origin_n, tuple snapshots, dict mapping, object log):
    """One embedding into one `(reaction, dedupe key)` pair, or None when it patches to nothing.

    Everything this raises is caught by the caller and logged (N11).  Nothing here touches `work`,
    which every candidate of the enumeration reads.
    """
    cdef str rule = smk_rule_id(t)
    cdef dict react_to_prod = {}
    cdef dict prod_to_new = {}
    cdef object number, pair, sid, key, part, u, v
    for number in t.mapped_pairs:
        pair = <tuple> t.mapped_pairs[number]
        react_to_prod[pair[0]] = pair[1]
        prod_to_new[pair[1]] = mapping[pair[0]]

    # DELETION BY ABSENCE.  `deleted_atoms` is already `M`-exempt: `read_smirks` skips a masked atom
    # when it works the set out, so an atom named purely as context is here on neither list.
    cdef set doomed = set()
    for sid in t.deleted_atoms:
        doomed.add(mapping[sid])
    cdef set remain = set()
    for sid in mapping.values():
        if sid not in doomed:
            remain.add(sid)
    if doomed:
        doomed = smk_closure(work, doomed, remain)

    cdef set touched = set()
    for sid in mapping.values():
        touched.add(origin[sid])

    cdef MoleculeContainer new = work.copy()
    # EVERY configured parity, snapshotted before the edit.  O(atoms), and correct without reasoning
    # about how far the patch reaches -- which is the point: the arena decides what survives, and a
    # cheaper snapshot would be this file guessing at that decision again.
    cdef set was_configured = set()
    for sid in new.atom_numbers:
        if new.parity_of(<uint32_t> sid):
            was_configured.add(sid)

    # atoms whose hydrogen count the patch invalidates.  A deleted atom takes its bonds with it, so
    # every surviving neighbour of one is on the list before the edit even opens.
    cdef set changed = set()
    for sid in doomed:
        for key in work.neighbors_of(<uint32_t> sid):
            if key not in doomed:
                changed.add(key)

    cdef set prod_bond_set = set(t.product_bonds)
    cdef tuple fields
    cdef uint32_t nid
    cdef int order
    with new.edit():
        for sid in sorted(doomed):
            new.delete_atom(<uint32_t> sid)

        for sid in sorted(t.created_atoms):
            fields = smk_atom_fields(<tuple> t.product_atom_build[sid], sid in t.product_radicals)
            # `implicit_h` is deliberately omitted, so the atom is born H_UNKNOWN and the recompute
            # below is what gives it a number -- once its bonds exist
            nid = new.add_atom(smk_element(fields), charge=<int> fields[1],
                               isotope=<int> fields[2], radical=<bint> fields[3])
            prod_to_new[sid] = nid
            changed.add(<object> nid)

        for sid in sorted(react_to_prod.values()):
            nid = <uint32_t> prod_to_new[sid]
            fields = smk_atom_fields(<tuple> t.product_atom_build[sid], sid in t.product_radicals)
            # written only where the value DIFFERS: an equal write would land this atom on the
            # hydrogen-recompute list, and recomputing turns a count the input left unknown into a
            # number the template never asked for
            if fields[0] and smk_element(fields) != work.element_of(nid):
                new.set_element(nid, smk_element(fields))
                changed.add(<object> nid)
            if <int> fields[1] != work.charge_of(nid):
                new.set_charge(nid, <int> fields[1])
                changed.add(<object> nid)
            if <int> fields[2] != work.isotope_of(nid):
                new.set_isotope(nid, <int> fields[2])
            if <bint> fields[3] != work.radical_of(nid):
                new.set_radical(nid, <bint> fields[3])
                changed.add(<object> nid)

        # a reactant bond survives when the product side RESTATES it between the partners of its two
        # endpoints.  An endpoint that pairs with nothing is deleted or masked, and either way the
        # bond is not this loop's business.
        for pair in t.reactant_bonds:
            u = react_to_prod.get((<tuple> pair)[0])
            v = react_to_prod.get((<tuple> pair)[1])
            if u is None or v is None:
                continue
            if ((u, v) if u < v else (v, u)) in prod_bond_set:
                continue
            u = mapping[(<tuple> pair)[0]]
            v = mapping[(<tuple> pair)[1]]
            if u in doomed or v in doomed:
                continue
            new.delete_bond(<uint32_t> u, <uint32_t> v)
            changed.add(u)
            changed.add(v)

        for pair in t.product_bonds:
            order = smk_bond_order(<tuple> t.product_bond_build[pair])
            u = prod_to_new[(<tuple> pair)[0]]
            v = prod_to_new[(<tuple> pair)[1]]
            # add-versus-set is decided from the PRE-EDIT state, because `add_bond` inside an open
            # scope does not check for an existing bond -- it cannot, the arena still holds the
            # pre-scope graph -- and a second add of one bond is a duplicate edge at apply time
            if ((u, v) if u < v else (v, u)) in work_bonds:
                if work.order_of(<uint32_t> u, <uint32_t> v) != order:
                    new.set_order(<uint32_t> u, <uint32_t> v, order)
                    changed.add(u)
                    changed.add(v)
            else:
                new.add_bond(<uint32_t> u, <uint32_t> v, order)
                changed.add(u)
                changed.add(v)

    cdef set alive = set(new.atom_numbers)

    # Hydrogens BEFORE the stereo directives, and the order is load-bearing: whether an anchor's
    # fourth direction is an atom or its implicit hydrogen is a fact about the count, so a directive
    # read against a stale count would be read against the wrong frame.
    cdef list writes = []
    for sid in sorted(changed):
        if sid in alive:
            writes.append((sid, smk_hydrogens(new, <uint32_t> sid)))
    if writes:
        with new.edit():
            for pair in writes:
                new.set_hydrogens(<uint32_t> (<tuple> pair)[0], <int> (<tuple> pair)[1])

    # `@~`: the other of the unit's two states
    cdef set stated = smk_stereo_invert(t, new, prod_to_new)
    # and a drawn set whose members were correlated in one group.  Separate pass rather than a branch
    # inside the first, because a correlated group is decided for all of its members at once and an
    # atom-at-a-time loop has nowhere to put that decision.
    cdef tuple correlated = smk_stereo_correlated(t, new, prod_to_new, rule, log)
    stated |= <set> correlated[0]
    # and the geometries `/` and `\` drew.  Third pass rather than a branch, because it is the only one
    # whose statement is absolute and the only one whose frame spans two atoms.
    stated |= smk_stereo_geometry(t, new, prod_to_new, rule, log)

    # N7: WHAT THE ARENA DECIDED, READ BACK.  A parity that was configured, whose atom is still here
    # and whose sign is now unset, is one `rebase_parity` refused to carry -- more than one leftover
    # position on a side of the frame, or the unit's kind changed under it.  A silently correct
    # RE-BASE gets no record, so a reaction that changes nothing at a labelled centre emits zero of
    # these, which is N7's negative control.  An atom a directive wrote is not a loss and not
    # reported: the template SAID what the configuration is, which supersedes whatever was carried.
    cdef list dropped = []
    for sid in sorted(was_configured):
        if sid in alive and sid not in stated and not new.parity_of(<uint32_t> sid):
            dropped.append(sid)

    # N3: the template's own groups, renumbered.  An atom the template groups is not a candidate for
    # the N8 clear -- the group it is about to get is the group it should have.
    cdef dict groups = smk_stereo_group_ids(t, new, prod_to_new, <set> correlated[1])
    if groups:
        groups = smk_group_anchors(new, groups, rule, log)

    # THE REACTION CENTRE'S OWN CONFIGURATION, dropped unless a directive spoke for it: `@=`, `@~`, a
    # group, or a correlated sign.  After N7's report, because these two drops are different facts
    # about the molecule and one message may not stand in for the other -- N7 says the frame stopped
    # reading, this says the reaction happened here.
    cdef set spared = set(stated)
    for sid in groups:
        spared.add(sid)
    for sid in t.product_stereo_keep:
        spared.add(prod_to_new[sid])
    for sid in t.product_stereo_invert:
        spared.add(prod_to_new[sid])
    cdef list rc = smk_rc_stereo(new, changed, spared)

    cdef list clears = []
    for sid in dropped:
        if sid not in groups and new.stereo_group_of(<uint32_t> sid)[0]:
            clears.append(sid)
    for sid in rc:
        if new.stereo_group_of(<uint32_t> sid)[0]:
            clears.append(sid)

    # N3: A GROUPED UNIT IS A CONFIGURED ONE.  `&<n>` on a unit the patch left unconfigured is the
    # racemize spelling, and a group written beside parity 0 would say "unconfigured, and in a mixture"
    # -- two claims that cannot both be true, and the exact conflation N3 exists against.  WHICH parity
    # is written does not matter and is not a choice this
    # file is making: a one-member AND group means this configuration and its mirror, an OR group
    # means one of the two and nobody knows which, and both are symmetric in the value.
    #
    # EVERY KIND, not just tetrahedral.  A group states that a unit's two states are both present, and
    # each kind chython models has exactly two -- so `&<n>` on a cis/trans anchor is an E/Z mixture for
    # the same reason it is a racemate on a centre, which is the answer a template with no facial or
    # geometric control has.  What gates it is `stereogenic` alone: a group on a CH2 or on a
    # 1,1-disubstituted alkene is a false claim about a mixture, and the group is dropped with it below.
    cdef list configure = []
    cdef dict unit
    for sid in sorted(groups):
        if new.parity_of(<uint32_t> sid) or (<tuple> groups[sid])[0] == SMI_SG_ABS:
            continue
        unit = <dict> new.unit_of(<uint32_t> sid)
        if unit is not None and <bint> unit['stereogenic']:
            configure.append(sid)
        else:
            del groups[sid]
            log.append(mc_record(rule, (sid,),
                                'the template puts atom %d in an enhanced-stereo group and the '
                                'patched molecule has no stereogenic unit there, so no group was '
                                'written' % sid))
    if clears or groups or rc:
        with new.edit():
            for sid in rc:
                new.set_parity(<uint32_t> sid, 0)
            for sid in configure:
                new.set_parity(<uint32_t> sid, 1)
            # N8: the group goes with the parity, in the same operation.  A group id left behind on an
            # atom whose parity is gone rejoins that atom to a mixture the next pass to set a parity
            # there has no reason to believe it belongs to.
            for sid in clears:
                new.set_stereo_group(<uint32_t> sid, 0)
            for sid in sorted(groups):
                pair = <tuple> groups[sid]
                new.set_stereo_group(<uint32_t> sid, <int> pair[0], <int> pair[1])

    for sid in dropped:
        log.append(mc_record(rule, (sid,),
                             'atom %d carried a configured parity and the patch changed its '
                             'neighbourhood past what the arena could re-base the sign against, so '
                             'the parity was dropped%s' % (sid, ' and its stereo group cleared with '
                                                           'it' if sid in clears else '')))
    for sid in rc:
        log.append(mc_record(rule, (sid,),
                             'atom %d holds the configuration of a stereo unit the patch wrote part '
                             'of, so the reaction happened at that unit and the configuration was '
                             'dropped%s.  A template whose reaction carries it through states `@=` on '
                             'the atom, and one that turns it over states `@~`'
                             % (sid, ' and its stereo group cleared with it' if sid in clears else '')))

    # N4: the check half of the product side, on the finished molecule.  Last of the passes and
    # before the split, because every primitive it reads is a property of the whole patched graph.
    cdef tuple refused = smk_post_filter(t, new, prod_to_new, alive)
    if refused is not None:
        log.append(mc_record(rule, <tuple> refused[0], <str> refused[1]))
        return None

    # THE MAPPING THIS FILE IMPOSES, before the split so both sides read one table.
    cdef tuple numbering = smk_numbering(origin, alive, touched)
    cdef dict numbers = <dict> numbering[0]
    if <Py_ssize_t> numbering[1]:
        log.append(mc_record(rule, (), 'the reaction pairs more than %d atoms, which is the ceiling '
                                        'on a map number, so %d of them came back unmapped'
                               % (MAP_NUMBER_MAX, <Py_ssize_t> numbering[1])))
    cdef MoleculeContainer numbered = smk_write_numbers(new, numbers)

    # ONLY TOUCHED INPUTS, AND NO COMPONENT DISAPPEARS.  A component of a touched input is a product
    # whether or not the reaction centre reached it, so a salt's counter-ion comes out as its own
    # product molecule -- that is a property of this filter rather than a check bolted onto it.  An
    # input the template never reached is on neither side of the reaction at all.  `origin` is keyed on
    # ATOM ids, which `numbered` shares with `new`; a map number is not one of them.
    cdef list products = []
    if alive:
        for part in numbered.split():
            for sid in (<MoleculeContainer> part).atom_numbers:
                if origin.get(sid) is None or origin[sid] in touched:
                    products.append(part)
                    break
    # THE SAME NUMBERING, RESTATED IN EACH INPUT'S OWN ATOM IDS.  A product atom IS the work atom its
    # reactant partner is, so the pairing is already in `numbers`; `origin_n` is only the change of
    # coordinates.  Copies, because one snapshot is shared by every candidate of the enumeration.
    cdef dict per_input = {}
    cdef list reactants = []
    for key in sorted(touched):
        per_input[key] = {}
    for sid in sorted(numbers):
        (<dict> per_input[<int> origin[sid]])[origin_n[sid]] = numbers[sid]
    for key in sorted(touched):
        reactants.append(smk_write_numbers(<MoleculeContainer> snapshots[<int> key],
                                           <dict> per_input[key]))

    # N9: the dedupe key is STRUCTURE.  `canonical_bytes` is the identity of a stereo-bearing molecule
    # and a SMILES string is not, so nothing here dedupes on `str(r)`.
    cdef list identities = []
    for part in products:
        identities.append((<MoleculeContainer> part).canonical_bytes)
    identities.sort()
    # `prod_to_new` is the only place the template's own numbering and the built atom ids are both in
    # scope; a caller that must edit "the atom the product side called :1" has no other way back.
    # Keyed OUT by map number, not by product sid: the sid is private to this template's parse.
    cdef dict where = {}
    for sid in t.product_map_numbers:
        where[t.product_map_numbers[sid]] = prod_to_new[sid]
    return (_SMK_REACTION(reactants=tuple(reactants), products=tuple(products)),
            (tuple(sorted(touched)), tuple(identities)), where)


def smk_apply(ReactionTemplate t, tuple molecules, bint automorphism_filter, object log, bint report):
    """Validate eagerly, then hand back the generator.

    The split is the point: a bad argument raises from the CALL, not from the first `next()`, so a
    caller who wrote `template(some_string)` finds out where they wrote it.
    """
    smk_lazy_imports()
    cdef object m
    for m in molecules:
        if not isinstance(m, MoleculeContainer):
            if hasattr(m, '__iter__'):
                raise TypeError('a template takes its molecules as separate arguments; unpack the %s '
                                'at the call: template(*molecules)' % type(m).__name__)
            raise TypeError('a template applies to molecules; this one was handed a %s' % type(m).__name__)
    if not molecules:
        raise ValueError('a template needs at least one molecule to apply to')
    return smk_enumerate(t, list(molecules), automorphism_filter, log, report)


def smk_enumerate(ReactionTemplate t, list inputs, bint automorphism_filter, object log, bint report):
    """Every distinct outcome of one template over one set of inputs.

    The inputs are unioned into one working container, which is why an intramolecular template needs
    no separate entry point: `A.B>>C` finds its two fragments wherever they are.  `union` carries both
    sides' stereo, so nothing is racemised on the way in, and both sides' map numbers, which is why the
    reaction's mapping is imposed and not inherited -- two inputs each numbered from 1 would arrive with
    every low number used twice.
    """
    cdef str rule = smk_rule_id(t)
    cdef MoleculeContainer work = (<MoleculeContainer> inputs[0]).copy()
    cdef list snapshots = [(<MoleculeContainer> inputs[0]).copy()]
    cdef dict origin = {}
    cdef dict origin_n = {}
    cdef set before
    cdef list added
    cdef Py_ssize_t k
    cdef object sid, other, mapping, made, key, exc
    for sid in work.atom_numbers:
        origin[sid] = 0
        origin_n[sid] = sid
    for k in range(1, len(inputs)):
        before = set(work.atom_numbers)
        work = work.union(<MoleculeContainer> inputs[k])
        added = []
        for sid in work.atom_numbers:
            if sid not in before:
                origin[sid] = <int> k
                added.append(sid)
        added.sort()
        # THE WAY BACK TO AN INPUT'S OWN ATOM IDS.  `union` appends `other`'s atoms in slot order with
        # fresh ascending ids and never reuses one, so zipping the two sequences pairs them.  Needed
        # because a reactant snapshot keeps its own ids while `work` renumbers every input after the
        # first.
        for sid, other in zip(added, (<MoleculeContainer> inputs[k]).atom_numbers):
            origin_n[sid] = other
        snapshots.append((<MoleculeContainer> inputs[k]).copy())
    cdef tuple snaps = tuple(snapshots)

    # the working container's bonds, once, low-first: `smk_one` needs the PRE-EDIT answer to
    # "does this bond exist" for every candidate and the graph it asks about never changes
    cdef set work_bonds = set()
    for sid in work.atom_numbers:
        for other in work.neighbors_of(<uint32_t> sid):
            if sid < other:
                work_bonds.add((sid, other))

    cdef set seen_keys = set()
    cdef Py_ssize_t start
    for mapping in t.reactants.get_mapping(work, automorphism_filter):
        # THE OUTCOME CARRIES WHAT ITS OWN PATCH REPORTED, on `rxn.log` under stage `react`, whether or
        # not the caller passed a `log=`.  A candidate that produced nothing has no container to write
        # to, so its line is only on the caller's list -- and a duplicate outcome is not yielded, so
        # nothing is stamped twice.
        start = len(log)
        # N11: THE WHOLE CANDIDATE IS INSIDE THE GUARD, dedupe key included.  A guard around the patch
        # alone lets an exception from the dedupe escape the generator, and that takes the tail of the
        # enumeration with it.
        try:
            made = smk_one(t, work, work_bonds, origin, origin_n, snaps, <dict> mapping, log)
        except Exception as exc:
            log.append(mc_record(rule, tuple(sorted((<dict> mapping).values())),
                                 'the patch raised %s: %s -- this match produced nothing and the '
                                 'enumeration continues' % (type(exc).__name__, exc)))
            continue
        if made is None:
            continue
        key = (<tuple> made)[1]
        if key in seen_keys:
            continue
        seen_keys.add(key)
        if len(log) > start:      # `_smiles_read.pxi:smi_one` states why the touch is guarded
            (<object> (<tuple> made)[0]).log.absorb('react', log[start:])
        if report:
            yield ((<tuple> made)[0], (<tuple> made)[2])
        else:
            yield (<tuple> made)[0]
