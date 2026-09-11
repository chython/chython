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
# THE CHEMISTRY VALENCE MODEL.  THERE IS A SECOND VALENCE MODEL AND MERGING THEM IS A BUG.
#
# The rules are DATA and they live in `chython/core/valence_rules.tsv`, which is the authority;
# the tables below are compiled from it by `chython/core/test/gen_valence_rules.py compile` and
# committed.  Read that file's header first -- it says what a row means.  This file is the query
# side: a key lookup, a multiset comparison against the neighbourhood, and three questions asked
# of the result.
#
# THE COLLECTION IS A POSITIVE LIST, NOT A THEORY.
#
# It exists to catch bad input.  It is neither complete nor ideal, it grows by adding rows with
# evidence behind them, and NOTHING IN CHYTHON MAY REJECT A STRUCTURE BECAUSE OF A VERDICT FROM
# THIS FILE.  A structure whose valence violates every rule is stored, reported and left to the
# repair pipeline; a reader that refused it would lose a real record, and a reader that silently
# "fixed" it would lose the evidence.  That is why the answer is a three-state verdict rather
# than a bool -- a bool invites `if not legal: raise`, and the third state makes the honest
# behaviour the easy one.
#
# THE OTHER MODEL.
#
# This file answers a question about a MOLECULE: is this a valence state chemistry is known to
# allow, and how many hydrogens does it come with?  `smv_default_h` (in the SMILES layer) answers
# a question about a NOTATION: what hydrogen count does a bracketless atom in a SMILES string
# imply, and where are brackets therefore mandatory?  Its authority is the Daylight/OpenSMILES
# specification.  The two disagree, on purpose:
#
#   * A bare `S` with six single-bonded carbons is legal SMILES implying zero hydrogens -- the
#     notation permits it because sulfur's wide valence set includes 6.  This file answers "no
#     rule" for that atom, while answering 0 for `CS(=O)(=O)C`.  Same element, same charge, same
#     valence, different neighbourhood: the chemistry model takes an environment and the notation
#     model has nowhere to put one, because a string's syntax cannot depend on what its atoms are
#     bonded to.
#   * Neutral five-valent nitrogen is spellable without brackets and has no row here at all --
#     the collection wants the charge-separated form.
#
# Merging the two therefore fails in both directions: either the SMILES reader starts rejecting
# strings RDKit emits, or MDL read starts accepting states this collection calls violations.
# `test_the_two_valence_models_answer_differently` pins one input where the answers must differ,
# in this suite AND in the SMILES writer's, so a future merge fails a test rather than passing a
# review.
#
# CONSUMERS -- named, because the wrong caller is the hazard.
#
# MDL read (for the hydrogen count only), the InChI bridge, and standardization's valence check.
# All of them genuinely ask a chemistry question.
#
# NOT the SMILES writer, and not any other writer.  "Strict out" means spec-conformant SYNTAX,
# never chemically validated CONTENT: a writer spells what is stored, including a molecule whose
# valence is impossible, because refusing to write is refusing to let a user see what they are
# holding.  Not `canonical_order` either, and not isomorphism: neither asks a chemistry question,
# and a valence check inside them would make an illegal molecule silently unmatchable instead of
# visibly wrong.
#
# THREE QUESTIONS, THREE ENTRY POINTS.
#
# "How many hydrogens does this atom get" (`val_implicit_h`), "what does the collection make of
# this state" (`val_check`), and "does the collection describe this valence at all"
# (`val_has_rules`).  The first two are NOT each other's inverse: several rows can sit at one key,
# the hydrogen count is the FIRST matching row's, and the verdict accepts ANY matching row's.  So an
# atom can legally carry a count that `val_implicit_h` would not have chosen, and collapsing the two
# would make one of those behaviours unreachable.
#
# TWO QUESTIONS, TWO TABLES, AND ONLY ONE OF THEM IS HOT.
#
# One collection, but the hydrogen question runs on every atom of every parsed molecule and the
# check runs when somebody asks.  Measured: of 1036 rows, 1031 answer the hydrogen question with no
# reference to a neighbourhood, and 583 exist only to be checked against.  Charging the parse path a
# binary search over 576 keys plus a multiset compare, to reach a row whose environment is `*` in
# 99.4% of cases, is paying the checker's price on the parser's traffic.
#
# So there are two generated artifacts and they are queried by different code:
#
#   VAL_H_PAT / VAL_H_ROW    hot.  `(z, charge, radical) -> pattern`, `[pattern][bonds] -> count`.
#                            One indexed load and a branch.  65 patterns, about 2.7 KB total.
#   VAL_KEY / VAL_H / VAL_ENV cold.  Every row, in scan order, with the environments interned.
#                            `val_check` and `val_has_rules` only.
#
# The hot table is a PROJECTION of the cold one and never a summary: `-1` means no row covers the
# state (not zero hydrogens), and `-2` means an environment decides, which sends that atom -- 0.6%
# of them, mostly sulfones -- through the cold table.  Dropping the `-2` state would answer
# "no rule" for every sulfone, nitro group and perchlorate, which is why the two tables are not
# simply "the 43 hydrogen rules" and "the rest": that split does not survive contact with the data.
# `test_the_dense_table_is_the_full_scan` sweeps the whole domain and fails if they ever disagree.
#
# "NO RULE" IS NOT ZERO, AND "NO RULE" IS NOT "WRONG".
#
# `val_implicit_h` returns VAL_NO_RULE and its Python wrapper returns None.  Returning 0 would
# conflate "this atom has no hydrogens" with "the collection has nothing to say about this atom",
# and the arena cannot store the distinction -- `atom_t.hydrogens` is two nibbles with no unset
# state -- so it has to survive in the answer rather than in the molecule.
#
# `val_check` splits the miss further, and the boundary is worth stating precisely: a state is a
# VIOLATION when the collection describes this element in this charge and radical state and no
# row accepts what you have; it is UNKNOWN when the collection says nothing about that element,
# charge and radical state at all.  Pentavalent neutral carbon is a violation -- neutral carbon is
# thoroughly described, so the absence of a row is a claim.  A radical lanthanide is unknown --
# nobody wrote anything about it, and a library should not invent a verdict.  UNKNOWN is the
# data-driven-development hook: `gen_valence_rules.py coverage` counts exactly those atoms, by
# element, and a `mined:<corpus>` row is how the count goes down.
#
# ROW ORDER IS OBSERVABLE, WHICH IS WHY THE FILE IS SORTED AND NOT A SET.
#
# Fifteen keys in the shipped collection hold rows that would answer the hydrogen question
# differently in a different order.  Eleven are bare atoms where a `common` row and a `curated`
# row disagree (`[C]` is methane's 4 hydrogens because the common row is scanned first, not the
# atomic-carbon row's 0); the other four are bonded states of phosphorus and sulfur.  See
# `test_row_order_within_a_key_is_observable`.
#
# ONE INCONSISTENCY IN THE SHIPPED COLLECTION, PRESERVED RATHER THAN SMOOTHED.
#
# The fifteenth is not a bare atom: phosphorus at three bonds with `-O =O` should by the curated row
# carry two hydrogens (phosphorous acid), but the common valence 3 row sits at the same key with zero
# and is scanned first, so the curated row's two-hydrogen state is unreachable.  The `common` and
# `curated` rows are derived from chython 2 and a test re-derives them, so changing this means
# deleting or re-keying a row in the TSV -- not a change to this file -- and it must be argued on
# chemistry.
#
# ELEMENT DATA HANGS OFF AN ATOMIC NUMBER, NOT OFF A PER-ELEMENT TYPE.
#
# There is no 118-element class hierarchy anywhere in the core.  Isotope masses and abundances belong
# to the InChI epic.
#
# AROMATIC BONDS ARE NOT THIS FILE'S BUSINESS.
#
# There are no aromatic rows -- the collection is about localised bond orders -- so a caller
# holding an order-4 bond must decide its own policy, and `valence_implicit_h` refuses order 4
# rather than guessing.  Hard-coding benzene-shaped neutral carbon and answering None for everything
# else aromatic is arithmetic a caller must own where the aromatic system is understood, not here.
#
DEF VAL_NO_RULE = -1     # `val_implicit_h`: no row covers this state.  Distinct from 0 hydrogens
DEF VAL_ANY_H = -1       # `val_scan`: report the first matching row's count, do not filter by it
DEF VAL_ENV_SHIFT = 8    # a VAL_ENV entry is (bond order << VAL_ENV_SHIFT) | atomic number
DEF VAL_ORDER_MAX = 3    # the highest bond order any row's environment mentions

# The dense hydrogen table's domain, and the third sentinel in it.  Mirrored in the generator as
# H_*, which fails the compile if a new row leaves these extents -- they are the collection's
# measured span and not a guess.  A state outside them has no row, so it needs no load.
#
# THE CHARGE DOMAIN HERE IS THE RULES' DOMAIN AND NOT THE ARENA'S STORAGE RANGE, and the two are
# close enough to mislead whoever keys a new row.  `atom_t.charge` stores CHARGE_MIN = -4 to
# CHARGE_MAX = 8 (`_molecule_arena.pxi`); the collection describes -4 to +4 and NOTHING above +4.
# Measured over `valence_rules()`, all 1036 rows across all 118 elements: -4:1, -3:10, -2:28,
# -1:60, 0:791, +1:69, +2:30, +3:39, +4:8.
#
# So `[S+6]` is STORABLE AND UNDESCRIBED, and it gets VAL_NO_RULE / VAL_UNKNOWN -- never a
# confident 0.  Same ruling as the sparse metal valence states: an absent row is a GAP, not a
# violation, it gets no opinion, and it shows up in the coverage report.  Answering 0 for a state
# nobody looked at invents chemistry and destroys the distinction the coverage report exists to
# preserve.  Do not widen the span to 13 by reflex to match the storage range: the four slots above
# +4 would be structurally empty, and the guard in `val_implicit_h` is what keeps them from being
# indexed at all.  `test_the_charge_domain_is_the_rules_span_not_the_storage_range` pins it.
DEF VAL_H_Z_MAX = 118
DEF VAL_H_CHARGE_MIN = -4
DEF VAL_H_CHARGE_MAX = 4
DEF VAL_H_CHARGE_SPAN = 9          # VAL_H_CHARGE_MAX - VAL_H_CHARGE_MIN + 1
DEF VAL_H_BONDS_MAX = 8            # the widest bond-order sum any row states
DEF VAL_H_STRIDE = 9               # VAL_H_BONDS_MAX + 1
DEF VAL_H_CONSULT = -2             # this state's hydrogen count depends on the neighbourhood

# The verdict.  A three-way state and NOT a bool: `if val_check(...)` is always a bug in C,
# because two of the three states are truthy.  Compare against these names.  The numbers are
# indices into VAL_VERDICT_NAMES and carry no ordering: there is no "worse" among them.
DEF VAL_UNKNOWN = 0      # the collection says nothing about this element, charge and radical
DEF VAL_VALID = 1        # a row accepts this state, hydrogen count included
DEF VAL_VIOLATION = 2    # the collection describes this element here, and no row accepts it

# The key domain, declared once (RULES.md §6).  `gen_valence_rules.py` emits VAL_KEY packed this
# way and nothing checks the two copies by inspection -- if they disagreed, the exhaustive sweep
# in test_valence.py would fail to find any row through its own key.
DEF VAL_KEY_Z_SHIFT = 16
DEF VAL_KEY_CHARGE_SHIFT = 12
DEF VAL_KEY_RADICAL_SHIFT = 11
DEF VAL_CHARGE_BIAS = 8            # so an unsigned sort of VAL_KEY is a signed sort of the charge
DEF VAL_BONDS_MAX = 2047           # the 11 bits below the radical bit
DEF VAL_NO_KEY = 0x7FFFFFFF        # an unpackable state.  Larger than any real key (the
                                   # highest is 118 << 16) and small enough to stay a C int, which
                                   # 0xFFFFFFFF is not -- a DEF that overflows int becomes a Python
                                   # object and drags the whole comparison out of nogil


# --- BEGIN GENERATED TABLES: python chython/core/test/gen_valence_rules.py compile ---
# Compiled from chython/core/valence_rules.tsv, which is the authority.  Do not edit by
# hand -- run the command in the marker above.  The TSV is in scan order and so is this,
# so the k-th entry here is the k-th row there.
#
# 1036 rules over 576 keys on 118 elements; 219 distinct environments, 732 entries.
# Provenance: 256 common, 780 curated.
#
# The hot artifact is separate and is a projection of the same rows: 65 distinct hydrogen
# patterns behind 2142 states.  See `hydrogen_tables`.
cdef extern from *:
    """
    /* sorted: (z << 16) | ((charge + 8) << 12) | (radical << 11) | bonds */
    static const unsigned int VAL_KEY[576] = {
    94208, 98305, 100352, 102400, 163840, 229376, 229377, 233472,
    294912, 294914, 303104, 356352, 356353, 356354, 356355, 356356,
    360448, 360449, 360450, 360451, 362496, 362497, 362498, 421888,
    421889, 421890, 421891, 425984, 425985, 425986, 425987, 425988,
    428032, 428033, 428034, 428035, 430080, 430081, 430082, 430083,
    487424, 487425, 487426, 491520, 491521, 491522, 491523, 493568,
    493569, 493570, 495616, 495617, 495618, 495619, 495620, 548864,
    552960, 552961, 557056, 557057, 557058, 559104, 559105, 561152,
    561153, 561154, 561155, 618496, 622592, 622593, 688128, 753664,
    753665, 757760, 819200, 819202, 823297, 827392, 872454, 880640,
    880641, 880642, 880643, 880644, 884736, 884737, 884738, 884739,
    888832, 888833, 888834, 892928, 892929, 897024, 942086, 950272,
    950273, 950274, 950275, 950276, 1011712, 1011713, 1011714, 1011718,
    1015808, 1015809, 1015810, 1015811, 1015812, 1015813, 1017856, 1017857,
    1017858, 1017859, 1017860, 1019904, 1019905, 1019906, 1019907, 1019908,
    1073152, 1077248, 1077249, 1081344, 1081345, 1081346, 1081348, 1081350,
    1083392, 1083393, 1083394, 1083395, 1085443, 1085445, 1142784, 1142786,
    1146880, 1146881, 1146883, 1146885, 1146887, 1212416, 1277952, 1277953,
    1282048, 1343488, 1343490, 1351680, 1396742, 1409024, 1409027, 1421312,
    1466374, 1474560, 1474562, 1474563, 1474564, 1482754, 1490944, 1540096,
    1540098, 1540099, 1540100, 1540101, 1548288, 1548290, 1552384, 1605632,
    1605634, 1605635, 1605636, 1605638, 1613824, 1617920, 1671168, 1671170,
    1671171, 1671172, 1671174, 1671175, 1679360, 1683456, 1736704, 1736706,
    1736707, 1744896, 1748992, 1785862, 1789957, 1789958, 1794052, 1794054,
    1798147, 1802240, 1802241, 1802242, 1802243, 1810432, 1810433, 1814528,
    1867776, 1867778, 1867779, 1871873, 1875968, 1921026, 1929218, 1933312,
    1933313, 1933314, 1937408, 1941504, 1990660, 1998848, 1998850, 2002945,
    2007040, 2060292, 2064384, 2064385, 2064387, 2076672, 2121734, 2129920,
    2129921, 2129922, 2129923, 2129924, 2191366, 2195456, 2195459, 2195461,
    2199552, 2199553, 2199554, 2199555, 2199556, 2252800, 2256896, 2256897,
    2260992, 2260993, 2260994, 2260996, 2260998, 2265091, 2322432, 2322434,
    2326528, 2326529, 2326531, 2326533, 2326535, 2392064, 2457600, 2457601,
    2461696, 2523136, 2523138, 2531328, 2588672, 2588675, 2600960, 2654208,
    2654210, 2654211, 2654212, 2719744, 2719746, 2719747, 2719748, 2719749,
    2785280, 2785284, 2785285, 2785286, 2850816, 2850820, 2916352, 2916356,
    2916357, 2916359, 2916360, 2969606, 2977796, 2981888, 2981889, 2981890,
    2981891, 2981892, 2981894, 3039236, 3047424, 3047426, 3051521, 3055616,
    3108866, 3112960, 3112961, 3112962, 3117056, 3170308, 3178496, 3178498,
    3186688, 3244032, 3244033, 3244035, 3256320, 3301382, 3309568, 3309570,
    3309571, 3309572, 3313667, 3317760, 3371014, 3375104, 3375107, 3375109,
    3379200, 3379201, 3379202, 3379203, 3379204, 3436549, 3440640, 3440641,
    3440642, 3440644, 3440646, 3444739, 3502080, 3502082, 3506176, 3506177,
    3506179, 3506181, 3506183, 3510274, 3571712, 3571714, 3571716, 3571718,
    3571720, 3637248, 3637249, 3641344, 3702784, 3702786, 3710976, 3768320,
    3768323, 3780608, 3833856, 3833859, 3833860, 3846144, 3899392, 3899395,
    3899396, 3911680, 3964928, 3964930, 3964931, 3964932, 3977216, 4030464,
    4030467, 4042752, 4096000, 4096002, 4096003, 4108288, 4161536, 4161538,
    4161539, 4173824, 4227072, 4227075, 4239360, 4292608, 4292611, 4292612,
    4304896, 4358144, 4358147, 4358148, 4370432, 4423680, 4423682, 4423683,
    4435968, 4489216, 4489219, 4501504, 4554752, 4554754, 4554755, 4567040,
    4620288, 4620290, 4620291, 4632576, 4685824, 4685827, 4698112, 4751360,
    4751361, 4751362, 4751363, 4751364, 4816896, 4816901, 4882432, 4882438,
    4947968, 5013504, 5013510, 5013512, 5066758, 5079040, 5079041, 5079042,
    5079043, 5079044, 5079045, 5079046, 5144576, 5144578, 5144580, 5144582,
    5152768, 5206020, 5210112, 5210115, 5214208, 5222400, 5275648, 5275650,
    5283840, 5328902, 5341184, 5341185, 5345280, 5345282, 5353472, 5398532,
    5406720, 5406722, 5406724, 5414912, 5472256, 5472257, 5472258, 5472259,
    5472260, 5472261, 5484544, 5537792, 5537794, 5537796, 5537798, 5599232,
    5603328, 5603329, 5603333, 5607424, 5668864, 5668866, 5672961, 5734400,
    5734401, 5738496, 5799936, 5799938, 5808128, 5865472, 5865475, 5877760,
    5931008, 5931010, 5931011, 5931012, 5947392, 5996544, 5996546, 5996547,
    5996548, 5996549, 6012928, 6062080, 6062083, 6062084, 6062085, 6062086,
    6070276, 6074368, 6078464, 6127616, 6127618, 6127619, 6127620, 6127621,
    6127622, 6127623, 6131716, 6135812, 6139904, 6144000, 6193152, 6193154,
    6193155, 6193156, 6193157, 6193158, 6197252, 6201348, 6205440, 6209536,
    6258688, 6258690, 6258691, 6258692, 6270976, 6324224, 6324226, 6324227,
    6324228, 6389760, 6389762, 6389763, 6389764, 6402048, 6406144, 6455296,
    6455298, 6455299, 6455300, 6467584, 6520832, 6520834, 6520835, 6533120,
    6586368, 6586370, 6586371, 6598656, 6651904, 6651906, 6651907, 6664192,
    6717440, 6717442, 6717443, 6725632, 6782976, 6782979, 6795264, 6848512,
    6848516, 6864896, 6914048, 6979584, 7045120, 7110656, 7176192, 7241728,
    7307264, 7372800, 7438336, 7503872, 7569408, 7634944, 7700480, 7766016
    };
    /* first rule of the k-th key */
    static const unsigned short VAL_KEY_OFF[576] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24,
    25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37,
    38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
    62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 78, 82, 83, 84, 85, 86, 87, 88,
    89, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
    103, 104, 105, 106, 107, 108, 109, 110, 112, 114, 115, 116,
    118, 121, 125, 127, 129, 131, 132, 133, 134, 135, 136, 137,
    138, 139, 140, 141, 143, 144, 145, 208, 255, 257, 259, 260,
    261, 268, 270, 271, 272, 273, 274, 276, 280, 281, 282, 283,
    284, 285, 286, 287, 288, 289, 290, 291, 292, 301, 302, 309,
    318, 319, 321, 322, 323, 324, 329, 334, 340, 341, 342, 343,
    344, 345, 346, 350, 354, 355, 356, 357, 359, 361, 363, 364,
    365, 366, 367, 368, 369, 370, 371, 372, 373, 375, 376, 381,
    382, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395,
    396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407,
    408, 409, 414, 415, 418, 419, 420, 421, 423, 424, 425, 426,
    427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438,
    439, 441, 442, 443, 458, 465, 467, 468, 473, 474, 475, 478,
    482, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
    495, 506, 513, 514, 515, 516, 518, 523, 529, 530, 535, 538,
    542, 543, 545, 546, 547, 548, 549, 550, 552, 553, 554, 556,
    560, 561, 562, 563, 566, 567, 568, 569, 570, 573, 574, 575,
    576, 577, 578, 579, 580, 581, 582, 586, 587, 588, 589, 590,
    594, 595, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606,
    607, 608, 609, 611, 612, 613, 625, 629, 631, 632, 635, 636,
    637, 654, 661, 668, 669, 670, 671, 672, 676, 679, 680, 681,
    682, 683, 684, 685, 686, 687, 688, 689, 690, 694, 695, 696,
    697, 699, 700, 701, 712, 713, 714, 715, 716, 717, 718, 719,
    730, 731, 732, 733, 744, 745, 746, 747, 748, 749, 750, 751,
    753, 754, 755, 756, 758, 759, 760, 771, 772, 773, 774, 775,
    776, 777, 788, 789, 790, 791, 802, 803, 804, 805, 806, 807,
    808, 809, 810, 816, 817, 818, 824, 825, 830, 831, 832, 834,
    836, 837, 838, 843, 848, 849, 850, 851, 852, 853, 854, 857,
    859, 860, 861, 862, 865, 866, 867, 868, 869, 870, 871, 872,
    873, 874, 875, 876, 878, 879, 880, 886, 887, 888, 890, 897,
    898, 900, 903, 904, 905, 906, 910, 912, 913, 914, 915, 916,
    917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928,
    929, 930, 933, 935, 936, 937, 938, 939, 940, 941, 942, 943,
    944, 945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 955,
    956, 957, 958, 959, 960, 961, 962, 963, 971, 972, 973, 974,
    975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 987,
    988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999,
    1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011,
    1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023,
    1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035
    };
    /* how many rules that key has, in scan order */
    static const unsigned char VAL_KEY_LEN[576] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1,
    1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
    1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 3, 4, 2, 2,
    2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 63, 47,
    2, 2, 1, 1, 7, 2, 1, 1, 1, 1, 2, 4, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 9, 1, 7, 9, 1, 2, 1, 1,
    1, 5, 5, 6, 1, 1, 1, 1, 1, 1, 4, 4, 1, 1, 1, 2,
    2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 5, 1,
    3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 3, 1, 1, 1, 2,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 1, 1, 15, 7, 2, 1, 5, 1, 1, 3, 4, 2, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 11, 7, 1, 1, 1, 2, 5, 6,
    1, 5, 3, 4, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 4,
    1, 1, 1, 3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 4, 1, 1, 1, 1, 4, 1, 2, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 2, 1, 1, 12, 4, 2, 1, 3, 1, 1,
    17, 7, 7, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 4, 1, 1, 1, 2, 1, 1, 11, 1, 1, 1, 1,
    1, 1, 1, 11, 1, 1, 1, 11, 1, 1, 1, 1, 1, 1, 1, 2,
    1, 1, 1, 2, 1, 1, 11, 1, 1, 1, 1, 1, 1, 11, 1, 1,
    1, 11, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 5,
    1, 1, 2, 2, 1, 1, 5, 5, 1, 1, 1, 1, 1, 1, 3, 2,
    1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
    1, 1, 6, 1, 1, 2, 7, 1, 2, 3, 1, 1, 1, 4, 2, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };
    /* implicit hydrogen count of the k-th rule */
    static const unsigned char VAL_H[1036] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 2, 1, 0, 3, 0, 2, 1, 0, 2, 1, 0,
    3, 2, 1, 0, 4, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 3, 2, 1,
    0, 2, 1, 0, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 1, 0, 3, 2, 1, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 2, 1, 0, 0, 3, 2, 1, 0, 2, 1,
    0, 1, 0, 0, 0, 4, 0, 3, 2, 1, 0, 2, 1, 0, 0, 0, 3, 0, 2, 1, 0, 2, 1, 1,
    1, 0, 0, 0, 0, 2, 4, 1, 3, 0, 2, 1, 0, 4, 3, 2, 1, 0, 0, 1, 0, 2, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 2, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 3, 2, 1, 0, 0, 0, 0, 0, 4,
    3, 2, 1, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 4, 3, 2, 1, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0
    };
    /* the k-th rule's environment, interned by content */
    static const unsigned short VAL_ENV_OFF[1036] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 4, 10, 0, 0, 0, 0, 0, 16, 18, 21,
    16, 0, 18, 21, 16, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 24, 26, 28, 30, 32, 34, 18, 36, 39, 42, 45,
    48, 51, 54, 57, 60, 63, 66, 69, 72, 21, 75, 78,
    81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114,
    117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150,
    153, 156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196,
    200, 204, 208, 212, 216, 219, 222, 226, 230, 234, 238, 242,
    246, 250, 254, 258, 262, 266, 270, 274, 278, 282, 286, 290,
    294, 298, 302, 306, 310, 314, 318, 322, 326, 330, 334, 338,
    342, 346, 350, 354, 358, 362, 366, 370, 374, 378, 382, 386,
    4, 390, 396, 0, 0, 0, 0, 0, 0, 402, 404, 406,
    408, 411, 414, 417, 420, 424, 0, 428, 0, 0, 16, 430,
    433, 436, 441, 445, 448, 0, 0, 0, 0, 0, 0, 0,
    4, 0, 0, 0, 4, 452, 458, 464, 470, 476, 222, 216,
    481, 0, 486, 488, 490, 492, 494, 496, 498, 430, 499, 502,
    505, 508, 16, 511, 514, 516, 0, 498, 496, 0, 0, 0,
    430, 499, 502, 16, 517, 168, 519, 523, 24, 66, 436, 433,
    527, 530, 441, 534, 0, 498, 0, 0, 0, 0, 24, 18,
    168, 208, 216, 222, 242, 302, 0, 0, 0, 0, 498, 16,
    508, 24, 24, 222, 448, 0, 0, 0, 0, 0, 0, 0,
    470, 530, 538, 470, 168, 519, 523, 543, 547, 4, 430, 499,
    508, 0, 551, 0, 0, 0, 552, 0, 0, 0, 16, 0,
    0, 553, 488, 0, 0, 0, 0, 0, 547, 0, 0, 0,
    0, 555, 168, 519, 523, 543, 0, 3, 2, 559, 0, 0,
    4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24,
    560, 28, 18, 66, 562, 21, 60, 51, 117, 138, 196, 565,
    569, 573, 222, 302, 274, 230, 330, 326, 4, 408, 402, 0,
    577, 490, 579, 488, 492, 0, 0, 16, 508, 430, 433, 436,
    441, 445, 448, 581, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 486, 488, 490, 492, 494, 585, 587, 589, 591,
    593, 595, 430, 499, 502, 505, 511, 514, 516, 0, 0, 498,
    16, 516, 168, 519, 523, 24, 560, 436, 538, 441, 534, 433,
    530, 0, 168, 519, 523, 24, 560, 436, 538, 596, 216, 222,
    4, 481, 0, 24, 560, 0, 601, 604, 448, 608, 452, 470,
    523, 0, 551, 3, 498, 496, 612, 553, 0, 0, 4, 547,
    168, 519, 0, 0, 0, 0, 488, 496, 553, 0, 0, 486,
    0, 547, 0, 0, 0, 0, 3, 2, 559, 1, 0, 0,
    470, 0, 498, 496, 612, 488, 408, 0, 408, 408, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0, 613, 0, 0, 0,
    0, 24, 18, 51, 21, 618, 196, 569, 565, 176, 208, 622,
    626, 222, 470, 4, 390, 408, 402, 0, 492, 577, 488, 0,
    0, 16, 630, 402, 406, 508, 430, 499, 502, 632, 414, 635,
    638, 641, 644, 647, 408, 417, 433, 436, 441, 445, 650, 653,
    657, 448, 662, 667, 673, 680, 686, 581, 591, 0, 486, 168,
    4, 216, 298, 481, 608, 691, 696, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 547, 168, 519, 24, 0, 0,
    0, 168, 24, 0, 0, 498, 496, 486, 488, 490, 492, 702,
    704, 706, 591, 494, 0, 168, 0, 0, 0, 0, 0, 496,
    486, 488, 490, 492, 702, 704, 706, 591, 494, 498, 0, 0,
    0, 496, 486, 488, 490, 492, 702, 704, 706, 591, 494, 498,
    0, 0, 0, 0, 0, 0, 0, 168, 24, 0, 0, 0,
    168, 24, 0, 0, 496, 486, 488, 490, 492, 702, 704, 706,
    591, 494, 498, 0, 0, 0, 0, 0, 0, 496, 486, 488,
    490, 492, 702, 704, 706, 591, 494, 498, 0, 0, 0, 496,
    486, 488, 490, 492, 702, 704, 706, 591, 494, 498, 0, 0,
    0, 0, 0, 0, 3, 490, 499, 502, 505, 511, 514, 516,
    0, 0, 436, 538, 441, 534, 433, 530, 0, 216, 222, 4,
    452, 481, 0, 0, 216, 222, 608, 708, 452, 0, 1, 713,
    3, 2, 559, 488, 490, 492, 553, 612, 0, 0, 436, 4,
    0, 0, 714, 718, 24, 4, 216, 0, 519, 0, 499, 722,
    502, 0, 0, 0, 0, 0, 452, 0, 0, 0, 591, 0,
    547, 18, 0, 0, 24, 547, 18, 168, 519, 208, 0, 0,
    3, 2, 488, 490, 492, 553, 612, 725, 727, 0, 519, 24,
    436, 433, 530, 0, 0, 0, 24, 519, 523, 543, 216, 4,
    0, 0, 0, 433, 0, 0, 486, 713, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 490, 492, 494, 502, 505, 0,
    0, 0, 498, 728, 0, 0, 0, 0, 0, 0, 0, 0,
    24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 24,
    0, 0, 0, 727, 612, 731, 498, 488, 490, 492, 494, 0,
    0, 0, 0, 24, 24, 0, 0, 0, 0, 0, 0, 0,
    0, 498, 494, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0
    };
    /* how many neighbours it demands; 0 is `env=*` */
    static const unsigned char VAL_ENV_LEN[1036] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 0, 2, 3, 3,
    2, 0, 3, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 0, 0, 0, 0, 0, 0, 2, 2, 2,
    3, 3, 3, 3, 4, 4, 0, 2, 0, 0, 2, 3, 3, 5, 4, 3, 4, 0, 0, 0, 0, 0, 0, 0,
    6, 0, 0, 0, 6, 6, 6, 6, 6, 5, 4, 3, 5, 0, 2, 2, 2, 2, 2, 2, 1, 3, 3, 3,
    3, 3, 2, 3, 2, 1, 0, 1, 2, 0, 0, 0, 3, 3, 3, 2, 2, 4, 4, 4, 2, 3, 5, 3,
    3, 4, 4, 4, 0, 1, 0, 0, 0, 0, 2, 3, 4, 4, 3, 4, 4, 4, 0, 0, 0, 0, 1, 2,
    3, 2, 2, 4, 4, 0, 0, 0, 0, 0, 0, 0, 6, 4, 5, 6, 4, 4, 4, 4, 4, 6, 3, 3,
    3, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 2, 2, 0, 0, 0, 0, 0, 4, 0, 0, 1,
    0, 4, 4, 4, 4, 4, 0, 1, 1, 1, 0, 0, 6, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 6, 3, 2, 0, 2, 2, 2, 2, 2, 0, 0, 2, 3, 3, 3, 5,
    4, 3, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 1, 3, 3, 3, 3, 3, 2, 1, 0, 0, 1, 2, 1, 4, 4, 4, 2, 2, 5, 5, 4, 4, 3,
    4, 0, 4, 4, 4, 2, 2, 5, 5, 5, 3, 4, 6, 5, 0, 2, 2, 0, 3, 4, 4, 4, 6, 6,
    4, 0, 1, 1, 1, 2, 1, 2, 0, 0, 6, 4, 4, 4, 0, 0, 1, 0, 2, 2, 2, 0, 0, 2,
    0, 4, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 6, 0, 1, 2, 1, 2, 3, 0, 3, 3, 0, 6,
    0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 6, 6, 6, 3, 2, 0, 2, 2, 2, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 5, 4, 3, 3, 4, 5, 4, 5, 6, 7, 6, 5, 4, 2, 0, 2, 4,
    6, 3, 4, 5, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 2, 0, 0,
    0, 4, 2, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 4, 0, 0, 0, 0, 0, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
    0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 4, 2, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 1, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3, 2, 1,
    0, 0, 5, 5, 4, 4, 3, 4, 0, 3, 4, 6, 6, 5, 0, 0, 3, 4, 4, 5, 6, 0, 1, 1,
    1, 1, 1, 2, 2, 2, 2, 1, 0, 0, 5, 6, 0, 0, 4, 4, 2, 6, 3, 0, 4, 0, 3, 3,
    3, 0, 0, 0, 0, 0, 6, 0, 0, 0, 2, 0, 4, 3, 0, 0, 2, 4, 3, 4, 4, 4, 0, 0,
    1, 1, 2, 2, 2, 2, 1, 2, 1, 0, 4, 2, 5, 3, 4, 0, 0, 0, 2, 4, 4, 4, 3, 6,
    0, 0, 0, 3, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 3, 3, 0,
    0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,
    0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0
    };
    /* (bond order << 8) | atomic number */
    static const unsigned short VAL_ENV[732] = {
    262, 264, 291, 273, 265, 265, 265, 265, 265, 265, 262, 262,
    262, 265, 265, 265, 264, 520, 264, 264, 520, 262, 262, 520,
    520, 520, 519, 520, 519, 519, 518, 520, 518, 518, 518, 519,
    264, 272, 520, 263, 264, 520, 264, 265, 520, 264, 273, 520,
    264, 291, 520, 262, 264, 520, 263, 263, 520, 263, 273, 520,
    262, 263, 520, 263, 272, 520, 273, 273, 520, 291, 291, 520,
    272, 272, 520, 262, 265, 520, 262, 273, 520, 262, 291, 520,
    262, 272, 520, 262, 271, 520, 265, 265, 519, 262, 262, 519,
    262, 264, 519, 262, 273, 519, 262, 272, 519, 262, 263, 519,
    263, 263, 519, 263, 264, 519, 264, 264, 519, 262, 262, 518,
    262, 265, 518, 262, 272, 518, 262, 263, 518, 272, 272, 518,
    263, 272, 518, 263, 263, 518, 264, 264, 518, 262, 262, 528,
    262, 264, 528, 264, 264, 528, 262, 263, 528, 262, 272, 528,
    263, 265, 265, 265, 262, 262, 262, 265, 262, 265, 265, 265,
    265, 265, 265, 265, 262, 262, 262, 264, 262, 262, 264, 264,
    262, 264, 264, 264, 263, 264, 264, 264, 262, 263, 264, 264,
    262, 262, 263, 264, 262, 262, 273, 273, 262, 262, 263, 263,
    262, 262, 262, 272, 262, 262, 262, 262, 262, 262, 262, 309,
    520, 520, 520, 518, 520, 520, 264, 264, 520, 520, 263, 264,
    520, 520, 262, 264, 520, 520, 264, 272, 520, 520, 264, 265,
    520, 520, 264, 273, 520, 520, 264, 291, 520, 520, 264, 309,
    520, 520, 263, 263, 520, 520, 262, 263, 520, 520, 263, 272,
    520, 520, 263, 265, 520, 520, 263, 273, 520, 520, 262, 262,
    520, 520, 262, 272, 520, 520, 262, 265, 520, 520, 262, 273,
    520, 520, 262, 291, 520, 520, 262, 309, 520, 520, 265, 265,
    520, 520, 273, 273, 520, 520, 265, 273, 520, 520, 263, 264,
    519, 520, 263, 263, 519, 520, 262, 263, 519, 520, 264, 264,
    519, 520, 262, 264, 519, 520, 262, 262, 519, 520, 262, 273,
    519, 520, 262, 265, 519, 520, 262, 262, 519, 519, 262, 262,
    518, 520, 262, 264, 518, 520, 263, 264, 518, 520, 264, 264,
    518, 520, 262, 263, 518, 520, 263, 263, 518, 520, 262, 264,
    518, 519, 264, 264, 520, 528, 262, 264, 520, 528, 262, 262,
    520, 528, 264, 264, 528, 528, 264, 265, 265, 265, 265, 265,
    262, 265, 265, 265, 265, 265, 262, 518, 272, 518, 262, 519,
    262, 262, 262, 261, 262, 262, 262, 262, 264, 262, 262, 263,
    262, 262, 262, 520, 262, 262, 263, 520, 273, 309, 265, 265,
    265, 264, 520, 520, 265, 265, 265, 265, 265, 265, 265, 265,
    520, 265, 520, 520, 264, 520, 520, 520, 273, 273, 273, 273,
    273, 273, 291, 291, 291, 291, 291, 291, 309, 309, 309, 309,
    309, 309, 264, 264, 264, 264, 264, 264, 264, 264, 264, 264,
    520, 265, 265, 265, 265, 520, 265, 265, 273, 273, 291, 291,
    309, 309, 257, 257, 264, 264, 520, 273, 273, 273, 291, 291,
    291, 309, 309, 309, 264, 264, 264, 263, 263, 263, 263, 519,
    775, 273, 520, 273, 273, 273, 273, 291, 291, 291, 291, 272,
    528, 528, 264, 264, 264, 520, 273, 273, 273, 520, 273, 273,
    273, 273, 273, 309, 309, 309, 309, 264, 264, 264, 264, 257,
    263, 272, 272, 257, 257, 257, 257, 309, 528, 528, 265, 265,
    520, 262, 262, 264, 273, 262, 262, 291, 291, 262, 262, 264,
    291, 291, 309, 273, 291, 265, 520, 520, 520, 257, 273, 257,
    291, 257, 309, 262, 262, 263, 263, 519, 291, 291, 291, 291,
    291, 273, 273, 518, 262, 273, 273, 518, 520, 520, 520, 520,
    528, 262, 264, 273, 273, 273, 262, 273, 273, 273, 262, 264,
    273, 273, 264, 264, 273, 273, 262, 520, 262, 264, 273, 262,
    264, 264, 262, 263, 264, 262, 265, 265, 262, 273, 273, 262,
    262, 273, 262, 520, 520, 262, 264, 264, 520, 262, 264, 264,
    264, 264, 264, 264, 264, 520, 520, 264, 264, 264, 264, 264,
    520, 265, 265, 265, 265, 265, 265, 265, 265, 265, 265, 265,
    265, 520, 265, 265, 265, 520, 520, 265, 265, 520, 520, 520,
    264, 264, 264, 264, 520, 520, 265, 273, 265, 291, 265, 309,
    264, 264, 520, 520, 520, 265, 263, 263, 273, 273, 263, 263,
    264, 264, 264, 273, 273, 290, 290, 546, 257, 257, 257, 564
    };
    /* hot: pattern of ((z * 9 + charge + 4) * 2 + radical) */
    static const unsigned char VAL_H_PAT[2142] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    5, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 6, 0, 5, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 7, 0, 6, 7, 5, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 8, 0, 7, 8, 6, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 8, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 0, 9, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0,
    5, 0, 11, 0, 7, 0, 8, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    10, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 12, 0, 13, 14, 5, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 8, 0, 15, 16, 17, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 19, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0,
    0, 0, 20, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    10, 0, 0, 0, 21, 0, 0, 0, 22, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 23, 0, 0, 0, 18, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 1, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 1, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26, 0, 0, 0,
    1, 0, 1, 0, 0, 0, 10, 0, 27, 0, 28, 0, 29, 0, 30, 0,
    0, 0, 31, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    32, 0, 9, 0, 1, 0, 0, 0, 0, 0, 0, 0, 22, 0, 0, 0,
    22, 0, 33, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    34, 0, 0, 0, 4, 0, 9, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 34, 0, 35, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 10, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 36, 0, 5, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 8, 0, 15, 0, 29, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 19, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 37, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 39, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0,
    34, 0, 41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    34, 0, 0, 0, 4, 0, 9, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 22, 0, 42, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 34, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 35, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 21, 0, 29, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 36, 0,
    5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 43, 0,
    15, 0, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    18, 0, 19, 0, 22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 45, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 45, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    46, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 20, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 47, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 47, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 45, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 45, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 47, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    20, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 47, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 47, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 48, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 49, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0,
    0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 53, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 34, 0, 54, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 3, 0, 18, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 34, 0, 0, 0, 55, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    53, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 57, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 18, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 58, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    59, 0, 0, 0, 34, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 60, 0, 34, 0, 34, 0, 1, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 61, 0, 34, 0, 34, 0, 1, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 62, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 63, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 0, 0, 0, 0, 0,
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    26, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 26, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 26, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 64, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    /* hot: [pattern][bonds] -> hydrogens, -1 no rule, -2 an environment decides */
    static const signed char VAL_H_ROW[585] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 0, -1, -1, -1, -1, -1, -1, -1,
    0, 0, -1, -1, -1, -1, -1, -1, -1,
    0, -1, 0, -1, -1, -1, -1, -1, -1,
    4, 3, 2, 1, 0, -1, -1, -1, -1,
    3, 2, 1, 0, -1, -1, -1, -1, -1,
    2, 1, 0, -1, -1, -1, -1, -1, -1,
    1, 0, -1, -1, -1, -1, -1, -1, -1,
    -1, -2, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -2, -1, -1,
    0, 2, 1, 0, -1, -1, -1, -1, -1,
    2, 1, 0, -1, -1, -1, -2, -1, -1,
    3, 2, 1, 0, -2, 0, -1, -1, -1,
    2, 1, 0, 1, 0, -1, -1, -1, -1,
    2, 1, 0, -1, -2, -1, -2, -1, -1,
    1, 0, 1, 0, -1, -1, -1, -1, -1,
    -1, -1, -1, -2, -1, -2, -1, -1, -1,
    0, -1, -2, -1, -1, -1, -1, -1, -1,
    1, 0, -1, -2, -1, -2, -1, -2, -1,
    0, -1, -1, 0, -1, -1, -1, -1, -1,
    0, -1, -2, -2, 0, -1, -1, -1, -1,
    -1, -1, -2, -1, -1, -1, -1, -1, -1,
    0, -1, 0, -2, -2, -2, -1, -1, -1,
    0, -1, 0, 0, -2, -1, -2, -1, -1,
    0, -1, 0, -2, -2, -1, -2, -2, -1,
    0, -1, 0, 0, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -2, -2, -1, -1,
    -1, -1, -1, -1, -2, -1, -2, -1, -1,
    -1, -1, -1, -2, -1, -1, -1, -1, -1,
    0, -2, 0, 0, -1, -1, -1, -1, -1,
    0, -2, -1, -1, -1, -1, -1, -1, -1,
    0, -1, 0, -2, -1, -1, -1, -1, -1,
    0, 0, 0, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -2, -1, -1, -1, -1,
    0, -2, -1, 0, -1, -1, -1, -1, -1,
    0, -1, -1, 0, -1, 0, -1, -1, -1,
    0, -1, -2, -2, -2, -2, -1, -1, -1,
    0, -1, -1, -1, -2, -2, -2, -1, -1,
    0, -1, -1, -1, -2, -1, -1, -1, -1,
    0, -1, -1, -1, -2, -2, -1, -2, -2,
    0, -2, -2, 0, 0, -1, -2, -1, -1,
    0, 0, -2, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -2, -1, -1, -1,
    0, -1, -2, -1, -2, -1, -2, -1, -2,
    0, -1, -1, 0, -2, -1, -1, -1, -1,
    0, -1, -2, 0, -2, -1, -1, -1, -1,
    0, -1, -2, 0, -1, -1, -1, -1, -1,
    0, -2, -2, -2, 0, -1, -1, -1, -1,
    0, -1, -1, -1, -1, -2, -1, -1, -1,
    0, -1, -1, -1, -1, -1, -2, -1, -1,
    0, -1, -1, -1, -1, -1, -2, -1, -2,
    0, -2, -2, 0, 0, -2, -2, -1, -1,
    0, -1, 0, -1, -2, -1, -2, -1, -1,
    0, -1, -1, -2, -1, -1, -1, -1, -1,
    0, -1, 0, -1, -2, -1, -1, -1, -1,
    0, -2, -2, 0, -2, -2, -1, -1, -1,
    0, 0, -1, -1, -1, -2, -1, -1, -1,
    0, -1, -2, -2, 0, 0, -1, -1, -1,
    0, -1, -1, 0, 0, 0, 0, -1, -1,
    0, -1, 0, 0, 0, 0, 0, 0, -1,
    0, -1, -2, 0, 0, 0, 0, -1, -1,
    0, -1, 0, 0, 0, -1, -1, -1, -1,
    0, -1, -2, 0, 0, -1, -1, -1, -1,
    0, -1, -1, -1, 0, -1, -1, -1, -1
    };
    """
    const uint32_t VAL_KEY[576]
    const uint16_t VAL_KEY_OFF[576]
    const uint8_t VAL_KEY_LEN[576]
    const uint8_t VAL_H[1036]
    const uint16_t VAL_ENV_OFF[1036]
    const uint8_t VAL_ENV_LEN[1036]
    const uint16_t VAL_ENV[732]
    const uint8_t VAL_H_PAT[2142]
    const int8_t VAL_H_ROW[585]

DEF VAL_KEY_COUNT = 576
DEF VAL_H_PATTERNS = 65
# --- END GENERATED TABLES ---


cdef uint32_t val_key(uint32_t z, int charge, bint radical, uint32_t bonds) noexcept nogil:
    """The packed lookup key, or VAL_NO_KEY for a state no key could express.

    The guard is not paranoia: a charge outside the biased field would carry into the atomic
    number's bits and find another element's rows, which is the one failure mode of a packed key
    that a test on real inputs would never see.
    """
    if z < 1 or z > 118:
        return VAL_NO_KEY
    if charge < -VAL_CHARGE_BIAS or charge > 15 - VAL_CHARGE_BIAS:
        return VAL_NO_KEY
    if bonds > VAL_BONDS_MAX:
        return VAL_NO_KEY
    return ((z << VAL_KEY_Z_SHIFT) |
            (<uint32_t> (charge + VAL_CHARGE_BIAS) << VAL_KEY_CHARGE_SHIFT) |
            (<uint32_t> radical << VAL_KEY_RADICAL_SHIFT) | bonds)


cdef uint32_t val_lower_bound(uint32_t key) noexcept nogil:
    """The first index in VAL_KEY whose key is >= `key`; VAL_KEY_COUNT if there is none."""
    cdef uint32_t lo = 0
    cdef uint32_t hi = VAL_KEY_COUNT
    cdef uint32_t mid
    while lo < hi:
        mid = (lo + hi) >> 1
        if VAL_KEY[mid] < key:
            lo = mid + 1
        else:
            hi = mid
    return lo


cdef int val_key_find(uint32_t key) noexcept nogil:
    """The index of `key` in VAL_KEY, or -1.  10 comparisons over 576 keys."""
    cdef uint32_t at
    if key == VAL_NO_KEY:
        return -1
    at = val_lower_bound(key)
    if at < VAL_KEY_COUNT and VAL_KEY[at] == key:
        return <int> at
    return -1


cdef bint val_described(uint32_t z, int charge, bint radical) noexcept nogil:
    """Does the collection say anything at all about this element, charge and radical state?

    The VIOLATION/UNKNOWN boundary of `val_check`.  A key range scan rather than a second table:
    VAL_KEY is sorted and the packing puts `bonds` in the low bits, so every key for one
    `(z, charge, radical)` is one contiguous run and its lower bound is one binary search.
    """
    cdef uint32_t key = val_key(z, charge, radical, 0)
    cdef uint32_t at
    if key == VAL_NO_KEY:
        return False
    at = val_lower_bound(key)
    return (at < VAL_KEY_COUNT and
            (VAL_KEY[at] >> VAL_KEY_RADICAL_SHIFT) == (key >> VAL_KEY_RADICAL_SHIFT))


cdef bint val_env_ok(uint32_t rule, const uint16_t *env, uint32_t env_len) noexcept nogil:
    """Does this atom's neighbourhood contain everything the row demands?

    A MULTISET COMPARISON, and over a handful of entries that is two loops rather than a set plus a
    count dict: a subset test adds nothing once the counts are compared, since every count is at least
    one.  The multiplicity is the load-bearing half -- it is what stops a sulfone's `=O =O` row from
    firing on a sulfoxide.

    The requirement is a LOWER BOUND, never an exact match: one single-bonded carbon covers
    methanol and dimethyl ether alike.  Cost is bounded by the longest environment in the
    collection (7 entries) times this atom's degree.
    """
    cdef uint32_t off = VAL_ENV_OFF[rule]
    cdef uint32_t n = VAL_ENV_LEN[rule]
    cdef uint32_t i, j, need, have
    cdef uint16_t token

    for i in range(n):
        token = VAL_ENV[off + i]
        # a repeated token (`=O =O`) is counted once per occurrence and checked twice with the
        # same answer, which costs an iteration and saves a dedup pass over the table
        need = 0
        for j in range(n):
            if VAL_ENV[off + j] == token:
                need += 1
        have = 0
        for j in range(env_len):
            if env[j] == token:
                have += 1
        if have < need:
            return False
    return True


cdef int val_scan(uint32_t z, int charge, bint radical, uint32_t order_sum,
                  const uint16_t *env, uint32_t env_len, int want_h) noexcept nogil:
    """The one scanner both hydrogen questions are asked through; VAL_NO_RULE if nothing matches.

    Rows at a key are in the TSV's order, and that order is observable -- see this file's header.
    `want_h` is VAL_ANY_H to report the first match's count, or a count to require: that is the whole
    difference between "how many hydrogens" and "is this state described".  The filter has to keep
    scanning past a row it rejects rather than stop, because the verdict accepts any matching row and
    not the first.
    """
    cdef int at = val_key_find(val_key(z, charge, radical, order_sum))
    cdef uint32_t off, n, i
    cdef int h

    if at < 0:
        return VAL_NO_RULE
    off = VAL_KEY_OFF[at]
    n = VAL_KEY_LEN[at]
    for i in range(n):
        h = VAL_H[off + i]
        if want_h != VAL_ANY_H and want_h != h:
            continue
        if val_env_ok(off + i, env, env_len):
            return h
    return VAL_NO_RULE


cdef int val_implicit_h(uint32_t z, int charge, bint radical, uint32_t order_sum,
                        const uint16_t *env, uint32_t env_len) noexcept nogil:
    """How many hydrogens the collection gives this atom, or VAL_NO_RULE if it has no row for it.

    THE HOT PATH.  This runs on every atom of every molecule read from a format with a hydrogen
    convention, and it is one indexed load and a branch: no search, no environment compare, and
    none of the 583 rows that exist only to be checked against.  The dense table is generated as a
    projection of the same TSV -- `hydrogen_tables` explains why precomputing is sound and
    `test_the_dense_table_is_the_full_scan` proves it over every state in the domain.

    `order_sum` is the sum of the bond orders to EXPLICIT neighbours -- hydrogen atoms included,
    the implicit count excluded, since that is the answer.  An explicit H is an ordinary neighbour in
    both the sum and the environment; a row written against a hydrogen neighbour is unreachable
    otherwise.

    `env` and `env_len` are consulted for the 0.6% of atoms whose answer an environment decides --
    a sulfone, a nitro group, a phosphonic acid.  A caller with no environment to offer may pass
    `NULL, 0` and will get VAL_NO_RULE there rather than a guess.
    """
    cdef int h
    if (z < 1 or z > VAL_H_Z_MAX or charge < VAL_H_CHARGE_MIN or charge > VAL_H_CHARGE_MAX
            or order_sum > VAL_H_BONDS_MAX):
        return VAL_NO_RULE
    h = VAL_H_ROW[<uint32_t> VAL_H_PAT[((z * VAL_H_CHARGE_SPAN + (charge - VAL_H_CHARGE_MIN)) << 1)
                                       | <uint32_t> radical] * VAL_H_STRIDE + order_sum]
    if h >= 0:
        return h
    if h == VAL_H_CONSULT:
        return val_scan(z, charge, radical, order_sum, env, env_len, VAL_ANY_H)
    return VAL_NO_RULE


cdef int val_check(uint32_t z, int charge, bint radical, uint32_t order_sum,
                   const uint16_t *env, uint32_t env_len, uint32_t implicit_h) noexcept nogil:
    """VAL_VALID, VAL_VIOLATION or VAL_UNKNOWN for this exact state.  Never a reason to reject.

    Not `val_implicit_h(...) == implicit_h`: see the scanner's note on why a legal count need not
    be the chosen one.
    """
    if implicit_h > H_NIBBLE_MAX:
        # unstorable rather than merely unknown, and no row could accept it
        return VAL_VIOLATION if val_described(z, charge, radical) else VAL_UNKNOWN
    if val_scan(z, charge, radical, order_sum, env, env_len, <int> implicit_h) != VAL_NO_RULE:
        return VAL_VALID
    if val_described(z, charge, radical):
        return VAL_VIOLATION
    return VAL_UNKNOWN


cdef bint val_has_rules(uint32_t z, int charge, bint radical, uint32_t bonds) noexcept nogil:
    """Does any row cover this exact `(element, charge, radical, bonds)`, environment ignored?

    THE ANSWER IGNORES THE ENVIRONMENT, and callers depend on that boundary being separate from "no
    row matched": a resonance pass refuses a charge shift outright when no row covers the target
    valence at all.  A caller that only wants "what does the collection make of this atom" should use
    `val_check`; this is for the ones that must distinguish an unreachable valence from an unusual
    neighbourhood.
    """
    return val_key_find(val_key(z, charge, radical, bonds)) >= 0


cdef uint16_t *val_env_from_python(environment, uint32_t *env_len) except NULL:
    """A caller's `[(order, element), ...]` as VAL_ENV tokens; the caller frees the block.

    Bond order 4 and order 8 are REFUSED rather than counted or skipped.  Skipping order 8 silently
    would decide for the caller; refusing makes it state its policy, since "a dative bond contributes
    no electron pair and so nothing to the valence" is a modelling decision about a complex and not a
    fact in this collection.  Callers who have made that decision -- the MDL reader has -- filter
    before calling, which is exactly the visibility wanted.
    """
    cdef uint32_t n = len(environment)
    cdef uint32_t i = 0
    cdef uint32_t number
    cdef int order
    cdef uint16_t *block
    cdef object entry           # declared because warn.undeclared is an error's younger sibling

    # one entry allocated even for an empty environment, so a NULL return means only failure
    block = <uint16_t *> PyMem_Malloc((n + 1) * sizeof(uint16_t))
    if block is NULL:
        raise MemoryError()
    try:
        for entry in environment:
            order = entry[0]
            number = _to_atomic_number(entry[1])
            if order == 4:
                raise ValueError('bond order 4 has no valence rule: kekulise first, or decide '
                                 'your own policy -- the collection has no aromatic rows')
            if order == 8:
                raise ValueError('bond order 8 has no valence rule: decide whether a dative bond '
                                 'contributes to valence and pass the environment you meant')
            if order < 1 or order > VAL_ORDER_MAX:
                raise ValueError(f'bond order {order} is outside 1..{VAL_ORDER_MAX}')
            block[i] = <uint16_t> ((order << VAL_ENV_SHIFT) | number)
            i += 1
    except:
        PyMem_Free(block)
        raise
    env_len[0] = n
    return block


# The verdict crosses into Python as a string, the way `kekule_classify`'s atom classes do.  Not
# three module-level integer constants: under this build's `warn.undeclared` a module-level Python
# binding is a warning, so the alternative would be an import-and-compare dance for a value that
# ends up in a log line as text anyway.
cdef tuple VAL_VERDICT_NAMES = ('unknown', 'valid', 'violation')


def valence_implicit_h(element, int charge, bint radical, uint32_t order_sum, environment=()):
    """How many hydrogens the collection gives this atom, or None when it has no row for it.

    `element` is a symbol or an atomic number; `order_sum` is the sum of the bond orders to
    explicit neighbours, hydrogens included; `environment` is `[(order, element), ...]` for those
    same neighbours, again hydrogens included.

    None is not zero: see this file's header.
    """
    cdef uint32_t z = _to_atomic_number(element)
    cdef uint32_t env_len = 0
    cdef uint16_t *env = val_env_from_python(environment, &env_len)
    cdef int h
    try:
        h = val_implicit_h(z, charge, radical, order_sum, env, env_len)
    finally:
        PyMem_Free(env)
    if h == VAL_NO_RULE:
        return None
    return h


def valence_check(element, int charge, bint radical, uint32_t order_sum, uint32_t implicit_h,
                  environment=()):
    """`'valid'`, `'violation'` or `'unknown'` for this exact state.  Never a reason to reject.

    `'violation'` means the collection describes this element in this charge and radical state and
    no row accepts what you have.  `'unknown'` means it describes nothing there at all, which is a
    gap in the collection and not a claim about the molecule -- report it, count it, and leave the
    atom alone.
    """
    cdef uint32_t z = _to_atomic_number(element)
    cdef uint32_t env_len = 0
    cdef uint16_t *env = val_env_from_python(environment, &env_len)
    cdef int verdict
    try:
        verdict = val_check(z, charge, radical, order_sum, env, env_len, implicit_h)
    finally:
        PyMem_Free(env)
    return VAL_VERDICT_NAMES[verdict]


def valence_has_rules(element, int charge, bint radical, uint32_t bonds):
    """Does any row cover this `(element, charge, radical, bonds)` at all, environment ignored?

    "This valence is unreachable" rather than "this neighbourhood matched no row" -- see
    `val_has_rules`.
    """
    return val_has_rules(_to_atomic_number(element), charge, radical, bonds)


def valence_rules():
    """The compiled collection back as rows, in scan order, for the tests that re-derive it.

    `[(atomic_number, charge, radical, bonds, implicit_h, ((order, atomic_number), ...)), ...]` --
    the TSV's rows minus the provenance, which the C tables do not carry because nothing at query
    time may behave differently for a mined row than for a curated one.
    """
    cdef uint32_t k, i, j, off, n, eoff, en, key
    cdef list out = []
    cdef list env

    for k in range(VAL_KEY_COUNT):
        key = VAL_KEY[k]
        off = VAL_KEY_OFF[k]
        n = VAL_KEY_LEN[k]
        for i in range(n):
            env = []
            eoff = VAL_ENV_OFF[off + i]
            en = VAL_ENV_LEN[off + i]
            for j in range(en):
                env.append((VAL_ENV[eoff + j] >> VAL_ENV_SHIFT,
                            VAL_ENV[eoff + j] & ((1 << VAL_ENV_SHIFT) - 1)))
            out.append((key >> VAL_KEY_Z_SHIFT,
                        <int> ((key >> VAL_KEY_CHARGE_SHIFT) & 0xF) - VAL_CHARGE_BIAS,
                        (key >> VAL_KEY_RADICAL_SHIFT) & 1 != 0,
                        key & VAL_BONDS_MAX,
                        VAL_H[off + i],
                        tuple(env)))
    return out
