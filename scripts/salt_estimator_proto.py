# -*- coding: utf-8 -*-
"""Prototype for `chython/chemistry/tables/ionization.tsv` + `chemistry/_ionization.py`.

Two questions get confused under one name, `salt_estimator`:

  Q1  the record already carries the counterion -- what is in it?
      chython answers today: `thiele()` + `strip_salts()` over `tables/salts.tsv`.
  Q2  the record is the free base -- how many equivalents will it pick up?
      that is what the RDKit code computes, and what this file prototypes.

Q2's knowledge is ONE table with a pKa column, not two hand-copied SMARTS lists.  `FG_FOR_FA_LIST` and
`FG_FOR_TFA_LIST` differ by exactly the rings whose pKaH falls between formic acid and TFA -- pyrazole
2.5, pyridazine 2.3, pyrimidine 1.3 -- so the second list is the first list read at a lower threshold.
A ROW COUNTS WHEN `site.pka - titrant.pka >= DELTA`, and `DELTA=1` reproduces both lists including their
shared exclusion of anilines, which land at 4.6 - 3.75 = 0.85.

THREE RULES KEEP A COUNT FROM DOUBLE-COUNTING, and each one is a case the differential harness found:

  * A SITE IS AN ATOM, NOT A MATCH.  A symmetric pattern maps both ways round and a guanidine matches
    three of its own nitrogens; a dict keyed by the atom that `:1` landed on makes that one site.
  * ROWS ARE ORDERED, SPECIFIC BEFORE GENERAL, AND THE FIRST TO CLAIM AN ATOM TYPES IT.  A pyrimidine
    nitrogen is 1.5 and never also pyridine's 5.2; a xanthine is 0.6 and never imidazole's 7.0.
  * ONE BASIC SITE PER FUSED AROMATIC SYSTEM.  A second protonation of the same delocalised ring is
    orders of magnitude weaker than the first, so adenine takes one equivalent and not two.  Saturated
    rings are exempt: piperazine genuinely takes two.
"""
from chython import smarts


#: a row counts when the pKa gap to the titrant is at least this.  1.0 is what reproduces both RDKit
#: lists; 2.0 is the classical criterion for an isolable 1:1 salt and drops pyridines.
DELTA = 1.0

#  id                 role  pka    smarts (`:1` is the site atom)
#  `pka` is the site's own for an acid, its conjugate acid's for a base.  ORDER IS PRECEDENCE.
IONIZATION = [
    # --- bases ----------------------------------------------------------------------------------- #
    ('guanidine',       'base', 13.6, '[N;x0:1]=[C;D3;!R](-[N;z1;x0])-[N;z1;x0]'),
    ('amidine',         'base', 11.6, '[N;x0:1]=[C;D3;z2](-[#6])-[N;z1;x0]'),
    ('amine_2',         'base', 11.0, '[N;h1;D2;z1;x0:1](-[C;z1])-[C;z1]'),
    ('amine_1',         'base', 10.6, '[N;h2;D1;z1;x0:1]-[C;z1]'),
    ('amine_3',         'base', 10.0, '[N;h0;D3;z1;x0:1](-[C;z1])(-[C;z1])-[C;z1]'),
    ('ammonia',         'base',  9.2, '[N;*;h3:1]'),
    ('imine',           'base',  7.0, '[N;D2;z2;x0:1](-[C;z1])=[C;z2]'),
    # a xanthine's imidazole sits between two amide carbonyls and is not basic; it must precede the
    # two imidazole rows, which would otherwise read caffeine as one formate equivalent.
    ('xanthine',        'base',  0.6, '[N;a;D2;h0;r5:1]:[C;a;r5]:[C;a;r5;r6]:[C;a;r6]=[O]'),
    ('imidazole',       'base',  7.0, '[N;a;D2;h0;r5;x0:1]:[C;a;h1]:[N;a;D2;h1;r5]'),
    ('n_alkyl_azole',   'base',  7.0, '[N;a;D2;h0;r5;x0:1]:[C;a;h1]:[N;a;D3;r5](-[C;z1])'),
    # the WHOLE C-N(H)-N-C spine, not just the N-N: it is what excludes a tetrazole, far less basic.
    ('pyrazole',        'base',  2.5, '[C;a;r5]:[N;a;D2;h1;r5]:[N;a;D2;h0;r5:1]:[C;a;r5]'),
    ('diazine',         'base',  1.5, '[N;a;D2;h0;r6;x0:1]:[C;a]:[N;a;D2;h0;r6]'),
    ('pyridine',        'base',  5.2, '[N;a;D2;h0;r6;x0:1](:[C;a]):[C;a]'),
    # --- acids ----------------------------------------------------------------------------------- #
    ('sulfonic_acid',   'acid', -1.0, '[O;D1;z1;h1:1]-[S;D4](=[O])=[O]'),
    ('tetrazole',       'acid',  4.9, '[N;a;D2;h1;r5:1]:[N;a;D2;h0;r5]:[N;a;D2;h0;r5]'),
    ('carboxylic_acid', 'acid',  4.5, '[O;D1;z1;h1;x0:1]-[C;z2;x2;D3]=[O]'),
    ('phenol',          'acid', 10.0, '[O;D1;z1;h1;x0:1]-[C;a]'),
]

#: a titrant is a compound with a measured pKa; its own role says which sites it reaches.
TITRANTS = {'TFA': ('acid', 0.23), 'HCOOH': ('acid', 3.75), 'AcOH': ('acid', 4.76),
            'HCl': ('acid', -6.0), 'NH3': ('base', 9.25), 'NaOH': ('base', 15.7)}


def _compile(pattern):
    """The query, and the stable id its `:1` sits on -- the site atom the row types."""
    query = smarts(pattern)
    return query, next(a for a, n in query.map_numbers().items() if n == 1)


_COMPILED = [(i, role, pka, *_compile(p)) for i, role, pka, p in IONIZATION]


def ionizable_sites(molecule):
    """`((atom, role, pka, rule_id), ...)` -- one entry per atom, the first matching row wins."""
    seen = {}
    for rule_id, role, pka, query, subject in _COMPILED:
        for mapping in query.get_mapping(molecule):
            seen.setdefault(mapping[subject], (role, pka, rule_id))
    return tuple((a, *v) for a, v in sorted(seen.items()))


def _fused_systems(molecule):
    """Aromatic ring systems as atom sets, fused rings merged."""
    systems = []
    for ring in molecule.aromatic_rings:
        ring = set(ring)
        touching = [s for s in systems if s & ring]
        for s in touching:
            systems.remove(s)
            ring |= s
        systems.append(ring)
    return systems


def salt_equivalents(molecule, titrant, *, delta=DELTA):
    """Equivalents of `titrant` this molecule takes up.  An acid titrant counts basic sites."""
    role, pka = TITRANTS[titrant] if isinstance(titrant, str) else titrant
    wanted, sign = ('base', 1) if role == 'acid' else ('acid', -1)
    sites = [(a, p) for a, r, p, _ in ionizable_sites(molecule)
             if r == wanted and sign * (p - pka) >= delta]
    if wanted == 'acid':
        return len(sites)
    systems = _fused_systems(molecule)
    count, claimed = 0, set()
    for atom, _ in sorted(sites, key=lambda s: -s[1]):          # strongest base claims its system
        system = next((i for i, s in enumerate(systems) if atom in s), None)
        if system is None:
            count += 1
        elif system not in claimed:
            claimed.add(system)
            count += 1
    return count
