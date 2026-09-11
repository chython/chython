Standardization
===============

A reader stores and logs what a file says; it does not repair it. An illegal valence, a nonsense
charge, an aromatic ring with no Kekule form are all *parsed*, and what the reader thought is on
``mol.log``. Repair is this page: an explicit pipeline you run afterwards.

Two rules run through all of it. **A pass either repairs by design or refuses by design, never
both** — ``kekule()`` will separate the charges on an aromatic N-oxide to find a Kekule form and say
so, while ``thiele()`` is single-purpose and declines rather than repairing. And **a pipeline is never
stronger than its parts**: ``canonicalize()`` is literally the stages below run by hand, so each one
heals what its own writes made derivable.


The Pipeline
------------

========================= ================================================ ==========================
pass                      writes                                           returns
========================= ================================================ ==========================
``kekule()``              definite bond orders; **repairs** to find a form ``KekuleResult``
``standardize()``         group and metal-ligand rule tables               ``bool``
``implicify_hydrogens()`` hydrogen atoms folded into counts                ``int``, the count removed
``neutralize()``          charges paired off                               ``bool``
``thiele()``              the aromatic form; **refuses**, never repairs    ``ThieleResult``
``standardize_isomers()`` canonical placement of a mobile hydrogen         ``bool``
``canonicalize()``        all six, to a fixed point                        ``bool``
========================= ================================================ ==========================

Every pass mutates in place and writes to ``mol.log``. No pass takes a ``log=`` argument and nothing is
conditional: only readers, writers and ``depict`` have one. Each is also a plain function in
``chython``, so ``mol.standardize()`` and ``standardize(mol)`` are the same pass.

``canonicalize()`` really is its own stages, which is the point of running it rather than a hand-rolled
sequence — get the order wrong and two drawings of one compound stop agreeing:

.. testcode::

    from chython import (smiles, standardize, implicify_hydrogens, neutralize,
                         standardize_isomers)

    a = smiles('Oc1ccccn1')
    a.canonicalize()

    b = smiles('Oc1ccccn1')
    b.kekule()                    # 1. definite orders, so the group rules can match
    standardize(b)                # 2. repair the drawing
    implicify_hydrogens(b)        # 3. hydrogen ATOMS are a key difference, so they go
    neutralize(b)                 # 4. pair off the charges acids.tsv can pair off
    b.thiele()                    # 5. back to the aromatic form
    standardize_isomers(b)        # 6. canonical mobile-hydrogen placement

    print(a == b, a)

.. testoutput::

    True C=1C(=O)NC=CC=1

Every picture on this page is one pass drawn as an arrow: the left side is what the record stored, the
right side is what the pass left. A pass being a plain function of one molecule is all it takes.

.. testcode::

    from chython import ReactionContainer, canonicalize, explicify_hydrogens, fix_resonance, kekule

    def before_after(mol, pass_):
        after = mol.copy()
        pass_(after)
        return ReactionContainer([mol], [after]).depict()


Canonicalize
------------

``canonicalize()`` brings a molecule to the representation two drawings of one compound share, so that
``canonical_bytes``, ``__hash__`` and ``__eq__`` answer "same compound" rather than "same drawing".
The stage order above is a correctness constraint, not taste, and steps 2 to 6 run **to a fixed point**
rather than once: the placement stage can unblock a repair, a hydroxy-azine whose mobile hydrogen sits
on the very ring nitrogen the repair rule needs free being the shape that needs the second round. Five
rounds is the cap, and a molecule still changing at five is logged ``canonicalize:rounds`` as ``LOST``
rather than hung on — a result, but not one two drawings are guaranteed to share.

.. testcode::

    mol = smiles('OC(=O)c1ccccc1N(=O)=O')
    print(mol.canonicalize(), mol)

.. testoutput::

    True C(O)(=O)c1ccccc1[N+]([O-])=O

.. testcode::

    svg = before_after(smiles('Oc1ccccn1'), canonicalize)

.. figure:: images/standardize-canonicalize.svg
   :width: 420px

   2-hydroxypyridine as the file drew it and as ``canonicalize()`` left it. The ring is kekulized,
   the mobile hydrogen ends up on the nitrogen, and the pyridone is the form ``O=c1cccc[nH]1`` and
   ``OC1=NC=CC=C1`` also reach — which is what makes them one key.

The return value answers "was the molecule changed", not "does ``str()`` differ" — a drawing that was
already canonical is ``False`` even though the output SMILES walks the atoms in its own order:

.. testcode::

    mol = smiles('O=C(O)c1ccccc1')
    print(mol.canonicalize(), mol)

.. testoutput::

    False c1(ccccc1)C(=O)O

Charges are **paired off, not preserved atom by atom**, because step 4 runs ``neutralize()``. Glycine's
zwitterion and its neutral drawing therefore share a key, while the NET charge is untouched — sodium
acetate stays sodium acetate, having no proton in it to move:

.. testcode::

    a, b = smiles('[NH3+]CC(=O)[O-]'), smiles('NCC(=O)O')
    a.canonicalize()
    b.canonicalize()
    print(a == b, a)

.. testoutput::

    True C(CN)(=O)O

What the pipeline did is on ``mol.log``, always. Each entry is a ``LogRecord(rule, atoms, message,
severity, stage, subject)`` tagged with the stage that wrote it, so one log reads as the trace:

.. testcode::

    mol = smiles('c1ccccc1N(=O)=O')
    mol.canonicalize()

    for record in mol.log:
        print(record.stage, record.rule, record.atoms)

.. testoutput::

    kekule kekule:kekulized (1, 2, 3, 4, 5, 6)
    standardize groups:12 (6, 7, 8, 9)
    thiele thiele:aromatized (1, 2, 3, 4, 5, 6)

The ring is kekulized, the nitro group is normalized by ``groups:12`` on its own four atoms, the ring
is aromatized again. ``mol.log`` accumulates across calls, so clear it (``del mol.log[:]``) to read one
call's records alone, and it has the three queries a pipeline is usually asked about:

.. testcode::

    print([r.rule for r in mol.log.by_stage('standardize')])
    print(mol.log.repaired(), mol.log.lost(), mol.log.refused())

.. testoutput::

    ['groups:12']
    [] [] []

Two options:

.. testcode::

    mol.canonicalize(
        fix_tautomers=True,    # forwarded to standardize(); False gives up part of the guarantee
        keep_kekule=False,     # True costs a second kekulisation rather than a skipped thiele()
    )

``fix_tautomers=False`` withholds the local repair rules — ``Oc1ccccn1`` and ``O=c1cccc[nH]1`` stop
hashing equal — and deliberately does **not** reach ``standardize_isomers()``, which picks between two
valid drawings rather than repairing a wrong one.


Aromaticity
-----------

``kekule()`` writes definite bond orders, ``thiele()`` writes the aromatic form. Both return a result
object rather than a bare bool, because both have a third thing to say: which systems they could not
do:

.. testcode::

    mol = smiles('c1ccccc1')

    r = mol.kekule()
    print(r.changed, r.unresolved, mol)

    r = mol.thiele()
    print(r.changed, r.refused, mol)

.. testoutput::

    True [] C1=CC=CC=C1
    True [] c1ccccc1

.. testcode::

    svg = before_after(smiles('c1ccccc1'), kekule)

.. figure:: images/standardize-kekule.svg
   :width: 340px

   The dashed inner bonds are the aromatic order the reader stored; ``kekule()`` writes the definite
   alternation. Benzene is the case where ``thiele()`` takes that straight back.

``kekule()`` repairs by design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ring written aromatic that has no Kekule form as drawn is a mis-drawing, and ``kekule()`` is the pass
that fixes it. An aromatic N-oxide gets its charges separated, because that is the form the ring needs:

.. testcode::

    mol = smiles('c1ccn(O)cc1')
    mol.kekule()
    print(mol)
    print(mol.log[0].rule, mol.log[0].severity)

.. testoutput::

    [N+]=1([O-])C=CC=CC=1
    kekule:cation-oxide repaired

.. testcode::

    svg = before_after(smiles('c1ccn(O)cc1'), kekule)

.. figure:: images/standardize-n-oxide.svg
   :width: 420px

   The drawing on the left has no Kekule form: a six-membered aromatic ring needs a double bond at
   every atom, and a neutral nitrogen already spending two ring bonds and a hydroxyl has no valence
   left for one. ``kekule()`` separates the charges to earn it, and logs the change ``repaired``.

An unspecified bond between two aromatic atoms is aromatic per OpenSMILES, so biphenyl arrives with an
aromatic bond in no ring at all. That is the reader being right, and ``kekule()`` coping:

.. testcode::

    mol = smiles('c1ccccc1c1ccccc1')
    mol.kekule()
    print(mol)
    print([(r.rule, r.severity) for r in mol.log])

.. testoutput::

    C=1C=C(C=CC=1)C1=CC=CC=C1
    [('kekule:acyclic-bond', 'repaired'), ('kekule:kekulized', 'info'), ('kekule:kekulized', 'info')]

A system with genuinely no Kekule form does **not** raise and does not stop the pipeline. It comes back
in ``.unresolved``, is logged ``LOST``, and the atoms left in a bad state are what ``check_valence()``
finds — cyclopentadiene written aromatic without the charge that would earn it:

.. testcode::

    mol = smiles('c1cccc1')
    r = mol.kekule()
    print(r.unresolved)
    print(mol.log[0].rule, mol.log[0].severity)
    print(mol.check_valence())

.. testoutput::

    [(1, 2, 3, 4, 5)]
    kekule:no-kekule-form lost
    [(3, 'violation')]

``thiele()`` refuses instead
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``thiele()`` writes bond orders and nothing else. A ring that looks like a Kekule aromatic and is
declined for a fixable reason comes back in ``.refused`` with the reason logged. A radical brings one pi
electron and not two:

.. testcode::

    mol = smiles('C1=CC=CC=C1 |^1:0|')     # phenyl radical, in CXSMILES
    r = mol.thiele()
    print(r.changed, r.refused)
    print(mol.log[-1].rule, mol.log[-1].severity)

.. testoutput::

    False [(1, 2, 3, 4, 5, 6)]
    thiele:radical-ring refused

The two are therefore **not a symmetric round trip**. ``kekule()`` on the N-oxide above changed the
drawing for good; ``thiele()`` takes that repaired form back to aromatic, and the charge separation
stays, because undoing a repair is not aromatization's job:

.. testcode::

    mol = smiles('c1ccn(O)cc1')
    mol.kekule()
    mol.thiele()
    print(mol)

.. testoutput::

    c1[n+](cccc1)[O-]


Functional Groups
-----------------

``standardize()`` normalizes functional groups. Two tables, each row a SMARTS with the atom and bond
edits to apply, and each recording under its own table-qualified id with the row's own ``why`` as the
message:

========================== ========== =========================================================
table                      ids        is
========================== ========== =========================================================
``standardize_groups.tsv`` ``groups`` a functional group drawn some other way than the one form
``standardize_metals.tsv`` ``metals`` a metal-ligand bond drawn covalent that is dative
========================== ========== =========================================================

.. testcode::

    mol = smiles('c1ccccc1N(=O)=O')   # nitro drawn as a pentavalent nitrogen
    print(mol.standardize(), mol.log[0].rule, mol.log[0].severity)
    print(mol.log[0].message.split(';')[0])

.. testoutput::

    True groups:12 info
    A pentavalent nitrogen carrying an N=O and a second double bond

.. testcode::

    svg = before_after(smiles('c1ccccc1N(=O)=O'), standardize)

.. figure:: images/standardize-nitro.svg
   :width: 460px

   ``groups:12`` on nitrobenzene. Both drawings are accepted input; the charge-separated one is the
   single form the rest of the library then matches against.

Both tables log ``INFO``: the group was drawn in one of its accepted forms and is now in the chosen one.
A ``metals`` row installs the **dative bond** (order 8) that records coordination, which is
exactly the bond ``split_salts()`` reads as its exemption signal — so the two passes compose, the second
declining what the first said to preserve:

.. testcode::

    mol = smiles('C(#N)[Fe]')             # a cyanide ligand drawn with a covalent Fe-C bond
    mol.standardize()
    print(mol, mol.order_of(1, 3))
    print(mol.split_salts(), [(r.rule, r.severity) for r in mol.log])

.. testoutput::

    [C-](#N)~[Fe+] 8
    False [('metals:01', 'info'), ('salts:metal', 'refused')]

``remove_coordinate_bonds()`` is the complement of that row: it deletes every order-8 bond in place and
answers how many went, for a consumer whose representation has no dative bond to give — a fingerprint
scheme, a format field, another toolkit. It is the deletion only; the charges the ``metals`` row wrote
stay written.

.. testcode::

    print(mol.remove_coordinate_bonds(), mol)

.. testoutput::

    1 [C-]#N.[Fe+]

A row's ``tautomer`` column marks it as one of the local repairs ``fix_tautomers=False`` withholds; the
rest run either way.


Hydrogens
---------

**A parsed molecule already has its implicit hydrogen counts**, filled at read time by one algorithm
every reader shares. Nothing here has to be called to get them:

.. testcode::

    mol = smiles('CCO')
    print({n: mol.implicit_h_of(n) for n in mol})

.. testoutput::

    {1: 3, 2: 2, 3: 1}

``explicify_hydrogens()`` turns those counts into atoms and ``implicify_hydrogens()`` folds them back.
Both answer how many hydrogens they moved:

.. testcode::

    print(mol.explicify_hydrogens(), mol)
    print(mol.implicify_hydrogens(), mol)

.. testoutput::

    6 O(C([H])(C([H])([H])[H])[H])[H]
    6 C(C)O

.. testcode::

    svg = before_after(smiles('CCO'), explicify_hydrogens)

.. figure:: images/standardize-explicit-h.svg
   :width: 420px

   The same six hydrogens, first as counts on three heavy atoms and then as six atoms with their own
   ids. Hydrogen atoms are a key difference, which is why ``canonicalize()`` folds them back.

Exactly one class of atom is left unknown at read time: the aromatic pnictogen whose class the ring
decides — pyrrole's nitrogen carries a hydrogen, pyridine's does not, and an aromatic ring alone does
not say which. The count is a registered *unknown* rather than a silent zero, so it reads back as
``None``, ``check_valence()`` calls it ``'unknown'`` rather than a violation, and ``thiele()`` refuses
the ring because that count **is** the atom's aromatic class:

.. testcode::

    from chython import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        ring = [mol.add_atom('N')] + [mol.add_atom('C') for _ in range(4)]
        for i in range(5):
            mol.add_bond(ring[i], ring[(i + 1) % 5], 4)      # 4 is the aromatic order

    print(mol.implicit_h_of(1), mol.check_valence())
    print(mol.thiele().refused, mol.log[0].rule)

.. testoutput::

    None [(1, 'unknown')]
    [(1, 2, 3, 4, 5)] thiele:unknown-h

``unknown_h_count`` asks the same question about the whole record, and it is the one check to run before
trusting anything a hydrogen count feeds: a formula, a mass (a lower bound while the count is non-zero),
or an ``h``/``H`` query primitive, which cannot match such an atom in either direction. A count and not a
list -- the atoms themselves are the ones whose ``implicit_h_of`` is ``None``. Five here and not one,
because ``add_atom`` states no count and a hand edit derives nothing: every reader runs the derivation,
and a caller building a molecule by hand runs it when the structure is finished. The same ring read from
a CTAB with order-4 bonds reports 1, the carbons derived and only the pnictogen left to the ring:

.. testcode::

    print(mol.unknown_h_count, mol.implicit_h_of(1))
    mol.derive_hydrogens(fill_only=True)      # the read-time pass; it answers {n: reason} for the rest
    print(mol.unknown_h_count, mol.implicit_h_of(2))

.. testoutput::

    5 None
    1 1

``kekule()`` closes that class itself — its definite orders are exactly what decides the count — which
is why ``canonicalize()`` needs no separate hydrogen-filling stage:

.. testcode::

    mol.kekule()
    print(mol.implicit_h_of(1), mol.unknown_h_count, mol)

.. testoutput::

    1 0 C=1C=CNC=1

``mol.calc_implicit(n)`` recomputes one atom's count and is the call to reach for after an edit that
changed a bond order — unlike ``derive_hydrogens(fill_only=True)`` it replaces a count already stored.
It answers ``None`` where nothing local derives one, storing ``H_UNKNOWN`` rather than a guessed zero,
and it never raises:

.. testcode::

    mol = smiles('CCO')
    mol.set_hydrogens(1, 0)
    print(mol.calc_implicit(1), mol.calc_implicit(3))

.. testoutput::

    3 1

``chython.calc_implicit(mol, n)`` and ``check_valence`` delegate to the same derivation, so there is no
second copy to disagree with this one.


Neutralize
----------

``neutralize()`` moves a proton from a cation that has one onto an anion that can take one, leaving both
ends neutral. Only charges and implicit hydrogen counts are written: no bond is cut and no component is
deleted, which is what separates it from ``split_salts()`` below.
``canonicalize()`` runs it as a stage, so a zwitterion and its neutral drawing share a key:

.. testcode::

    mol = smiles('[NH3+]CC(=O)[O-]')       # glycine, zwitterionic
    print(mol.neutralize(), mol)

    salt = smiles('C[NH3+].[Cl-]')         # methylammonium chloride
    salt.neutralize()
    print(salt)

.. testoutput::

    True C(CN)(=O)O
    CN.Cl

.. testcode::

    svg = before_after(smiles('[NH3+]CC(=O)[O-]'), neutralize)

.. figure:: images/standardize-neutralize.svg
   :width: 430px

   Glycine's zwitterion and its neutral drawing are one compound, and the proton moves within it: no
   bond is cut, no atom is added or deleted, and the net charge was zero on both sides.

Nothing overshoots zero: a component is only ever taken closer to neutral. Both oxygens of a nitrate
match, and only one is protonated, because the second would take the nitrate to +1; sulfate has two
charges to spend and takes two protons:

.. testcode::

    for s in ['[NH3+]CC[NH3+].[O-][N+](=O)[O-]', '[NH3+]CC[NH3+].[O-]S(=O)(=O)[O-]']:
        mol = smiles(s)
        mol.neutralize()
        print(mol)

.. testoutput::

    [N+](=O)(O)[O-].C([NH3+])CN
    O=S(O)(=O)O.C(N)CN

The default ``keep_charge=True`` moves protons in pairs, so the net charge is preserved exactly and a
record that cannot be balanced comes back partly neutral, as the nitrate above did. A cation with no
proton to give — a quaternary ammonium, a metal — keeps its counterion, and every site left charged is
logged with the row that recognized it:

.. testcode::

    mol = smiles('[NH3+]CC[NH3+].[Cl-]')
    mol.neutralize()
    print(mol.log.refused()[0].rule)

.. testoutput::

    acids:ammonium

``keep_charge=False`` lets one side act alone, as far as its own component's charge allows, which is how
a lone ion is neutralized:

.. testcode::

    mol = smiles('C[NH3+]')
    print(mol.neutralize(), mol.neutralize(keep_charge=False), mol)

.. testoutput::

    False True CN

The sites are ``chython/chemistry/tables/acids.tsv``, whose ``h`` primitive reads implicit hydrogens —
call ``implicify_hydrogens()`` first on a molecule that carries hydrogen atoms. Whether the neutral form
is allowed at all is asked of the shared valence collection rather than written into the table, so a
matched row can still refuse:

.. testcode::

    mol = smiles('[O-]')
    print(mol.neutralize(keep_charge=False), mol)
    print(mol.log[0].rule, mol.log[0].message)

.. testoutput::

    False [O-]
    acids:hydroxide atom 1 was left charged: neutral with 1 implicit hydrogen(s) is a valence violation


Mobile Hydrogens
----------------

``standardize_isomers()`` places a mobile hydrogen canonically, so two drawings of one compound store
the same molecule. It decides three shapes: an aromatic ring system, a ring system stored Kekule as a
lactam is, and an acyclic amidine or guanidine. Neither form is a precondition, and it has to be that
way: a pyridone comes out of the pipeline stored Kekule, ``thiele()`` declining to aromatize a lactam
ring, so requiring the aromatic form would put the whole lactam family out of reach:

.. testcode::

    for a, b in [('Cc1cnc[nH]1', 'Cc1c[nH]cn1'),           # aromatic: 4-methylimidazole
                 ('CC1=CC=NC(=O)N1', 'CC1=NC(=O)NC=C1'),   # Kekule: 4-methylpyrimidin-2-one
                 ('CC(N)=NC', 'CC(=N)NC')]:                # acyclic: N-methylacetamidine
        ma, mb = smiles(a), smiles(b)
        ma.standardize_isomers()
        mb.standardize_isomers()
        assert ma.canonical_bytes == mb.canonical_bytes, (a, b)

.. testcode::

    svg = before_after(smiles('Cc1cnc[nH]1'), standardize_isomers)

.. figure:: images/standardize-isomers.svg
   :width: 400px

   4-methylimidazole: the hydrogen moves to the ring nitrogen further from the methyl. Neither side is
   a mis-drawing, so the records are ``INFO`` and the choice is only about which one both drawings of
   this compound store.

It picks between two valid drawings rather than repairing a wrong one, so its log records are ``INFO``
and ``fix_tautomers=False`` deliberately does not switch it off.


Resonance
---------

``fix_resonance()`` moves charges and radicals along an alternating path so that a dipolar or biradical
drawing becomes the covalent one. It is not part of ``canonicalize()``: call it when a source is known
to draw amides as zwitterions. ``mol.fix_resonance()`` is the same pass as a method, as for every pass
on this page:

.. testcode::

    for s in ['[O-]C=[NH+]C',                  # N-methylformamide as its zwitterion
              '[CH2]C=C[CH2] |^1:0,3|']:       # a biradical four atoms apart
        mol = smiles(s)
        print(mol.fix_resonance(), mol, [r.rule for r in mol.log])

.. testoutput::

    True C(NC)=O ['resonance:donor_anion']
    True C(C=C)=C ['resonance:radical']

.. testcode::

    svg = before_after(smiles('[O-]C=[NH+]C'), fix_resonance)

.. figure:: images/standardize-resonance.svg
   :width: 410px

   N-methylformamide drawn as its zwitterion and left covalent. Two atoms lose a charge and the two
   bonds along the path swap orders — one path, written whole or not at all.

The walk is bounded by bond order alone: only orders 1 to 3 are crossed and only orders 1 to 3 are
produced, so an aromatic bond (4) and a dative one (8) are impassable and a resonance path never runs
through a ring the aromatizer has claimed. Every atom the path touches then has to pass the shared
valence collection, and one that does not refuses the whole path — nothing is written half-way.


Salts
-----

Two passes read ``chython/chemistry/tables/salts.tsv``, and both act on all 93 metals the core's ``[M]``
accepts. One edits and one reports; **neither deletes a component**.

``split_salts()`` cuts the ionic bond and moves the charge onto the two ends. The atom count does not
change:

.. testcode::

    mol = smiles('CC(=O)O[Na]')
    print(mol.split_salts(), mol)

.. testoutput::

    True C(C)([O-])=O.[Na+]

The test is **all-or-nothing per cation atom**, which is what makes generality over 93 metals safe: a
dative bond, a neighbour that matches no acceptor row, an untabulated resulting charge or an implicit
hydrogen on the cation refuses the whole atom and logs why. Cisplatin, ferrocene and the metal
carbonyls therefore come back intact rather than half-split. A dative bond is the *exemption* signal and
never the trigger — ``standardize()`` installs it to record coordination that must be preserved:

.. testcode::

    mol = smiles('N[Pt](N)(Cl)Cl')
    print(mol.split_salts(), mol)
    print(mol.log.refused()[0].rule)

.. testoutput::

    False N[Pt](Cl)(N)Cl
    salts:metal

``decompose_salts()`` reads the record as a compound plus what was drawn beside it, and **changes
nothing and logs nothing** — the return value is the whole report. Four fields: the compound, then the
counterions and solvates keyed by their row, and lone cations keyed by element symbol:

.. testcode::

    r = smiles('NCC(=O)O.OC(=O)C(F)(F)F.O').decompose_salts()   # glycine TFA salt, monohydrate
    print([str(p) for p in r.parents], r.counterions, r.solvates, r.cations)

.. testoutput::

    ['C(CN)(=O)O'] {'salts:tfa': 1} {'salts:water': 1} {}

The work runs on a copy that is hydrogen-implicified, salt-split, neutralized and aromatized, so the
counts compare across a corpus rather than across a drawing — the three ways of writing sodium acetate
agree, and neutralizing first is why no conjugate base needs a row of its own:

.. testcode::

    for s in ['CC(=O)O[Na]', 'CC(=O)[O-].[Na+]', 'CC(=O)O.[Na+]']:
        r = smiles(s).decompose_salts()
        print([str(p) for p in r.parents], r.counterions, r.cations)

.. testoutput::

    ['C(C)(=O)O'] {} {'Na': 1}
    ['C(C)(=O)O'] {} {'Na': 1}
    ['C(C)(=O)O'] {} {'Na': 1}

A tabulated species is a counterion only when **something else is there to be the compound**, which is
why sodium acetate answers acetic acid rather than an empty ``parents`` and one equivalent of
``salts:acetic``. Every component is a solvate row, a species row or unmatched, and the three decide
together: the unmatched components are the parents when there are any, failing that the species rows are
with the solvates counted, and failing that the solvates are.

.. testcode::

    r = smiles('CC(=O)O.O').decompose_salts()      # acetic acid monohydrate
    print([str(p) for p in r.parents], r.solvates)
    print([str(p) for p in smiles('O.O').decompose_salts().parents])

.. testoutput::

    ['C(C)(=O)O'] {'salts:water': 1}
    ['O']

``cations`` is keyed by element symbol and not by row id because all 93 metals share one row: keyed by
row, a sodium and a potassium salt would be indistinguishable. Parents dedup by canonical bytes, so two
drawn equivalents of one compound are one parent and two enantiomers are two.


Stereo and Isotopes
-------------------

``clean_stereo()`` drops every stereo annotation, leaving the constitution; ``clean_isotopes()`` drops
every isotope label. Both are how you ask a question about a constitution rather than about a specific
labelled stereoisomer:

.. testcode::

    mol = smiles('C[C@H](N)C(=O)O')
    mol.clean_stereo()
    print(mol)

    mol = smiles('[13CH4]')
    mol.clean_isotopes()
    print(mol)

.. testoutput::

    C(C(C)N)(=O)O
    C

.. testcode::

    svg = before_after(smiles('C[C@H](N)C(=O)O'), lambda m: m.clean_stereo())

.. figure:: images/standardize-clean-stereo.svg
   :width: 430px

   One enantiomer of alanine, then alanine. The wedge is gone because the parity behind it is gone,
   and the constitution is untouched.

Neither is a ``canonicalize()`` stage: a stereoisomer and a labelled compound are different compounds,
so discarding either is the caller's decision and never the pipeline's.


Bond Perception
---------------

A file that states coordinates returns a record of atoms and coordinates rather than a molecule, and two
passes turn one into chemistry. They answer different questions and are called in this order:

===================  ===============================================  ==========================
pass                 the question                                     what it reads
===================  ===============================================  ==========================
``perceive_bonds``   which pairs of atoms are bonded                  one stored 3D model
``saturate``         at what order the bonds already there are        the hydrogen counts
===================  ===============================================  ==========================

Neither runs on read, and neither has a method: their callers are the coordinate readers, and where
those hand their record over is still being designed; a method now would pin the shape first.

PDB and mmCIF state connectivity, so ``chython.formats.pdb.build_molecule()`` hands ``saturate()`` a
molecule that already has its bonds and ``perceive_bonds()`` is not in that path. XYZ states no bond at
all, so a frame goes through both.


``perceive_bonds()`` finds the pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One threshold and nothing else: a pair of atoms closer than ``radius_multiplier`` times the sum of their
covalent radii is bonded, and a pair further apart is not. No element special case, no ring closure and
no hydrogen rule, so the answer is a function of the geometry the file stated. The radii are
``chython/chemistry/tables/covalent_radii.tsv`` and cover elements 1–96; an atom whose element has no row
— the ``R`` marker among them — takes no bond and is named in a ``perceive:no-radius`` record.

Every bond it adds is **order 1**, the order no valence rule has to justify, and a bond the molecule
already holds is left exactly as it is, order included. ``model=`` picks which stored geometry is read,
and a molecule carrying none raises ``IndexError`` rather than being handed an invented one.

.. testcode::

    from chython import implicify_hydrogens, perceive_bonds, saturate, xyz
    from chython.formats.xyz import build_molecule

    acetonitrile = """6
    from a geometry optimization: C-C 1.458, C-N 1.157, C-H 1.087 A
    C    0.000000    0.000000    0.000000
    C    0.000000    0.000000    1.458000
    N    0.000000    0.000000    2.615000
    H    1.024000    0.000000   -0.362000
    H   -0.512000   -0.887000   -0.362000
    H   -0.512000    0.887000   -0.362000
    """
    mol = build_molecule(xyz(acetonitrile)[0])
    print(perceive_bonds(mol), mol)      # five bonds, every one single
    print(saturate(mol), mol)            # the zeroed hydrogen counts force the nitrile
    implicify_hydrogens(mol)             # XYZ writes hydrogens as atoms; a chemist does not
    print(mol == smiles('CC#N'))         # as a STRUCTURE: the writer emits the file's atom order

.. testoutput::

    True [C]([N])C([H])([H])[H]
    True C(#N)C([H])([H])[H]
    True

The default multiplier of ``1.25`` answers every bond length and every nonbonded contact in
``chython/chemistry/test/test_covalent_radii_tsv.py``, two corpora of measured geometries that leave it
the window **1.2386 to 1.2961** — the floor being fluorine's long F–F bond and the ceiling
cyclobutadiene's tight transannular C···C, which a looser threshold turns into a bicyclobutane. That is a
4.6% window, so the knob is not a free parameter; a caller reading a metal cluster or a stretched
transition state moves it and reads the log.

**One threshold cannot answer every geometry.** Bicyclo[1.1.1]pentane's bridgehead carbons are 1.845 Å
apart and not bonded, and no multiplier rejects that pair while still reaching F–F: the pass bonds them.
A single distance rule is what ``perceive_bonds()`` is, and where a structure is strained enough for a
contact and a bond to overlap, the ``perceive:bonds`` record is what says how many bonds the geometry
was read as having.


``saturate()`` finds the orders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It perceives **orders and never connectivity**: it raises existing bonds and never adds one. The stated
hydrogen counts are what it reads — a count is what *forces* a multiple bond:

.. testcode::

    from chython import saturate

    mol = smiles('CC=N')
    with mol.edit():                  # state that neither heavy atom carries a hydrogen
        mol.set_hydrogens(2, 0)
        mol.set_hydrogens(3, 0)

    svg = before_after(mol, saturate)  # the same pass, drawn on a copy
    print(saturate(mol), mol)

.. testoutput::

    True C(C)#N

.. figure:: images/standardize-saturate.svg
   :width: 340px

   Nothing in the left drawing says the nitrogen is a nitrile; the two zeroed hydrogen counts do, by
   leaving the order sum no other way to close. The double bond already there is a stated fact, so it
   is raised and never lowered.

Where nothing pins the sum every bond stays single, and the pass says which case that was. An
``H_UNKNOWN`` count is slack, not zero, so a ring of six carbons none of whose hydrogens are known is
cyclohexane and not benzene — the smallest answer, never invented unsaturation:

.. testcode::

    from chython import H_UNKNOWN, MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        ring = [mol.add_atom('C') for _ in range(6)]
        for i in range(6):
            mol.add_bond(ring[i], ring[(i + 1) % 6], 1)
        for n in ring:
            mol.set_hydrogens(n, H_UNKNOWN)

    print(saturate(mol), mol, mol.log[0].rule, mol.log[0].severity)

.. testoutput::

    True C1CCCCC1 saturate:unforced info

Assignment is **all-or-nothing per connected fragment**, and writes bond orders only — never a charge, a
radical, a hydrogen count or the set of bonds, and never a lowered order, a multiple bond already in the
molecule being a stated fact. An unsatisfiable fragment is ``False`` plus a log, never a silently
unbalanced answer:

.. testcode::

    mol = smiles('C1CCCC1')
    with mol.edit():                  # five atoms each needing one unit of extra order: five is odd
        for n in mol.atom_numbers:
            mol.set_hydrogens(n, 1)
    print(saturate(mol))
    print(mol.log.refused()[0].rule)

.. testoutput::

    False
    saturate:no-valence-state


Contracted Groups
-----------------

A drawing that writes ``OMe`` on one atom has stated a methoxy group, and the reader invents neither the
oxygen nor the methyl: the text is kept as the atom's alias and the atom itself is the marker, element 0.
``expand_abbreviations()`` turns the labels ``chython/chemistry/tables/abbreviations.tsv`` knows into
atoms. It is the third pass on this page with **no method** — what a drawing wrote on an atom is a fact
about a file, so the pass belongs beside the reader that stored the label and not on every molecule.

.. testcode::

    from chython import expand_abbreviations

    mol = smiles('c1ccc(*)cc1* |$;;;;OMe;;;NO2$|')
    print(mol.aliases)
    print(expand_abbreviations(mol), mol)
    print([r.rule for r in mol.log.repaired()])

.. testoutput::

    {5: b'OMe', 8: b'NO2'}
    True c1([N+](=O)[O-])cccc(OC)c1
    ['abbreviations:OMe', 'abbreviations:NO2']

The label survives as an alias until a row claims it, so a spelling the table does not hold costs nothing
and the table is safe to grow. The absences are deliberate: ``R``, ``X``, ``Pol``, ``Ar`` and ``PEG`` name
a marker or a family rather than a structure, and stay R atoms carrying their label — see :doc:`io` for
what a marker is.

.. testcode::

    mol = smiles('CC(C)C |$;;;Pol$|')
    print(expand_abbreviations(mol), mol, list(mol.log))

.. testoutput::

    False C(C)(C)C |$;Pol;;$| []

Each row's fragment is a molecule with one ``*``, and that marker is what makes the hydrogen counts come
out right: ``*NS(=O)(=O)c1ccc(C)cc1`` gives the nitrogen one hydrogen, being the count for an atom that
already carries the bond to the rest of the structure. The labelled atom is then **transmuted** into the
fragment's attachment atom and the rest grafted on, rather than deleted and rebuilt — its stable id stays,
so a neighbouring stereocentre's parity remains a statement about the frame it was written over:

.. testcode::

    mol = smiles('*[C@@H](N)CBr |$OMe;;;;$|')
    print(expand_abbreviations(mol), mol)

.. testoutput::

    True [C@@H](CBr)(N)OC

The group takes the place the label held, mirror image included. Per site the pass is all-or-nothing: a
label reached by anything but one single bond, or on an atom the record also gave a charge, a radical or
an isotope, keeps its alias and gets a ``REFUSED`` record — two statements about one atom, and nothing
here ranks them. Grafted atoms take the labelled atom's coordinates, so a record with a depiction needs
``clean2d()`` afterwards.


Valence Checking
----------------

``check_valence()`` returns ``(atom, verdict)`` pairs against the same valence collection every other
pass consults. Two verdicts, and the difference matters: ``'violation'`` is a state no rule admits,
``'unknown'`` is a state nothing has decided yet.

.. testcode::

    print(smiles('C=N=Cc1ccccc1').check_valence())    # a pentavalent nitrogen
    print(smiles('O(C)(C)C').check_valence())         # trimethyloxonium drawn neutral
    print(smiles('c1ccccc1').check_valence())         # nothing wrong
    print(smiles('c1cc[nH]c1').check_valence())       # pyrrole: the ring decided
    print(smiles('F[Xe]F').check_valence())           # a hypervalent noble gas is tabulated

.. testoutput::

    [(2, 'violation')]
    [(1, 'violation')]
    []
    []
    []

Only the aromatic pnictogen above answers ``'unknown'``, and it stops doing so as soon as ``kekule()``
has run. A ``'violation'`` is also how the atoms an unresolved aromatic system left behind are found,
``kekule()`` having logged the system itself as ``LOST`` rather than raising.


Deduplication
-------------

``__eq__``, ``__hash__`` and ``canonical_bytes`` answer for what the molecule is **storing**, so on
their own they mean "same drawing" — see :doc:`molecule`. ``canonicalize()`` is what turns them into
"same compound", which makes deduplicating a corpus one call per molecule and then a set:

.. testcode::

    mols = [smiles('Oc1ccccn1'), smiles('O=c1cccc[nH]1'), smiles('OC1=NC=CC=C1')]
    print(len(set(mols)))

    for m in mols:
        m.canonicalize()
    print(len(set(mols)), mols[0])

.. testoutput::

    3
    1 C=1C(=O)NC=CC=1

``str(mol)`` after the pass is that same key as a string, for a corpus that has to be keyed across
processes; ``canonical_bytes`` is it as bytes, and is what the container hashes.

The pass has to run on **every** side of the comparison. Canonicalizing one and not the other compares a
canonical form against a drawing, which is not a weaker answer but a wrong one.
