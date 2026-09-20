Salts
=====

A record with several components is the normal case, not the exception: a hydrate, an amine
hydrochloride, a sodium carboxylate, a solvate a crystallization left behind, and a reagent drawn beside
the compound it acted on all arrive as one record of several parts. This page is the surface that reads
them.

Three entry points, and **none of them deletes a component**:

=============================== ============================================================ ====================
surface                         does                                                         returns
=============================== ============================================================ ====================
``split_salts()``               cuts the ionic bond, moving the charge onto the two ends      ``bool``
``decompose_salts()``           reads the record as one featurized row per component          ``SaltComposition``
``SaltComposition.compose()``   puts a chosen selection back together                        ``MoleculeContainer``
=============================== ============================================================ ====================

The first edits, the second only reports — it **changes nothing and logs nothing**, the return value
being the whole report. A fourth surface is not on this page: the charges a salt was drawn without are
moved by a ``standardize()`` stage, and :doc:`standardize` documents it beside the other stages.


The Table
---------

``chython/chemistry/tables/salts.tsv`` holds all of it: 134 rows, of which 110 name a whole species. Two
columns carry two different questions — ``match`` is *how* a row is recognized, ``klass`` is *what* the
row is — and a pass selects rows by ``klass`` and never by ``match``, which is what makes a new class
invisible to a pass that did not name it.

========= ======== ====================================================================================
``match`` pattern  recognized by
========= ======== ====================================================================================
``embed`` SMARTS   embedding; the ``:1`` atom is the one the row is about
``whole`` SMILES   canonical-key equality against the whole component
========= ======== ====================================================================================

Four classes are about an atom and the other fourteen name a species:

=================== ======================================================================================
group               classes
=================== ======================================================================================
atom rows           ``metal_cation``, ``charge_acceptor``, ``protic_acid``, ``metal_protic``
acids               ``mineral_acid``, ``sulfonic_acid``, ``short_carboxylic_acid``, ``carboxylic_acid``,
                    ``aromatic_acid``, ``fatty_acid``
bases               ``amine_base``, ``amino_acid``
solvents            ``water``, ``alcohol``, ``hydrocarbon``, ``halo_solvent``, ``aprotic_solvent``
cations             ``quaternary_ammonium``
=================== ======================================================================================

``SALT_CLASSES`` and ``SALT_MATCHES`` are that vocabulary at runtime, so a caller comparing a row's
``klass`` against a name it typed can check the name exists.

``metal_protic`` is the one class neither pass on this page reads. Water and an alcohol are sites only for
a free s-block metal — the :doc:`standardize` stage's question and nobody else's — so ``decompose_salts()``
still reads ``CCN.CCO`` as a solvate, a hydrate stays droppable, and neither is an ``acid_salt``. Sulfur is
an acid site: a thiol, a thiophenol, a thioacid and a xanthic acid each have a ``protic_acid`` row, so an
amine beside one is the ``acid_salt`` it is on the shelf.

The acid, base and solvent classes are the **formers** — the species a salt is *made with*, which is the
word the tags below use. ``quaternary_ammonium`` and ``metal_cation`` are not: they are the cation half of
a pair rather than something added to it.

A species row is compared by the canonical SMILES with stereo disabled, because a row names a
**constitution**: whether the tartrate is L, D or meso does not change that it is the counterion. And a
**conjugate is not a row** — the probe below neutralizes before it matches, so ``[Cl-]`` is ``Cl`` and
``CC([O-])=O`` is acetic acid by the time the key is taken. Only a species that cannot be neutral gets a
row of its own: ``salts:choline`` is a ``quaternary_ammonium`` and has nowhere to put the charge.

One more consequence of ``klass`` being the selector: **a species that is both solvent and base is
tabulated as the base.** ``salts:pyridine`` is an ``amine_base``, because as a solvent class it would
read pyridine hydrochloride as hydrochloric acid, while as a base that record answers with both
components as candidate parents — which a caller can reject.


Splitting the Ionic Bond
------------------------

``split_salts()`` cuts a metal–heteroatom bond and moves the charge onto the two ends. The atom count
does not change, and all 93 metals the core's ``[M]`` accepts are in scope:

.. testcode::

    from chython import smiles, DEFAULT_STABILIZER_CLASSES

    mol = smiles('CC(=O)O[Na]')
    print(mol.split_salts(), mol)

.. testoutput::

    True C(C)([O-])=O.[Na+]

The test is **all-or-nothing per cation atom**, which is what makes generality over 93 metals safe: a
dative bond, a neighbour that matches no acceptor row, an untabulated resulting charge or an implicit
hydrogen on the cation refuses the whole atom and logs why. Cisplatin, ferrocene and the metal carbonyls
therefore come back intact rather than half-split. A dative bond is the *exemption* signal and never the
trigger — ``standardize()`` installs it to record coordination that must be preserved:

.. testcode::

    mol = smiles('N[Pt](N)(Cl)Cl')
    print(mol.split_salts(), mol)
    print(mol.log.refused()[0].rule)

.. testoutput::

    False N[Pt](Cl)(N)Cl
    salts:metal

**A drawn metal–oxygen bond states the alkoxide**, so ``salts:alkoxide`` reads the saturated case the
carboxylate and phenolate rows do not cover. Silicon and boron are outside ``[M]``, which is why a silyl
ether and a borate ester are not salts:

.. testcode::

    for s in ('CCO[Na]', 'CC(C)(C)O[K]', 'CCS[K]', 'C[Si](C)(C)O[Na]', 'CCO[Si](C)(C)C', 'CCOCC'):
        mol = smiles(s)
        print(mol.split_salts(), mol)

.. testoutput::

    True C(C)[O-].[Na+]
    True C(C)([O-])(C)C.[K+]
    True C(C)[S-].[K+]
    True C[Si]([O-])(C)C.[Na+]
    False C(C)O[Si](C)(C)C
    False C(C)OCC

That is a different question from ``CCO.[Na]``, where nothing is drawn between the two components: there
the :doc:`standardize` stage decides, and it decides on the metal.


Reading the Record
------------------

``decompose_salts()`` answers what the record is made of. It works on a copy that is
hydrogen-implicified, salt-split, neutralized and aromatized, which is what makes the answer comparable
across a corpus rather than across a drawing — the four ways of writing sodium acetate report the same
components in the same roles:

.. testcode::

    for s in ['CC(=O)O[Na]', 'CC(=O)[O-].[Na+]', 'CC(=O)O.[Na+]', 'CC(=O)O.[Na]']:
        r = smiles(s).decompose_salts()
        print(len(r.parents), [row.is_lone_metal for row in r.parents].count(True), r.stabilizers)

.. testoutput::

    2 1 ()
    2 1 ()
    2 1 ()
    2 1 ()

What converges is what each component *is* — molecule, species, class and role. ``charge`` is the charge
as drawn and keeps the drawing's own answer, and so do the tags.


The row
~~~~~~~

Each component is one ``ComponentRow``:

===================== =======================================================================
field                 is
===================== =======================================================================
``atoms``             the atom ids in the **caller's** molecule, a back-pointer
``molecule``          the **normalized** component
``species``           the row id it matched, ``salts:water``, or ``None``
``klass``             that row's class, or ``None``
``equivalents``       how many rows share this canonical key **and** this role
``heavy_atoms``       count
``carbon_count``      count
``ring_count``        count
``charge``            the charge as drawn
``residual_charge``   the charge left after ``neutralize(keep_charge=False)``
``is_lone_metal``     the component is one metal atom, whatever its charge
``is_organometallic`` a carbon bonded to a metal, or to boron
``role``              ``'parent'`` or ``'stabilizer'``
===================== =======================================================================

.. testcode::

    r = smiles('CC(=O)Oc1ccccc1C(=O)O.O').decompose_salts()     # aspirin monohydrate
    for row in r.components:
        print(row.role, row.species, row.klass, row.heavy_atoms, row.carbon_count, row.ring_count)

.. testoutput::

    parent None None 13 9 1
    stabilizer salts:water water 1 0 0

``molecule`` and ``atoms`` are the answer to "which one is it": the normalized copy costs chemistry and
not identity, so a caller that needs comparability has the form that compares, and one that needs its own
drawing has the ids to find it. **One drawn component can yield two rows**, because the probe splits the
ionic bond first, and then the two rows' ``atoms`` partition the one component the caller drew:

.. testcode::

    mol = smiles('CC(=O)O[Na]')                                 # one component as drawn
    print(len(mol.split()), [(row.atoms, str(row.molecule)) for row in mol.decompose_salts().components])

.. testoutput::

    1 [((1, 2, 3, 4), 'C(C)(=O)O'), ((5,), '[Na+]')]

``is_lone_metal`` is **charge-blind**, and that is load-bearing: a neutral ``[Na]`` has residual charge 0,
so keyed on charge it would read as droppable and sodium metal would be stripped off an acid.
``equivalents`` counts key **and** role, so a chloride on counter-ion duty and a free HCl copy in the same
record each count only their own kind:

.. testcode::

    r = smiles('C[N+](C)(C)C.[Cl-].Cl').decompose_salts()
    for row in r.components:
        print(row.role, row.species, row.equivalents, row.charge, row.residual_charge)
    print(r.equivalents_by_species())

.. testoutput::

    parent salts:hcl 1 -1 0
    stabilizer salts:hcl 1 0 0
    parent None 1 1 1
    {'salts:hcl': 1}

``equivalents_by_species()`` regroups the **stabilizer** rows only, which is the stoichiometry a name like
"hydrochloride monohydrate" states, and it omits an untabulated row rather than keying it on ``None``.


Roles
~~~~~

``components`` splits into ``parents`` and ``stabilizers``. **A stabilizer is eligible, never leftover:**
a component leaves only when its residual charge is zero — the charge it still carries after
``neutralize(keep_charge=False)``, which separates a protonation state from an intrinsic charge — and it
is neither a lone metal nor organometallic nor on counter-ion duty, and its class is one ``classes``
names. The default is narrow on purpose: water, the mineral acids, the sulfonic acids and the C1–C2
carboxylic acids, which are the species a hydrate or an ``-HCl`` is made of and nothing else.

.. testcode::

    r = smiles('CC(=O)Oc1ccccc1C(=O)O.O').decompose_salts()     # aspirin monohydrate
    print([str(row.molecule) for row in r.parents], r.equivalents_by_species())
    print([str(row.molecule) for row in smiles('CC(=O)O.O').decompose_salts().parents])
    print([str(row.molecule) for row in smiles('C[N+](C)(C)C.[Cl-]').decompose_salts().parents])

.. testoutput::

    ['O=C(Oc1c(C(O)=O)cccc1)C'] {'salts:water': 1}
    ['C(C)(=O)O', 'O']
    ['Cl', 'C[N+](C)(C)C']

The second line is the one guard: the parents must contain a row that is neither a lone metal nor a
recognized solvent, and when they do not, every component is a parent — acetic acid is itself one of the
C1–C2 acids, so acetic acid monohydrate answers both its components rather than nothing. The third line is
counter-ion duty, below.

Three keywords widen or narrow what may leave, and nothing else does:

================ ==================================================================================
keyword          is
================ ==================================================================================
``classes``      which ``klass`` values may be a stabilizer
``max_atoms``    heavy-atom ceiling on a stabilizer; ``None`` is no ceiling
``discardable``  species keys, row ids or class names that are eligible on top of ``classes``
================ ==================================================================================

.. testcode::

    salt = smiles('CCN.OC(=O)c1ccccc1')                         # ethylamine benzoate
    print([str(row.molecule) for row in salt.decompose_salts().parents])
    print([str(row.molecule) for row in
           salt.decompose_salts(classes=DEFAULT_STABILIZER_CLASSES + ('aromatic_acid',)).parents])
    print([str(row.molecule) for row in
           smiles('CC(=O)Oc1ccccc1C(=O)O.O').decompose_salts(max_atoms=0).parents])

.. testoutput::

    ['c1(ccccc1)C(=O)O', 'C(C)N']
    ['C(C)N']
    ['O=C(Oc1c(C(O)=O)cccc1)C', 'O']

A benzoate is tabulated, classified and a parent by default: which of an amine and an aromatic acid is
the compound is a question about a collection, so widening it is the caller's decision and never the
pass's.

**A widening needs no caller-side pre-check.** A registry that treats solvent of crystallization as
packaging rather than substance widens by the four solvent classes, and needs no "is this record nothing
but solvent?" test in front of it:

.. testcode::

    wide = DEFAULT_STABILIZER_CLASSES + ('alcohol', 'hydrocarbon', 'halo_solvent', 'aprotic_solvent')
    for s in ['CC(=O)Oc1ccccc1C(=O)O.Cc1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O.CS(C)=O.CCO',
              'Cc1ccccc1', 'CCO.O']:
        r = smiles(s).decompose_salts(classes=wide)
        print([str(row.molecule) for row in r.parents], sorted(r.tags))

.. testoutput::

    ['O=C(Oc1c(C(O)=O)cccc1)C'] ['acid_salt', 'solvate']
    ['O=C(Oc1c(C(O)=O)cccc1)C'] ['acid_salt', 'competing_formers', 'solvate']
    ['c1c(C)cccc1'] ['single', 'stabilizer_only']
    ['C(C)O', 'O'] ['competing_formers', 'hydrate', 'solvate', 'stabilizer_only']

The last two lines are the one guard again: a solvent name in a reagent field is a compound to whoever
wrote it, so a record left with no parent promotes every component back and says so with
``stabilizer_only``. The widening cannot empty a record, however wide it is.


Counter-ion duty
~~~~~~~~~~~~~~~~

Two clauses keep a component that would otherwise be eligible, both of them about a charge the rest of the
record needs.

**An intrinsic charge keeps its counter-ion.** A component whose residual charge is zero but whose drawn
charge opposes the record's intrinsic charge stays, up to that many equivalents of it. Residual charge
alone would let a chloride leave a quaternary ammonium; charge balance alone would take the water off a
sodium sulfonate monohydrate together with the sulfonate. That is the ``C[N+](C)(C)C.[Cl-]`` line above:
the record's intrinsic charge is ``+1``, the chloride is drawn ``-1``, so it stays and the pair is an
``ion_pair``.

**A lone metal that lacks its anion** puts every anion and every neutral acid in the record on that duty,
all of them and with no budget: the missing charge is the one quantity the record does not state, so there
is nothing to count equivalents against. The metal lacks it two ways — drawn neutral it states no charge
at all, and in a record whose drawn charges do not balance, in either direction, it carries one the
drawing does not account for. The second way is why the answer does not turn on whether ``standardize()``
ran first: ``[Mg+]`` beside two acetates is a magnesium the salt-charge stage itself calls under-charged.
Water is on no duty, being neither an anion nor an acidic site, so the hydrate below still loses its
water; and a balanced sodium chloride puts nothing on duty at all, so the acetic acid beside it is a
stabilizer:

.. testcode::

    for s in ['Cl.[Na].CCBr', 'CC(=O)O.[Na+].CCBr', 'CC(=O)[O-].CC(=O)[O-].[Mg+].CCBr',
              'CC(=O)O.[Na].O.CCBr', '[Na+].[Cl-].CC(=O)O.CCBr']:
        r = smiles(s).decompose_salts()
        print([str(row.molecule) for row in r.parents], r.equivalents_by_species())

.. testoutput::

    ['[Na]', 'Cl', 'C(C)Br'] {}
    ['C(C)(=O)O', 'C(C)Br', '[Na+]'] {}
    ['C(C)(=O)O', 'C(C)(=O)O', 'C(C)Br', '[Mg+]'] {}
    ['C(C)(=O)O', '[Na]', 'C(C)Br'] {'salts:water': 1}
    ['Cl', 'C(C)Br', '[Na+]'] {'salts:acetic': 1}


Tags
~~~~

``tags`` says what the record **is** — never what this pass did with it, which is what ``role`` says. A
record carries every tag that applies:

===================== =========================================================================
tag                   the record
===================== =========================================================================
``single``            has one component
``mixture``           has several parents, none a former, a lone metal or on duty
``hydrate``           has more than one component, one of them water drawn uncharged
``solvate``           has more than one component, one of them an uncharged non-water solvent
``acid_salt``         holds an acid or an acidic site beside a component that is not a lone metal
``base_salt``         holds a tabulated base beside an acidic site
``metal_salt``        pairs a lone-metal parent with an anion or a neutral acid
``elemental_metal``   holds a lone metal drawn neutral with nothing to pair with
``ion_pair``          has a parent still charged after neutralization: the charge is intrinsic
``competing_formers`` names more than one former species
``charge_unbalanced`` has drawn charges that do not sum to zero
``charges_undrawn``   draws a lone metal neutral beside its anion: the salt-charge stage repairs it
``stabilizer_only``   holds nothing that is a compound in its own right, so every component is a parent
===================== =========================================================================

.. testcode::

    for s in ['O', '[Na]', 'CCN.Cl', 'CCN(CC)CC.Cl', 'CC(=O)O[Na]', 'CC(=O)O.[Na]', 'CCBr.CCO',
              'CCBr.CCI', 'C[N+](C)(C)C.[Cl-]', '[Na+].[Cl-].[Cl-].CCBr', '[OH-].[Na+]',
              'CC(C)(C)[O-].[K+]']:
        print(s, sorted(smiles(s).decompose_salts().tags))

.. testoutput::

    O ['single', 'stabilizer_only']
    [Na] ['elemental_metal', 'single', 'stabilizer_only']
    CCN.Cl ['acid_salt']
    CCN(CC)CC.Cl ['acid_salt', 'base_salt', 'competing_formers']
    CC(=O)O[Na] ['ion_pair', 'metal_salt']
    CC(=O)O.[Na] ['charges_undrawn', 'metal_salt']
    CCBr.CCO ['solvate']
    CCBr.CCI ['mixture']
    C[N+](C)(C)C.[Cl-] ['acid_salt', 'ion_pair']
    [Na+].[Cl-].[Cl-].CCBr ['acid_salt', 'charge_unbalanced', 'ion_pair', 'metal_salt']
    [OH-].[Na+] ['ion_pair', 'metal_salt', 'stabilizer_only']
    CC(C)(C)[O-].[K+] ['ion_pair', 'metal_salt', 'stabilizer_only']

**The two solvent tags read the drawn charge, where ``klass`` reads the neutralized form.** The two last
lines are why: a conjugate is not a row, so hydroxide keys as ``water`` and an alkoxide as ``alcohol`` —
right for identity, and the wrong question for a tag that claims solvent of crystallization. Sodium
hydroxide and potassium *tert*-butoxide are ``metal_salt`` and nothing more.

Two more lines are worth reading twice. ``CCBr.CCO`` is the tag and the role saying different things,
which is the design: the record *is* a solvate, and ethanol is an ``alcohol``, a class the default does
not let leave — so both components are parents, and the tag is what tells a caller a solvate is what it is
looking at. And ``CCN.Cl`` is an
``acid_salt`` without being a ``base_salt``: the acid half is recognized from an acidic site, which every
drawing states, while the base half is recognized from a tabulated species, and ethylamine is not one
where triethylamine is.


Putting a Selection Back
------------------------

``compose()`` reassembles a chosen selection, over ``union()``, which preserves every parity and stereo
group the components carried — ``substructure()`` would drop them, a parity being stated in a frame that
cutting bonds destroys:

.. testcode::

    r = smiles('C[C@H](N)C(=O)O.O').decompose_salts()
    print(r.compose(r.parents))

.. testoutput::

    C([C@H](C)N)(=O)O

It takes rows of this composition and raises ``ValueError`` on an empty selection or a row from another
record, whose atom numbers this one never issued. What it returns is built from the **normalized**
components, so it is a key rather than a reproduction of the drawing: an ion pair recomposes with its
counter-ion neutral.


Linking Records Across Sources
------------------------------

Three sources describe one compound three ways: one delivers a single neutral molecule and states the
salt in metadata, the next delivers the salt, the third a hydrate of it. Linking them means deciding, per
record, which part is the compound — and the roles above are that decision already made, so the key is
two lines:

.. testcode::

    def parent_key(mol):
        r = mol.decompose_salts()
        parent = r.compose(r.parents)
        parent.canonicalize()
        return str(parent)

    for s in ['NCC(=O)O', 'NCC(=O)O.Cl', 'NCC(=O)O.OC(=O)C(F)(F)F.O', 'CC(=O)Oc1ccccc1C(=O)O.O']:
        print(parent_key(smiles(s)))

.. testoutput::

    C(CN)(=O)O
    C(CN)(=O)O
    C(CN)(=O)O
    O=C(Oc1c(C(O)=O)cccc1)C

Free glycine, its hydrochloride and its TFA salt monohydrate link; aspirin monohydrate links to aspirin.
The stoichiometry is not lost by linking on the parent — ``equivalents_by_species()`` is the rest of the
record, so a link can carry ``{'salts:tfa': 1, 'salts:water': 1}`` as an attribute of the association
rather than of the compound.

**A metal salt does not link to its free acid, and an ion pair does not link to its ion.** Neither needs a
clause: a lone metal is never eligible, so it is a parent and travels in the key, and an anion on
counter-ion duty travels with the cation it balances.

.. testcode::

    for s in ['CC(=O)[O-].[Na+]', 'CC(=O)O[Na]', 'CC(=O)O.[Na]', 'CC(=O)O.[Na+]']:
        print(parent_key(smiles(s)))
    print(parent_key(smiles('C[N+](C)(C)C.[Cl-]')))

.. testoutput::

    C(C)([O-])=O.[Na+]
    C(C)([O-])=O.[Na+]
    C(C)([O-])=O.[Na+]
    C(C)([O-])=O.[Na+]
    C[N+](C)(C)C.Cl

All four spellings of sodium acetate reach one key — the covalent one, the one that states no charge and
the one whose charges do not balance included — and tetramethylammonium chloride keeps both halves rather
than becoming a cation. The charge repair is
``canonicalize()``'s, not this page's: ``decompose_salts()`` reports a mis-drawn record faithfully and
tags it ``charges_undrawn``, and a pipeline is what repairs it.

The other question a corpus asks is which records are a **compound** at all and which are a reagent that
drove a reaction. ``tags`` and ``klass`` answer it, and the class names a caller counts as a driver are
the caller's list — a collection of amine syntheses and a collection of amine reactions do not draw that
line in the same place:

.. testcode::

    DRIVERS = frozenset({'mineral_acid', 'sulfonic_acid', 'short_carboxylic_acid', 'water', 'alcohol',
                         'hydrocarbon', 'halo_solvent', 'aprotic_solvent', 'amine_base'})

    def is_driver(mol):
        r = mol.decompose_salts()
        return 'elemental_metal' in r.tags or all(row.klass in DRIVERS for row in r.components)

    for s in ['Cl', 'CCOCC', 'CCN(CC)CC', '[Na]', 'ClCCl.Cl', 'NCC(=O)O.Cl',
              'CC(=O)Oc1ccccc1C(=O)O.O']:
        print(s, is_driver(smiles(s)))

.. testoutput::

    Cl True
    CCOCC True
    CCN(CC)CC True
    [Na] True
    ClCCl.Cl True
    NCC(=O)O.Cl False
    CC(=O)Oc1ccccc1C(=O)O.O False

``elemental_metal`` is the one tag that has to be named: sodium metal is a lone metal, so it is a parent
and no class covers it. Everything else falls out of a row being tabulated and classified — triethylamine
is an ``amine_base``, dichloromethane a ``halo_solvent``. Glycine hydrochloride and aspirin monohydrate
each hold a component **no row matches at all**, whose ``klass`` is ``None``, and that is the shape of a
compound: a table of 110 species names what a compound is beside, never the compound.

What is deliberately absent from both recipes is a threshold: an element allow-list, a cap on how many
components a record may have, a size above which a component stops being a counterion. Each is a claim
about one collection rather than a fact about a molecule, so each belongs in the caller's tree and not in
this pass. :doc:`substructure` and the screening section of :doc:`standardize` are where the rest of that
tree is built.
