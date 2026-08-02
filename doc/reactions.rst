Reactions & Templates
=====================

Parsing, CGR, atom-atom mapping, reaction templates, functional groups, and deprotection.


Parsing Reactions
-----------------

.. code-block:: python

    from chython import smiles

    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')

    # Access components
    rxn.reactants   # tuple of MoleculeContainers
    rxn.products    # tuple of MoleculeContainers
    rxn.reagents    # tuple of MoleculeContainers

    # Iterate all molecules
    for mol in rxn.molecules():
        print(str(mol))

    # Metadata
    rxn.name = 'Amination'
    rxn.meta['temperature'] = '100'


Reaction Signatures
-------------------

Canonical reaction SMILES (SMIRKS) with molecules sorted in canonical order:

.. code-block:: python

    str(rxn)             # canonical reaction SMILES
    format(rxn, 'm')     # with atom mapping

    # Format specifiers (same as molecule, plus extras)
    format(rxn, '!c')    # keep original molecule order (no sorting)
    format(rxn, '!C')    # skip CXSMILES fragment contraction

Reactions are hashable and comparable, same as molecules:

.. code-block:: python

    rxn1 == rxn2         # True if same canonical signature
    {rxn1, rxn2}         # set deduplication


Reaction Standardization
------------------------

ReactionContainer has the same standardization methods as molecules.
They are applied to all molecules in the reaction:

.. code-block:: python

    rxn.canonicalize()
    rxn.standardize()
    rxn.kekule()
    rxn.thiele()
    rxn.neutralize()
    rxn.explicify_hydrogens()   # tries to preserve atom-atom mapping
    rxn.implicify_hydrogens()

Reaction-specific methods:

.. code-block:: python

    # Move unchanged reactants to reagents (based on atom-atom mapping)
    rxn.remove_reagents(keep_reagents=True)

    # Merge ions into single multicomponent molecules
    rxn.contract_ions()

Example workflow:

.. code-block:: python

    rxn = smiles('[Na+:1].[OH-:2].[CH3:7][O:5][C:4]([CH3:3])=[O:6]>>[CH3:3][C:4]([OH:8])=[O:6]')
    rxn.contract_ions()       # merge Na+ and OH- into NaOH
    rxn.remove_reagents(keep_reagents=True)  # NaOH/MeOH become reagents


Condensed Graph of Reaction (CGR)
----------------------------------

CGR overlays reactant and product graphs, showing bond changes:

.. code-block:: python

    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')

    # Compose CGR
    cgr = ~rxn             # shorthand
    cgr = rxn.compose()    # same

    # From two mapped molecules
    mol_r = smiles('[CH3:1][OH:2]')
    mol_p = smiles('[CH3:1][NH2:3]')
    cgr = mol_r.compose(mol_p)

    # Reaction center
    center = cgr.center_atoms  # tuple of atom numbers where bonds change

    # Extract center substructure
    center_cgr = cgr.substructure(center)
    aug_cgr = cgr.augmented_substructure(center, deep=1)

    # CGR supports isomorphism search
    for mapping in cgr_query.get_mapping(cgr):
        print(mapping)


Atom-Atom Mapping
-----------------

Neural attention-based mapping (requires ``chython-rxnmap`` package):

.. code-block:: python

    from chython import smiles

    rxn = smiles('CCO.CC(=O)O>>CCOC(=O)C.O')

    # Neural attention mapping (ONNX-based, CPU only)
    rxn.attention_mapping()

    # With score (float, higher = more confident)
    score = rxn.attention_mapping(return_score=True)

    # Keep original reactant atom numbers
    rxn.attention_mapping(keep_reactants_numbering=True)

``attention_mapping`` loads the ONNX model once on first call.

Utility methods:

.. code-block:: python

    # Reset mapping: deduplicate atom numbers across components (no model)
    rxn.reset_mapping()

    # Rule-based fix for known mapping mistakes (called automatically by attention_mapping)
    rxn.fix_mapping()
    log = rxn.fix_mapping(logging=True)


Reactor (Multi-Reactant Templates)
------------------------------------

Reactor applies SMARTS-pattern transformations to one or more reactant molecules.
Atom numbers in query patterns and product templates must be mapped to each other.

.. code-block:: python

    from chython import smarts, smiles, Reactor

    # Define patterns using mapped SMARTS
    acid = smarts('[C:1]([O;D1:2])=[O:3]')    # carboxylic acid
    alco = smarts('[C:4][O;D1:5]')             # alcohol

    # Product template (reuses atom mapping numbers from patterns)
    ester = smarts('[C:1](=[O:3])[O:5][C:4]')

    # Create reactor (patterns and products are tuples)
    reactor = Reactor((acid, alco), (ester,),
                      delete_atoms=True,   # remove atoms not in product
                      one_shot=False)      # apply multiple times if possible

    # Apply to reactants (pass as positional args, not a list)
    acid_mol = smiles('CC(=O)O')
    alco_mol = smiles('CCO')
    for product in reactor(acid_mol, alco_mol):
        print(str(product))

Use ``one_shot=True`` to apply the template only once (first match).


Transformer (Single Molecule)
------------------------------

Transformer applies a pattern replacement to a single molecule:

.. code-block:: python

    from chython import smarts, smiles, Transformer

    pattern = smarts('[C:1]=[O:2]')
    replacement = smarts('[C:1][O:2]')

    t = Transformer(pattern, replacement)

    mol = smiles('CC=O')
    for result in t(mol):
        print(str(result))


Predefined Reactors
--------------------

The ``@`` operator enumerates possible reactions between molecules using predefined
Reactor templates built on functional group detection. It returns a generator
of ``(reaction_name, ReactionContainer)`` tuples:

.. code-block:: python

    from chython import smiles

    acid = smiles('CC(=O)O')
    amine = smiles('CCN')

    # Two-component reaction
    for name, rxn in acid @ amine:
        print(name, rxn)
    # amidation  CCN.O=C(O)C>>CCNC(=O)C

The operator is symmetric — order does not matter:

.. code-block:: python

    list(amine @ acid) == list(acid @ amine)  # True

Multi-component reactions use a list:

.. code-block:: python

    aldehyde = smiles('CC=O')
    amine = smiles('CCN')
    isocyanide = smiles('[C-]#[N+]C')

    for name, rxn in aldehyde @ [amine, isocyanide]:
        print(name, rxn)
    # ugi_3cr  C(=O)C.CCN.[C-]#[N+]C>>CCNC(C(=O)NC)C

Each molecule's ``functional_groups`` property is used as a pre-filter to select
applicable Reactor templates. Only matching combinations are attempted.

Selective reaction application with the ``reaction`` keyword:

.. code-block:: python

    arx = smiles('Brc1ccccc1')
    boronic = smiles('OB(O)c1ccccc1')

    # Only Suzuki coupling (skip other possible reactions)
    for name, rxn in arx.react(boronic, reaction='suzuki'):
        print(name, rxn)


Oxidation, Reduction & Transformation
---------------------------------------

Single-molecule transformations are available via dedicated methods:

.. code-block:: python

    mol = smiles('OCC')
    mol.canonicalize()

    # Oxidation products
    for name, rxn in mol.oxidize():
        print(name, rxn)
    # alcohol_to_aldehyde  CCO>>C=O

    # Reduction products
    ketone = smiles('CC(=O)c1ccccc1')
    ketone.canonicalize()
    for name, rxn in ketone.reduce():
        print(name, rxn)
    # ketone_to_alcohol  c1ccccc1C(C)=O>>OC(C)c1ccccc1

    # Functional group interconversions (Appel, borylation, ring closures, etc.)
    for name, rxn in mol.transform():
        print(name, rxn)
    # appel  CCO>>CCBr

    # All single-molecule transformations at once (~ operator)
    for name, rxn in ~mol:
        print(name, rxn)

All methods accept an optional ``reaction`` keyword to apply selectively:

.. code-block:: python

    mol = smiles('OC(C)c1ccccc1')
    mol.canonicalize()

    # Only oxidize to ketone (skip other possible transformations)
    for name, rxn in mol.oxidize(reaction='alcohol_to_ketone'):
        print(name, rxn)

    # Only Appel (alcohol → bromide)
    for name, rxn in mol.transform(reaction='appel'):
        print(name, rxn)


Functional & Protective Groups
-------------------------------

Detect functional groups and their counts:

.. code-block:: python

    mol = smiles('Clc1ccc(Br)cc1O')
    mol.canonicalize()

    mol.functional_groups
    # {'aryl_chloride': 1, 'aryl_bromide': 1, 'phenol': 1}

Detect protective groups:

.. code-block:: python

    mol = smiles('CC(NC(=O)OC(C)(C)C)CNC(=O)OC(C)(C)C')
    mol.canonicalize()

    mol.protective_groups
    # {'amine_boc': 2}

Remove protective groups (in-place). ``remove_protection`` takes a **list** of
group names (or ``None`` for all) and returns ``True`` if anything changed:

.. code-block:: python

    mol = smiles('c1ccccc1NC(=O)OC(C)(C)C')
    mol.canonicalize()

    # Remove a specific set of protective groups (pass a list, not a string)
    mol.remove_protection(['amine_boc'])  # returns True if changed
    str(mol)  # 'c1ccccc1N'

    # Remove all known protective groups
    mol2 = smiles('CC(NC(=O)OC(C)(C)C)COC(=O)OCC=C')
    mol2.canonicalize()
    mol2.remove_protection()  # None => removes all found PGs

Passing a list removes exactly that set, leaving other groups in place. Unknown
names raise ``ValueError``:

.. code-block:: python

    mol = smiles('CC(C)(C)OC(=O)Nc1ccc(OC(C)(C)C)cc1')
    mol.canonicalize()
    mol.protective_groups              # {'amine_boc': 1, 'hydroxyl_tbu': 1}

    mol.remove_protection(['amine_boc'])   # strip only the Boc
    mol.protective_groups              # {'hydroxyl_tbu': 1}  (tBu ether kept)

With ``logging=True`` it returns the list of freed (retained) atom numbers
instead of a bool — the heteroatoms revealed by deprotection (e.g. the unmasked
amine nitrogen). These feed the ``masked`` argument of the sticky enumerators
below:

.. code-block:: python

    mol = smiles('O=C(OC(C)(C)C)NCc1ccc(Br)cc1')
    mol.canonicalize()
    freed = mol.remove_protection(logging=True)   # [<amine N number>]
    # no-op deprotection returns [] (the logging analogue of False)


Sticky Fragments & Linkers
--------------------------

Cut a molecule at each detected coupling role and emit the pieces as
concatenation-ready SMILES for combinatorial reassembly. Roles are coupling
handles (``aryl_halide``, ``aryl_acyl``, ``alkyl_amine``, ...). The molecule must
be a single connected component.

``sticky_fragments()`` enumerates mono-attachment pieces — one open bond each:

.. code-block:: python

    from chython import smiles

    mol = smiles('Brc1ccc(cc1)C(=O)O')
    mol.canonicalize()

    for frag in mol.sticky_fragments():
        print(frag.role, frag.sticky_left, frag.sticky_right)
    # aryl_halide  -c1ccc(C(=O)O)cc1   O=C(O)c1ccc(cc1)-
    # aryl_acyl    -C(=O)c1ccc(Br)cc1  Brc1ccc(cc1)C(=O)-

    # Restrict to one role
    list(mol.sticky_fragments('aryl_halide'))

Each ``StickyFragment`` carries the attachment three ways: ``sticky_left``
(``-...``, glues onto the left of another piece), ``sticky_right`` (``...-``,
glues onto the right), and ``canonical_smiles`` (the ``[At]``-capped fragment,
the dedup key).

``sticky_linkers()`` enumerates bi-attachment pieces — open bonds at both ends.
The left role is capped ``[210At]``, the right ``[211At]``:

.. code-block:: python

    for lk in mol.sticky_linkers('aryl_halide', 'aryl_acyl'):
        print(lk.role_left, lk.role_right)
        print(lk.sticky_left)   # -c1ccc(cc1)C(=O)-   (role_left end first)
        print(lk.sticky_right)  # -C(=O)c1ccc(cc1)-   (flipped orientation)

    # No role arguments: all role pairs
    list(mol.sticky_linkers())

``sticky_left`` / ``sticky_right`` are the same linker traversed from opposite
ends; ``canonical_smiles`` (with ``[210At]``/``[211At]`` caps) is the dedup key.

Both enumerators accept ``masked`` — atom numbers barred from becoming a step-1
attachment point. Feed it the output of ``remove_protection(logging=True)`` so a
group revealed by an incidental deprotection never becomes a standalone handle:

.. code-block:: python

    mol = smiles('O=C(OC(C)(C)C)NCc1ccc(Br)cc1')
    mol.canonicalize()
    freed = mol.remove_protection(logging=True)   # reveal the amine

    # the freed amine is now barred as a fragment attachment point
    roles = {f.role for f in mol.sticky_fragments(masked=freed)}
    # 'alkyl_amine' not in roles

For linkers, ``masked`` bars the atom from the **left** (``[210At]``, step-1) end
only; the right (``[211At]``, deferred step-2) end is exempt. ``masked=None`` and
``masked=()`` behave identically to omitting it.


Reconstruct Mapping
-------------------

Annotate a reaction by trying to reconstruct the product from reactants using
predefined templates. If successful, sets atom-to-atom mapping and returns
matched reaction labels:

.. code-block:: python

    from chython import smiles

    rxn = smiles('Brc1ccccc1.OB(O)c1ccc(F)cc1>>Fc1ccc(-c2ccccc2)cc1')
    rxn.reset_mapping()

    labels = rxn.reconstruct_mapping()
    # ['react:suzuki']

    # Mapping is now set on the product
    format(rxn, 'm')  # reaction SMILES with atom mapping

Supports single-product reactions. Tries in order:

1. Standalone deprotection
2. Standalone protection (reverse)
3. Single-molecule transforms (oxidize/reduce/transform)
4. Deprotection + transform composition
5. Multi-component reactions (subset-based)

Returns an empty list if no template matches:

.. code-block:: python

    rxn = smiles('CC>>CCC')
    rxn.reconstruct_mapping()  # []
