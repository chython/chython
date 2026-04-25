Reactions & Templates
=====================

Parsing, CGR, atom-atom mapping, reaction templates, and deprotection.


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

Neural attention-based mapping (requires ``chytorch-rxnmap`` package):

.. code-block:: python

    import chython
    chython.torch_device = 'cpu'  # set before first use; 'cuda:0' for GPU

    rxn = smiles('CCO.CC(=O)O>>CCOC(=O)C.O')
    rxn.reset_mapping()

    # Rule-based fix for known mapping mistakes
    rxn.fix_mapping()
    log = rxn.fix_mapping(logging=True)

``reset_mapping`` loads the neural model once on first call.
To use GPU, set ``chython.torch_device`` before the first call.
For multiprocessing, call ``reset_mapping`` only inside workers to avoid a single-GPU bottleneck.


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


Deprotection
------------

Built-in templates for ~50+ protective group removals:

.. code-block:: python

    from chython.reactor.deprotection import apply_all, hydroxyl_benzyl, amine_boc

    mol = smiles('...')  # protected molecule

    # Remove specific protection group
    for result in hydroxyl_benzyl(mol):
        print(str(result))

    # Remove all known protection groups iteratively
    result = apply_all(mol)
