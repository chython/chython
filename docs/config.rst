Configuration & Integrations
============================

Global settings, RDKit interoperability, 3D conformers, and pandas support.


Configuration Reference
-----------------------

.. testcode::

    import chython

    # 2D layout engine — validated at assignment
    chython.clean2d_engine = 'smilesdrawer'  # default
    # Options: 'rdkit', 'smilesdrawer', 'cdk', 'obabel', 'indigo'

    # 3D conformer engine
    chython.conformer_engine = 'rdkit'  # default
    # Options: 'rdkit', 'cdpkit'

    # Java JAR paths (CDK, OPSIN); None means read CDK_PATH / OPSIN_PATH from environment
    chython.class_paths = ['/path/to/cdk.jar', '/path/to/opsin.jar']


RDKit Interoperability
-----------------------

``from chython.interop import rdkit`` gives a single callable that dispatches on its argument:
pass a chython molecule to export, pass an RDKit molecule to import.  Every toolkit has one such
callable — ``indigo``, ``openbabel``, ``cdk``, ``cdpkit`` and ``iupac`` beside it — and the export
direction of each is also a method on the container.

.. testcode::
    :skipif: __import__('importlib').util.find_spec('rdkit') is None

    from chython import smiles
    from chython.interop import rdkit

    mol = smiles('c1ccccc1')

    # Export to RDKit, as a function or as the method that calls it
    rdkit_mol = rdkit(mol)
    rdkit_mol = mol.to_rdkit()

    # Import from RDKit.  No method: there is no chython container to hang it on yet
    mol_back = rdkit(rdkit_mol)

An import's records land on the container it produced: ``mol_back.log`` says what the conversion
clamped or dropped, in stage ``'interop'``, with nothing passed in.  Each molecule of an imported
reaction keeps its own, and ``rxn.log.by_subject('products[0]')`` reads the same records from the
reaction end.  The export direction keeps a ``log=`` list instead, since it returns a foreign object
and there is no chython container for its records to go on.

A reaction converts as an ``rdChemReactions.ChemicalReaction``, with the three sides in the three
template lists and the atom–atom mapping carried both ways.  RDKit's atom map field holds one
integer, and the two things chython could put in it are separate flags: ``keep_mapping`` (on by
default) writes ``map_number``, so an unmapped molecule exports with an empty field, while
``keep_numbers`` writes chython's stable atom ids instead — a label to match results back on, not a
mapping.

.. testcode::
    :skipif: __import__('importlib').util.find_spec('rdkit') is None

    from rdkit.Chem import MolToSmiles
    from rdkit.Chem.rdChemReactions import ReactionToSmiles

    from chython import smiles

    print(MolToSmiles(smiles('c1ccccc1O').to_rdkit()))
    print(ReactionToSmiles(smiles('[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]').to_rdkit()))

.. testoutput::
    :skipif: __import__('importlib').util.find_spec('rdkit') is None

    Oc1ccccc1
    [CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]

The other five converters are molecule-only; ``mol.to_indigo()``, ``mol.to_openbabel()``,
``mol.to_cdk()``, ``mol.to_cdpkit()`` and the ``mol.iupac`` property (openclatura, Python ≥ 3.11) are
their methods.  ``rxn.to_rdkit()`` is the only reaction method, because RDKit is the only one of the
five with a reaction form.


3D Conformers
-------------

``generate_conformers`` stores the generated geometry as models of the molecule and returns how many
landed.  The generated set replaces any models the molecule carried; the layout is untouched.  The
engine is controlled by ``chython.conformer_engine``.

.. testcode::
    :skipif: __import__('importlib').util.find_spec('rdkit') is None

    from chython import smiles
    from chython.interop.conformers import generate_conformers

    mol = smiles('CCO')

    stored = generate_conformers(mol, limit=2)
    # one model per stored conformer, and nothing generated came from a file
    print(stored == len(mol.conformers), mol.conformer(0).ext_index)

.. testoutput::

    True None


Pandas Integration
------------------

.. testcode::
    :skipif: __import__('importlib').util.find_spec('pandas') is None

    import pandas as pd
    from chython import smiles, patch_pandas

    # Call once to enable molecule display in DataFrames
    patch_pandas()

    df = pd.DataFrame({
        'mol': [smiles('CCO'), smiles('c1ccccc1')],
        'name': ['ethanol', 'benzene'],
    })
    # Molecules display correctly in DataFrame
