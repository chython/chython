Configuration & Integrations
============================

Global settings, RDKit interoperability, 3D conformers, and pandas support.


Configuration Reference
-----------------------

.. code-block:: python

    import chython

    # 2D layout engine
    chython.clean2d_engine = 'smilesdrawer'  # default
    # Options: 'rdkit', 'smilesdrawer', 'cdk', 'obabel', 'indigo'

    # 3D conformer engine
    chython.conformer_engine = 'rdkit'  # default
    # Options: 'rdkit', 'cdpkit'

    # Neural AAM device (set before first reset_mapping call)
    chython.torch_device = 'cpu'  # default; 'cuda:0' for GPU

    # Java JAR paths (CDK, OPSIN)
    chython.class_paths = ['/path/to/cdk.jar', '/path/to/opsin.jar']
    # Or via environment variables: CDK_PATH, OPSIN_PATH


RDKit Interoperability
-----------------------

.. code-block:: python

    from chython import smiles, MoleculeContainer

    mol = smiles('c1ccccc1')

    # Convert to RDKit Mol
    rdkit_mol = mol.to_rdkit()

    # Convert from RDKit Mol
    mol_back = MoleculeContainer.from_rdkit(rdkit_mol)


3D Conformers
-------------

.. code-block:: python

    import chython
    chython.conformer_engine = 'rdkit'  # or 'cdpkit'

    mol = smiles('CCO')

    # Generate conformer (stored internally)
    # Access via mol._conformers after generation

    # Write 3D coordinates to SDF
    from chython import SDFWrite
    with SDFWrite('output_3d.sdf') as writer:
        writer.write(mol, write3d=0)  # conformer index


Pandas Integration
------------------

.. code-block:: python

    import pandas as pd
    from chython import smiles, patch_pandas

    # Call once to enable molecule display in DataFrames
    patch_pandas()

    df = pd.DataFrame({
        'mol': [smiles('CCO'), smiles('c1ccccc1')],
        'name': ['ethanol', 'benzene'],
    })
    # Molecules display correctly in DataFrame
