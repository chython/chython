Input / Output
==============

Reading and writing molecules and reactions in all supported formats.


String Parsers
--------------

SMILES
~~~~~~

.. code-block:: python

    from chython import smiles

    mol = smiles('CCO')          # ethanol
    mol = smiles('c1ccccc1')     # benzene (aromatic)
    mol = smiles('[Cu+2]')       # copper ion
    mol = smiles('C/C=C/C')     # trans-2-butene (with stereo)

    # Reaction SMILES
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')


InChI
~~~~~

.. code-block:: python

    from chython import inchi

    mol = inchi('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')


IUPAC Name
~~~~~~~~~~

Requires OPSIN JAR. Set path via ``OPSIN_PATH`` env variable or ``chython.class_paths[1]``.

.. code-block:: python

    from chython import iupac

    mol = iupac('acetic acid')
    mol = iupac('2-acetoxybenzoic acid')


MDL MOL Block
~~~~~~~~~~~~~

.. code-block:: python

    from chython import mdl_mol

    mol_block = """
      Mrv2211 03232310102D

      3  2  0  0  0  0            999 V2000
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        3.0800    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
      2  3  1  0  0  0  0
    M  END
    """
    mol = mdl_mol(mol_block)


XYZ Coordinates
~~~~~~~~~~~~~~~

.. code-block:: python

    from chython import xyz

    mol = xyz((('O', 0., 0., 0.), ('H', 1., 0., 0.), ('H', 0., 1., 0.)))
    mol.clean2d()


SMARTS
~~~~~~

Parses SMARTS into ``QueryContainer`` for substructure matching. See :doc:`substructure` for full SMARTS syntax.

.. code-block:: python

    from chython import smarts

    q = smarts('[C;r5,r6;a]-;!@[C;h0,h1]')
    print(q)  # canonical atom order


File Readers
------------

SDF / RDF
~~~~~~~~~

.. code-block:: python

    from chython import SDFRead, RDFRead

    # Iterate molecules from SDF
    with SDFRead('molecules.sdf') as reader:
        for mol in reader:
            print(mol.name, mol.meta)

    # Iterate reactions from RDF
    with RDFRead('reactions.rdf') as reader:
        for rxn in reader:
            print(rxn)

    # Read all at once
    with SDFRead('molecules.sdf') as reader:
        mols = reader.read()

    # Read a limited batch
    with SDFRead('molecules.sdf') as reader:
        first_100 = reader.read(amount=100)

    # Generator: get first record, then read rest
    with RDFRead('reactions.rdf') as f:
        first = next(f)
        rest = f.read()

    # Indexed random access (Unix only)
    with SDFRead('molecules.sdf', indexable=True) as reader:
        reader.reset_index()
        mol = reader[42]
        total = len(reader)

Pathlib supported:

.. code-block:: python

    from pathlib import Path

    with RDFRead(Path('reactions.rdf')) as r:
        rxn = next(r)

Opened file objects supported (text mode for all formats except MRV):

.. code-block:: python

    with open('reactions.rdf') as f, RDFRead(f) as r:
        rxn = next(r)


MRV
~~~

MRV files require **binary mode**:

.. code-block:: python

    from chython import MRVRead

    with MRVRead(open('structures.mrv', 'rb')) as reader:
        for mol in reader:
            print(mol)


Reading from Archives and Network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Readers accept any file-like object, enabling transparent reading from compressed sources:

.. code-block:: python

    # gzip
    from gzip import open as gzip_open
    with gzip_open('data.rdf.gz', 'rt') as f, RDFRead(f) as r:
        rxn = next(r)

    # zip
    from zipfile import ZipFile
    from io import TextIOWrapper
    with ZipFile('data.zip') as z, z.open('data.rdf') as c:
        with TextIOWrapper(c) as f, RDFRead(f) as r:
            rxn = next(r)

    # tar.gz
    from tarfile import open as tar_open
    with tar_open('data.tar.gz') as t:
        c = t.extractfile('data.rdf')
        with TextIOWrapper(c) as f, RDFRead(f) as r:
            rxn = next(r)

    # URL via requests
    from requests import get
    from io import StringIO
    with StringIO(get('https://example.com/data.rdf').text) as f, RDFRead(f) as r:
        rxn = next(r)


Other Readers
~~~~~~~~~~~~~

All readers share the same API (iteration, ``.read()``, context manager).

- **SDFRead** - MOL/SDF (V2000, V3000)
- **RDFRead** - RXN/RDF
- **MRVRead** - ChemAxon MRV (requires binary mode)
- **PDBRead** - PDB format (explicit hydrogens only)

For SMILES, InChI and XYZ, use the string parsers (``smiles()``, ``inchi()``, ``xyz()``) directly.
To process files with one record per line, iterate lines manually:

.. code-block:: python

    from chython import smiles

    with open('molecules.smi') as f:
        for line in f:
            mol = smiles(line)
            print(mol)


Reader Options
~~~~~~~~~~~~~~

MDL readers (SDFRead, RDFRead) accept these options:

.. code-block:: python

    with SDFRead('molecules.sdf',
                 ignore=True,              # try to fix/skip errors (default)
                 remap=False,              # renumber atoms from 1
                 ignore_stereo=False,      # discard stereochemistry
                 ignore_bad_isotopes=False, # reset invalid isotopes
                 calc_cis_trans=False,      # recalculate cis/trans from 2D
                 ) as reader:
        for mol in reader:
            pass


File Writers
------------

SDF / RDF
~~~~~~~~~

.. code-block:: python

    from chython import SDFWrite, RDFWrite, ESDFWrite, ERDFWrite

    # Write molecules to SDF (V2000)
    with SDFWrite('output.sdf') as writer:
        writer.write(mol)

    # V3000 extended format
    with ESDFWrite('output_v3000.sdf') as writer:
        writer.write(mol)

    # Write reactions to RDF (V2000)
    with RDFWrite('output.rdf') as writer:
        writer.write(rxn)

    # V3000 reactions
    with ERDFWrite('output_v3000.rdf') as writer:
        writer.write(rxn)

    # Append mode
    with SDFWrite('output.sdf', append=True) as writer:
        writer.write(mol)

    # Write with 3D coordinates (conformer index)
    with SDFWrite('output_3d.sdf') as writer:
        writer.write(mol, write3d=0)

    # Ongoing writing without context manager
    f = RDFWrite('output.rdf')
    for rxn in data:
        f.write(rxn)
    f.close()


MRV
~~~

.. code-block:: python

    from chython import MRVWrite

    with MRVWrite('output.mrv') as writer:
        writer.write(mol)


SMILES Strings
~~~~~~~~~~~~~~

.. code-block:: python

    from chython import smiles

    mol = smiles('CCO')

    # Canonical SMILES
    s = str(mol)           # or format(mol)

    # Format specifiers
    format(mol, 'm')   # include atom mapping numbers
    format(mol, 'h')   # show implicit hydrogens
    format(mol, 'r')   # random SMILES (non-canonical)
    format(mol, 'a')   # asymmetric closures
    format(mol, 'A')   # aromatic bonds (: notation) instead of lowercase atoms
    format(mol, '!s')  # without stereo
    format(mol, '!x')  # without CXSMILES extensions
    format(mol, '!z')  # without charges
    format(mol, '!b')  # without bond tokens
    format(mol, 'mh')  # combine multiple: mapping + hydrogens
    format(mol, 'h!b') # implicit H, no bond tokens

    # Works with f-strings, %-formatting, .format()
    print(f'{mol:A}')
    print('smiles: %s' % mol)


Serialization
-------------

Pickle
~~~~~~

Full pickle support for all containers. Faster than file formats for temporary storage:

.. code-block:: python

    from pickle import loads, dumps

    data = dumps(mol)
    mol = loads(data)

    # Works for reactions too
    data = dumps(rxn)
    rxn = loads(data)


Chython Binary Pack
~~~~~~~~~~~~~~~~~~~

Compact binary format. Stores 2D coordinates, stereo, charges, isotopes, radicals, atom numbers.
Size ~1.5-2x larger than SMILES. Parsing faster than pickle.

.. code-block:: python

    from chython import MoleculeContainer, ReactionContainer, unpack

    # Pack to bytes (zlib compressed by default)
    data = mol.pack()
    data = bytes(mol)          # same as pack()

    # Unpack (auto-detects molecule or reaction)
    restored = unpack(data)

    # Or unpack with explicit type
    restored = MoleculeContainer.unpack(data)

    # Uncompressed
    data = mol.pack(compressed=False)
    restored = unpack(data, compressed=False)

    # Peek at atom count without unpacking
    count = MoleculeContainer.pack_len(data, compressed=False)

    # Reactions
    rxn_data = rxn.pack()
    rxn_restored = ReactionContainer.unpack(rxn_data)


Metadata
--------

SDF/RDF files store metadata in molecule and reaction objects:

.. code-block:: python

    rxn = next(RDFRead('reactions.rdf'))
    rxn.meta           # dict of DTYPE/DATUM fields
    rxn.name           # reaction title from RDF

    mol = rxn.reactants[0]
    mol.name           # molecule title from MOL block
    mol.meta           # molecule metadata dict

    # Set metadata for writing
    mol.name = 'Ethanol'
    mol.meta['boiling_point'] = '78.37'