Input / Output
==============

Reading and writing molecules and reactions in all supported formats.


One Name, Both Directions
-------------------------

A format whose two directions are one operation has **one callable**, and the direction is the
argument's type: a container in exports, a string or a buffer in imports. The direction-stating names
stay beside it, and a loop calls those — they are where the keywords that serve one direction only
live, and where a reader that reports damage instead of raising lives.

=================  ==============  =====================================  ==================================
Format             Both ways       Import only                            Export only
=================  ==============  =====================================  ==================================
SMILES             ``smiles``      ``read_smiles``,                       ``write_smiles``,
                                   ``read_reaction_smiles``               ``write_reaction_smiles``,
                                                                          ``detached_smiles``
InChI              ``inchi``       ``inchi_to_molecule``                  ``molecule_to_inchi``,
                                                                          ``inchikey``
IUPAC name         ``iupac``       --                                     --
MDL molfile        ``mol``         --                                     --
MDL RXN            ``rxn``         --                                     --
MRV                ``mrv``         ``read_mrv``                           ``write_mrv``
CML                ``cml``         ``read_cml``                           ``write_cml``
pach               ``pach``        ``unpach``/``unpack``, ``pach_load``,  ``pach_dump``,
                                   ``reaction_pach_load``                 ``reaction_pach_dump``
SMARTS / SMIRKS    --              ``smarts``, ``read_smirks``            -- (a query has no string form)
MOL2               --              ``mol2``, ``read_mol2``, ``mol2_mol``  --
XYZ                --              ``xyz``                                --
PDB / mmCIF        --              ``pdb``, ``mmcif``, ``read_pdb``,      --
                                   ``read_mmcif``
SDF / RDF files    --              ``SDFRead``, ``RDFRead``               ``SDFWrite``, ``ESDFWrite``,
                                                                          ``RDFWrite``, ``ERDFWrite``
arena              --              ``MoleculeContainer.from_bytes``       ``mol.to_bytes()``
=================  ==============  =====================================  ==================================

``inchikey`` has no import direction at all and never will: nothing reads a hash back. XYZ, PDB and
mmCIF import a record of atoms and coordinates rather than a ``MoleculeContainer``. PDB and mmCIF
state connectivity and no orders, and ``chython.chemistry.saturate()`` is the explicit call that
perceives the orders once ``chython.formats.pdb.build_molecule()`` has built the bonds. XYZ states no
bond either, so a frame becomes a molecule through one call more: ``build_molecule()`` places the atoms,
``perceive_bonds()`` reads the connectivity out of the geometry, ``saturate()`` the orders. Every one of
them is the caller's, and no reader on this page perceives anything.


String Parsers
--------------

SMILES
~~~~~~

``smiles()`` goes both ways, the direction decided by the argument's type:

.. testcode::

    from chython import smiles

    mol = smiles('CCO')          # ethanol
    mol = smiles('c1ccccc1')     # benzene (aromatic)
    mol = smiles('[Cu+2]')       # copper ion
    mol = smiles('C/C=C/C')     # trans-2-butene (with stereo)

    # Reaction SMILES
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')

    # and back, for a molecule and for a reaction alike.  `m` writes the map numbers, which the
    # default spec does not: an atom-to-atom mapping is a claim about one record, not part of the
    # structure's identity.
    print(smiles(mol), smiles(rxn), smiles(rxn, spec='m'), sep='\n')

.. testoutput::

    C(/C)=C\C
    CO>>CN
    [CH3:1][OH:2]>>[CH3:1][NH2:3]

The **options live on the function**, as they do for InChI below: ``spec=`` is the ``format()`` spec, so
``smiles(mol, spec='a')`` and ``format(mol, 'a')`` are one string and ``mol.smiles`` never grows a
parameter.  ``read_smiles`` is the import-only name, and it is the one a file loop calls: it takes the
``log`` list that a door deciding a direction cannot.

.. testcode::

    from chython import read_smiles

    log = []
    for line in ('CCO', 'C-1CCCCC=1'):
        record = read_smiles(line, log)
    print(len(log), 'line(s) about the whole run')

.. testoutput::

    1 line(s) about the whole run

Two export methods write a string meant to be **concatenated** rather than read on its own, which is why
neither is canonical, neither is cached and neither is a ``format()`` spec: the result is a function of
the atoms the caller named, not of the molecule alone.

``sticky_smiles()`` writes a traversal that starts at one named atom and ends at another. ``remove_*``
drops that end's atom token and ``keep_bond_*`` keeps its bond token, so a piece can be left open at
either end:

.. testcode::

    ethanolamine = smiles('OCCN')
    print(ethanolamine.sticky_smiles(left=1, right=4))
    print(ethanolamine.sticky_smiles(left=1, remove_left=True, keep_bond_left=True))

.. testoutput::

    OCCN
    -CCN

``detached_smiles()`` cuts bonds instead, leaving each cut as an **open ring bond**: the halves are
separate ``.`` components of the text and the shared closure id is what makes them one molecule again.
``cuts`` is ``{attachment id: (keep, drop)}``, the pair ordered because ``C-C`` does not say which half
the caller wants, and an attachment id is 10..99 -- a two-digit ``%nn`` closure no ordinary written ring
uses:

.. testcode::

    from chython import DetachedSmiles

    amine = smiles('OCCN').detached_smiles({10: (2, 1)})       # keep the chain, drop the oxygen
    aryl = smiles('c1ccccc1Br').detached_smiles({10: (6, 7)})  # keep the ring, drop the bromine
    print(amine.text, aryl.text)
    print(smiles(DetachedSmiles.join(amine, aryl).text))

.. testoutput::

    C%10CN c1c%10cccc1
    c1c(cccc1)CCN

``join`` is a static call over the pieces: it concatenates the bodies as components and reports what is
left open, so an enumeration builds strings and parses once at the end. ``reserve`` withholds another
fragment's attachment ids, which is how several pieces are written for one join without colliding.
:doc:`reactions` documents the corpus-driven form of the same operation, where a row of ``roles.tsv``
decides where the cut goes.


InChI
~~~~~

Requires the InChI library (shipped with the wheel; available when
``chython.inchi_library_loaded()`` returns ``True``).  ``inchi()`` goes both ways, the direction
decided by the argument's type:

.. testcode::

    from chython import inchi, inchikey

    mol = inchi('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')   # ethanol from InChI
    text = inchi(mol)                                  # and back
    key = inchikey(mol)                                # one-way: nothing reads a key back

    print(mol.inchi)                                   # the same two, as properties
    print(mol.inchikey)

.. testoutput::

    InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3
    LFQSCWFLJHTTHZ-UHFFFAOYSA-N

``inchi_to_molecule``, ``molecule_to_inchi`` and ``molecule_to_inchikey`` are the same functions under
their direction-stating names. The **options live on the function** and never on the property — the same
split as ``mol.smiles`` and ``format(mol, spec)``, so a property never grows a parameter:

.. testcode::

    from chython import molecule_to_inchi

    butene = smiles('C/C=C/C')
    print(butene.inchi)
    print(molecule_to_inchi(butene, options='-SNon'))   # no stereo

.. testoutput::

    InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+
    InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3

Both properties raise ``ImportError`` when libinchi is not loaded, and neither is cached: a string
cached on a mutable container outlives the edit that makes it wrong.


IUPAC Name
~~~~~~~~~~

Requires OPSIN JAR. Set path via ``OPSIN_PATH`` env variable or ``chython.class_paths[1]``.

.. testcode::
   :skipif: not any(__import__('os').path.isfile(str(p)) for p in __import__('chython').class_paths if 'opsin' in str(p).lower())

    from chython import iupac

    mol = iupac('acetic acid')
    mol = iupac('2-acetoxybenzoic acid')


MDL MOL Block
~~~~~~~~~~~~~

``mol()`` goes both ways, the direction decided by the argument's type; ``rxn()`` is the same door for
an RXN block.  A record copied out of an SD file brings its ``$$$$`` along, which is not an error.

.. testcode::

    from chython import mol

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
    ethanol = mol(mol_block)
    block = mol(ethanol)                 # and back, V2000 unless the molecule needs V3000


XYZ Coordinates
~~~~~~~~~~~~~~~

``xyz()`` takes an XYZ-format string and returns a list of ``XYZFrame`` objects — not a
``MoleculeContainer``, because XYZ states no bonds and a reader that invented them would be
repairing on the read path:

.. testcode::

    from chython import xyz

    xyz_text = """3
    water
    O  0.000  0.000  0.000
    H  0.757  0.586  0.000
    H -0.757  0.586  0.000
    """
    frames = xyz(xyz_text)
    frame = frames[0]    # XYZFrame; atoms and coordinates, no bonds

Every frame in the string comes back, so a trajectory reads whole rather than first-frame-only.

``build_molecule()`` is the call that turns frames into a molecule — the caller's, never the reader's.
It lives at ``chython.formats.xyz`` and not on the façade, where ``build_molecule`` is the PDB record's
builder: the same question asked of a different record. The first frame states the atoms, every frame
becomes one conformer, and nothing arrives bonded:

.. testcode::

    from chython.formats.xyz import build_molecule

    frame_mol = build_molecule(frames)   # the water frames read above
    print(len(frame_mol), len(list(frame_mol.bonds())), frame_mol.has_3d)

.. testoutput::

    3 0 True

Every atom states **zero** implicit hydrogens, which is the format's own statement: an XYZ record lists
every atom, so a hydrogen not written is a hydrogen not there. Two explicit passes finish the job —
``perceive_bonds()`` reads the connectivity out of one stored model, and ``saturate()`` raises the orders
the hydrogen counts force. Both are documented in :doc:`standardize`.

.. testcode::

    from chython import perceive_bonds

    print(perceive_bonds(frame_mol), frame_mol)

.. testoutput::

    True [H]O[H]

``xyz_conformers()`` is the other direction of that rule: given a molecule the caller already has, it
stores the frames as its conformers and returns how many landed. Atoms match **positionally** — frame
atom *i* is ``molecule.atom_numbers[i]`` — and a frame whose length or element sequence disagrees with
the molecule is logged and skipped, the rest still landing. Frames append to whatever models the
molecule already carries, and each frame's ordinal in the sequence becomes its ``ext_index``.

.. testcode::

    from chython import smiles, xyz_conformers

    water = smiles('O([H])[H]')          # O, H, H -- the frame's order
    stored = xyz_conformers(water, frames)
    print(stored, len(water.conformers), water.conformer(0).ext_index)
    print(water.conformer(0).xyz_of(water.atom_numbers[1]))

.. testoutput::

    1 1 0
    (0.757, 0.586, 0.0)


SMARTS
~~~~~~

Parses SMARTS into ``QueryContainer`` for substructure matching. See :doc:`substructure` for full SMARTS syntax.

.. testcode::

    from chython import smarts

    q = smarts('[C;r5,r6;a]-;!@[C;h0,h1]')
    print(len(q), 'atoms')

.. testoutput::

    2 atoms

A ``QueryContainer`` is a sealed query, not a string: it has no SMARTS round trip, so print what you
want to know about it rather than the container itself.


File Readers
------------

.. testsetup::

    # Write the sample files that the file-reading examples below open.
    from chython import smiles as _s, SDFWrite as _SW, RDFWrite as _RW, ReactionContainer as _RC
    import gzip as _gz, zipfile as _zf, tarfile as _tf, io as _io

    _benzene = _s('c1ccccc1')
    _benzene.set_title(b'benzene')
    _aspirin = _s('CC(=O)Oc1ccccc1C(=O)O')
    _aspirin.set_title(b'aspirin')

    with _SW('molecules.sdf') as _f:
        _f.write(_benzene)
        _f.write(_aspirin)

    _rxn_eg = _RC(reactants=(_s('CCO'),), products=(_s('CC=O'),))

    with _RW('reactions.rdf') as _f:
        _f.write(_rxn_eg)
        _f.write(_rxn_eg)

    _rdf_buf = _io.StringIO()
    with _RW(_rdf_buf) as _f:
        _f.write(_rxn_eg)
    _rdf_text = _rdf_buf.getvalue()

    with _gz.open('data.rdf.gz', 'wt') as _f:
        _f.write(_rdf_text)

    with _zf.ZipFile('data.zip', 'w') as _f:
        _f.writestr('data.rdf', _rdf_text)

    _rdf_bytes = _rdf_text.encode()
    with _tf.open('data.tar.gz', 'w:gz') as _f:
        _info = _tf.TarInfo(name='data.rdf')
        _info.size = len(_rdf_bytes)
        _f.addfile(_info, _io.BytesIO(_rdf_bytes))

    with open('molecules.smi', 'w') as _f:
        _f.write('c1ccccc1 benzene\n')
        _f.write('CC(=O)Oc1ccccc1C(=O)O aspirin\n')


SDF / RDF
~~~~~~~~~

.. testcode::

    from chython import SDFRead, RDFRead

    # The name line and the data fields are ON THE MOLECULE
    with SDFRead('molecules.sdf') as reader:
        for mol in reader:
            print(mol.title, dict(mol.meta))

.. testoutput::

    benzene {}
    aspirin {}

.. note::

   ``mol.title``, ``mol.meta`` and ``mol.log`` travel with the molecule, so ``list(SDFRead(f))``
   keeps them and ``SDFWrite`` writes them back.  ``reader.title`` and ``reader.meta`` read the
   record just yielded -- the same two values off the molecule the reader last returned -- so they
   are a convenience inside the loop and say nothing once it moves on.  What no container holds is
   the framing: ``reader.version`` and the rest stay on the reader.

.. testcode::

    # Iterate reactions from RDF
    with RDFRead('reactions.rdf') as reader:
        for rxn in reader:
            print(rxn)

.. testoutput::

    C(C)O>>C(C)=O
    C(C)O>>C(C)=O

.. testcode::

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

Pathlib supported:

.. testcode::

    from pathlib import Path

    with RDFRead(Path('reactions.rdf')) as r:
        rxn = next(r)

Opened file objects supported (text mode):

.. testcode::

    with open('reactions.rdf') as f, RDFRead(f) as r:
        rxn = next(r)


MRV
~~~

``mrv()`` goes both ways over XML strings, the direction decided by the argument's type. A read always
answers a list, however many molecules the document holds:

.. testcode::

    from chython import smiles, mrv

    mol = smiles('c1ccccc1')

    mrv_string = mrv(mol)          # a molecule, or a list of them, out to MRV XML
    mols = mrv(mrv_string)         # XML text back to a list of MoleculeContainer

    # Read an MRV file: these take text, so hand them the file's contents
    with open('structures.mrv', 'w') as fh:
        fh.write(mrv_string)
    with open('structures.mrv') as fh:
        mols = mrv(fh.read())

``read_mrv`` and ``write_mrv`` are the same functions under their direction-stating names, and
``read_xml`` sniffs the dialect instead of being told it.


Reading from Archives
~~~~~~~~~~~~~~~~~~~~~

Readers accept any file-like object, enabling transparent reading from compressed sources:

.. testcode::

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


Other Readers
~~~~~~~~~~~~~

All readers share the same API (iteration, ``.read()``, context manager).

- **SDFRead** / **SDFWrite** — MOL/SDF (V2000, V3000)
- **RDFRead** / **RDFWrite** — RXN/RDF
- **mrv** (**read_mrv** / **write_mrv**) — ChemAxon MRV (XML strings)
- **cml** (**read_cml** / **write_cml**) — Chemical Markup Language (XML strings)
- **read_pdb** / **read_mmcif** — legacy PDB and PDBx/mmCIF; return ``PDBRecord`` rather than a molecule
- **mol2** / **read_mol2** / **mol2_mol** — Tripos MOL2

``chython.formats.pdb.build_molecule()`` takes one ``PDBRecord`` or a sequence of them, and a sequence
collapses to **one** molecule carrying a conformer per record — the shape of an NMR ensemble, where
each ``MODEL`` is the same molecule again. Atoms match by chain, residue sequence, insertion code,
residue name, atom name and alternate id, never by ``serial``, which no spec pins across models. A
model whose atom set differs from the first is logged and skipped whole; the first record's geometry is
model 0, and each ``MODEL`` serial is stored as that model's ``ext_index``.

For SMILES, InChI, XYZ and MOL2 use the string callables (``smiles()``, ``inchi()``, ``xyz()``,
``mol2()``) directly. ``read_mol2`` and ``read_pdb`` open *files*; the callables read *text*, so a
document is never mistaken for a filename. To process files with one record per line, iterate lines
manually:

.. note::

   A reader that returns a molecule puts that record's lines on ``mol.log`` as well, so the damage
   report survives a caller who passed no ``log=``. One that returns a record object instead —
   ``PDBRecord``, ``XYZFrame`` — keeps them on the record, there being no container yet;
   ``build_molecule()`` folds the record's log onto the molecule it builds. A ``log=`` list is always a
   second copy of the same lines, per record and never pooled: molecule 3's lines are on molecule 3.

.. testcode::

    from chython import smiles

    with open('molecules.smi') as f:
        for line in f:
            mol = smiles(line.split()[0])
            print(mol)

.. testoutput::

    c1ccccc1
    O=C(Oc1c(C(O)=O)cccc1)C


Reader Options
~~~~~~~~~~~~~~

MDL readers (SDFRead, RDFRead) accept this option:

.. testcode::

    with SDFRead('molecules.sdf',
                 ignore_stereo=False,      # skip the stereo step (constitution always read)
                 ) as reader:
        for mol in reader:
            pass


File Writers
------------

SDF / RDF
~~~~~~~~~

.. testcode::

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

    # Ongoing writing without context manager
    f = RDFWrite('output.rdf')
    for r in [rxn]:
        f.write(r)
    f.close()


MRV
~~~

Handed a container, ``mrv()`` returns an XML string:

.. testcode::

    from chython import mrv

    mrv_string = mrv(mol)
    with open('output.mrv', 'w') as fh:
        fh.write(mrv_string)


CML
~~~

Chemical Markup Language works the same way through ``cml()``, again with the direction taken from the
argument:

.. testcode::

    from chython import cml, read_xml

    cml_string = cml(mol)          # a container out to CML XML
    mols = cml(cml_string)         # XML text back to a list of MoleculeContainer

    # Dialect-sniffing entry point (CML or MRV, chosen from the document's root element)
    mols = read_xml(cml_string)

``read_cml`` and ``write_cml`` are the same functions under their direction-stating names.


SMILES Strings
~~~~~~~~~~~~~~

.. testcode::

    from chython import smiles

    mol = smiles('CCO')
    print(str(mol))            # canonical SMILES; same as format(mol)

.. testoutput::

    C(C)O

``format(mol, spec)`` takes these specifiers. A lowercase letter adds something, and ``!`` withdraws
something that is on by default; they combine in one string (``'mh'``, ``'h!b'``).

=========  ===============================================================
specifier  effect
=========  ===============================================================
``m``      include atom map numbers
``h``      write implicit hydrogen counts explicitly
``a``      asymmetric ring closures
``A``      aromaticity on the bonds (``:``) instead of lowercase atoms
``r``      a fresh random atom order instead of the canonical one
``!s``     without stereo
``!x``     without CXSMILES extensions
``!z``     without charges
``!b``     without bond tokens
=========  ===============================================================

Only the specifiers a given molecule can actually exercise change its string, so ethanol above is
unchanged by all of them but ``h``:

.. testcode::

    print(format(mol, 'h'))
    print(format(smiles('[CH3:1][OH:2]'), 'm'))
    print(format(smiles('c1ccccc1'), 'A'))

.. testoutput::

    [CH2]([CH3])[OH]
    [CH3:1][OH:2]
    [CH]:1:[CH]:[CH]:[CH]:[CH]:[CH]:1

``A`` brackets every aromatic atom and states its hydrogen count, as above. That is not verbosity:
``C:C`` is outside OpenSMILES — an aromatic bond between aliphatic atoms — so no rule says what
hydrogen count it implies and readers differ. Stating the count makes every reader agree.

``r`` writes the same molecule a different way each call, which is what a training-set augmenter wants:

.. testcode::

    from random import seed

    seed(0)
    ethanol = smiles('CCO')
    print(len({format(ethanol, 'r') for _ in range(20)}) > 1)
    print(all(smiles(format(ethanol, 'r')) == ethanol for _ in range(20)))

.. testoutput::

    True
    True

Two consequences. The draws come from the ``random`` module, so ``random.seed()`` makes a run
repeatable even though its strings are not predictable. And an ``r`` string is **not** canonical, so it
must never reach a hash or an equality test — ``==`` compares ``canonical_bytes`` and never a string.
``r`` and the stored-order key ``i`` both name where the atom order comes from and raise together.

Any spelling that reaches ``__format__`` works — f-strings, ``%``, ``str.format``:

.. testcode::

    print(f'{mol:h}')
    print('smiles: %s' % mol)

.. testoutput::

    [CH2]([CH3])[OH]
    smiles: C(C)O


Data Labels (S-groups)
----------------------

A CTfile ``DAT`` S-group is text attached to atoms or bonds -- a stereo descriptor beside a centre, a
note beside a fragment.  ``mol.add_data_sgroup()`` attaches one in a single call and both emitters write
it, so the label survives a round trip through either CTAB version:

.. testcode::

    from chython import mol as mol_facade

    block = """
      chython

      4  3  0  0  0  0            999 V2000
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        3.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
      2  3  1  0  0  0  0
      3  4  1  0  0  0  0
    M  END
    """
    butane = mol_facade(block)          # 2D coordinates, so an anchor can be computed
    n, m = butane.atom_numbers[1], butane.atom_numbers[2]

    record = butane.add_data_sgroup('StereoLabel', '(R)', atoms=[n])
    print(record.field_data, record.disp[:2])

    # both versions carry it; `data_sgroups` reads them back by FIELDNAME
    back, = mol_facade(mol_facade(butane, version=3000)).data_sgroups('StereoLabel')
    print(back.field_data)

.. testoutput::

    (R) (1.0, 0.0)
    (R)

The anchor is the ``FIELDDISP`` display position.  By default it is the mean of the referenced atoms'
coordinates, which puts a bond label between its endpoints when the bond's atoms are listed too:

.. testcode::

    label = butane.add_data_sgroup('StereoLabel', '(Z)', atoms=[n, m], bonds=[(n, m)])
    print(label.disp[:2])

.. testoutput::

    (1.5, 0.5)

``disp=(x, y)`` states the position instead, and ``disp=False`` writes no anchor at all.  A molecule
with no coordinates has none to give, so the record is written without one and the reason is logged
rather than raised.  ``mol.data_sgroups()`` with no name answers every ``DAT`` record, and each call
appends -- unlike ``set_sgroups``, which replaces the whole set:

.. testcode::

    butane.add_data_sgroup('NOTE', ['first', 'second'], atoms=[n])    # a multi-value datum
    print(sorted(r.name for r in butane.data_sgroups()))

.. testoutput::

    ['NOTE', 'StereoLabel', 'StereoLabel']

``add_data_sgroup(molecule, ...)`` and ``data_sgroups(molecule)`` are the same two as functions, for a
caller holding the molecule at arm's length. ``mol.sgroups`` beside them is the arena's own dicts for
every S-group kind there is; these two are the parsed ``DAT`` records.


Atom Labels and Markers
-----------------------

A file may name an atom with something that is not an element: a drawn abbreviation (``OMe``), a
polymer bead (``Pol``), a resin support, a registry identifier. Every one of them reads as **the
marker** — element 0, matching nothing — carrying that text as its **alias**, so the record is read
and its own words are kept. What a label abbreviates is not in the structure until an explicit pass
puts it there.

.. testcode::

    from chython import smiles

    labelled = smiles('[Pol]CC')
    marker = next(a for a in labelled.atoms() if a.is_r)
    print(marker.element, labelled.aliases[marker.n])

.. testoutput::

    0 b'Pol'

``*`` is that same marker and ``[R7]`` is one carrying an index. The CXSMILES ``$...$`` field is read
and written in the same index space as ``^1:``, which is how a label survives a round trip through the
string; ``_R<n>`` there is ChemAxon's spelling of the index and ``_AP<n>`` an attachment point, whose
ordinal nothing in the arena holds.

.. testcode::

    labelled = smiles('CC* |$;;_R1$|')        # the marker is in the body, its index in the tail
    print(next(a for a in labelled.atoms() if a.is_r).r_index)

    labelled = smiles('CCC |$;;OMe$|')        # the body writes a carbon, so a carbon is what is stored
    tagged = next(a for a in labelled.atoms() if a.n in labelled.aliases)
    print(tagged.element, tagged.is_r, labelled.aliases[tagged.n])
    print(smiles(labelled))

.. testoutput::

    1
    6 False b'OMe'
    C(C)C |$;OMe;$|

**The body decides the element and the tail only adds a word to it.** ``[OMe]`` names no element, so
that atom is the marker; ``CCC |$;;OMe$|`` names a carbon, so the atom is a carbon carrying the label,
and ``CCO |$;;OMe$|`` an oxygen carrying it. A ``_R<n>`` on an atom the body spells as an element is
dropped with a log line for the same reason — an R index belongs to a marker. What the label
abbreviates is put in by ``chython.chemistry.expand_abbreviations``, which reads the alias in every one
of those cases and gives ``C(C)OC``.

A character the field cannot hold — ``;`` ends an entry, ``$`` the field, ``|`` the block, ``&`` a
reference — travels as a numeric character reference, in both directions:

.. testcode::

    print(smiles(smiles('CC |$a&#59;b;$|')))

.. testoutput::

    CC |$;a&#59;b$|


Serialization
-------------

Pickle
~~~~~~

Full pickle support for all containers. Faster than file formats for temporary storage:

.. testcode::

    from pickle import loads, dumps

    data = dumps(mol)
    mol = loads(data)

    # Works for reactions too
    data = dumps(rxn)
    rxn = loads(data)


Chython Binary Pack
~~~~~~~~~~~~~~~~~~~

Compact binary format, and the one to store when something other than this release will read the
bytes back.  :doc:`pach` is its wire specification, every version byte of it; this is the API.

``pach()`` is the short door and goes both ways; ``unpach()`` is its import half under its own name,
with ``unpack`` as chython 2's spelling of the same function.  Nothing has to be declared about a
buffer: byte 0 says which era wrote it -- molecule pach, reaction pach, or the arena -- and whether
anybody compressed it.

.. testcode::

    from chython import pach, unpach

    data = pach(mol)                      # zlib compressed, current layout
    restored = unpach(data)               # molecule or reaction, raw or compressed
    restored = unpach(pach(rxn))          # the reaction shape through the same door

    raw = pach(mol, compressed=False)     # the record itself
    old = pach(mol, version=2, drop='*')  # the legacy layout, losing what it cannot carry

``unpach`` chooses between two error policies, and ``log=`` is the choice. Without it, it is an answer
boundary and raises ``ValueError`` for a damaged record — including damage it recovered from, since a
caller who cannot have a structure is owed the whole story. With a list, it is the loop-safe door:

.. testcode::

    log = []
    store = [pach(mol), b'not a record at all', pach(rxn)]
    good = [x for x in (unpach(record, log=log) for record in store) if x is not None]
    print(len(good), 'of', len(store), 'read;', len(log), 'complaint(s)')

.. testoutput::

    2 of 3 read; 1 complaint(s)

``pach_load`` and ``reaction_pach_load`` are the same reporting behaviour under direction-stating
names, returning ``(structure, problems)`` per shape; ``pack``/``unpack`` on the containers are the
methods, and ``pach``/``unpach`` on them chython 2's names for those.

.. testcode::

    from chython import MoleculeContainer, ReactionContainer, pach_load, pach_record_length

    data = mol.pack()
    restored, problems = pach_load(data)
    restored = MoleculeContainer.unpack(data)
    rxn_restored = ReactionContainer.unpack(rxn.pack())

    # Byte length of the first record in a buffer, for walking a stream of them
    n = pach_record_length(data)

chython 2's ``check=``, ``order=`` and ``skip_labels_calculation=`` are a ``TypeError`` rather than
accepted and ignored. ``drop=`` is what waives a refusal here, and it names the field it waives, so a
silently swallowed ``check=False`` would promise a waiver and then let the encoder raise anyway.


Metadata
--------

``MoleculeContainer`` and ``ReactionContainer`` both carry ``meta``, a plain ``dict`` of the
record's data fields, and ``title``, a ``str``.  A reader fills them, so a writer needs no keyword
to put them back:

.. testcode::

    rxn = next(RDFRead('reactions.rdf'))
    rxn.meta           # dict of DTYPE/DATUM fields
    rxn.title          # reaction title from RDF (str)

    mol = rxn.reactants[0]
    mol.meta           # dict of the SDF data fields the record carried
    mol.title          # molecule title from MOL block (str)

    # Set the title on a molecule or reaction
    mol.set_title('Ethanol')

    # Metadata is a dict on the container; write to it directly
    rxn.meta['boiling_point'] = '78.37'
    mol.meta['boiling_point'] = '78.37'

    # The writer reads it off the molecule; `meta=` overrides, and `meta={}` writes none
    with SDFWrite('output.sdf') as writer:
        writer.write(mol)

The XML formats carry the same dict as a ``<propertyList>`` of ``<property>``/``<scalar>`` inside
``<molecule>`` -- where ``molconvert`` puts an SD field converting an SDF -- so a field survives a
round trip through either dialect, newlines and markup characters included:

.. testcode::

    from chython import cml, mrv

    mol.meta['note'] = 'R at C2'
    back, = mrv(mrv(mol))      # a read answers a list, however many molecules the document holds
    assert back.meta['note'] == 'R at C2'
    back, = cml(cml(mol))
    assert back.meta['note'] == 'R at C2'
