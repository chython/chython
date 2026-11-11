HELM and BILN
=============

``helm()`` and ``biln()`` read a biopolymer notation into a molecule and write a molecule as its
**canonical** notation: one molecule, one string, whatever order its atoms were drawn in and however the
notation it came from was spelt.

.. testcode::

    from chython import helm, biln, smiles

    oxytocin = helm('PEPTIDE1{C.Y.I.Q.N.C.P.L.G.[am]}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$')
    print(helm(oxytocin))
    print(biln(oxytocin))

    shuffled = smiles(format(oxytocin, 'r'))       # the same molecule, atoms in another order
    assert helm(shuffled) == helm(oxytocin)

.. testoutput::

    PEPTIDE1{C.Y.I.Q.N.C.P.L.G.[am]}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$
    C(1,3)-Y-I-Q-N-C(1,3)-P-L-G-am


The Monomer Library
-------------------

Monomers come from a library in the Pistoia HELM JSON format, each record a SMILES whose leaving atoms
carry their R number as a map number (``C[C@H](N[H:1])C([OH:2])=O``). The shipped default is the Pistoia
HELMCoreLibrary (MIT). ``chython.monomer_library`` stacks more libraries on top of it, a later file's
symbol overriding an earlier one's, and ``()`` restores the default; the set is compiled once and cached:

.. code-block:: text

    chython.monomer_library = 'my_monomers.json'          # one file
    chython.monomer_library = ['public.json', 'mine.json']
    chython.monomer_library = ()                          # the shipped library alone

A HELM string is relative to the library it was written with: adding a monomer can turn an inline SMILES
piece into a symbol.


What Is Not a Monomer
---------------------

Nothing is guessed and nothing happens silently:

=====================================================  ============================================
Input                                                  Result
=====================================================  ============================================
a HELM or BILN string that does not parse              ``HelmError`` / ``BilnError``
a bond to an R the monomer lacks, or one already used  ``HelmError`` / ``BilnError``
a HELM position past the end of its chain              ``HelmError``
a symbol the library does not hold                     a marker atom carrying the symbol, logged
a piece matching several monomers                      one written, logged: PEPTIDE before RNA,
                                                       DNA, CHEM, then the latest library's
a BLOB, a repeat range, ``(A,G)``, ``?``                a marker atom carrying the text, logged
a hydrogen pair, an annotation, a polymer group        not kept, logged
a piece of a molecule matching no monomer              an inline SMILES monomer, logged
a piece matching a monomer by constitution only        ``StereoMismatch``
enhanced stereo groups                                 not written, logged
a molecule BILN cannot spell                           ``BilnError``
=====================================================  ============================================

.. testcode::

    log = []
    text = helm(smiles('C[C@H](N)C(=O)NCC(=O)OCC1CCCCC1'), log=log)
    print(text.split('{')[0], sorted({r.rule for r in log}))

.. testoutput::

    PEPTIDE1 ['polymers:inline-monomer']


How a Molecule Is Cut
---------------------

The library decides where a molecule is cut: a bond whose two ends are a pair of attachment atoms two
library monomers can join by, or whose far side is a terminal cap. The pieces are looked up by canonical
key, stereo included, and the covering leaving the fewest atoms unmatched and then using the fewest
monomers wins -- asparagine is one monomer, not an aspartic acid and an amide. Every choice follows the
canonical atom order, which is what makes the string canonical.

=====================================================  ============================================
Rule                                                   Instance
=====================================================  ============================================
never cut between two saturated carbons                an alkyl R3 does not split a side chain
never cut inside a ring of eight atoms or fewer        Tic stays Tic
a monomer is found whole, however many pieces           a PEG monomer cut at every ether bond
a non-ring aromatic bond is single                     ``c1ccccc1c1ccccc1`` is biphenyl, Bip intact
=====================================================  ============================================

A leaving hydrogen in a library record is read as the attachment, so a configuration it references
holds: ``[H:1]/C=C\[H:2]`` is ethylene as a molecule and the cis linker as a monomer.
