Machine learning views
======================

Two views turn a container into ``int32`` numpy arrays: ``state_view`` for a molecule, and
``transition_view`` for a molecule or a mapped reaction. chython emits arrays; the framework layer
repacks them into torch or jax tensors and owns batching, prompts and attention.

One array per column, every one ``int32`` and C-contiguous, so ``torch.from_numpy(view.elements)`` wraps
the buffer chython filled rather than copying it. The views hold arrays and no methods: anything a
framework wants to do to them it does faster in its own tensor library.

numpy is an optional dependency — ``pip install chython[ml]``. Without it these methods raise an
``ImportError`` naming the extra.

Three methods in all, and the third is a repackaging of the second rather than a third kind of view:

======================== ======================== ========================================
method                   container                output
======================== ======================== ========================================
``state_view()``         molecule                 arrays, one side
``transition_view()``    molecule and reaction    arrays, a before and an after side
``modeling_view()``      reaction                 dicts keyed by map number
======================== ======================== ========================================

``transition_view`` is on **both** containers, which is what lets one vocabulary serve both.
``modeling_view`` is the same union assembled into dicts, for a consumer taking one record at a time and
looking an atom up across the arrow; it costs about twice the arrays and reports an underivable hydrogen
count as the physical sentinel rather than as ``0``.

Atom state
----------

``state_view()`` is four per-atom columns, one pairwise matrix, and the token column when the encoding
carries a vocabulary. ``n`` below is the molecule's atom count.

=============== ============ ==========================================================================
attribute       shape        contents
=============== ============ ==========================================================================
``elements``    ``(n,)``     atomic number
``hydrogens``   ``(n,)``     implicit hydrogen count
``neighbors``   ``(n,)``     heavy-atom degree
``distances``   ``(n, n)``   bond counts; ``0`` on the diagonal, ``-1`` for a pair with no path
``tokens``      ``(n,)``     vocabulary id, or ``None`` when the encoding carries no vocabulary
=============== ============ ==========================================================================

.. testcode::

    from chython import smiles

    mol = smiles('CCO')
    view = mol.state_view()
    print(view.elements.tolist())
    print(view.hydrogens.tolist())
    print(view.neighbors.tolist())
    print(view.distances.tolist())
    print(view.elements.shape, view.elements.dtype, view.distances.shape, view.tokens)

.. testoutput::

    [6, 6, 8]
    [3, 2, 1]
    [1, 2, 1]
    [[0, 1, 2], [1, 0, 1], [2, 1, 0]]
    (3,) int32 (3, 3) None

``hydrogens`` is the implicit count and ``neighbors`` is the heavy-atom degree. They stay separate
columns: a single ``degree + hydrogens`` number cannot say which half was unstated.

**A row is a position in container order, never an atom id.** A view carries no ids at all;
``mol.atom_numbers[i]`` is the id of row ``i``, and after a deletion the two stop coinciding:

.. testcode::

    from chython import smiles

    mol = smiles('CCON')
    with mol.edit() as e:
        e.delete_atom(2)

    print(mol.atom_numbers)
    print(mol.state_view().elements.tolist())

.. testoutput::

    [1, 3, 4]
    [6, 8, 7]

``distances`` counts bonds along the shortest path, so it is a graph distance and never a geometry —
a conformer is not consulted. A pair in two different components has no path, which is ``-1``:

.. testcode::

    from chython import smiles

    print(smiles('CCO.O').state_view().distances.tolist())

.. testoutput::

    [[0, 1, 2, -1], [1, 0, 1, -1], [2, 1, 0, -1], [-1, -1, -1, 0]]

Encoding
--------

``TensorEncoding`` is built once and reused. Every default is identity except ``unknown_h``, so a bare
encoding — or none — yields physical values.

.. testcode::

    from chython import TensorEncoding, smiles

    encoding = TensorEncoding(element_shift=2, neighbor_shift=2, distance_shift=2,
                              disconnected=1, max_distance=10)
    view = smiles('[Na+].[Cl-]').state_view(encoding)
    print(view.distances.tolist())

.. testoutput::

    [[2, 1], [1, 2]]

The clamp applies before the shift, and ``disconnected`` is written verbatim: a shifted sentinel would
collide with a real distance, which no consumer can detect.

============================ =====================================================================
argument                     meaning
============================ =====================================================================
``element_shift``            added to the element column
``hydrogen_shift``           added to the hydrogen columns
``neighbor_shift``           added to the degree columns
``distance_shift``           added to a real distance, never to ``disconnected``
``disconnected``             the value for a pair with no path; ``-1`` physically
``unknown_h``                what an underivable hydrogen count reports; ``0`` by default
``max_distance``             clamp on a distance, before the shift; ``0`` is off
``max_neighbors``            clamp on a heavy degree, before the shift; ``0`` is off
``width``                    pad every array to this atom count; ``0`` is off
``pad``                      the value in a padded cell
``pad_diagonal``             the value on a padded row's diagonal; ``0`` is off
``vocabulary``               ``{(element, h_before, n_before, h_after, n_after): token}``
``unknown``                  the token for a key the vocabulary does not hold; ``-1`` by default
============================ =====================================================================

``unknown_h`` is that exception: physically an underivable count is ``H_UNKNOWN``, and the default reports
it as ``0``. The value reaches the vocabulary key as well as the column, so a vocabulary trained under one
``unknown_h`` and read under another misses every atom whose count was unstated. Pass ``unknown_h=15`` to
see the physical sentinel — ``ReactionContainer.modeling_view()`` does exactly that.

.. testcode::

    from chython import TensorEncoding, smiles

    # a nitrogen no valence rule covers, so its hydrogen count is underivable
    mol = smiles('N(F)(F)(F)F')
    print(mol.state_view().hydrogens.tolist())
    print(mol.state_view(TensorEncoding(unknown_h=15)).hydrogens.tolist())

.. testoutput::

    [0, 0, 0, 0, 0]
    [15, 0, 0, 0, 0]

The four fluorines report ``0`` under either setting because their count is derivable and genuinely
zero — which is the reason the sentinel is worth a knob rather than a fixed value.

.. testcode::

    from chython import TensorEncoding, smiles

    view = smiles('CCO').state_view(TensorEncoding(width=5, pad=0, pad_diagonal=1))
    print(view.elements.tolist())
    print(view.distances[4].tolist())

.. testoutput::

    [6, 6, 8, 0, 0]
    [0, 0, 0, 0, 1]

``width`` makes output stackable: every structure the same shape, so ``numpy.stack`` is a
concatenation. ``pad_diagonal`` leaves one unmasked cell on a fully padded row, which an attention
softmax needs. A structure larger than ``width`` raises ``ValueError`` naming its atom count —
truncating to fit is a wrong training example nothing downstream can detect.

Batching
--------

One ``width`` and ``numpy.stack`` is the whole batch path — chython emits one structure at a time and
owns none of it:

.. testcode::

    from numpy import stack
    from chython import TensorEncoding, smiles

    encoding = TensorEncoding(width=6, pad=0, pad_diagonal=1)
    batch = [smiles(s).state_view(encoding) for s in ('CCO', 'c1ccccc1', 'CC(=O)O')]

    elements = stack([v.elements for v in batch])
    distances = stack([v.distances for v in batch])
    print(elements.shape, distances.shape)
    print(elements.tolist())
    print((elements != encoding.pad).sum(axis=1).tolist())

.. testoutput::

    (3, 6) (3, 6, 6)
    [[6, 6, 8, 0, 0, 0], [6, 6, 6, 6, 6, 6], [6, 6, 8, 8, 0, 0]]
    [3, 6, 4]

``elements != pad`` is the attention mask, and that is what ``pad`` is a knob for: the padded value has
to be one no real atom carries. **Element 0 is the R marker**, not an impossibility, so a corpus holding
attachment points wants a ``pad`` outside the element range and a mask taken against that value.

Tokens
------

``tokens`` is the column an embedding layer indexes: one integer per atom in place of the four an atom's
state otherwise takes. A vocabulary is what turns that state into the integer, and the state is five
numbers — every one of them a column the views above already report:

================ =================================================================================
key position     value
================ =================================================================================
``element``      atomic number, ``0..127``
``h_before``     implicit hydrogen count before the transformation, ``0..15``
``n_before``     heavy-atom degree before, ``0..255``
``h_after``      hydrogen count after
``n_after``      degree after
================ =================================================================================

A molecule has no transformation, so its halves are equal: the methyl carbon of ``CCO`` is
``(6, 3, 1, 3, 1)`` — carbon, three hydrogens, one heavy neighbour, unchanged. The key is five numbers
rather than three so that **one vocabulary serves molecules and reactions**; ``transition_view`` is what
fills the two halves differently.

A vocabulary is the caller's — chython compiles the table it is handed and never holds ids — so building
one is a pass over a corpus:

.. testcode::

    from chython import TensorEncoding, smiles

    vocabulary = {}
    for mol in (smiles('CCO'), smiles('CC=O'), smiles('c1ccccc1')):
        view = mol.state_view()
        for e, h, n in zip(view.elements.tolist(), view.hydrogens.tolist(), view.neighbors.tolist()):
            vocabulary.setdefault((e, h, n, h, n), len(vocabulary) + 1)

    for key, token in vocabulary.items():
        print(key, token)

.. testoutput::

    (6, 3, 1, 3, 1) 1
    (6, 2, 2, 2, 2) 2
    (8, 1, 1, 1, 1) 3
    (6, 1, 2, 1, 2) 4
    (8, 0, 1, 0, 1) 5

Five states cover three molecules, and benzene contributes none of them: its aromatic CH is
``(6, 1, 2, 1, 2)``, the state acetaldehyde's carbonyl CH already claimed. **The key is counts and
nothing else** — no aromaticity, no bond orders, no ring membership. A model that has to tell an arene
from an aldehyde reads that off ``distances``, where a ring is a ring, and not off the token.

Read a structure back through that table and ``tokens`` is a column like any other, with ``unknown``
where the corpus held no such state:

.. testcode::

    encoding = TensorEncoding(vocabulary=vocabulary, unknown=999)
    print(smiles('CCO').state_view(encoding).tokens.tolist())
    print(smiles('CC(C)O').state_view(encoding).tokens.tolist())

.. testoutput::

    [1, 2, 3]
    [1, 999, 1, 3]

Isopropanol's central carbon is ``(6, 1, 3, 1, 3)``: three heavy neighbours, a state no molecule in that
three-molecule corpus carries. Its two methyls and its oxygen are states the corpus does carry, so only
one atom of the four falls to ``unknown``.

In a reaction the halves differ, which is the whole reason there are five numbers:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view()
    print(list(zip(view.elements.tolist(), view.h_before.tolist(), view.n_before.tolist(),
                   view.h_after.tolist(), view.n_after.tolist())))

.. testoutput::

    [(6, 3, 1, 3, 1), (35, 0, 1, 0, 0), (8, 1, 0, 1, 1)]

The bromine goes from one heavy neighbour to none and the hydroxide oxygen from none to one, so neither
key is one a molecule can produce: methanol's oxygen is ``(8, 1, 1, 1, 1)`` and the same oxygen in this
transition state is ``(8, 1, 0, 1, 1)``. **A vocabulary trained on molecules holds no key for a reacting
atom** — the halves are equal in every key it learned.

Three rules are worth knowing before a table is compiled, because each one is a wrong training example
with no symptom otherwise. A key outside its field would collide with another key, so the ranges above
are checked and a violation raises ``ValueError`` naming the key. ``(0, 0, 0, 0, 0)`` is refused
outright: it packs to the value that marks an empty slot, and it is a real state — a lone R marker — so
that atom always reads ``unknown``. A token must fit ``int32``.

Tokens are never shifted; the vocabulary owns its own id space. The key is built from physical values
with ``unknown_h`` already applied and **before** the clamps, so an atom that trips ``max_neighbors`` or
``max_distance`` still finds its token — building the key after the clamps would hand back ``unknown``
for every atom in a structure that tripped one.

Transition state
----------------

``transition_view`` reports each atom on both sides of a transformation. For a molecule the two sides
are equal, so molecules and reactions share one layout and one vocabulary. ``n`` below is the union's
atom count and ``m`` its bond count.

================================ ============ ============================================================
attribute                        shape        contents
================================ ============ ============================================================
``elements``                     ``(n,)``     atomic number over the union
``h_before`` / ``h_after``       ``(n,)``     implicit hydrogen count on each side
``n_before`` / ``n_after``       ``(n,)``     heavy-atom degree on each side
``map_numbers``                  ``(n,)``     the map number of each row
``distances``                    ``(n, n)``   bond counts over the union graph
``bonds``                        ``(m, 2)``   the two **row indices** a bond joins, never map numbers
``bond_before`` / ``bond_after`` ``(m,)``     bond order on each side; ``0`` is absent there
``tokens``                       ``(n,)``     vocabulary id, or ``None`` without a vocabulary
================================ ============ ============================================================

.. testcode::

    from chython import smiles

    rxn = smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]')
    view = rxn.transition_view()
    print(view.map_numbers.tolist())
    print(view.n_before.tolist(), view.n_after.tolist())
    print(view.bonds.tolist())
    print(view.bond_before.tolist(), view.bond_after.tolist())

.. testoutput::

    [1, 2, 3]
    [1, 1, 0] [1, 0, 1]
    [[0, 1], [0, 2]]
    [1, 0] [0, 1]

``bond_before`` and ``bond_after`` are the reaction centre: ``(1, 0)`` a broken bond, ``(0, 1)`` a
formed one, ``(1, 2)`` an order change. ``0`` is a bond absent on that side.

A ``bonds`` row holds **row indices into the columns above**, so reading it as map numbers is the one
misreading to avoid — ``[0, 1]`` joins the first two rows, which happen to carry map numbers ``1`` and
``2`` here. ``map_numbers`` is the translation:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view()
    maps = view.map_numbers.tolist()
    for (i, j), before, after in zip(view.bonds.tolist(), view.bond_before.tolist(),
                                     view.bond_after.tolist()):
        print(f'{maps[i]}-{maps[j]}  {before} -> {after}')

.. testoutput::

    1-2  1 -> 0
    1-3  0 -> 1

Two bond rows and one reaction centre: the C-Br bond breaks and the C-O bond forms. A bond present on
neither side does not occupy a row at all.

Atom order over the union is reactant atoms in container order, then product-only atoms, so components
stay contiguous. Agents contribute nothing. ``distances`` spans the union graph, so a bond formed on the
product side shortens a path — here the ``1-3`` bond that only exists after the reaction puts the oxygen
one bond from the carbon:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').transition_view()
    print(view.distances.tolist())

.. testoutput::

    [[0, 1, 1], [1, 0, 2], [1, 2, 0]]

An atom that bonds to nothing on either side is ``disconnected`` from every other row — a spectator ion
is one matrix row of ``-1``, and one matrix covers the union rather than a per-side pair of them:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][Br:2].[Na+:4].[OH-:3]>>[CH3:1][OH:3].[Br-:2].[Na+:4]').transition_view()
    print(view.map_numbers.tolist())
    print(view.distances.tolist())

.. testoutput::

    [1, 2, 4, 3]
    [[0, 1, -1, 1], [1, 0, -1, 2], [-1, -1, 0, -1], [1, 2, -1, 0]]

The sodium is the third row, ``-1`` everywhere off its own diagonal, and the map numbers show the union
order: reactant atoms in container order — the sodium among them — then whatever only the products carry.

A badly mapped record is reported and not refused. An atom the record does not number is **placed on the
side it came from** — reactant-only or product-only, which is all an unmapped atom can be — and counted
in ``unmapped``. ``map_numbers`` carries ``0`` for its row, that column being the record's own numbering:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][C:2](=[O:3])[OH:4].CO>>[CH3:1][C:2](=[O:3])[O:4]C.O').transition_view()
    print(view.unmapped)
    print(view.map_numbers.tolist())
    print(view.n_before.tolist())
    print(view.n_after.tolist())

.. testoutput::

    {'reactants': 2, 'products': 2}
    [1, 2, 3, 4, 0, 0, 0, 0]
    [1, 3, 1, 1, 1, 1, 0, 0]
    [1, 3, 1, 2, 1, 1, 1, 0]

Map ``4``, the acid's hydroxyl oxygen, reads one heavy neighbour before and two after: the arriving methyl
is a row, so its bond is a bond. Dropping the unnumbered atoms instead would leave that ``2`` a ``1`` and
would make an aryl chloride, bromide and iodide the same transition state.

The price is on the other side of the same record. Nothing says methanol's carbon and oxygen are the
product's methyl and water, so they are four rows rather than two — the reactant pair leaves and the
product pair arrives. ``unmapped`` is what a consumer that cannot accept that reads, and it is nonzero
here for exactly that reason.

``collisions`` lists a map number claimed twice on one side. The first claim keeps its union row; the
second contributes no state and no bond, because a union cannot hold two atoms at one key.

``ReactionContainer.modeling_view()`` hands this same union back as two dicts, with ``unmapped`` and
``collisions`` unchanged:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]').modeling_view()
    print(view.states)
    print(view.union_bonds)

.. testoutput::

    {1: (6, 3, 1, 3, 1), 2: (35, 0, 1, 0, 0), 3: (8, 1, 0, 1, 1)}
    {(1, 2): (1, 0), (1, 3): (0, 1)}

``states`` is ``{map_number: (element, h_before, n_before, h_after, n_after)}`` — the transition-view key,
one per atom — and ``union_bonds`` is ``{(n, m): (order_before, order_after)}`` with ``n < m``. It takes no
encoding: the hydrogen count is ``H_UNKNOWN`` where the record does not state one, so nothing here is
shifted, clamped or padded.

A dict needs one key per atom, so an unmapped one is keyed **negatively**, ``-1`` onwards in union order.
Zero would put every such atom on one entry and a positive key would be indistinguishable from a number
the record stated:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1]Br>>[CH3:1]O').modeling_view()
    print(view.states)
    print(view.union_bonds)

.. testoutput::

    {1: (6, 3, 1, 3, 1), -1: (35, 0, 1, 0, 0), -2: (8, 1, 0, 1, 1)}
    {(-1, 1): (1, 0), (-2, 1): (0, 1)}

The record names one atom and the transition state is still whole: the C–Br bond breaks, the C–O bond
forms, and the carbon counts one heavy neighbour on each side.

What a record can state
-----------------------

Every difference between the two sides lands in one of these six places, and reading a transition state
is reading those columns against each other:

======================== ===============================================================================
the record states        where the view puts it
======================== ===============================================================================
an atom on both sides    one row; the halves differ wherever that atom changed
a leaving atom           one row, present on the reactant side only; each of its bonds reads
                         ``(order, 0)``
an arriving atom         the mirror, each bond reading ``(0, order)``
a hydrogen count change  ``h_before != h_after`` on that row
a degree change          ``n_before != n_after`` on that row
a bond order change      one bond row carrying two orders, neither of them ``0``
======================== ===============================================================================

Being on one side only is not the same thing as being unnumbered. An atom is reactant-only when the
products do not carry it — which is most records, since a record is free to omit its byproducts — and
whether it was numbered changes its key and nothing else. The three columns below are identical for
``[CH3:1][Br:2]>>[CH4:1]`` and ``[CH3:1]Br>>[CH4:1]``; only ``map_numbers`` and ``unmapped`` differ.

An oxidation moves hydrogens and one bond order and nothing else:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]').transition_view()
    print(view.h_before.tolist(), view.h_after.tolist())
    print(view.n_before.tolist(), view.n_after.tolist())
    print(view.bond_before.tolist(), view.bond_after.tolist())
    print(view.unmapped)

.. testoutput::

    [3, 2, 1] [3, 1, 0]
    [1, 2, 1] [1, 2, 1]
    [1, 1] [1, 2]
    {'reactants': 0, 'products': 0}

The carbon and the oxygen each lose a hydrogen and the bond between them goes single to double. No
degree moves, because nothing left and nothing arrived.

A leaving atom keeps the state it had in the fragment it left with: ``n_after`` counts the neighbours
that left beside it, while those same bonds read ``(order, 0)``. The degree columns describe the atom as
its own side drew it — the byproduct is still a molecule — and the bond columns describe the union:

.. testcode::

    from chython import smiles

    view = smiles('[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])([CH3:7])[CH3:8].[OH2:9]'
                  '>>[CH3:1][C:2](=[O:3])[OH:9]').transition_view()
    print(view.map_numbers.tolist())
    print(view.n_before.tolist())
    print(view.n_after.tolist())

.. testoutput::

    [1, 2, 3, 4, 5, 6, 7, 8, 9]
    [1, 3, 1, 2, 4, 1, 1, 1, 0]
    [1, 3, 1, 1, 4, 1, 1, 1, 1]

Map ``4``, the ester oxygen, drops from two heavy neighbours to one: the bond to the carbonyl carbon
breaks and the bond to the tert-butyl carbon leaves with it. Map ``5`` keeps all four, all four of its
neighbours having left with it. Map ``9``, the water oxygen, gains one and loses a hydrogen.

One record reaches all of it — a leaving atom, an arriving one, a changing one, a bond broken and a bond
formed:

.. testcode::

    from chython import smiles

    view = smiles('[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1Br.[NH:7]1[CH2:8][CH2:9]1'
                  '>>[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[N:7]1[CH2:8][CH2:9]1.Br').modeling_view()
    print(view.unmapped)
    print(view.states[6], view.states[7])
    print({k: v for k, v in view.union_bonds.items() if v[0] != v[1]})

.. testoutput::

    {'reactants': 1, 'products': 1}
    (6, 0, 3, 0, 3) (7, 1, 2, 0, 3)
    {(-1, 6): (1, 0), (6, 7): (0, 1)}

The bromide leaves, hydrogen bromide arrives, and the aziridine nitrogen loses its hydrogen and gains a
neighbour. The aryl carbon is the row to check: three heavy neighbours before and three after, the
bromine counted on one side and the nitrogen on the other. Neither bromine is numbered, so nothing pairs
them — ``-1`` leaves and ``-2`` arrives, and ``unmapped`` counts one on each side.

Measured cost
-------------

``python -m chython.core.test.bench_ml``, single thread, public compounds.  2000 molecules, median 14
atoms:

======================================== ============
path                                     time
======================================== ============
``distance_matrix()`` alone              0.80 µs/mol
``state_view()``, physical               1.60 µs/mol
``state_view()``, shifted and clamped    1.57 µs/mol
``state_view()``, padded to 64           3.16 µs/mol
``state_view()``, with a vocabulary      1.86 µs/mol
``transition_view()``, molecule          3.27 µs/mol
======================================== ============

500 mapped reactions, 12 atoms per reaction over both sides:

======================================== ============
path                                     time
======================================== ============
``transition_view()``, reaction          3.18 µs/rxn
``modeling_view()``                      7.09 µs/rxn
======================================== ============

Starting from a stored pach record rather than a live container, by wire version:

============ ============== ============
wire version median record  time
============ ============== ============
``2``        173 B          5.18 µs/mol
``3``        130 B          4.76 µs/mol
``4``        130 B          4.65 µs/mol
============ ============== ============

Version 3 shortened the record by a quarter and the three decode at the same speed within measurement
noise. On this corpus ``unpach`` accounts for ~3.0 µs of that ~4.7 µs/mol, against ~1.6 µs for the array
derivation — so the byte path is dominated by decoding, and caching arrays pays where re-deriving them
from stored records does not.

Distance is the dominant term in every path: an every-source BFS over the CSR, ``O(V·E)``.
