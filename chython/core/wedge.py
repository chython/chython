# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""Geometry to parity, and parity back to wedges.

A drawn structure does not store configurations, it stores a *drawing*: 2D coordinates plus a wedge
or hash on one bond per stereocentre, or 3D coordinates and no wedges at all.  Turning that into the
core's parity is this module's read direction; choosing a drawing that reproduces a stored parity is
its write direction.  Both go through :func:`tetrahedral_parity`, which is why a round trip cannot
invert a sign: the writer picks the wedge by asking the reader's own function which one it would read
back.

**IN ``core`` AND NOT IN A FORMAT, because a wedge is not a file's idea.**  Putting it in
``formats/ctfile/`` -- on the reasoning that a CTfile is where wedges are written down -- makes every
other caller cross the layering: ``formats/xml/_cml.py`` and ``_mrv.py`` reach sideways into another
dialect's private module, and ``depict/svg.py`` reaches *up* into ``formats``, a renderer importing a
file format to find out which way to draw a triangle.  What the module contains is signed volumes,
permutation parities and the core's own stereo-unit frame.  There is nothing in it a CTfile knows and a
CML document or an SVG canvas does not.  So it sits below all of them, where a format and a renderer may
both see it, and the file formats keep only what is genuinely theirs: MDL's integer atom-parity field,
MRV's and CML's ``<bondStereo>`` letters, and each version's field layout.

It is Python and not a ``.pxi``: nothing here is on a measured hot path, and a pure-Python leaf in
``core`` -- ``reaction.py`` is the precedent -- buys the layering without adding a compile step to
the geometry.  Its tests stay in ``formats/ctfile/test/test_wedge.py``, which is deliberate: they
exercise the chooser through real drawn records, and a corpus of molfiles cannot be read from
``core/test/`` without ``core`` importing a format to test itself.

TWO CONVENTIONS ARE PINNED HERE, and neither is a free choice.

**Frame (ruling F26).**  A parity is meaningless without an ordered list of directions to measure it
against, and the stored byte is relative to the order the core itself reports in
``stereo_units()[k]['refs']`` -- heavy neighbours in CSR-slot ascending order, then explicit
hydrogens ascending, then a hole for each direction with no atom of its own.  Nothing here invents
an order, sorts an order, or assumes that creation order, file order and CSR order coincide.  They
often do, which is exactly why assuming it survives testing and then fails on a real file.

**Sign.**  Even (1) or odd (2) is a label, and which geometric handedness gets which label is fixed
by the core's own legacy contract: ``set_stereo(True)`` is documented as parity 2, and chython 2's
``_stereo`` boolean is True for an *anticlockwise* environment -- SMILES ``@``.  Measured, not
assumed: for ``F[C@](Cl)(Br)I`` chython 2 reports ``_stereo is True`` with environment
``(F, Cl, Br, I)``, and the signed volume of that environment with F wedged toward the viewer is
negative.  So::

    signed volume < 0  ==  anticlockwise  ==  SMILES @   ==  parity 2 (odd)
    signed volume > 0  ==  clockwise      ==  SMILES @@  ==  parity 1 (even)

``test_wedge.py`` pins both directions against the chython 2 stack over a real file, because a
reader-writer-reader round trip passes with the sign globally inverted and is therefore not evidence.
"""

from ._core import WEDGE_DOWN, WEDGE_EITHER, WEDGE_NONE, WEDGE_UP
from ._log import LogRecord, LOST, REFUSED, REPAIRED


__all__ = ['assign_parities', 'tetrahedral_parity', 'cis_trans_parity', 'cis_trans_frame',
           'cis_trans_letter', 'cis_trans_for_write', 'allene_parity', 'atropisomer_parity',
           'stated_parity', 'stated_cis_trans', 'wedges_for_write', 'wedge_in_file_order',
           'signed_volume', 'SU_TETRA', 'SU_CIS_TRANS', 'SU_ALLENE', 'SU_ATROPISOMER']


# The core's stereo-unit kinds.  Mirrored rather than imported because they are `DEF` constants in
# `_stereo.pxi` and so exist only at Cython compile time; `test_wedge.py` asserts the values against
# the kinds the core actually emits for a tetrahedron, an alkene, an allene and a biaryl.
SU_TETRA = 0
SU_CIS_TRANS = 1
SU_ALLENE = 2
SU_ATROPISOMER = 3

# The out-of-plane displacement a wedge stands for.  The magnitude is arbitrary -- the determinant is
# linear in it, so the sign of the answer does not depend on the number -- but 1.0 is what chython 2
# uses, and keeping it identical means the oracle test compares two computations of the same
# quantity rather than two quantities that happen to agree.
_WEDGE_Z = {WEDGE_UP: 1.0, WEDGE_DOWN: -1.0, WEDGE_NONE: 0.0}


def signed_volume(v0, v1, v2, v3):
    """``det(v1 - v0, v2 - v0, v3 - v0)`` -- six times the signed volume of the tetrahedron.

    Positive for a clockwise ``(v1, v2, v3)`` seen from ``v0``, negative for anticlockwise.
    """
    ax, ay, az = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
    bx, by, bz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
    cx, cy, cz = v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]
    return ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)


def _plane_sign(a, b, c):
    """z of ``(b - a) x (c - b)`` -- which side of the directed line ``a->b`` the point `c` lies on."""
    return (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])


def tetrahedral_parity(mol, unit, z=None, wedge_of=None, log=None):
    """Parity of a tetrahedral `unit` from the drawing: 0 unset, 1 even, 2 odd.

    `z` maps stable id to a third coordinate, for a genuinely 3D file; when it is given and not flat
    the wedges are ignored, because a 3D drawing states the configuration outright and a wedge on top
    of it is at best redundant.  `wedge_of` overrides where the wedge codes come from, which is what
    lets the writer ask "what would I read back if I drew *this*".

    Returns 0 -- deliberately, not an exception -- when the drawing says nothing: no wedge at all, a
    wavy "either" bond, a degenerate layout, or two undrawn directions.  Every one of those is a real
    statement by the file and none of them is an error in it.
    """
    anchor = unit['anchor']
    refs = unit['refs']
    if wedge_of is None:
        wedge_of = mol.wedge_of
    z = z or {}
    flat = not any(z.values())

    cx, cy = mol.xy_of(anchor)
    cz = z.get(anchor, 0.0)
    vectors = []
    holes = []
    drawn = False
    for i, r in enumerate(refs):
        if r is None:
            vectors.append(None)
            holes.append(i)
            continue
        x, y = mol.xy_of(r)
        if flat:
            w = wedge_of(anchor, r)
            if w == WEDGE_EITHER:
                if log is not None:
                    log.append(LogRecord('wedge:either-bond', (anchor, r),
                                         f'atom {anchor}: bond to {r} drawn as either, configuration '
                                         f'left unset', LOST))
                return 0
            dz = _WEDGE_Z.get(w, 0.0)
            if dz:
                drawn = True
        else:
            dz = z.get(r, 0.0) - cz
            drawn = True
        vectors.append((x - cx, y - cy, dz))

    if not drawn:
        return 0
    if len(holes) > 1:
        if log is not None:
            log.append(LogRecord('wedge:undrawn-directions', (anchor,),
                                 f'atom {anchor}: {len(holes)} undrawn directions, configuration left unset', LOST))
        return 0
    if holes:
        # The one direction with no atom of its own -- an implicit hydrogen, or a lone pair on a
        # sulfoxide sulfur -- sits opposite the sum of the drawn ones.  That is the same statement as
        # "use the stereocentre itself as the fourth vertex", up to a positive factor: for vectors
        # relative to the centre, det(v0, v1, v2) and the four-point volume with the phantom at
        # -(v0 + v1 + v2) differ by exactly -4, so they never disagree about the sign.
        sx = sy = sz = 0.0
        for v in vectors:
            if v is not None:
                sx += v[0]
                sy += v[1]
                sz += v[2]
        vectors[holes[0]] = (-sx, -sy, -sz)

    volume = signed_volume(*vectors)
    if volume == 0.0:
        if log is not None:
            log.append(LogRecord('wedge:degenerate-geometry', (anchor,),
                                 f'atom {anchor}: degenerate geometry, configuration left unset', LOST))
        return 0
    return 2 if volume < 0.0 else 1


def cis_trans_parity(mol, unit, log=None, plane=None):
    """Parity of a cis/trans `unit` from 2D coordinates: 0 unset, 1 even, 2 odd.

    A CTfile has no field for double-bond geometry -- bond stereo 3 says only "cis or trans, unknown
    which" -- so the coordinates are the whole statement, and a file drawn without meaningful
    coordinates has not stated one.

    ``refs`` for a bond kind is two pairs, the anchor's end first, and a hole is always the second
    slot of its pair, so ``refs[0]`` and ``refs[2]`` are the two named substituents that define the
    frame.  ``refs[0]`` and ``refs[2]`` on the same side is parity 2, again by the legacy alignment:
    chython 2 sets its bond ``_stereo`` True for a same-side first pair.

    `plane` is ``{stable id: (x, y)}``, and it is the same parameter :func:`wedges_for_write` takes,
    for the same reason its paragraph gives: a caller holding a layout of its own -- a renderer -- is
    asking about THAT drawing and not about whatever coordinates the molecule happens to store.  The
    question matters here because a plain double bond asserts whatever geometry the drawing shows, so a
    renderer has to know which configuration its own plane reads back before it draws one; and it
    matters that the answer come from HERE, because a renderer that works the sign out for itself is a
    second copy of the convention this module exists to keep single.
    """
    if plane is not None:
        mol = _Planar(mol, plane)
    # The frame comes from `cis_trans_frame`, which is also what a writer names in an `atomRefs4`: the
    # partner terminal is the anchor's neighbour that also neighbours the far substituent, found rather
    # than assumed, because the double bond's order is not a reliable marker here (an atropisomeric or
    # cumulated system can put more than one candidate bond on the anchor).  Shared so that the frame a
    # parity is read in and the frame a descriptor is written in are one expression.
    frame = cis_trans_frame(mol, unit)
    if frame is None:
        anchor = unit['anchor']
        refs = unit['refs']
        if log is not None and refs[0] is not None and refs[2] is not None:
            # Only the missing partner earns a line.  An unnamed substituent is an ordinary unnamed
            # substituent and the caller can see it in `refs`; no partner terminal on a unit the core
            # emitted as cis/trans is a graph the reader did not expect.
            log.append(LogRecord('wedge:no-partner-terminal', (anchor,),
                                 f'bond stereo at atom {anchor}: partner terminal not found, left unset', LOST))
        return 0
    near, anchor, partner, far = frame

    a = mol.xy_of(near)
    b = mol.xy_of(anchor)
    c = mol.xy_of(partner)
    d = mol.xy_of(far)
    s1 = _plane_sign(a, b, c)
    s2 = _plane_sign(b, c, d)
    if s1 == 0.0 or s2 == 0.0:
        if log is not None:
            log.append(LogRecord('wedge:collinear', (anchor,),
                                 f'bond stereo at atom {anchor}: collinear layout, left unset', LOST))
        return 0
    return 2 if s1 * s2 > 0.0 else 1


def cis_trans_frame(mol, unit):
    """``(near, anchor, partner, far)`` -- the four atoms a cis/trans parity is measured over.

    ``None`` when the unit names no frame, which is the same condition :func:`cis_trans_parity` returns
    0 for: an unnamed substituent on either end, or no partner terminal.  Factored out of that function
    rather than duplicated beside it, so the frame a parity was *read* in and the frame a descriptor is
    *written* in cannot come apart -- which is the failure this module exists to prevent, one kind
    further along than the wedges.
    """
    refs = unit['refs']
    near, far = refs[0], refs[2]
    if near is None or far is None:
        return None
    anchor = unit['anchor']
    for m in mol.neighbors_of(anchor):
        if far in mol.neighbors_of(m):
            return near, anchor, m, far
    return None


def cis_trans_for_write(mol, framed=True):
    """``[(anchor, letter, frame), ...]`` -- every cis/trans configuration a document can state as a letter.

    The molecule-side half of writing ``<bondStereo>C``/``T``, so that CML and MRV do not each walk the
    stereo units and each decide what a writable descriptor is.  What is left to a dialect is the element
    name and whether it spells the frame out, which is `framed`: see :func:`cis_trans_letter` for why a
    bare letter is not writable everywhere a framed one is.

    The anchors are also what a writer hands :func:`wedges_for_write` as ``cis_trans_stated``, which is
    the other reason this is one list and not a loop in each writer -- the descriptor it wrote and the
    loss it suppressed have to be the same set.
    """
    out = []
    for unit in mol.stereo_units():
        if unit['kind'] != SU_CIS_TRANS:
            continue
        said = cis_trans_letter(mol, unit, framed=framed)
        if said is not None:
            out.append((unit['anchor'], said[0], said[1]))
    return out


def cis_trans_letter(mol, unit, framed=True):
    """``('C'|'T', frame)`` for a configured cis/trans `unit`, or ``None``.

    THE ONE PLACE A PARITY BECOMES A NON-GEOMETRIC DESCRIPTOR, and it is in ``core`` for the same
    reason the wedge chooser is: CML and MRV both spell a double-bond configuration as a letter, a
    third dialect will, and two of them working it out separately is two chances to disagree about a
    sign.  What a *format* owns is the element name and whether it states the frame; the sign is
    chemistry.

    ``C`` is parity 2 and ``T`` is parity 1, which is :func:`cis_trans_parity`'s own convention read
    backwards -- ``refs[0]`` and ``refs[2]`` on the same side is parity 2, and same-side is cis.
    Measured rather than taken on trust, because a sign convention asserted in a docstring is how the
    wedges went wrong: over this repository's V2000/SDF corpus, 122 drawn cis/trans bonds agree with
    this mapping and 0 disagree, comparing against the plane-sign arithmetic directly.  The negation
    would score 0 of 122, so the corpus decides rather than tolerating both.  Independently anchored
    at the other end by the SMILES reader, which has no coordinates at all: ``C/C=C/C`` is E and comes
    back parity 1, so ``T``.

    **The frame is returned with the letter and is not optional.**  A bare letter states a
    configuration without saying of what -- ChemAxon writes one, and gets away with it because the
    drawing is in the same record -- so a caller with no drawing to fall back on must write the frame
    or write nothing.  Handing them back together is what stops a writer from emitting the first and
    forgetting the second.

    `framed` is ``False`` for a dialect whose spelling has nowhere to put the frame -- MRV is one -- and
    then ``None`` is returned for a bond where a bare letter would be ambiguous.  **The test is
    :func:`stated_cis_trans`'s own**, and it is here rather than in the writer for the reason this whole
    module exists: what a reader can take back and what a writer may state have to be one condition, or
    the round trip loses a descriptor that looked writable.
    """
    parity = mol.parity_of(unit['anchor'])
    if not parity:
        return None
    frame = cis_trans_frame(mol, unit)
    if frame is None:
        return None
    if not framed and (mol.degree_of(frame[1]) > 2 or mol.degree_of(frame[2]) > 2):
        return None
    return ('C' if parity == 2 else 'T'), frame


def stated_cis_trans(mol, unit, letter, refs=(), log=None):
    """Core parity from a **non-geometric** double-bond descriptor, or 0 when it states nothing usable.

    :func:`cis_trans_letter` read backwards, and it is the cis/trans counterpart of
    :func:`stated_parity`: the second, coordinate-free statement a *document* can make about a double
    bond.  No CTAB can make it -- bond stereo 3 says only "cis or trans, unknown which" -- but CML and
    MRV both spell it ``<bondStereo>C``/``T``, so a record with no meaningful layout still states a
    configuration there, and refusing to read it would drop a fact the file states plainly.

    `letter` is ``'C'`` or ``'T'``; `refs` is the four atoms the document measured it in, as stable ids
    in the document's own order ``(outer, terminal, terminal, outer)``, or empty when the document names
    none.

    **A BARE LETTER IS READ ONLY WHERE IT CANNOT BE AMBIGUOUS.**  ``C`` says two substituents are on one
    side without saying which two, and where either terminal carries a second substituent there are two
    answers and the file has chosen neither.  ChemAxon writes bare letters and gets away with it because
    the drawing is in the same record; a reader with no drawing cannot.  So the letter is taken when both
    terminals carry exactly one substituent apiece -- one frame exists, so naming it adds nothing -- and
    reported and dropped otherwise.

    Where `refs` is given the sign is **translated, not assumed**: a document is free to measure ``C``
    over the *other* substituent at either end, and each such swap inverts the parity.  A frame that does
    not name this double bond, or names an atom that is not a substituent of the terminal it is written
    beside, is a frame describing something else and is dropped with a line rather than applied to the
    bond that happens to be here.
    """
    frame = cis_trans_frame(mol, unit)
    if frame is None:
        return 0
    near, anchor, partner, far = frame
    parity = 2 if letter == 'C' else 1 if letter == 'T' else 0
    if not parity:
        if log is not None:
            log.append(LogRecord('wedge:bad-letter', (anchor,),
                                 f'bond stereo at atom {anchor}: configuration {letter!r} is not a '
                                 f'configuration, left unset', LOST))
        return 0

    if not refs:
        if mol.degree_of(anchor) > 2 or mol.degree_of(partner) > 2:
            if log is not None:
                log.append(LogRecord('wedge:ambiguous-bare-letter', (anchor,),
                                     f'bond stereo at atom {anchor}: configuration {letter} names no reference '
                                     f'atoms and a terminal carries two substituents, so which pair is {letter} '
                                     f'is not stated; left unset', LOST))
            return 0
        return parity

    if len(refs) != 4:
        if log is not None:
            log.append(LogRecord('wedge:wrong-ref-count', (anchor,),
                                 f'bond stereo at atom {anchor}: configuration {letter} measured over '
                                 f'{len(refs)} atoms, expected 4; left unset', LOST))
        return 0
    o1, t1, t2, o2 = refs
    if (t1, t2) == (partner, anchor):
        o1, t1, t2, o2 = o2, t2, t1, o1   # the document wrote the frame from the other terminal
    if (t1, t2) != (anchor, partner):
        if log is not None:
            log.append(LogRecord('wedge:wrong-bond', (anchor,),
                                 f'bond stereo at atom {anchor}: configuration {letter} is measured over the '
                                 f'bond {t1}={t2}, which is not this one; left unset', LOST))
        return 0
    for outer, terminal, ours in ((o1, anchor, near), (o2, partner, far)):
        if outer == ours:
            continue
        if outer == (partner if terminal == anchor else anchor) or \
                outer not in mol.neighbors_of(terminal):
            if log is not None:
                log.append(LogRecord('wedge:wrong-substituent', (anchor, outer, terminal),
                                     f'bond stereo at atom {anchor}: configuration {letter} names atom {outer} '
                                     f'at the {terminal} end, which is not a substituent there; left unset', LOST))
            return 0
        parity = 3 - parity   # the other substituent at one end is the opposite configuration
    return parity


def _chain_ends(mol, anchor):
    """The two terminals of the cumulene chain through `anchor`, as ``[(terminal, inward), ...]``.

    WALKED, NOT LOOKED UP, and that is the whole point of this function.  An allene's unit is anchored
    on the *centre* of its chain, and the core emits ``SU_ALLENE`` for every odd chain length -- five
    and seven cumulated carbons included -- so on anything longer than three the terminals are not the
    anchor's neighbours.  A reader that took ``neighbors_of(anchor)`` for the terminals is correct on
    every allene and silently wrong on every longer cumulene, which is a bug an allene-only test file
    cannot see.

    ``inward`` is the terminal's chain neighbour, needed because a terminal's undrawn direction sits
    opposite the sum of its drawn ones and the chain bond is one of those.

    ``order_of(...) == 2`` is the same predicate the core walks with: the arena stores an aromatic bond
    as order 4 with its aromatic flag set, so a Kekule double bond is exactly order 2 here.
    """
    ends = []
    for start in mol.neighbors_of(anchor):
        if mol.order_of(anchor, start) != 2:
            continue
        prev, cur = anchor, start
        while True:
            nxt = [m for m in mol.neighbors_of(cur) if m != prev and mol.order_of(cur, m) == 2]
            if len(nxt) != 1:  # a terminal, or a branch this pass has no business resolving
                break
            prev, cur = cur, nxt[0]
        ends.append((cur, prev))
    return ends


def allene_parity(mol, unit, z=None, wedge_of=None, log=None):
    """Parity of an allene `unit` from the drawing: 0 unset, 1 even, 2 odd.

    An allene is *axially* chiral: the four positions that make up its frame sit on two different
    atoms, one pair per terminal, and the two pairs are mutually perpendicular.  So the sign is the
    signed volume of those four positions in the core's ``refs`` order -- the identical rule
    :func:`tetrahedral_parity` applies to a tetrahedron's four -- and no separate convention is
    introduced here.  ``core/_inchi.pxi`` states that equivalence from the other side, deriving both of
    its frames from one signed-volume rule.

    WHAT THE WEDGE MEANS ON AN ALLENE, and why it is not the tetrahedral reading.  The wedge sits on a
    bond from a *terminal*, while the unit is anchored on the *centre*, so the two are never on the
    same atom.  In a flat drawing of an allene the wedged terminal's two substituents lie in a plane
    perpendicular to the paper: one toward the viewer and one away, and their drawn in-plane offsets
    are the drawing's lie about a molecule that is not planar.  This function therefore *discards*
    those offsets, collapsing both substituents onto the terminal's own position with ``z = +1`` and
    ``z = -1``, and keeps the far terminal's drawn positions with ``z = 0``.  Two wedges on that
    terminal -- up on one bond, down on the other, which is how drawing packages state
    "perpendicular to the paper" -- are consistent with a single one and produce the same two z values;
    two pointing the *same* way are not, and are reported.

    Analytically the determinant then collapses to ``(m1 - m0) * cross_z(f0 - t, f1 - t)``, which
    depends on the wedged terminal only through the sign of its wedge -- the same quantity chython 2's
    ``_allene_sign`` computes from the mark, the axis and one far substituent.  The two stacks are
    compared centre by centre over this repository's corpus in ``test_wedge.py``.

    `z` maps stable id to a third coordinate; when it is given and not flat the wedges are ignored and
    the real positions are used, exactly as in the tetrahedral case.

    Returns 0 -- never raises -- for every drawing that states no configuration: no wedge, a wavy
    "either" bond, wedges on both terminals (two statements about a frame that takes one, the second
    of them measured against a plane the first has already tilted), two wedges pointing the same way,
    or a layout whose determinant vanishes.  There is no second source to fall back on the way a
    tetrahedral centre has the atom parity field: **no CTfile field states an axial configuration**, so
    an unreadable drawing means the molecule has none, and saying so is the whole of the contract.
    """
    anchor = unit['anchor']
    refs = unit['refs']
    if wedge_of is None:
        wedge_of = mol.wedge_of
    z = z or {}

    ends = _chain_ends(mol, anchor)
    if len(ends) != 2:  # pragma: no cover - the core anchors SU_ALLENE on a two-ended chain
        if log is not None:
            log.append(LogRecord('wedge:allene-chain-ends', (anchor,),
                                 f'atom {anchor}: allene chain has {len(ends)} ends, configuration left unset', LOST))
        return 0
    # Which end owns which ref pair is asked of the constitution rather than assumed from the chain
    # walk's order: `refs[0:2]` belongs to the terminal nearer the anchor's lower CSR slot, and the
    # walk returns the two ends in neighbour order, which is not the same statement.
    inward = dict(ends)
    near = [t for t in inward if refs[0] in mol.neighbors_of(t)]
    far = [t for t in inward if refs[2] in mol.neighbors_of(t)]
    if len(near) != 1 or len(far) != 1 or near[0] == far[0]:  # pragma: no cover - F47 forbids it
        if log is not None:
            log.append(LogRecord('wedge:allene-terminals', (anchor,),
                                 f'atom {anchor}: allene terminals not identified, configuration left unset', LOST))
        return 0
    pairs = ((near[0], 0), (far[0], 2))
    for _, base in pairs:
        # A hole -- an implicit hydrogen or a lone pair -- has a position this can compute; an EMPTY
        # slot is a direction that does not exist, and three points define no volume.  The unnamed mask
        # is the only thing that tells the two apart.
        # Untested by construction rather than by omission: an sp2 terminal has three directions, so
        # its ref pair is always full or holed, never short.  The guard stays because reaching
        # `signed_volume` with a None would be a crash rather than a report.
        if refs[base + 1] is None and not (unit['unnamed_mask'] >> (base + 1)) & 1:
            if log is not None:
                log.append(LogRecord('wedge:allene-one-direction', (anchor,),
                                     f'atom {anchor}: an allene terminal has one direction, left unset', LOST))
            return 0

    pos, reason, named = _axial_frame(mol, unit, pairs, inward, z, wedge_of)
    if reason == 'either':
        if log is not None:
            t, r = named
            log.append(LogRecord('wedge:either-bond', (anchor, t, r),
                                 f'atom {anchor}: bond {t}-{r} drawn as either, configuration '
                                 f'left unset', LOST))
        return 0
    if reason == 'flat':
        return 0  # a flat allene: unspecified, and the common case in real files
    if reason == 'both-ends':
        if log is not None:
            log.append(LogRecord('wedge:allene-both-terminals', (anchor,),
                                 f'atom {anchor}: wedges on both terminals of the allene, configuration '
                                 f'left unset', LOST))
        return 0
    if reason == 'same-direction':
        if log is not None:
            log.append(LogRecord('wedge:allene-same-direction', (anchor,),
                                 f'atom {anchor}: both bonds of one allene terminal wedged the '
                                 f'same way, configuration left unset', LOST))
        return 0

    volume = signed_volume(*pos)
    if volume == 0.0:
        if log is not None:
            log.append(LogRecord('wedge:degenerate-geometry', (anchor,),
                                 f'atom {anchor}: degenerate allene layout, configuration left unset', LOST))
        return 0
    return 2 if volume < 0.0 else 1


def _axial_frame(mol, unit, pairs, inward, z, wedge_of):
    """The four positions an axial parity is measured over: ``(pos, reason, atoms)``.

    ONE GEOMETRY FOR THE TWO AXIAL KINDS, and it reports rather than logs.  An allene's frame sits on
    two chain terminals and an atropisomer's on two ring pivots; the construction is the same signed
    volume over the same ``refs`` order, but a log line has to name the noun it is about, so `reason`
    -- ``'either'``, ``'flat'``, ``'both-ends'``, ``'same-direction'`` -- is worded by the caller and
    `atoms` carries what its message needs.  `pos` is ``None`` whenever `reason` is not.

    `pairs` is ``((end, base), (end, base))`` with the anchor's pair first, and `inward` maps each end
    to its neighbour along the axis.  A hole in a pair is filled from the drawing by
    :func:`_fill_terminal`; a wedged end is COLLAPSED onto its own position with the wedge's sign,
    because in a flat drawing of an axial system that end's two directions are perpendicular to the
    paper and their drawn offsets are the drawing's lie about them.
    """
    refs = unit['refs']
    pos = [None] * 4
    if any(z.values()):
        for t, base in pairs:
            _fill_terminal(mol, refs, base, t, inward[t], pos, z)
        return pos, None, ()

    marks = {}
    for t, base in pairs:
        for j in (base, base + 1):
            r = refs[j]
            if r is None:
                continue
            w = wedge_of(t, r)
            if w == WEDGE_EITHER:
                return None, 'either', (t, r)
            if _WEDGE_Z.get(w, 0.0):
                marks[j] = _WEDGE_Z[w]
    if not marks:
        return None, 'flat', ()
    if len({0 if j < 2 else 2 for j in marks}) != 1:
        return None, 'both-ends', ()
    wedged = 0 if next(iter(marks)) < 2 else 2
    for t, base in pairs:
        if base == wedged:
            m0, m1 = marks.get(base), marks.get(base + 1)
            if m0 is None:
                m0 = -m1
            elif m1 is None:
                m1 = -m0
            elif m0 == m1:
                return None, 'same-direction', ()
            x, y = mol.xy_of(t)
            # The collapse: the drawn offsets at this end say nothing, the wedge says everything.
            pos[base] = (x, y, m0)
            pos[base + 1] = (x, y, m1)
        else:
            _fill_terminal(mol, refs, base, t, inward[t], pos, None)
    return pos, None, ()


def atropisomer_parity(mol, unit, z=None, wedge_of=None, log=None):
    """Parity of an atropisomer `unit` from the drawing: 0 unset, 1 even, 2 odd.

    A biaryl axis is chiral for the same reason an allene is -- four positions on two atoms, two pairs
    that cannot become coplanar -- so it is read the same way: the signed volume of the four in the
    core's ``refs`` order, which for this kind is each pivot's two RING directions in CSR ascending
    order, the anchor pivot's pair first.  No new sign convention is introduced, and none could be: the
    module's one rule (negative volume is parity 2) fixes it.

    WHAT THE WEDGE MEANS ON AN AXIS.  Both rings are drawn flat in the paper, which is the drawing's
    lie -- their planes are nearly perpendicular in fact.  A wedge on a pivot-to-ortho RING bond says
    that ortho comes toward the viewer, so that pivot's ring is the one standing out of the paper and
    its other ortho goes behind: exactly the allene terminal's situation, and handled by the same
    collapse onto the pivot's own position with ``z = +1`` and ``z = -1``.  The collapse discards a real
    along-axis component of those orthos, and cannot change the answer: it moves the determinant's apex
    along the axis line, while both of the far pivot's orthos lie on the same side along it.

    Returns 0, never raises, for every drawing that states nothing: no wedge on the axis -- the common
    case, and 4280 of the 4368 axes in a 119,534-molfile sample of a production corpus -- a wavy "either"
    bond, wedges on both pivots, both of one pivot's ring bonds wedged the same way, or a vanishing
    determinant.  There is no second channel to fall back on: **no CTfile field states an axial
    configuration**, in either version.

    `z` and `wedge_of` are :func:`tetrahedral_parity`'s, for the same purposes.
    """
    anchor = unit['anchor']
    if wedge_of is None:
        wedge_of = mol.wedge_of
    z = z or {}

    pivots = _atropisomer_pivots(mol, unit)
    if pivots is None:  # pragma: no cover - the core emits this kind only for a bond it perceived
        if log is not None:
            log.append(LogRecord('wedge:atropisomer-pivots', (anchor,),
                                 f'atom {anchor}: the axis partner pivot was not found, configuration left '
                                 f'unset', LOST))
        return 0
    pairs = ((pivots[0], 0), (pivots[1], 2))
    # No hole check and no missing-direction check, unlike the allene: `_atropisomer_end` admits a pivot
    # only at degree three with no hydrogen of any kind, so both of its ring directions are named atoms
    # and `unnamed_mask` is 0 for this kind.
    inward = {pivots[0]: pivots[1], pivots[1]: pivots[0]}
    pos, reason, named = _axial_frame(mol, unit, pairs, inward, z, wedge_of)
    if reason == 'either':
        if log is not None:
            t, r = named
            log.append(LogRecord('wedge:either-bond', (anchor, t, r),
                                 f'atom {anchor}: bond {t}-{r} drawn as either, configuration left unset', LOST))
        return 0
    if reason == 'flat':
        return 0  # a biaryl drawn without an axial statement, which is most of them
    if reason == 'both-ends':
        if log is not None:
            log.append(LogRecord('wedge:atropisomer-both-pivots', (anchor,),
                                 f'atom {anchor}: wedges on both pivots of the axis, configuration left '
                                 f'unset', LOST))
        return 0
    if reason == 'same-direction':
        if log is not None:
            log.append(LogRecord('wedge:atropisomer-same-direction', (anchor,),
                                 f'atom {anchor}: both ring bonds of one pivot wedged the same way, '
                                 f'configuration left unset', LOST))
        return 0

    volume = signed_volume(*pos)
    if volume == 0.0:
        if log is not None:
            log.append(LogRecord('wedge:degenerate-geometry', (anchor,),
                                 f'atom {anchor}: degenerate atropisomer layout, configuration left unset', LOST))
        return 0
    return 2 if volume < 0.0 else 1


def _atropisomer_pivots(mol, unit):
    """``(anchor, partner)`` -- the two pivots of an axis, or ``None``.

    ASKED OF THE CONSTITUTION, not taken as "the anchor's third neighbour": the partner is the
    neighbour that carries the far ref pair, which is the same question :func:`cis_trans_frame` asks
    about its partner terminal and the same reason -- a frame found by elimination is a frame that
    still looks right when the refs are not what this function assumed.
    """
    refs = unit['refs']
    anchor = unit['anchor']
    for m in mol.neighbors_of(anchor):
        neighbors = mol.neighbors_of(m)
        if refs[2] in neighbors and refs[3] in neighbors:
            return anchor, m
    return None


def _fill_terminal(mol, refs, base, terminal, inward, pos, z):
    """Positions for one terminal's ref pair, from the coordinates as drawn.

    The hole -- the direction with no atom -- goes opposite the sum of the terminal's other directions,
    the chain bond included.  Same construction as the tetrahedral hole, one direction shorter: an sp2
    terminal has three directions, not four.  `z` of None means a flat drawing, where this side of the
    frame is the one that lies in the paper.
    """
    xt, yt = mol.xy_of(terminal)
    zt = 0.0 if z is None else z.get(terminal, 0.0)
    named = [j for j in (base, base + 1) if refs[j] is not None]
    p = {j: (*mol.xy_of(refs[j]), 0.0 if z is None else z.get(refs[j], 0.0)) for j in named}
    xi, yi = mol.xy_of(inward)
    zi = 0.0 if z is None else z.get(inward, 0.0)
    for j in (base, base + 1):
        if refs[j] is not None:
            pos[j] = p[j]
        else:
            sx = (xi - xt) + sum(p[k][0] - xt for k in named)
            sy = (yi - yt) + sum(p[k][1] - yt for k in named)
            sz = (zi - zt) + sum(p[k][2] - zt for k in named)
            pos[j] = (xt - sx, yt - sy, zt - sz)


def _permutation_is_odd(keys):
    """Whether sorting `keys` into ascending order takes an odd number of transpositions.

    A parity is a statement about an *ordered* frame, so the same configuration reported against two
    different orders differs by exactly the parity of the permutation between them.  This is the whole
    of the arithmetic that translates the file's frame into the core's.  Bubble counting rather than
    cycle decomposition because a tetrahedron has four neighbours and clarity is worth more than the
    asymptotics.
    """
    n = len(keys)
    swaps = 0
    for i in range(n):
        for j in range(i + 1, n):
            if keys[i] > keys[j]:
                swaps += 1
    return bool(swaps & 1)


def stated_parity(unit, field, positions, log=None):
    """Core parity from the file's own atom parity field, or 0 when it cannot be translated.

    `field` is the raw value -- V2000 ``sss``, V3000 ``CFG=`` -- where 1 is odd and 2 is even; 0 means
    unstated and 3 means "either", and neither is a configuration.  `positions` maps stable id to
    0-based atom-block position, which is the frame the field is measured in.

    THE FRAME AND THE SIGN ARE BOTH MEASURED, not read off the specification.  The field's neighbour
    order is the atoms in **ascending atom-block position**, with the direction that has no atom of its
    own -- an implicit hydrogen, a lone pair -- ranking last, and ``1`` maps to core parity 1 directly
    once the permutation into the core's frame is applied.  The evidence is the repository corpus:
    on the 1053 stereocentres where both the field and the wedge drawing state a configuration, that
    reading agrees with the drawing on 1046.  The competing reading -- the specification's own
    "viewing the centre from behind the highest-numbered neighbour", taken as putting that neighbour
    first in the frame -- agrees on 7, necessarily, because a cyclic rotation of four elements is an
    odd permutation and the two readings are therefore exact opposites.  The corpus is choosing
    between them rather than tolerating both.

    ``0`` is returned, not raised, for a field that is not a configuration and for a centre with more
    than one undrawn direction, where the frame has two indistinguishable slots and the permutation
    into it is not defined.
    """
    if field not in (1, 2):
        return 0
    # The undrawn direction ranks above every real atom.  Any value past the largest position does,
    # and there is at most one such slot for the parity to be defined at all.
    high = len(positions) + 1
    keys = [high if r is None else positions[r] for r in unit['refs']]
    if keys.count(high) > 1:
        if log is not None:
            log.append(LogRecord('wedge:ambiguous-parity-frame', (unit["anchor"],),
                                 f'atom {unit["anchor"]}: stated parity {field} not read, {keys.count(high)} '
                                 f'undrawn directions leave the field\'s frame ambiguous', LOST))
        return 0
    return (3 - field) if _permutation_is_odd(keys) else field


def assign_parities(mol, z=None, log=None, stated=None, positions=None, configurations=None):
    """Read every stereo unit's configuration out of the drawing and store it.  Returns log lines.

    Only units the core calls stereogenic get a parity.  A wedge on a non-stereogenic atom is the
    file over-specifying -- common, and harmless -- and is reported rather than stored, because
    storing it would claim a configuration the constitution cannot distinguish.

    `stated` maps stable id to the file's atom parity field and `positions` to atom-block position;
    together they are the second, coordinate-free statement a CTfile can make about a tetrahedral
    centre.  **The drawing wins where both speak, and the disagreement is logged.**  Two reasons, and
    the second is why this is not simply deference to the specification's "ignored by readers":

    * the drawing is what a person looked at and approved, while the parity field is derived and the
      specification tells readers to ignore it precisely because writers get it wrong;
    * where the corpus's 1053 doubly-stated centres disagree, they disagree on records that are
      pathological *as drawings* -- a wedge hub with three bonds pointing at one atom, a centre
      carrying both an up and a down wedge -- so the seven contradictions are a signal that something
      is wrong with the record, not evidence that the field is the better source.

    Where the drawing is silent the field is used.  **On this repository's corpus that recovers
    nothing**: the fixed-point search below reaches the same centres from the drawing, which is the
    better source for them.  What the field is kept for is the case the corpus does not contain -- a
    record with no meaningful layout, where the field is the only statement there is -- and
    ``test_wedge.py`` therefore covers it with a hand-written fixture rather than a corpus sweep.
    chython 2 reads the field too and prefers it over the drawing; this is the one place the two
    deliberately differ, and it differs on 7 centres out of 1064.

    `configurations` is the same arrangement for **double** bonds, one kind further along: it maps the
    bond's two stable ids, low first, to ``(letter, refs)`` as :func:`stated_cis_trans` takes them --
    the non-geometric ``<bondStereo>C``/``T`` an XML dialect can state and no CTAB can.  It is a
    separate argument rather than a fifth key in `stated` because the two are keyed differently, an
    atom against a bond, and because a dialect can supply either without the other.  **The drawing wins
    here too**, and for the stronger of the two reasons above: the drawing is what a person looked at.
    """
    out = [] if log is None else log
    stated = stated or {}
    positions = positions or {}
    configurations = configurations or {}
    geometry = mol.has_coordinates or bool(z and any(z.values()))
    if not geometry:
        if any(mol.wedges()) and any(u['stereogenic'] for u in mol.stereo_units()):
            out.append(LogRecord('wedge:no-coords-read', (),
                                 'wedges present but no coordinates; stereo read from the atom parity field '
                                 'only'))
        if not stated and not configurations:
            return out
    else:
        # BEFORE the search, because it repairs the DRAWING and every read below asks the drawing.  A
        # wedge pointing the wrong way along its bond is otherwise not read at all, and the centre it
        # was drawn for comes out flat.
        _reanchor_wedges(mol, z, out)

    # ONE PASS IS NOT ENOUGH, and this loop is not an optimisation but a correctness requirement.
    # Whether an atom is stereogenic depends on whether its neighbours are *distinguishable*, and
    # configuring a neighbouring centre is what makes two otherwise identical branches differ.  So a
    # chain of centres resolves from the outside in: in `test/stereo.sdf` record 198 the middle atom
    # of a five-centre chain is not stereogenic until the four around it carry parities, and a
    # single-pass reader drops it -- with the actively misleading log line "wedge drawn on a
    # non-stereogenic centre" on an atom that is one.
    #
    # THE KINDS NEED NO ORDERING RELATIVE TO EACH OTHER, and it is worth saying why, because "read the
    # tetrahedral centres first, then the axial ones" is the obvious wrong answer.  The dependency runs
    # both ways: in `test/stereo.sdf` record 242 an allene is stereogenic only once the two tetrahedral
    # centres on one of its terminals carry parities, and nothing forbids the mirror case.  What makes
    # a single loop sufficient is that it re-reads `stereo_units()` on every pass, so whichever kind
    # becomes readable first is read first, by the constitution rather than by a hardcoded order.
    #
    # Logging is suppressed during the search and the reasons are collected in a final pass, because
    # an atom that stays undetermined would otherwise be reported once per iteration.
    # At worst one centre resolves per pass, so the number of stereo units bounds the search.
    for _ in range(sum(1 for _ in mol.stereo_units()) + 1):
        parities = {}
        # Collected, then written in one scope.  Every geometry read needs a clean arena, so a
        # `set_parity` in the middle of the loop would have to apply its journal before the next
        # read -- which the core does do, silently, at one buffer copy per centre.
        for unit in mol.stereo_units():
            anchor = unit['anchor']
            if not unit['stereogenic'] or mol.parity_of(anchor):
                continue
            parity = _unit_parity(mol, unit, z, stated, positions, configurations, geometry, None)
            if parity:
                parities[anchor] = parity
        if not parities:
            break
        with mol.edit():
            for anchor, parity in parities.items():
                mol.set_parity(anchor, parity)
    else:  # pragma: no cover - one centre resolves per pass at worst, so the bound cannot be reached
        out.append(LogRecord('wedge:unsettled', (),
                             'stereo resolution did not settle; some configurations may be unread', LOST))

    # The final pass exists to explain what is left, against the settled constitution.  Anything it
    # says about an atom is true of the molecule the caller gets, which is not something a line
    # emitted mid-search can promise.
    wedged = {narrow for narrow, _, _ in mol.wedges()}
    unsupported = set()
    for unit in mol.stereo_units():
        anchor = unit['anchor']
        kind = unit['kind']
        if not unit['stereogenic']:
            if anchor in wedged:
                out.append(LogRecord('wedge:non-stereogenic', (anchor,),
                                     f'atom {anchor}: wedge drawn on a non-stereogenic centre, ignored', LOST))
            continue
        if kind not in (SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER):  # pragma: no cover
            # THE GUARD FOR A KIND THIS MODULE DOES NOT READ.  No such kind exists -- the core emits
            # four and all four are read here -- so this is unreachable today and stays as the thing a
            # fifth kind meets instead of being silently measured in another kind's frame.
            unsupported.add(kind)
        elif not mol.parity_of(anchor):
            _unit_parity(mol, unit, z, stated, positions, configurations, geometry, out)
        elif kind == SU_CIS_TRANS and configurations:
            # The same contradiction check the tetrahedral branch below runs, for the other channel a
            # document can state a configuration in.  Reached only when the drawing already settled the
            # bond, so agreement is silent and disagreement is the whole content of the line.
            spelled = _stated_for_bond(mol, unit, configurations, out)
            drawn = mol.parity_of(anchor)
            if spelled and spelled != drawn:
                out.append(LogRecord('wedge:drawing-field-disagree', (anchor,),
                                     f'bond stereo at atom {anchor}: the drawing and the stated configuration '
                                     f'disagree (drawn {"C" if drawn == 2 else "T"}, the document says '
                                     f'{"C" if spelled == 2 else "T"}); keeping the drawing', REFUSED))
        elif kind == SU_TETRA and anchor in stated:
            said = stated_parity(unit, stated[anchor], positions, out)
            drawn = mol.parity_of(anchor)
            if said and said != drawn:
                # A record that contradicts itself.  Logged loudly and per atom, because the caller
                # cannot see it any other way and it usually means the drawing is damaged.
                out.append(LogRecord('wedge:drawing-field-disagree', (anchor,),
                                     f'atom {anchor}: the drawing and the atom parity field disagree (drawn '
                                     f'{drawn}, field says {said}); keeping the drawing. Check the wedges on '
                                     f'this centre -- two wedges on one centre, or a wedge whose narrow end is '
                                     f'elsewhere, produce exactly this', REFUSED))
    for kind in sorted(unsupported):  # pragma: no cover - see the guard above
        out.append(LogRecord('wedge:unsupported-kind', (),
                             f'kind {kind} stereo not read from the drawing', LOST))
    return out


def _reanchor_wedges(mol, z, out):
    """Move a wedge whose narrow end cannot hold a configuration onto the end that can.  A REPAIR.

    A CTfile wedge points at the atom it is a statement about: its narrow end is the stereocentre, and
    the reader looks only at the half-edges leaving the centre.  A drawing that puts the point at the
    other end -- ``11 13  1  6`` where 11 is an amine nitrogen with nothing to configure and 13 is the
    tetrahedral centre -- therefore states a configuration the reader cannot see, and the centre comes
    out flat.

    The bond and the direction are the file's; only which end the triangle's point sits on moves, so the
    signed volume the reader computes is the one the drawing shows.  Under the input posture that is a
    repair with a log line, not a loss.

    Narrow, deliberately: the destination is a stereogenic tetrahedral centre the bond already reaches,
    the source end holds no stereo unit of its own -- an allene TERMINAL holds one, since
    :func:`allene_parity` reads the terminals' wedges rather than the anchor's -- the centre carries no
    wedge of its own to contradict, and the moved wedge has to read back as a configuration. Anything
    looser starts inventing stereochemistry out of a decorative triangle.
    """
    if z and any(z.values()):     # a 3D drawing states the configuration outright; wedges are ignored
        return
    holders = set()
    tetra = {}
    for unit in mol.stereo_units():
        if not unit['stereogenic']:
            continue
        if unit['kind'] == SU_ALLENE:
            holders.update(t for t, _ in _chain_ends(mol, unit['anchor']))
        elif unit['kind'] == SU_ATROPISOMER:
            # BOTH pivots, for the allene's reason: `atropisomer_parity` reads the wedges of either
            # end, so a wedge at the far pivot is a statement this module now understands and must not
            # be re-anchored onto a neighbouring tetrahedral centre.
            holders.update(_atropisomer_pivots(mol, unit) or (unit['anchor'],))
        else:
            holders.add(unit['anchor'])
            if unit['kind'] == SU_TETRA:
                tetra[unit['anchor']] = unit
    wedged = {narrow for narrow, _, _ in mol.wedges()}
    moves = []
    for narrow, wide, code in mol.wedges():
        if narrow in holders or wide in wedged or wide not in tetra:
            continue
        unit = tetra[wide]
        if narrow not in unit['refs']:
            continue

        def probe(a, b, _w=wide, _n=narrow, _c=code):
            return _c if (a, b) == (_w, _n) else mol.wedge_of(a, b)

        parity = tetrahedral_parity(mol, unit, z, probe)
        if parity:
            moves.append((narrow, wide, code, parity))
    if not moves:
        return
    with mol.edit():
        for narrow, wide, code, _ in moves:
            mol.set_wedge(narrow, wide, WEDGE_NONE)
            mol.set_wedge(wide, narrow, code)
    for narrow, wide, code, parity in moves:
        out.append(LogRecord('wedge:reanchored', (wide, narrow),
                             f'atom {wide}: the wedge on bond {narrow}-{wide} had its narrow end at {narrow}, '
                             f'which holds no configuration; re-anchored at {wide}, which reads parity '
                             f'{parity}', REPAIRED))


def _stated_for_bond(mol, unit, configurations, log):
    """A cis/trans unit's stated configuration, looked up by the bond rather than by the anchor.

    ONE LOOKUP FOR TWO CALLERS -- the fixed-point search and the contradiction report -- for the same
    reason :func:`_unit_parity` exists: a second copy of "which bond is this unit about" is a second
    chance for the log to describe a bond other than the one that was read.  A document keys the letter
    to a bond and the core keys a parity to an anchor, and the frame is what translates between them.
    """
    frame = cis_trans_frame(mol, unit)
    if frame is None:
        return 0
    _, anchor, partner, _ = frame
    said = configurations.get((anchor, partner) if anchor < partner else (partner, anchor))
    if said is None:
        return 0
    return stated_cis_trans(mol, unit, said[0], said[1], log)


def _unit_parity(mol, unit, z, stated, positions, configurations, geometry, log):
    """One unit's parity from every source, drawing first.  0 when nothing determines it.

    Split out of :func:`assign_parities` because it is called from two places that must not disagree:
    the fixed-point search, with logging off, and the final reporting pass, with logging on.  Two
    copies of this decision would be two chances for the log to describe something other than what
    was stored.
    """
    anchor = unit['anchor']
    if unit['kind'] == SU_TETRA:
        parity = tetrahedral_parity(mol, unit, z, log=log) if geometry else 0
        if not parity and anchor in stated:
            parity = stated_parity(unit, stated[anchor], positions, log)
        return parity
    if unit['kind'] == SU_CIS_TRANS:
        # No parity field exists for a double bond in either MDL version -- bond stereo 3 says only
        # "cis or trans, unknown which" -- so for a molfile the coordinates are the only source there
        # is, and `configurations` arrives empty.  An XML dialect states the letter as well, and it is
        # consulted second for the same reason the atom parity field is: the drawing wins.
        parity = cis_trans_parity(mol, unit, log=log) if geometry else 0
        if not parity and configurations:
            parity = _stated_for_bond(mol, unit, configurations, log)
        return parity
    if unit['kind'] == SU_ALLENE:
        # Nor for an axial configuration, in either version: the drawing is the only statement.
        return allene_parity(mol, unit, z, log=log) if geometry else 0
    if unit['kind'] == SU_ATROPISOMER:
        return atropisomer_parity(mol, unit, z, log=log) if geometry else 0
    return 0  # pragma: no cover - every kind the core emits has a branch above


class _Planar:
    """`mol` with :meth:`xy_of` answering out of a supplied layout instead of the stored coordinates.

    A PROXY RATHER THAN A ``plane`` ARGUMENT THREADED THROUGH THE READ PATH, deliberately.  Every
    geometric question the chooser asks is asked through :func:`tetrahedral_parity`, which is the
    reader's own function -- that is what makes a round trip incapable of inverting a sign, and it is
    worth more than the indirection costs.  Giving that function a second way to obtain coordinates
    would put a plane parameter on the read path too, where no reader ever wants one, so the
    substitution happens at the object instead.  Everything other than the coordinates is forwarded,
    so the proxy cannot drift from the molecule it stands for.
    """
    __slots__ = ('_mol', '_plane')

    def __init__(self, mol, plane):
        self._mol = mol
        self._plane = plane

    def xy_of(self, n):
        return self._plane[n]

    @property
    def has_coordinates(self):
        return True

    def __getattr__(self, name):
        return getattr(self._mol, name)


def wedges_for_write(mol, log=None, plane=None, cis_trans_stated=()):
    """Wedges to emit, as ``[(narrow, wide, wedge), ...]``.  Returns ``(wedges, log)``.

    `plane` is ``{stable id: (x, y)}`` -- the layout to draw for, defaulting to the molecule's own
    coordinates.  It is not a convenience: **which bond carries a wedge and which way it points is a
    property of the drawing, not of the molecule**, so a renderer that has just computed a temporary
    layout must ask for the wedges of *that* layout.  Asking for the stored ones instead states the
    configuration of a picture nobody is drawing, and for a molecule with no stored coordinates at
    all it states nothing and drops the stereochemistry silently.

    A molecule read from a CTfile already carries the wedges it was drawn with, and those are
    returned untouched -- re-deriving them would silently redraw a file the user asked to round-trip.
    **That shortcut is only taken when `plane` is not given**, because stored wedges are only valid
    for the coordinates they were read with: reflect the layout and the same codes on the same bonds
    state the opposite configuration.  A caller who supplies a layout is asking about that layout.

    Otherwise a wedge is *chosen* per configured unit, by trying up and down on a candidate bond and
    keeping whichever one the matching reader -- :func:`tetrahedral_parity`, :func:`allene_parity` or
    :func:`atropisomer_parity` -- reads back as the stored parity.  The sign is therefore never computed
    twice in two places, which is the only way to make a round trip structurally incapable of inverting
    it.

    Every kind this module reads is written, and an axial kind needs writing for the same reason it
    needs reading: a molecule built in code or parsed from a pach record carries an axial parity and no
    wedges, and emitting it flat would drop the descriptor silently.  Where the wedge goes differs by
    kind -- a tetrahedral centre wedges one of *its own* bonds, an allene one of a *terminal's*, an
    atropisomer one of a *pivot's ring* bonds -- and the anchor of either axial kind has no bond of its
    own to wedge at all.

    A cis/trans unit needs nothing here: its configuration is in the coordinates, which are written
    whatever this returns.

    **A WRITER WITH NO LAYOUT LOSES EVERY DOUBLE-BOND CONFIGURATION THE MOLECULE STATES, and this
    function is the one place that says so.**  A molecule from SMILES or built in code holds the fact
    perfectly well (``C/C=C/C`` carries ``kind: 1, parity: 1`` with no coordinates anywhere), and a CTAB
    has nowhere to put it: the coordinates are its only channel.  The loss is the format's; the silence
    would have been ours.

    `cis_trans_stated` is the anchors the *caller* is about to state some other way, and those are not
    reported.  CML and MRV have :func:`stated_cis_trans`'s non-geometric channel and use it, so the line
    would be false for them -- per anchor rather than per format, because MRV's spelling names no
    reference atoms and is therefore only writable where it cannot be ambiguous.  A caller with a channel
    passes what it wrote; the CTAB versions pass nothing and keep the line.
    """
    out = [] if log is None else log
    # THE FORMAT CANNOT HOLD THIS AND SAYS SO, reported before every branch below so that no shortcut
    # can skip it: a molecule whose only stated stereo is a double-bond configuration produces no
    # `configured` units at all and would otherwise leave at the early return three lines down, in
    # silence.
    #
    # ONE EMITTER FOR FIVE CALLERS, and that is the whole point of it living here.  A copy of this loop
    # per format -- same condition, same wording, differently shaped -- leaves the two CTAB versions
    # silent while the two XML dialects speak, and gives the next author a fourth chance to word it
    # differently.  The fact is molecule-side ("a configuration is stated and there is no drawing to put
    # it in") and identical for every coordinate-carrying format, so it belongs beside the chooser every
    # one of them calls, not replicated in each.  Nothing here is MDL-specific, and the message must not say MDL.
    #
    # A SEPARATE STATEMENT RATHER THAN A WIDER `configured`, and that is not a stylistic choice.  That
    # list drives wedge *selection*, and a cis/trans unit wants no wedge -- putting `SU_CIS_TRANS` in it
    # would hand the allocator a unit it has no bond to serve and make one noun count two things.  What
    # the two cases share is the writer's knowledge that it has no layout, which is here.
    #
    # `plane is not None` is a caller drawing a layout it computed, and a laid-out molecule expresses
    # its double bonds in the coordinates like any other, so nothing is lost and nothing is said.  That
    # is also why `chython.depict` never sees this line: it always supplies the layout it drew.
    #
    # Three things are deliberately NOT done.  It does not raise -- a writer that refuses produces no
    # file at all, which is worse for every consumer than a file missing one descriptor a layout would
    # regenerate.  It does not write MDL bond stereo 3, "cis or trans, unknown which": we know which, so
    # that flag is a false statement and a reader taking it back would gain a WRONG fact where it
    # currently gains none.  And it does not generate coordinates: inventing a drawing to preserve a
    # descriptor invents the descriptor's evidence, and no read or write path in this tree calls a layout
    # engine.
    #
    # Writing CML's and MRV's `<bondStereo>C`/`T` is not on that list any more, and it never belonged
    # HERE: a wedge chooser returns wedges, and a non-geometric descriptor is a field on the record.  It
    # is `CtabBond.configuration`, the two XML writers fill it from `cis_trans_letter`, and what reaches
    # this function is the resulting `cis_trans_stated` -- the anchors that already have a home, which
    # this line must not claim were dropped.
    #
    # One line per unit, naming the anchor, because that is the atom the caller has to look at; an
    # unconfigured stereogenic double bond is not a loss and is not counted -- nothing was stated, so
    # nothing is dropped, and counting it would fire the line on most molecules in any corpus.
    if plane is None and not mol.has_coordinates:
        for unit in mol.stereo_units():
            if unit['kind'] == SU_CIS_TRANS and mol.parity_of(unit['anchor']) \
                    and unit['anchor'] not in cis_trans_stated:
                out.append(LogRecord('wedge:double-bond-no-coords', (unit["anchor"],),
                                     f'unsupported: stereo: atom {unit["anchor"]}: a double-bond configuration '
                                     f'on a molecule with no coordinates is not written', LOST))
    if plane is None:
        existing = list(mol.wedges())
        if existing:
            return existing, out
    configured = [u for u in mol.stereo_units()
                  if u['kind'] in (SU_TETRA, SU_ALLENE, SU_ATROPISOMER) and mol.parity_of(u['anchor'])]
    if not configured:
        return [], out
    if plane is not None:
        mol = _Planar(mol, plane)
    elif not mol.has_coordinates:
        out.append(LogRecord('wedge:no-coords-write', (),
                             f'{len(configured)} configured stereocentre(s) but no coordinates; '
                             f'no wedges written', LOST))
        return [], out

    # THE AXIAL KINDS FIRST, and each one reserves every bond its own reader looks at.  An axial unit
    # reads the wedges of BOTH ends, so a wedge chosen afterwards for some other centre on one of those
    # bonds can change what this one reads back -- a tetrahedral centre cannot be disturbed that way,
    # since it only ever reads bonds incident to itself.  Reserving makes the choice order-independent
    # in the one direction where order would otherwise decide the answer.
    #
    # It is also why the tetrahedral prioritisation below runs SECOND and over `taken`: an axial unit has
    # very few bonds it can possibly use, while a tetrahedral centre usually has three, so letting the
    # prettier-drawing heuristic bid against an allene's only option would trade a configuration for an
    # aesthetic.  Correctness is the gate; the heuristic works on what is left.
    axial_wedges = {}
    taken = set()
    for unit in sorted((u for u in configured if u['kind'] in (SU_ALLENE, SU_ATROPISOMER)),
                       key=lambda u: (u['kind'], u['anchor'])):
        anchor = unit['anchor']
        target = mol.parity_of(anchor)
        if unit['kind'] == SU_ALLENE:
            reader = allene_parity
            kind_name = 'allene'
            # Every single bond from either terminal to a named direction.  `order_of == 1` is not
            # decoration: the chain bonds are double, and MDL has no wedge for a double bond.
            candidates = [(t, r) for t, _ in _chain_ends(mol, anchor) for r in unit['refs']
                          if r is not None and r in mol.neighbors_of(t) and mol.order_of(t, r) == 1]
        else:
            reader = atropisomer_parity
            kind_name = 'atropisomer'
            # A pivot's two RING bonds, and no order filter: an axial statement about a biaryl is drawn
            # on the ring bond by every package that draws one, and the ring bond is aromatic or double
            # as often as it is single.  Both refs of each pair are named atoms for this kind.
            pivots = _atropisomer_pivots(mol, unit)
            candidates = [] if pivots is None else \
                [(t, r) for t, base in zip(pivots, (0, 2)) for r in unit['refs'][base:base + 2]]
        reserved = list(candidates)
        # An explicit hydrogen first because wedging it is what every drawing package does, then a
        # terminal atom, then lowest id for determinism.  Left as it was rather than routed through
        # `_draw_cost`: that cost is built out of comparisons between bonds at one tetrahedral centre,
        # and an axial unit's candidates come from two different atoms, so its ring and degree terms
        # would be comparing things the reasoning behind them does not cover.
        candidates.sort(key=lambda ab: (mol.element_of(ab[1]) != 1, mol.degree_of(ab[1]) != 1, ab))
        placed = False
        for narrow, r in candidates:
            if (narrow, r) in taken or (r, narrow) in taken:
                continue
            for wedge in (WEDGE_UP, WEDGE_DOWN):
                trial = dict(axial_wedges)
                trial[(narrow, r)] = wedge

                def probe(a, b, _t=trial):
                    return _t.get((a, b), WEDGE_NONE)

                if reader(mol, unit, None, probe) == target:
                    axial_wedges[(narrow, r)] = wedge
                    taken.add((narrow, r))
                    taken.update(reserved)
                    placed = True
                    break
            if placed:
                break
        if not placed:
            out.append(LogRecord('wedge:no-writable-bond', (anchor,),
                                 f'atom {anchor}: no bond can carry a wedge that reproduces its '
                                 f'configuration; {kind_name} written flat', LOST))

    # Every configured anchor, both kinds, for cost term 3 -- a wedge should not point at any of them.
    centres = {u['anchor'] for u in configured}
    options = {}
    for unit in configured:
        if unit['kind'] != SU_TETRA:
            continue
        anchor = unit['anchor']
        feasible = [(r, code) for r, code in _feasible_wedges(mol, unit, mol.parity_of(anchor))
                    if (anchor, r) not in taken and (r, anchor) not in taken]
        if not feasible:
            out.append(LogRecord('wedge:no-writable-bond', (anchor,),
                                 f'atom {anchor}: no bond can carry a wedge that reproduces its '
                                 f'configuration; centre written flat', LOST))
            continue
        options[anchor] = sorted((_draw_cost(mol, anchor, r, centres), r, code)
                                 for r, code in feasible)

    wedges = [(narrow, wide, code) for (narrow, wide), code in axial_wedges.items()]
    if options:
        # An axial unit's wedge is as crowding as any other, so its endpoints count as drawn for the
        # adjacency preference -- but through `busy`, which is advice, not through `taken`, which is a
        # veto.  A tetrahedral centre may end up sharing an atom with an axial wedge; it may not steal
        # the bond.
        wedges += _assign(mol, options, out, taken,
                          {a for pair in axial_wedges for a in pair})
    # Sorted, so the emitted order is the molecule's own and not the chooser's cost ranking.
    return sorted(wedges), out


def wedge_in_file_order(wedge_of, n, m):
    """`(a, b, code)` -- the bond's endpoints in the order CTfile wants, and its wedge code or None.

    In CTfile the point of a wedge is at the FIRST atom of the bond line, and that atom is the one the
    wedge is a statement about, so a wedge whose narrow end is `m` reverses the bond as written.  Both
    writers had this, and it does not come from `Bond.wedge` or `mol.wedge_between`: `wedges_for_write`
    returns wedges it may have chosen for a molecule that carries none, so the arena is the wrong
    place to ask.
    """
    code = wedge_of.get((n, m))
    if code is not None:
        return n, m, code
    code = wedge_of.get((m, n))
    if code is not None:
        return m, n, code
    return n, m, None


def _feasible_wedges(mol, unit, target):
    """``[(ref, code), ...]`` -- every bond out of this centre that can state `target`, unordered.

    WHY THIS CAN BE ANSWERED ONE CENTRE AT A TIME, which is what makes choosing tractable: a wedge
    lives on a single half-edge, the one leaving its narrow end, and :func:`tetrahedral_parity` at an
    anchor reads only the half-edges leaving that anchor.  So no wedge chosen for one TETRAHEDRAL centre
    can change what another reads back, and feasibility never has to be re-tested once the search starts
    trading bonds between centres.  A CTfile wedge always has its narrow end at the stereocentre, so
    that is not an accident of the storage but the convention itself.

    **This is a property of tetrahedral units only, and the difference is a trap.**
    :func:`allene_parity` reads the wedges of both of its terminals -- neither of which is its anchor --
    so an allene's reading *is* disturbed by a wedge chosen elsewhere.  That is why the caller settles
    every allene first and hands the bonds they depend on to :func:`_assign` as a veto, and why this
    function is only ever asked about tetrahedral units.

    At most one of up and down can be right -- the determinant is linear in the out-of-plane
    displacement, so flipping the code flips the sign -- and neither is right when the centre's
    geometry is degenerate along that bond, which is why a bond can fail to be a candidate at all.
    A multiple bond is excluded outright, for a reason that is about the file and not the geometry;
    see below.
    """
    anchor = unit['anchor']
    out = []
    for r in unit['refs']:
        # `None` is a direction with no atom of its own: an implicit hydrogen, or a lone pair. It is a
        # perfectly good frame direction and a hopeless wedge, because a CTfile wedge is a line in the
        # bond block and there is no second atom to name. This is the whole reason a ring-fusion carbon
        # has nothing but ring bonds to offer.
        if r is None:
            continue
        # THE BOND STEREO FIELD IS OVERLOADED BY BOND ORDER, so a multiple bond has no wedge to give.
        # On a single bond `sss` is the wedge -- 1 up, 6 down, 4 either -- and on a double bond it is
        # cis/trans instead: 0 "use the coordinates", 3 "either". Writing 1 there does not state a
        # wedge to any conforming reader; it puts an out-of-domain value in the double bond's own
        # stereo field and loses the configuration it was meant to state. chython's reader accepts it
        # coming back, which is exactly what would hide the defect in a round-trip test.
        #
        # The allene half of `wedges_for_write` has always known this and filters its candidates on
        # `order_of == 1`. This half needs it for the same reason and did not have it: a sulfoxide,
        # sulfimide or phosphine oxide is a tetrahedral centre with a double bond among its refs, and
        # that bond is a *terminal, acyclic* one -- so it wins on every term of `_draw_cost` and is
        # exactly the bond the chooser reaches for first.
        if mol.order_of(anchor, r) != 1:
            continue
        for code in (WEDGE_UP, WEDGE_DOWN):
            def probe(a, b, _r=r, _c=code):
                return _c if (a, b) == (anchor, _r) else WEDGE_NONE

            if tetrahedral_parity(mol, unit, None, probe) == target:
                out.append((r, code))
                break
    return out


#: ``cos(30 degrees) ** 2``.  Two bonds at one atom separated by less than 30 degrees are treated as
#: too close for a wedge to be attributed to one of them rather than the other.
#:
#: THIRTY IS READ OFF THE CORPUS AND NOT TUNED.  Over the 877 bond directions at the configured
#: centres of ``test/wedge_stereo.sdf`` the nearest-sibling separation is bimodal, with a mode at
#: 110-120 degrees, 4 directions below 30 and **nothing at all between 20 and 30** -- so every cut in
#: that empty band classifies the corpus identically, and the measured effect of the term is the same
#: at 20 as at 30.  It is a gap in the data rather than a knob.  The lower bound is geometric: the
#: drawn triangle's own half-angle is ``atan(wedge_space / bond length)``, about 4.6 degrees at
#: ``depict``'s default width, so a sibling inside 30 degrees is well inside the region where the two
#: marks are read together.
_COLLINEAR_COS2 = 0.75


def _near_collinear(mol, anchor, r):
    """Is another bond at `anchor` drawn within 30 degrees of the bond to `r`?

    Compared as a squared cosine so the test is exact arithmetic on the stored coordinates: for the
    angle to be under the threshold the dot product must be positive *and* its square must exceed
    ``cos(30) ** 2`` times the two squared lengths.  Both guards are needed -- squaring loses the sign,
    and without the first test a bond at 150 degrees would read as one at 30.
    """
    ax, ay = mol.xy_of(anchor)
    ux, uy = mol.xy_of(r)
    ux -= ax
    uy -= ay
    u2 = ux * ux + uy * uy
    for other in mol.neighbors_of(anchor):
        if other == r:
            continue
        vx, vy = mol.xy_of(other)
        vx -= ax
        vy -= ay
        dot = ux * vx + uy * vy
        if dot > 0.0 and dot * dot > _COLLINEAR_COS2 * u2 * (vx * vx + vy * vy):
            return True
    return False


def _draw_cost(mol, anchor, r, centres):
    """How badly a wedge from `anchor` to `r` reads.  Lower is better; a total order, so deterministic.

    Derived from what the drawing asks a viewer to do, in decreasing order of how badly it misleads
    them.  Every term is a comparison between two bonds at the SAME centre, so nothing here is a
    judgement about the molecule.

    1. **A RING BOND LAST.**  The bonds of a drawn ring are what a viewer reads as the ring's plane, so
       a wedge on one of them asks for a ring atom to be lifted out of a plane the same picture asserts
       is flat.  The two readings are both available and the drawing does not say which is meant.  On a
       fused or bridged bond, which belongs to two rings, there are two planes to contradict.  An
       acyclic bond has no such second reading, so it wins outright -- this term dominates the rest.

    2. **NOT ON TOP OF A SIBLING BOND.**  A wedge is attributed to a bond by lying along it, so a
       second bond leaving the same atom within 30 degrees leaves the reader unable to say which of the
       two the triangle belongs to -- and the two state opposite things about the centre.  This is the
       one term that is about the *mark* rather than about what lies past it, which is why it ranks
       above the terms that are: a wedge whose bond cannot be identified says nothing at all, while a
       wedge with a branch point at its wide end says something imprecise.  It ranks below the ring test
       because that one is measured to dominate and because a near-collinear pair is rare -- 1 centre in
       282 on the corpus -- so promoting it above ring would trade 21 ring wedges for one.

    3. **A LEAF WIDE END NEXT.**  The wedge widens away from the centre and claims that the neighbour
       at the wide end is toward the viewer.  If that neighbour is a leaf -- hydroxyl oxygen, halogen,
       methyl -- nothing is drawn past it and the claim ends there.  If it is a branch point, every
       atom drawn beyond it is implicitly lifted too, and the viewer has to guess how far the claim
       reaches.

    4. **NOT AT ANOTHER STEREOCENTRE.**  A wedge whose wide end is itself a configured centre reads as
       if it also said something about that centre, which it does not: a reader keys on the narrow end.
       It is also the direct cause of the adjacent pairs :func:`_assign` then has to untangle, so
       discouraging it here is cheaper than repairing it there.  This term only ever discriminates
       among branch points, since a stereocentre has three heavy neighbours at least and so can never
       be a leaf -- it does not compete with term 3, it refines it.

    5. **THE HYDROGEN, BY CONVENTION.**  A wedge to the hydrogen is how a stereocentre is drawn
       everywhere, and it is the prior behaviour of this function.  It ranks here rather than first
       because it cannot lose by ranking here: an explicit hydrogen is a leaf and no hydrogen is in a
       ring, so it is already in the best class of every term above and this only breaks a tie inside
       it.

    6. **THE LONGER BOND.**  Pure geometry, and a tie-break rather than a preference: the wedge is a
       triangle whose base sits at the wide end, so on a short bond it is stubby and its base crowds
       whatever line art meets the neighbour.  A longer bond gives it room.  Read from the layout being
       drawn, which is the only place the answer exists -- and the reason this is a tie-break at all is
       that a normalised layout makes most bonds the same length, so it decides only where the terms
       above genuinely cannot.

    7. **THE LOWEST STABLE ID**, so that the order is total and two runs cannot differ.
    """
    ax, ay = mol.xy_of(anchor)
    rx, ry = mol.xy_of(r)
    # Term 1 stays first because `_assign` tiers on `cost[0]` -- adding anything ahead of the ring test
    # would silently retier the pass on whatever was put there.
    return (mol.bond_in_ring(anchor, r),
            _near_collinear(mol, anchor, r),
            mol.degree_of(r) != 1,
            r in centres,
            mol.element_of(r) != 1,
            -((rx - ax) ** 2 + (ry - ay) ** 2),
            r)


def _assign(mol, options, out, reserved, busy):
    """One wedge per tetrahedral centre out of `options`, as cheap as a greedy pass can make it.

    `reserved` is bonds an allene already owns and `busy` the atoms its wedges touch.  The asymmetry
    between them is the point: `reserved` is a veto, because taking one of those bonds would change what
    the allene reads back, and `busy` is only advice, because sharing an atom with an allene's wedge is
    merely crowded.

    TWO CONSTRAINTS AND ONE PREFERENCE, and which of them is allowed to yield is the whole design:

    * **a bond carries at most one wedge, ever.**  Never relaxed.  A bond drawn as a wedge from both
      ends is not a picture of anything, and since two centres bonded to each other are the only way to
      want it, refusing costs almost nothing.
    * **no two wedges share an atom.**  A shared atom is lifted by one wedge and is the base of
      another, and a viewer cannot hold both claims at once; a run of them along a chain is the defect
      that makes a drawing unreadable rather than merely ugly.  Relaxed for a centre that would
      otherwise have no wedge at all, because losing a configuration is a worse outcome than an ugly
      one -- correctness is the gate and beauty is the goal.
    * **the cost.**  Advice throughout.

    CHEAPEST-FIRST ACROSS ALL CENTRES, not centre by centre.  A centre holding a clean terminal bond
    should take it before a centre with nothing but ring bonds is asked to choose; walking centres in
    stable-id order does the opposite, and lets whichever comes first take the bond a later centre
    needed more.  That is a greedy pass over a global order and not an optimum -- the exact problem is
    a minimum-cost matching -- but the fallback means the difference is measured in ugliness rather
    than in lost stereochemistry, and the corpus metric in ``test_wedge.py`` is what says whether it is
    worth more machinery.

    RING-OR-NOT IS A TIER AND NOT MERELY THE FIRST TERM OF THE COST, which is the one subtlety here.
    The cost is lexicographic, so within a single pass the ring term already dominates the others --
    but it has to dominate the adjacency constraint too, and that constraint lives in the pass and not
    in the cost.  Measured, on the corpus of ``test_wedge.py``: relaxing adjacency only after both
    tiers had been tried left six centres on a ring bond that had an acyclic bond going spare, because
    their acyclic bond happened to touch a wedge already drawn.  Exhausting the acyclic tier first --
    adjacency and all -- removes all six and takes the corpus's shared-atom count from 2 to 8, which is
    the trade this makes explicitly: a ring wedge is *ambiguous*, while a shared atom is merely
    crowded, and an ambiguous drawing is the worse failure.  83 ring wedges is then exactly the number
    of centres in that corpus with no acyclic bond to offer, so nothing is left on the table.
    """
    chosen = {}
    used_bonds = set(reserved)
    used_atoms = set(busy)
    crowded = []

    def free_bond(anchor, r):
        return (anchor, r) not in used_bonds and (r, anchor) not in used_bonds

    def take(anchor, r, code):
        chosen[anchor] = (r, code)
        used_bonds.add((anchor, r))
        used_atoms.update((anchor, r))

    # Cost terms 1 and 2 -- ring, then near-collinear -- and both are tiers for the same reason: each
    # has to outrank the adjacency constraint, which lives in this pass and not in the cost.  Their
    # order between themselves is the cost's own, so an acyclic collinear bond is still preferred to a
    # clean ring bond; what the tiering adds is that a whole class is exhausted, sharing an atom
    # included, before the next is asked.
    for tier in ((False, False), (False, True), (True, False), (True, True)):
        tiered = {anchor: [o for o in opts if o[0][:2] == tier] for anchor, opts in options.items()}
        for _, anchor, r, code in sorted((cost, a, r, code) for a, opts in tiered.items()
                                         for cost, r, code in opts):
            if anchor in chosen or anchor in used_atoms or r in used_atoms:
                continue
            if free_bond(anchor, r):
                take(anchor, r, code)
        # Whoever this tier could not place cleanly shares an atom rather than dropping to the next
        # tier. Still preferring a wide end nothing else has claimed, since one shared atom reads
        # better than two.
        for anchor in sorted(tiered):
            if anchor in chosen or not tiered[anchor]:
                continue
            for _, r, code in sorted(tiered[anchor], key=lambda o: (o[1] in used_atoms, o[0])):
                if free_bond(anchor, r):
                    take(anchor, r, code)
                    crowded.append(anchor)
                    break

    for anchor in sorted(options):
        if anchor not in chosen:
            out.append(LogRecord('wedge:no-free-bond', (anchor,),
                                 f'atom {anchor}: every bond that could state its configuration is already '
                                 f'drawn as a wedge; centre written flat', LOST))

    # Reported rather than silent, and as one line each rather than one per atom: these are the two
    # compromises the drawing contains, a reader of the log can do nothing about either, and a steroid
    # would otherwise contribute a dozen lines saying the same thing.
    ring = sorted(a for a, (r, _) in chosen.items() if mol.bond_in_ring(a, r))
    if ring:
        out.append(LogRecord('wedge:ring-bond-wedge', tuple(ring),
                             f'{len(ring)} stereocentre(s) could only be drawn with a wedge on a ring bond '
                             f'(atoms {", ".join(map(str, ring))}); every other bond at them is a ring bond too'))
    if crowded:
        out.append(LogRecord('wedge:crowded', tuple(crowded),
                             f'{len(crowded)} stereocentre(s) share an atom with another wedge '
                             f'(atoms {", ".join(map(str, crowded))}); no independent bond was left for them'))

    # `out` is appended to in place, so only the wedges come back -- the caller has allene wedges of its
    # own to merge with these and is the one place that decides the emitted order.
    return [(a, r, code) for a, (r, code) in chosen.items()]
