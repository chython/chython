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
# MoleculeContainer: the Python-facing molecule, and the journal it is built through.
#
# Every mutation is appended to the journal as a 16-byte `journal_t` record and applied to the
# arena in one pass (`_apply`), so the arena is never observed half-edited.  `_EditScope`
# batches a run of edits into a single apply.


DEF JOURNAL_MIN_CAP = 64


# The replay gate in _apply range-checks op between the first and the last member, so these values
# must stay contiguous.  A new op goes last, updates `OP_HIGHEST`, and needs an arm in the replay
# chain, whose closing else raises if one is missing.
#
# `OP_HIGHEST` IS THE ONE UPPER BOUND, declared beside the values it bounds (§6.1) and read by the
# gate.  A bound pinned in a comment as well is pinned in neither: an op the gate does not admit
# raises `NotImplementedError` from the range check while the arm that would handle it sits in the
# replay chain, so the number is stated once and nowhere else.
cdef enum:
    OP_ADD_ATOM = 1
    OP_DELETE_ATOM = 2
    OP_ADD_BOND = 3
    OP_DELETE_BOND = 4
    OP_SET_ORDER = 5
    OP_SET_CHARGE = 6
    OP_SET_ISOTOPE = 7
    OP_SET_RADICAL = 8
    OP_SET_MAP_NUMBER = 9
    OP_SET_HYDROGENS = 10
    OP_SET_STEREO = 11
    OP_SET_XY = 12
    OP_SET_WEDGE = 13
    OP_SET_STEREO_GROUP = 14
    OP_SET_ATOM_CIP = 15
    OP_SET_BOND_CIP = 16
    OP_SET_ELEMENT = 17
    OP_SET_XYZ = 18
    OP_SET_R_INDEX = 19
    # NAMES THE PARITY SEGMENT AND NOTHING ELSE, so it states no atom and no value.  The one op whose
    # whole effect is on the layout `_apply` allocates.
    OP_WANT_PARITY = 20
    # THE TWO OPS THAT NAME A MODEL AND NO ATOM.  `OP_ADD_CONFORMER` carries the stated `ext_index` in
    # `a` and its destination in `model`; `OP_DROP_CONFORMER` carries the source model in `model`.
    # Both are intercepted in the replay before `_work_index`, which would read `a` as a stable id.
    OP_ADD_CONFORMER = 21
    OP_DROP_CONFORMER = 22
    OP_HIGHEST = OP_DROP_CONFORMER


# deliberately NOT packed: natural alignment gives exactly 20 bytes and aligned loads
cdef struct journal_t:
    uint8_t op
    # IN THE PADDING AFTER `op`, WHERE THE LAYOUT ALREADY HAD ROOM.  Natural alignment puts the next
    # word at offset 4 either way, so the model index costs nothing and no fifth payload word appears.
    # `CONF_MAX_MODELS` is 0xFFFF and the highest index is one below it, so 16 bits is the domain and
    # not a truncation.  Read by `OP_SET_XYZ`, `OP_ADD_CONFORMER` and `OP_DROP_CONFORMER`.
    uint16_t model
    uint32_t a
    uint32_t b
    int32_t v
    # A FOURTH PAYLOAD WORD FOR ONE OP, AND THE ALTERNATIVE WAS WORSE.  `OP_SET_XYZ` carries an atom and
    # three coordinates, which is four values where `a`/`b`/`v` are three.  The other way to fit it is a
    # PAIR of adjacent journal entries written together and read together -- and that makes correctness
    # depend on an ordering invariant nothing in the struct can state, in a buffer that is scanned by
    # index in four separate loops.  A field costs 4 bytes on a TRANSIENT scratch buffer that is freed
    # at seal and never serialised; the invariant would have cost a class of bug.
    #
    # Only `OP_SET_XYZ` reads it.
    int32_t w


cdef struct wedge_edit_t:
    uint32_t src
    uint32_t dst
    uint8_t code


# A BOND'S CIP DESCRIPTOR HAS TO BE CARRIED, WHERE AN ATOM'S DOES NOT, and the asymmetry is the
# layout's and not the chemistry's.  An atom descriptor lives in `atom_t.reserved` and `_apply` copies
# `atom_t` whole (`fresh_atoms[ni] = work[i]`), so it survives with no help.  A bond descriptor lives
# in `halfedge_t.flags`, and `_emit_half` writes that word FRESH from the order on every apply -- so a
# descriptor left to itself is erased by the next `edit()`, silently, which is worse than not storing
# one.  This is the wedge's problem and it gets the wedge's answer: a scratch array seeded from the old
# structure, replayed through the journal, and re-applied after `csr_build`.
#
# Same three fields as `wedge_edit_t` and DELIBERATELY NOT THE SAME ARRAY, though it was tried.  Their
# survival rules differ where it matters: a wedge whose bond was deleted is skipped only when it came
# from the seed (`k < seed_wcount`) and is a KeyError when the caller just set it, whereas a descriptor
# on a molecule whose graph changed is dropped wholesale by the rule below.  Sharing one array would
# have put two different lifetimes behind one `if`.
cdef struct cip_edit_t:
    uint32_t src
    uint32_t dst
    uint8_t code


cdef struct apply_scratch_t:
    void *block
    atom_t *work
    uint8_t *live
    int32_t *newidx
    edge_edit_t *edits
    wedge_edit_t *wedge
    xy_t *wxy
    xyz_t *wxyz
    uint32_t *wext        # one `ext_index` per model of the SEALED molecule, in its order
    uint8_t *wsg
    uint8_t *wpar
    cip_edit_t *cip


cdef struct apply_sizes_t:
    size_t work
    size_t live
    size_t newidx
    size_t edits
    size_t wedge
    size_t wxy
    size_t wxyz
    size_t wext
    size_t wsg
    size_t wpar
    size_t cip


cdef apply_sizes_t _scratch_sizes(uint32_t n_max, uint32_t e_max, uint32_t w_max,
                                   bint want_xy, bint want_sg,
                                   uint32_t conf_models = 0, bint want_parity = False) noexcept nogil:
    """Compute 8-aligned region sizes for the _apply scratch block.

    Both _apply and _apply_scratch_probe call this helper, so they cannot drift apart.
    """
    cdef apply_sizes_t s
    s.work   = align8((n_max if n_max else 1) * sizeof(atom_t))
    s.live   = align8(n_max if n_max else 1)
    s.newidx = align8((n_max if n_max else 1) * sizeof(int32_t))
    s.edits  = align8((e_max if e_max else 1) * sizeof(edge_edit_t))
    s.wedge  = align8((w_max if w_max else 1) * sizeof(wedge_edit_t))
    s.wxy    = align8((n_max if n_max else 1) * sizeof(xy_t)) if want_xy else 0
    # ONE COLUMN PER MODEL, sized here and not at the call site so the arithmetic has one home.  A
    # count rather than a flag: `True` is 1, which is what a one-model caller passes either way.
    s.wxyz   = align8(<size_t> conf_models * (n_max if n_max else 1) * sizeof(xyz_t))
    s.wext   = align8(<size_t> conf_models * sizeof(uint32_t))
    s.wsg    = align8(n_max if n_max else 1) if want_sg else 0
    s.wpar   = align8(n_max if n_max else 1) if want_parity else 0
    # SIZED BY BONDS AND NOT BY A COUNT OF ITS OWN, so this region needs no parameter: one descriptor
    # per bond is the ceiling, and `e_max` already bounds the bonds.  It is unconditional rather than
    # gated on a `want_cip`, because a gate would have to be computed from the journal AND from the old
    # structure's half-edges, and a region that is sometimes absent is the shape that produced the
    # v3-compat segfault in the S-group work -- an index that is right for one version and wrong for
    # the other.  The cost of always having it is 8 bytes per bond of scratch, freed at the end of the
    # apply.
    s.cip    = align8((e_max if e_max else 1) * sizeof(cip_edit_t))
    return s


def _apply_scratch_probe(uint32_t n_max, uint32_t e_max, uint32_t w_max,
                          bint want_xy, bint want_sg, uint32_t conf_models = 0,
                          bint want_parity = False):
    """Return offsets, element sizes, and total for the _apply scratch block.

    Uses _scratch_sizes so the layout matches _apply exactly.  The element
    sizes (work_esz, live_esz, ...) come from C sizeof — tests use them to
    verify sufficiency independently of the alignment arithmetic.
    """
    cdef apply_sizes_t s = _scratch_sizes(n_max, e_max, w_max, want_xy, want_sg, conf_models,
                                          want_parity)
    cdef size_t o_work   = 0
    cdef size_t o_live   = s.work
    cdef size_t o_newidx = s.work + s.live
    cdef size_t o_edits  = s.work + s.live + s.newidx
    cdef size_t o_wedge  = s.work + s.live + s.newidx + s.edits
    cdef size_t o_wxy    = s.work + s.live + s.newidx + s.edits + s.wedge
    cdef size_t o_wxyz   = o_wxy + s.wxy
    cdef size_t o_wext   = o_wxyz + s.wxyz
    cdef size_t o_wsg    = o_wext + s.wext
    cdef size_t o_wpar   = o_wsg + s.wsg
    cdef size_t o_cip    = o_wpar + s.wpar
    cdef size_t total    = o_cip + s.cip
    return {
        'work_off':   o_work,
        'live_off':   o_live,
        'newidx_off': o_newidx,
        'edits_off':  o_edits,
        'wedge_off':  o_wedge,
        'wxy_off':    o_wxy if want_xy else None,
        'wxyz_off':   o_wxyz if conf_models else None,
        'wext_off':   o_wext if conf_models else None,
        'wsg_off':    o_wsg if want_sg else None,
        'wpar_off':   o_wpar if want_parity else None,
        'cip_off':    o_cip,
        'total':      total,
        # Element sizes from C sizeof — independent of alignment arithmetic.
        # Tests compute required_bytes = (count if count else 1) * esz and
        # assert off + required_bytes <= next_off (sufficiency check).
        'work_esz':   sizeof(atom_t),
        'live_esz':   sizeof(uint8_t),
        'newidx_esz': sizeof(int32_t),
        'edits_esz':  sizeof(edge_edit_t),
        'wedge_esz':  sizeof(wedge_edit_t),
        'wxy_esz':    sizeof(xy_t) if want_xy else None,
        'wxyz_esz':   sizeof(xyz_t) if conf_models else None,
        'wext_esz':   sizeof(uint32_t) if conf_models else None,
        'wsg_esz':    sizeof(uint8_t) if want_sg else None,
        'wpar_esz':   sizeof(uint8_t) if want_parity else None,
        'cip_esz':    sizeof(cip_edit_t),
    }


# CIP DESCRIPTORS AT THE PYTHON BOUNDARY.  Index IS the stored code, so these tables define the
# encoding rather than describe it, and the reverse maps below are BUILT FROM THEM -- there is no
# second literal to disagree with the first.  Slot 0 is None, which is why a `.index()` on these is
# also the encoder.
#
# THE TWO TABLES SHARE TWO LETTERS AND ARE STILL TWO TABLES.  M and P are axial descriptors on both
# sides, so one merged table would encode the same letter as two different codes depending on which
# field it was headed for -- and a caller passing an atom's 'M' to a bond would get silence instead of
# a refusal.  Separate tables make that a ValueError.
#
# LOWERCASE IS NOT A SPELLING VARIANT.  r/s are the pseudo-asymmetric descriptors of CIP's auxiliary
# rules: a different determination about a different kind of centre.  Nothing on this path calls
# `.upper()` or `.lower()`, and 'e' is refused rather than read as 'E' -- a caller whose input has been
# case-folded needs to learn that here, because the same folding turns an 'R' into an 'r'.
cdef tuple ATOM_CIP_CODES = (None, 'R', 'S', 'r', 's', 'M', 'P', 'm', 'p')
cdef tuple BOND_CIP_CODES = (None, 'E', 'Z', 'M', 'P')
cdef dict _ATOM_CIP_BY_NAME = {}
cdef dict _BOND_CIP_BY_NAME = {}
cdef uint8_t _cip_i
for _cip_i in range(1, len(ATOM_CIP_CODES)):
    _ATOM_CIP_BY_NAME[ATOM_CIP_CODES[_cip_i]] = _cip_i
for _cip_i in range(1, len(BOND_CIP_CODES)):
    _BOND_CIP_BY_NAME[BOND_CIP_CODES[_cip_i]] = _cip_i


cdef uint8_t _cip_code(object value, dict table, tuple names, str what) except? 0:
    """A descriptor string to its stored code.  `None` is 0 -- "no descriptor" is a value, not an error.

    Refuses anything else, INCLUDING a correct letter for the other kind: 'E' on an atom and 'R' on a
    bond are both mistakes worth hearing about, and the message names the domain it was checked
    against so the caller can see which of the two tables it hit.
    """
    if value is None:
        return 0
    if not isinstance(value, str):
        raise TypeError('a CIP descriptor for %s must be a str or None, not %s'
                        % (what, type(value).__name__))
    cdef object code = table.get(value)
    cdef list known
    cdef Py_ssize_t i
    if code is None:
        # An explicit loop and not a comprehension: in a `cdef` function Cython gives the loop
        # variable an implicit declaration and warns, and the gate for this file is zero warnings.
        known = []
        for i in range(1, len(names)):
            known.append(repr(names[i]))
        raise ValueError('%r is not a CIP descriptor for %s; case is significant, and the accepted '
                         'ones are %s' % (value, what, ', '.join(known)))
    return <uint8_t> <int> code


# The keys an S-group dict may carry.  `_alias` is here because `sgroups` emits it and a caller that
# round-trips what it was given must not be told its own key is unknown -- the read shape and the write
# shape are ONE shape, and a read-only field is how they quietly stop being that.
cdef frozenset _SGROUP_KEYS = frozenset((
    'type', 'subtype', 'name', 'disp', 'disp_tail', 'index', 'ext_index', 'parent',
    'atoms', 'patoms', 'bonds', 'cstates', 'data', 'fields', 'log', '_alias'))


cdef bytes _as_bytes(object value, str what):
    """`bytes` unchanged, `str` encoded UTF-8, anything else refused.

    ACCEPTS `str` FOR CONVENIENCE AND STORES `bytes` ALWAYS, and the asymmetry is deliberate: a caller
    with a name out of a UTF-8 file should not have to encode it, but nothing here may DECODE, because
    a byte that is not valid UTF-8 is exactly what an SDF from another vendor contains and losing it is
    the one thing this segment exists to prevent.  So the conversion is one-way by construction.
    """
    if isinstance(value, bytes):
        return <bytes> value
    if isinstance(value, str):
        return (<str> value).encode('utf8')
    raise TypeError('sgroup %s must be bytes or str, not %s' % (what, type(value).__name__))


cdef bytes _as_title_bytes(object value):
    """A title as the bytes the blob stores: `str` encoded with `surrogateescape`, a buffer copied.

    SEPARATE FROM `_as_bytes` so an S-group name keeps refusing exactly what it always refused.  The
    handler is the whole trick: `title` hands out a `str` decoded the same way, so the round trip is
    exact for every byte a file can hold.
    """
    if isinstance(value, str):
        return (<str> value).encode('utf8', 'surrogateescape')
    if isinstance(value, bytes):
        return <bytes> value
    if isinstance(value, (bytearray, memoryview)):
        return bytes(value)
    raise TypeError('title must be str, bytes, bytearray or memoryview, not %s'
                    % type(value).__name__)


cdef uint16_t _sgroup_number(object value, str what) except? 0:
    """An Sgroup number, or SGROUP_NO_INDEX for none.

    THE BOUND CITED IS SGROUP_INDEX_MAX AND NOT THE FIELD'S WIDTH.  0xFFFF is the sentinel, so the
    greatest number a file may carry is 0xFFFE, and a validator that let 0xFFFF through would store a
    numbered record that reads back as unnumbered.
    """
    cdef long long n = value
    if n == <long long> SGROUP_NO_INDEX:
        return <uint16_t> SGROUP_NO_INDEX
    if n < 0 or n > <long long> SGROUP_INDEX_MAX:
        raise ValueError('sgroup %s must be 0..%d or %d for none, not %d'
                         % (what, SGROUP_INDEX_MAX, SGROUP_NO_INDEX, n))
    return <uint16_t> n


cdef int32_t _fixed_point(object value) except? -1:
    """A display coordinate as the arena's x10000 fixed point.

    `xy_t` IS EXACTLY F10.4 AND THAT IS WHY IT IS THE RIGHT TYPE HERE, not merely an available one: a
    FIELDDISP anchor is written by MDL as F10.4, an int32 of tenths of a thousandth spans +-214748.3647,
    and a float32 would have lost the fourth decimal on a five-digit coordinate -- silently, and only
    for large drawings.
    """
    cdef double x = value
    if x != x or x < -214748.0 or x > 214748.0:
        raise ValueError('display coordinate %r is outside the +-214748 this format holds' % value)
    return <int32_t> round(x * XY_SCALE)


cdef dict _sgroup_normalise(object src, bint alias):
    """Normalise one caller-supplied S-group dict: fill defaults, reject unknown keys, check ranges.

    Module level and not a method because it touches no molecule -- and that is not tidiness: a
    per-molecule version would be reachable only through a container, and the range checks below are the
    kind of thing a test wants to drive directly.
    """
    cdef dict d = dict(src)
    cdef set unknown = set(d) - _SGROUP_KEYS
    if unknown:
        raise ValueError('unknown sgroup key(s): %s' % ', '.join(sorted(unknown)))
    cdef list data = []
    cdef list fields = []
    cdef list log = []
    cdef object x, key, value
    for x in d.get('data', ()):
        data.append(_as_bytes(x, 'data'))
    for key, value in d.get('fields', ()):
        fields.append((_as_bytes(key, 'field key'), _as_bytes(value, 'field value')))
    for x in d.get('log', ()):
        log.append(_as_bytes(x, 'log'))
    cdef dict out = {'type': _as_bytes(d.get('type', b''), 'type'),
                     'subtype': _as_bytes(d.get('subtype', b''), 'subtype'),
                     'name': _as_bytes(d.get('name', b''), 'name'),
                     'disp_tail': _as_bytes(d.get('disp_tail', b''), 'disp_tail'),
                     'disp': d.get('disp'),
                     'atoms': tuple(d.get('atoms', ())),
                     'patoms': tuple(d.get('patoms', ())),
                     'bonds': tuple(d.get('bonds', ())),
                     'cstates': tuple(d.get('cstates', ())),
                     'data': tuple(data), 'fields': tuple(fields), 'log': tuple(log),
                     '_alias': alias or bool(d.get('_alias'))}
    cdef object name
    for name in ('index', 'ext_index', 'parent'):
        out[name] = _sgroup_number(d.get(name, SGROUP_NO_INDEX), name)
    if out['_alias'] and len(out['atoms']) != 1:
        raise ValueError('an atom alias labels exactly one atom, not %d' % len(out['atoms']))
    for name in ('atoms', 'patoms'):
        if len(out[name]) > SGROUP_LIST_MAX:
            raise ValueError('sgroup %s list of %d references exceeds the %d this format holds'
                             % (name, len(out[name]), SGROUP_LIST_MAX))
    # REFUSED AND NOT TRUNCATED.  A record too long for the format is a record this build cannot store,
    # and storing a prefix of it would be the silent loss the segment exists to prevent.
    if 2 * len(out['bonds']) > SGROUP_LIST_MAX or 2 * len(out['cstates']) > SGROUP_LIST_MAX:
        raise ValueError('sgroup pair list exceeds the %d slots this format holds' % SGROUP_LIST_MAX)
    if (len(out['data']) > SGROUP_LIST_MAX or 2 * len(out['fields']) > SGROUP_LIST_MAX
            or len(out['log']) > SGROUP_LIST_MAX):
        raise ValueError('sgroup string list exceeds the %d handles this format holds'
                         % SGROUP_LIST_MAX)
    return out


# --- the standardization hook ------------------------------------------------------------------- #
#
# `MoleculeContainer.standardize()` is a method on a cdef class, so its body cannot be attached from
# outside; and its body is 101 chemistry rules read out of `chython/chemistry/tables/`, so it cannot
# live here either.  The package layout is `core <- chemistry`, and this file importing
# `chython.chemistry` at module scope would invert it.  So the core owns the NAME and
# `chython.chemistry` registers the BODY, exactly as the InChI writer's kekuliser is registered by
# `_ich_set_kekule_fn`.
#
# NOTHING in this file names `chython.chemistry` in an import, at any scope.  A runtime-only
# fallback import was tried and removed: see `_standardize_fn` for why a lazy cycle is still a cycle.

cdef object _standardize_impl = None


def _set_standardize_fn(fn):
    """Register the standardization pass.  `chython.chemistry` calls this at its own init.

    `fn` must have signature `fn(molecule, *, fix_hydrogens=True, fix_tautomers=True) -> bool` and must
    mutate `molecule` in place, recording onto `molecule.log`.
    """
    global _standardize_impl
    _standardize_impl = fn


cdef object _standardize_fn():
    """The registered pass, or an error that names the package which registers it.

    NO FALLBACK IMPORT HERE, and the omission is the point.  Importing `chython.chemistry` on
    demand would work, and it was written that way first, but it makes the bottom layer name the
    package above it -- a cycle that is merely lazy rather than absent.  It also breaks the core's
    hard invariant that `chython/core/` imports nothing from this distribution but `chython.core`
    (`core/test/test_no_chython_two_imports.py`); an exception carved out for one sibling is an
    exception someone later widens.

    The cost is this message, for the narrow case of importing `chython.core` directly and calling
    `standardize()` without the package that implements it -- and the message teaches the layering
    instead of hiding it.  `import chython` registers the pass, so no ordinary caller sees this.
    """
    global _standardize_impl
    if _standardize_impl is None:
        raise ImportError('standardization is implemented in `chython.chemistry`, which registers '
                          'it on import; `import chython.chemistry` (or `import chython`) first. '
                          'The core owns the method name, not the 101 chemistry rules behind it.')
    return _standardize_impl


cdef object _canonicalize_impl = None


def _set_canonicalize_fn(fn):
    """Register the canonicalization façade.  `chython.chemistry` calls this at its own init.

    `fn` must have signature
    `fn(molecule, *, fix_tautomers=True, keep_kekule=False) -> bool` and must mutate `molecule` in
    place, recording onto `molecule.log`.  Registered rather than imported for the reason `_set_standardize_fn` gives:
    the direction is `core <- chemistry`, and it orchestrates the 101 rules the core does not know
    about.
    """
    global _canonicalize_impl
    _canonicalize_impl = fn


cdef object _canonicalize_fn():
    global _canonicalize_impl
    if _canonicalize_impl is None:
        raise ImportError('canonicalization is implemented in `chython.chemistry`, which registers '
                          'it on import; `import chython.chemistry` (or `import chython`) first.')
    return _canonicalize_impl


# --- the depiction hook ------------------------------------------------------------------------- #
#
# Same shape and the same reason as the two above, one layer further out.  A 2D layout is a JavaScript
# bundle on QuickJS or a third-party toolkit; the core neither imports one nor knows what an engine
# name means, and `chython.depict` is above it in the order `core <- ... <- depict`.  So the core owns
# the method NAMES and `chython.depict` registers the BODIES at its own import time.
#
# ELEVEN SLOTS IN THREE GROUPS, ONE SETTER, ONE GUARD PER GROUP.  The setter is one function because a
# group is one registration -- see `_set_depict_fns` -- and each guard is one function because its
# message is one fact stated once, per RULES.md 6.

cdef object _clean2d_impl = None
cdef object _layout2d_impl = None
cdef object _rescale2d_impl = None
cdef object _reaction_clean2d_impl = None
cdef object _reaction_layout2d_impl = None

cdef object _depict_impl = None
cdef object _scene_impl = None
cdef object _reaction_depict_impl = None
cdef object _reaction_scene_impl = None

cdef object _depict3d_impl = None
cdef object _view3d_impl = None


def _set_depict_fns(*, clean2d=None, layout2d=None, rescale2d=None, reaction_clean2d=None,
                    reaction_layout2d=None, depict=None, scene=None, reaction_depict=None,
                    reaction_scene=None, depict3d=None, view3d=None):
    """Register the depiction entry points onto the sealed containers.

    ONE function for a whole group, rather than a setter per name, because a group is one registration:
    a `depict` that had registered two of the layout four would leave `rxn.clean2d()` raising while
    `mol.clean2d()` worked, and no caller could tell why.  `chython.depict` calls this at its own import
    time.

    THE LAYOUT GROUP.  The molecule pair takes `(molecule, *, engine=None, force=False)`; `clean2d`
    stores and returns None, `layout2d` returns `{n: (x, y)}` and stores nothing.  `rescale2d` takes
    the molecule alone, rewrites the stored plane to average bond length 0.825 and answers whether it
    rescaled.  The reaction pair takes `(reaction, *, engine=None, force=False)`; `layout2d` returns
    `(planes, arrow, signs)` and `clean2d` stores the planes and returns `(arrow, signs)`.

    THE 3D GROUP.  `depict3d` takes `(molecule, index)` and returns an X3DOM document; `view3d` takes
    `(molecule, index, width, height)` and returns a notebook widget.  A molecule side only -- a
    reaction has no conformer -- and no layout member, because a conformer is read, never computed
    here: `chython.interop.conformers.generate_conformers` is the call that makes one.

    THE DRAWING GROUP.  `scene`/`depict` take `(molecule, *, style=None, plane=None, log=None)` and
    return a `Scene` and an SVG document respectively; `reaction_scene`/`reaction_depict` take
    `(reaction, *, style=None, log=None)`.  None of the four stores anything -- a molecule with no
    layout is drawn against a temporary and says so in `log`.

    ENFORCED HERE, AT THE CALL, and BOTH-OR-NEITHER PER GROUP.  Offering any member of a group obliges
    all of it.  The rule is per group and not across every argument this function will ever take: a call
    that offers nothing from a group leaves that group alone, which is what lets the two register
    separately rather than making one of them a silent half-registration.  Each group has its own
    all-or-nothing check below, and a new group adds a third beside them.

    The argument for checking it here is the argument the paragraph above makes for one setter: read
    time is too late.  `_reaction_depict_fns` also refuses a half-registered pair, but it refuses at the
    first `rxn.clean2d()`, which is somewhere else entirely from the hook that caused it.
    """
    global _clean2d_impl, _layout2d_impl, _rescale2d_impl
    global _reaction_clean2d_impl, _reaction_layout2d_impl
    global _depict_impl, _scene_impl, _reaction_depict_impl, _reaction_scene_impl
    global _depict3d_impl, _view3d_impl
    cdef list missing

    if (clean2d is not None or layout2d is not None or rescale2d is not None
            or reaction_clean2d is not None or reaction_layout2d is not None):
        missing = []
        if clean2d is None:
            missing.append('clean2d')
        if layout2d is None:
            missing.append('layout2d')
        if rescale2d is None:
            missing.append('rescale2d')
        if reaction_clean2d is None:
            missing.append('reaction_clean2d')
        if reaction_layout2d is None:
            missing.append('reaction_layout2d')
        if missing:
            raise ValueError('the five layout functions are ONE registration and this call gave only '
                             'some of them; missing: %s. A half-registered layout leaves one of '
                             '`mol.clean2d()` / `rxn.clean2d()` raising while the other works, and the '
                             'caller cannot see why.' % ', '.join(missing))
        _clean2d_impl = clean2d
        _layout2d_impl = layout2d
        _rescale2d_impl = rescale2d
        _reaction_clean2d_impl = reaction_clean2d
        _reaction_layout2d_impl = reaction_layout2d

    if (depict is not None or scene is not None or reaction_depict is not None
            or reaction_scene is not None):
        missing = []
        if depict is None:
            missing.append('depict')
        if scene is None:
            missing.append('scene')
        if reaction_depict is None:
            missing.append('reaction_depict')
        if reaction_scene is None:
            missing.append('reaction_scene')
        if missing:
            raise ValueError('the four drawing functions are ONE registration and this call gave only '
                             'some of them; missing: %s. A half-registered drawing leaves one of '
                             '`mol.depict()` / `rxn.depict()` raising while the other works, and the '
                             'caller cannot see why.' % ', '.join(missing))
        _depict_impl = depict
        _scene_impl = scene
        _reaction_depict_impl = reaction_depict
        _reaction_scene_impl = reaction_scene

    if depict3d is not None or view3d is not None:
        missing = []
        if depict3d is None:
            missing.append('depict3d')
        if view3d is None:
            missing.append('view3d')
        if missing:
            raise ValueError('the two 3D functions are ONE registration and this call gave only some '
                             'of them; missing: %s. `view3d` IS `depict3d` in a widget, so half of the '
                             'pair is a notebook that renders nothing.' % ', '.join(missing))
        _depict3d_impl = depict3d
        _view3d_impl = view3d


cdef object _depict_fn(object impl, str name):
    """The registered body, or an error that names the package which registers it.

    NO FALLBACK IMPORT, for the reason `_standardize_fn` states at length: a lazy cycle is still a
    cycle, and `core/test/test_no_chython_two_imports.py` holds the core to importing nothing from
    this distribution but `chython.core`.

    `ImportError` and not `NotImplementedError`: the method is not unimplemented, it is unregistered,
    and the one caller who can see this -- somebody who imported `chython.core` alone -- fixes it with
    an import.  The message says what implements it and what to type.
    """
    if impl is None:
        raise ImportError('`%s` is implemented in `chython.depict`, which registers it on import; '
                          '`import chython.depict` (or `import chython`) first. The core owns the '
                          'method name, not the layout engine behind it.' % name)
    return impl


def _reaction_depict_fns():
    """The reaction-side pair, for `chython/core/reaction.py`, which cannot see a `cdef` global.

    A Python-visible `def` in the extension read by a Python module in the same package -- the mirror
    of `_set_reaction_factory(ReactionContainer)` at the foot of that file, and the same reason.

    BOTH OR NEITHER, which is the read-time half of `_set_depict_fns`' one-registration argument: a
    caller who can lay a reaction out but not store it has a half-registered `depict`, and finding
    that out at the second call is worse than at the first.
    """
    return (_depict_fn(_reaction_clean2d_impl, 'ReactionContainer.clean2d'),
            _depict_fn(_reaction_layout2d_impl, 'ReactionContainer.layout2d'))


def _reaction_draw_fns():
    """The reaction-side DRAWING pair, for `chython/core/reaction.py`.  Same door, same reason.

    BOTH OR NEITHER again, and read-time as well as call-time: a caller who can build a reaction's
    `Scene` but cannot serialize it has a half-registered `depict`.
    """
    return (_depict_fn(_reaction_depict_impl, 'ReactionContainer.depict'),
            _depict_fn(_reaction_scene_impl, 'ReactionContainer.scene'))


cdef object _isomers_impl = None


def _set_isomers_fn(fn):
    """Register the isomer-placement pass.  `chython.chemistry` calls this at its own init.

    `fn` must have signature `fn(molecule) -> bool` and must mutate `molecule` in place.
    A SETTER OF ITS OWN rather than a second argument to `_set_canonicalize_fn`, unlike the hydrogen
    pair below: those two share one setter because a tree where `implicify_hydrogens()` exists and
    `explicify_hydrogens()` raises is a worse diagnostic than neither existing.  These two do not
    share that property -- `canonicalize()` is a pipeline and this is one stage of it, and a caller
    may reasonably want the stage without the pipeline: a corpus already kekulised and standardized,
    wanting only the placement unified.
    """
    global _isomers_impl
    _isomers_impl = fn


cdef object _isomers_fn():
    global _isomers_impl
    if _isomers_impl is None:
        raise ImportError('isomer placement is implemented in `chython.chemistry`, which registers '
                          'it on import; `import chython.chemistry` (or `import chython`) first.')
    return _isomers_impl


cdef object _reconstruct_impl = None


def _set_reconstruct_fn(fn):
    """Register template-driven mapping reconstruction.  `chython.reactions` calls this at its own init.

    `fn` must have signature
    `fn(reaction, *, max_size_ratio=5., min_filter_size=42) -> tuple[str, ...]`
    and must mutate `reaction` in place -- it canonicalizes both sides and writes map numbers.

    ITS OWN SETTER, and NOT another key in `_set_reactions_fns`, because the two hooks serve different
    classes.  That dict exists for names on a `cdef class`, which cannot be extended from outside;
    this one is read by `chython/core/reaction.py`, a plain Python module, and the accessor below is
    the Python-visible `def` it reads -- the same arrangement as `_reaction_depict_fns` and for the
    same reason.  Accepts `None` so a test can exercise the unregistered branch and put the body back.
    """
    global _reconstruct_impl
    _reconstruct_impl = fn


def _reaction_reconstruct_fn():
    """The reconstruction body, for `chython/core/reaction.py`, which cannot see a `cdef` global.

    NO FALLBACK IMPORT HERE, for the reason `_standardize_fn` states at length: a lazy import inside an
    accessor is still an import cycle, deferred to the first call, where the diagnostic is worse.  The
    core owns the method name; `chython.reactions` owns the corpus of templates the answer comes from.
    """
    if _reconstruct_impl is None:
        raise ImportError('mapping reconstruction is implemented in `chython.reactions`, which '
                          'registers it on import; `import chython.reactions` (or `import chython`) '
                          'first.  The core owns the method name, not the corpus of reaction '
                          'templates the reconstruction searches.')
    return _reconstruct_impl


cdef object _attention_impl = None


def _set_attention_fn(fn):
    """Register the neural atom-atom mapper.  `chython.reactions` calls this at its own init.

    `fn` must have signature
    `fn(reaction, *, multiplier=1.75, keep_reactant_mapping=False, threads=None) -> MappingResult`
    and must mutate `reaction` in place -- it writes map numbers and touches nothing else.

    ITS OWN SETTER, for the reason `_set_reconstruct_fn` gives: this hook is read by
    `chython/core/reaction.py`, a plain Python module, and not by a `cdef class`.  Accepts `None` so a
    test can exercise the unregistered branch and put the body back.
    """
    global _attention_impl
    _attention_impl = fn


def _reaction_attention_fn():
    """The mapper's body, for `chython/core/reaction.py`, which cannot see a `cdef` global.

    NO FALLBACK IMPORT, for the reason `_reaction_reconstruct_fn` states.  The message names the
    package that registers a body; the body's own message names the `chython[mapping]` extra when the
    runtime and the weights are what is missing.  A caller who imported chython and skipped the extra
    must read the second, so this branch never mentions it.
    """
    if _attention_impl is None:
        raise ImportError('attention mapping is implemented in `chython.reactions`, which registers '
                          'it on import; `import chython.reactions` (or `import chython`) first.  The '
                          'core owns the method name, not the model the mapping comes from.')
    return _attention_impl


cdef object _implicify_impl = None
cdef object _explicify_impl = None


def _set_hydrogens_fns(implicify, explicify):
    """Register the two hydrogen passes together.  `chython.chemistry` calls this at its own init.

    ONE SETTER FOR THE PAIR, unlike the two above, because they are one decision: a tree where
    `implicify_hydrogens()` exists and `explicify_hydrogens()` raises `ImportError` is a worse
    diagnostic than neither existing.  Both signatures are
    `fn(molecule) -> int` and both mutate `molecule` in place.

    THEY ARE REGISTERED DESPITE READING NO RULE TABLE, unlike the two passes above, and that is not the
    usual reason for the split.  The line is "does this belong to a pack": implicify is a stage of
    `canonicalize()`, shares its recording conventions and its `hydrogens:` rule-id namespace with the 101
    rules, and would have to be imported by `chemistry` wherever it lived -- so it belongs to the
    standardization pack whose home is `chemistry`, and explicify follows its inverse.  A table is the
    usual EVIDENCE that something belongs to a pack rather than the criterion for it.
    """
    global _implicify_impl, _explicify_impl
    _implicify_impl = implicify
    _explicify_impl = explicify


cdef object _hydrogens_fn(bint explicit):
    global _implicify_impl, _explicify_impl
    cdef object fn = _explicify_impl if explicit else _implicify_impl
    if fn is None:
        raise ImportError('the hydrogen passes are implemented in `chython.chemistry`, which '
                          'registers them on import; `import chython.chemistry` (or `import chython`) '
                          'first.')
    return fn


cdef object _valence_impl = None


def _set_valence_fn(fn):
    """Register the valence report.  `chython.chemistry` calls this at its own init.

    `fn` must have signature `fn(molecule) -> list[tuple[int, str]]` and must NOT edit `molecule`.

    THE ONLY REGISTERED SURFACE THAT IS A QUESTION RATHER THAN A PASS, and that is why it gets a setter
    of its own rather than joining one of the pairs above.  Every other hook here mutates and answers
    "did anything change"; this one answers "what does the valence collection say about what is already
    stored" and answering it must leave the molecule alone.  A shared setter would put a read-only
    function and an in-place one behind the same name, and the next reader would have to check which.

    Registered rather than compiled in for the reason the whole hook block gives, plus one specific to
    it: the verdicts come from `valence_check` in the core's generated tables, but the DECISION to report
    an aromatic atom as `'unknown'` rather than checking it is standardization policy, and policy lives
    in `chemistry`.
    """
    global _valence_impl
    _valence_impl = fn


cdef object _valence_fn():
    global _valence_impl
    if _valence_impl is None:
        raise ImportError('the valence report is implemented in `chython.chemistry`, which registers '
                          'it on import; `import chython.chemistry` (or `import chython`) first.')
    return _valence_impl


cdef dict _salts_impl = {}


def _set_salts_fns(**fns):
    """Register the salt passes.  `chython.chemistry` calls this at its own init.

    Two names: `split_salts` edits and answers `bool`, `decompose_salts` reports and edits nothing.
    Keyword-only and dict-backed rather than positional, so a third pass added here cannot silently take
    an earlier one's slot -- a keyword cannot be swapped in the way a positional pair can.
    """
    global _salts_impl
    _salts_impl = dict(fns)


cdef object _salts_fn(str name):
    global _salts_impl
    cdef object fn = _salts_impl.get(name)
    if fn is None:
        raise ImportError('the salt passes are implemented in `chython.chemistry`, which registers '
                          'them on import; `import chython.chemistry` (or `import chython`) first. '
                          'The core owns the method name, not the rows of salt knowledge behind it.')
    return fn


cdef object _resonance_impl = None


def _set_resonance_fn(fn):
    """Register the charge-separation repair.  `chython.chemistry` calls this at its own init.

    One function and one setter: `fn(molecule) -> bool` is the whole resonance surface the core owns.

    WHY IT IS HOOKED AT ALL, `standardize()` not calling it and nothing in the library doing so either:
    a caller holds a molecule, not a package, and having to remember which pass is a method and which is
    a module function is the API describing the layer boundary rather than the chemistry.  `bool` and
    `molecule.log`, the contract every pass above has.

    `saturate` is deliberately NOT here.  It is bond perception for a file that gave connectivity and no
    orders, so its callers are the coordinate readers, and where those hand their record over is still
    being designed -- a method now would pin the shape before the design picks one.
    """
    global _resonance_impl
    _resonance_impl = fn


cdef object _resonance_fn():
    global _resonance_impl
    if _resonance_impl is None:
        raise ImportError('resonance repair is implemented in `chython.chemistry`, which registers it '
                          'on import; `import chython.chemistry` (or `import chython`) first.  The core '
                          'owns the method name, not the charge-separation rules behind it.')
    return _resonance_impl


cdef dict _sgroup_impl = {}


def _set_sgroup_fns(**fns):
    """Register the CTfile data-label helpers.  `chython.formats` calls this at its own init.

    THE ONLY HOOK REGISTERED BY `formats`, and it exists because the storage and the convention live in
    different layers.  The core owns S-group STORAGE -- `sgroups()` and `set_sgroups()` read and replace
    the arena's own segment for every S-group kind there is -- while what a `DAT` record means, how a
    multi-value `FIELDDATA` is spelled and where `FIELDDISP` puts the anchor is CTfile knowledge, which
    may not move down here.  Without the hook the two halves of one job are a method and a module
    function, which tells the caller about the layer boundary instead of about S-groups.

    Keyword-only and dict-backed for `_set_salts_fns`' reason.  `add_data_sgroup` appends and answers the
    record; `data_sgroups` reads and answers a list -- so unlike the pairs above these two are not
    interchangeable in type, and the keyword is here for uniformity rather than to catch a swap.
    """
    global _sgroup_impl
    _sgroup_impl = dict(fns)


cdef object _sgroup_fn(str name):
    global _sgroup_impl
    cdef object fn = _sgroup_impl.get(name)
    if fn is None:
        raise ImportError('the data-label helpers are implemented in `chython.formats`, which registers '
                          'them on import; `import chython.formats` (or `import chython`) first.  The '
                          'core owns S-group storage, not the CTfile `DAT` convention.')
    return fn


cdef object _protomers_impl = None


def _set_protomers_fn(fn):
    """Register the acid/base pass.  `chython.chemistry` calls this at its own init.

    One function and one setter, not a dict: `neutralize` is the whole protomer surface the core owns, and
    `fn(molecule, *, keep_charge=True) -> bool` has no sibling that could be swapped for it.
    """
    global _protomers_impl
    _protomers_impl = fn


cdef object _protomers_fn():
    global _protomers_impl
    if _protomers_impl is None:
        raise ImportError('`neutralize` is implemented in `chython.chemistry`, which registers it on '
                          'import; `import chython.chemistry` (or `import chython`) first.  The core owns '
                          'the method name, not the acid/base table behind it.')
    return _protomers_impl


cdef dict _reactions_impl = {}


def _set_reactions_fns(**fns):
    """Register the reaction-enumeration surface.  `chython.reactions` calls this at its own init.

    Keyword-only and dict-backed, the shape `_set_salts_fns` argues for: both of these take a molecule
    and answer an iterable, so a swapped positional pair would be no type error anywhere and would
    silently report functional groups where the caller asked for reactions.

    The core owns `react`, `functional_groups`, `functional_group_hits`, `protective_groups`,
    `protective_group_hits`, `deprotect`, `sticky_fragments`, `sticky_linkers` and `@` -- names on a
    `cdef class`, which cannot be extended from outside -- and `chython.reactions` owns the corpus behind
    them.  A special method is the reason this hook exists at all rather than the package simply exporting
    functions: `mol @ other` resolves through the type's slot, so `__matmul__` has to be defined here or
    `@` is a TypeError however many templates are loaded.

    EIGHT NAMES, THREE QUESTIONS.  `react`, `functional_groups` and `functional_group_hits` read the
    reaction corpus -- one enumerator, no per-scope variants of it -- `protective_groups`,
    `protective_group_hits` and `deprotect` are the protecting-group half, a different question over a
    different table and not another view of the same one -- and `sticky_fragments`/`sticky_linkers` cut a
    coupling handle and cap it with an R, a third question over `roles.tsv`.  A `*_hits` name is the same
    scan as the dict beside it with each row's id kept, and folds to it.
    """
    global _reactions_impl
    _reactions_impl = dict(fns)


cdef object _reactions_fn(str name):
    global _reactions_impl
    cdef object fn = _reactions_impl.get(name)
    if fn is None:
        raise ImportError('reaction enumeration is implemented in `chython.reactions`, which '
                          'registers it on import; `import chython.reactions` (or `import chython`) '
                          'first.  The core owns the method name and the SMIRKS reader, not the '
                          'corpus of templates behind them.')
    return fn


cdef dict _featurizer_impl = {}


def _set_featurizer_fns(**fns):
    """Register the F3 descriptor bodies.  Called by `chython.chemistry` on import.

    Keyword-only and dict-backed: the ten functions all share a signature that takes a molecule and
    returns a value, so a positional list would give no type error on a swap and a transposed pair
    would return the wrong descriptor for the right name with no visible failure until a caller
    cross-checked two properties.
    """
    global _featurizer_impl
    _featurizer_impl = dict(fns)


cdef object _featurizer_fn(str name):
    global _featurizer_impl
    cdef object fn = _featurizer_impl.get(name)
    if fn is None:
        raise ImportError(f'`{name}` is implemented in `chython.chemistry`, which registers it '
                          'on import; `import chython.chemistry` (or `import chython`) first.')
    return fn


def _featurizer_fn_for_test(str name):
    """Reach `_featurizer_fn` from Python.  Used only by the injection test."""
    return _featurizer_fn(name)


# --- the interop hook --------------------------------------------------------------------------- #
#
# Same shape and the same reason as every hook above, at the outermost layer: a converter imports RDKit
# or starts a JVM, `chython.interop` is a peer of `depict` above the core, and `to_rdkit` is a name on a
# `cdef class`.  So the core owns the five NAMES and `chython.interop` registers the BODIES.
#
# WHAT IS REGISTERED IS THE PACKAGE'S OWN DISPATCHER -- `interop.rdkit` and not `interop._rdkit.to_rdkit`
# -- because that function already chooses the direction by testing for a container, and `self` is one.
# The method is therefore the export half of the published callable and cannot drift from it; the lazy
# import of the toolkit also stays where it already is, inside the dispatcher.

cdef dict _interop_impl = {}


def _set_interop_fns(**fns):
    """Register the toolkit converters.  `chython.interop` calls this at its own init.

    Keyword-only and dict-backed, the shape `_set_salts_fns` argues for and with its own reason: all
    five take a container and answer a foreign object, so a positional list would give no type error on
    a swap -- `mol.to_rdkit()` would hand back an Indigo object, and the caller would find out several
    toolkit calls later.

    THE METHODS ARE A SECOND DOOR ONTO ONE IMPLEMENTATION, not a second implementation.  A method on a
    `cdef class` needs its body attached from outside, which is what this hook is for, and the whole of
    it is the export direction -- reading a foreign object stays a call to
    `interop.rdkit(x)`, because a method on a chython container is a strange door for an object that is
    not one yet.
    """
    global _interop_impl
    _interop_impl = dict(fns)


cdef object _interop_fn(str name):
    global _interop_impl
    cdef object fn = _interop_impl.get(name)
    if fn is None:
        raise ImportError(f'`{name}` conversion is implemented in `chython.interop`, which registers '
                          'it on import; `import chython.interop` (or `import chython`) first.  The '
                          'core owns the method name, not the toolkit behind it.')
    return fn


def _reaction_interop_fn(str name):
    """The converter, for `chython/core/reaction.py`, which cannot see a `cdef` global.

    The same arrangement as `_reaction_reconstruct_fn` and for the same reason: a Python-visible `def`
    in the extension, read by a Python module in the same package.
    """
    return _interop_fn(name)


# `Log`, `LogRecord` and the severities, bound on first use.  `chython/core/__init__.py` imports `._core`
# FIRST, so a module-level import here would run halfway through the package's own initialisation -- the
# reason `_smirks_patch.pxi:76` gives.  `_log.py` imports nothing itself; reaching it still goes through
# the package that is still executing.
#
# One translation unit, and this fragment is included before `_kekule`, `_thiele` and the four readers --
# so those five reach `LogRecord` through these globals rather than each growing a lazy import of its own.
cdef object _MC_LOG = None
cdef object _MC_RECORD = None
cdef object _MC_INFO = None
cdef object _MC_LOST = None
cdef object _MC_REPAIRED = None
cdef object _MC_REFUSED = None


cdef int mc_lazy_log_imports() except -1:
    global _MC_LOG, _MC_RECORD, _MC_INFO, _MC_LOST, _MC_REPAIRED, _MC_REFUSED
    # declared before they are imported: an import statement binds a name Cython never saw declared,
    # and `warn.undeclared` is on
    cdef object Log
    cdef object LogRecord
    cdef object INFO
    cdef object LOST
    cdef object REPAIRED
    cdef object REFUSED
    if _MC_LOG is None:
        from chython.core._log import Log, LogRecord, INFO, LOST, REPAIRED, REFUSED
        _MC_LOG = Log
        _MC_RECORD = LogRecord
        _MC_INFO = INFO
        _MC_LOST = LOST
        _MC_REPAIRED = REPAIRED
        _MC_REFUSED = REFUSED
    return 0


cdef inline object mc_lost():
    """`LOST`, importing it if nobody has yet.  A call site reading `_MC_LOST` directly would evaluate
    the argument BEFORE `_log_event` runs the import, and pass `None` on the first event of the run."""
    mc_lazy_log_imports()
    return _MC_LOST


cdef inline object mc_info():
    """`INFO`, on the same terms as `mc_lost`.  For a field this reader read and deliberately did not
    store because it states nothing storable -- so a count of LOSSES does not see it."""
    mc_lazy_log_imports()
    return _MC_INFO


cdef inline object mc_repaired():
    """`REPAIRED`, on the same terms as `mc_lost`."""
    mc_lazy_log_imports()
    return _MC_REPAIRED


cdef inline object mc_refused():
    """`REFUSED`, on the same terms as `mc_lost`."""
    mc_lazy_log_imports()
    return _MC_REFUSED


cdef inline object mc_record(str rule, tuple atoms, str message, object severity=None):
    """One `LogRecord`, with the lazy import done for the caller.  `stage` is deliberately absent here:
    the pass names it once, where it folds these records onto `mol.log` (`absorb('kekule', ...)`), rather
    than at every emit site where one of them could be misspelled."""
    mc_lazy_log_imports()
    return _MC_RECORD(rule, atoms, message, _MC_INFO if severity is None else severity)


cdef class MoleculeContainer:
    cdef Structure _structure
    cdef journal_t *_journal
    cdef uint32_t _journal_len
    cdef uint32_t _journal_cap
    cdef uint32_t _next_id        # next stable id to hand out; never decreases
    cdef uint32_t _first_pending  # lowest stable id not yet in the arena
    cdef uint32_t _gen
    cdef uint32_t _scope_depth    # nesting depth of open edit() scopes
    cdef dict _index_of           # n -> index
    cdef list _numbers         # index -> n
    cdef dict _order_cache        # atoms_order result, valid while _order_gen == _gen
    cdef uint32_t _order_gen
    cdef bytes _identity_cache    # canonical_bytes result, valid while _identity_gen == _gen
    cdef uint32_t _identity_gen
    cdef str _smiles_cache        # str(self), the EMPTY spec only, valid while _smiles_gen == _gen
    cdef uint32_t _smiles_gen
    cdef object _log              # this handle's Log, created on first use; see the `log` property
    cdef dict _meta               # record metadata, created on first access; see the `meta` property
    cdef bint _representation_change   # set ONLY by kekule/thiele around their own edit scope
    cdef uint32_t _conf_adds      # add_conformer ops in the open journal
    cdef uint32_t _conf_drops     # drop_conformer ops in the open journal
    # THE SOURCE MODELS THIS SESSION DROPPED, None until one is.  Counters and a set rather than a
    # scan of the journal: `set_xyz` is called once per atom per model, so an exclusivity rule that
    # re-scanned would be quadratic in the one op that is written most.
    cdef set _conf_dropped

    def __cinit__(self):
        self._journal = NULL
        self._journal_len = 0
        self._journal_cap = 0
        self._next_id = 1
        self._first_pending = 1
        self._gen = 0
        self._scope_depth = 0
        self._index_of = {}
        self._numbers = []
        self._order_cache = None
        self._order_gen = 0
        self._identity_cache = None
        self._identity_gen = 0
        self._smiles_cache = None
        self._smiles_gen = 0
        self._log = None
        self._meta = None
        self._representation_change = False
        self._conf_adds = 0
        self._conf_drops = 0
        self._conf_dropped = None
        self._structure = structure_alloc(0, 0, False)

    def __dealloc__(self):
        PyMem_Free(self._journal)
        self._journal = NULL

    cdef int _append(self, uint8_t op, uint32_t a, uint32_t b, int32_t v, int32_t w = 0,
                     uint16_t model = 0) except -1:
        cdef uint32_t cap
        cdef journal_t *grown
        if self._journal_len == self._journal_cap:
            cap = JOURNAL_MIN_CAP if self._journal_cap == 0 else self._journal_cap * 2
            grown = <journal_t *> PyMem_Realloc(self._journal, cap * sizeof(journal_t))
            if grown is NULL:
                raise MemoryError('journal reallocation failed')
            self._journal = grown
            self._journal_cap = cap
        cdef journal_t *rec = self._journal + self._journal_len
        rec.op = op
        rec.a = a
        rec.b = b
        rec.v = v
        rec.w = w
        rec.model = model
        self._journal_len += 1
        return 0

    cdef bint _has(self, uint32_t n) noexcept:
        # live if the arena knows it, or if it is a pending id from this journal
        if self._first_pending <= n and n < self._next_id:
            return True
        return n in self._index_of

    cdef int _require(self, uint32_t n) except -1:
        if not self._has(n):
            raise KeyError(n)
        return 0

    cdef int _discard(self) except -1:
        self._journal_len = 0
        self._first_pending = self._next_id
        self._conf_adds = 0
        self._conf_drops = 0
        self._conf_dropped = None
        return 0

    cdef int _maybe_apply(self) except -1:
        if self._scope_depth == 0:
            self._apply()
        return 0

    cdef int _require_clean(self) except -1:
        if self._journal_len:
            raise RuntimeError('the container has pending edits; the arena still holds the '
                               'pre-scope state, so this read would answer from stale data')
        return 0

    # THERE IS NO `_require_kekule` HERE, AND THE ABSENCE IS DELIBERATE.  No core reader of a bond
    # order needs to refuse order 4: hybridization answers 4, the feature words separate order 4 from
    # a dative bond on the aromatic bit, the canonical bond word folds the flag in, and the stereo
    # surface takes a five-line electron-budget correction (`_stereo.pxi`, `spent`) rather than a gate
    # -- axis detection tests degree, hydrogen count and ring membership and never an order.  A gate
    # there would refuse most drug-like molecules.
    #
    # NOR IS THERE ONE FOR A CALLER OUTSIDE THE CORE THAT CANNOT REPRESENT AN AROMATIC BOND.  A CTfile
    # bond block spells one as bond type 4, both MDL writers write it and the readers accept it, so
    # such a refusal makes `write(read(x))` fail on input chython itself read without complaint.  An
    # unused refusal reads as a policy the core has, and this one is not one.
    #
    # `is_kekule` and `aromatic_bond_count` are the whole surface, and the honest shape of it: a
    # caller that cannot represent order 4 asks which representation it holds and decides for itself.
    cdef inline uint32_t _work_index(self, uint32_t n, uint32_t n_atoms_old) except? 0xFFFFFFFF:
        if n >= self._first_pending:
            return n_atoms_old + (n - self._first_pending)
        return <uint32_t> self._index_of[n]

    cdef list _harvest_parities(self, set touched):
        """Snapshot every CONFIGURED stereo unit as
        `(anchor n, kind, parity, refs as numbers, unnamed mask)`.

        Returns None when there is nothing to carry, which is the common case: the gate below is the
        presence of SEG_PARITY, so a molecule that states no parity never perceives units.

        THE MASK IS PART OF THE FRAME, not a convenience (ruling F69).  `refs` collapses an unnamed
        direction and a slot that is no direction at all to the same `None`, and the re-basing
        correspondence has to tell them apart -- an empty slot may only correspond to an empty slot.
        So `spare >> SU_UNNAMED_SHIFT` travels with the refs it explains.

        BUILT UNMARKED.  This reads `kind`, `refs`, `anchor` and the unnamed mask, all of them pure
        constitution, so it perceives through `ensure_stereo_units_unmarked` and never pays the
        budgeted stereogenicity search (ruling F70).  A reader that wants `stereogenic` calls
        `ensure_stereo_units` later and the marking pass runs then, once.

        UNITS, NOT ATOMS, and that distinction is the drop rule (see `_replay_parities`).  A unit
        that is skipped here is a unit the apply will not touch AT ALL: the `touched` skip suppresses
        the DROP just as much as the re-base, and that is the intended reading rather than a leak.  A
        caller who states a parity inside the very edit that destroys the frame is stating a sign
        against the NEW molecule, not asking for the old one to be carried over, so the apply has
        nothing to carry and nothing to clear -- and the bit it leaves behind is exactly ruling F66's
        second case, a sign whose frame has not yet existed.  Such a bit is readable and survives
        `to_bytes`; reporting it, and clearing it if the consumer wants that, is `validate_stereo`'s
        business and never the apply's
        (`test_a_parity_stated_in_the_edit_that_destroys_the_frame_is_the_callers_own`).
        """
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef stereo_unit_t *units
        cdef stereo_unit_t *u
        cdef uint8_t *par
        cdef uint32_t i, k, count, r
        cdef list numbers = self._numbers
        cdef list snapshots
        cdef list refs
        cdef object n
        if n_atoms == 0:
            return None
        if not structure_has(self._structure, SEG_PARITY):
            return None
        # REALLOCATES THE ARENA (ruling F60): pointers into it are dead from here, and are
        # re-fetched below rather than carried across.  UNMARKED (ruling F70): nothing below reads a mark.
        ensure_stereo_units_unmarked(self._structure)
        count = structure_stereo_unit_count(self._structure)
        if count == 0:
            return None
        units = structure_stereo_units(self._structure)
        par = structure_parities(self._structure)
        snapshots = []
        for k in range(count):
            u = &units[k]
            if not par[u.anchor]:
                continue
            n = numbers[u.anchor]
            if touched is not None and n in touched:
                continue
            refs = []
            for i in range(4):
                r = u.refs[i]
                refs.append(None if r == SU_NO_REF else numbers[r])
            snapshots.append((n, <int> u.kind, <int> par[u.anchor], tuple(refs),
                              <int> (u.spare >> SU_UNNAMED_SHIFT)))
        return snapshots if snapshots else None

    cdef int _replay_parities(self, list snapshots) except -1:
        """Re-base or drop each harvested sign against the rebuilt arena.

        WHERE THE DROP LINE RUNS, and it is the whole reason the harvest snapshots units:

          * a frame that EXISTED AND WAS DESTROYED loses its sign.  The bit is only meaningful
            against the anchor's CSR row and direction count, so once those are gone, keeping it
            would silently re-read the sign against a frame nobody wrote it in.
          * a frame that HAS NOT YET EXISTED keeps its bit, because it is not in `snapshots` at all.
            A parity on an atom that never anchored a unit was never interpreted against anything,
            and a container mid-edit legitimately carries one -- a lone carbon that will be a
            stereocentre once its neighbours arrive.  Reporting and clearing that is
            `validate_stereo`'s business, on a consumer's demand, and never the apply's: a global
            "clear every parity whose atom anchors no unit" sweep would pass every test in this
            file and delete the datum at apply time.
        """
        cdef uint32_t anchor_slot
        cdef uint32_t old_refs[4]
        cdef uint32_t i
        cdef int rc
        cdef bint wrote = False
        cdef bint has_parity = structure_has(self._structure, SEG_PARITY)
        cdef dict index_of = self._index_of
        cdef tuple snap
        cdef object n, ref
        # REALLOCATES THE ARENA (ruling F60): every atom pointer below is taken after this call.
        # UNMARKED (ruling F70): `rebase_parity` reads kind, refs and the unnamed nibble only.
        #
        # THE TABLE HERE IS NECESSARILY ABSENT (ruling F71), and the reason is structural rather than
        # sampled.  `_apply` does not edit an arena, it ALLOCATES A FRESH ONE (`structure_alloc_full`)
        # and rebinds `self._structure` to it, and the only thing that fills derived segments on that
        # path is `rebuild_derived`, which builds six of them and `SEG_STEREO_UNIT` is not among them.
        # So no stereo table can exist at this point, let alone a marked one -- not even by way of
        # `copy()`, which shares an arena but cannot be reached by an apply, because the apply rebinds
        # instead of mutating.  Corroborated rather than established by measurement: with
        # `not structure_has(SEG_STEREO_UNIT)` asserted before this call and `[2] == 0` after it, the
        # whole tree passes (1691 / 1 skipped / 2 xfailed) and neither assertion fires.  That is why this
        # call is the cheap `_unmarked` one, and why nothing below has to invalidate anything.
        ensure_stereo_units_unmarked(self._structure)
        for snap in snapshots:
            n = snap[0]
            if n not in index_of:
                continue            # the anchor was deleted; its bit went with its atom
            anchor_slot = <uint32_t> index_of[n]
            for i in range(4):
                ref = snap[3][i]
                if ref is None:
                    old_refs[i] = SU_NO_REF
                elif ref in index_of:
                    old_refs[i] = <uint32_t> index_of[ref]
                else:
                    # the direction's ATOM is gone.  Not a drop by itself: one direction leaving is
                    # the hydrogen that stopped being drawn, and `rebase_parity` decides.
                    old_refs[i] = RB_GONE
            rc = rebase_parity(self._structure, anchor_slot, <uint8_t> snap[1],
                               <uint8_t> snap[2], old_refs, <uint8_t> snap[4])
            if rc == snap[2]:
                continue            # the permutation was even, or the identity: nothing to write
            wrote = True
            # Only where a segment exists: an arena that states no parity has no segment,
            # and `structure_set_parity` refuses one rather than storing into the shared zero page.
            if has_parity:
                structure_set_parity(self._structure, anchor_slot, 0 if rc == RB_DROP else <uint8_t> rc)
        # THE FEATURE WORDS ARE RE-BASED (ruling F78) and only when something above wrote.  This
        # runs after `rebuild_derived`, so a drop or a flip leaves word IV's bit 6 stating the sign
        # the arena no longer holds -- publicly, through `features_of()` and `_union_feature_words`, on every
        # edit that re-bases or drops a parity.  Same helper as `validate_stereo`'s clear, which is
        # the point: the two writers forgot this independently once already.
        if wrote:
            refresh_parity_features(self._structure)
        # NOTHING IS INVALIDATED HERE (rulings F70, F71 and F74).  The table this replay wrote into is
        # the one the rebuild just derived, and it is unmarked -- measured by assertion over the whole
        # core suite -- so there are no SU_STEREOGENIC marks that could have been computed against the
        # parities just changed.  `structure_invalidate_stereo_units` survives for
        # `validate_stereo`, which clears stereogenicity marks in a CLONE, and `structure_clone` copies
        # derived segments verbatim.
        return 0

    cdef int _apply(self) except -1:
        if self._journal_len == 0:
            return 0

        cdef xy_t *fxy = NULL
        cdef uint8_t *fsg = NULL
        cdef uint8_t *fpar = NULL
        cdef uint32_t sg_var[3]
        cdef int sg_lost
        cdef uint32_t alias_lost
        cdef uint32_t wcount = 0
        cdef uint32_t ccount = 0
        cdef uint32_t n_dropped = 0
        cdef bint cip_stale = False
        cdef bint want_xy = structure_has(self._structure, SEG_XY)
        cdef bint want_sg = structure_has(self._structure, SEG_STEREO_GROUPS)
        cdef uint32_t src_models = structure_conformer_count(self._structure)
        cdef uint32_t conf_adds = 0
        cdef uint32_t conf_drops = 0
        # A PARITY MAY ARRIVE AFTER THE SEAL, which is why the segment is not simply implied by the
        # journal: the SMILES reader's stereo pass reads a perceived frame, so it cannot state its
        # parities until the arena it perceives in exists.  `request_parity` is how it asks.
        cdef bint want_parity = structure_has(self._structure, SEG_PARITY)

        cdef uint32_t n_atoms_old = self._structure.header.atom_count
        cdef uint32_t jn = self._journal_len
        cdef journal_t *jr = self._journal
        cdef uint32_t i, k, wi, wj, n_add = 0, e_add = 0, w_add = 0
        cdef uint8_t jop
        # The anchors the journal states a parity for itself.  Built here rather than in the
        # harvest because this loop is already reading every record, and left None when there is no
        # such record so that the common edit allocates nothing.
        cdef set touched = None
        cdef list snapshots = None
        for i in range(jn):
            jop = jr[i].op
            if jop == OP_ADD_ATOM:
                n_add += 1
                cip_stale = True      # see the CIP invalidation note below
            elif jop == OP_ADD_BOND:
                e_add += 1
                cip_stale = True      # see the CIP invalidation note below
            elif jop == OP_SET_XY:
                want_xy = True
            elif jop == OP_ADD_CONFORMER:
                conf_adds += 1
            elif jop == OP_DROP_CONFORMER:
                conf_drops += 1
            elif jop == OP_SET_WEDGE:
                w_add += 1
            elif jop == OP_SET_STEREO_GROUP:
                want_sg = True
            elif jop == OP_WANT_PARITY:
                want_parity = True
            elif jop == OP_SET_STEREO:
                want_parity = True
                if touched is None:
                    touched = set()
                touched.add(jr[i].a)
            # WHAT INVALIDATES A STORED CIP DESCRIPTOR, decided here, in one place, from the journal.
            #
            # The rule is about the MOLECULE and not about the fields the op happens to write: a
            # descriptor is an assertion the INPUT made about the molecule, so an operation that leaves
            # the molecule the same cannot make the assertion false, and one that changes which atoms
            # exist, which are bonded, or what a bond's order is can.  Ranking at a centre depends on
            # all three, and it depends on them ANYWHERE in the molecule -- so a deletion far from the
            # centre invalidates just as thoroughly as one at it.  Being loudly conservative is the
            # right side to err on: a dropped descriptor is logged and can be recomputed once an
            # assignment algorithm exists, while a carried wrong 'R' is a different molecule to a
            # chemist and nothing downstream can tell.
            #
            # `kekule`/`thiele` ARE EXEMPT, and the exemption is not a special case invented here --
            # RULES and `_kekule.pxi` already name those two as the only operations in the library
            # allowed to change a representation.  They reach this loop as ordinary OP_SET_ORDER
            # records, which is why they cannot be recognised from the journal and must announce
            # themselves instead (`_representation_change`).
            #
            # NOT INVALIDATING, and each is deliberate: charge, isotope, radical, map number, hydrogen
            # count, coordinates, wedges, stereo flags, stereo groups.  Isotope is the interesting one,
            # because CIP Rule 2 ranks by mass and so an isotope edit CAN change a computed descriptor.
            # It is still not dropped, because this layer is not the one that computed it -- and an
            # assignment algorithm that trusted a stored descriptor rather than recomputing would be
            # wrong for a reason no drop rule here could fix.
            #
            # THE ELEMENT IS INVALIDATING AND THE ISOTOPE IS NOT, which is a line worth drawing where
            # it can be read.  Both feed a CIP ranking, but atomic number is Rule 1 -- the primary
            # criterion, ahead of mass -- and an atom whose element changed is not the atom the input
            # made its assertion about.  `set_element` is nearer to swapping one atom for another than
            # to relabelling one, so it is grouped with the ops that change which atoms exist.
            # THE ADD ARMS SET THE FLAG WHERE THEY COUNT, not here: the two adds are matched by earlier
            # arms of the same chain, so a rule stated only here would never run for them.  A dead
            # `elif` in an if/elif chain is not a warning in any language here -- the only thing that
            # catches it is a test per op, which is why there is one.
            elif jop == OP_DELETE_ATOM or jop == OP_DELETE_BOND or jop == OP_SET_ELEMENT:
                cip_stale = True
            elif jop == OP_SET_ORDER:
                if not self._representation_change:
                    cip_stale = True
        cdef uint32_t n_max = n_atoms_old + n_add
        cdef uint32_t e_max = self._structure.header.bond_count + e_add
        cdef uint32_t w_max = 2 * self._structure.header.bond_count + w_add
        # A THIRD COORDINATE IS INDEPENDENT OF THE FIRST TWO, which is D2 of the design and not an
        # oversight: a 3D file fills both segments, so `clean2d()` may replace the depiction without
        # touching the geometry, and a 2D file that never had a z does not grow one here.
        #
        # THE COUNT IS DERIVED, NEVER READ OFF THE SOURCE.  A session that dropped every model must
        # seal with no segment, and the source having one is exactly the state a `structure_has` would
        # read as "keep it" -- so the source count is the starting point and the journal decides.  The
        # first model of a flat molecule is journalled by `set_xyz` itself, so an add is the only thing
        # here that creates one and no case is derived from a coordinate.
        cdef uint32_t conf_models = src_models + conf_adds - conf_drops
        cdef bint want_xyz = conf_models != 0

        # one allocation: 8-aligned regions carved from a single block;
        # wxy and wsg contribute zero bytes when their segment is not wanted.
        # align8() and _scratch_sizes() are the reference — see structure_alloc_full for
        # the same pattern applied to arena segments.
        cdef apply_sizes_t _sz = _scratch_sizes(n_max, e_max, w_max, want_xy, want_sg, conf_models,
                                                want_parity)
        cdef apply_scratch_t scratch
        scratch.block = PyMem_Malloc(_sz.work + _sz.live + _sz.newidx + _sz.edits + _sz.wedge
                                     + _sz.wxy + _sz.wxyz + _sz.wext + _sz.wsg + _sz.wpar
                                     + _sz.cip)
        if scratch.block is NULL:
            self._discard()
            raise MemoryError('journal apply scratch allocation failed')
        # Each pointer is paired with its own size field by hand below.  That pairing is the one
        # property _apply_scratch_probe cannot check -- it calls _scratch_sizes too, so a swap
        # here (e.g. _bp += _sz.wedge after scratch.edits) would leave the probe unchanged.
        # An edit to this carve needs re-derivation by hand.
        cdef char *_bp = <char *> scratch.block
        scratch.work   = <atom_t *>       _bp;  _bp += _sz.work
        scratch.live   = <uint8_t *>      _bp;  _bp += _sz.live
        scratch.newidx = <int32_t *>      _bp;  _bp += _sz.newidx
        scratch.edits  = <edge_edit_t *>  _bp;  _bp += _sz.edits
        scratch.wedge  = <wedge_edit_t *> _bp;  _bp += _sz.wedge
        if want_xy:
            scratch.wxy = <xy_t *> _bp
        else:
            scratch.wxy = NULL
        _bp += _sz.wxy
        if want_xyz:
            scratch.wxyz = <xyz_t *> _bp
        else:
            scratch.wxyz = NULL
        _bp += _sz.wxyz
        if want_xyz:
            scratch.wext = <uint32_t *> _bp
        else:
            scratch.wext = NULL
        _bp += _sz.wext
        if want_sg:
            scratch.wsg = <uint8_t *> _bp
        else:
            scratch.wsg = NULL
        _bp += _sz.wsg
        if want_parity:
            scratch.wpar = <uint8_t *> _bp
        else:
            scratch.wpar = NULL
        _bp += _sz.wpar
        scratch.cip = <cip_edit_t *> _bp

        cdef atom_t *work = scratch.work
        cdef uint8_t *live = scratch.live
        cdef int32_t *newidx = scratch.newidx
        cdef edge_edit_t *edits = scratch.edits
        cdef wedge_edit_t *wedge = scratch.wedge
        cdef xy_t *wxy = scratch.wxy
        cdef xyz_t *wxyz = scratch.wxyz
        cdef uint32_t *wext = scratch.wext
        cdef uint8_t *wsg = scratch.wsg
        cdef uint8_t *wpar = scratch.wpar
        cdef cip_edit_t *cip = scratch.cip
        cdef uint32_t ecount = 0, e_new = 0, n_new = 0
        cdef uint32_t seed_wcount = 0
        cdef uint32_t *ptr
        cdef halfedge_t *edges
        cdef journal_t *rec
        cdef halfedge_t e
        cdef halfedge_t *he
        cdef atom_t *fresh_atoms
        cdef atom_t *wa
        cdef edge_edit_t *ed
        cdef edge_edit_t *enew
        cdef wedge_edit_t *wg
        cdef cip_edit_t *cg
        cdef xy_t *wxy_i
        cdef xyz_t *wxyz_i
        cdef xyz_t *fxyz = NULL
        cdef conformer_t *frec = NULL
        cdef conformer_t *src_rec
        cdef uint32_t model
        cdef int32_t ni, ni_src, ni_dst
        cdef Structure fresh
        cdef list numbers
        cdef dict index_of
        cdef uint8_t op
        cdef bint found
        cdef int rc
        cdef uint32_t d
        try:
            # FIRST, and before any pointer into the old arena is taken: the harvest perceives the
            # OLD molecule's units, which reallocates it (ruling F60).  Inside the `try` so that a
            # failure here still discards the journal, like every other failure in the apply.
            snapshots = self._harvest_parities(touched)
            if n_atoms_old:
                memcpy(work, self._structure.atoms(), n_atoms_old * sizeof(atom_t))
            if n_add:
                memset(<void *> (work + n_atoms_old), 0, n_add * sizeof(atom_t))
            if want_xy:
                if n_atoms_old:
                    if structure_has(self._structure, SEG_XY):
                        memcpy(wxy, structure_xy(self._structure), n_atoms_old * sizeof(xy_t))
                    else:
                        memset(<void *> wxy, 0, n_atoms_old * sizeof(xy_t))
                if n_add:
                    memset(<void *> (wxy + n_atoms_old), 0, n_add * sizeof(xy_t))
            if want_xyz:
                # MEMSET FIRST AND COPY OVER IT, which collapses two origin cases into one: an atom
                # added this session (design D7 -- there is no honest z for an atom the caller placed
                # by connectivity alone, and a per-atom validity bitmap was declined, so `xyz_of`
                # answers (0, 0, 0) rather than None) and a model added this session, whose every
                # atom is at the origin for the same reason.  Neither needs an arm of its own.
                memset(<void *> wxyz, 0, <size_t> conf_models * n_max * sizeof(xyz_t))
                for model in range(conf_models):
                    wext[model] = <uint32_t> CONF_NO_INDEX
                # THE SURVIVING SOURCE MODELS, COMPACTED.  `d` is the destination cursor, so a dropped
                # model is skipped rather than blanked and everything above it moves down one -- which
                # is what makes an index a list position and not a stable name.
                if src_models:
                    src_rec = structure_conformer_records(self._structure)
                    d = 0
                    for model in range(src_models):
                        if self._conf_dropped is not None and model in self._conf_dropped:
                            continue
                        if n_atoms_old:
                            memcpy(wxyz + <size_t> d * n_max,
                                   structure_conformer_xyz(self._structure, model),
                                   n_atoms_old * sizeof(xyz_t))
                        wext[d] = src_rec[model].ext_index
                        d += 1
            if want_sg:
                if n_atoms_old and structure_has(self._structure, SEG_STEREO_GROUPS):
                    memcpy(wsg, structure_stereo_groups(self._structure), n_atoms_old)
                elif n_atoms_old:
                    memset(<void *> wsg, 0, n_atoms_old)
                if n_add:
                    memset(<void *> (wsg + n_atoms_old), 0, n_add)
            if want_parity:
                # SEEDED FROM THE OLD ARENA'S SEGMENT when it has one, and otherwise zero because
                # an arena with no segment states no parity.
                if n_atoms_old and structure_has(self._structure, SEG_PARITY):
                    memcpy(wpar, structure_parities(self._structure), n_atoms_old)
                elif n_atoms_old:
                    memset(<void *> wpar, 0, n_atoms_old)
                if n_add:
                    memset(<void *> (wpar + n_atoms_old), 0, n_add)
            memset(live, 1, n_max)

            # existing bonds, canonical half-edges only; seed wedges from all half-edges
            ptr = csr_ptr(self._structure)
            edges = csr_edges(self._structure)
            for i in range(n_atoms_old):
                for k in range(ptr[i], ptr[i + 1]):
                    e = edges[k]
                    if e.to > i:
                        ed = &edits[ecount]
                        ed.src = i
                        ed.dst = e.to
                        ed.order = e.order
                        ecount += 1
                        # CANONICAL HALF ONLY, unlike the wedge below, because a wedge is directional
                        # and a descriptor is not.  Seeding both halves would write the same bond
                        # twice and make `ccount` exceed the `e_max` this region is sized by.
                        if he_cip(&edges[k]):
                            cg = &cip[ccount]
                            cg.src = i
                            cg.dst = e.to
                            cg.code = he_cip(&edges[k])
                            ccount += 1
                    if e.wedge:
                        wg = &wedge[wcount]
                        wg.src = i
                        wg.dst = e.to
                        wg.code = e.wedge
                        wcount += 1
            seed_wcount = wcount
            # THE DROP HAPPENS HERE, AFTER THE SEED AND BEFORE THE REPLAY, so that a descriptor set in
            # the SAME scope as the invalidating edit still wins.  `with mol.edit(): delete_atom(x);
            # set_bond_cip(a, b, 'E')` is a caller stating a descriptor about the molecule it just
            # made, and the journal's order is what says so.  Clearing after the replay instead would
            # throw that away, and clearing before the seed would leave the old values to be re-seeded.
            if cip_stale:
                if ccount:
                    self._log_event('container:cip-dropped', 'edit:cip',
                                    '%d bond CIP descriptor(s) dropped: the molecule changed' % ccount,
                                    mc_lost())
                    ccount = 0
                n_dropped = 0
                for i in range(n_atoms_old):
                    if at_cip(&work[i]):
                        at_set_cip(&work[i], 0)
                        n_dropped += 1
                if n_dropped:
                    self._log_event('container:cip-dropped', 'edit:cip',
                                    '%d atom CIP descriptor(s) dropped: the molecule changed'
                                    % n_dropped, mc_lost())

            # replay in emission order, so a later record simply overwrites an earlier one
            for i in range(jn):
                rec = jr + i
                op = rec.op
                if op < OP_ADD_ATOM or op > OP_HIGHEST:
                    raise NotImplementedError('journal op %d has no apply rule yet' % op)
                if op == OP_WANT_PARITY or op == OP_DROP_CONFORMER:
                    # NAMES NO ATOM AND WRITES NO FIELD HERE: the first was read in the pre-pass, the
                    # second by the seed, whose compaction is the whole of its effect.
                    continue
                if op == OP_ADD_CONFORMER:
                    # NAMES NO ATOM EITHER, but it does write a field: the stated number, at the
                    # destination the op was given.  The coordinates are the memset's already.
                    wext[rec.model] = rec.a
                    continue
                wi = self._work_index(rec.a, n_atoms_old)
                wa = &work[wi]
                if op == OP_ADD_ATOM:
                    wa.element = <uint8_t> rec.v
                    wa.n = rec.a
                    # A NEW ATOM'S HYDROGEN COUNT IS UNKNOWN, NOT ZERO, and this line is the whole
                    # of that rule.  The slot arrives zeroed from the memset above, and a zero
                    # nibble is a STATEMENT -- "this atom has no hydrogens" -- which the builder is
                    # in no position to make on the caller's behalf.  `kekule()` reading a stored 0
                    # on a two-coordinate aromatic N pins pyrrole to zero hydrogens and makes it
                    # unkekulisable; reading H_UNKNOWN it correctly treats the count as unstated
                    # and chooses.  An OP_SET_HYDROGENS emitted later in the same journal
                    # overwrites this, so `add_atom(implicit_h=0)` still stores a real zero.
                    at_set_h(wa, H_UNKNOWN, 0)
                elif op == OP_DELETE_ATOM:
                    live[wi] = 0
                elif op == OP_ADD_BOND:
                    ed = &edits[ecount]
                    ed.src = wi
                    ed.dst = self._work_index(rec.b, n_atoms_old)
                    ed.order = <uint8_t> rec.v
                    ecount += 1
                elif op == OP_DELETE_BOND or op == OP_SET_ORDER:
                    wj = self._work_index(rec.b, n_atoms_old)
                    found = False
                    for k in range(ecount):
                        ed = &edits[k]
                        if ed.order and ((ed.src == wi and ed.dst == wj)
                                         or (ed.src == wj and ed.dst == wi)):
                            ed.order = 0 if op == OP_DELETE_BOND else <uint8_t> rec.v
                            found = True
                            break
                    if not found:
                        raise KeyError((rec.a, rec.b))
                elif op == OP_SET_ELEMENT:
                    # The HYDROGEN COUNT IS LEFT ALONE, and that is the whole of this arm's
                    # subtlety: a count derived for the old element is almost certainly wrong for
                    # the new one, but this layer does not derive counts -- `calc_implicit` does,
                    # and it lives in `chemistry`.  Writing H_UNKNOWN here would throw away a count
                    # the caller may be about to write correctly; guessing one is what the core
                    # refuses everywhere else.  `set_element`'s docstring says so out loud.
                    wa.element = <uint8_t> rec.v
                elif op == OP_SET_CHARGE:
                    wa.charge = <int8_t> rec.v
                elif op == OP_SET_ISOTOPE:
                    wa.isotope = <uint16_t> rec.v
                elif op == OP_SET_MAP_NUMBER:
                    wa.map_number = <uint16_t> rec.v
                elif op == OP_SET_RADICAL:
                    at_set_radical(wa, rec.v != 0)
                elif op == OP_SET_STEREO:
                    # `want_parity` is guaranteed true here: the pre-pass set it from this very record.
                    wpar[wi] = <uint8_t> rec.v
                elif op == OP_SET_HYDROGENS:
                    at_set_h(wa, <uint8_t> rec.v, at_explicit_h(wa))
                    at_set_h_pinned(wa, True)
                elif op == OP_SET_XY:
                    wxy_i = &wxy[wi]
                    wxy_i.x = <int32_t> rec.b
                    wxy_i.y = rec.v
                elif op == OP_SET_XYZ:
                    # The one op that spends `rec.w` -- see the note on `journal_t`.
                    wxyz_i = wxyz + <size_t> rec.model * n_max + wi
                    wxyz_i.x = <int32_t> rec.b
                    wxyz_i.y = rec.v
                    wxyz_i.z = rec.w
                elif op == OP_SET_WEDGE:
                    wg = &wedge[wcount]
                    wg.src = wi
                    wg.dst = self._work_index(rec.b, n_atoms_old)
                    wg.code = <uint8_t> rec.v
                    wcount += 1
                    # A wedge is geometry only (Ruling F54); a parity is not written here.
                    # `wedge_of()` / `wedges()` read the wedge segment directly.  The parity
                    # will be derived from (wedge code, coordinates, reference order) by the
                    # wedge-ingestion work.
                elif op == OP_SET_STEREO_GROUP:
                    wsg[wi] = sg_pack(<uint8_t> rec.v, <uint8_t> rec.b)
                elif op == OP_SET_ATOM_CIP:
                    # `wi` came from `_work_index`, which raises for an atom the journal has not added
                    # yet -- so a descriptor op that precedes its atom's ADD_ATOM fails here rather
                    # than landing on whatever slot that index happens to name.
                    at_set_cip(wa, <uint8_t> rec.v)
                elif op == OP_SET_BOND_CIP:
                    # Overwrite in place if this bond already has an entry, so a scope that sets the
                    # same bond twice does not spend two slots -- `ccount` is bounded by `e_max`, and
                    # a bond set twice would otherwise overrun that bound on a molecule where every
                    # bond is set once already.
                    wj = self._work_index(rec.b, n_atoms_old)
                    found = False
                    for k in range(ccount):
                        cg = &cip[k]
                        if (cg.src == wi and cg.dst == wj) or (cg.src == wj and cg.dst == wi):
                            cg.code = <uint8_t> rec.v
                            found = True
                            break
                    if not found:
                        cg = &cip[ccount]
                        cg.src = wi
                        cg.dst = wj
                        cg.code = <uint8_t> rec.v
                        ccount += 1
                elif op == OP_SET_R_INDEX:
                    # THE CROSS-FIELD RULE, AT ITS ONE ENFORCEMENT POINT.  Element is final here and
                    # not at queue time: `_atom` refuses a read while the journal is dirty, so
                    # `add_atom('C')` followed by `set_r_index` in one scope can only be caught on
                    # replay.
                    if wa.element != 0:
                        raise ValueError('atom %d is not an R; the R index is only meaningful on '
                                         'element 0' % rec.a)
                    at_set_r_index(wa, <uint8_t> rec.v)
                else:
                    raise AssertionError('journal op %d is in range but has no replay arm' % op)

            # compact the survivors
            for i in range(n_max):
                if live[i]:
                    newidx[i] = <int32_t> n_new
                    n_new += 1
                else:
                    newidx[i] = -1
            for k in range(ecount):
                ed = &edits[k]
                if not ed.order:
                    continue                       # a deleted bond, marked by order 0
                ni_src = newidx[ed.src]
                ni_dst = newidx[ed.dst]
                if ni_src < 0 or ni_dst < 0:
                    continue                       # an endpoint was deleted
                # Order 4 passes through UNCHANGED, because input fidelity is the invariant: a caller
                # who states an aromatic bond gets an aromatic bond stored, and `kekule()` is the
                # explicit operation that converts one. `_emit_half` sets HE_AROMATIC from the order,
                # so the flag and the order cannot be written apart.
                # e_new <= k, so enew may alias ed on the pass-through case.  Both index
                # values were read out above, and `order` is copied last, so the aliasing
                # write order is harmless.
                enew = &edits[e_new]
                enew.src = <uint32_t> ni_src
                enew.dst = <uint32_t> ni_dst
                enew.order = ed.order
                e_new += 1

            # S-GROUPS ARE SIZED AT THE OLD MOLECULE'S LENGTHS, AND THAT IS ALWAYS ENOUGH.  The record
            # count never changes (an emptied record stays present and empty), the blob is copied byte
            # for byte, and only the index run can shrink -- so `structure_sgroup_var_len` reads three
            # lengths off the source and the carry compacts into them.  The slack that leaves at the end
            # of the index segment is zeroed and unreferenced; it is not worth a second pass to reclaim,
            # because reclaiming it would mean sizing before the carry has decided what survives.
            if structure_has(self._structure, SEG_OPAQUE_BLOB):
                structure_sgroup_var_len(self._structure, sg_var)
                fresh = structure_alloc_full(n_new, e_new, n_new > 65535,
                                             (SEG_MASK_XY if want_xy else 0)
                                             | (SEG_MASK_STEREO if want_sg else 0)
                                             | (SEG_MASK_PARITY if want_parity else 0), sg_var,
                                             conf_models)
            else:
                fresh = structure_alloc_full(n_new, e_new, n_new > 65535,
                                             (SEG_MASK_XY if want_xy else 0)
                                             | (SEG_MASK_STEREO if want_sg else 0)
                                             | (SEG_MASK_PARITY if want_parity else 0), NULL,
                                             conf_models)
            if structure_has(self._structure, SEG_OPAQUE_BLOB):
                # The blob first and verbatim: it is the OPAQUE half by definition, and no edit to the
                # chemistry can invalidate a byte of it.
                memcpy(structure_blob(fresh), structure_blob(self._structure), sg_var[2])
                alias_lost = 0
                sg_lost = structure_carry_sgroups(fresh, self._structure, newidx, &alias_lost)
                # TWO LINES AND NOT ONE.  An alias is stored as a record like any other, but a writer
                # emits it as a display label on a single atom while an S-group gets its own block, so a
                # single count would send a reader to the wrong half of the file.
                #
                # ALIASES FIRST, and the order is part of the contract rather than an accident of which
                # `if` came first.  A V2000 record puts the `A`/`V` alias lines in the atom-adjacent part
                # and the `M  ST*` properties block after them, so this order is the one a reader
                # scanning the file already holds -- diffing this log against a file walks both in the
                # same direction.  Measured and requested by the format epic, which builds against it.
                if alias_lost:
                    self._log_event('container:alias-lost', 'edit:sgroup',
                                    '%d alias(es) lost the atom they label' % alias_lost, mc_lost())
                if sg_lost:
                    self._log_event('container:sgroup-lost', 'edit:sgroup',
                                    '%d sgroup record(s) lost a reference to a deleted atom' % sg_lost,
                                    mc_lost())
            fresh_atoms = fresh.atoms()
            if want_xy:
                fxy = structure_xy(fresh)
            if want_xyz:
                fxyz = structure_conformer_xyz(fresh, 0)
                frec = structure_conformer_records(fresh)
                for model in range(conf_models):
                    frec[model].ext_index = wext[model]
            if want_sg:
                fsg = structure_stereo_groups(fresh)
            if want_parity:
                fpar = structure_parities(fresh)
            numbers = []
            for i in range(n_max):
                if live[i]:
                    ni = newidx[i]
                    fresh_atoms[ni] = work[i]
                    if want_xy:
                        fxy[ni] = wxy[i]
                    if want_xyz:
                        # A DELETED ATOM DROPS ITS COLUMN IN EVERY MODEL, in the same pass and by the
                        # same `newidx` that moves the atom record.  That is the whole of design D7's
                        # first half: geometry is per atom, so it follows the atom's slot and cannot
                        # desynchronise from it.  A bond edit reaches this loop with every `live[i]`
                        # set and every `newidx[i] == i`, so it copies the geometry through unchanged.
                        # The models are contiguous, so `structure_conformer_xyz(fresh, m)` is
                        # `fxyz + m * n_new` and one base pointer indexes them all.
                        for model in range(conf_models):
                            fxyz[<size_t> model * n_new + ni] = wxyz[<size_t> model * n_max + i]
                    if want_sg:
                        fsg[ni] = wsg[i]
                    if want_parity:
                        fpar[ni] = wpar[i]
                    numbers.append(work[i].n)
            with nogil:
                rc = csr_build(fresh, edits, e_new)
            if rc:
                raise MemoryError('csr scratch allocation failed')
            ptr = csr_ptr(fresh)
            edges = csr_edges(fresh)
            for i in range(n_new):
                d = ptr[i + 1] - ptr[i]
                fresh_atoms[i].degree = <uint8_t> (d if d < 255 else 255)
                for k in range(ptr[i] + 1, ptr[i + 1]):
                    e = edges[k]
                    if e.to == edges[k - 1].to:
                        raise ValueError(f'duplicate bond between stable ids {numbers[i]} '
                                         f'and {numbers[e.to]}')

            for k in range(wcount):
                wg = &wedge[k]
                ni_src = newidx[wg.src]
                ni_dst = newidx[wg.dst]
                if ni_src < 0 or ni_dst < 0:
                    continue          # an endpoint was deleted; its wedge goes with it
                he = csr_find(fresh, <uint32_t> ni_src, <uint32_t> ni_dst)
                if he is NULL:
                    if k < seed_wcount:
                        continue      # bond deleted after the wedge was set; the wedge goes with it
                    raise KeyError((numbers[ni_src], numbers[ni_dst]))
                he.wedge = wg.code

            # BOTH HALVES, FROM ONE CALL SITE.  This is the only place in the library that writes a
            # bond descriptor, which is what makes "a descriptor is not direction-dependent" a
            # property of the code rather than a rule someone has to remember.  `he_set_cip` takes a
            # single half-edge deliberately, so that no caller can reach for a one-sided shortcut.
            for k in range(ccount):
                cg = &cip[k]
                ni_src = newidx[cg.src]
                ni_dst = newidx[cg.dst]
                if ni_src < 0 or ni_dst < 0:
                    continue          # unreachable while `cip_stale` covers deletion; cheap insurance
                he = csr_find(fresh, <uint32_t> ni_src, <uint32_t> ni_dst)
                if he is NULL:
                    raise KeyError((numbers[ni_src], numbers[ni_dst]))
                he_set_cip(he, cg.code)
                he = csr_find(fresh, <uint32_t> ni_dst, <uint32_t> ni_src)
                he_set_cip(he, cg.code)

            rebuild_derived(fresh)

            index_of = {}
            for i in range(n_new):
                index_of[numbers[i]] = i
            self._structure = fresh
            self._numbers = numbers
            self._index_of = index_of
            self._gen += 1
        finally:
            PyMem_Free(scratch.block)
            # on both paths: a journal is applied exactly once, and a failed apply is dropped
            # rather than left to re-raise on every later read
            self._discard()
        # AFTER the finally, deliberately: the replay perceives the NEW molecule, which reads
        # through the container's own accessors and so needs the journal already discarded and the
        # arena already committed.  The bits it writes were carried through the memcpy unchanged, so
        # a re-base is a rewrite of a value that is already in place -- there is no window in which
        # a parity is missing, only one in which it is not yet re-based.
        if snapshots is not None:
            self._replay_parities(snapshots)
        return 0

    cdef bint _has_bond(self, uint32_t n, uint32_t m) except -1:
        if n not in self._index_of or m not in self._index_of:
            return False
        return csr_find(self._structure, <uint32_t> self._index_of[n],
                        <uint32_t> self._index_of[m]) is not NULL

    cdef inline atom_t *_atom(self, uint32_t n) except NULL:
        self._require_clean()
        return self._structure.atoms() + <uint32_t> self._index_of[n]

    cdef inline uint32_t _slot(self, uint32_t n) except *:
        """The arena slot of stable id `n`.  `_atom`'s sibling, for the segments that are indexed by
        slot rather than reached through an `atom_t *`."""
        self._require_clean()
        return <uint32_t> self._index_of[n]

    @property
    def journal_length(self):
        return self._journal_len

    @property
    def generation(self):
        return self._gen

    def journal_record(self, uint32_t i):
        if i >= self._journal_len:
            raise IndexError(i)
        cdef journal_t *rec = self._journal + i
        # SIX FIELDS AND NOT FOUR, since `w` and then `model` landed.  A probe that hid a payload word
        # would let a wrong `w` through in exactly the tests written to catch one, so every op's
        # record is reported whole and the ops that do not spend a field report the zero they carry.
        # `model` goes LAST rather than beside `op`, because a tuple index is what the tests read:
        # `[3]` is `v` in both spellings.
        return (rec.op, rec.a, rec.b, rec.v, rec.w, rec.model)

    @property
    def atom_count(self):
        self._require_clean()
        return self._structure.header.atom_count

    @property
    def bond_count(self):
        self._require_clean()
        return self._structure.header.bond_count

    @property
    def unknown_h_count(self):
        """How many atoms carry NO implicit hydrogen count -- the sentinel, not a zero.

        Zero on every molecule built from a record that stated its hydrogens, so `if
        mol.unknown_h_count:` is the one test a caller needs before trusting anything derived from
        hydrogen counts: `float(mol)` (a lower bound while this is non-zero), a formula, a valence
        check, any `h`/`H` query primitive (which cannot match such an atom in either direction).

        A COUNT AND NOT A LIST, on purpose.  The question a caller actually has is "is this record
        complete"; the atoms themselves are reachable with `implicit_h_of` returning None, and a
        property that built a list would allocate one on every molecule to answer a question that
        is almost always "none".

        AT AN OUTPUT BOUNDARY THIS IS A LOSS TO REPORT, and how big a loss depends on whether the
        format can spell "no count".  Two cases, and a writer must know which one it is in:

        * THE FORMAT RESERVES A VALUE FOR "NOT STATED", so omission is honest.  CTfile does: V2000's
          `vvv` valence field reads 0 as "default" (and 15 as ZERO valence, not as fifteen), and
          V3000's `VAL=` reads 0 the same way.  The loss is still real but it is smaller -- what
          goes missing is that a downstream reader will re-derive the count and may land somewhere
          else -- so the writer omits the field AND logs it.
        * THE FORMAT'S OMISSION MEANS ZERO, so omission is a false statement.  An absent hydrogen
          term inside SMILES brackets means ZERO in OpenSMILES, not "unspecified": a writer that
          drops the term has stated a number the molecule never had.  Omission is right there only
          where the language supplies the count back, which is a bare `C` outside brackets and
          nothing else.

        Never a silent zero in either case.  Same rule as a valence table with no row for a charge:
        the third state is worth having only while it stays distinguishable from a real answer.
        """
        self._require_clean()
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t i
        cdef uint32_t acc = 0
        for i in range(self._structure.header.atom_count):
            if at_implicit_h_unknown(&atoms[i]):
                acc += 1
        return acc

    @property
    def atom_numbers(self):
        # The atom numbers in arena order. An atom number IS its stable id and survives every edit,
        # which is why atom-atom mapping lives in the separate `map_number` field.
        self._require_clean()
        return list(self._numbers)

    def index_of(self, uint32_t n):
        self._require_clean()
        return self._index_of[n]

    def number_of(self, uint32_t index):
        """The atom number at arena position `index`; the inverse of `index_of`."""
        self._require_clean()
        if index >= <uint32_t> len(self._numbers):
            raise IndexError(index)
        return self._numbers[index]

    def element_of(self, uint32_t n):
        return self._atom(n).element

    def radius_of(self, uint32_t n):
        """The calculated atomic radius of atom `n` in angstroms, and 0.0 for the R marker.

        Element data rather than per-atom state, so it answers the same for every atom of an element;
        `el_atomic_radius` states which radius and where the published set stops.
        """
        return el_atomic_radius(self._atom(n).element)

    def charge_of(self, uint32_t n):
        return self._atom(n).charge

    def isotope_of(self, uint32_t n):
        return self._atom(n).isotope

    def map_number_of(self, uint32_t n):
        return self._atom(n).map_number

    def degree_of(self, uint32_t n):
        return self._atom(n).degree

    def implicit_h_of(self, uint32_t n):
        """Implicit hydrogens on atom `n`, or None when the record does not say.

        None and not 15: the sentinel is a STORAGE spelling and does not belong on this surface,
        where every caller would have to know the number to avoid adding it to something -- see
        H_UNKNOWN in _molecule_arena.pxi.

        None is not zero.  A record that omits a hydrogen count and a record that states none are
        different records, and only the second one is a methane carbon.
        """
        cdef atom_t *a = self._atom(n)
        if at_implicit_h_unknown(a):
            return None
        return at_implicit_h(a)

    def explicit_h_of(self, uint32_t n):
        return at_explicit_h(self._atom(n))

    def heteroatoms_of(self, uint32_t n):
        return self._atom(n).heteroatoms

    def hybridization_of(self, uint32_t n):
        return at_hybridization(self._atom(n))

    def total_h_of(self, uint32_t n):
        """Implicit plus explicit hydrogens, or None when the implicit count is unknown.

        A sum with an unknown term is unknown, so this refuses rather than reporting the explicit
        count alone -- which would read as a total and be wrong by however many the record left out.
        `explicit_h_of` is always a number and is the way to ask for the part that IS known.
        """
        cdef atom_t *a = self._atom(n)
        if at_implicit_h_unknown(a):
            return None
        return at_implicit_h(a) + at_explicit_h(a)

    def radical_of(self, uint32_t n):
        return at_radical(self._atom(n))

    def stereo_of(self, uint32_t n):
        """True when atom `n`'s configured parity is odd.  `parity_of` is the three-state read."""
        return structure_parity_at(self._structure, self._slot(n)) == 2

    def parity_of(self, uint32_t n):
        """Three-state parity of atom `n`: 0 = no parity configured, 1 = even, 2 = odd.

        A wedge does NOT configure one (Ruling F54): `set_wedge` writes no parity, so a wedge-drawn
        centre reads 0 exactly as an undrawn one does.  Parity is derived from (wedge code, coordinates,
        reference order); `wedge_of` answers whether a wedge exists.
        """
        return structure_parity_at(self._structure, self._slot(n))

    def in_ring_of(self, uint32_t n):
        return at_in_ring(self._atom(n))

    def bond_in_ring(self, uint32_t n, uint32_t m):
        self._require_clean()
        cdef halfedge_t *e = csr_find(self._structure, <uint32_t> self._index_of[n],
                                      <uint32_t> self._index_of[m])
        if e is NULL:
            raise KeyError((n, m))
        return (e.flags & HE_IN_RING) != 0

    def order_of(self, uint32_t n, uint32_t m):
        self._require_clean()
        cdef halfedge_t *e = csr_find(self._structure, <uint32_t> self._index_of[n],
                                      <uint32_t> self._index_of[m])
        return None if e is NULL else e.order

    def neighbors_of(self, uint32_t n):
        self._require_clean()
        if n not in self._index_of:
            raise KeyError(n)
        cdef uint32_t i = <uint32_t> self._index_of[n]
        cdef uint32_t *ptr = csr_ptr(self._structure)
        cdef halfedge_t *edges = csr_edges(self._structure)
        cdef list numbers = self._numbers
        cdef uint32_t k
        cdef list out = []
        for k in range(ptr[i], ptr[i + 1]):
            out.append(numbers[edges[k].to])
        return out

    def edge_words_of(self, uint32_t n):
        """The isomorphism edge word of each half-edge out of this atom, in CSR order."""
        self._require_clean()
        cdef uint32_t i = self._index_of[n]
        cdef uint32_t *ptr = csr_ptr(self._structure)
        cdef uint64_t *words = structure_edge_words(self._structure)
        cdef uint32_t k
        cdef list out = []
        for k in range(ptr[i], ptr[i + 1]):
            out.append(words[k])
        return out

    def component_labels(self):
        """{stable id: 0-based connected-component label}."""
        self._require_clean()
        ensure_component_labels(self._structure)
        cdef uint32_t *label = structure_component_labels(self._structure)
        cdef list numbers = self._numbers
        cdef uint32_t i
        cdef dict out = {}
        for i in range(self._structure.header.atom_count):
            out[numbers[i]] = label[i]
        return out

    def edit(self):
        return _EditScope(self)

    def __enter__(self):
        # `with mol:` is `with mol.edit():` -- same counter, so the two nest in either order.
        # An unapplied journal never reached the arena in the first place, so __exit__ drops it
        # and there is nothing to undo.
        self._scope_depth += 1
        return self

    @cython.warn.unused_arg(False)
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._scope_depth -= 1
        if self._scope_depth == 0:
            if exc_type is None:
                self._apply()
            else:
                self._discard()
        return False

    def remap(self, dict mapping not None):
        """
        Relabel stable ids in place. `mapping` is {old id: new id}; ids absent from it keep theirs.

        Atom order, bonds and every derived descriptor are untouched -- only the labels move, so
        this costs one buffer copy and no re-perception. For a relabelled copy: mol.copy().remap().
        """
        self._require_clean()
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef uint32_t i, old, new, high = 0
        cdef object key, value
        cdef list numbers = []
        cdef dict index_of = {}
        for key in mapping:
            if key not in self._index_of:
                raise KeyError(key)
        # resolve the whole relabelling before touching anything: a rejected mapping must
        # leave the container exactly as it was
        for i in range(n_atoms):
            old = <uint32_t> self._numbers[i]
            value = mapping.get(old, old)
            if not isinstance(value, int) or isinstance(value, bool):
                raise TypeError(f'stable id must be an int, got {value!r}')
            if value < 1 or value > 0xFFFFFFFE:
                raise ValueError(f'stable id {value!r} is out of range 1..4294967294')
            new = <uint32_t> value
            if new in index_of:
                raise ValueError(f'remap is not injective: two atoms would both get id {new}')
            index_of[new] = i
            numbers.append(new)
            if new > high:
                high = new

        cdef Structure fresh = structure_clone(self._structure)
        cdef atom_t *atoms = fresh.atoms()
        for i in range(n_atoms):
            atoms[i].n = <uint32_t> numbers[i]

        self._structure = fresh
        self._numbers = numbers
        self._index_of = index_of
        if high >= self._next_id:
            self._next_id = high + 1
        self._first_pending = self._next_id
        self._gen += 1

    def copy(self):
        # O(1): the arena is immutable, and _apply rebinds _numbers/_index_of rather than
        # mutating them, so all three are safe to share with a second container
        self._require_clean()
        cdef MoleculeContainer mol = MoleculeContainer.__new__(MoleculeContainer)
        mol._structure = self._structure
        mol._numbers = self._numbers
        mol._index_of = self._index_of
        mol._next_id = self._next_id
        mol._first_pending = self._next_id
        mol._gen = self._gen
        # The canonical form travels with the copy: it is a function of the arena, the arena is
        # immutable and shared, so recomputing it would buy nothing and cost the whole
        # canonicalisation.  `_order_cache` is left behind only because it is keyed by stable id and
        # cheap; this one is neither.
        #
        # ONLY WHEN IT IS CURRENTLY VALID, and the guard is the whole point rather than a belt.  A
        # cache row here is a PAIR -- bytes plus the generation they belong to -- and copying the
        # bytes while stamping them with `self._gen` re-validates a row the source itself would have
        # rejected.  Every edit bumps `_gen` and leaves `_identity_gen` behind, so
        # `m.canonical_bytes` / edit / `m.copy()` would hand the copy the PRE-EDIT molecule's identity,
        # sworn to be the post-edit arena's, and `==`, `hash()` and `str()` would answer for a molecule
        # that does not exist.  `split()` inherits this path through `copy()`, so a patched product
        # would report the identity of the input it came from.
        if self._identity_cache is not None and self._identity_gen == self._gen:
            mol._identity_cache = self._identity_cache
            mol._identity_gen = self._gen
        # and the canonical SMILES for the same reason, word for word: it is a function of the arena,
        # the arena is shared and immutable, and it is the same canonicalisation being paid for
        if self._smiles_cache is not None and self._smiles_gen == self._gen:
            mol._smiles_cache = self._smiles_cache
            mol._smiles_gen = self._gen
        # Shallow: the values are the caller's objects and copying them would be a second policy
        # nobody asked for.
        mol._meta = None if self._meta is None else dict(self._meta)
        return mol

    def shares_arena_with(self, MoleculeContainer other not None):
        return self._structure is other._structure

    def atom(self, uint32_t n):
        self._require_clean()
        if n not in self._index_of:
            raise KeyError(n)
        cdef Atom a = Atom.__new__(Atom)
        a._molecule = self
        a._n = n
        a._gen = self._gen
        return a

    def bond(self, uint32_t n, uint32_t m):
        self._require_clean()
        if not self._has_bond(n, m):
            raise KeyError((n, m))
        cdef Bond bd = Bond.__new__(Bond)
        bd._molecule = self
        bd._n = n
        bd._m = m
        bd._gen = self._gen
        return bd

    def conformer(self, uint32_t index):
        """Model `index` as a `Conformer` view; IndexError when there is no such model.

        IndexError and not KeyError, where `atom()` raises KeyError: an index is a list position, so
        the exception a sequence raises is the one a caller expects.
        """
        self._require_clean()
        if index >= structure_conformer_count(self._structure):
            raise IndexError(index)
        cdef Conformer c = Conformer.__new__(Conformer)
        c._molecule = self
        c._index = index
        c._gen = self._gen
        return c

    @property
    def conformers(self):
        """Every model, as a tuple of `Conformer` views.

        A tuple and not a lazy view: it already is a read-only sequence with `len`, iteration and
        indexing, the count is small, and `edit()` is the one way to write a model.  Empty for a
        molecule with no geometry.
        """
        self._require_clean()
        cdef uint32_t i
        cdef list out = []
        for i in range(structure_conformer_count(self._structure)):
            out.append(self.conformer(i))
        return tuple(out)

    def atoms(self):
        self._require_clean()
        cdef uint32_t n
        for n in self._numbers:
            yield self.atom(n)

    def bonds(self):
        # Yields each bond once, as a Bond view; its endpoints are bond.n and bond.m, and `n` is the
        # EARLIER atom -- the `to > i` test below is what makes that true, and a caller building a
        # {(n, m): ...} map from this depends on it.  Pinned by
        # test_bonds_yields_each_bond_once_with_the_lower_stable_id_FIRST.
        self._require_clean()
        # Hold the arena for the generator's whole life: a mutation mid-iteration repoints
        # self._structure, and without this reference the buffer would be freed under us.
        #
        # `pinned` protects the Structure OBJECT, not its buffer.  Building a lazy derived segment
        # -- component_labels(), stereo_units() -- calls structure_append, which reallocates
        # through PyMem_Realloc and may MOVE the buffer, and the body of a `for b in m.bonds():`
        # loop is free to do exactly that.  So ptr/edges are re-read on every resume rather than
        # cached across the yield; the CSR itself cannot change while the arena is clean, so
        # re-reading costs two loads and changes no answer.
        cdef Structure pinned = self._structure
        cdef list numbers = self._numbers
        cdef uint32_t gen = self._gen
        cdef uint32_t i, k, to, begin, end
        cdef Bond bd
        for i in range(pinned.header.atom_count):
            begin = csr_ptr(pinned)[i]
            end = csr_ptr(pinned)[i + 1]
            for k in range(begin, end):
                to = csr_edges(pinned)[k].to
                if to > i:
                    bd = Bond.__new__(Bond)
                    bd._molecule = self
                    bd._n = <uint32_t> numbers[i]
                    bd._m = <uint32_t> numbers[to]
                    bd._gen = gen
                    yield bd

    # ================================================================================
    # THE OPERATOR SURFACE.
    #
    # The dunders and the `*_of(n)` accessors -- `atom_count`, `order_of`, `element_of`, `to_bytes`
    # -- coexist by design.  An accessor takes no view object and allocates nothing, so it is the
    # fast path a tight loop wants, and it sits BESIDE the pretty spelling rather than instead of
    # it.  `atom(n)` and `bond(n, m)` return real views with real properties; the dunders are the
    # same surface at the container level.
    #
    # `__eq__` AND `__hash__` REST ON `canonical_bytes`, AND ON NOTHING ELSE.  Three candidate
    # substrates were measured; two of them are wrong and the record of why belongs here, because
    # both wrong ones are one line long and the right one is not:
    #
    #   * `_union_feature_words` is the OR of every atom's feature words -- a prefilter, lossy by
    #     construction.  Propane, butane and pentane share one value, so it would report
    #     `smiles('CCC') == smiles('CCCC')` as True.  A permissive `==` silently merges rows in
    #     every set, dict and dedup that touches it, which is the worst available failure.  It is
    #     private now precisely so this cannot be reached for by accident.
    #   * a canonical SMILES string is right in KIND and sound only when the labelling behind it is
    #     the extremal one: a labelling that stops at the first discrete leaf answers many strings
    #     for one compound, measured as forty over sixty relabelings of one cubane skeleton.
    #   * `canonical_bytes` -- the extremal labelling's certificate plus one parity digit per
    #     position -- is the sound one.  Measured: cis and trans 2-butene give two values; the four
    #     parity combinations of hexa-2,4-diene collapse to THREE, with (2E,4Z) and (2Z,4E) sharing
    #     one because they are one compound; and 720 creation orders of each configuration give one
    #     value each.
    #
    # `__eq__` REJECTS ON THE CHEAP DISCRIMINATORS FIRST -- atom count, bond count, then the union
    # words -- and only then pays for a canonical form.  All three are exact rejections: two
    # molecules that differ in any of them cannot be equal, so a mismatch is a decision and not a
    # guess.  The union words earn their keep exactly here, as a screen, which is what they are.
    #
    # WHICH PUTS A CONDITION ON THE WORDS THAT ONE BIT DID NOT MEET: to be an exact rejection a screen
    # feature must be FRAME-FREE, a property of the compound rather than of the atom order it happens
    # to be stored in.  Word IV's bit 6 is the raw SEG_PARITY sign and is not, so `==` returned False
    # on pairs whose canonical forms were equal -- measured on the two spellings of meso-2,3-butanediol
    # and of cis-cyclohexane-1,2-diol.  It is masked out with W4_FRAME_FREE_MASK; the screen keeps the
    # frame-free "a parity is configured" pair in bits 7 and 8, which is the part that discriminates.
    # ================================================================================

    def __len__(self):
        """Atom count."""
        self._require_clean()
        return self._structure.header.atom_count

    def __bool__(self):
        """True when the molecule has atoms. An empty container is falsy; note that `len(mol) == 0`
        and `not mol` therefore agree, which is what makes `if mol:` mean what it reads as."""
        self._require_clean()
        return self._structure.header.atom_count != 0

    def __iter__(self):
        """Iterate atom NUMBERS, not `Atom` views.

        Numbers even though views would look friendlier: `for n in mol` beside `mol.atom(n)` is one
        idiom, while yielding views would make `for n in mol: mol.atom(n)` a type error and would
        silently change what `list(mol)`, `set(mol)` and `dict.fromkeys(mol)` mean.  Iterate
        `mol.atoms()` for the views.

        A copy of the id list, so mutating the molecule mid-loop cannot corrupt the iteration --
        the loop then walks the atoms as they were when it started, which is the only behaviour a
        list can honestly offer.
        """
        self._require_clean()
        return iter(list(self._numbers))

    cdef bint _any_r(self, uint32_t index, bint any_index) noexcept:
        """Does the R bucket hold a marker -- any marker when `any_index`, else one whose index is
        `index`?  The bucket is contiguous, so the walk is bounded by the molecule's R count."""
        cdef uint32_t begin = element_bucket_begin(self._structure, 0)
        cdef uint32_t end = element_bucket_end(self._structure, 0)
        if any_index:
            return end > begin
        cdef atom_t *atoms = self._structure.atoms()
        # The bucket bounds index the element index list, which lives 120 words past its own header.
        cdef uint32_t *idx = structure_element_index(self._structure) + 120
        cdef uint32_t k
        for k in range(begin, end):
            if <uint32_t> at_r_index(&atoms[idx[k]]) == index:
                return True
        return False

    def __contains__(self, item):
        """Four questions behind one operator, dispatched on the type of `item`:

            n in mol          an int is an ATOM NUMBER: is that atom present
            'C' in mol        a str is an ELEMENT SYMBOL: does any atom have it
            query in mol      a QueryContainer: does it match anywhere in this molecule
            mol2 in mol       another molecule: is it a substructure of this one

        Both the int and the str reading are wanted often enough to be worth the dispatch, and no atom
        number is ever a string, so nothing is ambiguous.  An int is NOT read as an atomic number:
        `6 in mol` asks about atom 6, and `'C' in mol` is how you ask about carbon.

        The symbol question goes through the element index, so it costs a bucket lookup rather than a
        scan.  An unknown symbol is False rather than an error: `'Xx' in mol` is a question with an
        answer, and so is `'R0'` -- an unindexed R spells `R`.  `'R'` asks about ANY R and costs the
        same lookup; `'R7'` asks about one index and walks the R bucket, which the R count bounds.
        """
        cdef uint32_t number
        self._require_clean()
        if isinstance(item, str):
            if item == 'R' or (item.startswith('R') and item[1:].isdigit() and item != 'R0'):
                if item == 'R':
                    return self._any_r(0, True)
                if int(item[1:]) > R_INDEX_MAX:
                    return False
                return self._any_r(<uint32_t> int(item[1:]), False)
            number = SYMBOL_TO_NUMBER.get(item, NOT_AN_ELEMENT)
            if number == NOT_AN_ELEMENT:
                return False
            return element_bucket_end(self._structure, number) > \
                element_bucket_begin(self._structure, number)
        if isinstance(item, int) and not isinstance(item, bool):
            return item in self._index_of
        if isinstance(item, QueryContainer):
            return (<QueryContainer> item).is_substructure(self)
        if isinstance(item, MoleculeContainer):
            return (<MoleculeContainer> item).is_substructure(self)
        raise TypeError('membership test wants an atom number, an element symbol, a '
                        'QueryContainer or a MoleculeContainer, not %s'
                        % type(item).__name__)

    def __int__(self):
        """Total formal charge.

        A sum over the atoms rather than a stored total: nothing in the arena maintains one, and a
        stored total is a second truth that a `set_charge` can contradict.
        """
        self._require_clean()
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t i
        cdef int acc = 0
        for i in range(self._structure.header.atom_count):
            acc += <int> atoms[i].charge
        return acc

    def __float__(self):
        """Molecular mass in daltons.

        Each atom contributes the exact mass of its isotope, or the abundance-weighted average over
        the natural isotopes when no isotope is stated, plus one average hydrogen mass per IMPLICIT
        hydrogen.  Explicit hydrogens are atoms and are already in the sum.

        AN ATOM WHOSE IMPLICIT COUNT IS UNKNOWN CONTRIBUTES NO HYDROGEN MASS, so the result comes out
        LIGHT by however many hydrogens the record failed to state, and says nothing about it.  That
        is deliberate: a mass is not the place to raise, and there is no better number to use.  Ask
        `unknown_h_count` first if the answer has to be trusted; it is non-zero on exactly the records
        where this mass is a lower bound rather than a mass.
        """
        self._require_clean()
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t i
        cdef double acc = 0.0
        cdef double h = element_mass(1, 0)
        for i in range(self._structure.header.atom_count):
            acc += element_mass(atoms[i].element, atoms[i].isotope)
            if not at_implicit_h_unknown(&atoms[i]):
                acc += h * <double> at_implicit_h(&atoms[i])
        return acc

    def __bytes__(self):
        """The arena's persistent prefix, i.e. `to_bytes()`.

        THE ARENA AND NOT THE PACH RECORD.  `bytes(mol)` is lossless and this build's own format;
        `pack()` is the small, lossy, frozen pach record another chython reads.  V2 answers `pack()`
        here, so a consumer that stored `bytes(mol)` gets a buffer of the other kind.
        """
        return self.to_bytes()

    def __copy__(self):
        """`copy.copy(mol)` is `mol.copy()`: a new container over the SAME immutable arena, O(1)."""
        return self.copy()

    cdef str _smiles(self):
        """`write_smiles(self, '')`, CACHED, keyed on the same generation counter `atoms_order` uses.

        THE DEFAULT SPEC ONLY, and the narrowness is the design rather than a shortcut.  Measured
        2026-09-03 on 4999 public records: a canonical write is 6.14 µs/molecule, of which
        `canonical_order()` is 4.62 -- so `str(mol)` in a loop over a hundred thousand molecules is
        two thirds of a second of relabelling the same graphs.  This cache is what removes that.

        One field and not a `{spec: string}` dict, for three reasons that all point the same way.  The
        empty spec is what `str`, `repr`, `f'{mol}'` and `.smiles` all ask for, so it is the call that
        happens in a loop; a dict would need `normalize_smiles_spec` to avoid storing `'sm'` and `'ms'`
        separately, and that parse costs a measurable fraction of the write it is trying to save; and a
        per-spec cache grows with whatever specs a caller happens to pass, which is a footprint decision
        this class should not be making on their behalf.  A caller who writes one molecule under many
        specs in a hot loop wants `write_smiles` and their own dict keyed on `normalize_smiles_spec`,
        which exists and whose docstring states exactly the invariant that makes such a cache sound.

        INVALIDATION IS `_gen` AND NOTHING ELSE, the same rule `_identity` documents at length: `_gen`
        is bumped by `_apply` and by every in-place writer, `_require_clean` refuses a read while a
        journal is pending, and `kekule`/`thiele` bump it too -- which they must, since they change the
        string.  It over-invalidates, because a coordinate or wedge change bumps `_gen` and reaches
        nothing here, and that is the safe direction.

        THE TWO STORES ARE ORDERED and for the reason `_identity` gives: the string goes down BEFORE
        the generation it belongs to, so a torn read misses the cache and recomputes rather than
        pairing a fresh generation with a stale string.

        `_require_clean` COMES BEFORE THE CACHE READ, and getting that backwards is a bug this method
        shipped for about ten minutes.  `_gen` is bumped by `_apply`, i.e. when an edit scope CLOSES,
        so inside an open scope the counter still matches the row that was stored before the scope
        opened -- and a cache read placed first happily returns the pre-scope string for a molecule the
        caller has just added three atoms to.  A plain `write_smiles` would have refused.  So the guard
        that every other reader on this class runs is run here first, and the cache is consulted only
        once the container has agreed to be read at all.
        """
        self._require_clean()
        if self._smiles_cache is not None and self._smiles_gen == self._gen:
            return self._smiles_cache
        cdef str out = write_smiles(self, '')
        self._smiles_cache = out
        self._smiles_gen = self._gen
        return out

    def __str__(self):
        """Canonical SMILES."""
        return self._smiles()

    def __repr__(self):
        """`smiles('...')` -- the canonical SMILES as an expression that rebuilds this molecule.

        Eval-able rather than pretty, because that is what `repr` is for and because the default
        `<chython.core._core.MoleculeContainer object at 0x...>` tells a chemist at a REPL nothing at
        all about the molecule in hand.  `repr` and not `__str__` carries the wrapper: `str(mol)` is
        the string you paste into a file, `repr(mol)` is the line you paste back into Python.

        IT DOES NOT RAISE, and that is the whole reason this is not a one-liner.  A `repr` is called
        by debuggers, by tracebacks, by `%r` in the error messages of other code, and by every REPL --
        i.e. exactly when something is already wrong -- so a molecule with a pending journal, or one
        whose configuration the writer refuses, must still be printable.  What comes out then names
        the container and its size, which is the honest answer, and the exception is re-spelled in it
        rather than swallowed: a `repr` that hid the reason would send the reader looking in the wrong
        place.  Nothing else in this class is allowed to be lenient like this.

        The empty molecule is `MoleculeContainer()` and not `smiles('')`, because the point of this
        method is an expression that rebuilds the object and the SMILES reader does not accept an empty
        string -- so the one case where the general form would not round-trip is spelled out.

        The atom count in the fallback is the ARENA's, which under a pending journal is the pre-scope
        count and not what the caller has just added.  That is the same state every other reader on
        this class sees and refuses to answer from, so it is the state a diagnostic should show.
        """
        # declared rather than left to `except ... as`, which Cython reports as an implicit
        # declaration and this tree treats a Cython warning as a build failure
        cdef object e
        if not self._structure.header.atom_count and not self._journal_len:
            return 'MoleculeContainer()'
        try:
            return 'smiles(%r)' % self._smiles()
        except Exception as e:
            return '<MoleculeContainer, arena holds %d atoms, unwritable: %s: %s>' % (
                self._structure.header.atom_count, type(e).__name__, e)

    @property
    def smiles(self):
        """Canonical SMILES, i.e. `str(self)`.

        A property and not only `str`, because it reads better in a comprehension and it is the
        spelling every caller in the wild uses.  For anything other than the default spec use
        `format(mol, spec)`; there is deliberately no `smiles_without_stereo`-style property per key.

        Cached with `str` and `repr`, so reading it twice costs one write -- which also means it is a
        property that can be read in a comprehension without a performance surprise, the thing that
        makes a property the wrong shape when it is not true.
        """
        return self._smiles()

    def __format__(self, str format_spec not None):
        """Canonical SMILES with a format spec, so `f'{mol:a}'` and `format(mol, 'a')` work.

        The spec goes straight to `write_smiles`, which owns what the letters mean; this is a
        forwarding method and deliberately does not document them a second time.  Three things about
        it are not the writer's business but a caller's, so they are here:

        THE SPEC CAN CHANGE THE ATOM ORDER, NOT ONLY THE TOKENS.  With the stereo-seeded canonical
        order, `format(mol, 's')` and `format(mol, '!s')` can order the ATOMS differently, not merely
        add and remove stereo signs, because the stereo configuration is part of what the canonical
        order is computed from.  Anyone diffing the two strings positionally, or reusing an atom
        order fetched under one spec while writing under another, gets silently wrong answers.  Fetch
        the order under the same spec you write under, every time.

        `format(mol, 'i')` IS STORED SLOT ORDER AND IS EXPLICITLY NOT CANONICAL.  It is a debugging
        view of the arena's own layout, it changes when the arena is rebuilt, and two records of one
        compound built in different orders give different strings.  It must never reach a hash, a
        dict key or an equality test -- `==` and `hash` use `canonical_bytes` and do not go through
        this method at all, which is the point.  `format(mol, 'r')` is the same warning drawn fresh per
        call: a random order, for augmentation, and two calls are two strings of one molecule.

        AN UNKNOWN KEY RAISES `ValueError` and this method does not soften it: a silently ignored key
        means the caller asked for one thing and shipped another.

        THE EMPTY SPEC IS CACHED and every other one is not, because `format(mol)` and `f'{mol}'` are
        the same call as `str(mol)` and must not be the slow spelling of it.  A caller who writes one
        molecule under several specs repeatedly wants `write_smiles` with their own dict keyed on
        `normalize_smiles_spec`; `_smiles` says why this class does not keep that dict for them.
        """
        if not format_spec:
            return self._smiles()
        return write_smiles(self, format_spec)

    def sticky_smiles(self, left=None, right=None, *, bint remove_left=False,
                      bint remove_right=False, tries=10, bint keep_bond_left=False,
                      bint keep_bond_right=False, bint hydrogens=False):
        """A SMILES that STARTS at atom `left` and ENDS at atom `right`, for a caller who glues strings.

        THE SIGNATURE IS FIXED, because consumers outside this repository call it: `left` and
        `right` are atom ids, `remove_*` drops that end's atom token and `keep_bond_*` keeps its bond
        token, and `hydrogens` shows every implicit count.  `sticky_smiles(left=n, remove_left=True,
        keep_bond_left=True)` gives `-CCO`, and `A + B` where A ends open and B starts open is a
        molecule.  Everything the letters mean is at the core `sticky_smiles`, which this forwards to.

        `tries` IS ACCEPTED AND IGNORED.  The core constrains the traversal and proves it ends where it
        was asked to (`smw_sticky_traverse`), so there is nothing to retry and no failure to retry it
        for.  Kept in the signature so that existing calls -- both in-repo callers pass it -- do not
        have to change; it will go when the signature is next allowed to.

        NOT CANONICAL, NOT CACHED and no CXSMILES tail: the order depends on the atoms named, and a tail
        index counts atoms from the start of a string a caller is about to prepend to.  Radicals and
        stereo groups are therefore NOT carried by the returned text.
        """
        tries       # ignored, and NAMED so that Cython's unused-argument warning does not fire on it:
                    # this module's build gate is zero Cython warnings, and a suppression comment would
                    # be a promise the compiler does not check.
        return sticky_smiles(self, left, right, 'h' if hydrogens else '', remove_left=remove_left,
                             remove_right=remove_right, keep_bond_left=keep_bond_left,
                             keep_bond_right=keep_bond_right)

    def detached_smiles(self, cuts not None, str spec='', reserve=None):
        """This molecule minus the dropped side of every cut, the cut bonds left as open ring bonds.

        `cuts` is `{attachment_id: (keep_n, drop_n)}` -- ORDERED pairs, nothing in `C-C` saying which
        half the caller wants -- and `reserve` withholds the other fragments' attachment ids so several
        fragments can be written for one join.  Answers a `DetachedSmiles`.  The module-level
        `detached_smiles` in `_smiles_write.pxi` owns the rules and the four refusals; this forwards.

        NOT CANONICAL AND NOT CACHED, for `sticky_smiles`' reason and one of its own: the cuts and the
        reserved ids are the caller's, so no spec identifies the result.  That is also why it is a
        method rather than a `format()` key -- `format(mol, spec)` promises a pure function of the
        molecule and the spec, and this takes two arguments neither of them can carry.
        """
        return detached_smiles(self, cuts, spec, reserve)

    def as_query(self):
        """This molecule as a `QueryContainer` that matches it and its supergraphs.

        WHAT IS DEMANDED, and every omission below is deliberate: element, isotope when one is
        stated, formal charge and the radical flag per atom, and the order per bond.

        WHAT IS NOT DEMANDED, because a substructure test must survive embedding.  Not degree --
        the whole point is that a fragment may have neighbours the query never mentioned.  Not
        hydrogen counts, and that omission is exactly what makes Ramil's `smiles('CO') <
        smiles('COC')` come out True: methanol's oxygen carries one hydrogen and the ether's
        carries none, so a query demanding `implicit_h 1` would refuse the very comparison the
        operator exists for.  Not heteroatom count, hybridization, ring membership or ring size,
        all of which are properties of the whole molecule rather than of the fragment.  Not the map
        number: AAM is annotation.  Not stereo -- see the note on the comparison operators.

        An atom with no stated isotope demands nothing about isotopes, so `smiles('C')` matches
        `smiles('[13CH4]')`.  That is the permissive reading, and it is the right one for the same
        reason as the hydrogens: an unstated isotope is the absence of a demand, not a demand for
        the absence.  `no_isotope` is the primitive for a caller who means the other thing.

        An aromatic stored bond is demanded as aromatic, so it matches an aromatic stored bond and
        NOT its Kekule twin: `benzene <= cyclohexane` is False in both directions and so is
        `aromatic benzene <= Kekule benzene`.  A caller comparing two differently-written records
        must kekulise both first; nothing here normalises a representation behind its back.

        A molecule carrying an R marker is REFUSED rather than converted.  The marker matches nothing,
        so the query would be structurally incapable of matching and every comparison through it would
        answer False for a reason the caller never asked about.  The reverse question is fine and
        answers False: an R-bearing molecule is a legitimate substructure TARGET.
        """
        self._require_clean()
        # An R matches nothing, so a query demanding one could never match: refused here rather than
        # returning a query that is silently always False.  `as_query()` is the one door -- the
        # comparison operators and `is_substructure` all route through it, and `QueryContainer` itself
        # refuses element 0 at `atom_primitive`, so no hand-written query can hold one.
        if element_bucket_end(self._structure, 0) > element_bucket_begin(self._structure, 0):
            raise ValueError('this molecule carries an R marker, which matches nothing, so it cannot '
                             'be used as a query. Compare fragments by their canonical SMILES, or '
                             'write the pattern in SMARTS, where an attachment point is unwritten.')
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t *ptr = csr_ptr(self._structure)
        cdef halfedge_t *edges = csr_edges(self._structure)
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef uint32_t i, k, qn
        cdef list qids = []
        cdef QueryContainer q = QueryContainer()
        with q.edit():
            for i in range(n_atoms):
                qn = q.add_atom()
                qids.append(qn)
                # `and_low` BETWEEN the terms, not after: the journal is a token stream and an
                # operator is the juxtaposition, so a run of bare primitives is a malformed term
                # and query_seal says so at the first match rather than here.
                q.atom_primitive(qn, 'element', <int32_t> atoms[i].element)
                q.atom_operator(qn, 'and_low')
                q.atom_primitive(qn, 'charge', <int32_t> atoms[i].charge)
                q.atom_operator(qn, 'and_low')
                # `radical` IGNORES its value: the positive term means "is a radical" and the
                # NEGATED term means "is not one".  Passing `value=0` and expecting "not a radical"
                # compiles, seals and then matches nothing at all, which is how this was found.
                q.atom_primitive(qn, 'radical', 0, not at_radical(&atoms[i]))
                if atoms[i].isotope:
                    q.atom_operator(qn, 'and_low')
                    q.atom_primitive(qn, 'isotope', <int32_t> atoms[i].isotope)
            for i in range(n_atoms):
                for k in range(ptr[i], ptr[i + 1]):
                    if edges[k].to > i:
                        q.add_bond(<uint32_t> qids[i], <uint32_t> qids[edges[k].to])
                        # A STORED ORDER-4 BOND IS SPELLED `bond_aromatic`, NOT `bond_order 4`:
                        # the order primitive accepts 1, 2, 3 and 8 only and refuses 4 at seal.
                        # That is not an oversight to route around -- 4 and 8 share one feature
                        # bit and are separated by the aromatic bit alone, so the aromatic demand
                        # IS the order-4 demand, spelled as the thing it actually tests.
                        if edges[k].order == 4:
                            q.bond_primitive(<uint32_t> qids[i], <uint32_t> qids[edges[k].to],
                                             'bond_aromatic', 0)
                        else:
                            q.bond_primitive(<uint32_t> qids[i], <uint32_t> qids[edges[k].to],
                                             'bond_order', <int32_t> edges[k].order)
        return q

    def is_substructure(self, other):
        """Is this molecule a substructure of `other`?  See `as_query` for what that demands.

        `other` may be a `MoleculeContainer`.  Costs one query compilation per call and caches
        nothing, so a caller testing one fragment against many molecules should build the query
        once with `as_query()` and reuse it.
        """
        if not isinstance(other, MoleculeContainer):
            raise TypeError('is_substructure wants a MoleculeContainer, not %s'
                            % type(other).__name__)
        return bool(self.as_query().is_substructure(<MoleculeContainer> other))

    def __eq__(self, other):
        """Same compound?  Equal iff the two canonical forms are equal.

        WHAT COUNTS AS THE SAME COMPOUND: the graph, elements, isotopes, charges, radicals, implicit
        hydrogen counts, bond orders, aromatic bits and one stereo parity per centre.  Atom NUMBERS
        do not count, nor does creation order, nor do coordinates, wedges or map numbers -- `==` is a
        question about the compound and those four are annotation.  A REPRESENTATION DOES count: an
        aromatic ring and its Kekule twin are two different records of two different inputs and they
        compare unequal, which is what storing an input faithfully means.  Kekulise both first if
        that is not the question you meant to ask.

        Enhanced stereo (ABS / AND / OR groups) is NOT part of this yet, so a racemate and a single
        enantiomer of the same skeleton compare equal.  That is a real gap and it is stated rather
        than hidden; the groups are stored, they are simply not in the canonical form.

        Raises `AutomorphismBudgetExceeded` if a canonical search is truncated, and it is allowed to:
        an exception from `==` is a caller's problem to see, while a fallback answer would be a wrong
        one that nobody sees.
        """
        if not isinstance(other, MoleculeContainer):
            return NotImplemented
        cdef MoleculeContainer o = <MoleculeContainer> other
        if self is o:
            return True
        self._require_clean()
        o._require_clean()
        if self._structure.header.atom_count != o._structure.header.atom_count:
            return False
        if self._structure.header.bond_count != o._structure.header.bond_count:
            return False
        # `copy()` shares the arena, and the arena is immutable, so this is not an optimisation of a
        # rare case: it is the case every `copy()`, `substructure()` and dict round trip produces.
        if self._structure is o._structure:
            return True
        cdef uint64_t *fa = structure_features(self._structure)
        cdef uint64_t *fb = structure_features(o._structure)
        cdef uint32_t i
        # WORDS I..III ONLY, AND WORD IV WITH ITS FRAME-RELATIVE BIT MASKED OFF.  Word IV bit 6 is the
        # raw SEG_PARITY sign, a statement in each molecule's own slot frame and not a property of the
        # compound -- so two spellings of one meso compound differ there while their canonical forms
        # agree, and an unmasked screen rejects them.  See W4_FRAME_FREE_MASK.
        for i in range(3):
            if fa[i] != fb[i]:
                return False
        if fa[3] & W4_FRAME_FREE_MASK != fb[3] & W4_FRAME_FREE_MASK:
            return False
        return self._identity() == o._identity()

    def __hash__(self):
        """Hash of the canonical form, so equal molecules hash alike and a molecule is a usable dict
        key and set member.

        NOT built from the union feature words even though they are 32 bytes already sitting in the
        arena: they are lossy, a hash built on them would put propane and pentane in one bucket, and
        while a colliding hash is legal it would make every `set` of molecules degenerate into a
        linear scan of `__eq__` calls.  A hash should discriminate as well as equality does.

        Costs one canonical form on first call and nothing afterwards -- see `_identity` for the
        cache and its invalidation -- so hashing a molecule twice is cheap and hashing a MUTATED
        molecule is correct.  A molecule is mutable, which by the usual Python rule argues against
        hashing it at all; hashing it anyway is deliberate, because the alternative is that no
        molecule can be a dict key, and the generation counter makes the answer follow the mutation
        rather than go stale.  A caller who mutates a molecule while it sits in a set still gets what that always
        gets in Python: a member that can no longer be found.
        """
        return hash(self._identity())

    # THE FOUR COMPARISONS ARE CONTAINMENT TESTS.
    #
    #     mol1 <= mol2     mol1 is a substructure of mol2
    #     mol1 <  mol2     ... and is strictly smaller
    #     mol1 >= mol2     mol2 is a substructure of mol1
    #     mol1 >  mol2     ... and mol1 is strictly larger
    #
    # THE `len()` GUARDS ON THE STRICT PAIR ARE NOT DECORATION: they make `<` antisymmetric by
    # construction, so
    # `a < b` and `b < a` cannot both hold, and they reject before the isomorphism kernel runs, which
    # is the same free short-circuit `__eq__` opens with.  All four are written out explicitly and
    # `functools.total_ordering` is deliberately NOT used: deriving `<=` from `<` is valid only for a
    # total order, and substructure containment is a poset -- benzene and cyclohexane are simply
    # incomparable, both `<` and `>` False.
    #
    # TO PUT MOLECULES IN A DETERMINISTIC SEQUENCE, sort with an explicit key, and pick the one that
    # matches what the sequence is for:
    #
    #   * `sorted(mols, key=str)` -- determinism.  Lexicographic on the canonical SMILES, which is
    #     stereo-seeded and order-stable, so it is exact and ties only between identical compounds.
    #     This is the one genuine TOTAL order on molecules.  Use it for stable file output and dedup.
    #   * `sorted(mols, key=float)` -- chemistry.  Molecular mass; ties across every isomer set, and
    #     a float, so the ties are unstable near-ties.  Fine for "roughly smallest first", useless
    #     for reproducibility.
    #   * a formula -- grouping.  A key that collects isomers together; it is not an ordering.
    #
    # NOT STEREO-AWARE, on both sides: `as_query` demands no parity, so one enantiomer is a
    # substructure of the other.  That is the honest answer until a query can carry a parity demand
    # built from a molecule.
    #
    # A `QueryContainer` IS ACCEPTED, on the CONTAINED side of all four: `mol >= q` and `mol > q` ask
    # whether the query matches inside the molecule, and `q <= mol` / `q < mol` are the same question
    # spelled the other way round and live on `QueryContainer`.  The other four combinations --
    # `mol <= q`, `mol < q`, `q >= mol`, `q > mol` -- would ask whether a molecule embeds in a
    # pattern, which the kernel cannot answer in that direction, and they raise `TypeError` rather
    # than quietly answering the question that was not asked.

    def __le__(self, other):
        if isinstance(other, MoleculeContainer):
            return self.is_substructure(<MoleculeContainer> other)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, MoleculeContainer):
            if len(self) >= len(<MoleculeContainer> other):
                return False
            return self.is_substructure(<MoleculeContainer> other)
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, MoleculeContainer):
            return (<MoleculeContainer> other).is_substructure(self)
        if isinstance(other, QueryContainer):
            return bool((<QueryContainer> other).is_substructure(self))
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, MoleculeContainer):
            if len(self) <= len(<MoleculeContainer> other):
                return False
            return (<MoleculeContainer> other).is_substructure(self)
        if isinstance(other, QueryContainer):
            # a query's `len` is its atom count, so the same guard applies and means the same thing
            if len(self) <= len(<QueryContainer> other):
                return False
            return bool((<QueryContainer> other).is_substructure(self))
        return NotImplemented

    def substructure(self, atoms):
        """A new molecule holding just `atoms` and the bonds between them, numbers preserved.

        WHAT SURVIVES: element, isotope, charge, radical, map number, the implicit hydrogen count
        as stored, and coordinates -- 2D and 3D alike -- when the source has them.  Atom numbers are the source's, so
        `mol.substructure([3, 7])` returns a molecule whose atoms are 3 and 7, which is why this builds
        and then `remap`s rather than handing back 1 and 2.

        WHAT DOES NOT SURVIVE, and is not silently approximated: stereo parities, wedges, stereo
        groups and CIP descriptors.  A parity is a statement about a frame of NEIGHBOURS, and cutting
        bonds destroys the frame it was stated in; retranslating one needs `translate_stereo`'s
        machinery and is stereo-epic work.  Until then a substructure is stereo-free and says so,
        which is the failure a caller can see.

        A CIP DESCRIPTOR IS LOST HERE WITHOUT A `cip_log` LINE, and that is the one place the log is
        not the record of a loss.  An edit does not ask to lose one, so an edit reports it; a cut is a
        caller asking for a smaller molecule, and this paragraph is where the answer lives -- the same
        standing this operation already gives parities, which are not logged either.  A descriptor is
        a ranking at a centre and a ranking reads the WHOLE molecule, so a cut anywhere can change it.

        THE HYDROGEN COUNTS ARE NOW WRONG, DELIBERATELY.  Cutting a bond leaves the atom that lost
        it with the count it had, because the core does not derive hydrogens at all -- deriving them
        is `calc_implicit`'s job in the standardization pass, and a container that guessed here
        would be inventing chemistry inside a graph operation.  Run the repair pipeline on the
        result if you need a valid molecule; this returns a faithful cut.
        """
        cdef set keep = set()
        cdef object x
        self._require_clean()
        for x in atoms:
            if x not in self._index_of:
                raise KeyError(x)
            keep.add(x)
        if not keep:
            raise ValueError('a substructure of no atoms is not a molecule')
        cdef atom_t *src = self._structure.atoms()
        cdef uint32_t *ptr = csr_ptr(self._structure)
        cdef halfedge_t *edges = csr_edges(self._structure)
        cdef bint want_xy = self.has_coordinates
        cdef uint32_t models = structure_conformer_count(self._structure)
        cdef uint32_t model
        cdef Conformer conf
        cdef uint32_t i, k, n
        cdef dict back = {}          # new id -> source id, which is exactly remap's argument
        cdef dict new_of_slot = {}   # source slot -> new id, for the bond pass
        cdef MoleculeContainer mol = MoleculeContainer()
        with mol.edit():
            for i in range(self._structure.header.atom_count):
                n = <uint32_t> self._numbers[i]
                if n not in keep:
                    continue
                # `at_implicit_h` is the RAW nibble, so an unknown count copies across as the
                # sentinel rather than as a zero -- an atom nobody had a count for must not acquire
                # one by being cut out of a bigger molecule.
                new_of_slot[i] = mol.add_atom(<int> src[i].element, charge=<int> src[i].charge,
                                            isotope=<int> src[i].isotope,
                                            radical=at_radical(&src[i]),
                                            map_number=<int> src[i].map_number,
                                            implicit_h=<int> at_implicit_h(&src[i]))
                back[new_of_slot[i]] = n
            # a second pass for the bonds: every endpoint must exist before any bond is journalled
            for i in new_of_slot:
                for k in range(ptr[i], ptr[i + 1]):
                    if edges[k].to > i and edges[k].to in new_of_slot:
                        mol.add_bond(<uint32_t> new_of_slot[i],
                                   <uint32_t> new_of_slot[edges[k].to],
                                   <int> edges[k].order)
        mol.remap(back)
        if want_xy or models:
            with mol.edit():
                if want_xy:
                    for n in mol._numbers:
                        mol.set_xy(n, *self.xy_of(n))
                # ADDED FOR MODEL 0 TOO, rather than letting the first `set_xyz` create it: one
                # uniform loop covers every model, and only an explicit add carries the stated number.
                for model in range(models):
                    conf = self.conformer(model)
                    mol.add_conformer(ext_index=conf.ext_index)
                    for n in mol._numbers:
                        mol.set_xyz(n, *conf.xyz_of(n), model=model)
        return mol

    def split(self):
        """This molecule's connected components, one `MoleculeContainer` each.

        A LIST, always -- a one-component molecule comes back as a list holding one copy of itself,
        never as a bare molecule, so a caller never has to test the return type.  An empty molecule
        gives an empty list.  Components come in first-seen atom order and each keeps the source's
        atom numbers, so `union` of the pieces reproduces the numbering.

        THIS REPORTS THE COMPONENTS THE MOLECULE ALREADY HAS, AND A SALT OFTEN IS NOT DRAWN AS ONE.
        `CC(=O)O[Na]` is ONE component here, because the record drew a covalent Na-O bond and deciding
        that a drawn bond is wrong is a chemistry judgement rather than a graph one.  A caller who splits
        that and concludes the record carries no counterion has a silently wrong answer, which is the
        expensive kind.  Run `split_salts()` first if that matters -- on a `copy()` when the original
        must survive -- and then `split()` gives the two pieces with their charges on them.

        UNLIKE `substructure`, THIS PRESERVES STEREO, and the reason is that it cuts no bond.  A
        parity is a statement about a frame of neighbours; splitting components leaves every atom
        with exactly the neighbours it had, so every frame survives and is re-based rather than
        dropped.  Parities, wedges and stereo groups all come through.  Implemented as `copy()` plus
        the other components' `delete_atom`s for precisely that reason -- routing it through
        `substructure` would throw away configurations no atom lost.

        THE ONE LOSS IS CIP: a deletion anywhere drops every stored descriptor, because a ranking
        reads the whole molecule and the arena will not reason about which deletions could not have
        changed one.  The drop is counted in `cip_log`, as it is for any edit.
        """
        return mol_split(self)

    def augmented_substructure(self, atoms, int deep=1):
        """The environment of `atoms` out to `deep` bonds, as a molecule.

        `deep` COUNTS BONDS: 0 is `substructure(atoms)` itself, 1 adds the seed's direct
        neighbours, 2 their neighbours, and so on.  The seed is always included.  A radius that
        exceeds the seed's own components is not an error -- it saturates, and the answer is those
        components entire.

        A MOLECULE AND NOT A PROJECTION, built by `substructure`, so everything that docstring says
        applies here word for word: numbers are the source's, bonds that leave the selection are
        gone, hydrogen counts are the source's and therefore wrong on any atom that lost a bond, and
        stereo does not survive.  That last one is the difference from `split`, and it is not a
        shortcut: this DOES cut bonds, so the frame a parity was stated in is genuinely destroyed.
        """
        cdef list levels = mol_augmented_levels(self, atoms, deep)
        # `levels[-1]` and NOT: this translation unit compiles with `wraparound=False`, which turns a
        # negative list index into an unchecked read past the front of the array -- a segfault, found
        # by writing exactly that line first.
        return self.substructure(levels[len(levels) - 1])

    def augmented_substructures(self, atoms, int deep=1):
        """Every shell of `augmented_substructure`, from the seed alone outwards, as a list.

        `[0]` is the seed, `[1]` the seed plus its neighbours, and so on.  THE LIST MAY BE SHORTER
        THAN `deep + 1`: growth stops as soon as a shell adds no atom, so `deep=99` on a small
        molecule returns one entry per shell that exists rather than a hundred copies of the
        component.  Read `len()`, never `deep`.
        """
        cdef list out = []
        cdef object level
        for level in mol_augmented_levels(self, atoms, deep):
            out.append(self.substructure(level))
        return out

    def adjacency_matrix(self, bint set_bonds=False):
        """The `(n, n)` uint32 adjacency matrix: 1 where a bond exists, 0 where none does.

        `set_bonds=True` writes the bond's stored order instead of 1, so a double bond is a 2 and an
        aromatic bond a 4.  It is positional or keyword.

        ROWS AND COLUMNS ARE POSITIONS, NOT ATOM NUMBERS: row `i` is `self.atom_numbers[i]`.  A matrix
        keyed by atom number is not expressible once numbers are sparse, which they are after any
        deletion.  The matrix is symmetric, and a dative bond (order 8) is a bond here like any other.
        """
        self._require_clean()
        return mol_adjacency_matrix(self._structure, set_bonds)

    def distance_matrix(self):
        """The `(n, n)` int32 matrix of topological distances, counted in bonds.

        0 ON THE DIAGONAL, -1 WHERE THERE IS NO PATH.  The convention is the one chytorch's
        `graph_distances` needs: it adds 2 and documents "1 marks a pair in different components, 2
        an atom with itself, 3 neighbours", leaving 0 free for padding, so -1 / 0 / 1 shift onto
        exactly those.  Values are exact and unclamped -- a consumer that wants a cutoff applies it.

        Rows are positions, like `adjacency_matrix`: row `i` is `self.atom_numbers[i]`.  Symmetric, and
        a dative bond is a step like any other -- ring perception excludes order 8 because a dative
        bond closes no ring, but a walk has no such argument.
        """
        self._require_clean()
        return mol_distance_matrix(self._structure)

    def state_view(self, TensorEncoding encoding=None):
        """Per-atom int32 arrays for a model: element, hydrogens, heavy degree, distances.

        See `docs/ml.rst`.  numpy is an optional dependency, `chython[ml]`.
        """
        self._require_clean()
        return mol_state_view(self, encoding)

    def transition_view(self, TensorEncoding encoding=None):
        """Per-atom int32 arrays with a before and an after side; for a molecule the two are equal.

        See `docs/ml.rst`.  numpy is an optional dependency, `chython[ml]`.
        """
        self._require_clean()
        return mol_transition_view(self, encoding)

    # --- fingerprints ----------------------------------------------------------------------------

    def atom_invariants(self):
        """`(n,)` uint32: the default featurization label of every atom, in this molecule's order.

        Row `i` is `self.atom_numbers[i]`, the same indexing `adjacency_matrix` and `distance_matrix`
        use.  Pass a vector of this shape and dtype back as `invariants=` to any fingerprint method
        to fingerprint over a different atom typing.
        """
        self._require_clean()
        return fp_atom_invariants(self._structure)

    def morgan_hash_counts(self, int min_radius=1, int max_radius=4, *, invariants=None):
        """Circular fragments of radius `min_radius..max_radius` as `{hash: count}`, unfolded.

        Radius 1 is the atom by itself and radius `r` reaches `r - 1` bonds out, so the counts of a
        single radius sum to `atom_count`.  Comparable to ECFP / RDKit Morgan, and the bit values are
        chython's own -- no other toolkit's hashes are reproduced here.

        `invariants` takes a `(atom_count,)` uint32 vector to fingerprint over a different atom
        typing; `atom_invariants()` returns the default one.
        """
        self._require_clean()
        _fp_check_radii(min_radius, max_radius)
        return fp_morgan_counts(self._structure, <uint32_t> min_radius, <uint32_t> max_radius,
                                invariants)

    def morgan_hash_set(self, int min_radius=1, int max_radius=4, *, invariants=None):
        """The distinct circular fragment hashes, unfolded.  See `morgan_hash_counts`."""
        return set(self.morgan_hash_counts(min_radius, max_radius, invariants=invariants))

    def morgan_bit_set(self, int min_radius=1, int max_radius=4, int length=1024,
                       int number_active_bits=2, *, invariants=None):
        """The folded bit positions, `0 <= p < length`.  See `morgan_hash_counts`.

        Cheaper than `morgan_fingerprint` when the answer is a Tanimoto: `len(a & b) / len(a | b)`.
        """
        _fp_check_folding(length, number_active_bits)
        return fp_fold_bit_set(self.morgan_hash_counts(min_radius, max_radius,
                                                       invariants=invariants),
                               <uint32_t> length, <uint32_t> number_active_bits)

    def morgan_fingerprint(self, int min_radius=1, int max_radius=4, int length=1024,
                           int number_active_bits=2, *, invariants=None):
        """`(length,)` uint8 of 0 and 1, the folded binary fingerprint.  See `morgan_hash_counts`."""
        _fp_check_folding(length, number_active_bits)
        return fp_fold_binary(self.morgan_hash_counts(min_radius, max_radius,
                                                      invariants=invariants),
                              <uint32_t> length, <uint32_t> number_active_bits)

    def morgan_count_vector(self, int min_radius=1, int max_radius=4, int length=1024,
                            int number_active_bits=2, *, invariants=None):
        """`(length,)` uint32 of folded fragment counts.  See `morgan_hash_counts`.

        The counted analogue of `morgan_fingerprint`: a real count per bit, not one extra bit per
        repeat.
        """
        _fp_check_folding(length, number_active_bits)
        return fp_fold_counted(self.morgan_hash_counts(min_radius, max_radius,
                                                       invariants=invariants),
                               <uint32_t> length, <uint32_t> number_active_bits)

    def linear_hash_counts(self, int min_radius=1, int max_radius=4, *, invariants=None):
        """Simple paths of `min_radius..max_radius` ATOMS as `{hash: count}`, unfolded.

        The radii count atoms and not bonds: length 1 is a lone atom, length 2 is a bond, so the
        counts at length 2 sum to `bond_count`.  Comparable to RDKit's RDKFingerprint in spirit; the
        bit values are chython's own.

        `invariants` takes a `(atom_count,)` uint32 vector, exactly as in `morgan_hash_counts`.

        PATH ENUMERATION IS EXPONENTIAL IN `max_radius` ON DENSELY FUSED RINGS -- there is
        deliberately no cap, because truncating would return a wrong fingerprint rather than a slow
        one.  The default of 4 is cheap everywhere.
        """
        self._require_clean()
        _fp_check_radii(min_radius, max_radius)
        return fp_linear_counts(self._structure, <uint32_t> min_radius, <uint32_t> max_radius,
                                invariants)

    def linear_hash_set(self, int min_radius=1, int max_radius=4, *, invariants=None):
        """The distinct path fragment hashes, unfolded.  See `linear_hash_counts`.

        Path enumeration is exponential in `max_radius` on densely fused rings; no cap is applied.
        """
        return set(self.linear_hash_counts(min_radius, max_radius, invariants=invariants))

    def linear_bit_set(self, int min_radius=1, int max_radius=4, int length=1024,
                       int number_active_bits=2, *, invariants=None):
        """The folded bit positions, `0 <= p < length`.  See `linear_hash_counts`.

        Path enumeration is exponential in `max_radius` on densely fused rings; no cap is applied.
        """
        _fp_check_folding(length, number_active_bits)
        return fp_fold_bit_set(self.linear_hash_counts(min_radius, max_radius,
                                                       invariants=invariants),
                               <uint32_t> length, <uint32_t> number_active_bits)

    def linear_fingerprint(self, int min_radius=1, int max_radius=4, int length=1024,
                           int number_active_bits=2, *, invariants=None):
        """`(length,)` uint8 of 0 and 1, the folded binary fingerprint.  See `linear_hash_counts`.

        Path enumeration is exponential in `max_radius` on densely fused rings; no cap is applied.
        """
        _fp_check_folding(length, number_active_bits)
        return fp_fold_binary(self.linear_hash_counts(min_radius, max_radius,
                                                      invariants=invariants),
                              <uint32_t> length, <uint32_t> number_active_bits)

    def linear_count_vector(self, int min_radius=1, int max_radius=4, int length=1024,
                            int number_active_bits=2, *, invariants=None):
        """`(length,)` uint32 of folded path counts.  See `linear_hash_counts`.

        A real count per bit, not one extra bit per repeat.

        Path enumeration is exponential in `max_radius` on densely fused rings; no cap is applied.
        """
        _fp_check_folding(length, number_active_bits)
        return fp_fold_counted(self.linear_hash_counts(min_radius, max_radius,
                                                       invariants=invariants),
                               <uint32_t> length, <uint32_t> number_active_bits)

    # --- graph descriptors -----------------------------------------------------------------------
    #
    # Thin by policy: validation and one forward call, with the definitions, the papers and the
    # disconnected-molecule behaviour in `_descriptors.pxi` beside the arithmetic they describe.
    # NOTHING HERE IS CACHED, for the reason `rings_count` gives -- a cached derived number is a
    # second truth an edit can contradict, and `_require_clean()` is the only gate needed.

    @property
    def carbon_count(self):
        """How many carbon atoms.  See `carbon_sp3_fraction`, whose denominator this is."""
        self._require_clean()
        return desc_carbon_count(self._structure)

    @property
    def carbon_sp3_count(self):
        """How many carbons store hybridization 1 -- sp3, and not an aromatic or allenic carbon."""
        self._require_clean()
        return desc_carbon_sp3_count(self._structure)

    @property
    def carbon_sp3_fraction(self):
        """`carbon_sp3_count / carbon_count`, and 0.0 for a molecule with no carbon.

        0.0 rather than nan or a refusal: this number goes into a descriptor vector where one nan
        poisons the row.
        """
        self._require_clean()
        cdef uint32_t total = desc_carbon_count(self._structure)
        if total == 0:
            return 0.0
        return <double> desc_carbon_sp3_count(self._structure) / <double> total

    @property
    def heteroatoms_count(self):
        """How many atoms are neither carbon, hydrogen, nor R (element 0).

        A count of ATOMS.  `heteroatoms_of(n)` is the per-atom count of heteroatom NEIGHBOURS and
        summing it answers a different question.
        """
        self._require_clean()
        return desc_heteroatoms_count(self._structure)

    @property
    def valence_electrons_count(self):
        """Sum of Zv - charge + implicit hydrogens over every atom.

        Zv is the group number convention stated in the header of `elements.tsv`.  Raises ValueError
        on an f-block atom, which states no count, and on an atom whose implicit hydrogen count is
        unknown, which makes the sum underivable -- run `chython.chemistry.calc_implicit` first.
        """
        self._require_clean()
        return desc_valence_electrons(self)

    @property
    def aromatic_rings_count(self):
        """How many rings of the basis have every bond stored order 4 -- `len(aromatic_rings)`.

        A kekulized molecule answers 0, which is the representation state and not a perception
        failure; `thiele()` is what changes it.
        """
        self._require_clean()
        cdef uint32_t counts[5]
        desc_ring_classes(self._structure, counts)
        return counts[0]

    @property
    def aliphatic_rings_count(self):
        """How many rings of the basis are not aromatic.

        The complement of `aromatic_rings_count`, so the two sum to `rings_count`.  NOT the same as
        `saturated_rings_count`: tetralin's carbocycle is aliphatic and unsaturated at once, because
        it shares an order-4 bond with the arene.
        """
        self._require_clean()
        cdef uint32_t counts[5]
        desc_ring_classes(self._structure, counts)
        return counts[1]

    @property
    def saturated_rings_count(self):
        """How many rings of the basis have every bond stored order 1."""
        self._require_clean()
        cdef uint32_t counts[5]
        desc_ring_classes(self._structure, counts)
        return counts[2]

    @property
    def heterocycles_count(self):
        """How many rings of the basis hold an atom that is neither carbon, hydrogen, nor R (element 0)."""
        self._require_clean()
        cdef uint32_t counts[5]
        desc_ring_classes(self._structure, counts)
        return counts[3]

    @property
    def aromatic_heterocycles_count(self):
        """How many rings of the basis are aromatic and heterocyclic at once."""
        self._require_clean()
        cdef uint32_t counts[5]
        desc_ring_classes(self._structure, counts)
        return counts[4]

    @property
    def spiro_atoms_count(self):
        """How many atoms are shared by two rings of the basis that share nothing else."""
        self._require_clean()
        cdef uint32_t counts[3]
        desc_ring_atoms(self._structure, counts)
        return counts[0]

    @property
    def bridgehead_atoms_count(self):
        """How many atoms bridge two rings of the basis that share a path of at least two bonds.

        A bridgehead carries at least three ring bonds, which is what excludes the middle of the
        shared path -- norbornane has two, not three.  Fused rings share a single bond and have none.
        """
        self._require_clean()
        cdef uint32_t counts[3]
        desc_ring_atoms(self._structure, counts)
        return counts[1]

    @property
    def fused_ring_systems_count(self):
        """How many connected components the subgraph of ring bonds has.

        An isolated ring is one system, so benzene answers 1 and biphenyl 2; a spiro atom merges its
        two rings into one.  0 for an acyclic molecule.
        """
        self._require_clean()
        cdef uint32_t counts[3]
        desc_ring_atoms(self._structure, counts)
        return counts[2]

    def eccentricities(self):
        """A `(n,)` int32 array of eccentricities: the largest distance from each atom to an atom it
        can reach.

        A METHOD, not a property, because it allocates an array on every call -- the same reason
        `distance_matrix` is one.  Indexed like `distance_matrix`'s rows: entry `i` belongs to
        `atom_numbers[i]`.  An atom that can reach nothing gets 0.
        """
        self._require_clean()
        return desc_eccentricities(self._structure)

    @property
    def wiener_index(self):
        """The Wiener index: the sum of topological distances over unordered pairs of atoms.

        Wiener, JACS 69 (1947) 17.  Defined on the hydrogen-suppressed graph -- explicit hydrogen
        atoms are vertices here and raise it.  A pair in different components is skipped, so the index
        is additive over components.
        """
        self._require_clean()
        return desc_wiener(self._structure)

    @property
    def graph_radius(self):
        """The smallest eccentricity.  0 for an empty molecule, and 0 for any molecule holding an
        isolated atom -- see `eccentricities`.  Use `split()` to ask per-component.
        """
        self._require_clean()
        cdef int32_t rd[2]
        desc_radius_diameter(self._structure, rd)
        return rd[0]

    @property
    def graph_diameter(self):
        """The largest eccentricity: the longest shortest path in the molecule.  For a disconnected
        molecule that is the widest component's diameter.
        """
        self._require_clean()
        cdef int32_t rd[2]
        desc_radius_diameter(self._structure, rd)
        return rd[1]

    def zagreb_index(self, uint32_t order=1):
        """The first Zagreb index (`order=1`, the default) or the second (`order=2`).

        M1 is the sum of squared degrees, M2 the sum over bonds of the degree product; Gutman and
        Trinajstic, Chem. Phys. Lett. 17 (1972) 535.  A METHOD because it takes the order, and no
        other order is defined -- 0 or 3 raises ValueError rather than answering.

        THIS IS THE ONLY PLACE THE 1-OR-2 DOMAIN IS STATED.  `desc_zagreb` takes a `bint`, so there is
        no third value for it to answer wrongly under `nogil`, where a refusal is impossible.
        """
        self._require_clean()
        if order != 1 and order != 2:
            raise ValueError('order must be 1 or 2; Gutman and Trinajstic define no others')
        return desc_zagreb(self._structure, order == 2)

    @property
    def randic_index(self):
        """The Randic branching index: sum over bonds of 1 / sqrt(deg(u) * deg(v)).

        Randic, JACS 97 (1975) 6609.  Equal to `chi(1)` -- the first-order connectivity index is the
        same sum -- and 0.0 for a molecule with no bonds.
        """
        self._require_clean()
        return desc_randic(self._structure)

    @property
    def balaban_j(self):
        """Balaban's average distance sum connectivity index J.

        `q / (mu + 1) * sum over bonds of 1 / sqrt(s(u) * s(v))`; Balaban, Chem. Phys. Lett. 89 (1982)
        399.  RAISES ValueError ON A DISCONNECTED MOLECULE -- a vertex distance sum needs a path to
        every atom -- and names `split()`, which yields parts that each answer.  A one-atom molecule is
        connected and answers 0.0.
        """
        self._require_clean()
        return desc_balaban_j(self)

    @property
    def bertz_ct(self):
        """Bertz's molecular complexity index CT.

        Bertz, JACS 103 (1981) 3599.  An information content over the molecule's CONNECTIONS -- pairs
        of bonds sharing an atom -- plus an element diversity term.  Two connections are equivalent
        when their central atoms share an `atoms_order` class and their outer atoms' classes agree as
        an unordered pair; that partition is chython's stated reading of the paper, which leaves it
        open, so this number is not comparable with another toolkit's CT.

        Reads symmetry, not bond orders, so benzene and cyclohexane agree.  Defined for a disconnected
        molecule, and not additive over its components.
        """
        self._require_clean()
        return desc_bertz_ct(self._structure)

    def chi(self, uint32_t order, bint valence=False):
        """The Kier-Hall connectivity index of the given order; `valence=True` for the delta-v variant.

        Order 0 sums 1/sqrt(delta) over atoms, order m sums 1/sqrt(the delta product) over simple paths
        of m bonds.  `chi(1)` is `randic_index`.  Kier and Hall, Rev. Comput. Chem. 2 (1991) 367-422.

        delta is the CSR row length -- what `degree_of()` reports -- so a dative bond is an edge and
        an explicit hydrogen vertex is a vertex, both counted; the SMARTS `D` primitive answers a
        different question.  delta-v is Zv - h where h is the IMPLICIT hydrogen count only: an explicit
        hydrogen is a vertex in this graph with its own delta-v of 1, so the plain and valence variants
        agree on a hydrogen-suppressed graph and part on one carrying explicit H atoms.  The formal
        charge is NOT subtracted.  An atom with delta 0 contributes nothing and no path through it
        contributes either -- 1/sqrt(0) is not a number.  `valence=True` raises ValueError on an
        f-block atom and on an unknown implicit hydrogen count; the plain variant never refuses.

        Orders 0 to 4 only: past that the path enumeration grows exponentially and no published index
        uses it.
        """
        self._require_clean()
        if order > 4:
            raise ValueError('order must be 0-4; Kier and Hall tabulate no higher connectivity index')
        return desc_chi_of(self, order, valence)

    def estate_intrinsic_states(self):
        """`(n,)` float64: Kier and Hall's intrinsic state per atom, in this molecule's atom order.

        `I = ((2/N)**2 * dv + 1) / d` with N the period, dv the valence delta and d the degree; Kier and
        Hall, Pharm. Res. 1990, 7, 801.  Public because it is the other number the paper prints, and
        because the sum of `estate_indices()` equals the sum of these.

        `nan` FOR AN ATOM WITH NO HEAVY NEIGHBOUR, never 0.0, which is an ordinary intrinsic state --
        `split()` yields parts whose atoms each have a neighbour if that is what the caller wants.
        Raises ValueError on an f-block atom and on an unknown implicit hydrogen count, the same two
        refusals `chi(valence=True)` has and for the same reason.
        """
        self._require_clean()
        return desc_estate_of(self, True)

    def estate_indices(self):
        """`(n,)` float64: the electrotopological state per atom, in this molecule's atom order.

        `S_i = I_i + sum_j (I_i - I_j) / (dist_ij + 1)**2` over every other atom, with I the intrinsic
        state; Kier and Hall, Pharm. Res. 1990, 7, 801.  Row `i` is `self.atom_numbers[i]`, the same
        indexing `atom_invariants` and `distance_matrix` use.

        `nan` FOR AN ATOM WITH NO HEAVY NEIGHBOUR, never 0.0, and such an atom is left out of every
        other atom's sum as well -- it has no intrinsic state to perturb with.  A pair in different
        components is skipped, so a salt's organic part answers what it answers alone.  Refuses what
        `estate_intrinsic_states` refuses.
        """
        self._require_clean()
        return desc_estate_of(self, False)

    @property
    def hall_kier_alpha(self):
        """The Hall-Kier alpha: the sum of each atom's covalent radius relative to an sp3 carbon's,
        less one.

        Hall and Kier, Rev. Comput. Chem. 2 (1991) 367-422.  0.0 for a saturated hydrocarbon, -0.78 for
        benzene.  An element the paper's table omits contributes 0.0 -- a reference, not a measured
        value -- so a metal adds nothing rather than making the index refuse.
        """
        self._require_clean()
        return desc_hall_kier_alpha(self._structure)

    def kappa(self, uint32_t order, bint alpha=False):
        """Kier's kappa shape index of order 1, 2 or 3; `alpha=True` applies the Hall-Kier correction.

        Each compares the molecule's count of `order`-bond paths against the counts of the extremal
        graphs of the same size, so it reads as a linearity.  kappa1 cannot see branching (it reads only
        the bond count); kappa2 is where n-butane's 3.0 parts from isobutane's 1.333.  Kier and Hall,
        Rev. Comput. Chem. 2 (1991) 367-422.

        `alpha=True` replaces n by n + `hall_kier_alpha` and P by P + `hall_kier_alpha`.  A molecule
        with no path of that length answers 0.0 -- there is no shape of that length to report.
        """
        self._require_clean()
        if order < 1 or order > 3:
            raise ValueError('order must be 1, 2 or 3; Kier defines three shape indices')
        return desc_kappa(self._structure, order, alpha)

    def union(self, MoleculeContainer other not None, bint remap=True):
        """Both molecules in one container, as separate components.

        `remap=True` renumbers `other`'s atoms above this molecule's highest number so the two
        cannot collide; `remap=False` refuses when any number is shared, rather than quietly
        merging two different atoms into one.

        NOTHING OF EITHER SIDE'S STEREO IS LOST.  This rebuilds `other` through `add_atom`/`add_bond`,
        which carries no stereo by itself, so `other`'s parities, stereo groups, coordinates and wedges
        are copied across explicitly -- without that, `mol1.union(mol2)` returns a molecule whose
        SECOND side has been quietly racemised.  Carrying a parity VERBATIM is
        valid here for the reason `split()` gives: a parity is a statement about a frame of
        neighbours, this cuts no bond, and the atoms are appended in slot order, so `other`'s
        ascending-neighbour order -- the order decision D1 stores the sign against -- is preserved
        by a monotone remap.  Nothing has to be re-based.

        OR AND AND GROUP IDS ARE RENUMBERED, because they are the one piece of stereo whose meaning
        is not local to an atom.  Both sides number from 1, so carrying `other`'s ids verbatim would
        merge its OR 1 with this molecule's OR 1 -- two independent mixtures declared to be one.
        Each of `other`'s groups gets an id this molecule does not spend, and the memberships are
        otherwise untouched; `canonical_stereo_groups` owns the ids a caller should compare.  ABS is
        one bucket rather than a numbered group and needs no renumbering.  When the two together
        would need more than the 63 ids one kind has, this raises rather than dropping a group.

        The result has coordinates when EITHER side does, which was already reachable (this molecule
        with them, `other` without) and is now symmetric.  It carries the LARGER of the two model
        counts, and in a model index one side does not have, that side's atoms sit at the origin.
        """
        self._require_clean()
        other._require_clean()
        cdef set mine = set(self._numbers)
        cdef set theirs = set(other._numbers)
        if not remap and (mine & theirs):
            raise ValueError('the two molecules share atom numbers %s; pass remap=True to '
                             'renumber the second' % sorted(mine & theirs))
        cdef MoleculeContainer mol = self.copy()
        cdef atom_t *src = other._structure.atoms()
        cdef uint32_t *ptr = csr_ptr(other._structure)
        cdef halfedge_t *edges = csr_edges(other._structure)
        # `other`'s arena is never the one being edited below -- `mol.edit()` copies on write, and it
        # copies MOL's arena -- so these four pointers survive the scope even when `other is self`.
        cdef uint8_t *osg = (structure_stereo_groups(other._structure)
                             if structure_has(other._structure, SEG_STEREO_GROUPS) else NULL)
        # TAKEN BEFORE THE EDIT, and the models are contiguous, so one base pointer plus
        # `model * o_atoms` indexes them all.
        cdef uint32_t my_models = structure_conformer_count(self._structure)
        cdef uint32_t their_models = structure_conformer_count(other._structure)
        cdef uint32_t o_atoms = other._structure.header.atom_count
        cdef xyz_t *oxyz = (structure_conformer_xyz(other._structure, 0)
                            if their_models else NULL)
        cdef xyz_t *p_xyz
        cdef uint32_t model
        cdef xy_t *oxy = (structure_xy(other._structure)
                          if structure_has(other._structure, SEG_XY) else NULL)
        cdef uint32_t i, k, kind, group
        cdef uint8_t p
        cdef dict new_of_slot = {}
        cdef dict back = {}
        # (kind, other's id) -> the id it is given here.  Built before the edit, from the groups this
        # molecule already spends, so a free id is free against both sides at once.
        cdef dict group_map = {}
        cdef set used = set()
        cdef object key
        if osg is not NULL:
            for key in mol.stereo_groups():
                used.add(key)
            for i in range(other._structure.header.atom_count):
                if not osg[i]:
                    continue
                kind = sg_kind(osg[i])
                if kind != 2 and kind != 3:
                    continue
                key = (<int> kind, <int> sg_group(osg[i]))
                if key in group_map:
                    continue
                for group in range(1, STEREO_GROUP_MAX + 1):
                    if (<int> kind, <int> group) not in used:
                        break
                else:
                    raise ValueError('the two molecules together need more than %d %s stereo groups'
                                     % (STEREO_GROUP_MAX, 'OR' if kind == 2 else 'AND'))
                used.add((<int> kind, <int> group))
                group_map[key] = <int> group
        with mol.edit():
            # THE RESULT CARRIES THE WIDER SIDE'S COUNT.  The copy already holds this molecule's
            # models; the extra ones exist for `other`'s atoms and leave this molecule's at the
            # origin, which is the answer a union with a flat partner already gives the other way
            # round.  The extras come first because `set_xyz` names a model that exists.
            for model in range(my_models, their_models):
                mol.add_conformer(ext_index=other.conformer(model).ext_index)
            for i in range(other._structure.header.atom_count):
                new_of_slot[i] = mol.add_atom(<int> src[i].element, charge=<int> src[i].charge,
                                            isotope=<int> src[i].isotope,
                                            radical=at_radical(&src[i]),
                                            map_number=<int> src[i].map_number,
                                            implicit_h=<int> at_implicit_h(&src[i]))
                back[new_of_slot[i]] = <uint32_t> other._numbers[i]
            for i in range(other._structure.header.atom_count):
                for k in range(ptr[i], ptr[i + 1]):
                    if edges[k].to > i:
                        mol.add_bond(<uint32_t> new_of_slot[i],
                                   <uint32_t> new_of_slot[edges[k].to],
                                   <int> edges[k].order)
            # AFTER THE BONDS, and that is the whole reason it is a second loop: a parity stated
            # before its frame exists is a sign about nothing, and `_apply`'s harvest would have
            # nothing to hold it against.  By here every direction is in place.
            for i in range(other._structure.header.atom_count):
                p = structure_parity_at(other._structure, i)
                if p:
                    mol.set_parity(<uint32_t> new_of_slot[i], <int> p)
                if osg is not NULL and osg[i]:
                    kind = sg_kind(osg[i])
                    group = <uint32_t> <int> group_map.get((<int> kind, <int> sg_group(osg[i])), 0)
                    mol.set_stereo_group(<uint32_t> new_of_slot[i], <int> kind, <int> group)
                if oxy is not NULL:
                    mol.set_xy(<uint32_t> new_of_slot[i], xy_read_x(oxy + i), xy_read_y(oxy + i))
                # CARRIED SEPARATELY FROM `oxy` AND NOT INSTEAD OF IT.  A union of a 3D molecule with
                # a 2D one gives a result whose depiction is complete and whose geometry is complete
                # only where the sources had one; that is the honest answer, and collapsing the two
                # segments into one would have made it unrepresentable.
                if oxyz is not NULL:
                    for model in range(their_models):
                        p_xyz = oxyz + <size_t> model * o_atoms + i
                        mol.set_xyz(<uint32_t> new_of_slot[i], xyz_read_x(p_xyz),
                                    xyz_read_y(p_xyz), xyz_read_z(p_xyz), model=model)
            for i in range(other._structure.header.atom_count):
                for k in range(ptr[i], ptr[i + 1]):
                    if edges[k].wedge:
                        mol.set_wedge(<uint32_t> new_of_slot[i],
                                      <uint32_t> new_of_slot[edges[k].to], <int> edges[k].wedge)
        if not remap:
            # numbers were checked disjoint above, so putting `other`'s back is injective
            mol.remap(back)
        return mol

    def __and__(self, other):
        """`mol & atoms` is the substructure on those atom numbers."""
        return self.substructure(other)

    def __sub__(self, other):
        """`mol - atoms` is the substructure on everything EXCEPT those atom numbers.

        Raises `ValueError` for a number this molecule does not have, because subtracting an atom that
        is not there is a caller bug that would otherwise return a wrong answer quietly.  Raises too
        when nothing would be left, for the same reason
        `substructure` does.
        """
        cdef set drop = set(other)
        self._require_clean()
        if drop - set(self._numbers):
            raise ValueError('invalid atom numbers %s' % sorted(drop - set(self._numbers)))
        return self.substructure(set(self._numbers) - drop)

    def __or__(self, other):
        """`mol | other` is the union of two molecules, renumbering the second where they collide."""
        if not isinstance(other, MoleculeContainer):
            return NotImplemented
        return self.union(<MoleculeContainer> other)

    def add_atom(self, element, *, int charge=0, int isotope=0, bint radical=False,
                 int map_number=0, implicit_h=None, bint stereo=False):
        """Append one atom and return its stable id.

        `implicit_h` takes three kinds of value and there are only TWO statements among them:

        * an int in 0..14 -- the record states this many implicit hydrogens.  0 states zero.
        * `H_UNKNOWN` -- the record does not say, and nothing can derive it.  Stored as the
          sentinel; `implicit_h_of` will answer None; `float(mol)` will leave the mass light.
        * `None` -- the DEFAULT, and it now stores `H_UNKNOWN` too.  Omitting the argument says
          nothing about hydrogens, and the sentinel is how the arena spells "nothing was said".

        THE DEFAULT MAY NOT BE A STORED ZERO.  A caller who omits the argument has not made a
        statement, so storing zero puts words in their mouth -- and leaves `kekule()` unable to tell a
        builder-made aromatic N that nobody has counted from one stated to carry no hydrogen, which is
        the difference between pyrrole and pyridine.  Every consumer that does arithmetic on the count
        asks `at_implicit_h_unknown` first, so a caller -- or a test -- that means zero says
        `implicit_h=0`.

        The value carries this and not a flag.  `at_h_pinned` looks like it records the same fact,
        but it does not survive `copy()`: the copy path passes `implicit_h=at_implicit_h(src)`
        unconditionally, so every copied atom comes out pinned whatever the source was.  A sentinel
        in the nibble copies as itself.
        """
        cdef uint32_t number = _to_atomic_number(element)
        cdef uint8_t r_index = 0
        if number == 0 and isinstance(element, str) and len(element) > 1:
            r_index = _parse_r_index(element)
        if charge < CHARGE_MIN or charge > CHARGE_MAX:
            raise ValueError('charge must be in %d..%d' % (CHARGE_MIN, CHARGE_MAX))
        if isotope < 0 or isotope > ISOTOPE_MAX:
            raise ValueError('isotope must be an absolute mass number, 0 for unset')
        if map_number < 0 or map_number > MAP_NUMBER_MAX:
            raise ValueError('map_number must be in 0..%d' % MAP_NUMBER_MAX)
        cdef int hydrogens = 0
        if implicit_h is not None:
            hydrogens = implicit_h
            if (hydrogens < 0 or hydrogens > H_IMPLICIT_MAX) and hydrogens != H_UNKNOWN:
                raise ValueError('implicit_h must be in 0..%d, H_UNKNOWN (%d) for a record that '
                                 'does not state one, or None to state nothing (also H_UNKNOWN)'
                                 % (H_IMPLICIT_MAX, H_UNKNOWN))

        if self._next_id >= 0xFFFFFFFF:
            raise OverflowError('stable id space is exhausted; ids are never reused')
        cdef uint32_t n = self._next_id
        self._append(OP_ADD_ATOM, n, 0, <int32_t> number)
        self._next_id += 1
        if r_index:
            self._append(OP_SET_R_INDEX, n, 0, r_index)
        if charge:
            self._append(OP_SET_CHARGE, n, 0, charge)
        if isotope:
            self._append(OP_SET_ISOTOPE, n, 0, isotope)
        if radical:
            self._append(OP_SET_RADICAL, n, 0, 1)
        if map_number:
            self._append(OP_SET_MAP_NUMBER, n, 0, map_number)
        if implicit_h is not None:
            self._append(OP_SET_HYDROGENS, n, 0, hydrogens)
        if stereo:
            self._append(OP_SET_STEREO, n, 0, 2)  # True -> known-odd (parity 2)
        self._maybe_apply()
        return n

    def delete_atom(self, uint32_t n):
        self._require(n)
        self._append(OP_DELETE_ATOM, n, 0, 0)
        self._maybe_apply()

    def add_bond(self, uint32_t n, uint32_t m, int order=1):
        self._require(n)
        self._require(m)
        if n == m:
            raise ValueError('self loops are not allowed')
        if order not in ALLOWED_ORDERS:
            raise ValueError(ALLOWED_ORDERS_MSG)
        if self._scope_depth == 0 and self._has_bond(n, m):
            raise ValueError(f'bond {n}-{m} already exists; use set_order to change its order')
        self._append(OP_ADD_BOND, n, m, order)
        self._maybe_apply()

    def delete_bond(self, uint32_t n, uint32_t m):
        self._require(n)
        self._require(m)
        if self._scope_depth == 0 and not self._has_bond(n, m):
            raise KeyError((n, m))
        self._append(OP_DELETE_BOND, n, m, 0)
        self._maybe_apply()

    def set_order(self, uint32_t n, uint32_t m, int order):
        self._require(n)
        self._require(m)
        if order not in ALLOWED_ORDERS:
            raise ValueError(ALLOWED_ORDERS_MSG)
        if self._scope_depth == 0 and not self._has_bond(n, m):
            raise KeyError((n, m))
        self._append(OP_SET_ORDER, n, m, order)
        self._maybe_apply()

    def set_element(self, uint32_t n, element):
        """Change what atom `n` IS, keeping its bonds, its charge and its place in the arena.

        WHY THE MUTATION SURFACE CARRIES THIS FIELD.  Without it the only way to turn a carbon into a
        nitrogen is to delete the atom and rebuild it with its bonds -- which allocates a NEW stable
        id, and so silently breaks every mapping the caller was holding.
        A reaction patcher is the first consumer (`_smirks_patch.pxi`) and it needs the id to survive,
        because the id is what pairs a product atom with the reactant atom it came from.

        `element` is what `add_atom` takes -- an atomic number in 1..118 or an element symbol.

        WHAT THIS DOES NOT TOUCH: the implicit hydrogen count.  A count derived for the old element
        is very probably wrong for the new one, and this layer does not derive counts -- `calc_implicit`
        does, from `chemistry`, which `core` cannot import.  So a caller that changes an element and
        cares about hydrogens writes the count itself, or hands the molecule to `calc_implicit`
        afterwards; nothing here will guess one, and nothing here will erase the one that is stored.

        Stored CIP descriptors ARE dropped, unlike on an isotope edit -- see the invalidation note in
        `_apply` for the argument.  Parities are not: a parity is a statement about a frame of
        neighbours, and the frame is exactly what an element change leaves alone.
        """
        self._require(n)
        cdef uint32_t number = _to_atomic_number(element)
        self._append(OP_SET_ELEMENT, n, 0, <int32_t> number)
        self._maybe_apply()

    def set_charge(self, uint32_t n, int charge):
        self._require(n)
        if charge < CHARGE_MIN or charge > CHARGE_MAX:
            raise ValueError('charge must be in %d..%d' % (CHARGE_MIN, CHARGE_MAX))
        self._append(OP_SET_CHARGE, n, 0, charge)
        self._maybe_apply()

    def set_isotope(self, uint32_t n, int isotope):
        self._require(n)
        if isotope < 0 or isotope > ISOTOPE_MAX:
            raise ValueError('isotope must be an absolute mass number, 0 for unset')
        self._append(OP_SET_ISOTOPE, n, 0, isotope)
        self._maybe_apply()

    def set_radical(self, uint32_t n, bint radical):
        self._require(n)
        self._append(OP_SET_RADICAL, n, 0, 1 if radical else 0)
        self._maybe_apply()

    def set_map_number(self, uint32_t n, int map_number):
        self._require(n)
        if map_number < 0 or map_number > MAP_NUMBER_MAX:
            raise ValueError('map_number must be in 0..%d' % MAP_NUMBER_MAX)
        self._append(OP_SET_MAP_NUMBER, n, 0, map_number)
        self._maybe_apply()

    def set_hydrogens(self, uint32_t n, int hydrogens):
        """Write the implicit hydrogen count of atom `n`, or `H_UNKNOWN` to record that the source
        did not state one.  0..14 are counts; 15 is the sentinel and is spelled `H_UNKNOWN`.

        There is no None here, unlike `add_atom`: this method exists to write a value, and a caller
        with no value to write has nothing to call.  Writing 0 states zero hydrogens; writing
        `H_UNKNOWN` states that the number is unknown.  Both are writes and neither is the other.
        """
        self._require(n)
        if (hydrogens < 0 or hydrogens > H_IMPLICIT_MAX) and hydrogens != H_UNKNOWN:
            raise ValueError('implicit_h must be in 0..%d, or H_UNKNOWN (%d) for a record that '
                             'does not state one' % (H_IMPLICIT_MAX, H_UNKNOWN))
        self._append(OP_SET_HYDROGENS, n, 0, hydrogens)
        self._maybe_apply()

    def set_r_index(self, uint32_t n, int index):
        """Set the R index of the R at `n`.  Refused at apply if `n` is not element 0."""
        self._require(n)
        if index < 0 or index > R_INDEX_MAX:
            raise ValueError('R index must be in 0-%d' % R_INDEX_MAX)
        self._append(OP_SET_R_INDEX, n, 0, index)
        self._maybe_apply()

    def set_stereo(self, uint32_t n, bint stereo):
        """Set the stereo flag on atom `n`.  A legacy API that writes a CONFIGURED parity:
        True -> parity 2 (known-odd), False -> parity 1 (known-even).  Use `set_parity` for
        the three-state write (0 = unset, 1 = even, 2 = odd)."""
        self._require(n)
        self._append(OP_SET_STEREO, n, 0, 2 if stereo else 1)
        self._maybe_apply()

    def set_parity(self, uint32_t n, int parity):
        """Set the three-state parity of atom `n`.

        `parity` must be 0 (unset -- no wedge drawn), 1 (even), or 2 (odd).
        The same `OP_SET_STEREO` journal op is used; the three values map directly to the
        apply handler's 0/1/2 dispatch.  Raises `ValueError` for any value outside 0..2.
        """
        self._require(n)
        if parity < 0 or parity > 2:
            raise ValueError('parity must be 0 (unset), 1 (even), or 2 (odd)')
        self._append(OP_SET_STEREO, n, 0, parity)
        self._maybe_apply()

    def request_parity(self):
        """Lay out the parity segment even though this edit states no parity.

        For a writer that states one AFTER the seal: a parity is read against a perceived frame, so a
        reader deriving one needs the sealed arena's CSR first, and by then the persistent block has
        been laid out and cannot grow.  A caller using `set_parity` needs nothing here.
        """
        self._append(OP_WANT_PARITY, 0, 0, 0)
        self._maybe_apply()

    def set_xy(self, uint32_t n, double x, double y):
        self._require(n)
        if not (-214748.0 <= x <= 214748.0) or not (-214748.0 <= y <= 214748.0):
            raise ValueError('coordinate out of fixed point range')
        cdef int32_t ix = <int32_t> round(x * XY_SCALE)
        cdef int32_t iy = <int32_t> round(y * XY_SCALE)
        self._append(OP_SET_XY, n, <uint32_t> ix, iy)
        self._maybe_apply()

    @property
    def has_coordinates(self):
        self._require_clean()
        return structure_has(self._structure, SEG_XY)

    def xy_of(self, uint32_t n):
        self._require_clean()
        if not structure_has(self._structure, SEG_XY):
            return None
        cdef xy_t *p = structure_xy(self._structure) + <uint32_t> self._index_of[n]
        return (xy_read_x(p), xy_read_y(p))

    def add_conformer(self, ext_index=None):
        """Append a model with every atom at the origin; return its index.

        `ext_index` is the file's own MODEL or frame number, stored VERBATIM and never interpreted;
        None is `CONF_NO_INDEX`, which is what a generated conformer states.  The returned index is
        the source count plus this session's adds, so a `set_xyz` naming it is stable for the rest of
        the scope.  A bare `set_xyz` on a molecule with no model appends one and counts as an add, so a
        session that fills model 0 that way and then calls this gets index 1.

        A session either adds or drops, never both: a drop shifts the models above it while an add
        counts from the unshifted source count, so allowing the pair would make a returned index mean
        one thing before the seal and another after.  Two sessions cost one respan each.
        """
        if self._conf_drops:
            raise ValueError('this edit session already dropped a conformer; an add would count '
                             'from a source count the drop has not shifted yet, so use a second '
                             'session')
        cdef uint32_t models = structure_conformer_count(self._structure) + self._conf_adds
        if models >= <uint32_t> CONF_MAX_MODELS:
            raise ValueError('a molecule holds at most %d conformers' % int(CONF_MAX_MODELS))
        cdef uint32_t ext = <uint32_t> CONF_NO_INDEX
        if ext_index is not None:
            if ext_index < 0 or ext_index > <uint32_t> CONF_EXT_INDEX_MAX:
                raise ValueError('ext_index %r is outside 0..%d; the value above that range is '
                                 'CONF_NO_INDEX, which is how "no number" is stored'
                                 % (ext_index, int(CONF_EXT_INDEX_MAX)))
            ext = <uint32_t> ext_index
        self._append(OP_ADD_CONFORMER, ext, 0, 0, 0, <uint16_t> models)
        self._conf_adds += 1
        self._maybe_apply()
        return models

    def drop_conformer(self, uint32_t index):
        """Drop model `index`; the seal compacts the survivors down.

        Indices are list positions in the molecule as it stands, so dropping model 1 of three leaves
        the old model 2 at index 1 after the seal, carrying its own `ext_index`.  A session that drops
        may neither add nor `set_xyz`: both name an index the compaction moves.  Dropping every model
        leaves no segment at all, which is how a molecule with no conformers is stored.
        """
        if self._conf_adds:
            raise ValueError('this edit session already added a conformer; a drop would shift the '
                             'index that add returned, so use a second session')
        if index >= structure_conformer_count(self._structure):
            raise IndexError(index)
        if self._conf_dropped is None:
            self._conf_dropped = set()
        if index in self._conf_dropped:
            raise ValueError('conformer %d is already dropped in this edit session' % int(index))
        self._conf_dropped.add(index)
        self._conf_drops += 1
        self._append(OP_DROP_CONFORMER, 0, 0, 0, 0, <uint16_t> index)
        self._maybe_apply()

    def set_xyz(self, uint32_t n, double x, double y, double z, uint32_t model=0):
        """Set atom `n`'s three-dimensional coordinates in model `model`.

        Independent of `set_xy`, deliberately: a reader of a 3D file calls BOTH, so the depiction and
        the geometry are two facts about one atom rather than one fact read two ways.  `clean2d()`
        rewrites the first and must never touch the second.

        `model` names a model that exists, counting this session's `add_conformer` calls; the one
        exception is `model=0` on a molecule that carries none, which appends the first model.

        The range check is on the SLOT and not on any file's field.  `int32_t` at `XY_SCALE` reaches
        +/-214748.0, which is far wider than the ten-column `F10.4` a V3000 line spends on a
        coordinate -- and RULES.md §1.4 says a slot's domain is the slot's, so a format that wants a
        narrower one enforces it at its own boundary.
        """
        self._require(n)
        if self._conf_drops:
            raise ValueError('this edit session dropped a conformer and the compaction moves every '
                             'index above it, so set coordinates in a second session')
        cdef uint32_t models = structure_conformer_count(self._structure) + self._conf_adds
        if models == 0 and model == 0:
            # THE IMPLICIT FIRST MODEL IS A REAL ADD, journalled here rather than derived at the seal: a
            # later `add_conformer` in this same session counts from `_conf_adds`, so a derived one
            # would have it return 0 and overwrite the coordinates this call is placing.
            self._append(OP_ADD_CONFORMER, <uint32_t> CONF_NO_INDEX, 0, 0, 0, 0)
            self._conf_adds += 1
        elif model >= models:
            raise ValueError('this molecule holds %d conformer(s), so model %d does not exist; '
                             'add_conformer() appends one' % (int(models), int(model)))
        if not (-214748.0 <= x <= 214748.0) or not (-214748.0 <= y <= 214748.0) \
                or not (-214748.0 <= z <= 214748.0):
            raise ValueError('coordinate out of fixed point range')
        cdef int32_t ix = <int32_t> round(x * XY_SCALE)
        cdef int32_t iy = <int32_t> round(y * XY_SCALE)
        cdef int32_t iz = <int32_t> round(z * XY_SCALE)
        self._append(OP_SET_XYZ, n, <uint32_t> ix, iy, iz, <uint16_t> model)
        self._maybe_apply()

    @property
    def has_3d(self):
        """True when the arena carries a conformer, exactly as `has_coordinates` answers for `SEG_XY`.

        Says nothing about whether the coordinates are any good -- a molecule whose every atom sits at
        the origin answers True, because the segment is there.  That is the same honesty
        `has_coordinates` offers and for the same reason: the arena knows what it stores, not what it
        means.
        """
        self._require_clean()
        return structure_conformer_count(self._structure) != 0

    def xyz_of(self, uint32_t n):
        """Atom `n`'s `(x, y, z)`, or None when the molecule carries no conformer at all.

        None is a statement about the MOLECULE and never about the atom: within a conformer every atom
        has a coordinate, because the segment is one dense column per model.  An atom added by an edit
        to a 3D molecule therefore reads (0.0, 0.0, 0.0) and not None -- see design D7, which records
        that as the slice's known defect.
        """
        self._require_clean()
        if not structure_conformer_count(self._structure):
            return None
        cdef xyz_t *p = structure_conformer_xyz(self._structure, 0) \
            + <uint32_t> self._index_of[n]
        return (xyz_read_x(p), xyz_read_y(p), xyz_read_z(p))

    def set_wedge(self, uint32_t narrow, uint32_t wide, int wedge):
        self._require(narrow)
        self._require(wide)
        if narrow == wide:
            raise ValueError('a wedge needs two distinct atoms')
        if wedge < 0 or wedge > 3:
            raise ValueError('wedge must be 0 none, 1 up, 2 down, 3 either')
        self._append(OP_SET_WEDGE, narrow, wide, wedge)
        self._maybe_apply()

    def wedge_of(self, uint32_t narrow, uint32_t wide):
        self._require_clean()
        cdef halfedge_t *he = csr_find(self._structure, <uint32_t> self._index_of[narrow],
                                       <uint32_t> self._index_of[wide])
        if he is NULL:
            raise KeyError((narrow, wide))
        return he.wedge

    def wedge_between(self, uint32_t n, uint32_t m):
        """The wedge drawn on the bond between `n` and `m` as `(narrow_id, code)`, or None.

        `wedge_of(narrow, wide)` is the DIRECTIONAL read and stays: it answers "is there a wedge
        pointing this way", which is what a round trip needs.  This one answers "is there a wedge, and
        which way does it point", which is what everything else needs -- and which three call sites
        were assembling by looking the pair up, then looking it up reversed.  Not spelled
        `wedge_from(n, m)`: a name reading "from n" that can answer "the narrow end is m" is a trap,
        and RULES.md 9.4 keeps `set_wedge`'s narrow-first pair for the same reason.

        Precondition: at most one of `set_wedge(n, m, ...)` and `set_wedge(m, n, ...)` is set.
        If both half-edges carry a wedge code (reachable from a caller but not from a CTfile reader),
        the `n`-forward half-edge wins.

        Raises KeyError((n, m)) when the two atoms are not bonded, as `wedge_of` and `bond_cip_of` do.
        """
        self._require_clean()
        cdef uint32_t i = <uint32_t> self._index_of[n]
        cdef uint32_t j = <uint32_t> self._index_of[m]
        cdef halfedge_t *forward = csr_find(self._structure, i, j)
        if forward is NULL:
            raise KeyError((n, m))
        if forward.wedge:
            return (n, forward.wedge)
        cdef halfedge_t *back = csr_find(self._structure, j, i)
        if back is not NULL and back.wedge:
            return (m, back.wedge)
        return None

    @property
    def aromatic_rings(self):
        """`rings` filtered to the rings whose every bond has order 4.

        A FILTER over the ring set and not a rename of `rings` -- which is why the chython 2 names
        block declines to alias it.  Kekulized molecules have no order-4 bonds and get an empty list,
        which is the honest answer rather than a perception fallback.
        """
        self._require_clean()
        cdef list out = []
        cdef tuple ring          # `rings` yields tuples, and the filter hands the same objects back
        cdef Py_ssize_t i, size
        cdef uint32_t prev, cur
        cdef halfedge_t *e
        for ring in self.rings:
            size = len(ring)
            # The closing bond walked first, then every other one -- the wrap-around without a
            # modulo.  Spelled with a carried `prev` and NOT as `ring[i - 1]`: this translation unit
            # compiles with `wraparound=False` (`_core.pyx:8`), so `ring[-1]` is an unchecked read one
            # slot before the tuple rather than its last element.  It segfaults, and RULES.md 7.6 is
            # the family -- nothing in the pipeline objects.
            prev = <uint32_t> self._index_of[ring[size - 1]]
            for i in range(size):
                cur = <uint32_t> self._index_of[ring[i]]
                e = csr_find(self._structure, prev, cur)
                if e is NULL or e.order != 4:
                    break
                prev = cur
            else:
                out.append(ring)
        return out

    @property
    def has_layout(self):
        """True when the stored coordinates are a plane something can be drawn from.

        `has_coordinates` answers only that the arena carries an XY segment, and a molecule that got a
        segment from a writer has one with every atom at the origin.  The test a renderer needs is
        whether the plane is non-degenerate -- a bounding box with span in at least one axis -- so it
        lives here rather than at each of the five call sites, one of them per atom.  The threshold is
        .01 molecule units, in the arena's fixed point.

        One atom needs no layout, so a one-atom molecule with a segment is True; an empty one is False.
        """
        self._require_clean()
        if not structure_has(self._structure, SEG_XY):
            return False
        cdef uint32_t count = self._structure.header.atom_count
        if count == 0:
            return False
        if count == 1:
            return True
        cdef xy_t *xy = structure_xy(self._structure)
        cdef xy_t *p = xy
        # int32_t bounds: accumulate as int32 (no overflow on an individual coordinate), but
        # compute the span in int64_t — two int32 values at opposite ends of the int32 range
        # differ by up to 2×INT32_MAX ≈ 4.3×10⁹, which overflows a signed int32 subtraction.
        cdef int32_t min_x = p.x, max_x = p.x, min_y = p.y, max_y = p.y
        cdef uint32_t i
        cdef int64_t span_x, span_y
        for i in range(1, count):
            p = xy + i
            if p.x < min_x:
                min_x = p.x
            elif p.x > max_x:
                max_x = p.x
            if p.y < min_y:
                min_y = p.y
            elif p.y > max_y:
                max_y = p.y
        span_x = <int64_t> max_x - <int64_t> min_x
        span_y = <int64_t> max_y - <int64_t> min_y
        return span_x >= 100 or span_y >= 100

    def coordinates(self):
        """`{n: (x, y)}` for the whole plane, or `{}` when none is stated.

        Eight call sites built this dict one `xy_of` at a time, each paying a segment check, a dict
        lookup, two divisions and a tuple.  Keys come out in arena order, so the result can be zipped
        against `iter(mol)`.

        `{}` rather than a plane of zeros when the segment is absent, because "no coordinates" and
        "every atom at the origin" are different facts and the arena distinguishes them.  A caller that
        wants origins writes `dict.fromkeys(mol, (0., 0.))` and says so.
        """
        self._require_clean()
        cdef dict out = {}
        if not structure_has(self._structure, SEG_XY):
            return out
        cdef xy_t *xy = structure_xy(self._structure)
        cdef xy_t *p
        cdef list ids = self._numbers
        cdef uint32_t i
        for i in range(self._structure.header.atom_count):
            p = xy + i
            out[ids[i]] = (xy_read_x(p), xy_read_y(p))
        return out

    def set_atom_cip(self, uint32_t n, object descriptor):
        """State the CIP descriptor of atom `n`: 'R' 'S' 'r' 's' 'M' 'P' 'm' 'p', or None to clear.

        STORAGE ONLY.  Nothing here computes, checks or corrects a descriptor -- this records what an
        input said.  A descriptor that disagrees with the structure is stored exactly as given, in
        keeping with the arena's posture everywhere else: the caller can see what it holds, and only an
        explicit operation changes it.

        Lowercase r/s are the pseudo-asymmetric descriptors and are NOT accepted as spellings of R/S.
        """
        self._require(n)
        self._append(OP_SET_ATOM_CIP, n, 0,
                     <int32_t> _cip_code(descriptor, _ATOM_CIP_BY_NAME, ATOM_CIP_CODES, 'an atom'))
        self._maybe_apply()

    def atom_cip_of(self, uint32_t n):
        return ATOM_CIP_CODES[at_cip(self._atom(n))]

    def set_bond_cip(self, uint32_t n, uint32_t m, object descriptor):
        """State the CIP descriptor of the bond between `n` and `m`: 'E' 'Z' 'M' 'P', or None to clear.

        Not directional: `set_bond_cip(a, b, 'E')` and `set_bond_cip(b, a, 'E')` are the same statement,
        and reading it back from either end gives the same answer.  See the re-apply loop in `_apply`,
        which is the single site that writes both half-edges.
        """
        self._require(n)
        self._require(m)
        if n == m:
            raise ValueError('a bond needs two distinct atoms')
        self._append(OP_SET_BOND_CIP, n, m,
                     <int32_t> _cip_code(descriptor, _BOND_CIP_BY_NAME, BOND_CIP_CODES, 'a bond'))
        self._maybe_apply()

    def bond_cip_of(self, uint32_t n, uint32_t m):
        self._require_clean()
        cdef halfedge_t *he = csr_find(self._structure, <uint32_t> self._index_of[n],
                                       <uint32_t> self._index_of[m])
        if he is NULL:
            raise KeyError((n, m))
        return BOND_CIP_CODES[he_cip(he)]

    @property
    def cip_log(self):
        """Descriptors dropped because the molecule changed, newest last.  See `_apply`'s drop rule.

        A DROP IS ONLY RECOVERABLE FROM HERE.  Storage cannot tell an atom that never carried a
        descriptor from one whose descriptor was dropped -- both hold code 0 -- so if the event is not
        read from this list it is not reachable from the bytes afterwards.  A READ-ONLY VIEW over `log`,
        filtered to the `edit:cip` stage -- a tuple, so a caller cannot edit the record of a loss out of
        the container that suffered it.
        """
        if self._log is None:
            return ()
        # an explicit loop and not a genexpr: a comprehension's loop variable is a name Cython never saw
        # declared, and `warn.undeclared` reports it
        cdef object rec
        cdef list out = []
        for rec in self._log.by_stage('edit:cip'):
            out.append(str(rec))
        return tuple(out)

    def atom_cips(self):
        """`{n: descriptor}` for the atoms carrying one.  Absent means no descriptor stated."""
        self._require_clean()
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t i
        cdef uint8_t code
        cdef dict out = {}
        for i in range(self._structure.header.atom_count):
            code = at_cip(&atoms[i])
            if code:
                out[self._numbers[i]] = ATOM_CIP_CODES[code]
        return out

    def bond_cips(self):
        """`{(n, m): descriptor}` for the bonds carrying one, each pair once."""
        self._require_clean()
        cdef uint32_t *ptr = csr_ptr(self._structure)
        cdef halfedge_t *edges = csr_edges(self._structure)
        cdef uint32_t i, k
        cdef uint8_t code
        cdef dict out = {}
        for i in range(self._structure.header.atom_count):
            for k in range(ptr[i], ptr[i + 1]):
                if edges[k].to > i:
                    code = he_cip(&edges[k])
                    if code:
                        out[(self._numbers[i], self._numbers[edges[k].to])] = \
                            BOND_CIP_CODES[code]
        return out

    def wedges(self):
        self._require_clean()
        cdef uint32_t *ptr = csr_ptr(self._structure)
        cdef halfedge_t *edges = csr_edges(self._structure)
        cdef halfedge_t *e
        cdef uint32_t i, k
        cdef list out = []
        for i in range(self._structure.header.atom_count):
            for k in range(ptr[i], ptr[i + 1]):
                e = &edges[k]
                if e.wedge:
                    out.append((self._numbers[i], self._numbers[e.to], e.wedge))
        return out

    def set_stereo_group(self, uint32_t n, int kind, int group=0):
        self._require(n)
        if kind < 0 or kind > 3:
            raise ValueError('kind must be 0 unspecified, 1 abs, 2 or, 3 and')
        if kind == 2 or kind == 3:
            if group < 1 or group > STEREO_GROUP_MAX:
                raise ValueError('OR and AND groups must have a group id in 1..%d' % STEREO_GROUP_MAX)
        else:
            group = 0
        self._append(OP_SET_STEREO_GROUP, n, <uint32_t> group, kind)
        self._maybe_apply()

    @property
    def has_stereo_groups(self):
        self._require_clean()
        return structure_has(self._structure, SEG_STEREO_GROUPS)

    def stereo_group_of(self, uint32_t n):
        self._require_clean()
        cdef uint32_t i = <uint32_t> self._index_of[n]
        if not structure_has(self._structure, SEG_STEREO_GROUPS):
            return (STEREO_UNSPECIFIED, 0)
        cdef uint8_t v = structure_stereo_groups(self._structure)[i]
        return (sg_kind(v), sg_group(v))

    def stereo_groups(self):
        self._require_clean()
        if not structure_has(self._structure, SEG_STEREO_GROUPS):
            return {}
        cdef uint8_t *sg = structure_stereo_groups(self._structure)
        cdef uint32_t i
        cdef dict out = {}
        for i in range(self._structure.header.atom_count):
            if sg[i]:
                out.setdefault((sg_kind(sg[i]), sg_group(sg[i])), []).append(
                    self._numbers[i])
        return out

    def canonical_stereo_groups(self):
        """`stereo_groups()` with the opaque stored ids replaced by canonical ones.

        Returns {(kind, canonical_id): [n, ...]}, the same shape `stereo_groups()` returns
        and with the same memberships -- only the second half of each key moves.  Two molecules that
        are the same molecule with the same stereo groups get the same dict here whatever ids their
        input files happened to use AND whatever order their atoms were created in, UP TO the
        permutation `canonical_stereo_group_ambiguities()` reports: where that tuple is empty -- the
        common case, and the only case in which this dict may be compared or hashed key by key -- the
        two dicts are equal outright, and where it is not, the ids inside one of its classes may be
        exchanged between the two readings.  Measured identical in processes started with different
        PYTHONHASHSEED, since nothing on the path from the arena to the ids is Python-hash-ordered.
        It is NOT a persistence format: the ids ride the refinement class numbering, so a future
        change to the atom invariant may renumber them.

        OR and AND groups are numbered densely from 1, each kind independently, ordered by member
        count first and then by the members' symmetry classes -- so a two-member group precedes a
        three-member one -- with the canonical atom order breaking ties only between groups nothing
        else can separate.  ABS keeps its stored number (0), because ABS is one bucket rather than a
        numbered group -- `set_stereo_group` will not accept anything else for it -- so every
        canonical id this returns is a legal input to `set_stereo_group`.

        WHERE THE MOLECULE'S OWN SYMMETRY CAN EXCHANGE TWO GROUPS, which id each one got is
        arbitrary and `canonical_stereo_group_ambiguities()` says so -- read it before hashing or
        comparing this dict key by key (ruling F89).  That report errs in ONE direction: it may call a
        pinned pair ambiguous, it never calls an exchangeable pair pinned, so an empty tuple is a
        guarantee and a non-empty one is an upper bound on what moves.

        COMPARE MEMBERS AS GROUP-RELATIVE DATA -- their count, their parities READ IN A FRAME YOU
        NAME, their symmetry classes -- and NOT as `canonical_order()` positions.  The frame is part
        of the advice and not a refinement of it: `parity_of` returns the STORED byte, whose reference
        order is the order the atoms were created in (ruling F26), so it is not a property of the
        molecule, and two encodings of one molecule differ in it legitimately.  Take a parity through
        `translate_stereo(atom, refs)` with `refs` you chose yourself -- for a ring centre, say, the
        next ring atom, the previous one, the substituent, the hydrogen.  Measured on the OR pair
        {C1,C2}, {C3,C4} of 1,2,3,4-tetrachlorocyclobutane with ring-frame parities 1, 2, 2, 2: a gate
        joining members on `(count, sorted(parity_of(...)))` sees SIX different shapes over the 96
        creation orders swept below, and the same join through `translate_stereo` in the ring frame
        sees one.

        `canonical_order()` positions fail for a second reason -- that order is seeded stereo-blind,
        so it moves with the creation order even where these ids do not.  On the same ring with a
        three-member and a one-member OR group this method and `canonical_stereo_group_ambiguities()`
        certify every id pinned, and joining the two reads takes twelve distinct values over the 96
        creation orders (every permutation of the four ring carbons, each with the chlorines created
        before them, after them, and interleaved with them; 72 of the 96 distinct).

        A READ (ruling F79).  Nothing is renumbered in the arena, so a `copy()` that shares it is
        unaffected and `stereo_groups()` keeps reporting the stored ids.  Raises
        `AutomorphismBudgetExceeded` for the same reason `canonical_order` does, and by the same
        route: there is no degraded canonical id.
        """
        self._require_clean()
        if not structure_has(self._structure, SEG_STEREO_GROUPS):
            return {}
        cdef uint8_t ids[256]
        canonical_stereo_group_ids(self._structure, ids, NULL)
        cdef uint8_t *sg = structure_stereo_groups(self._structure)
        cdef uint32_t i
        cdef dict out = {}
        for i in range(self._structure.header.atom_count):
            if sg[i]:
                out.setdefault((sg_kind(sg[i]), ids[sg[i]]), []).append(self._numbers[i])
        return out

    def canonical_stereo_group_ambiguities(self):
        """Which keys of `canonical_stereo_groups()` may be interchangeable under this molecule's
        symmetry.

        Returns a tuple of frozensets of `(kind, canonical_id)` keys.  Keys in one frozenset belong to
        groups that nothing invariant here tells apart -- same kind, same size, same everything the
        rules below can see -- so the ids inside a frozenset are handed out in an arbitrary order while
        the SET of member sets they cover is not.  That is a weaker statement than "the molecule maps
        one onto the other", and deliberately so: see the next paragraph.  A tuple with no entries, the
        common case, means every id is pinned and the whole view may be compared key by key.

        THE REPORT ERRS IN ONE DIRECTION ONLY, and a caller may rely on that.  Two groups are reported
        together when no invariant rule available here can tell them apart, which is a weaker test than
        "some symmetry of the molecule exchanges them" -- colour refinement is incomplete, so a pair
        that is really pinned can land in a class.  The converse cannot happen: a pair some symmetry
        exchanges is invisible to every invariant rule and so is never split out.  An empty tuple is
        therefore a guarantee that the whole view compares, and a non-empty one is an upper bound on
        what may move -- never a claim that it does.

        The ids in one frozenset are consecutive, and a key outside every frozenset is pinned no
        matter how many groups are tied elsewhere -- an ambiguity cannot move an id it does not
        cover.  So a caller comparing two molecules key by key needs to special-case only these
        classes, and each of them only as an unordered set.

        THE TUPLE ITSELF COMPARES: the classes are ordered by their smallest `(kind, canonical_id)`,
        which no ambiguity can move because the ambiguity is confined inside a class (ruling F92), so
        `==` between two forms of one molecule is True and does not need `set()` around it.  As with
        `canonical_stereo_groups()`, compare the members as group-relative data -- and read a parity in
        a frame you name, never as `parity_of`'s stored byte -- and NOT as `canonical_order()`
        positions: that order is stereo-blind, and a molecule this method certifies pinned can still
        take twelve different position joins over the 96 creation orders named there.

        Ruling F89.  1,2,3,4-tetrachlorocyclobutane with its centres in OR groups (1, 1) and (2, 2) is
        the smallest witness: read each centre's parity in the RING FRAME -- the next ring atom, the
        previous one, the chlorine, the hydrogen, which is the frame `translate_stereo()` will give you
        and NOT the stored byte, whose frame follows the creation order -- and let those parities
        alternate 1, 2, 1, 2.  The ring's rotation by two then carries one group onto the other and
        each parity onto an equal one (a rotation preserves the ring frame, so it preserves parities),
        so both id assignments describe the same mixture.  This does not merge them -- two OR groups
        describe four stereoisomers where one describes two -- and it does not raise: a caller comparing
        two molecules should compare the keys outside these classes directly, and each class only as
        its unordered set of member sets.

        A read, like `canonical_stereo_groups()`, and it does the same work: call one or the other,
        not both, if the cost matters.  Raises `AutomorphismBudgetExceeded` on the same path.
        """
        self._require_clean()
        if not structure_has(self._structure, SEG_STEREO_GROUPS):
            return ()
        cdef uint8_t ids[256]
        cdef uint8_t amb[256]
        canonical_stereo_group_ids(self._structure, ids, amb)
        cdef uint8_t *sg = structure_stereo_groups(self._structure)
        cdef uint32_t i
        cdef dict classes = {}
        cdef list out = []
        cdef object number
        for i in range(self._structure.header.atom_count):
            if sg[i] and amb[sg[i]]:
                classes.setdefault(amb[sg[i]], set()).add((sg_kind(sg[i]), ids[sg[i]]))
        # The C side numbers the classes by ascending smallest canonical id (ruling F92), which is
        # invariant because an ambiguity permutes ids only inside its own class and so cannot move a
        # class's minimum -- so this sort is the invariant order and the tuple compares with `==`.
        for number in sorted(classes):
            out.append(frozenset(classes[number]))
        return tuple(out)

    cdef dict _unit_dict(self, stereo_unit_t *u):
        cdef list numbers = self._numbers
        cdef uint32_t k, r
        cdef list refs = []
        for k in range(4):
            r = u.refs[k]
            refs.append(None if r == SU_NO_REF else numbers[r])
        cdef int live_parity
        # Parity is read from SEG_PARITY at the anchor's slot.
        live_parity = structure_parity_at(self._structure, u.anchor)
        return {'kind': u.kind, 'parity': live_parity, 'n_refs': u.n_refs,
                'anchor': numbers[u.anchor], 'refs': tuple(refs),
                # The high nibble only, not the flag bits sharing the byte with it.  A MASK of which
                # `refs` slots hold an unnamed direction, not a count (ruling F41) -- renamed from
                # `unnamed_directions` so that no consumer reads the new value as the old one.
                # `bin(mask).count('1')` is the count.
                'unnamed_mask': u.spare >> SU_UNNAMED_SHIFT,
                # The flag nibble's answer to "can this unit hold two configurations": exact, set by
                # `mark_stereogenic` while the table is being built, so it is never absent and never
                # stale.  False here means no molecule information is lost by leaving the unit
                # unconfigured.
                'stereogenic': (u.spare & SU_STEREOGENIC) != 0,
                }

    def stereo_units(self):
        """Every place in this molecule that COULD carry a configuration, as a list of dicts.

        A unit is `{kind, parity, n_refs, anchor, refs, unnamed_mask, stereogenic}`.  `anchor` is the stable id
        the parity is stored against.  `refs` is always a 4-tuple naming the anchor's directions,
        laid out per kind: an ATOM kind (tetrahedral) packs four directions in one order -- heavy
        neighbours in CSR ascending order, then explicit hydrogens ascending, then `None` for each
        direction with no atom of its own -- while a BOND kind (cis/trans, allene, atropisomer) packs
        two such pairs, the anchor's end first, so a `None` can appear in the middle and the entries
        need not ascend across the pair boundary.  Within every list a hydrogen direction sits after
        every heavy one, so the tuple's shape does not change when a hydrogen starts or stops being
        drawn and a stored parity is never re-based.

        `unnamed_mask` says WHICH slots those `None`s are directions rather than absences: bit i set
        means `refs[i]` is a direction with no atom of its own (an implicit hydrogen, or a sulfur
        lone pair), bit i clear with `refs[i] is None` means slot i holds no direction at all.
        `CH3CH=NOH` is why that distinction is in the record: its `refs` are `(CH3, None, O, None)`
        and only slot 1 is a direction, so a scalar count could not tell it from a molecule whose
        unnamed direction is on the other terminal.  The popcount is the count, and two unnamed
        directions in ONE direction list mean that list cannot be told apart, which is how
        `stereogenic_units()` rejects a methyl group -- `CH2=CHCH3` (mask `0b1011`) against
        `CH3CH=CHCH3` (mask `0b1010`) is exactly that test, and the reason the mask has to be per
        slot rather than a count.  For a cumulene the first of those two never reaches a record,
        because perception already refuses a terminal all of whose directions are unnamed; for an
        atom kind, and for whatever kinds later tasks add, the count would be ambiguous.

        Being in this list is not being a stereocentre; `stereogenic` is.  Every record here is a
        CANDIDATE, and the flag says whether flipping it would give a different molecule -- so
        `stereo_units()` is the perception layer's answer and `stereogenic_units()` is the chemist's.
        `parity` is whatever has been stored, which for an unconfigured molecule is 0; a stereogenic
        unit with parity 0 is a real stereocentre nobody has stated the configuration of.
        """
        self._require_clean()
        ensure_stereo_units(self._structure)
        # both taken after ensure_stereo_units, which reallocates the arena
        cdef uint32_t count = structure_stereo_unit_count(self._structure)
        cdef stereo_unit_t *units = structure_stereo_units(self._structure)
        cdef uint32_t k
        cdef list out = []
        for k in range(count):
            out.append(self._unit_dict(&units[k]))
        return out

    def stereogenic_units(self):
        """The subset of `stereo_units()` that can really hold two configurations.

        The same dicts, filtered on `stereogenic`.  This is the list a stereo-aware consumer wants:
        every record in it names a place where the molecule's identity depends on a configuration,
        whether or not one has been stated.  Which of them HAVE been stated is `parity != 0`, and the
        unsigned subset is a comprehension over this list rather than an accessor of its own -- one
        filter is a filter, two are an API to keep in step.
        """
        cdef list out = []
        cdef dict unit
        for unit in self.stereo_units():
            if unit['stereogenic']:
                out.append(unit)
        return out

    def chiral_atoms(self):
        """`{n: unit}` for every stereogenic ATOM-kind unit.

        Tetrahedral centres and allene/cumulene axes, which are named on their central atom (spec
        3.2) and so belong with the atoms rather than with the bonds.  The value is the same dict
        `stereo_units()` yields, so `mol.chiral_atoms()[n]['parity']` is a one-liner and the keys
        alone are the set of centres.

        WHICH QUESTION THIS ANSWERS, in the terms other toolkits use.  It is RDKit's
        `FindPotentialStereo`: every site whose configuration this record's identity depends on,
        labelled or not.  It is NOT `FindMolChiralCenters(includeUnassigned=False)`, which returns
        only the sites that carry a label, and it is NOT V2's `__chiral_centers`, which returns the
        stereogenic sites MINUS the labelled ones.  The one place it parts company with
        `FindPotentialStereo` is the
        epic's per-record convention: this predicate asks about THIS record, so a pseudo-asymmetric
        centre whose neighbours are unlabelled is not stereogenic here and is "potential" there --
        `test_trimethylcyclohexane_answers_per_diastereomer` is that difference, measured.
        `parity == 0` in the value means the configuration is NOT CONFIGURED (ruling F54); it never
        means "no wedge was drawn", which is a question the core does not yet answer at all.

        NOT CACHED, and deliberately.  The value depends on the stereo unit table, which is a lazy
        derived segment that any edit invalidates -- and an edit is exactly what a caller does
        between two reads of this.  A `cached_property` here would have to be dropped by every
        mutating path in the container, which is a list nobody can keep complete; recomputing is a
        table scan over a segment that is already built, so the cache would buy a scan and cost a
        class of stale-answer bug.
        """
        self._require_clean()
        ensure_stereo_units(self._structure)
        cdef uint32_t count = structure_stereo_unit_count(self._structure)
        cdef stereo_unit_t *units = structure_stereo_units(self._structure)
        cdef list numbers = self._numbers
        cdef uint32_t k
        cdef dict out = {}
        for k in range(count):
            if not (units[k].spare & SU_STEREOGENIC):
                continue
            if units[k].kind == SU_TETRA or units[k].kind == SU_ALLENE:
                out[numbers[units[k].anchor]] = self._unit_dict(&units[k])
        return out

    def chiral_bonds(self):
        """Stereogenic BOND-kind units, as a dict keyed on the bond.

        KEYED ON THE BOND, not on the anchor: `{(n, m): unit}` with the pair sorted.  Which
        end anchors a cis/trans or atropisomer unit is a function of slot order, not of chemistry --
        perception takes the lower-indexed end and relocates when that end is already claimed
        (ruling F45) -- so an anchor-keyed answer would move when the same molecule is read with its
        atoms in another order.  The bond does not move.  For a cumulene longer than one bond the key
        is the pair of TERMINALS, which is what the unit is named on.

        Not cached, for the reason `chiral_atoms` gives.  A fresh dict every call, so a caller may
        mutate the result.  It answers the same question `chiral_atoms` does -- RDKit's
        `FindPotentialStereo`, for the bond kinds -- with the same per-record convention and the same
        reading of `parity == 0` as "not configured" (ruling F54).
        """
        self._require_clean()
        ensure_stereo_units(self._structure)
        cdef uint32_t count = structure_stereo_unit_count(self._structure)
        cdef stereo_unit_t *units = structure_stereo_units(self._structure)
        cdef list numbers = self._numbers
        cdef uint32_t k, partner
        cdef dict out = {}
        for k in range(count):
            if not (units[k].spare & SU_STEREOGENIC):
                continue
            partner = stereo_unit_partner(self._structure, &units[k])
            if partner == SU_NO_REF:
                continue
            out[tuple(sorted((numbers[units[k].anchor], numbers[partner])))] = self._unit_dict(&units[k])
        return out

    def is_chiral(self, uint32_t n):
        """Does the atom `n` anchor a stereogenic unit?

        TRUE FOR A LABELLED SITE AS MUCH AS AN UNLABELLED ONE.  It answers "is this site
        stereogenic" and nothing else.  V2's `__chiral_centers` answers the narrower "which sites still
        need a sign", which is a comprehension over this one:
        `[s for s in mol.chiral_atoms() if mol.parity_of(s) == 0]`.  And `parity_of(s) == 0` there
        means the site is NOT CONFIGURED (ruling F54) -- not that no wedge was drawn on it, which the
        core cannot yet be asked, since deriving parity from wedges is a later task.

        Anchored, so this is the atom-side question: a cis/trans unit answers True at whichever
        terminal happens to be keyed in SEG_PARITY, which is why `chiral_bonds()` keys on the bond
        instead.  An atom that is only a DIRECTION of some unit anchors none and answers False.
        """
        self._require_clean()
        if n not in self._index_of:
            raise KeyError(n)
        ensure_stereo_units(self._structure)
        cdef stereo_unit_t *u = stereo_unit_of(self._structure,
                                               <uint32_t> self._index_of[n])
        return u is not NULL and (u.spare & SU_STEREOGENIC) != 0

    @property
    def stereo_truncated(self):
        """Was this molecule's stereogenicity marking taken conservatively?

        True when the symmetry search that decides which candidate units are really stereogenic ran
        out of budget, so every unit it had not settled was marked rather than dropped.  Then
        `chiral_atoms`, `chiral_bonds`, `is_chiral`, `stereo_units` and `stereogenic_units` are an
        OVER-approximation: a site may be reported whose refuting symmetry was never reached, and no
        site the molecule really has can be missing.  That is the sound direction for "could this be a
        stereocentre", which is why these readers answer instead of raising (ruling F62) -- and it is
        the opposite choice from `automorphism_orbits`, which raises, because a truncated canonical
        labelling is not a conservative answer but a wrong one.

        PER MOLECULE, not per unit: the flag cannot tell you WHICH marks are unproven, because the
        search abandons the record as a whole once the node budget is gone.  A caller that needs
        certainty has to re-ask a smaller record -- in practice one component at a time, since it is
        multi-component symmetry that exhausts the budget.

        Reads the stereo unit table, so it builds it on first call like every other stereo accessor,
        and an edit that changes the record clears both.  It never raises for truncation, which is the
        whole point of it.
        """
        self._require_clean()
        ensure_stereo_units(self._structure)
        return structure_stereo_truncated(self._structure)

    def unit_of(self, uint32_t n):
        """The stereo unit anchored at `n`, as `stereo_units()` describes it, or None.

        An atom anchors at most one unit -- see the invariant `_stereo.pxi` asserts -- so this is
        a function of the atom and not a choice among candidates.  An atom that is only a
        DIRECTION of some unit has none of its own and returns None.
        """
        self._require_clean()
        if n not in self._index_of:
            raise KeyError(n)
        ensure_stereo_units(self._structure)
        cdef stereo_unit_t *u = stereo_unit_of(self._structure,
                                               <uint32_t> self._index_of[n])
        return None if u is NULL else self._unit_dict(u)

    def validate_stereo(self):
        """The stated parities this molecule cannot justify, ASCENDING BY STABLE ID -- and clear them.

        A parity is accepted into the arena unconditionally when it is stated (`set_parity` checks the
        range 0..2 and nothing else), because the centre that justifies it need not exist yet.  This is
        where the question is asked, once, on a molecule that is finished -- asking it inside
        `add_atom_stereo` instead loses a configuration whose centre the next bond creates.  The answer
        is a REPORT rather than a raise, because a
        container mid-edit legitimately carries a sign nothing justifies yet and raising would make it
        unusable while it is being built.  A consumer that needs strictness calls this and acts on a
        non-empty result.

        Two shapes reach the list: an atom that anchors no unit at all, and one that anchors a unit
        `mark_stereogenic` did not mark.  Both are readable -- the bits mean what they always meant --
        which is why nothing before this point clears them (ruling F66).  The apply drops a sign
        whose FRAME was destroyed, because that one is no longer readable at all; this one preserves a
        sign whose frame has not yet existed, and only a caller asking gets it removed.

        A TRUNCATED SEARCH CANNOT MAKE THIS CLEAR REAL INPUT, and the instinct runs the other way.
        `mark_stereogenic` marks every undecided unit of a record whose search ran out of budget
        (decision 5), so a parity on a truncated candidate sits on a MARKED unit and survives; the
        conservative direction of that approximation keeps input rather than discarding it, and
        `stereo_truncated` stays a report and never a raise (ruling F62).

        THE CLEAR GOES INTO A CLONE (ruling F65).  `copy()` shares the arena outright on the grounds
        that it is immutable, so a clear written straight into SEG_PARITY would strip the parity
        from every container sharing it, including ones the caller never named.  So an empty report is
        a PURE READ -- no clone, no `_gen` bump, `shares_arena_with` still True -- and a non-empty one
        clones, clears in the clone and rebinds, exactly as `remap` does.

        The clone carries the derived stereo table verbatim (`structure_clone` copies the whole
        buffer), and its SU_STEREOGENIC marks were computed against the parities just removed:
        `mark_stereogenic` pins every CONFIGURED unit to its stored value, so dropping a parity
        loosens a pin and can only make a refuting witness easier to find.  The table is therefore
        invalidated in the clone (ruling F74) and the next reader derives it again.  That the table is
        dropped rather than carried is measured by the arena growing on the next read
        (`test_a_reported_parity_is_cleared`), not by the marks: no fixture I could build moves a mark
        across this clear, and the report says why.

        Idempotent, and only for the reason in the paragraph just above: the second call has nothing
        left to report BECAUSE clearing a parity removes a pin from `_stereo_consistent` and so can only
        loosen the system -- a unit that was marked stays marked -- and no record has been found
        where a clear moves a mark (the counts are with the assertion in
        `test_a_reported_parity_is_cleared`).  Should such a record turn up, this is a FIXPOINT
        rather than one-shot idempotence: a clear could then unmark a unit whose parity is still
        stated, the next call would report that one too, and a caller who needs the property must
        loop until the report comes back empty.
        """
        cdef uint32_t *slots = NULL
        cdef uint32_t count = 0
        cdef uint32_t k
        cdef uint32_t slot
        cdef list report = []
        cdef Structure fresh
        cdef uint8_t *fpar
        cdef object n
        self._require_clean()
        collect_stereo_rejections(self._structure, &slots, &count)
        try:
            for k in range(count):
                report.append(self._numbers[slots[k]])
        finally:
            PyMem_Free(slots)
        # SLOT order out of the collector, STABLE ID order out of here: the two agree until somebody
        # remaps, and a caller comparing this list against ids of its own is entitled to one order.
        report.sort()
        if not report:
            return report
        fresh = structure_clone(self._structure)
        # the marks in the copied table were decided against the parities about to go; drop the table
        structure_invalidate_stereo_units(fresh)
        # THE GUARDED FORM, not `structure_set_parity`: a molecule that states no parity has no
        # segment, so clearing a parity that is not stored is a no-op.
        fpar = structure_parities(fresh) if structure_has(fresh, SEG_PARITY) else NULL
        for n in report:
            slot = <uint32_t> self._index_of[n]
            if fpar is not NULL:
                fpar[slot] = 0
        # The clone carried SEG_FEATURES verbatim and word IV screens the SEG_PARITY byte just
        # cleared, so the derived words are re-based here (ruling F78).  Not optional bookkeeping:
        # `features_of()` and `_union_feature_words` hand those words to Python verbatim, so without
        # this the accessor disagrees with `parity_of` and with a `from_bytes` round trip of the same
        # molecule.  The features are re-derived rather than dropped because, unlike the stereo
        # table above, they are cheap to keep and every isomorphism call would pay for the drop.
        refresh_parity_features(fresh)
        self._structure = fresh
        self._gen += 1
        return report

    def clean_stereo(self):
        """Wipe EVERY kind of stereo state, unconditionally.  Returns what was wiped.

        THE OPPOSITE POSTURE FROM `validate_stereo`, and that is why both exist.  `validate_stereo`
        asks a question -- which stated parities this constitution can justify -- and clears only the
        answers it cannot.  This asks nothing and judges nothing: a caller reaching for it has already
        decided that what the molecule says about configuration is not to be kept.

        FOUR KINDS OF STATE, and clearing only the parities is the trap:

        * the atom parities -- ALL of them, justified or not;
        * the WEDGES on the edges.  Not optional: `chython/formats/ctfile` writes back the wedges a
          molecule carries rather than re-deriving them, and its reader derives parities FROM wedges,
          so a wipe that left them behind would be undone by one molfile round trip;
        * the ABS/AND/OR STEREO GROUP membership.  An AND-group membership with no parity inside it
          names a configuration that no longer exists;
        * the stored CIP DESCRIPTORS on atoms and bonds.  An `(R)` on an atom with no parity is
          actively wrong, and stored CIP exists for external consumers, so it would be trusted.

        COORDINATES ARE NOT DROPPED.  A layout is not a configuration, and a caller who wanted the
        drawing gone would have said so -- `clean2d()` is that call.  Dropping them would also destroy
        the only thing a depiction has to work with on a molecule the caller asked to flatten.

        THE WIPE GOES INTO A CLONE (ruling F65), for `validate_stereo`'s reason word for word:
        `copy()` shares the arena outright on the grounds that it is immutable, so a clear written
        straight into SEG_PARITY would strip the stereo out of every container sharing it,
        including ones the caller never named.  So an empty report is a PURE READ -- no clone, no
        `_gen` bump, `shares_arena_with` still True -- and a non-empty one clones, clears in the clone
        and rebinds.  The stereo unit table is invalidated in the clone (ruling F74) because
        `structure_clone` carries it verbatim and its SU_STEREOGENIC marks were computed against the
        parities just removed; feature word IV is re-based (ruling F78) because the clone carried
        SEG_FEATURES verbatim and that word screens the SEG_PARITY byte.

        THE REPORT IS KEYED BY READER, one key per kind, and a key is ABSENT when its reader was
        empty -- so `{}` means "this molecule had no stereo at all" and the return value is falsy
        exactly then.  `validate_stereo`'s flat list of stable ids is not extended because it cannot
        be: that method reports one kind of state, where a list of ids says everything there is to
        say, while this one touches five readers, and a union list would tell a caller that atom 2
        "had something" without saying what.  Each value is the corresponding reader's own answer,
        taken before the wipe and unchanged in shape:

            {'parities':      [n, ...],                    # ascending stable ids, as `validate_stereo`
             'wedges':        [(narrow, wide, wedge), ...],  # `wedges()`
             'stereo_groups': {(kind, group): [n, ...]},     # `stereo_groups()`
             'atom_cips':     {n: descriptor},               # `atom_cips()`
             'bond_cips':     {(n, m): descriptor}}          # `bond_cips()`

        NOTHING IS APPENDED TO `cip_log`.  That log exists for a descriptor lost as a SIDE EFFECT of
        an edit, where the event is unrecoverable from the bytes afterwards; here the drop is the
        thing the caller asked for and the descriptors are handed straight back.

        THE SEGMENT SURVIVES THE MEMBERSHIPS.  `has_stereo_groups` reports whether the arena carries
        SEG_STEREO_GROUPS, and this zeroes the payload rather than dropping the segment -- so it can
        still answer True on a molecule with no groups left.  That is not a new state: the public
        `set_stereo_group(n, 0)` leaves exactly it, and both CTfile writers read the memberships
        through `canonical_stereo_groups()`, which comes back empty either way.  Dropping the segment
        needs a whole new buffer (`structure_respan`), which is a cost this call has no reason to pay.
        """
        cdef Structure fresh
        cdef atom_t *atoms
        cdef uint32_t *ptr
        cdef halfedge_t *edges
        cdef halfedge_t *e
        cdef uint32_t i, k
        cdef uint32_t n_atoms
        cdef list parities
        cdef dict report
        cdef object key, value
        self._require_clean()
        n_atoms = self._structure.header.atom_count
        parities = []
        for i in range(n_atoms):
            if structure_parity_at(self._structure, i):
                parities.append(self._numbers[i])
        # SLOT order out of the walk, STABLE ID order out of here -- `validate_stereo`'s promise, and
        # a caller comparing the two lists is entitled to one order from both.
        parities.sort()
        report = {}
        if parities:
            report['parities'] = parities
        for key, value in (('wedges', self.wedges()), ('stereo_groups', self.stereo_groups()),
                           ('atom_cips', self.atom_cips()), ('bond_cips', self.bond_cips())):
            if value:
                report[key] = value
        if not report:
            return report
        fresh = structure_clone(self._structure)
        # the marks in the copied table were decided against the parities about to go; drop the table
        structure_invalidate_stereo_units(fresh)
        # EVERY POINTER BELOW IS TAKEN AFTER THE INVALIDATE (RULES section 2): it retires a segment,
        # and a pointer bound before it names the arena as it was.
        # The clone carried the segment; every stated parity in it goes.
        structure_clear_parities(fresh)
        atoms = fresh.atoms()
        for i in range(n_atoms):
            at_set_cip(&atoms[i], 0)
        ptr = csr_ptr(fresh)
        edges = csr_edges(fresh)
        for i in range(n_atoms):
            for k in range(ptr[i], ptr[i + 1]):
                e = &edges[k]
                # A wedge lives on ONE half-edge (the narrow end's) and a CIP descriptor on both;
                # walking every half-edge covers both cases without having to know which is which.
                e.wedge = 0
                he_set_cip(e, 0)
        if structure_has(fresh, SEG_STEREO_GROUPS):
            memset(structure_stereo_groups(fresh), 0,
                   structure_seg_len(fresh, SEG_STEREO_GROUPS))
        # word IV screens the SEG_PARITY byte, just cleared by structure_clear_parities above
        refresh_parity_features(fresh)
        self._structure = fresh
        self._gen += 1
        return report

    def translate_stereo(self, uint32_t n, tuple order):
        """Translate the stored parity of the unit anchored at `n` to the caller's direction order.

        `order` is a tuple of stable ids (with `None` for a direction that has no atom of its
        own, or for a pinned slot that is no direction at all) representing the caller's desired
        ordering of the unit's directions.  Returns the parity in that order: 0 when the unit's
        parity is unset, 1 when even, 2 when odd.

        For bond kinds (cis/trans, allene, atropisomer) the two pairs in `refs` must be
        preserved.  The ``unnamed_mask`` (``u.spare >> SU_UNNAMED_SHIFT``) distinguishes real
        unnamed directions (mask bit set -- e.g. an implicit hydrogen, movable within its pair)
        from pinned slots (mask bit clear with refs[i] == SU_NO_REF -- e.g. a nitrogen lone pair,
        not a direction at all).  ``order[0:2]`` must map entirely to one stored pair and
        ``order[2:4]`` to the other; a cross-pair mix is rejected with ``ValueError``.

        A pinned slot is frozen at **its offset within its pair**, not at its absolute index
        (Ruling F55).  The wholesale pair exchange legally moves a pinned slot at absolute index 1
        to absolute index 3 (carrying its pair along), and that is accepted.  An absolute-index
        reading (F51, superseded) would make every reversed oxime raise -- which
        ``test_acetaldoxime_on_pair_exchange_does_not_flip`` explicitly refutes.  What is
        forbidden is swapping a pinned slot *within* its pair; that raises ``ValueError``.

        Raises `KeyError` when `n` is not in the molecule or does not anchor any unit.
        Raises `ValueError` when `order` is not a valid permutation of the unit's refs.
        """
        # All cdef declarations at the top of the function scope (Cython rule: no cdef inside
        # conditional blocks).
        cdef uint32_t anchor_slot
        cdef stereo_unit_t *u
        cdef uint32_t n_refs
        cdef uint32_t want[4]
        cdef uint32_t i, j, k
        cdef uint32_t refs_used        # bitmask: bit j = refs[j] already matched (phases 1 and 3)
        cdef uint32_t n_no_ref_want    # SU_NO_REF count in want (phase 1)
        cdef uint32_t n_no_ref_refs    # SU_NO_REF count in refs (phase 1)
        cdef bint found
        cdef uint32_t perm[4]          # perm[i] = j means want[i] comes from refs[j] (atom kind)
        cdef uint32_t refs_no_ref[4]   # positions of unnamed slots in refs (atom kind)
        cdef uint32_t n_refs_no_ref    # count of above
        cdef uint32_t nr_used          # how many refs_no_ref slots consumed so far
        cdef uint32_t unnamed_mask     # (u.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK
        cdef uint32_t src_pair         # which stored pair want[0:2] maps to (bond kind)
        cdef uint32_t base             # first index of the k-th pair in the per-pair loop
        cdef uint32_t wo               # want offset for the k-th pair (0 or 2)
        cdef uint32_t swap[2]          # within-pair swap flags: swap[k]=0 or 1
        cdef uint32_t pp               # XOR of swap[0] and swap[1] (bond kind parity contribution)
        cdef uint32_t named            # named atom of the current pair (slot 0 by F26)
        cdef uint32_t other            # other-slot of the current pair (slot 1)
        cdef uint8_t parity
        cdef object elem

        # Initialise perm to zeros -- seatbelt: an unmatched slot (unreachable if the
        # validations below all pass) reads PERM_PARITY_4[0] = 0 rather than from garbage.
        perm[0] = 0; perm[1] = 0; perm[2] = 0; perm[3] = 0
        # Initialise swap and pp; they are always set before use on the bond path, but Cython
        # cannot prove that from the two-check structure.
        swap[0] = 0u; swap[1] = 0u
        pp = 0u

        self._require_clean()
        if n not in self._index_of:
            raise KeyError(n)
        # Builds the table for the refs and the unnamed mask, and NEVER READS SU_STEREOGENIC: this is
        # arithmetic on a stored parity in the caller's direction order, so a conservatively-marked
        # table (`stereo_truncated`) changes nothing here -- an unproven mark is not one of its inputs.
        ensure_stereo_units(self._structure)
        anchor_slot = <uint32_t> self._index_of[n]
        u = stereo_unit_of(self._structure, anchor_slot)
        if u is NULL:
            raise KeyError(n)

        n_refs = u.n_refs
        unnamed_mask = (u.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK
        if len(order) != n_refs:
            raise ValueError('order must have exactly %d elements' % n_refs)

        # Convert stable ids to slot values; None -> SU_NO_REF
        for i in range(n_refs):
            elem = order[i]
            if elem is None:
                want[i] = SU_NO_REF
            else:
                if elem not in self._index_of:
                    raise ValueError('stable id %d is not in this molecule' % elem)
                want[i] = <uint32_t> self._index_of[elem]

        # ================================================================
        # Phase 1 — kind-independent validation (above the kind split).
        # These two checks run for BOTH bond and atom kinds.  No validation
        # lives below the kind split except checks that are genuinely about
        # that kind's pair structure; any permutation-level invariant belongs
        # here so that no future restructure of the computation can drop it.
        # ================================================================

        # Named-atom uniqueness: every named want[i] must appear in refs exactly once.
        refs_used = 0
        for i in range(n_refs):
            if want[i] != SU_NO_REF:
                found = False
                for j in range(n_refs):
                    if u.refs[j] == want[i] and not (refs_used & (1u << j)):
                        refs_used |= 1u << j
                        found = True
                        break
                if not found:
                    raise ValueError(
                        'order is not a permutation of the unit refs: '
                        'atom not found or duplicated in order')

        # None count: want's SU_NO_REF count must equal refs' SU_NO_REF count.
        # For bond kinds a pinned slot appears as None in want (just as an unnamed direction does)
        # so the count is still one per SU_NO_REF slot regardless of kind.
        n_no_ref_want = 0
        n_no_ref_refs = 0
        for j in range(n_refs):
            if u.refs[j] == SU_NO_REF:
                n_no_ref_refs += 1
        for i in range(n_refs):
            if want[i] == SU_NO_REF:
                n_no_ref_want += 1
        if n_no_ref_want != n_no_ref_refs:
            raise ValueError(
                'order is not a permutation of the unit refs: '
                'None count does not match the number of unnamed directions')

        # ================================================================
        # Phase 2 — kind-specific structural validation.
        # Bond: pair assignment, cross-pair membership, None correspondence,
        #        and pin check -- one loop over the two pairs.
        # Atom: phase 1 is sufficient; no further structural checks.
        # ================================================================
        if u.kind != SU_TETRA:
            # ----------------------------------------------------------------
            # BOND KIND: direct pair decomposition (Rulings F55, F56).
            # ----------------------------------------------------------------
            # refs[0:2] = P0 (anchor terminal's pair), refs[2:4] = P1 (other terminal).
            # By Ruling F26, each pair's named atom is always at offset 0 within the pair;
            # offset 1 is the "other-slot": a named atom, a real unnamed direction (SU_NO_REF
            # with unnamed_mask bit set), or a pinned slot (SU_NO_REF + mask bit clear).
            #
            # The per-pair checks and swap computation are one loop over k in range(2).
            # swap[k] = 0 means the k-th assigned pair's slots appear in refs order;
            # swap[k] = 1 means they are reversed.  A pinned other-slot forbids swap.
            #
            # PERM_PARITY_4 / _perm_index4 / permutation_parity_of serve SU_TETRA only; not
            # used here.

            # -- pair assignment: find which stored pair's named atom (at slot 0 by F26)
            #    appears in want[0:2].  The search includes both want[0] and want[1] so that
            #    a want where the named atom is in the second slot is still resolved correctly
            #    (the pair's own slot-0 atom is at want[1]).
            src_pair = 0u
            found = False
            for i in range(2):
                if want[i] != SU_NO_REF:
                    for j in range(2):
                        if u.refs[j * 2u] == want[i]:
                            src_pair = <uint32_t> j
                            found = True
                            break
                    if found:
                        break
            if not found:
                # want[0:2] contains no named atom from any pair.  By F47 every terminal has
                # at least one named substituent, so a legal want[0:2] always contains one.
                #
                # REACHABLE, but never the SOLE violation: measured over 2,229 orders on nine
                # bond fixtures, every order that reaches this raise is also rejected by the
                # k=0 cross-pair check or by None correspondence, both of which run later.  So
                # disabling this raise alone fails no test -- that is a fact about check
                # ORDERING, not about reachability, and it is not a licence to delete the
                # branch: without it an illegal order would be reported against the wrong pair.
                raise ValueError(
                    'want[0:2] contains no named atom from either stored pair; the order is '
                    'illegal (every pair has exactly one named atom per Ruling F26 / F47)')

            # Per-pair loop: k=0 handles want[0:2] (pair `src_pair`),
            #                k=1 handles want[2:4] (the other pair).
            for k in range(2):
                if k == 0:
                    base = src_pair * 2u       # first index of the pair for want[0:2]
                    wo = 0u                    # want offset
                else:
                    base = (1u - src_pair) * 2u
                    wo = 2u

                # Cross-pair: every named atom in want[wo:wo+2] must belong to refs[base:base+2].
                for i in range(2):
                    if want[wo + i] != SU_NO_REF:
                        if u.refs[base] != want[wo + i] and u.refs[base + 1u] != want[wo + i]:
                            raise ValueError(
                                'order mixes directions across pair boundary; bond-kind units '
                                'have two ordered pairs and order[%d:%d] must come from the '
                                'same stored pair (base %d)' % (wo, wo + 2u, base))

                # Swap: slot 0 of the pair is always the named atom (F26).
                # If want[wo] != refs[base], the pair is presented in reversed order (swap=1).
                named = u.refs[base]
                if named == SU_NO_REF:
                    # Unreachable: perception (F47) refuses any terminal with no named
                    # substituent, so refs[base] is always named.  Guard is a future-proof raise.
                    raise ValueError(
                        'pair at base %d has SU_NO_REF at slot 0; perception invariant '
                        'F26/F47 violated -- this is a bug in the perception layer' % base)
                swap[k] = 0u if want[wo] == named else 1u

                # None correspondence: each None in want[wo:wo+2] must map to a SU_NO_REF slot.
                #
                # REACHABLE whenever the two pairs hold DIFFERENT numbers of SU_NO_REF slots --
                # the commonest real E/Z shape, one terminal disubstituted and the other bearing
                # a hydrogen.  (Z)-2-chlorobut-2-ene, refs = (Cl, CH3, CH3', None), reaches it
                # on order (Cl, None, CH3', CH3): four such orders per spelling.  When the two
                # pairs hold EQUAL SU_NO_REF counts it is pre-empted -- with none, by phase 1's
                # None count; with one each, by uniqueness / None count / src_pair / the k=0
                # cross-pair check.
                #
                # It is never the SOLE violation, though: a misplaced None implies a misplaced
                # named atom, which the OTHER pair's cross-pair check also rejects -- one
                # iteration later.  So disabling this branch alone changes the message and not
                # the outcome, which is why only a message-discriminating test can pin it
                # (test_chlorobut_2_ene_none_does_not_correspond_raises).  That is a fact about
                # check ordering, not about reachability.
                for i in range(2):
                    if want[wo + i] == SU_NO_REF:
                        if u.refs[base + (i ^ swap[k])] != SU_NO_REF:
                            raise ValueError(
                                'None in order at position %d does not correspond to a '
                                'SU_NO_REF slot in pair at base %d' % (wo + i, base))

                # Pin check: if the other-slot is pinned (SU_NO_REF, mask bit clear),
                # the within-pair order is frozen -- swap is forbidden (Ruling F55).
                other = u.refs[base + 1u]
                if other == SU_NO_REF and not (unnamed_mask & (1u << (base + 1u))):
                    if swap[k]:
                        raise ValueError(
                            'order swaps a pinned slot in pair at base %d; '
                            'the within-pair order is frozen (Ruling F55)' % base)

            # Ruling F56 -- why the answer is `stored XOR swap[0] XOR swap[1]` and nothing else.
            # Any legal bond-kind order is a composition of at most three permutations of the
            # four slots: the wholesale pair exchange, and one within-pair transposition per
            # pair.  Parity is a homomorphism to {even, odd}, so their contributions XOR:
            #
            #   * the pair exchange is (0 2)(1 3) -- TWO transpositions, an EVEN permutation --
            #     so it contributes NOTHING.  That is why `src_pair` only selects which pair
            #     want[0:2] is compared against and never enters the arithmetic, and why
            #     test_but_2_ene_pair_exchange_does_not_flip expects an unflipped parity.
            #   * each within-pair swap is ONE transposition, odd, contributing 1.
            #
            # Hence pp = swap[0] ^ swap[1], and the SU_TETRA permutation table is not needed
            # (and would be wrong here: it cannot express the pair constraint).
            pp = swap[0] ^ swap[1]

        # ================================================================
        # Phase 3 — computation.
        # ================================================================

        # Parity is read from SEG_PARITY at the anchor's slot.
        parity = structure_parity_at(self._structure, anchor_slot)

        if u.kind != SU_TETRA:
            # Bond kind: XOR the two within-pair swaps; no perm table needed.
            if parity == 0u:
                return 0
            return <int> (((parity - 1u) ^ pp) + 1u)

        # ATOM KIND (SU_TETRA): build perm from phase-1-validated want, then use perm table.
        # All SU_NO_REF in TETRA are real unnamed directions (mask bit set); no pinned slots.
        n_refs_no_ref = 0
        for j in range(n_refs):
            if u.refs[j] == SU_NO_REF:
                refs_no_ref[n_refs_no_ref] = j
                n_refs_no_ref += 1
        nr_used = 0
        refs_used = 0
        for i in range(n_refs):
            if want[i] == SU_NO_REF:
                perm[i] = refs_no_ref[nr_used]
                nr_used += 1
            else:
                for j in range(n_refs):
                    if u.refs[j] == want[i] and not (refs_used & (1u << j)):
                        perm[i] = j
                        refs_used |= 1u << j
                        break
        return translate_parity(parity, perm)

    # --------------------------------------------------------------------------------------------
    # S-GROUPS, ALIASES AND THE TITLE.  Three names for one storage decision: everything a CTfile
    # attaches to a SET OF ATOMS rather than to an atom.
    #
    # REFERENCES ARE STABLE IDS ON THIS SIDE AND INDICES ON THE OTHER, and the asymmetry is the point
    # of having a boundary here at all.  The arena stores atom INDICES because that makes the carry
    # across an edit one loop over one array against one map, with no dict; a caller cannot use indices
    # because they are re-assigned by every deletion.  So the translation happens exactly here, in both
    # directions, and no layer above the core ever sees an index.
    #
    # THE WHOLE SET IS REPLACED, NEVER AMENDED, and that is a storage fact and not a taste: the three
    # S-group segments are persistent, their sizes are a function of the data, and a persistent segment
    # never grows in place (ruling F60).  `structure_respan` says the rest.
    # --------------------------------------------------------------------------------------------

    cdef inline bytes _raw_title(self):
        """The name line as the blob stores it.  `set_sgroups` and `set_aliases` re-serialise the blob
        and must not encode a title they never decoded."""
        return blob_bytes(self._structure, SEG_OPAQUE_BLOB, 0)

    @property
    def title(self):
        """The molecule's name line, as `str`.  `''` for a molecule never given one -- absent and empty
        are the same answer for a title, unlike for a coordinate.

        Handle 0 of the blob, which is why a molecule with a title and no S-groups still carries the
        segment.  THE BLOB STORES BYTES AND THIS DECODES WITH `surrogateescape`: an SDF name line is not
        required to be UTF-8, and an undecodable byte becomes a lone surrogate that `set_title` and this
        library's writers re-encode to the same byte.  The fidelity promise is therefore kept, and kept
        as something a test can state rather than as a type every caller has to decode.

        The cost, stated because it is real: such a title raises `UnicodeEncodeError` if it is encoded as
        UTF-8 WITHOUT the handler -- `json.dumps` on it, or a stream opened with no `errors=`.  That
        needs an undecodable byte in the source file, and it is the trade `os.fsdecode` makes for the
        same reason.  XML cannot take that trade at all, which is why the CML and MRV writers replace
        instead (`chython/formats/xml/_dialect.py::xml_text`).
        """
        self._require_clean()
        return self._raw_title().decode('utf8', 'surrogateescape')

    def set_title(self, title):
        """Replace the name line.  `str`, `bytes`, `bytearray` or `memoryview`; a raw name line is
        honoured, and comes back as the `str` `surrogateescape` makes of it."""
        self._require_clean()
        self._rebuild_blob(_as_title_bytes(title), self._sgroup_records())
        return None

    @property
    def meta(self):
        """Record metadata -- an SDF's data fields, an RDfile's DTYPE/DATUM pairs.  Created on first
        access.

        A PLAIN DICT, and the same property `ReactionContainer.meta` is.  One implementation for both
        containers is the whole point: the reason `CtfileRecord` and `FieldsView` existed was that this
        was missing, and a second mapping class here would have re-created them in the core.

        NOT in the arena and not in `to_bytes`: the arena is chemistry plus what a file drew, and a
        boiling point is neither.  `__reduce__` carries it, `copy()` copies it shallow, `substructure`
        and `split` start empty -- a part is not the record the metadata described.
        """
        if self._meta is None:
            self._meta = {}
        return self._meta

    cdef inline int _log_event(self, str rule, str stage, str message, str severity) except -1:
        """The extension's only door into `log`.  `stage` is what the filtered views select on."""
        mc_lazy_log_imports()
        if self._log is None:
            self._log = _MC_LOG()
        self._log.append(_MC_RECORD(rule, (), message, severity, stage, ''))
        return 0

    @property
    def log(self):
        """What was read, repaired or lost about THIS molecule, oldest first.  Created on first use.

        THE STORAGE, AND IT IS NEVER CONDITIONAL.  Every reader, every pass and every edit session
        writes here, whether or not anyone asked: `mol.canonicalize()` with no arguments fills this, and
        there is deliberately no flag, no verbosity level and no `log is not None` branch anywhere on the
        path that could turn recording off.  NO PASS TAKES A `log=` -- reading it is `mol.log`, and a
        second sequence to write to is how a repair with nobody watching stops being recorded.  A reader
        still takes one, having no container to write to until it has produced one.

        It accumulates: nothing here clears it, so a caller wanting one call's records alone brackets
        the call with `len(mol.log)` or `del mol.log[:]`.

        A `chython.core.Log`, so `by_stage`, `by_severity`, `repaired()` and `lost()` all work.  ONE
        STORAGE: `sgroup_log` and `cip_log` are read-only views over this, filtered by stage.

        PER HANDLE.  `copy()` does not carry it, because `cip_log` records what THIS handle's edits
        lost -- the standing `core/test/test_cip_storage.py:306` states, now that both live here.  It is
        in neither `to_bytes` nor `pack`: derived diagnostics, not data.
        """
        if self._log is None:
            mc_lazy_log_imports()
            self._log = _MC_LOG()
        return self._log

    @property
    def sgroup_log(self):
        """Diagnostics this molecule accumulated about its own S-groups, oldest first.

        WHY A CONTAINER-LEVEL LIST AND NOT THE RECORDS' OWN `log` RUN.  A record's log is in the blob,
        the blob is copied byte for byte by every carry, and a carry is exactly when a reference is
        lost -- so the one event that most needs recording is the one event that cannot be written
        where the others live.  A Python list can grow; a persistent segment cannot.

        A READ-ONLY VIEW over `log`, filtered to the `edit:sgroup` stage; the storage is one list.
        """
        self._require_clean()
        if self._log is None:
            return ()
        # an explicit loop and not a genexpr: a comprehension's loop variable is a name Cython never saw
        # declared, and `warn.undeclared` reports it
        cdef object rec
        cdef list out = []
        for rec in self._log.by_stage('edit:sgroup'):
            out.append(str(rec))
        return tuple(out)

    @property
    def sgroups(self):
        """Every S-group record, in the order the file gave them, as dicts.

        Atom aliases are NOT here even though they share the storage -- see :attr:`aliases`.  A record
        whose atoms have all been deleted IS here, empty: an emptied record still says a field was
        attached to something, and dropping it is the silent loss the whole segment exists to prevent.
        """
        self._require_clean()
        cdef sgroup_t *rec = structure_sgroups(self._structure)
        cdef uint32_t n = structure_sgroup_count(self._structure)
        cdef uint32_t r
        cdef list out = []
        for r in range(n):
            if rec[r].flags & SGROUP_FLAG_ALIAS:
                continue
            out.append(self._sgroup_dict(&rec[r]))
        return tuple(out)

    @property
    def aliases(self):
        """Atom aliases as {n: bytes} -- V2000 `A  <n>` lines, MRV mrvAlias.

        STORED AS S-GROUP RECORDS, WHICH IS THE PROPOSAL THIS ACCESSOR IS.  An alias is not an S-group
        by any reading of the spec, but its STORAGE REQUIREMENT is an S-group's exactly and in full: a
        label bound to a set of atoms (of size one), which must survive a remap, must follow its atom
        into a substructure, and must be dropped AND REPORTED when its atom dies.  Every one of those is
        already implemented once for records; a second mechanism would be a second thing to get wrong,
        and the arena would need a fourth segment to hold it.  So an alias is a record with
        SGROUP_FLAG_ALIAS, exactly one atom, its text in the name slot and an empty `type`.
        """
        self._require_clean()
        cdef sgroup_t *rec = structure_sgroups(self._structure)
        cdef uint32_t *idx = structure_sgroup_index(self._structure)
        cdef uint32_t n = structure_sgroup_count(self._structure)
        cdef uint32_t r
        cdef dict out = {}
        for r in range(n):
            if not rec[r].flags & SGROUP_FLAG_ALIAS or rec[r].atoms_len != 1:
                # atoms_len 0 is an alias whose atom was deleted. It stays in storage (invariant 1)
                # and disappears from this view, because a label with no atom has nothing to label.
                continue
            out[self._numbers[idx[rec[r].refs_off]]] = blob_bytes(
                self._structure, SEG_OPAQUE_BLOB, rec[r].strings_off + 2)
        return out

    def set_sgroups(self, records):
        """Replace every S-group record, leaving aliases and the title alone.

        `records` is an iterable of dicts shaped like the ones :attr:`sgroups` returns; every key is
        optional except that an unknown key is an ERROR rather than ignored, because a misspelled
        `patoms` that silently did nothing is a lost reference set with no diagnostic.
        """
        self._require_clean()
        cdef list out = []
        cdef object r
        for r in records:
            out.append(_sgroup_normalise(r, False))
        for r in self._sgroup_records():
            if r['_alias']:
                out.append(r)
        self._rebuild_blob(self._raw_title(), out)
        return None

    def set_aliases(self, mapping):
        """Replace every atom alias, leaving S-group records and the title alone."""
        self._require_clean()
        cdef list out = []
        cdef object r, n, text
        for r in self._sgroup_records():
            if not r['_alias']:
                out.append(r)
        for n, text in dict(mapping).items():
            self._require(<uint32_t> n)
            out.append(_sgroup_normalise({'atoms': (n,), 'name': _as_bytes(text, 'alias')}, True))
        self._rebuild_blob(self._raw_title(), out)
        return None

    def add_data_sgroup(self, name, data, *, atoms=(), bonds=(), disp=None, log=None):
        """Attach one CTfile `DAT` S-group -- a text label on atoms or bonds -- and return the record.

        APPENDS, unlike `set_sgroups` above, which replaces the whole set: two calls give two labels.
        `data` is a string or a list of them for a multi-value `FIELDDATA`.  `atoms` and `bonds` are the
        references, a bond being an `(n, m)` pair, and an unknown atom or an unbonded pair raises.

        `disp=` is the `FIELDDISP` anchor: `None` computes the mean of the referenced atoms'
        coordinates, `(x, y)` states one, and `False` writes none.  A molecule with no coordinates has
        none to give, so the record is written anchorless and the reason is logged.

        Registered by `chython.formats`, not implemented here.  The core owns S-group STORAGE; what a
        `DAT` record MEANS is CTfile knowledge -- see `_set_sgroup_fns`.
        """
        return _sgroup_fn('add_data_sgroup')(self, name, data, atoms=atoms, bonds=bonds, disp=disp,
                                            log=log)

    def data_sgroups(self, name=None):
        """Every `DAT` S-group on this molecule as a list, or only those with this `FIELDNAME`.

        The read beside `add_data_sgroup` above, and the same layer split: `sgroups` hands back raw
        dicts from storage, this hands back the parsed records.  Registered by `chython.formats`.
        """
        return _sgroup_fn('data_sgroups')(self, name)

    cdef dict _sgroup_dict(self, sgroup_t *rec):
        """One record as a dict, with every atom index turned back into a stable id."""
        cdef uint32_t *idx = structure_sgroup_index(self._structure)
        cdef list numbers = self._numbers
        cdef Structure s = self._structure
        cdef uint32_t base = rec.refs_off
        cdef uint32_t sbase = rec.strings_off
        cdef uint32_t k, a, b
        cdef list atoms = [], patoms = [], bonds = [], cstates = [], data = [], fields = [], log = []
        for k in range(rec.atoms_len):
            atoms.append(numbers[idx[base + k]])
        base += rec.atoms_len
        for k in range(rec.patoms_len):
            patoms.append(numbers[idx[base + k]])
        base += rec.patoms_len
        for k in range(0, rec.bonds_len, 2):
            bonds.append((numbers[idx[base + k]], numbers[idx[base + k + 1]]))
        base += rec.bonds_len
        # The CSTATE tails are the run AFTER data, fields and log; `sgroup_strings_len` fixes that
        # order and this arithmetic is the only place that reads it.
        cdef uint32_t tails = sbase + 4 + rec.data_len + rec.fields_len + rec.log_len
        for k in range(0, rec.cstates_len, 2):
            a = idx[base + k]
            b = idx[base + k + 1]
            cstates.append(((None if a == SGROUP_NO_REF else (numbers[a], numbers[b])),
                            blob_bytes(s, SEG_OPAQUE_BLOB, tails + (k >> 1))))
        for k in range(rec.data_len):
            data.append(blob_bytes(s, SEG_OPAQUE_BLOB, sbase + 4 + k))
        for k in range(0, rec.fields_len, 2):
            fields.append((blob_bytes(s, SEG_OPAQUE_BLOB, sbase + 4 + rec.data_len + k),
                           blob_bytes(s, SEG_OPAQUE_BLOB, sbase + 4 + rec.data_len + k + 1)))
        for k in range(rec.log_len):
            log.append(blob_bytes(s, SEG_OPAQUE_BLOB,
                                  sbase + 4 + rec.data_len + rec.fields_len + k))
        return {'type': blob_bytes(s, SEG_OPAQUE_BLOB, sbase),
                'subtype': blob_bytes(s, SEG_OPAQUE_BLOB, sbase + 1),
                'name': blob_bytes(s, SEG_OPAQUE_BLOB, sbase + 2),
                # disp is a member of packed sgroup_t (line 391), so &rec.disp is forbidden by RULES.md
                # §2.3; this divides in place rather than calling xy_read_x/xy_read_y
                'disp': ((<double> rec.disp.x / XY_SCALE, <double> rec.disp.y / XY_SCALE)
                         if rec.flags & SGROUP_FLAG_DISP else None),
                'disp_tail': blob_bytes(s, SEG_OPAQUE_BLOB, sbase + 3),
                'index': rec.index, 'ext_index': rec.ext_index, 'parent': rec.parent,
                'atoms': tuple(atoms), 'patoms': tuple(patoms), 'bonds': tuple(bonds),
                'cstates': tuple(cstates), 'data': tuple(data), 'fields': tuple(fields),
                'log': tuple(log), '_alias': (rec.flags & SGROUP_FLAG_ALIAS) != 0}

    cdef list _sgroup_records(self):
        """Every record including aliases, as dicts -- the input shape of `_rebuild_blob`."""
        cdef sgroup_t *rec = structure_sgroups(self._structure)
        cdef uint32_t n = structure_sgroup_count(self._structure)
        cdef uint32_t r
        cdef list out = []
        for r in range(n):
            out.append(self._sgroup_dict(&rec[r]))
        return out

    cdef int _rebuild_blob(self, bytes title, list records) except -1:
        """Serialise `title` and `records` into a NEW arena, replacing this molecule's buffer.

        The order of the two runs is the layout, stated once:
        index run   -- atoms | patoms | bonds | cstates, each atom index one slot;
        blob run    -- type, subtype, name, disp_tail, data*, fields*, log*, cstate tails.
        """
        cdef dict index_of = self._index_of
        cdef list items = [title]
        cdef list slots = []
        cdef list recs = []
        cdef object d, a, b, pair, tail, key, value
        cdef dict rec
        cdef dict numbered = {}
        for d in records:
            if d['index'] != SGROUP_NO_INDEX:
                numbered[d['index']] = True
        for d in records:
            if d['parent'] != SGROUP_NO_INDEX and d['parent'] not in numbered:
                raise ValueError('sgroup names parent %d, which no record in this set carries as '
                                 'its index' % d['parent'])
        for d in records:
            rec = {'refs_off': len(slots), 'strings_off': len(items)}
            for a in d['atoms']:
                slots.append(index_of[a])
            for a in d['patoms']:
                slots.append(index_of[a])
            for a, b in d['bonds']:
                slots.append(index_of[a])
                slots.append(index_of[b])
            for pair, tail in d['cstates']:
                if pair is None:
                    slots.append(SGROUP_NO_REF)
                    slots.append(SGROUP_NO_REF)
                else:
                    slots.append(index_of[pair[0]])
                    slots.append(index_of[pair[1]])
            items.append(d['type'])
            items.append(d['subtype'])
            items.append(d['name'])
            items.append(d['disp_tail'])
            items.extend(d['data'])
            for key, value in d['fields']:
                items.append(key)
                items.append(value)
            items.extend(d['log'])
            for pair, tail in d['cstates']:
                items.append(tail)
            recs.append((rec, d))

        cdef uint32_t var_len[3]
        var_len[0] = <uint32_t> (len(recs) * sizeof(sgroup_t))
        var_len[1] = <uint32_t> (len(slots) * sizeof(uint32_t))
        var_len[2] = <uint32_t> blob_size_for(items)
        # A BLOB IS ALLOCATED EVEN FOR AN EMPTY TITLE AND NO RECORDS, because `items` always holds the
        # title handle -- `blob_size_for([b''])` is 24, never 0 -- and a zero-length segment would read
        # as absent.  The cost is 24 bytes on a molecule that has neither, which is why parsers only
        # call this when there is something to store.
        cdef Structure fresh = structure_respan(self._structure, var_len)
        cdef sgroup_t *out = structure_sgroups(fresh)
        cdef uint32_t *oidx = structure_sgroup_index(fresh)
        cdef uint32_t i
        for i in range(<uint32_t> len(slots)):
            oidx[i] = <uint32_t> slots[i]
        for i in range(<uint32_t> len(recs)):
            rec, d = recs[i]
            memset(&out[i], 0, sizeof(sgroup_t))
            out[i].refs_off = <uint32_t> rec['refs_off']
            out[i].strings_off = <uint32_t> rec['strings_off']
            out[i].atoms_len = <uint32_t> len(d['atoms'])
            out[i].patoms_len = <uint32_t> len(d['patoms'])
            out[i].bonds_len = <uint32_t> (2 * len(d['bonds']))
            out[i].cstates_len = <uint32_t> (2 * len(d['cstates']))
            out[i].index = <uint16_t> d['index']
            out[i].ext_index = <uint16_t> d['ext_index']
            out[i].parent = <uint16_t> d['parent']
            out[i].data_len = <uint16_t> len(d['data'])
            out[i].fields_len = <uint16_t> (2 * len(d['fields']))
            out[i].log_len = <uint16_t> len(d['log'])
            if d['_alias']:
                out[i].flags |= SGROUP_FLAG_ALIAS
            if d['disp'] is not None:
                out[i].flags |= SGROUP_FLAG_DISP
                out[i].disp.x = _fixed_point(d['disp'][0])
                out[i].disp.y = _fixed_point(d['disp'][1])
        structure_put_blob(fresh, SEG_OPAQUE_BLOB, items)
        rebuild_derived(fresh)
        self._structure = fresh
        self._gen += 1
        # Every cache keyed on `_gen` invalidates itself, but `_order_cache` holds a dict rather than a
        # generation-stamped scalar, so it is dropped by hand exactly as `_apply` does.
        self._order_cache = None
        return 0

    @property
    def connected_components_count(self):
        # isolated components: salts as ion pairs, and anything a reaction record glued together
        self._require_clean()
        cdef uint32_t n_refs = self._structure.header.atom_count
        if n_refs == 0:
            return 0
        cdef uint32_t *label = <uint32_t *> PyMem_Malloc(<size_t> n_refs * sizeof(uint32_t))
        if label is NULL:
            raise MemoryError('component labelling allocation failed')
        cdef Py_ssize_t comps
        try:
            with nogil:
                comps = label_components(self._structure, label)
        finally:
            PyMem_Free(label)
        if comps < 0:
            raise MemoryError('component labelling scratch allocation failed')
        return comps

    @property
    def connected_components(self):
        # one tuple of stable ids per component, each ascending by arena index
        self._require_clean()
        cdef list out = []
        cdef uint32_t n_atoms = self._structure.header.atom_count
        if n_atoms == 0:
            return out
        cdef uint32_t *label = <uint32_t *> PyMem_Malloc(<size_t> n_atoms * sizeof(uint32_t))
        if label is NULL:
            raise MemoryError('component labelling allocation failed')
        cdef Py_ssize_t comps, c
        cdef uint32_t i
        try:
            with nogil:
                comps = label_components(self._structure, label)
            if comps < 0:
                raise MemoryError('component labelling scratch allocation failed')
            for c in range(comps):
                out.append([])
            for i in range(n_atoms):
                (<list> out[label[i]]).append(self._numbers[i])
        finally:
            PyMem_Free(label)
        for c in range(comps):
            out[c] = tuple(out[c])
        return out

    @property
    def rings_count(self):
        # the circuit rank of the non-dative subgraph, since `rings` is its minimum cycle basis
        self._require_clean()
        if not structure_has(self._structure, SEG_RELEVANT_RINGS):
            return 0
        return structure_rings(self._structure)[0]

    @property
    def sssr(self):
        # the smallest set of smallest rings is exactly a minimum cycle basis
        return self.rings

    @property
    def rings(self):
        # A minimum cycle basis: the shortest independent cycles, one per unit of circuit rank.
        # Not the relevant-cycle set -- that is exponential, and the per-atom descriptors
        # (ring_sizes_of, ring_count_of, shares_ring) are the ones that carry its full
        # semantics, derived from prototypes without ever materialising the cycles.
        #
        # Order-8 bonds are excluded, so ferrocene is two five-rings and its iron is in none of
        # them. `mark_bridges` is where that happens; see its docstring for why it has to.
        self._require_clean()
        if not structure_has(self._structure, SEG_RELEVANT_RINGS):
            return []
        cdef uint32_t *r = structure_rings(self._structure)
        cdef uint32_t count = r[0]
        cdef uint32_t base = 2 + count
        cdef uint32_t i, k
        cdef uint32_t n_atoms = <uint32_t> len(self._numbers)
        cdef list out = []
        cdef list row
        for i in range(count):
            row = []
            for k in range(r[1 + i], r[2 + i]):
                # Defence in depth (Ruling F60), not a reachable error: with `_numbers`
                # correct, every ring member is a valid arena index, and nothing in the suite
                # reaches this raise.  It exists because this module compiles with
                # boundscheck=False, so a C uint32_t indexing a Python list reads a stale
                # PyObject* instead of raising IndexError -- which is precisely how a corrupt
                # `_numbers` presented as garbage integers inside `rings` rather than as an
                # exception.  A corrupt stable id must be loud.
                if r[base + k] >= n_atoms:
                    raise AssertionError(
                        'ring member index %d is out of range for %d atoms; the stable id '
                        'table is corrupt' % (r[base + k], n_atoms))
                row.append(self._numbers[r[base + k]])
            out.append(tuple(row))
        return out

    def ring_count_of(self, uint32_t n):
        return at_ring_count(self._atom(n))

    def ring_sizes_word_of(self, uint32_t n):
        return self._atom(n).ring_sizes

    def ring_sizes_of(self, uint32_t n):
        cdef uint32_t w = self._atom(n).ring_sizes
        cdef uint32_t size
        cdef list sizes = []
        for size in range(3, 25):
            if w >> size & 1:
                sizes.append(size)
        return frozenset(sizes)

    def macrocycle_of(self, uint32_t n):
        # ring_sizes is one uint32: bits 3-24 are exact sizes, and the three low bits are all
        # that is left for anything bigger. They exist so an atom on a 30-membered macrolactone
        # is not indistinguishable from an acyclic one -- the exact size is not recoverable from
        # them, so this reports the fact and `rings` carries the size.
        return (self._atom(n).ring_sizes & 7) != 0

    def shares_ring(self, uint32_t n, uint32_t m):
        self._require_clean()
        cdef uint32_t ia = self._index_of[n]
        cdef uint32_t ib = self._index_of[m]
        return structure_shares_ring(self._structure, ia, ib)

    @property
    def atoms_order(self):
        """Symmetry classes as {n: rank}, 1-based. Equivalent atoms share a rank.

        Refined to a fixed point, so two atoms share a rank only when no walk out of either one
        can tell them apart -- benzene is one class, toluene's ring is four. The ranks are ordered by
        the invariant itself (element first, so carbon precedes nitrogen), so the value carries meaning
        and not only the class. V2 ordered them by Python hash value, so an output order derived from
        these ranks does not agree numerically with a V2 one.
        """
        self._require_clean()
        if self._order_cache is not None and self._order_gen == self._gen:
            return self._order_cache
        cdef dict out = self._order_dict(NULL)
        self._order_cache = out
        self._order_gen = self._gen
        return out

    cdef dict _order_dict(self, uint32_t *seed):
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef list numbers = self._numbers
        cdef dict out = {}
        cdef uint32_t *rank
        cdef Py_ssize_t classes
        cdef uint32_t i
        if n_atoms == 0:
            return out
        rank = <uint32_t *> PyMem_Malloc(<size_t> n_atoms * sizeof(uint32_t))
        if rank is NULL:
            raise MemoryError()
        try:
            with nogil:
                classes = compute_atoms_order(self._structure, rank, seed)
            if classes < 0:
                raise MemoryError('atom order refinement failed to allocate')
            for i in range(n_atoms):
                out[numbers[i]] = rank[i]
        finally:
            PyMem_Free(rank)
        return out

    def refined_order(self, dict seed):
        """Refine caller-supplied starting classes to a fixed point, as {n: rank}.

        `atoms_order` starts from the atom records; this starts from whatever distinctions the
        caller already has, which is what stereo perception needs -- it differentiates
        stereocentres, feeds the result back in, and repeats. Labels are compared for equality
        only, so they need not be dense or 1-based, but they must be non-negative and cover every
        atom.
        """
        self._require_clean()
        cdef uint32_t *buf
        cdef dict out
        if self._structure.header.atom_count == 0:
            return {}
        buf = self._seed_labels(seed)
        try:
            out = self._order_dict(buf)
        finally:
            PyMem_Free(buf)
        return out

    cdef uint32_t *_seed_labels(self, dict seed) except NULL:
        """A {n: label} seed as the per-atom array the refinement wants.

        The caller owns the block and frees it with PyMem_Free. One implementation, shared by
        `refined_order` and the automorphism helpers, so the label lookup -- and the KeyError a
        seed that misses an atom raises, and the OverflowError a negative label raises -- cannot
        drift between them.
        """
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef list numbers = self._numbers
        cdef uint32_t i
        cdef uint32_t *buf = <uint32_t *> PyMem_Malloc(_alloc_at_least(n_atoms) * sizeof(uint32_t))
        if buf is NULL:
            raise MemoryError()
        try:
            for i in range(n_atoms):
                buf[i] = <uint32_t> seed[numbers[i]]
        except:
            PyMem_Free(buf)
            raise
        return buf

    def automorphism_orbits(self, dict seed=None):
        """Symmetry orbits as {n: orbit_id}, 1-based, from the automorphism group.

        Two atoms share an orbit when some automorphism of the molecule maps one onto the other.
        That is strictly stronger than sharing an `atoms_order` rank: refinement can call two
        atoms alike that no automorphism relates (regular graphs do this), so orbits are a
        subdivision of the ranks and never the other way round. Stereo perception needs the
        orbits: two ligands are interchangeable only if a symmetry really swaps them.

        Orbit numbers are dense and 1-based; beyond telling which atoms share an orbit they carry
        no meaning, and which orbit got which number is not stable. `seed` seeds the underlying
        colouring exactly as `refined_order` does, which lets a caller declare atoms distinct by
        hand -- the group then has to respect that.

        Raises `AutomorphismBudgetExceeded` when the search ran out of nodes. A truncated search
        returns orbits that may be FINER than the truth, which reads exactly like a right answer
        and would invent stereocentres, so there is nothing safe to degrade to.
        """
        self._require_clean()
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef list numbers = self._numbers
        cdef uint32_t *labels = NULL
        cdef uint32_t *orbits = NULL
        cdef uint32_t flags = 0
        cdef uint32_t i
        cdef dict out = {}
        if n_atoms == 0:
            return out
        if seed is not None:
            labels = self._seed_labels(seed)
        orbits = <uint32_t *> PyMem_Malloc(<size_t> n_atoms * sizeof(uint32_t))
        if orbits is NULL:
            PyMem_Free(labels)
            raise MemoryError()
        try:
            # The partition comes back from the search itself, which is the only thing the search
            # returns. The permutations it found are not available here and must not be: a bounded
            # sample of the group unioned into a partition under-reports symmetry.
            mol_automorphisms(self._structure, labels, orbits, &flags)
            if flags & CANON_BUDGET_EXCEEDED:
                raise AutomorphismBudgetExceeded(
                    'symmetry search exceeded its node budget (%d per pair, %d per call); the '
                    'orbits it reached may be finer than the true ones'
                    % (CANON_MAX_NODES_SEARCH, CANON_MAX_NODES_CALL))
            for i in range(n_atoms):
                out[numbers[i]] = orbits[i]
        finally:
            PyMem_Free(orbits)
            PyMem_Free(labels)
        return out

    def is_asymmetric(self, dict seed=None):
        """Does the molecule have no symmetry at all -- is its automorphism group trivial?

        True only when that is known. False means "not known to be asymmetric": a search that ran
        out of budget reports False rather than claim an asymmetry it did not prove, the consumer
        then does the full work it would have skipped, and that is the safe direction -- which is
        why this returns a bool where `automorphism_orbits` raises. `seed` behaves as in
        `automorphism_orbits`.
        """
        self._require_clean()
        cdef uint32_t *labels = NULL
        cdef uint32_t flags = 0
        if seed is not None and self._structure.header.atom_count:
            labels = self._seed_labels(seed)
        try:
            # NULL orbits: flags only, so the search stops at the first automorphism it finds. One
            # bit is all this reads.
            mol_automorphisms(self._structure, labels, NULL, &flags)
        finally:
            PyMem_Free(labels)
        return (flags & CANON_ASYMMETRIC) != 0

    def canonical_order(self, dict seed=None, *, uint32_t _node_budget=0):
        """The canonical atom order as {n: position}, 0-based and a permutation.

        Positions are the extremal labelling of the refinement tree, so they are a function of the
        molecule and not of the order its atoms were added in: rebuild the same molecule with its
        atoms in any order and the graph read back through these positions is the same graph. That
        is what `atoms_order` cannot do -- it stops at a partition, and on a symmetric molecule a
        partition leaves several atoms sharing a rank with nothing to say which comes first.

        Two atoms that some automorphism swaps have no position of their own: the pair of
        positions they occupy is fixed, which of them takes which is not. So the canonical FORM is
        unique -- hash it, compare it, serialise it -- while the labelling is unique only up to the
        automorphism group. `seed` seeds the underlying colouring exactly as `refined_order` and
        `automorphism_orbits` do, and any distinction it makes is respected here, which is how a
        caller pins a labelling down further than the structure alone can.

        Raises `AutomorphismBudgetExceeded` when the search ran out of nodes. There is deliberately
        no degraded answer: a truncated extremal search returns SOME labelling in place of THE
        labelling, which reads exactly like a right answer and would corrupt every hash built on
        it. `automorphism_orbits` raises for the mirror-image reason.

        `_node_budget` caps the refinement tree at that many nodes instead of the shipped
        1,000,000, and exists ONLY so that the paragraph above has a test -- no record small enough
        to run in a test suite can exhaust the real budget. Callers must not set it.
        """
        self._require_clean()
        cdef uint32_t n_atoms = self._structure.header.atom_count
        cdef list numbers = self._numbers
        cdef uint32_t *labels = NULL
        cdef uint32_t *order = NULL
        cdef uint32_t flags = 0
        cdef uint32_t i
        cdef dict out = {}
        if n_atoms == 0:
            return out
        if seed is not None:
            labels = self._seed_labels(seed)
        order = <uint32_t *> PyMem_Malloc(<size_t> n_atoms * sizeof(uint32_t))
        if order is NULL:
            PyMem_Free(labels)
            raise MemoryError()
        try:
            # Raises on a truncated search, and leaves `order` untouched when it does -- the flag
            # is read by no one here because the exception is the answer.
            # `stereo` True: the same search `canonical_bytes` runs, so the order this reports and
            # the order the canonical form is built on cannot be two different orders.
            if _node_budget:
                _canon_order(self._structure, labels, order, &flags, _node_budget, True)
            else:
                mol_canonical_order(self._structure, labels, order, &flags, True)
            for i in range(n_atoms):
                out[numbers[i]] = order[i]
        finally:
            PyMem_Free(order)
            PyMem_Free(labels)
        return out

    @property
    def atoms_order_classes(self):
        """How many distinct symmetry classes `atoms_order` found. Equal to atom_count when the
        molecule has no symmetry at all."""
        cdef dict order = self.atoms_order
        if not order:
            return 0
        return max(order.values())

    def kekule(self, aromatic_bonds=None, stated_h=None):
        """Rewrite this molecule's aromatic bonds as Kekule orders 1 and 2, in place.

        One of exactly two operations allowed to change a molecule's representation (`thiele` is
        the other), and it is always the caller's decision: nothing in the library kekulises
        behind your back, and no reader normalises what its input said.  Returns a
        `KekuleResult`: `.changed` is False on a molecule with no aromatic bonds and False on a
        second call.  The body is the module-level `kekule` in `_kekule.pxi` -- one operation, one
        name, reachable as a method on the thing it mutates and as a function for a caller holding
        the molecule at arm's length.
        """
        # the bare name is the module-level function and not this method: an unqualified lookup
        # inside a method body goes to module globals, never back through `self`
        return kekule(self, aromatic_bonds, stated_h)

    def derive_hydrogens(self, stated=None, *, fill_only=False):
        """Fill every derivable implicit hydrogen count from what this molecule is storing.

        THE READ-TIME PASS, and it is not one of the two representation-changing operations above: it
        changes no element, no charge, no bond and no count that was STATED.  It fills in the ones
        nobody has answered yet, so that a molecule which has just been parsed already has its
        hydrogens before `kekule()`, `standardize()` or `canonicalize()` is asked for anything.

        One algorithm, and the same one for every format -- the whole point of `_hydrogens.pxi`.
        `stated` is an iterable of the atoms whose count the record gave outright (an MDL
        `MRV_IMPLICIT_H`, a SMILES bracket, MRV's `hydrogenCount`); those are left exactly as they
        are.  Returns `{n: reason}` for the atoms it could not settle -- see
        `derive_implicit_hydrogens`, whose `HYD_*` codes those are.

        The atoms it cannot settle keep `H_UNKNOWN`, and for the aromatic pnictogen that is the
        pyrrole-versus-pyridine choice, which is `kekule()`'s to make and not a table's.  YOU DO NOT
        HAVE TO COME BACK FOR THOSE: `kekule()` runs this pass itself, in `fill_only` mode, on the
        atoms its own orders made derivable.  The mode is public anyway, for the one caller that has
        to do it by hand -- one holding an open edit scope across the kekulisation, where the orders
        are still pending and there is nothing to derive from until the scope closes.  Fill-only
        writes only where nothing is claimed, so it cannot undo a count the record gave or one a
        reader derived from something this pass cannot see.
        """
        return derive_implicit_hydrogens(self, stated, fill_only)

    def calc_implicit(self, uint32_t n):
        """Recompute and write ONE atom's implicit hydrogen count.  Returns what was derived.

        `derive_hydrogens()` above for a single atom, and the one to reach for after an edit that
        changed a bond order: it recomputes unconditionally, where the sweep's `fill_only` mode does
        not, so a count already stored is replaced rather than kept.

        `None` when nothing local can derive the count, and `H_UNKNOWN` is then what gets STORED --
        never zero, an atom whose hydrogens nobody can derive not being an atom with no hydrogens.
        That is also why this never raises: a metal the valence collection says nothing about has to
        survive a repair pass rather than stop it.  Two reasons for `None`, no row for this element in
        this charge and radical state, and the aromatic pnictogen whose count the ring decides;
        `kekule()` settles the second.
        """
        cdef object hydrogens
        hydrogens, _ = derive_implicit_hydrogen(self, n)
        self.set_hydrogens(n, H_UNKNOWN if hydrogens is None else hydrogens)
        return hydrogens

    def thiele(self):
        """Rewrite this molecule's Kekule bond orders as aromatic ones, in place.

        The other of the two operations allowed to change a molecule's representation, and the exact
        inverse spelling of `kekule` above -- same shape, same in-place contract, same reason to be a
        method as well as a function.  Returns a `ThieleResult`: `.changed`, and `.refused` listing the
        candidate systems declined with a `.log` line naming why.  The module-level `thiele` in
        `_thiele.pxi` documents the rule and its four proving molecules.

        NO WRITER CALLS IT.  A molecule holding Kekule orders is written Kekule, so `format(mol, 'A')`
        asks for aromatic BONDS in the string and does not aromatise the molecule to get them.  The
        aromatic spelling of a Kekule molecule is `mol.thiele()` first, and that is a mutation the
        caller performed rather than one a write performed behind them.
        """
        return thiele(self)

    def canonicalize(self, *, fix_tautomers=True, keep_kekule=False):
        """Bring this molecule to the representation two drawings of one compound share.  Changed?

        THE PASS TO RUN BEFORE DEDUPLICATING.  `canonical_bytes`, `__hash__` and `__eq__` are computed
        from what the molecule is STORING, so on their own they answer "same drawing" and not "same
        compound": `c1ccccc1O` and `C1=CC=CC=C1O` hash differently, and so do `[CH4]` and
        `[H]C([H])([H])[H]`.  This runs kekule, the group repairs, hydrogen implicification, thiele and
        the canonical placement of mobile hydrogens in the one order that is not a matter of taste, and
        after it those pairs agree.  A corpus deduplicated by hash WITHOUT it silently keeps duplicates.

        THE LAST GAP IS CLOSED.  Local shifts are unified by the group rules -- `Oc1ccccn1` and
        `O=c1cccc[nH]1` agree -- and a prototropic shift around a RING is unified by
        `standardize_isomers` below, so the two N-H forms of `Cc1cnc[nH]1` now agree too.  What remains
        outside this pass is what is outside the WORD: a ring-chain tautomer and a tautomer that moves a
        hydrogen between two separate molecules are different compounds by graph and stay different.

        Like `standardize()`, the body lives in `chython.chemistry` and arrives by registration
        rather than import, the return is a bool, and nothing in the library calls this for you.  Every
        record it wrote is on `self.log`, tagged with the stage that wrote it -- `log.by_stage('kekule')`
        through `log.by_stage('isomers')`.
        """
        return _canonicalize_fn()(self, fix_tautomers=fix_tautomers, keep_kekule=keep_kekule)

    def standardize_isomers(self):
        """Put every mobile hydrogen and charge where the canonical order says it goes.  Moved?

        THE STAGE THAT MAKES TWO ANNULAR TAUTOMERS OF ONE COMPOUND HASH EQUAL, and the last thing that
        stood between `canonical_bytes` and a sound compound identity.  The two N-H forms of
        `Cc1cnc[nH]1` are one compound stored two ways; after this they are stored one way.  It is
        narrower than the word "tautomer" suggests: a hydrogen moving between two heavy atoms ONE BOND
        apart is a local repair and `standardize()` already owns it, while a prototropic shift AROUND A
        RING moves the hydrogen and the double bonds together and no local rule sees both ends.

        REQUIRES THE AROMATIC FORM, which is why `canonicalize()` runs `thiele()` before this and not
        after: a mobile hydrogen is a property of the aromatic form, and in `CC1=CN=CN1` the double
        bonds have already said where it is.  A Kekule molecule is not refused here -- it simply has no
        mobile sites, and the honest answer to "did anything move" is `False`.

        WHAT IT CHOOSES BETWEEN IS PROVED, NOT SCORED.  A placement is admissible when the complete
        backtracking kekuliser finds a Kekule form for it, and among the admissible ones the winner is
        the one lowest in the canonical order of a skeleton with every candidate's hydrogen and charge
        STRIPPED.  That stripped frame is the whole trick: the ranks of the molecule as written depend
        on where the hydrogen already is, so they cannot be used to decide where it should go.

        Like `standardize()`, the body lives in `chython.chemistry` and arrives by registration rather
        than import, and the return is a bool.  Records land on `self.log` as `INFO`: both forms were
        valid molecules and neither was a defect to repair.
        """
        return _isomers_fn()(self)

    def check_valence(self):
        """`[(atom, verdict)]` for every atom the valence collection does not call valid.  A REPORT.

        Reads, never edits, and never raises -- so it is the way to ask whether a repair pass left an
        honest molecule behind, and the way to triage a corpus without dropping a row of it.  An empty
        list means every atom checked out.

        TWO VERDICTS, AND THEY ARE NOT THE SAME CLAIM.  `'violation'` says the collection describes this
        element in this charge and radical state and no row accepts what the molecule has -- a statement
        about the molecule.  `'unknown'` says the collection describes nothing there at all -- a gap in
        the collection, and no claim about the molecule whatsoever.  Reporting both under one word is how
        a coverage hole gets mistaken for bad input, so a caller filtering for defects filters for
        `'violation'`.

        KEKULISE FIRST IF YOU WANT AN ANSWER ABOUT A RING.  An atom carrying an aromatic bond is reported
        `'unknown'`: no valence row admits order 4, so there is nothing to check against and claiming a
        violation would be inventing one.  `mol.kekule()` then `mol.check_valence()` is the honest
        sequence, and it is what the repair passes' own tests use.

        Like `standardize()`, the body lives in `chython.chemistry` and arrives by registration rather
        than import -- the verdicts come from the core's generated valence tables, but the policy of
        which atoms get checked at all is standardization's.
        """
        return _valence_fn()(self)

    def implicify_hydrogens(self):
        """Fold ordinary hydrogen ATOMS into their neighbours' implicit COUNTS.  How many atoms went?

        A stage of `canonicalize()` above, exposed on its own because a caller may want the folding
        without the kekulisation and the group repairs -- reading an MDL record that spelled its
        hydrogens out is the common case.

        RETURNS AN INT, NOT A BOOL, and it is the number of hydrogen atoms REMOVED: `4` for
        `[H]C([H])([H])[H]`, not `1` for the one carbon touched and not `True`.  `0` is false, so
        `if mol.implicify_hydrogens():` reads as before.

        Five kinds of hydrogen are not a count and stay atoms -- an isotope, a charge, a radical, a
        bridging hydride, and an H2 whose neighbour is the other H2 -- and a bridging hydride is a
        `REFUSED` record rather than a raise.  A tetrahedral centre keeps
        its parity: `delete_atom` alone would racemise it.  `chython.chemistry._hydrogens` has the
        whole list with the reason for each.
        """
        return _hydrogens_fn(False)(self)

    def explicify_hydrogens(self):
        """Write every implicit hydrogen COUNT out as hydrogen ATOMS.  How many atoms arrived?

        The inverse of the above and NOT on the canonical path -- nothing canonical wants five atoms
        where one will do.  It is for consumers that need every hydrogen to be an addressable vertex:
        a depiction that labels them, a coordinate generator that places them.

        RETURNS AN INT: the number of hydrogen atoms ADDED.  The new atoms are UNMAPPED, and there is
        no keyword to number them -- a hydrogen this pass invented has no counterpart anywhere, so a
        map number would assert a correspondence that does not exist.  A caller who needs them mapped
        knows what it is mapping them to and assigns them itself.

        An atom whose implicit count is unknown gets nothing and a `LOST` record; there is no number to
        expand and inventing zero would answer a question the record never answered.
        """
        return _hydrogens_fn(True)(self)

    def split_salts(self, *, keep=()):
        """Cut every ionic cation-acceptor bond and put the charge on the two ends.  Anything cut?

        `CC(=O)O[Na]` becomes `CC(=O)[O-].[Na+]`.  The ATOM COUNT DOES NOT CHANGE and no component is
        deleted; nothing in the salt surface deletes one -- `decompose_salts()` reports instead.

        GENERAL OVER ALL 93 METALS `[M]` ACCEPTS, so a zinc or iron carboxylate splits exactly as an
        s-block one does.  What makes that safe is that the test is ALL-OR-NOTHING per cation atom: a
        dative bond, a neighbour that
        matches no acceptor row, an untabulated resulting charge or an implicit hydrogen on the cation
        refuses the whole atom and logs why.  `N[Pt](N)(Cl)Cl`, ferrocene and the metal carbonyls
        therefore come back intact rather than half-split.

        `keep=` takes row ids (`'salts:metal'`) and element symbols (`'Na'`).  The knowledge is
        `chython/chemistry/tables/salts.tsv`; `chython.chemistry._salts` has the reason for each
        refusal.
        """
        return _salts_fn('split_salts')(self, keep=keep)

    def decompose_salts(self):
        """What this record is, and what was drawn beside it.  A `SaltComposition` named tuple.

        `parents` is the compound, `counterions` and `solvates` are `{row id: equivalents}`, and
        `cations` counts lone cation atoms per element symbol -- separate because all 93 metals share
        one row id::

            smiles('NCC(=O)O.OC(=O)C(F)(F)F.O').decompose_salts()
            # parents=(smiles('C(CN)(=O)O'),) counterions={'salts:tfa': 1}
            # solvates={'salts:water': 1} cations={}

        NOTHING IS CHANGED AND NOTHING IS LOGGED: the work runs on a copy that is hydrogen-implicified,
        salt-split, neutralized and aromatized, so `CC(=O)O[Na]`, `CC(=O)[O-].[Na+]` and `CC(=O)O.[Na+]`
        all report one `Na`.  `parents` holds that form rather than the caller's drawing.

        A TABULATED SPECIES IS ONLY A COUNTERION WHEN SOMETHING ELSE IS THERE TO BE THE COMPOUND, so
        `smiles('CC(=O)O')` answers acetic acid as its own parent and `[Na+].[Cl-]` answers hydrochloric
        acid beside one `Na`.  `parents` is empty only for an empty molecule.
        """
        return _salts_fn('decompose_salts')(self)

    def neutralize(self, *, keep_charge=True):
        """Move every proton the acid/base table can from a cation onto an anion.  Anything moved?

        `[NH3+]CC(=O)[O-]` becomes `NCC(=O)O` and `C[NH3+].[Cl-]` becomes `CN.Cl`.  ONLY CHARGES AND
        IMPLICIT COUNTS ARE WRITTEN: no bond is cut and no component is deleted, which is what separates
        this from `split_salts` above.  `canonicalize()` runs it as a stage, so a
        zwitterion and its neutral drawing share a key; call it by hand on a molecule you are not
        canonicalizing.

        `keep_charge=True` moves protons in PAIRS and preserves the total charge exactly; a record that
        cannot be balanced comes back partly neutral rather than not at all.  NOTHING OVERSHOOTS ZERO in
        either mode: a component is only ever taken closer to neutral, so nitrate takes one proton and
        sulfate two.

        The knowledge is `chython/chemistry/tables/acids.tsv`, whose `h` primitive reads implicit
        hydrogens -- `implicify_hydrogens()` first if the molecule carries hydrogen atoms.
        """
        return _protomers_fn()(self, keep_charge=keep_charge)

    def functional_groups(self):
        """`{name: count}` for every functional group of `chython/reactions/tables/functional.tsv`.

        A METHOD AND NOT A CACHED PROPERTY.  A cached property on a mutable container survives an edit
        session that adds a nitrogen, and the enumerators below read this, so a stale cache would
        silently change which templates are tried.  The count is a question about the graph as
        it is now, and a call says so.  The enumerators compute it once per call and reuse it
        internally, so nothing pays for it twice in one enumeration.
        """
        return _reactions_fn('functional_groups')(self)

    def functional_group_hits(self):
        """Every functional group this molecule carries, in `functional.tsv`'s order, each with its id.

        A `GroupHit(id, name, count)` per group.  `functional_groups()` is this folded to `{name: count}`;
        the id is here because a consumer that stores group membership stores ids, an id being the
        identity a corpus edit does not move.
        """
        return _reactions_fn('functional_group_hits')(self)

    def react(self, *others, reaction=None):
        """Enumerate reactions between this molecule and the others, one per distinct outcome.

        Yields `EnumeratedReaction(name, reaction, rule_id)`.  THE WHOLE ENUMERATION SURFACE: one method
        over one `tables/reactions.tsv`, oxidations, reductions and interconversions included.  A row's
        slot count is not an arity filter, so a single-reactant row asks for no method of its own.

        A ROW IS OFFERED EVERY INPUT AT ONCE, so the argument order does not pick which molecule fills
        which slot: `a.react(b)` finds a coupling whichever of the two is the electrophile, and a mixture
        handed in as one container works.  What an outcome must satisfy instead is that EVERY input was
        touched, so `acid.react(amine, toluene)` yields nothing rather than an answer ignoring the
        toluene.  An untouched *component* of a touched input is a different thing and survives, which is
        what keeps a counter-ion from vanishing.

        `mol.react()` WITH NO PARTNER is the single-molecule question -- every oxidation, reduction and
        interconversion the corpus has.

        `reaction=` restricts to rows of that name (`'suzuki'`), and an unknown name raises rather than
        enumerating nothing.  The knowledge is `chython/reactions/tables/`.
        """
        return _reactions_fn('react')(self, others, reaction)

    def __matmul__(self, other):
        """`mol @ other` is `mol.react(other)`; `mol @ [b, c]` is `mol.react(b, c)`.

        Defined here rather than in `chython.reactions` because a special method resolves through the
        type's slot and a `cdef class` cannot be extended from outside -- see `_set_reactions_fns`.

        `~mol` IS NOT DEFINED: the single-molecule question is `mol.react()` with no partner.  An
        operator whose docstring can only say "same as the method" is a second spelling, not a shorthand.
        """
        if isinstance(other, (list, tuple)):
            return self.react(*other)
        return self.react(other)

    def protective_groups(self):
        """`{name: count}` for every protecting group of `chython/reactions/tables/protective.tsv`.

        A METHOD AND NOT A CACHED PROPERTY, for `functional_groups()`'s reason: a cache on a mutable
        container goes stale, and `deprotect()` reads this.

        THE COUNT IS CLAIMS, NOT MATCHES, and that is what makes the answer agree with `deprotect()`.
        The rules are tried most specific first and each site goes to the first rule that reaches it, so
        a Boc-protected alcohol reports one `hydroxyl_boc` and NOT the `hydroxyl_tbu` whose pattern also
        fits inside it.  A bis-Boc diamine reports two, because those are two sites and not two readings
        of one.
        """
        return _reactions_fn('protective_groups')(self)

    def protective_group_hits(self):
        """Every protecting group this molecule carries, most specific first, each with its id.

        A `GroupHit(id, name, count)` per group.  `protective_groups()` is this folded to `{name: count}`;
        the id is here because a consumer that stores group membership stores ids, an id being the
        identity a corpus edit does not move.
        """
        return _reactions_fn('protective_group_hits')(self)

    def deprotect(self, *names, protects=None, partial=False):
        """Enumerate deprotections of this molecule.  ALWAYS AN ITERATOR.

        Yields `EnumeratedDeprotection(names, reaction, rule_ids)`, where `reaction` has this molecule as
        its single reactant and the stripped molecule as its products.  A deprotection is a reaction, so
        it is reported as one and nothing here touches the caller's container.

        By default exactly one outcome, the full strip -- so `next(mol.deprotect(), None)` is the one-shot
        call and `None` means nothing was protected.  `partial=True` yields every non-empty subset of the SITES
        found, LARGEST FIRST, so element 0 is still the full strip and `2^k - 1` for *k* sites follow it.
        The unit is the site and not the rule because chemistry is not deterministic: a reagent that could
        cleave every Boc does not therefore cleave every Boc, incomplete cleavage is ordinary, and an
        N,N-di-Boc amine's second Boc is harder than its first -- so `R-N(Boc)2` answers `R-NHBoc`.  Two
        subsets giving the same products are one outcome, so a symmetry answers once.  Generated lazily, so
        taking the first few costs the first few.

        `*names` selects rules by name (`mol.deprotect('amine_boc')`) and `protects=` by what they reveal
        (`'amine'`, or several).  Neither changes what SHADOWS what -- specificity is computed over the
        whole table on every pass -- so asking for tert-butyl ethers on a Boc-protected alcohol yields
        nothing rather than cleaving the Boc down to a carbonate.

        HYDROGENS ON AN AROMATIC ATOM THE PATCH REACHED ARE `H_UNKNOWN`, as for any template application:
        the valence collection has no row for an atom holding aromatic bonds, so an underivable count is
        stored unknown rather than guessed.  `kekule()` then `chython.chemistry.calc_implicit` is the
        repair and it is the caller's, which is also why a deprotected N-aryl product does not compare
        `==` to a hand-written SMILES until that has run.
        """
        return _reactions_fn('deprotect')(self, names, protects=protects, partial=partial)

    def sticky_fragments(self, role=None, *, masked=None, bint hydrogens=False):
        """Every mono-attachment fragment this molecule exposes.  See `chython.reactions`."""
        return _reactions_fn('sticky_fragments')(self, role, masked=masked, hydrogens=hydrogens)

    def sticky_linkers(self, role_left=None, role_right=None, *, masked=None, bint hydrogens=False):
        """Every bi-attachment linker this molecule exposes.  See `chython.reactions`."""
        return _reactions_fn('sticky_linkers')(self, role_left, role_right, masked=masked,
                                              hydrogens=hydrogens)

    def standardize(self, *, fix_hydrogens=True, fix_tautomers=True):
        """Repair mis-drawn functional groups and metal-organic bonding in place.  Did it change?

        The body lives in `chython.chemistry`, which owns the 82 functional-group rules and the
        19 metal-organic ones and reads them out of `chython/chemistry/tables/`.  It is reached
        through the registration hook rather than an import because the dependency direction is
        `core <- chemistry`, and the core importing the package that imports it would invert that.
        The same shape as the kekuliser `molecule_to_inchi` reaches for.

        A REPAIR THE CALLER ASKED FOR.  Nothing in the library calls this: no reader standardizes
        what its input said and no writer standardizes what it is about to emit.  A molecule that
        parsed from a badly drawn record keeps that record's bonding until someone asks for this.

        `self.log` gets one record per patch applied and one per patch refused, whether or not anyone
        asked.  The return type is a bool, always.

        `fix_tautomers=False` withholds the 28 rules whose repair moves a hydrogen from one heavy
        atom to another -- enol to ketone, hydroxy-azine to amide, the ring-amidation family.  The
        rule table marks them; `chython.chemistry.standardize` documents why that has to be a column
        rather than something inferred, and this signature only forwards it.
        """
        return _standardize_fn()(self, fix_hydrogens=fix_hydrogens, fix_tautomers=fix_tautomers)

    def fix_resonance(self):
        """Collapse a separated charge pair back onto one atom where a neutral form exists.  Changed?

        `[O-][N+](=O)C` stays a nitro group -- no neutral form of pentavalent nitrogen exists -- while
        `C=[N+]([O-])C` written for an amide comes back `CC(=O)N` shaped.  The pass walks the ALTERNATING
        PATH between a `+` and a `-` and shifts the orders along it; `chython.chemistry._resonance` states
        which pairs qualify.

        NOT A STAGE OF `standardize()` OR `canonicalize()` and not called anywhere in the library.  A
        drawing's charge separation may be what the author meant, so removing it is a decision, and the
        decision is the caller's.  `bool` and one `self.log` record per shift, like every pass above.

        Registered by `chython.chemistry` rather than compiled in -- see `_set_resonance_fn`.
        """
        return _resonance_fn()(self)

    def layout2d(self, *, engine=None, force=False):
        """Compute a 2D layout and RETURN it as `{n: (x, y)}`, storing nothing.

        THE FORM A RENDERER WANTS: drawing must not change what it draws, so a layout computed for a
        picture is handed back rather than stored.  `clean2d()` is this plus the decision to keep it.

        `force=False` on a molecule that already `has_layout` hands back the STORED plane and computes
        nothing; `force=True` always recomputes.  `engine` overrides `chython.clean2d_engine` for this
        one call.

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_layout2d_impl, 'layout2d')(self, engine=engine, force=force)

    def clean2d(self, *, engine=None, force=False):
        """Compute a 2D layout and STORE it on this molecule.

        NOT ALWAYS A RECOMPUTATION: a molecule that already carries a non-degenerate plane
        (`has_layout`) is left alone, and `force=True` relays it unconditionally.  A no-op is the
        correct answer to "make sure this molecule has a layout".

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_clean2d_impl, 'clean2d')(self, engine=engine, force=force)

    def rescale2d(self):
        """Rescale the STORED coordinates to average bond length 0.825.  Did it rescale?

        0.825 is the scale the rest of the depiction geometry assumes, so a plane that came from a
        drawing editor is normalized by this rather than redrawn -- a redraw would throw the drawing
        away.  NOTHING BUT THE COORDINATES IS TOUCHED, and the plane is scaled about the origin, so
        every atom keeps its position relative to every other.

        False, and nothing stored, when there is no scale to read: a molecule with no coordinates, one
        with no bonds, or a plane collapsed tightly enough that dividing by its mean would be a
        singularity.

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_rescale2d_impl, 'rescale2d')(self)

    def scene(self, *, style=None, plane=None, overlays=(), log=None):
        """The drawing of this molecule as a backend-independent `Scene`, in molecule coordinates.

        The object a caller holds when it wants the geometry rather than a document -- to compose a
        grid, to place it beside another figure, or to serialize it more than once.  `depict()` is this
        plus `Scene.to_svg()`.

        `plane` draws a `{n: (x, y)}` layout the molecule does not carry; `style` defaults to
        the process default (`chython.depict.get_depict_style`), which is a default for this entry point
        and never a value a drawing function reads for itself.

        `overlays` is a sequence of overlay objects (`Highlight`, `AtomHalo`, `AtomField`, `BondScale`,
        `ValueLabels`) rendered in order onto the scene.  A colorbar is added when any overlay carries
        a colour scale.

        DRAWING STORES NOTHING.  A molecule with no layout is drawn against a temporary one and a
        `LogRecord` says so; `clean2d()` is the call for a layout that is kept.

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_scene_impl, 'scene')(self, style=style, plane=plane, overlays=overlays,
                                                log=log)

    def depict(self, *, style=None, plane=None, overlays=(), log=None):
        """This molecule as an SVG document.  `scene()` plus serialization, and nothing else.

        NOT CACHED, here or in `chython.depict`: a cached picture is a picture at ONE style, so the
        second style would silently get the first one's output.  Same arguments as `scene()`.

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_depict_impl, 'depict')(self, style=style, plane=plane, overlays=overlays,
                                                  log=log)

    def _repr_svg_(self):
        """Jupyter's hook.  The PROCESS DEFAULT style, because a notebook cell states none."""
        return _depict_fn(_depict_impl, 'depict')(self)

    def depict3d(self, uint32_t index=0):
        """Model `index` as an X3DOM document -- spheres for atoms, cylinders for bonds.

        A CONFORMER IS READ, NEVER GUESSED.  A molecule carrying no model is refused, where `depict()`
        falls back to a temporary plane and logs it: a 2D layout can be computed from the graph alone
        and a conformer cannot, so the fallback would be an invented geometry.  `has_3d` is the test
        and `chython.interop.conformers.generate_conformers` is the way to get one.

        Its rendering parameters are the module-level defaults in `chython/depict/x3dom.py` and are
        NOT `DepictStyle`, which is 2D throughout.

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_depict3d_impl, 'depict3d')(self, index)

    def view3d(self, uint32_t index=0, width='600px', height='400px'):
        """Model `index` as a Jupyter widget: `depict3d()` in a sized div that loads X3DOM.

        The viewer is the browser's, so the widget references x3dom.org and a notebook opened offline
        shows an empty box -- the document itself is complete and `depict3d()` is what to save.

        Registered by `chython.depict`, not implemented here -- see `_depict_fn`.
        """
        return _depict_fn(_view3d_impl, 'view3d')(self, index, width, height)

    # ------------------------------------------------------------------ toolkit conversion
    # The export half of `chython.interop`, whose dispatchers are the bodies -- see `_set_interop_fns`.
    # EXPORT ONLY, and the asymmetry is deliberate: `interop.rdkit(rd_mol)` reads a foreign object, and
    # there is no `self` to hang that on.  Keywords are forwarded rather than spelled out, because they
    # differ per toolkit and only `rdkit` takes any beyond `log`; the converter's own signature is the
    # one place they are documented, so a wrong one raises there and names the direction it reached.

    def to_rdkit(self, **kwargs):
        """This molecule as an RDKit `Mol`.  `chython.interop.rdkit` is the function behind it.

        `keep_mapping=True` by default and it means the ATOM-ATOM MAPPING: an unmapped molecule exports
        with a clean map field, so `MolToSmiles` of the result is a plain SMILES.  `keep_numbers=True`
        asks for chython's stable ids in that field instead, which is the label to match results back
        on.  `keep_hydrogens`, `keep_coordinates` and `absolute` are documented on `interop._rdkit`.
        """
        return _interop_fn('rdkit')(self, **kwargs)

    def to_indigo(self, **kwargs):
        """This molecule as an Indigo object.  `chython.interop.indigo` is the function behind it."""
        return _interop_fn('indigo')(self, **kwargs)

    def to_openbabel(self, **kwargs):
        """This molecule as an OpenBabel `OBMol`.  `chython.interop.openbabel` is behind it."""
        return _interop_fn('openbabel')(self, **kwargs)

    def to_cdk(self, **kwargs):
        """This molecule as a CDK `IAtomContainer`.  Starts a JVM through JPype on first use."""
        return _interop_fn('cdk')(self, **kwargs)

    def to_cdpkit(self, **kwargs):
        """This molecule as a CDPKit molecule.  `chython.interop.cdpkit` is behind it."""
        return _interop_fn('cdpkit')(self, **kwargs)

    @property
    def iupac(self):
        """The IUPAC name openclatura writes for this molecule, or `None` when it cannot name it.

        A PROPERTY AND NOT CACHED: `functools.cached_property`
        needs a `__dict__` a `cdef class` does not have, and a name cached on a mutable container
        outlives the edit session that makes it wrong.  It is not free -- an RDKit export plus
        openclatura per read -- so a caller naming a corpus should hold the string.

        openclatura is optional and requires Python >= 3.11; its absence raises `ImportError`.
        """
        return _interop_fn('iupac')(self)

    @property
    def inchi(self):
        """The standard InChI of this molecule.

        `mol.smiles` for the other identifier, and the same division of labour: the property is the
        plain answer, and the options live on the function.  `molecule_to_inchi(mol, standard=False,
        options=...)` is where a non-standard string or an InChI flag comes from, exactly as
        `format(mol, spec)` is where a non-default SMILES does.

        NOT CACHED, for `iupac`'s reason -- a `cdef class` has no `__dict__` for `cached_property`, and
        a string cached on a mutable container outlives the edit that makes it wrong.  Raises
        `ImportError` when libinchi is not loaded; `chython.inchi_library_loaded()` is the pre-check.
        """
        return molecule_to_inchi(self)

    @property
    def inchikey(self):
        """The standard InChIKey of this molecule: 27 characters, and nothing reads one back.

        `self.inchi` hashed, so the same caveats hold -- not cached, and `ImportError` without
        libinchi.  A key is a lossy fingerprint of a string and not the string, so it deduplicates a
        corpus and answers no question about structure.
        """
        return molecule_to_inchikey(self)

    def clean_isotopes(self):
        """Drop every isotope label, in place.  Did the molecule carry one?

        HERE AND NOT IN `chemistry` BECAUSE IT NEEDS NO RULE TABLE.  `standardize` above asks what a
        drawing meant and answers out of 101 rows; this asks nothing.  The dividing line for the
        standardization pack is knowledge, not mutability -- an operation that reads no table is a
        container method, which is also why `clean_stereo` is one.

        DEUTERIUM AND TRITIUM GO TOO.  A caller who wanted the heavy hydrogens kept has not asked for
        this; `implicify_hydrogens` is where that distinction lives, and it keeps an isotopic hydrogen
        as an explicit atom precisely because dropping the label would lose the fact.  Here losing the
        fact is the request.

        THE STEREO CONSEQUENCE IS REAL AND IS HANDLED, and it is the one place a mutation cannot
        leave stereo to `_apply`.  The journal's apply re-bases each configured parity against the new
        FRAME -- the anchor's directions -- and an isotope is in no frame, so a parity whose whole
        justification was the label survives the drop and starts lying.  Measured:
        `C[C@H](F)[13CH3]` differs at the two methyls only by the label, and dropping it leaves
        `[C@H](C)(F)C`, a configuration on an atom that has none.  So a drop is followed by
        `validate_stereo()`, which is the one implementation of that question this tree has; writing
        the narrow version here -- clear only what this call stranded -- is not even correct, since
        the label that justified the parity above sits on a NEIGHBOUR.

        The cost of borrowing it is that it also clears a parity that was already unjustifiable before
        the call, on a molecule that arrived that way.  That is stated rather than fixed: the caller
        asked for a repair, and the alternative is a second, worse `validate_stereo`.

        Returns False on a molecule with no isotope anywhere, and that answer is a PURE READ -- no
        clone, no journal, no generation bump, `shares_arena_with` still True.
        """
        cdef atom_t *atoms
        cdef uint32_t i
        cdef list labelled = []
        self._require_clean()
        atoms = self._structure.atoms()
        for i in range(self._structure.header.atom_count):
            if atoms[i].isotope:
                labelled.append(self._numbers[i])
        if not labelled:
            return False
        with self.edit():
            for i in labelled:
                self.set_isotope(i, 0)
        self.validate_stereo()
        return True

    def remove_coordinate_bonds(self, *, keep_stranded_hydrogens=True):
        """Delete every dative bond (order 8), in place.  How many went?

        THE COMPLEMENT OF WHAT `standardize` DOES.  The metal-organic half of the rule table CREATES
        these bonds -- a phosphine ligand becomes `C[P](~[Fe])(C)C` rather than a phosphonium -- so
        this is the call for a caller who wants the coordination sphere taken apart again: a metal
        carbonyl becomes an iron and three carbon monoxides, in one container, as separate components.
        STRANDING THE METAL IS THE POINT and not an accident to guard against, which is why the guard
        below is about hydrogen and nothing else.

        `keep_stranded_hydrogens=True` keeps the order-8 bonds of a hydrogen that has NO covalent bond
        at all -- a bridging hydride, a hydrogen held only by coordination -- because deleting them
        leaves a disconnected `[H]` that names nothing.  An ordinary hydrogen-bond donor is untouched
        by the guard and its contact IS deleted: that hydrogen keeps the covalent bond it came with.
        The test is exactly "has this hydrogen any non-order-8 neighbour", not "is it terminal".

        STEREO NEEDS NOTHING HERE, unlike `clean_isotopes` above: a bond is part of its anchor's frame,
        so the journal's apply re-bases or drops each configured parity against the frames the deletion
        left.  Measured on `C[C@](F)(Cl)~[Fe]`, whose centre the READER already declines.

        Returns 0 on a molecule with no dative bond, as a pure read.
        """
        cdef Structure s
        cdef atom_t *atoms
        cdef uint32_t *ptr
        cdef halfedge_t *edges
        cdef uint32_t i, k, to, n_atoms
        cdef bint covalent
        cdef set stranded = set()
        cdef list doomed = []
        self._require_clean()
        s = self._structure
        atoms = s.atoms()
        ptr = csr_ptr(s)
        edges = csr_edges(s)
        n_atoms = s.header.atom_count
        if keep_stranded_hydrogens:
            for i in range(n_atoms):
                if atoms[i].element != 1:
                    continue
                covalent = False
                for k in range(ptr[i], ptr[i + 1]):
                    if edges[k].order != 8:
                        covalent = True
                        break
                if not covalent:
                    stranded.add(i)
        for i in range(n_atoms):
            for k in range(ptr[i], ptr[i + 1]):
                to = edges[k].to
                # each bond once, and `i` is the lower slot -- `bonds()`' promise, same test
                if to > i and edges[k].order == 8 and i not in stranded and to not in stranded:
                    doomed.append((self._numbers[i], self._numbers[to]))
        if not doomed:
            return 0
        with self.edit():
            for i, k in doomed:
                self.delete_bond(i, k)
        return len(doomed)

    cdef bytes _identity(self):
        """`canonical_bytes`, CACHED, keyed on the same generation counter `atoms_order` uses.

        `__hash__` and `__eq__` both come through here and the canonical order is ~97% of the cost
        of the whole computation, so an uncached hash would make `set(molecules)` over a hundred
        thousand records a twenty-second mystery in a profile.  The cache is invalidated by the
        container's generation counter rather than by a second notion of staleness: `_gen` is bumped
        by `_apply` and by every in-place writer, and `_require_clean` refuses a read while a journal
        is pending, so there is no window where a stale row can be returned.  It over-invalidates --
        a coordinate or a map-number change bumps `_gen` and neither reaches the canonical form -- and
        that is the safe direction to be wrong in.

        THE TWO STORES ARE ORDERED, and the order is load-bearing under `freethreading_compatible`:
        the bytes go down BEFORE the generation they belong to.  Torn the other way -- generation
        first -- a second thread would see a fresh generation beside a stale row and return the wrong
        molecule's canonical form.  Torn this way the worst case is a reader that sees the new bytes
        beside the old generation, misses the cache and recomputes, which is a wasted 200
        microseconds and not a wrong answer.  Two threads racing to fill it compute the same value
        from the same immutable arena, so whichever store wins is the same store.  A thread mutating
        this molecule while another reads it is unsafe for reasons that have nothing to do with this
        cache -- the arena itself is being replaced -- and is not made safe here.
        """
        self._require_clean()
        if self._identity_cache is not None and self._identity_gen == self._gen:
            return self._identity_cache
        cdef bytes out = mol_identity_bytes(self._structure)
        self._identity_cache = out
        self._identity_gen = self._gen
        return out

    @property
    def canonical_bytes(self):
        """The canonical form of this molecule as bytes: equal bytes mean the same compound.

        THE IDENTITY, as opposed to `_union_feature_words`, which is the (private, lossy) screen.
        Built from the extremal canonical labelling -- so it is a function of the molecule and not of
        the order its atoms were added in -- and it carries the graph, the elements, isotopes,
        charges, radicals, implicit hydrogen counts, bond orders, aromatic bits, and one parity digit
        per position.  `mol_identity_bytes` documents the construction and its one residual.

        This is what `__eq__` and `__hash__` are built on, so the three answer alike by construction
        and the basis is inspectable rather than implied.  Cached; see `_identity`.

        Raises `AutomorphismBudgetExceeded` on a truncated canonical search, like `canonical_order`.
        """
        return self._identity()

    @property
    def _union_feature_words(self):
        """LOSSY. The OR of every atom's four feature words -- a 32-byte screen, and PRIVATE because
        nothing outside a test should read it.

        NOT AN IDENTITY AND NOT A "SIGNATURE" IN THE SENSE THE REST OF CHYTHON USES THAT WORD.
        `docs/reactions.rst` writes `rxn1 == rxn2  # True if same canonical signature`, meaning a
        canonical STRING, and that collision is exactly why this property kept being read as an
        identity.  It is not one.  Measured: propane, butane and pentane all return
        `(866942928268820480, 4611686018427387904, 9223372174432275457, 72198331526283521)`, because
        an OR over atoms cannot count them.  Enhanced stereo is absent from the words entirely --
        ABS, AND1 and OR1 on the same centre give one identical value, so a racemate and a single
        enantiomer are indistinguishable here.  `canonical_bytes` is the identity.

        Underscored rather than renamed, because a public name invites a caller and there should not
        be one: the words' whole job is to be the right-hand side of `query_may_match`'s gate, which
        reads `structure_features()` in C and never comes through this property.  Zero non-test
        consumers exist.  The legitimate test use is ruling F78's invariant -- that a to_bytes round
        trip and `refresh_parity_features` both reproduce what a fresh derivation builds -- where the
        union row is genuinely the object under test.  `query.header.signature[4]` on the C side
        keeps its name; there the meaning is unambiguous from its two lines of context.
        """
        self._require_clean()
        cdef uint64_t *f = structure_features(self._structure)
        return (f[0], f[1], f[2], f[3])

    def features_of(self, uint32_t n):
        """Four screening words for atom `n`: (element/bonds, heavy-element/radical,
        counts/charge/isotope, hybridization/stereo/rings)."""
        self._require_clean()
        cdef uint32_t i = self._index_of[n]
        cdef uint64_t *f = structure_features(self._structure)
        cdef uint32_t base = 4 + 4 * i
        return (f[base], f[base + 1], f[base + 2], f[base + 3])

    @property
    def aromatic_bond_count(self):
        """How many bonds are stored with order 4. The representation state, exactly.

        Not a flag, and deliberately not a three-way KEKULE/AROMATIC/MIXED enum. A stored flag can
        contradict the bond orders and then there are two truths; a count is recomputed from the
        bonds at every point where a bond can have changed, so it cannot. And a molecule may
        genuinely hold one type-4 ring beside one alternating ring -- two differently-written inputs,
        faithfully kept -- whereas calling that state "mixed" would require deciding whether the
        alternating ring OUGHT to be aromatic, which is a perception question no stored state can
        answer honestly.

        THE INVARIANT THIS EXISTS TO MAKE TESTABLE: no operation in the core changes this number
        except `kekule()` and `thiele()`. Nothing sanitises, normalises, kekulises or aromatises
        behind the caller's back -- not a parser, not a writer, not a perception pass, not InChI.
        """
        self._require_clean()
        return self._structure.aromatic_bond_count

    @property
    def is_kekule(self):
        """True when no bond is stored aromatic. Exactly `aromatic_bond_count == 0`."""
        self._require_clean()
        return self._structure.aromatic_bond_count == 0

    @property
    def persistent_len(self):
        self._require_clean()
        return self._structure.header.persistent_len

    @property
    def total_len(self):
        self._require_clean()
        return self._structure.total_len

    @property
    def persistent_view(self):
        self._require_clean()
        cdef uint32_t n_persistent = self._structure.header.persistent_len
        return memoryview(<bytes> self._structure.buffer[:n_persistent])

    def pack(self, *, bint compressed=True, drop=None, version=None):
        """One pach record.  See `chython.core.pach_dump`, which this forwards to.

        `version` is None for the current layout -- 3 with coordinates, 4 without -- 3 or 4 to state
        that layout outright, or 2 for the legacy one.  `3` asks for the coordinates the molecule has,
        so an undrawn molecule and a `drop=['coordinates']` caller both get version 4; `4` writes no
        coordinate block even for a drawn molecule.  pach is small and lossy; `to_bytes` is the arena
        verbatim and lossless.
        """
        return pach_dump(self, compressed=compressed, drop=drop, version=version)

    @staticmethod
    def unpack(data, *, compressed=None):
        """A molecule from either serialised form: a legacy pach record, or `to_bytes` output.

        NOTHING HAS TO BE DECLARED ABOUT THE BUFFER, because all three things it can be are
        distinguishable by their first byte. A pach record's is its format version, one of 0, 2, 3 and
        4 for a molecule and 1 or 5 for a reaction; the arena format's is 0x33, the low byte of its
        little-endian magic; and a zlib stream's low nibble is its compression method, always 8, so
        0, 2, 3, 4, 5 and 0x33 are six bytes no zlib header can spell. A caller holding bytes out of a
        store therefore does not have to know which era wrote them or whether anybody compressed them.

        `compressed` is for a caller who wants to be TOLD rather than accommodated: `True` insists the
        buffer be compressed and `False` insists it not be, and either is a ValueError when the buffer
        disagrees. The default sniffs.

        THIS IS AN ANSWER BOUNDARY AND IT RAISES. `ValueError` names what was wrong with the record,
        with every problem the decoder found appended -- including the ones it recovered from, since a
        caller who cannot have a molecule is owed the whole story. A caller who wants the recovered
        molecule and the complaints instead of an exception wants `chython.core.pach_load`.
        """
        cdef bytes raw = bytes(data)
        cdef object mol
        cdef list problems
        cdef bint looks_raw
        if not len(raw):
            raise ValueError('the buffer is empty')
        # 1 and 5 are the reaction pach versions, recognised here only so that a caller who fed a
        # stored reaction buffer to the molecule door is told which door it wanted. See `reaction.py`.
        looks_raw = (raw[0] == 0 or raw[0] == 1 or raw[0] == 2 or raw[0] == 3
                     or raw[0] == 4 or raw[0] == 5 or raw[0] == 0x33)
        if compressed is True and looks_raw:
            raise ValueError('compressed=True was stated and the buffer begins with %d, which is a '
                             'raw record and not a zlib header' % raw[0])
        if compressed is False and not looks_raw:
            raise ValueError('compressed=False was stated and the buffer begins with %d, which is '
                             'neither a pach version nor the arena magic' % raw[0])
        if not looks_raw:
            try:
                raw = zlib.decompress(raw)
            except Exception as err:
                raise ValueError('the buffer begins with %d, so it is neither a raw record nor a '
                                 'readable zlib stream: %s' % (raw[0], err))
            if not len(raw):
                raise ValueError('the buffer decompressed to nothing')
        if raw[0] == 1 or raw[0] == 5:
            raise ValueError('byte 0 is %d, which is a REACTION pach record and not a molecule; use '
                             'ReactionContainer.unpack' % raw[0])
        if raw[0] == 0 or raw[0] == 2 or raw[0] == 3 or raw[0] == 4:
            mol, problems = pach_load(raw, compressed=False)
            if mol is None:
                raise ValueError('this is not a readable pach record: %s' % '; '.join(problems))
            if problems:
                raise ValueError('this pach record is damaged: %s. pach_load() returns the molecule '
                                 'that could be recovered from it along with these problems'
                                 % '; '.join(problems))
            return mol
        return MoleculeContainer.from_bytes(raw)

    def pach(self, *, bint compressed=True, drop=None, version=None):
        """chython 2's name for `pack`, and the same record byte for byte.

        The three keywords chython 2 also took are gone rather than accepted and ignored: `check=` is
        not a choice here (this release refuses by field name, which `drop=` waives), `order=` states
        an atom order pach has never carried, and `skip_labels_calculation=` names a step this arena
        does not have.  Each is a `TypeError`, because a silently ignored `check=False` would promise
        a refusal was waived and let `pack` raise anyway.
        """
        return pach_dump(self, compressed=compressed, drop=drop, version=version)

    @staticmethod
    def unpach(data, *, compressed=None):
        """chython 2's name for `unpack`, with its behaviour: an answer boundary that raises."""
        return MoleculeContainer.unpack(data, compressed=compressed)

    def to_bytes(self):
        """The arena's persistent prefix, verbatim. Not the pach format -- see pack()."""
        self._require_clean()
        cdef uint32_t n_persistent = self._structure.header.persistent_len
        return <bytes> self._structure.buffer[:n_persistent]

    def __reduce__(self):
        # THE CANONICAL-FORM CACHE IS NOT PICKLED, and must never be: it goes through `to_bytes`,
        # which carries the arena and nothing derived.  A pickled cache surviving into a process where
        # the canonicalisation has changed would be a silent wrong answer with no way to notice it,
        # so it is recomputed on unpickle -- 200 microseconds once, against an unfindable defect.
        return (_from_bytes, (self.to_bytes(), self._meta))

    @staticmethod
    def from_bytes(data):
        cdef const unsigned char[::1] view = data
        cdef Structure structure
        cdef atom_t *atoms
        cdef uint8_t *par
        cdef uint32_t *norm_ptr
        cdef halfedge_t *norm_edges
        cdef uint32_t i, k
        cdef bint is_wedge_narrow
        cdef bint adopt = False
        # READ BEFORE THE CALL.  `structure_from_bytes` normalises an older buffer's header into the
        # current version in place, so afterwards the version byte no longer says where the parities
        # are.  Guarded on the length because boundscheck is off in this module and the call below is
        # what refuses a buffer too short to hold a header.
        cdef uint16_t src_version = 0
        if view.shape[0] >= 6:
            src_version = <uint16_t> view[4] | (<uint16_t> view[5] << 8)
        structure = structure_from_bytes(<const char *> &view[0], view.shape[0])

        # AN OLDER BUFFER KEEPS ITS PARITIES IN `atom_t.flags`, bit 1 the value and bit 7 the
        # "configured" bit; from STRUCT_VERSION 5 they are a byte per atom in SEG_PARITY.  Both steps
        # below read version-4 storage rather than this build's, which is why the bits are spelled as
        # literals instead of through an accessor.
        #
        # First, bit 1 set with bit 7 clear is a writer that spent bit 1 on a drawing flag rather than
        # on a parity value, and the wedge segment discriminates the two readings -- `halfedge_t.wedge`
        # lives in SEG_CSR_EDGE, part of the persistent prefix, so it is valid before rebuild_derived:
        #
        #   (a) the atom is NOT the narrow end of any nonzero wedge, so bit 1 was a stored parity
        #       direction: set bit 7;
        #   (b) the atom IS the narrow end of a nonzero wedge, so the wedge set bit 1: clear it,
        #       leaving parity_of = 0 and stereo_of = False.
        #
        # The wedge segment is sound as the discriminator because a wedge does not configure a parity,
        # so OP_SET_WEDGE writes neither parity bit and bit-1-alone is unreachable from this build's
        # writers.
        if src_version < STRUCT_VERSION_V5:
            atoms = structure.atoms()
            norm_ptr = csr_ptr(structure)
            norm_edges = csr_edges(structure)
            for i in range(structure.header.atom_count):
                if (atoms[i].flags & 0x02u) and not (atoms[i].flags & 0x80u):
                    is_wedge_narrow = False
                    for k in range(norm_ptr[i], norm_ptr[i + 1]):
                        if norm_edges[k].wedge:
                            is_wedge_narrow = True
                            break
                    if is_wedge_narrow:
                        atoms[i].flags &= <uint8_t> ~0x02u   # (b): clear geometry bit
                    else:
                        atoms[i].flags |= <uint8_t> 0x80u    # (a): promote to configured parity

            # Then the flags are MOVED into the segment, and that move is one-directional: the flags are
            # the only copy the buffer has, so an atom missed here loses its configuration outright.
            # Case (a) above sets bit 7, so the scan has to run after it.
            #
            # A FRESH BLOCK, not a patch: the incoming buffer has no room for the segment, the
            # persistent block is laid out once and `structure_from_bytes` copied that buffer verbatim.
            # Only for a buffer that states something -- a stereo-free record pays one scan of the flags
            # and no allocation.
            for i in range(structure.header.atom_count):
                if atoms[i].flags & 0x80u:
                    adopt = True
                    break
            if adopt:
                structure = structure_with_parity(structure)
                atoms = structure.atoms()
                par = structure_parities(structure)
                for i in range(structure.header.atom_count):
                    if atoms[i].flags & 0x80u:
                        par[i] = 2 if atoms[i].flags & 0x02u else 1
                    # THE MOVE IS COMPLETED HERE.  From version 5 these two bits are reserved and a
                    # reader refuses a buffer that sets them, so leaving them would make this
                    # molecule's own `to_bytes` unreadable by its own `from_bytes`.
                    atoms[i].flags &= <uint8_t> ~ATOM_FLAGS_RESERVED

        rebuild_derived(structure)
        # Ruling F60: rebuild_derived appends derived segments and so REALLOCATES the arena.
        # `atoms` above points into the pre-rebuild buffer, which may now be freed, so it is
        # re-fetched here rather than reused.  Reading stable ids through the stale pointer
        # silently corrupted `_numbers` on roughly three quarters of the round trips of a
        # 100-atom molecule; see test_pack.test_wedge_round_trip_keeps_stable_ids.
        atoms = structure.atoms()

        cdef uint32_t n
        cdef uint32_t high = 0
        cdef list numbers = []
        cdef dict index_of = {}
        for i in range(structure.header.atom_count):
            n = atoms[i].n
            numbers.append(n)
            index_of[n] = i
            if n > high:
                high = n

        cdef MoleculeContainer mol = MoleculeContainer.__new__(MoleculeContainer)
        mol._structure = structure
        mol._numbers = numbers
        mol._index_of = index_of
        mol._next_id = high + 1
        mol._first_pending = high + 1
        # THE NARROWING HAPPENED IN THE ARENA, WHICH HAS NO LOG.  `structure_from_bytes` works in
        # bytes and the version byte no longer says what came in, so this is the one place that can
        # both tell and be heard: `src_version` was read before the call for exactly this reason.
        if src_version < STRUCT_VERSION and structure_conformer_count(structure):
            mol._log_event('container:conformer-narrowed', 'read',
                           '%d conformer(s) came from a version-%d buffer, whose record carries '
                           'three words this build does not model; each model keeps its coordinates '
                           'and the number the file gave it'
                           % (int(structure_conformer_count(structure)), int(src_version)),
                           mc_lost())
        return mol

    def may_contain(self, MoleculeContainer other not None):
        """
        Cheap necessary-condition screen for `other` being a substructure of `self`.

        False is definitive: no substructure mapping of `other` into `self` exists, and
        a caller may skip the search. True means only that the search is worth running.
        Compares element composition, formal charge, isotope, radical state, bond orders
        and the presence of ring bonds -- never degree, hydrogen counts, heteroatom
        count, hybridization or ring descriptors, none of which survive embedding.
        """
        self._require_clean()
        other._require_clean()
        return bool(sig_contains(structure_features(self._structure),
                                 structure_features(other._structure)))

    def atoms_of_element(self, uint32_t number):
        """Return a tuple of stable ids for all atoms with the given atomic number, in index order."""
        self._require_clean()
        if number < 1 or number > 118:
            raise ValueError('number must be an atomic number in 1-118')
        cdef uint32_t begin = element_bucket_begin(self._structure, number)
        cdef uint32_t end = element_bucket_end(self._structure, number)
        cdef uint32_t *idx = structure_element_index(self._structure) + 120
        cdef uint32_t k
        cdef list out = []
        for k in range(begin, end):
            out.append(self._numbers[idx[k]])
        return tuple(out)

    @property
    def is_radical(self):
        """True when ANY atom carries the radical flag.

        A fold over the atom table and not a stored total, for the reason `__int__` gives: a stored
        flag is a second truth that a `set_radical` can contradict.  An empty molecule is False.

        NOT A COUNT AND NOT A MULTIPLICITY.  A biradical answers True exactly as a monoradical does;
        the arena stores one bit per atom and nothing about spin pairing, so the question this can
        answer honestly is "does this molecule have an unpaired electron somewhere".  `radical_of`
        per atom is the way to ask which, and how many.
        """
        self._require_clean()
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t i
        for i in range(self._structure.header.atom_count):
            if at_radical(&atoms[i]):
                return True
        return False

    @property
    def element_counts(self):
        """Return a dict mapping atomic number to atom count, for elements present in the molecule.

        The keys ascend by atomic number, so an R marker's key 0 comes first; `brutto` answers the
        same question in its own order, keyed by symbol.
        """
        self._require_clean()
        cdef uint32_t e, begin, end
        cdef dict out = {}
        begin = element_bucket_begin(self._structure, 0)
        end = element_bucket_end(self._structure, 0)
        if end > begin:
            out[0] = end - begin
        for e in range(1, 119):
            begin = element_bucket_begin(self._structure, e)
            end = element_bucket_end(self._structure, e)
            if end > begin:
                out[e] = end - begin
        return out

    @property
    def brutto(self):
        """The molecular formula as a dict from element symbol to count.

        HYDROGENS ARE FOLDED: 'H' counts the hydrogen ATOMS plus every IMPLICIT hydrogen, so ethanol
        is `{'C': 2, 'H': 6, 'O': 1}` however its six hydrogens are spelt.  That is the difference
        from `element_counts`, which is keyed by atomic number and counts atoms only -- the reason
        neither is an alias of the other.

        ISOTOPES AND CHARGES ARE NOT IN A FORMULA.  Heavy water is `{'H': 2, 'O': 1}` and the
        ammonium ion is `{'H': 4, 'N': 1}`; a formula counts elements, and the isotope and the charge
        belong to the atom.

        THE ORDER IS NOT HILL'S: C, H, O, N, B first -- in that sequence, with
        whichever of them the molecule lacks simply absent -- and then everything else by ascending
        atomic number.  Hill would put F before O; this puts O before F, because O is in the lead.
        The order is part of the answer, because `brutto_formula` is a join over it.

        AN ATOM WHOSE IMPLICIT COUNT IS UNKNOWN CONTRIBUTES NO HYDROGEN, the same silence
        `float(mol)` keeps, and for the same reason: a formula is not the place to raise, and there is
        no better number. `unknown_h_count` is non-zero on exactly the molecules where this formula is
        a lower bound on hydrogen rather than a formula.

        R MARKERS ARE COUNTED LAST under the key `'R'`.  Every R is one marker regardless of its index
        (R1, R7, …) — a formula counts atom kinds, and a fragment attachment point is one kind.  The
        `'R'` key sorts after every element, so `brutto_formula` ends with the marker count.
        """
        self._require_clean()
        cdef atom_t *atoms = self._structure.atoms()
        cdef uint32_t i, e, begin, end
        cdef uint32_t hydrogens = 0
        for i in range(self._structure.header.atom_count):
            if not at_implicit_h_unknown(&atoms[i]):
                hydrogens += <uint32_t> at_implicit_h(&atoms[i])
        # a dict preserves insertion order, so the lead is written first and the tail skips it
        cdef dict out = {}
        cdef tuple lead = (6, 1, 8, 7, 5)
        for e in lead:
            begin = element_bucket_begin(self._structure, e)
            end = element_bucket_end(self._structure, e)
            if e == 1:
                if end - begin + hydrogens:
                    out['H'] = end - begin + hydrogens
            elif end > begin:
                out[SYMBOLS[e - 1]] = end - begin
        for e in range(1, 119):
            if e in lead:
                continue
            begin = element_bucket_begin(self._structure, e)
            end = element_bucket_end(self._structure, e)
            if end > begin:
                out[SYMBOLS[e - 1]] = end - begin
        begin = element_bucket_begin(self._structure, 0)
        end = element_bucket_end(self._structure, 0)
        if end > begin:
            out['R'] = end - begin
        return out

    @property
    def brutto_formula(self):
        """`brutto` as a string, a count of 1 written as nothing: aspirin is `'C9H8O4'`.

        In `brutto`'s order and therefore not in Hill's -- triflic acid comes out `'CHO3F3S'`, not
        `'CHF3O3S'`.  An empty molecule gives an empty string.
        """
        cdef list out = []
        cdef object symbol, count
        for symbol, count in self.brutto.items():
            out.append(symbol if count == 1 else '%s%d' % (symbol, count))
        return ''.join(out)

    @property
    def brutto_formula_html(self):
        """`brutto_formula` with each count above one in a `<sub>`: aspirin is
        `'C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>'`.

        In `brutto`'s order, like `brutto_formula`, and a count of 1 is written as nothing rather than
        as `<sub>1</sub>`.  Nothing is escaped because an element symbol and a decimal count are the
        whole alphabet here.
        """
        cdef list out = []
        cdef object symbol, count
        for symbol, count in self.brutto.items():
            out.append(symbol if count == 1 else '%s<sub>%d</sub>' % (symbol, count))
        return ''.join(out)

    # --- alternative spellings -------------------------------------------------------------------
    #
    # Four names for questions this class also answers another way.  They are SPELLINGS AND NOT
    # DEPRECATIONS: `atoms_count` and `len(mol)` are the same question, so neither is the one true
    # form and neither warns -- a library must not print into an application's output for using an
    # API that works.
    #
    # ONLY AN EXACT EQUIVALENCE BELONGS HERE, checked against chython 2 over compounds covering
    # charge, radicals, isotopes, implicit hydrogens and multiple components.  A name whose meaning
    # differs from anything this class computes is not a spelling and gets no entry, because one that
    # is nearly right ports a consumer silently and wrongly.  `element_counts` versus `brutto` is the
    # example: one is keyed by atomic number and counts atoms, the other is keyed by symbol, folds
    # implicit hydrogens into 'H' and imposes a C/H/O/N/B ordering.  Two questions, so `brutto` is a
    # method of its own above, and so is `is_radical` -- an any-atom fold that `radical_of` answers
    # per atom.  `aromatic_rings` is a filter over the ring set, not a rename of `rings`.
    #
    # NOR IS A NAME THIS CLASS ALREADY ANSWERS TO: a second spelling of a live name SHADOWS it, so
    # the live name starts answering this body.  `atoms_numbers` is deliberately absent -- the list
    # is `atom_numbers`, and the plural spelling raises AttributeError, which is the intended signal.
    #
    # `has_atom` and `has_bond` were shims rather than spellings and are gone: `has_bond` differed
    # in the exception it raised, so a caller had to be edited either way.  `n in mol` and
    # `order_of(n, m) is not None` are the spellings.

    @property
    def atoms_count(self):
        """How many atoms.  `len(mol)` and `atom_count` are the other spellings."""
        return self.atom_count

    @property
    def bonds_count(self):
        """How many bonds.  `bond_count` is the other spelling."""
        return self.bond_count

    @property
    def molecular_charge(self):
        """Net charge.  `int(mol)` is the other spelling."""
        return self.__int__()

    @property
    def molecular_mass(self):
        """Average molecular mass in daltons.  `float(mol)` is the other spelling."""
        return self.__float__()

    # ------------------------------------------------------------------ F3 descriptors
    # Each body lives in `chython.chemistry` and is registered via `_set_featurizer_fns`.  The
    # core owns the ten names (they are slots on a `cdef class`, which cannot be extended from
    # outside), and every later task in the F3 phase fills in one body without touching the
    # extension.  Plain `@property` throughout: `functools.cached_property` needs a `__dict__`
    # that a `cdef class` does not have, and a descriptor on a mutable container must not
    # outlive an edit session anyway.

    @property
    def rotatable_bonds_count(self):
        """Number of rotatable bonds by chython's own definition (`tables/rotatable.tsv`).

        Not Lipinski's, not Veber's, not RDKit's strict or non-strict rule.  chython 2 counted
        MAPPINGS of a symmetric query and so returned twice this number; there is no flag to get
        the old answer back.
        """
        return _featurizer_fn('rotatable_bonds_count')(self)

    @property
    def hydrogen_bond_donors_count(self):
        """Hydrogen bond donor count over `tables/hbond.tsv`."""
        return _featurizer_fn('hydrogen_bond_donors_count')(self)

    @property
    def hydrogen_bond_acceptors_count(self):
        """Hydrogen bond acceptor count over `tables/hbond.tsv`."""
        return _featurizer_fn('hydrogen_bond_acceptors_count')(self)

    @property
    def tpsa(self):
        """Topological polar surface area in A^2.

        Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714.  N and O only, as published; the
        paper's S and P contributions are in `tables/tpsa.tsv` and are off by default.
        """
        return _featurizer_fn('tpsa')(self)

    @property
    def crippen_logp(self):
        """Wildman-Crippen atomic-contribution logP.

        Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868.
        """
        return _featurizer_fn('crippen_logp')(self)

    @property
    def crippen_mr(self):
        """Wildman-Crippen molar refractivity.

        Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868.
        """
        return _featurizer_fn('crippen_mr')(self)

    @property
    def qed(self):
        """Quantitative estimate of drug-likeness, QED_w,mo in the paper's notation.

        Bickerton, Paolini, Besnard, Muresan, Hopkins, Nat. Chem. 2012, 4, 90.  The mean-weight
        variant, which is the paper's recommended default; `chython.chemistry.qed(m, weights=...)`
        reaches the other two published weight sets.
        """
        return _featurizer_fn('qed')(self)

    def maccs_keys(self):
        """The 166 published MACCS structural keys as `ndarray(167) uint8`.

        Durant, Leland, Henry, Nourse, J. Chem. Inf. Comput. Sci. 2002, 42, 1273.  ONE-BASED:
        index 0 is permanently zero so that `keys[n]` is published key `n`.  Widely used
        implementations knowingly differ from the published keys; these are the published keys.
        """
        return _featurizer_fn('maccs_keys')(self)

    def maccs_bit_set(self):
        """The set of published MACCS key numbers this molecule sets.  1..166, never 0."""
        return _featurizer_fn('maccs_bit_set')(self)

    def pharmacophore_invariants(self):
        """Per-atom pharmacophore feature type as `ndarray(n) uint32`, in `atom_numbers` order.

        Six 2D types after Kutlushina, Khakimova, Madzhidov, Polishchuk, Molecules 2018, 23, 3094.
        Suitable as `invariants=` to any `morgan_*` or `linear_*` method.
        """
        return _featurizer_fn('pharmacophore_invariants')(self)


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    JOURNAL_OPS = {'add_atom': OP_ADD_ATOM, 'delete_atom': OP_DELETE_ATOM,
                   'add_bond': OP_ADD_BOND, 'delete_bond': OP_DELETE_BOND,
                   'set_order': OP_SET_ORDER, 'set_charge': OP_SET_CHARGE,
                   'set_isotope': OP_SET_ISOTOPE, 'set_radical': OP_SET_RADICAL,
                   'set_map_number': OP_SET_MAP_NUMBER, 'set_hydrogens': OP_SET_HYDROGENS,
                   'set_stereo': OP_SET_STEREO, 'set_xy': OP_SET_XY,
                   'set_xyz': OP_SET_XYZ,
                   'set_wedge': OP_SET_WEDGE, 'set_stereo_group': OP_SET_STEREO_GROUP,
                   'set_element': OP_SET_ELEMENT,
                   'set_r_index': OP_SET_R_INDEX, 'add_conformer': OP_ADD_CONFORMER,
                   'drop_conformer': OP_DROP_CONFORMER}
    WEDGE_NONE = 0
    WEDGE_UP = 1
    WEDGE_DOWN = 2
    WEDGE_EITHER = 3
    STEREO_UNSPECIFIED = 0
    STEREO_ABS = 1
    STEREO_OR = 2
    STEREO_AND = 3
    # The implicit-hydrogen sentinel, on the surface so a parser can WRITE it.  Readers get None
    # (`implicit_h_of`, `Atom.implicit_h`, `Atom.total_h`); only a writer needs the number, and it
    # is published from the C DEF so the two cannot drift apart -- this is the same 15 the nibble
    # holds, not a second constant that happens to agree today.
    #
    # THROUGH `globals()` BECAUSE `H_UNKNOWN = H_UNKNOWN` DOES NOT COMPILE.  A `DEF` name is
    # substituted textually wherever it appears as a name, assignment targets included, so the
    # obvious line becomes `15 = 15`.  The alternative was to give the Python surface a second
    # spelling, which is what this constant had before and what the project ruled against: the C
    # name and the Python name must be the same word.  A string key is the one place a DEF name
    # survives unsubstituted.
    globals()['H_UNKNOWN'] = H_UNKNOWN
    # And the bound that goes with it, for the same reason and by the same route.  A reader
    # validating a count it parsed out of a file needs the number 14, and the CTfile reader had
    # already reached for it and written a literal `H_MAX = 15` -- which admitted the sentinel as a
    # count on three separate write paths.  Exported so no caller has to restate it: a bound
    # duplicated as a literal is a bound that drifts, and this one drifting turns a stated count
    # into "nobody knows".
    globals()['H_IMPLICIT_MAX'] = H_IMPLICIT_MAX
    # The R index's domain, out for the same reason as the hydrogen bound above: the readers and
    # writers validate an index they parsed out of a file, and a bound restated as a literal there
    # would be a bound that drifts.  Two decimal digits is a promise the CTfile symbol column and the
    # depiction label both rely on.
    globals()['R_INDEX_MAX'] = R_INDEX_MAX


def journal_record_size():
    return sizeof(journal_t)


def _unit_parity_raw_probe(MoleculeContainer mol not None, uint32_t n):
    """Read the raw `u.parity` byte of the stereo unit anchored at `n`.

    Used in tests to verify that `translate_stereo` does not write the unit record.
    `u.parity` must always be 0; the field is reserved -- the parity is a byte per atom in
    SEG_PARITY, keyed by anchor slot.  `_stereo_emit` is the sole writer of the unit record --
    this probe reads the field so a test can assert the invariant is upheld.
    """
    mol._require_clean()
    if n not in mol._index_of:
        raise KeyError(n)
    ensure_stereo_units(mol._structure)
    cdef uint32_t slot = <uint32_t> mol._index_of[n]
    cdef stereo_unit_t *u = stereo_unit_of(mol._structure, slot)
    if u is NULL:
        raise KeyError(n)
    return <int> u.parity


def _stereo_forge_truncation(MoleculeContainer mol not None):
    """FORGE the stereo table's truncation word on an already-built table: the SECOND-READ path.

    Molecules do reach truncation on their own, and the example is a CONNECTED one:
    cyclo[CH(CH3)CH2]6 -- 1,3,5,7,9,11-hexamethylcyclododecane -- with every hydrogen written as an
    explicit atom, 54 atoms in one component, ~13 ms, which is a real test of its own
    (`test_a_connected_record_can_still_exhaust_the_budget`).  DISCONNECTED COPIES ARE NOT THE
    EXAMPLE: the witness search is restricted to the anchor's own component, so five copies of
    1,3,5,7-tetramethylcyclooctane decide in under a millisecond and report `stereo_truncated` FALSE.
    Build that record expecting a truncation and you will not get one; that is the restriction
    working, not this probe testing an unreachable state.

    What a really-truncating record cannot test is the state a table is in on every read AFTER the one
    that built it, because there the word is all that is left of the search: perception is not re-run.
    So this probe writes the word directly into a built table and a test reads it back
    (`test_a_forged_truncation_word_survives_every_later_read`) -- which is the ONLY remaining coverage
    of `chiral_bonds()`, `is_chiral()` and `unit_of()` on a record whose truncation word is set, so
    deleting this probe deletes that.  It is a forged record and only that: it does not simulate the
    search, and it is not evidence about which molecules truncate.
    """
    mol._require_clean()
    ensure_stereo_units(mol._structure)
    (<uint32_t *> mol._structure.segment(SEG_STEREO_UNIT))[1] = 1


def _rebase_parity_probe(MoleculeContainer mol not None, uint32_t n, int kind, int parity,
                         tuple old_refs, int old_unnamed_mask):
    """`rebase_parity` against one molecule's live table, with the old frame supplied by hand.

    Two of that function's refusals are not reachable by editing a molecule -- a kind change under
    an unmoved anchor, and a bond frame whose pairs cross -- because an edit that could produce
    either also changes the directions, and the leftover budget refuses first.  A probe is
    how those branches get measured rather than assumed; it is a forged frame and only that.

    `old_unnamed_mask` has NO DEFAULT on purpose (ruling F69): `old_refs` spells an unnamed direction
    and an empty slot the same way, the mask is the only thing that separates them, and a default
    would let a call site pass a frame whose flavours are silently wrong and still read a number back.
    Every caller states it.
    """
    cdef uint32_t refs[4]
    cdef int i
    mol._require_clean()
    if n not in mol._index_of:
        raise KeyError(n)
    # Range-checked here rather than trusted: the probe is the only Python-visible door to this
    # arithmetic, and unchecked an out-of-range `parity` comes back out as itself.
    if kind < 0 or kind > SU_HELICAL:
        raise ValueError(f'kind must be a stereo unit kind in 0..{<int> SU_HELICAL}, got {kind}')
    if parity < 0 or parity > 2:
        raise ValueError(f'parity must be 0 (none), 1 (even) or 2 (odd), got {parity}')
    if old_unnamed_mask < 0 or old_unnamed_mask > SU_UNNAMED_MASK:
        raise ValueError(f'old_unnamed_mask must be a 4-bit slot mask, got {old_unnamed_mask}')
    ensure_stereo_units(mol._structure)
    for i in range(4):
        refs[i] = SU_NO_REF if old_refs[i] is None else <uint32_t> mol._index_of[old_refs[i]]
    return rebase_parity(mol._structure, <uint32_t> mol._index_of[n], <uint8_t> kind,
                         <uint8_t> parity, refs, <uint8_t> old_unnamed_mask)


def _from_bytes(data, meta=None):
    # `meta` defaults to None so a pickle written before it was carried still loads.
    cdef MoleculeContainer mol = MoleculeContainer.from_bytes(data)
    if meta:
        mol._meta = dict(meta)
    return mol


def _segment_span(MoleculeContainer mol not None, int seg):
    """(offset, length) of a live molecule's persistent segment, for testing only.

    Returns (0, 0) for a segment absent from the table (id past `seg_count`).
    The offset is into the bytes that `mol.to_bytes()` returns; a caller that serialises first
    and then edits that buffer must use this probe BEFORE the serialisation, or the offsets will
    not match.
    """
    cdef uint16_t sc = mol._structure.header.seg_count
    cdef uint32_t off, ln
    if seg >= sc:
        return (0, 0)
    off = mol._structure.header.segments[seg].offset
    ln = mol._structure.header.segments[seg].length
    return (<uint32_t> off, <uint32_t> ln)


def _parity_bytes(MoleculeContainer mol not None):
    """The parity segment verbatim, `b''` when absent.  Test-only: it reads storage rather than the
    answer `parity_of` gives, which is what lets a test compare the two."""
    cdef Structure s = mol._structure
    if structure_seg_len(s, SEG_PARITY) == 0:
        return b''
    return <bytes> (<char *> structure_parities(s))[:s.header.atom_count]


cdef class _EditScope:
    cdef MoleculeContainer _molecule

    def __cinit__(self, MoleculeContainer molecule not None):
        self._molecule = molecule

    def __enter__(self):
        self._molecule._scope_depth += 1
        return self._molecule

    @cython.warn.unused_arg(False)
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._molecule._scope_depth -= 1
        if self._molecule._scope_depth == 0:
            if exc_type is None:
                self._molecule._apply()
            else:
                self._molecule._discard()
        return False
