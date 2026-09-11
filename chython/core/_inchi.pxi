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
# Native libinchi binding — both directions.
#
# THREAD SAFETY.  libinchi's main entry points (GetStdINCHI, GetINCHI, GetStructFromINCHI,
# GetINCHIKeyFromINCHI) allocate all working state on the stack.  The only module-level
# mutable is `bInterrupted`, which is written only by a signal handler and is never written
# during a normal call.  The bundled macOS .dylib was tested with 16 concurrent Python threads
# (mixed forward and reverse, 16,000 calls total) and produced zero disagreements relative to
# the single-threaded reference — with and without a serialising lock.  No lock is used.
# If a future platform's libinchi build is found to have thread-local state or global caches,
# re-introduce the lock at that point with a comment citing the specific build artifact.
#
# MEMORY. Forward direction: caller allocs ICH_Atom[] and ICH_Stereo0D[] with PyMem_Malloc and
# frees them in finally; library allocs ICH_Output strings, freed by FreeStdINCHI/FreeINCHI in
# finally.  Reverse direction: caller allocs nothing (input is a Python bytes object); library
# allocs ICH_OutputStruct.atom and .stereo0D, freed by FreeStructFromINCHI in finally.
#
# STEREO. Ruling F26 applies: parity_of() returns a frame-relative byte; always go through
# translate_stereo(anchor, refs) with refs in the order handed to libinchi.
# See §6 of the spec (docs/superpowers/specs/2026-09-02-inchi-native-design.md) for the
# full convention mapping between inchi_Stereo0D parities and the core's 1/2 encoding.
#
# KEKULÉ. InChI always outputs Kekulé bond orders, so the reverse path never encounters
# order 4 in practice.  Forward path: if the arena stores Kekulé bonds (HE_AROMATIC not
# set, no order-4 half-edges) the atoms are filled directly.  If aromatic bonds are stored
# (HE_AROMATIC flag or half-edge order 4 — live once the arom module lands) the forward
# path calls the registered kekuliser (_ich_kekule_fn) on a copy, then fills that copy.
# If no kekuliser is registered, molecule_to_inchi raises ValueError.  The caller's
# molecule is NEVER mutated — input fidelity is an invariant of molecule_to_inchi.
# Register a kekuliser: _ich_set_kekule_fn(fn).  The arom module calls this at init.


cdef extern from *:
    """
    #include <stdlib.h>
    #include <string.h>

    #ifdef _WIN32
    #  include <windows.h>
    static void *ich_open_lib(const char *path) { return (void *)LoadLibraryA(path); }
    static void *ich_get_sym(void *h, const char *n) { return (void *)GetProcAddress((HMODULE)h, n); }
    #else
    #  include <dlfcn.h>
    static void *ich_open_lib(const char *path) { return dlopen(path, RTLD_LAZY | RTLD_LOCAL); }
    static void *ich_get_sym(void *h, const char *n) { return dlsym(h, n); }
    #endif

    /* ---- subset of inchi_api.h needed here ----------------------------------- */
    typedef signed short   ICH_AT_NUM;
    typedef signed char    ICH_S_CHAR;

    typedef struct {
        double       x, y, z;
        ICH_AT_NUM   neighbor[20];
        ICH_S_CHAR   bond_type[20];
        ICH_S_CHAR   bond_stereo[20];
        char         elname[6];
        ICH_AT_NUM   num_bonds;
        ICH_S_CHAR   num_iso_H[4];
        ICH_AT_NUM   isotopic_mass;
        ICH_S_CHAR   radical;
        ICH_S_CHAR   charge;
    } ICH_Atom;

    typedef struct {
        ICH_AT_NUM   neighbor[4];
        ICH_AT_NUM   central_atom;
        ICH_S_CHAR   type;
        ICH_S_CHAR   parity;
    } ICH_Stereo0D;

    typedef struct {
        ICH_Atom     *atom;
        ICH_Stereo0D *stereo0D;
        char         *szOptions;
        ICH_AT_NUM    num_atoms;
        ICH_AT_NUM    num_stereo0D;
    } ICH_Input;

    typedef struct {
        char *szInChI;
        char *szAuxInfo;
        char *szMessage;
        char *szLog;
    } ICH_Output;

    typedef struct {
        char *szInChI;
        char *szOptions;
    } ICH_InputINCHI;

    typedef struct {
        ICH_Atom       *atom;
        ICH_Stereo0D   *stereo0D;
        ICH_AT_NUM      num_atoms;
        ICH_AT_NUM      num_stereo0D;
        char           *szMessage;
        char           *szLog;
        unsigned long   WarningFlags[2][2];
    } ICH_OutputStruct;

    /* parity and stereo type constants */
    #define ICH_PARITY_ODD     1
    #define ICH_PARITY_EVEN    2
    #define ICH_PARITY_UNKNOWN 3
    #define ICH_STEREO_DOUBLEBOND  1
    #define ICH_STEREO_TETRAHEDRAL 2
    #define ICH_STEREO_ALLENE      3
    #define ICH_NO_ATOM           -1
    #define ICH_ISOTOPIC_SHIFT_FLAG 10000

    /* Rounded standard atomic masses used by InChI for isotope delta encoding.
       index 0 is a dummy; valid range is [1..118].
       Source: chython.periodictable Element.mdl_isotope for each atomic number. */
    static const int ICH_STD_MASS[119] = {
        0,   /* 0 dummy */
        1,   /* H  */  4,   /* He */  7,   /* Li */  9,   /* Be */  11,  /* B  */
        12,  /* C  */  14,  /* N  */  16,  /* O  */  19,  /* F  */  20,  /* Ne */
        23,  /* Na */  24,  /* Mg */  27,  /* Al */  28,  /* Si */  31,  /* P  */
        32,  /* S  */  35,  /* Cl */  40,  /* Ar */  39,  /* K  */  40,  /* Ca */
        45,  /* Sc */  48,  /* Ti */  51,  /* V  */  52,  /* Cr */  55,  /* Mn */
        56,  /* Fe */  59,  /* Co */  59,  /* Ni */  64,  /* Cu */  65,  /* Zn */
        70,  /* Ga */  73,  /* Ge */  75,  /* As */  79,  /* Se */  80,  /* Br */
        84,  /* Kr */  85,  /* Rb */  88,  /* Sr */  89,  /* Y  */  91,  /* Zr */
        93,  /* Nb */  96,  /* Mo */  98,  /* Tc */  101, /* Ru */  103, /* Rh */
        106, /* Pd */  108, /* Ag */  112, /* Cd */  115, /* In */  119, /* Sn */
        122, /* Sb */  128, /* Te */  127, /* I  */  131, /* Xe */  133, /* Cs */
        137, /* Ba */  139, /* La */  140, /* Ce */  141, /* Pr */  144, /* Nd */
        145, /* Pm */  150, /* Sm */  152, /* Eu */  157, /* Gd */  159, /* Tb */
        163, /* Dy */  165, /* Ho */  167, /* Er */  169, /* Tm */  173, /* Yb */
        175, /* Lu */  178, /* Hf */  181, /* Ta */  184, /* W  */  186, /* Re */
        190, /* Os */  192, /* Ir */  195, /* Pt */  197, /* Au */  201, /* Hg */
        204, /* Tl */  207, /* Pb */  209, /* Bi */  209, /* Po */  210, /* At */
        222, /* Rn */  223, /* Fr */  226, /* Ra */  227, /* Ac */  232, /* Th */
        231, /* Pa */  238, /* U  */  237, /* Np */  244, /* Pu */  243, /* Am */
        247, /* Cm */  247, /* Bk */  251, /* Cf */  252, /* Es */  257, /* Fm */
        258, /* Md */  259, /* No */  260, /* Lr */  261, /* Rf */  270, /* Db */
        269, /* Sg */  270, /* Bh */  270, /* Hs */  278, /* Mt */  281, /* Ds */
        281, /* Rg */  285, /* Cn */  278, /* Nh */  289, /* Fl */  289, /* Mc */
        293, /* Lv */  297, /* Ts */  294  /* Og */
    };

    /* function pointer types */
    typedef int  (*ich_fn_GetStdINCHI_t)(ICH_Input *, ICH_Output *);
    typedef void (*ich_fn_FreeStdINCHI_t)(ICH_Output *);
    typedef int  (*ich_fn_GetINCHI_t)(ICH_Input *, ICH_Output *);
    typedef void (*ich_fn_FreeINCHI_t)(ICH_Output *);
    typedef int  (*ich_fn_GetStructFromINCHI_t)(ICH_InputINCHI *, ICH_OutputStruct *);
    typedef void (*ich_fn_FreeStructFromINCHI_t)(ICH_OutputStruct *);
    typedef int  (*ich_fn_GetINCHIKeyFromINCHI_t)(const char *, int, int, char *, char *, char *);
    """
    # `signed char` AND NOT `char` WHEREVER THE HEADER SAYS `S_CHAR`, which is `signed char`: plain
    # `char` is unsigned on Linux ARM, and this block is what Cython's coercions are generated from.
    # The C struct above is `ICH_S_CHAR` either way, so the LAYOUT never differed -- what differed is
    # the conversion `int(a_ptr.charge)` compiles to, which on aarch64 read a chloride's -1 as 255.
    # `elname` and every `sz*` stay plain `char`, being what the header says and being strings.
    ctypedef struct ICH_Atom:
        double x
        double y
        double z
        short neighbor[20]
        signed char bond_type[20]
        signed char bond_stereo[20]
        char elname[6]
        short num_bonds
        signed char num_iso_H[4]
        short isotopic_mass
        signed char radical
        signed char charge

    ctypedef struct ICH_Stereo0D:
        short neighbor[4]
        short central_atom
        signed char type
        signed char parity

    ctypedef struct ICH_Input:
        ICH_Atom *atom
        ICH_Stereo0D *stereo0D
        char *szOptions
        short num_atoms
        short num_stereo0D

    ctypedef struct ICH_Output:
        char *szInChI
        char *szAuxInfo
        char *szMessage
        char *szLog

    ctypedef struct ICH_InputINCHI:
        char *szInChI
        char *szOptions

    ctypedef struct ICH_OutputStruct:
        ICH_Atom *atom
        ICH_Stereo0D *stereo0D
        short num_atoms
        short num_stereo0D
        char *szMessage
        char *szLog

    int ICH_PARITY_ODD
    int ICH_PARITY_EVEN
    int ICH_PARITY_UNKNOWN
    int ICH_STEREO_DOUBLEBOND
    int ICH_STEREO_TETRAHEDRAL
    int ICH_STEREO_ALLENE
    int ICH_NO_ATOM
    int ICH_ISOTOPIC_SHIFT_FLAG
    int ICH_STD_MASS[119]

    ctypedef int  (*ich_fn_GetStdINCHI_t)(ICH_Input *, ICH_Output *)
    ctypedef void (*ich_fn_FreeStdINCHI_t)(ICH_Output *)
    ctypedef int  (*ich_fn_GetINCHI_t)(ICH_Input *, ICH_Output *)
    ctypedef void (*ich_fn_FreeINCHI_t)(ICH_Output *)
    ctypedef int  (*ich_fn_GetStructFromINCHI_t)(ICH_InputINCHI *, ICH_OutputStruct *)
    ctypedef void (*ich_fn_FreeStructFromINCHI_t)(ICH_OutputStruct *)
    ctypedef int  (*ich_fn_GetINCHIKeyFromINCHI_t)(const char *, int, int, char *, char *, char *)

    void *ich_open_lib(const char *path)
    void *ich_get_sym(void *h, const char *name)


# ------------------------------------------------------------------------------------------------
# BOND-KIND STEREO: THE TWO FRAMES, AND THE TWO SIGNS THAT RELATE THEM.
#
# A bond-kind unit (SU_CIS_TRANS, SU_ALLENE) stores its parity over FOUR SUBSTITUENT directions in
# two ruling-F26 pairs, `refs = (a1, a2 | b1, b2)`: `a1, a2` on one terminal, `b1, b2` on the other,
# and slot 0 of each pair is the NAMED atom (F26/F47, so `refs[0]` and `refs[2]` are never
# SU_NO_REF).  `translate_stereo` reads exactly that frame and enforces it: `order[0:2]` must map
# entirely to one stored pair and `order[2:4]` to the other, or it raises.
#
# InChI's `neighbor[4]` is a DIFFERENT frame, and for the allene a genuinely MIXED one:
# `{X, A, B, Y}` where `A` and `B` are the two CHAIN atoms and `X`, `Y` are ONE substituent each,
# on `A` and on `B` respectively (inchi_api.h ~lines 200-300).  So the record names two chain atoms
# that are not directions at all, and only half the substituents.  Two consequences:
#
#   * the code that FILLS `rec.neighbor` is right to put the chain atoms in slots 1 and 2, and must
#     keep doing so.  It is only the PARITY that has to be computed in the unit's own frame, by
#     handing `translate_stereo` the pair-grouped `(X, other-of-X's-pair, Y, other-of-Y's-pair)`.
#   * InChI names only ONE substituent per terminal, so if it names `a2` where we would name `a1`,
#     that is one WITHIN-PAIR transposition and the parity flips.  `translate_stereo` accounts for
#     that automatically once the order is pair-grouped -- which is why both directions below locate
#     `X` and `Y` inside `refs` rather than assuming `X == refs[0]`.
#
# WHAT REMAINS after the frames agree is one fixed sign per kind, and the two have DIFFERENT
# provenance.  Both were measured on 2026-09-02 against libinchi built from the bundled
# INCHI submodule (11a8798), which PERCEIVES both kinds from coordinates and is therefore an
# absolute reference -- the only one available (RDKit clears the allene tag on SanitizeMol,
# OpenBabel refuses the syntax, Indigo drops the layer on its own InChI export, and chython 2's
# axial frame is unstated, so asking V2 is circular).
#
# ALLENE -- DERIVED, not chosen.  The core's tetrahedral convention is documented in the SMILES
# writer and was measured there against RDKit from a hand-built conformer: parity 1 (even) is a
# POSITIVE signed volume `(p1-p0).((p2-p0)x(p3-p0))` over the ref order, i.e. "the remaining three
# clockwise seen from the first".  InChI's rule is the SAME sentence in ITS frame ("if A, B, Y are
# clockwise when seen from X then parity is 'e'"), and libinchi confirms it: for
# 1,3-dibromo-1,3-difluoroallene the two mirror conformers perceive as `/t1-/m1/s1` and
# `/t1-/m0/s1`, and a zero-coordinate 0D record of EVEN / ODD reproduces those two strings exactly.
# Applying one signed-volume rule to both frames is then pure geometry, and the answer is a FLIP:
# over 3,000 random non-degenerate allene geometries, `sign(V over (a1,a2,b1,b2))` was ALWAYS the
# opposite of `sign(V over (a1,A,B,b1))`.  So chython even (1) == InChI ODD.  This is not a free
# choice; changing it contradicts the measurement.
#
# CIS/TRANS -- CHOSEN, and this file is the first consumer to state it.  The signed-volume rule is
# DEGENERATE here: the four substituents are coplanar, V is zero, and it does not merely vanish at
# the planar limit but keeps the SAME sign on both sides of it (measured across the twisted-cumulene
# continuum: V stays negative at 5 deg and at 175 deg alike).  So the allene result does NOT
# transfer and nothing else in the tree states a meaning -- the SMILES writer does not yet write
# cis/trans at all (`test_cis_trans_configuration_is_not_yet_written`).  InChI's own rule is
# unambiguous and was confirmed against its 2D perception: EVEN == `X` and `Y` on OPPOSITE sides
# (but-2-ene, `/b4-3+` for the trans conformer and for a 0D EVEN record; `/b4-3-` for both cis
# spellings).  We adopt it unchanged, so the core's parity gains a stated geometric meaning:
#
#     chython parity 1 (even) == refs[0] and refs[2] are TRANS (opposite sides).
#
# If a later consumer (the SMILES writer's `/`-`\` work) needs the other polarity, THIS is the line
# to change, and the round trip will not notice -- only the absolute tests will.
#
# Both constants have the SAME meaning -- 1 flips, 0 does not -- and each is read in exactly two
# places, the export flip point and the import one.  They are the only sign decisions in this file.
# ------------------------------------------------------------------------------------------------
DEF ICH_ALLENE_FLIP    = 1     # 1: chython even -> InChI ODD.   DERIVED from libinchi's 3D sign.
DEF ICH_CIS_TRANS_FLIP = 0     # 0: chython even -> InChI EVEN, i.e. even == trans.  CHOSEN.


# ---- module-level state ----------------------------------------------------- #

cdef void *ich_lib_handle = NULL
cdef ich_fn_GetStdINCHI_t ich_fn_GetStdINCHI = NULL
cdef ich_fn_FreeStdINCHI_t ich_fn_FreeStdINCHI = NULL
cdef ich_fn_GetINCHI_t ich_fn_GetINCHI = NULL
cdef ich_fn_FreeINCHI_t ich_fn_FreeINCHI = NULL
cdef ich_fn_GetStructFromINCHI_t ich_fn_GetStructFromINCHI = NULL
cdef ich_fn_FreeStructFromINCHI_t ich_fn_FreeStructFromINCHI = NULL
cdef ich_fn_GetINCHIKeyFromINCHI_t ich_fn_GetINCHIKeyFromINCHI = NULL

# No serialising lock (see thread-safety note in the module header).
# A plain None sentinel is kept so the module's namespace stays auditable and can
# be patched if a specific build turns out to need serialisation.
cdef object _ich_lock
_ich_lock = None

# Kekuliser hook: set by the arom module when it initialises.
# molecule_to_inchi calls _ich_kekule_fn(mol) if aromatic bonds are detected.
# The function must return a Kekule copy without mutating its input.
cdef object _ich_kekule_fn
_ich_kekule_fn = None


def _ich_set_kekule_fn(fn):
    """Register the kekuliser used by molecule_to_inchi for aromatic-bond molecules.

    ``fn`` must have signature ``fn(mol: MoleculeContainer) -> MoleculeContainer``.
    It must return a new Kekule copy and must NOT mutate the input.

    Called by the arom module at its own initialisation time.
    """
    global _ich_kekule_fn
    _ich_kekule_fn = fn

# ---- element symbol → atomic number lookup ---------------------------------- #
# Populated at module init time from the SYMBOLS tuple defined in _elements.pxi.

cdef dict _ich_sym_to_z = {}
cdef uint32_t _ich_i
for _ich_i in range(len(SYMBOLS)):
    _ich_sym_to_z[SYMBOLS[_ich_i]] = _ich_i + 1
# _ich_i is a C uint32_t; it cannot be deleted. Its value is len(SYMBOLS) after the loop.


# ---- library load ----------------------------------------------------------- #

def ich_load_library(path):
    """Load libinchi from `path` (str or bytes) and resolve all required symbols.

    Returns True on success, False on any failure (library not found, symbol missing).
    Called once at module import time from chython/core/__init__.py.
    """
    cdef void *h
    cdef void *sym
    cdef bytes path_b

    global ich_lib_handle
    global ich_fn_GetStdINCHI, ich_fn_FreeStdINCHI
    global ich_fn_GetINCHI, ich_fn_FreeINCHI
    global ich_fn_GetStructFromINCHI, ich_fn_FreeStructFromINCHI
    global ich_fn_GetINCHIKeyFromINCHI

    if isinstance(path, str):
        path_b = path.encode()
    else:
        path_b = path

    h = ich_open_lib(<const char *> path_b)
    if h is NULL:
        return False

    sym = ich_get_sym(h, b'GetStdINCHI')
    if sym is NULL:
        return False
    ich_fn_GetStdINCHI = <ich_fn_GetStdINCHI_t> sym

    sym = ich_get_sym(h, b'FreeStdINCHI')
    if sym is NULL:
        return False
    ich_fn_FreeStdINCHI = <ich_fn_FreeStdINCHI_t> sym

    sym = ich_get_sym(h, b'GetINCHI')
    if sym is NULL:
        return False
    ich_fn_GetINCHI = <ich_fn_GetINCHI_t> sym

    sym = ich_get_sym(h, b'FreeINCHI')
    if sym is NULL:
        return False
    ich_fn_FreeINCHI = <ich_fn_FreeINCHI_t> sym

    sym = ich_get_sym(h, b'GetStructFromINCHI')
    if sym is NULL:
        return False
    ich_fn_GetStructFromINCHI = <ich_fn_GetStructFromINCHI_t> sym

    sym = ich_get_sym(h, b'FreeStructFromINCHI')
    if sym is NULL:
        return False
    ich_fn_FreeStructFromINCHI = <ich_fn_FreeStructFromINCHI_t> sym

    sym = ich_get_sym(h, b'GetINCHIKeyFromINCHI')
    if sym is NULL:
        return False
    ich_fn_GetINCHIKeyFromINCHI = <ich_fn_GetINCHIKeyFromINCHI_t> sym

    ich_lib_handle = h
    return True


def inchi_library_loaded():
    """Return True if libinchi has been loaded successfully."""
    return ich_lib_handle is not NULL


cdef MoleculeContainer _ich_ensure_kekule(MoleculeContainer mol):
    """Return mol unchanged if all bonds are Kekulé.

    If any half-edge has the HE_AROMATIC flag set or carries order 4, call the
    registered kekuliser (_ich_kekule_fn) on a copy and return that copy.  Raises
    ValueError if aromatic bonds are present but no kekuliser is registered.

    The caller's molecule is never mutated.
    """
    cdef halfedge_t *edges = csr_edges(mol._structure)
    cdef uint32_t *ptr = csr_ptr(mol._structure)
    cdef uint32_t n_atoms = mol._structure.header.atom_count
    cdef uint32_t i, k
    cdef bint has_arom = False
    cdef MoleculeContainer result

    for i in range(n_atoms):
        for k in range(ptr[i], ptr[i + 1]):
            if (edges + k).flags & HE_AROMATIC or (edges + k).order == 4:
                has_arom = True
                break
        if has_arom:
            break

    if has_arom:
        if _ich_kekule_fn is None:
            raise ValueError(
                'molecule contains aromatic bonds (HE_AROMATIC flag or order 4); '
                'the InChI path requires Kekule form. '
                'No kekuliser is registered — ensure the SMILES/arom module is '
                'compiled into this build, or call _ich_set_kekule_fn(fn) manually.'
            )
        result = _ich_kekule_fn(mol)
    else:
        result = mol
    return result


# ---- forward direction: MoleculeContainer → InChI string ------------------- #

cdef inline void _refuse_r_marker(MoleculeContainer mol) except *:
    """A molecule carrying an R has no InChI, whichever shared object was staged.

    Stated here rather than at each entry point so the two public functions cannot drift, and called
    before the `ich_lib_handle` check so the answer does not depend on the build.
    """
    if element_bucket_end(mol._structure, 0) > element_bucket_begin(mol._structure, 0):
        raise ValueError('this molecule carries an R marker, which InChI has no representation for. '
                         'Strip the markers, or use the canonical SMILES as the identity key.')


def _ich_platform_options(str options) -> str:
    """`options` with each flag's prefix rewritten to the one this platform's libinchi reads.

    `inchi_api.h` documents szOptions as "each is preceded by '/' or '-' depending on OS and compiler":
    `mode.h` sets `INCHI_OPTION_PREFX` to `/` under `_WIN32` and `-` otherwise, and `ichiparm.c` tests
    only that one character -- a token carrying the other prefix falls through to the input-path branch,
    so it is silently taken as a file name and `options='-SNon'` on Windows returns an InChI that still
    has its `/t` layer.  (`INCHI_ALT_OPT_PREFIX` is defined beside it and referenced nowhere.)  Both
    spellings are accepted here and rewritten, because a caller writes one string and the docstring,
    `docs/io.rst` and `test_container_methods.py` all write `-SNon`.
    """
    cdef str prefix

    from sys import platform

    prefix = '/' if platform == 'win32' else '-'
    return ' '.join(prefix + token[1:] if token[0] in '-/' else token for token in options.split())


def molecule_to_inchi(MoleculeContainer mol not None, *,
                      bint standard=True, str options=None) -> str:
    """Generate an InChI string from a MoleculeContainer.

    Parameters
    ----------
    mol:      the molecule.
    standard: True (default) calls GetStdINCHI; False calls GetINCHI.
    options:  option string, e.g. ``'-SNon'`` (no stereo).  None means no options.  Either prefix
              works -- `_ich_platform_options` rewrites it to the one this platform's library reads.

    Returns the InChI string, e.g. ``'InChI=1S/...'``.

    Raises ImportError if libinchi was not loaded.
    Raises ValueError on InChI generation failure, or if the molecule carries an R marker.

    What is discarded: atom-atom map numbers, enhanced stereo groups, atropisomers,
    XY coordinates.  See the spec for the full list.
    """
    cdef ICH_Atom *atoms
    cdef ICH_Stereo0D *s0d
    cdef ICH_Input inp
    cdef ICH_Output out
    cdef uint32_t n_atoms
    cdef int n_stereo
    cdef int rc = 0
    cdef bytes opt_bytes
    cdef str result = ''
    cdef str msg = ''

    _refuse_r_marker(mol)
    if ich_lib_handle is NULL:
        raise ImportError('libinchi not loaded; cannot generate InChI')

    mol._require_clean()
    # Kekulise a copy if needed (input fidelity: never mutate the caller's molecule).
    mol = _ich_ensure_kekule(mol)
    n_atoms = mol._structure.header.atom_count

    # Allocate atom array.
    atoms = <ICH_Atom *> PyMem_Malloc(<size_t> n_atoms * sizeof(ICH_Atom))
    if atoms is NULL:
        raise MemoryError('InChI atom array allocation failed')

    # Stereo0D records: worst case one per stereo unit; n is a safe upper bound.
    s0d = <ICH_Stereo0D *> PyMem_Malloc((<size_t> n_atoms + 1) * sizeof(ICH_Stereo0D))
    if s0d is NULL:
        PyMem_Free(atoms)
        raise MemoryError('InChI stereo0D array allocation failed')

    memset(&out, 0, sizeof(ICH_Output))

    try:
        _ich_fill_atoms(mol, atoms, n_atoms)
        n_stereo = _ich_fill_stereo(mol, s0d)

        if options is not None:
            opt_bytes = _ich_platform_options(options).encode()
        else:
            opt_bytes = b''

        inp.atom = atoms
        inp.stereo0D = s0d if n_stereo > 0 else NULL
        inp.szOptions = <char *> opt_bytes if opt_bytes else NULL
        inp.num_atoms = <short> n_atoms
        inp.num_stereo0D = <short> n_stereo

        if standard:
            rc = ich_fn_GetStdINCHI(&inp, &out)
        else:
            rc = ich_fn_GetINCHI(&inp, &out)
        try:
            if rc > 1:
                msg = out.szMessage.decode() if out.szMessage is not NULL else 'unknown error'
                raise ValueError(f'InChI generation failed (rc={rc}): {msg}')
            if out.szInChI is NULL:
                raise ValueError('InChI generation returned NULL string')
            result = out.szInChI.decode()
        finally:
            if standard:
                ich_fn_FreeStdINCHI(&out)
            else:
                ich_fn_FreeINCHI(&out)

        return result
    finally:
        PyMem_Free(atoms)
        PyMem_Free(s0d)


cdef void _ich_fill_atoms(MoleculeContainer mol, ICH_Atom *atoms, uint32_t n_atoms) except *:
    """Fill ICH_Atom[] from the arena's atom_t[] and CSR edge list.

    XY coordinates are set to zero; all stereo information goes through the 0D path.
    """
    cdef atom_t *src = mol._structure.atoms()
    cdef uint32_t *ptr = csr_ptr(mol._structure)
    cdef halfedge_t *edges = csr_edges(mol._structure)
    cdef atom_t *a
    cdef halfedge_t *he
    cdef ICH_Atom *dst
    cdef uint32_t i, k, nb
    cdef bytes sym_b
    cdef const char *sym_ptr
    cdef Py_ssize_t sym_len

    for i in range(n_atoms):
        a = src + i
        dst = atoms + i
        memset(dst, 0, sizeof(ICH_Atom))

        # element symbol: SYMBOLS is 0-indexed, element is 1-indexed
        sym_b = SYMBOLS[a.element - 1].encode('ascii')
        sym_ptr = sym_b
        sym_len = len(sym_b)
        if sym_len > 5:
            sym_len = 5
        memcpy(dst.elname, sym_ptr, sym_len)
        dst.elname[sym_len] = 0

        dst.charge = a.charge
        # isotope: 0 → not isotopic; absolute mass number, same in InChI
        dst.isotopic_mass = <short> a.isotope

        # radical: InChI DOUBLET=2 for any radical; our one-bit flag covers that case
        dst.radical = 2 if at_radical(a) else 0

        # implicit H: num_iso_H[0] = non-isotopic; isotopic H ([1..3] = D/T) not in arena.
        #
        # AN UNKNOWN COUNT GOES OUT AS InChI'S OWN -1 ("auto"), NOT AS THE SENTINEL.  Casting
        # `at_implicit_h(a)` unconditionally hands libinchi the sentinel as a literal 15: alanine comes
        # back as `C3H37NO2` with an `h` layer reading `3,6H15` and a moved InChIKey.  -1 asks libinchi
        # to fill the count from its own valence rules, which is the honest translation of "the record
        # does not say" into a format that has a way to say it -- and it is the exact value the import
        # path at `impl_h < 0` reads back as H_UNKNOWN, so the two directions agree.
        #
        # Reachable only from an explicit `implicit_h=H_UNKNOWN` or an MDL record with an undeterminable
        # count: both string readers state every count, and a builder atom takes one from the derivation.
        dst.num_iso_H[0] = -1 if at_implicit_h_unknown(a) else <signed char> at_implicit_h(a)

        # connectivity: emit all half-edges for atom i
        nb = 0
        for k in range(ptr[i], ptr[i + 1]):
            he = edges + k
            dst.neighbor[nb] = <short> he.to
            # arena bond orders: 1=single, 2=double, 3=triple, 8=dative
            # InChI bond orders: 1=single, 2=double, 3=triple, 4=aromatic
            # dative (8) → single; order 4 (aromatic) should have been kekulised by
            # _ich_ensure_kekule before reaching here, so it must not appear.
            if he.order == 8:
                dst.bond_type[nb] = 1
            elif he.order == 4 or he.flags & HE_AROMATIC:
                raise AssertionError(
                    'aromatic half-edge (order %d, flags 0x%x) reached _ich_fill_atoms; '
                    '_ich_ensure_kekule should have kekulised it' % (he.order, he.flags)
                )
            else:
                dst.bond_type[nb] = <signed char> he.order
            dst.bond_stereo[nb] = 0  # use 0D stereo, not wedge-based 2D
            nb += 1
        dst.num_bonds = <short> nb


cdef bint _ich_is_bonded(MoleculeContainer mol, int a_idx, int b_idx) noexcept:
    """True when slots `a_idx` and `b_idx` share a half-edge."""
    cdef uint32_t *ptr = csr_ptr(mol._structure)
    cdef halfedge_t *edges = csr_edges(mol._structure)
    cdef uint32_t k

    for k in range(ptr[a_idx], ptr[a_idx + 1]):
        if <int> (edges + k).to == b_idx:
            return True
    return False


cdef tuple _ich_bond_order_from_refs(tuple refs, object X_ref, object Y_ref):
    """Pair-group `refs` as `(X, other-of-X's-pair, Y, other-of-Y's-pair)`.

    This is the ONLY shape `translate_stereo` accepts for a bond kind: `order[0:2]` from one stored
    pair, `order[2:4]` from the other.  `X_ref`/`Y_ref` need not be the pairs' slot-0 atoms -- when
    InChI names the other member of a pair, the partner lookup (`i ^ 1`) puts it in the right slot
    and `translate_stereo` charges the within-pair transposition to the parity, which is precisely
    the flip that "InChI names only one substituent per terminal" implies.

    An element of the returned tuple is `None` where `refs` holds an unnamed direction (an implicit
    hydrogen) or a pinned slot; `translate_stereo` tells those apart from the unit's own
    ``unnamed_mask`` and enforces Ruling F55 on the pinned case, so nothing here has to.

    Returns `None` when `X_ref` and `Y_ref` are not both named refs of two DIFFERENT stored pairs --
    in which case the record cannot be expressed in this frame and the caller must skip it.
    """
    cdef int ix = -1, iy = -1, i

    for i in range(4):
        if refs[i] is not None and refs[i] == X_ref:
            ix = i
            break
    for i in range(4):
        if refs[i] is not None and refs[i] == Y_ref:
            iy = i
            break
    if ix < 0 or iy < 0 or (ix // 2) == (iy // 2):
        return None
    return (refs[ix], refs[ix ^ 1], refs[iy], refs[iy ^ 1])


cdef inline int _ich_chython_to_inchi(int p, int flip) noexcept nogil:
    """Core parity (1 even, 2 odd) -> InChI parity.  THE FLIP POINT, export direction."""
    if flip:
        return ICH_PARITY_ODD if p == 1 else ICH_PARITY_EVEN
    return ICH_PARITY_EVEN if p == 1 else ICH_PARITY_ODD


cdef inline int _ich_inchi_to_chython(int parity_ich, int flip) noexcept nogil:
    """InChI parity -> core parity (1 even, 2 odd).  THE FLIP POINT, import direction.

    Exact inverse of `_ich_chython_to_inchi` for the same `flip`, so the round trip is insensitive
    to the constants' values and only the absolute tests can pin them.
    """
    if flip:
        return 2 if parity_ich == ICH_PARITY_EVEN else 1
    return 1 if parity_ich == ICH_PARITY_EVEN else 2


cdef int _ich_fill_stereo(MoleculeContainer mol, ICH_Stereo0D *s0d) except -1:
    """Fill ICH_Stereo0D[] from the molecule's stereo units.

    Returns the number of 0D records written.  Skips units whose parity is unset (0)
    and SU_ATROPISOMER (no InChI layer for that kind).
    """
    cdef list units = mol.stereo_units()
    cdef dict u
    cdef int kind, n_stereo = 0, n_unnamed, p, _bit
    cdef uint32_t anchor_n
    cdef object X_ref, Y_ref, ref
    cdef tuple refs, order
    cdef list named_refs, idx_order
    cdef int partner_idx, A_idx, B_idx
    cdef ICH_Stereo0D *rec
    cdef int unnamed_mask, parity

    for u in units:
        parity = u['parity']
        if parity == 0:
            continue
        kind = u['kind']
        anchor_n = u['anchor']
        refs = u['refs']
        unnamed_mask = u['unnamed_mask']

        if kind == SU_TETRA:
            n_unnamed = bin(unnamed_mask).count('1')

            rec = s0d + n_stereo
            rec.type = ICH_STEREO_TETRAHEDRAL
            rec.central_atom = <short> mol._index_of[anchor_n]

            if n_unnamed == 0:
                # Four named neighbours: F26 order == InChI WXYZ order.
                idx_order = []
                for ref in refs:
                    if ref is None:
                        idx_order.append(-1)
                    else:
                        idx_order.append(<int> mol._index_of[ref])
                rec.neighbor[0] = <short> idx_order[0]
                rec.neighbor[1] = <short> idx_order[1]
                rec.neighbor[2] = <short> idx_order[2]
                rec.neighbor[3] = <short> idx_order[3]
                p = mol.translate_stereo(anchor_n, refs)
                rec.parity = ICH_PARITY_EVEN if p == 1 else ICH_PARITY_ODD

            elif n_unnamed == 1:
                # One implicit H.  InChI uses center as neighbor[0] proxy.
                # F26 order: (n0, n1, n2, None); InChI order: (center, n0, n1, n2).
                # Cyclic rotation of 4 = odd permutation → parity flips.
                named_refs = []
                for _bit in range(4):
                    if not (unnamed_mask >> _bit & 1) and refs[_bit] is not None:
                        named_refs.append(refs[_bit])
                if len(named_refs) != 3:
                    continue  # unexpected layout
                p = mol.translate_stereo(anchor_n, refs)
                rec.neighbor[0] = <short> mol._index_of[anchor_n]
                rec.neighbor[1] = <short> mol._index_of[named_refs[0]]
                rec.neighbor[2] = <short> mol._index_of[named_refs[1]]
                rec.neighbor[3] = <short> mol._index_of[named_refs[2]]
                # Flip: F26 parity 1 (even) → InChI parity ODD (the cyclic rotation)
                rec.parity = ICH_PARITY_ODD if p == 1 else ICH_PARITY_EVEN

            else:
                # Two or more unnamed directions: InChI cannot represent this uniquely.
                continue

            n_stereo += 1

        elif kind == SU_CIS_TRANS:
            # InChI DOUBLEBOND: neighbor = {X, A, B, Y}; central_atom = NO_ATOM.
            # refs layout from stereo_units(): pair 0 = (A_side_ref0, A_side_ref1),
            #                                 pair 1 = (B_side_ref0, B_side_ref1).
            X_ref = refs[0] if not (unnamed_mask & 1) else None
            Y_ref = refs[2] if not (unnamed_mask & 4) else None

            if X_ref is None or Y_ref is None:
                continue

            partner_idx = _ich_find_cis_trans_partner(mol, <int> mol._index_of[anchor_n])
            if partner_idx < 0:
                continue

            # Parity in the UNIT's frame, not InChI's: the chain atoms in neighbor[1:3] are not
            # directions and must never reach translate_stereo.  refs[0:2] is the anchor's own pair
            # (the anchor terminal is neighbor[1]), so X_ref = refs[0] is bonded to it as InChI
            # requires, and the grouped order is refs itself.
            order = _ich_bond_order_from_refs(refs, X_ref, Y_ref)
            if order is None:
                continue

            rec = s0d + n_stereo
            rec.type = ICH_STEREO_DOUBLEBOND
            rec.central_atom = ICH_NO_ATOM
            rec.neighbor[0] = <short> mol._index_of[X_ref]
            rec.neighbor[1] = <short> mol._index_of[anchor_n]
            rec.neighbor[2] = <short> partner_idx
            rec.neighbor[3] = <short> mol._index_of[Y_ref]

            p = mol.translate_stereo(anchor_n, order)
            rec.parity = <signed char> _ich_chython_to_inchi(p, ICH_CIS_TRANS_FLIP)
            n_stereo += 1

        elif kind == SU_ALLENE:
            # InChI ALLENE: neighbor = {X, A, B, Y}; central_atom = center index.
            X_ref = refs[0] if not (unnamed_mask & 1) else None
            Y_ref = refs[2] if not (unnamed_mask & 4) else None

            if X_ref is None or Y_ref is None:
                continue

            A_idx, B_idx = _ich_find_allene_terminals(mol, <int> mol._index_of[anchor_n])
            if A_idx < 0:
                continue

            # ORIENT THE TERMINALS.  `_ich_find_allene_terminals` returns them in CSR edge order,
            # which is unrelated to which pair perception put first -- the anchor here is the CENTRE,
            # so refs[0:2] is simply one terminal's pair, not "the anchor's".  InChI requires X to be
            # bonded to A (neighbor[1]), so swap unless A already carries X_ref.  Without this the
            # record claims X sits on the far terminal and the parity is meaningless.
            if not _ich_is_bonded(mol, A_idx, <int> mol._index_of[X_ref]):
                A_idx, B_idx = B_idx, A_idx

            order = _ich_bond_order_from_refs(refs, X_ref, Y_ref)
            if order is None:
                continue

            rec = s0d + n_stereo
            rec.type = ICH_STEREO_ALLENE
            rec.central_atom = <short> mol._index_of[anchor_n]
            rec.neighbor[0] = <short> mol._index_of[X_ref]
            rec.neighbor[1] = <short> A_idx
            rec.neighbor[2] = <short> B_idx
            rec.neighbor[3] = <short> mol._index_of[Y_ref]

            p = mol.translate_stereo(anchor_n, order)
            rec.parity = <signed char> _ich_chython_to_inchi(p, ICH_ALLENE_FLIP)
            n_stereo += 1

        # SU_ATROPISOMER: no InChI layer; skip silently.

    return n_stereo


cdef int _ich_find_cis_trans_partner(MoleculeContainer mol, int anchor_idx) except -1:
    """Walk the cumulated double-bond chain from anchor_idx and return the partner terminal index.

    Returns -1 if no double-bond neighbour is found.
    """
    cdef uint32_t *ptr = csr_ptr(mol._structure)
    cdef halfedge_t *edges = csr_edges(mol._structure)
    cdef uint32_t k, remaining = mol._structure.header.atom_count + 1
    cdef halfedge_t *he
    cdef int cur = anchor_idx
    cdef int prev = -1
    cdef int nxt

    # Walk the chain: follow double bonds, never turning back.
    # Bounded by atom_count: a cumulene chain cannot exceed the molecule's atom count.
    while remaining > 0:
        remaining -= 1
        nxt = -1
        for k in range(ptr[cur], ptr[cur + 1]):
            he = edges + k
            if he.order == 2 and <int> he.to != prev:
                nxt = <int> he.to
                break
        if nxt < 0:
            return -1 if prev < 0 else cur
        prev = cur
        cur = nxt
    return -1  # budget exhausted (should not occur in valid molecules)


cdef tuple _ich_find_allene_terminals(MoleculeContainer mol, int center_idx):
    """Return (A_idx, B_idx) for the two terminals of an allene through center_idx.

    Returns (-1, -1) if fewer or more than two double-bond neighbours are found.
    """
    cdef uint32_t *ptr = csr_ptr(mol._structure)
    cdef halfedge_t *edges = csr_edges(mol._structure)
    cdef uint32_t k
    cdef halfedge_t *he
    cdef list terminals = []

    for k in range(ptr[center_idx], ptr[center_idx + 1]):
        he = edges + k
        if he.order == 2:
            terminals.append(<int> he.to)

    if len(terminals) == 2:
        return (<int> terminals[0], <int> terminals[1])
    return (-1, -1)


# ---- InChIKey --------------------------------------------------------------- #

def molecule_to_inchikey(MoleculeContainer mol not None) -> str:
    """Generate the standard InChIKey for `mol`.

    Calls GetStdINCHI then GetINCHIKeyFromINCHI.
    Raises ImportError if libinchi was not loaded.
    Raises ValueError on failure, or if the molecule carries an R marker.
    """
    cdef bytes inchi_b
    cdef bytes kb
    cdef char key_buf[28]
    cdef int rc = 0
    cdef str inchi_str

    _refuse_r_marker(mol)
    if ich_lib_handle is NULL:
        raise ImportError('libinchi not loaded; cannot generate InChIKey')

    inchi_str = molecule_to_inchi(mol, standard=True)
    inchi_b = inchi_str.encode()
    memset(key_buf, 0, 28)

    rc = ich_fn_GetINCHIKeyFromINCHI(
        <const char *> inchi_b, 0, 0, key_buf, NULL, NULL)

    if rc != 0:
        raise ValueError(f'InChIKey generation failed (rc={rc})')
    # InChIKey is always 27 ASCII characters; decode via explicit bytes slice to
    # avoid Cython's per-char overflow guard on C arrays (flagged unreachable by clang).
    kb = key_buf[:27]
    return kb.decode('ascii')


# ---- reverse direction: InChI string → MoleculeContainer ------------------- #

def inchi_to_molecule(str inchi_string not None) -> MoleculeContainer:
    """Parse an InChI string into a MoleculeContainer.

    Returns a new MoleculeContainer with element, charge, isotope, radical, implicit
    hydrogen count, and (for standard InChI) tetrahedral/double-bond/allene stereo.

    What is NOT preserved (by InChI design): atom-atom map numbers, enhanced stereo
    groups, atropisomers, XY coordinates, exact tautomer form.

    Raises ImportError if libinchi was not loaded.
    Raises ValueError on parse failure.
    """
    cdef bytes inchi_b = inchi_string.encode()
    cdef ICH_InputINCHI inp_inchi
    cdef ICH_OutputStruct out_struct
    cdef ICH_Atom *a_ptr
    cdef int na, i, j, order, nbr, rc = 0
    cdef int charge, isotope, atom_z
    cdef bint radical
    cdef int impl_h, iso_kind, iso_count, iso_remaining
    cdef str sym
    cdef str msg = ''
    cdef set seen
    cdef tuple bond_key
    cdef list idx_to_n
    cdef list iso_h_jobs  # list of (parent_n, isotope_kind, count) for isotopic H
    cdef object iso_h_item
    cdef uint32_t n, parent_n, h_n
    cdef MoleculeContainer mol

    if ich_lib_handle is NULL:
        raise ImportError('libinchi not loaded; cannot parse InChI')

    memset(&inp_inchi, 0, sizeof(ICH_InputINCHI))
    memset(&out_struct, 0, sizeof(ICH_OutputStruct))

    inp_inchi.szInChI = <char *> inchi_b
    inp_inchi.szOptions = NULL

    rc = ich_fn_GetStructFromINCHI(&inp_inchi, &out_struct)

    try:
        if rc > 1:
            msg = out_struct.szMessage.decode() if out_struct.szMessage is not NULL else 'unknown'
            raise ValueError(f'InChI parse failed (rc={rc}): {msg}')

        na = out_struct.num_atoms
        if na <= 0:
            raise ValueError('InChI parse returned empty structure')

        mol = MoleculeContainer()
        idx_to_n = []
        iso_h_jobs = []

        with mol:
            for i in range(na):
                a_ptr = &out_struct.atom[i]
                sym = a_ptr.elname.decode().rstrip('\x00')
                if sym not in _ich_sym_to_z:
                    raise ValueError(f'unknown element symbol {sym!r} in InChI output')
                atom_z = _ich_sym_to_z[sym]
                charge = int(a_ptr.charge)
                # InChI isotopic_mass encoding:
                #   0            → not isotopic
                #   1..8999      → absolute mass number (used only on input; rare in output)
                #   >= 10000     → ICH_ISOTOPIC_SHIFT_FLAG + delta_from_std_mass
                # GetStructFromINCHI always uses the >= 10000 format for isotopic atoms.
                if a_ptr.isotopic_mass > 0 and a_ptr.isotopic_mass < ICH_ISOTOPIC_SHIFT_FLAG:
                    isotope = int(a_ptr.isotopic_mass)
                elif a_ptr.isotopic_mass >= ICH_ISOTOPIC_SHIFT_FLAG:
                    # delta is signed; reconstruct absolute mass from rounded std mass
                    if atom_z >= 1 and atom_z <= 118:
                        isotope = ICH_STD_MASS[atom_z] + (int(a_ptr.isotopic_mass) - ICH_ISOTOPIC_SHIFT_FLAG)
                    else:
                        isotope = 0  # fallback: unknown element
                else:
                    isotope = 0
                radical = (a_ptr.radical != 0)
                # num_iso_H[0] = non-isotopic implicit H count; -1 means "auto"
                impl_h = int(a_ptr.num_iso_H[0])

                if impl_h < 0:
                    # InChI's -1 means "work it out from the valence rules".  Stored as the sentinel
                    # HERE and derived once the whole structure exists -- see the sweep after the
                    # bond block, which is the layer that can do what InChI is asking for.  It cannot
                    # be done in this loop: the answer depends on bonds that have not been added yet.
                    n = mol.add_atom(atom_z, charge=charge, isotope=isotope,
                                       radical=radical, implicit_h=H_UNKNOWN)
                else:
                    n = mol.add_atom(atom_z, charge=charge, isotope=isotope,
                                       radical=radical, implicit_h=impl_h)
                idx_to_n.append(n)

                # Collect isotopic H (protium=1, deuterium=2, tritium=3) to add as
                # explicit bonded atoms after all heavy atoms are committed.
                for iso_kind in range(1, 4):
                    iso_count = int(a_ptr.num_iso_H[iso_kind])
                    if iso_count > 0:
                        iso_h_jobs.append((n, iso_kind, iso_count))

            # Emit each bond once (lower index first).
            seen = set()
            for i in range(na):
                a_ptr = &out_struct.atom[i]
                for j in range(a_ptr.num_bonds):
                    nbr = int(a_ptr.neighbor[j])
                    if nbr <= i:
                        continue  # already emitted from nbr's side
                    bond_key = (i, nbr)
                    if bond_key in seen:
                        continue
                    seen.add(bond_key)
                    order = int(a_ptr.bond_type[j])
                    # InChI output is always Kekule (1/2/3); order 4 and 0 should not appear
                    if order < 1 or order > 3:
                        continue
                    mol.add_bond(idx_to_n[i], idx_to_n[nbr], order)

        # Add isotopic hydrogen atoms (protium, deuterium, tritium) as explicit bonded H.
        # InChI stores these in num_iso_H[1..3]; they must be explicit atoms in chython.
        if iso_h_jobs:
            with mol:
                for iso_h_item in iso_h_jobs:
                    parent_n = <uint32_t> iso_h_item[0]
                    iso_kind   = <int>      iso_h_item[1]
                    iso_count  = <int>      iso_h_item[2]
                    iso_remaining = iso_count
                    while iso_remaining > 0:
                        iso_remaining -= 1
                        h_n = mol.add_atom(1, isotope=iso_kind, implicit_h=0)
                        mol.add_bond(parent_n, h_n, 1)

        # DO WHAT InChI ASKED.  `num_iso_H[0] == -1` is not "unknown", it is an instruction: derive
        # this count from the valence rules.  Storing the sentinel and stopping there would make an
        # InChI the one supported input whose hydrogens are left underivable ON PURPOSE -- and it is
        # the input where the derivation is at its most reliable, because InChI output is
        # always Kekule (order 4 never appears; see the bond loop above), so no atom here is
        # pyrrole-versus-pyridine ambiguous and the ordinary rows answer or nothing does.
        #
        # AFTER the isotopic hydrogens, not before: a deuterium is an explicit neighbour and belongs
        # in the bond-order sum, so deriving first would answer for a skeleton the molecule no longer
        # has.  `fill_only=True` is what keeps this to the atoms InChI declined to state -- every count
        # libinchi did give is already in the arena and is left exactly as it came.
        #
        # Reached through the module globals rather than as a `cdef` call, because `_hydrogens.pxi` is
        # included after this file; that works because this is a `def`, the same route
        # `_molecule_container.pxi` documents and `kekule()` uses for its own heal.
        derive_implicit_hydrogens(mol, fill_only=True)

        # The with-block(s) have called _apply; atoms and bonds are now in the arena.
        # Apply stereo parities if the InChI carried any 0D stereo records.
        if out_struct.num_stereo0D > 0:
            _ich_apply_stereo(mol, &out_struct, idx_to_n)

        return mol

    finally:
        ich_fn_FreeStructFromINCHI(&out_struct)


cdef void _ich_apply_stereo(MoleculeContainer mol,
                             ICH_OutputStruct *out_struct,
                             list idx_to_n) except *:
    """Set stereo parities on `mol` from InChI's 0D stereo records.

    Called after atoms and bonds have been committed to the arena.
    For each record we try parity=1, call translate_stereo, and flip to parity=2 if needed.
    """
    cdef ICH_Stereo0D *rec
    cdef int i, typ, parity_ich, central_idx, n0, n1, n2, n3
    cdef int desired_chython, p_try, flipped
    cdef uint32_t anchor_n
    cdef tuple order
    cdef dict unit_by_anchor = {}
    cdef dict u
    cdef list units

    # Perceive stereo units; needed before any parity write.
    units = mol.stereo_units()
    for u in units:
        unit_by_anchor[u['anchor']] = u

    for i in range(out_struct.num_stereo0D):
        rec = &out_struct.stereo0D[i]
        typ = rec.type
        parity_ich = rec.parity & 0x07  # low 3 bits = connected-structure parity

        if parity_ich != ICH_PARITY_ODD and parity_ich != ICH_PARITY_EVEN:
            continue

        # InChI EVEN (2) ↔ chython EVEN (1); InChI ODD (1) ↔ chython ODD (2).  This is the
        # TETRAHEDRAL mapping, where both frames are the same four directions and no sign is in
        # question.  The two bond kinds recompute it through `_ich_inchi_to_chython` with their own
        # constant, because their frames differ from InChI's -- see the convention block above.
        desired_chython = 1 if parity_ich == ICH_PARITY_EVEN else 2

        if typ == ICH_STEREO_TETRAHEDRAL:
            central_idx = int(rec.central_atom)
            if central_idx < 0 or central_idx >= len(idx_to_n):
                continue
            anchor_n = <uint32_t> idx_to_n[central_idx]
            if anchor_n not in unit_by_anchor:
                continue

            n0 = int(rec.neighbor[0])
            n1 = int(rec.neighbor[1])
            n2 = int(rec.neighbor[2])
            n3 = int(rec.neighbor[3])

            mol.set_parity(anchor_n, 1)

            if n0 == central_idx:
                # 3-neighbour form: neighbor[0] == center (proxy for implicit H).
                # F26 sees (n1, n2, n3, None); InChI sends (center, n1, n2, n3).
                # The cyclic rotation is an odd permutation → parity flips relative to F26.
                p_try = mol.translate_stereo(anchor_n, (
                    <uint32_t> idx_to_n[n1],
                    <uint32_t> idx_to_n[n2],
                    <uint32_t> idx_to_n[n3],
                    None))
                # p_try is in F26 order; the InChI 3-neighbour parity is flipped vs F26
                flipped = 2 if p_try == 1 else 1
                if flipped != desired_chython:
                    mol.set_parity(anchor_n, 2)
            else:
                # 4-neighbour form: InChI order IS the F26 direction order.
                p_try = mol.translate_stereo(anchor_n, (
                    <uint32_t> idx_to_n[n0],
                    <uint32_t> idx_to_n[n1],
                    <uint32_t> idx_to_n[n2],
                    <uint32_t> idx_to_n[n3]))
                if p_try != desired_chython:
                    mol.set_parity(anchor_n, 2)

        elif typ == ICH_STEREO_DOUBLEBOND:
            # neighbor = {X, A, B, Y}; central_atom = NO_ATOM
            n0 = int(rec.neighbor[0])
            n1 = int(rec.neighbor[1])  # anchor A
            n2 = int(rec.neighbor[2])  # partner B
            n3 = int(rec.neighbor[3])

            if n1 < 0 or n1 >= len(idx_to_n):
                continue
            anchor_n = <uint32_t> idx_to_n[n1]
            if anchor_n not in unit_by_anchor:
                continue

            # X and Y are SUBSTITUENTS (on the two chain atoms n1, n2); the chain atoms themselves
            # are not directions.  Group them into the unit's own refs frame.
            if n0 < 0 or n3 < 0 or n0 >= len(idx_to_n) or n3 >= len(idx_to_n):
                continue
            order = _ich_bond_order_from_refs(
                unit_by_anchor[anchor_n]['refs'],
                idx_to_n[n0], idx_to_n[n3])
            if order is None:
                continue

            desired_chython = _ich_inchi_to_chython(parity_ich, ICH_CIS_TRANS_FLIP)
            mol.set_parity(anchor_n, 1)
            if mol.translate_stereo(anchor_n, order) != desired_chython:
                mol.set_parity(anchor_n, 2)

        elif typ == ICH_STEREO_ALLENE:
            central_idx = int(rec.central_atom)
            if central_idx < 0 or central_idx >= len(idx_to_n):
                continue
            anchor_n = <uint32_t> idx_to_n[central_idx]
            if anchor_n not in unit_by_anchor:
                continue

            n0 = int(rec.neighbor[0])
            n1 = int(rec.neighbor[1])  # terminal A
            n2 = int(rec.neighbor[2])  # terminal B
            n3 = int(rec.neighbor[3])

            # Same as the double bond: neighbor[1:3] are the chain terminals, neighbor[0] and [3]
            # are the two named substituents.  Only the latter are directions.
            if n0 < 0 or n3 < 0 or n0 >= len(idx_to_n) or n3 >= len(idx_to_n):
                continue
            order = _ich_bond_order_from_refs(
                unit_by_anchor[anchor_n]['refs'],
                idx_to_n[n0], idx_to_n[n3])
            if order is None:
                continue

            desired_chython = _ich_inchi_to_chython(parity_ich, ICH_ALLENE_FLIP)
            mol.set_parity(anchor_n, 1)
            if mol.translate_stereo(anchor_n, order) != desired_chython:
                mol.set_parity(anchor_n, 2)


# ---- the facade: one name per format, both directions ----------------------- #

def inchi(data, *, bint standard=True, str options=None):
    """A molecule from an InChI string, or an InChI string from a molecule.

    :param data: an `InChI=`-prefixed string, or a `MoleculeContainer`.
    :param standard: export only; `False` asks for a non-standard InChI.
    :param options: export only; extra flags, e.g. `'-SNon'`, with their prefix rewritten to the one
        this platform's libinchi reads.

    The InChIKey has its own name, :func:`inchikey`, because it is one-way: nothing reads a key back
    into a structure, so it is not a direction of this call.
    """
    if isinstance(data, MoleculeContainer):
        return molecule_to_inchi(data, standard=standard, options=options)
    elif isinstance(data, str):
        if not data.startswith('InChI='):
            raise ValueError("an InChI string starts with 'InChI='; an InChIKey cannot be read back "
                             'into a structure')
        return inchi_to_molecule(data)
    raise TypeError(f'inchi() takes a molecule or an InChI string, not {type(data).__name__}')


def inchikey(molecule):
    """`molecule`'s InChIKey.

    One way by construction -- the key is a hash of the InChI, so there is no reader.
    """
    return molecule_to_inchikey(molecule)
