# -*- coding: utf-8 -*-
#
#  Copyright 2018-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ctypes import c_char, c_double, c_short, c_long, c_char_p, c_byte, POINTER, Structure, cdll, byref
from itertools import count
from sysconfig import get_platform
from warnings import warn
from .._convert import create_molecule
from ...containers import MoleculeContainer
from ...containers.bonds import Bond
from ...exceptions import ValenceError, IsChiral, NotChiral
from ...periodictable import H as _H


try:
    from importlib.resources import files, as_file
except ImportError:  # python3.8
    from importlib_resources import files, as_file


H = 1


def inchi(data, /, *, ignore_stereo: bool = False, _cls=MoleculeContainer) -> MoleculeContainer:
    """
    INCHI string parser
    """
    if lib is None:
        raise ImportError('libINCHI not found')

    structure = INCHIStructure()
    if lib.GetStructFromINCHI(byref(InputINCHI(data)), byref(structure)):
        lib.FreeStructFromINCHI(byref(structure))
        raise ValueError('invalid INCHI')

    atoms, bonds = [], []
    protium = {}
    deuterium = {}
    tritium = {}
    seen = set()
    for n in range(structure.num_atoms):
        seen.add(n)
        atom = structure.atom[n]

        atoms.append({'element': atom.atomic_symbol, 'charge': atom.charge, 'x': atom.x, 'y': atom.y,
                      'z': atom.z, 'isotope': atom.isotope, 'is_radical': atom.is_radical,
                      'implicit_hydrogens': atom.implicit_hydrogens, 'delta_isotope': atom.delta_isotope})
        if atom.implicit_protium:
            protium[n] = atom.implicit_protium
        if atom.implicit_deuterium:
            deuterium[n] = atom.implicit_deuterium
        if atom.implicit_tritium:
            tritium[n] = atom.implicit_tritium

        for k in range(atom.num_bonds):
            m = atom.neighbor[k]
            if m in seen:
                continue
            order = atom.bond_type[k]
            if order:
                bonds.append((n, m, order))

    stereo_atoms = []
    stereo_allenes = []
    stereo_cumulenes = []
    for i in range(structure.num_stereo0D):
        stereo = structure.stereo0D[i]
        sign = stereo.sign
        if sign is not None:
            if stereo.is_tetrahedral:
                stereo_atoms.append((stereo.central_atom, stereo.neighbors, sign))
            elif stereo.is_allene:
                nn, *_, nm = stereo.neighbors
                stereo_allenes.append((stereo.central_atom, nn, nm, sign))
            elif stereo.is_cumulene:
                nn, n, m, nm = stereo.neighbors
                stereo_cumulenes.append((n, m, nn, nm, sign))

    lib.FreeStructFromINCHI(byref(structure))

    tmp = {'atoms': atoms, 'bonds': bonds, 'stereo_atoms': stereo_atoms, 'stereo_allenes': stereo_allenes,
           'stereo_cumulenes': stereo_cumulenes, 'mapping': list(range(1, len(atoms) + 1)),
           'protium': protium, 'deuterium': deuterium, 'tritium': tritium}
    mol = create_molecule(tmp, skip_calc_implicit=True, _cls=_cls)
    postprocess_molecule(mol, tmp, ignore_stereo=ignore_stereo)
    return mol


def postprocess_molecule(molecule, data, *, ignore_stereo=False):
    atoms = molecule._atoms
    bonds = molecule._bonds

    # set hydrogen atoms. INCHI designed for hydrogens handling. hope correctly.
    free = count(len(atoms) + 1)
    to_add = []
    for n, atom in atoms.items():
        # in chython hydrogens never have implicit H. convert to explicit
        if atom == H and atom.implicit_hydrogens:
            for _ in range(atom.implicit_hydrogens):
                to_add.append((n, next(free), _H(implicit_hydrogens=0)))
            atom._implicit_hydrogens = 0

    for n, p in data['protium'].items():
        to_add.append((n + 1, next(free), _H(isotope=1, implicit_hydrogens=0)))
    for n, p in data['deuterium'].items():
        to_add.append((n + 1, next(free), _H(isotope=2, implicit_hydrogens=0)))
    for n, p in data['tritium'].items():
        to_add.append((n + 1, next(free), _H(isotope=3, implicit_hydrogens=0)))

    if to_add:
        for n, m, a in to_add:
            atoms[m] = a
            bonds[n][m] = b = Bond(1)
            bonds[m] = {n: b}
        molecule.calc_labels()  # reset labels

    if ignore_stereo or not data['stereo_atoms'] and not data['stereo_cumulenes'] and not data['stereo_allenes']:
        return

    st = molecule.stereogenic_tetrahedrons
    sa = molecule.stereogenic_allenes
    ctc = molecule._stereo_cis_trans_counterpart

    stereo = []
    for n, ngb, s in data['stereo_atoms']:
        n += 1
        if n in st:
            stereo.append((molecule.add_atom_stereo, n, [x + 1 for x in ngb], s))
    for n, nn, mn, s in data['stereo_allenes']:
        n += 1
        if n in sa:
            stereo.append((molecule.add_atom_stereo, n, nn + 1, mn + 1, s))
    for n, m, nn, nm, s in data['stereo_cumulenes']:
        n += 1
        if n in ctc:
            stereo.append((molecule.add_cis_trans_stereo, n, m + 1, nn + 1, nm + 1, s))

    while stereo:
        fail_stereo = []
        old_stereo = len(stereo)
        for f, *args in stereo:
            try:
                f(*args, clean_cache=False)
            except NotChiral:
                fail_stereo.append((f, *args))
            except IsChiral:
                pass
            except ValenceError:
                if 'chython_parsing_log' not in molecule.meta:
                    molecule.meta['chython_parsing_log'] = []
                molecule.meta['chython_parsing_log'].append('structure has errors, stereo data skipped')
                molecule.flush_cache()
                break
        else:
            stereo = fail_stereo
            if len(stereo) == old_stereo:
                break
            molecule.flush_stereo_cache()
            continue
        break


class InputINCHI(Structure):
    def __init__(self, string, options=None):
        if options is None:
            options = b''
        else:
            options = ' '.join(f'{opt_flag}{x}' for x in options).encode()
        super().__init__(string.encode(), options)

    _fields_ = [('szInChI', c_char_p),   # InChI ASCII string to be converted to a strucure
                ('szOptions', c_char_p)  # InChI options: space-delimited; each is preceded
                                         # by '/' or '-' depending on OS and compiler
                ]


class Atom(Structure):
    @property
    def atomic_symbol(self):
        return self.elname.decode()

    @property
    def isotope(self):
        if 0 < self.isotopic_mass < 9000:  # OVER NINE THOUSANDS!
            return self.isotopic_mass

    @property
    def delta_isotope(self):
        if self.isotopic_mass > 9000:
            return self.isotopic_mass - 10_000

    @property
    def is_radical(self):
        return bool(self.radical)

    @property
    def implicit_hydrogens(self):
        return self.num_iso_H[0]

    @property
    def implicit_protium(self):
        return self.num_iso_H[1]

    @property
    def implicit_deuterium(self):
        return self.num_iso_H[2]

    @property
    def implicit_tritium(self):
        return self.num_iso_H[3]

    _fields_ = [('x', c_double), ('y', c_double), ('z', c_double),  # atom coordinates
                ('neighbor', c_short * 20),    # adjacency list: ordering numbers of the adjacent atoms, >= 0
                ('bond_type', c_byte * 20),    # inchi_BondType
                ('bond_stereo', c_byte * 20),  # inchi_BondStereo2D; negative if the sharp end points to opposite atom
                ('elname', c_char * 6),        # zero-terminated chemical element name: "H", "Si", etc.
                ('num_bonds', c_short),        # number of neighbors, bond types and bond stereo in the adjacency list
                ('num_iso_H', c_byte * 4),     # implicit hydrogen atoms
                                               # [0]: number of implicit non-isotopic H
                                               # (exception: num_iso_H[0]=-1 means INCHI adds implicit H automatically),
                                               # [1]: number of implicit isotopic 1H (protium),
                                               # [2]: number of implicit 2H (deuterium),
                                               # [3]: number of implicit 3H (tritium)
                ('isotopic_mass', c_short),    # 0 => non-isotopic; isotopic mass or 10000 + mass - average atomic mass
                ('radical', c_byte),           # inchi_Radical,
                ('charge', c_byte)]            # positive or negative; 0 => no charge


class Stereo0D(Structure):
    @property
    def is_tetrahedral(self):
        return self.type == 2

    @property
    def is_allene(self):
        return self.type == 3

    @property
    def is_cumulene(self):
        return self.type == 1

    @property
    def neighbors(self):
        return tuple(self.neighbor)

    @property
    def sign(self):
        if self.parity == 1:
            return True
        elif self.parity == 2:
            return False

    _fields_ = [('neighbor', c_short * 4),  # 4 atoms always
                ('central_atom', c_short),  # central tetrahedral atom or a central atom of allene; otherwise NO_ATOM
                ('type', c_byte),           # inchi_StereoType0D
                ('parity', c_byte)]         # inchi_StereoParity0D


class INCHIStructure(Structure):
    _fields_ = [('atom', POINTER(Atom)),          # array of num_atoms elements
                ('stereo0D', POINTER(Stereo0D)),  # array of num_stereo0D 0D stereo elements or NULL
                ('num_atoms', c_short),           # number of atoms in the structure
                ('num_stereo0D', c_short),        # number of 0D stereo elements
                ('szMessage', c_char_p),          # Error/warning ASCII message
                ('szLog', c_char_p),              # log-file ASCII string, contains a human-readable list
                                                  # of recognized options and possibly an Error/warn message
                ('WarningFlags', (c_long * 2) * 2)]

# copy-pasted from INCHI-API
#  * Notes: 1. Atom ordering numbers (i, k, and atom[i].neighbor[j] below)
#  *           start from zero; max. ordering number is (num_atoms-1).
#  *        2. inchi_Atom atom[i] is connected to the atom[atom[i].neighbor[j]]
#  *           by a bond that has type atom[i].bond_type[j] and 2D stereo type
#  *           atom[i].bond_stereo[j] (in case of no stereo
#  *           atom[i].bond_stereo[j] = INCHI_BOND_STEREO_NONE)
#  *           Index j is in the range 0 <= j <= (atom[i].num_bonds-1)
#  *        3. Any connection (represented by atom[i].neighbor[j],
#  *           atom[i].bond_type[j], and atom[i].bond_stereo[j])
#  *           should be present in one or both adjacency list:
#  *             if k = atom[i].neighbor[j] then i may or may not be present in
#  *           atom[k].neighbor[] list. For example, the adjacency lists may be
#  *           populated with only such neighbors that atom[i].neighbor[j] < i
#  *           All elements of an adjacency list must be different, that is,
#  *           a bond must be specified in an adjacency list only once.
#  *        4. in Molfiles usually
#  *           (number of implicit H) = Valence - SUM(bond_type[])
#  *        5. Seemingly illogical order of the inchi_Atom members was
#  *           chosen in an attempt to avoid alignment problems when
#  *           accessing inchi_Atom from unrelated to C programming
#  *           languages such as Visual Basic.
#  *******************************************************************/
#
# /*******************************************************************
#     0D Stereo Parity and Type definitions
#  *******************************************************************
#             Note:
#             =====
#             o Below #A is the ordering number of atom A, starting from 0
#             o See parity values corresponding to 'o', 'e', and 'u' in
#               inchi_StereoParity0D definition below)
#
#            =============================================
#             stereogenic bond >A=B< or cumulene >A=C=C=B<
#            =============================================
#
#                                  neighbor[4]  : {#X,#A,#B,#Y} in this order
#      X                           central_atom : NO_ATOM
#       \            X      Y      type         : INCHI_StereoType_DoubleBond
#        A==B         \    /
#            \         A==B
#             Y
#
#     parity= 'e'    parity= 'o'   unknown parity = 'u'
#
#     Limitations:
#     ============
#     o Atoms A and B in cumulenes MUST be connected by a chain of double bonds;
#       atoms A and B in a stereogenic 'double bond' may be connected by a double,
#       single, or alternating bond.
#     o One atom may belong to up to 3 stereogenic bonds (i.g. in a fused
#       aromatic structure).
#     o Multiple stereogenic bonds incident to any given atom should
#       either all except possibly one have (possibly different) defined
#       parities ('o' or 'e') or should all have an unknown parity 'u'.
#
#       Note on parities of alternating stereobonds
#       ===========================================
#                                                      D--E
#       In large rings  (see Fig. 1, all              //   \\
#       atoms are C) all alternating bonds         B--C      F--G
#       are treated as stereogenic.              //              \\
#       To avoid "undefined" bond parities      A                  H
#       for bonds BC, DE, FG, HI, JK, LM, AN     \               /
#       it is recommended to mark them with       N==M       J==I
#       parities.                                     \     /
#                                                       L==K    Fig. 1
#       Such a marking will make
#       the stereochemical layer unambiguous
#       and it will be different from the          B--C      F--G
#       stereochemical layer of the second       //   \\   //    \\
#       structure (Fig. 2).                     A      D--E        H
#                                                \               /
#                                                 N==M       J==I
#       By default, double and alternating            \     /
#       bonds in 8-member and greater rings             L==K    Fig. 2
#       are treated by InChI as stereogenic.
#
#
#            =============================================
#             tetrahedral atom
#            =============================================
#
#    4 neighbors
#
#             X                    neighbor[4] : {#W, #X, #Y, #Z}
#             |                    central_atom: #A
#          W--A--Y                 type        : INCHI_StereoType_Tetrahedral
#             |
#             Z
#    parity: if (X,Y,Z) are clockwize when seen from W then parity is 'e' otherwise 'o'
#    Example (see AXYZW above): if W is above the plane XYZ then parity = 'e'
#
#    3 neighbors
#
#               Y          Y       neighbor[4] : {#A, #X, #Y, #Z}
#              /          /        central_atom: #A
#          X--A  (e.g. O=S   )     type        : INCHI_StereoType_Tetrahedral
#              \          \
#               Z          Z
#
#    parity: if (X,Y,Z) are clockwize when seen from A then parity is 'e',
#                                                           otherwise 'o'
#    unknown parity = 'u'
#    Example (see AXYZ above): if A is above the plane XYZ then parity = 'e'
#    This approach may be used also in case of an implicit H attached to A.
#
#            =============================================
#             allene
#            =============================================
#
#        X       Y                 neighbor[4]  : {#X,#A,#B,#Y}
#         \     /                  central_atom : #C
#          A=C=B                   type         : INCHI_StereoType_Allene
#
#                                       Y      X
#                                       |      |
#      when seen from A along A=C=B:  X-A    Y-A
#
#                           parity:   'e'    'o'
#
#    parity: if A, B, Y are clockwise when seen from X then parity is 'e',
#                                                           otherwise 'o'
#    unknown parity = 'u'
#    Example (see XACBY above): if X on the diagram is above the plane ABY
#                                                       then parity is 'o'
#
#    Limitations
#    ===========
#    o Atoms A and B in allenes MUST be connected by a chain of double bonds;
#
#
#    How InChI uses 0D parities
#    ==========================
#
#    1. 0D parities are used if all atom coordinates are zeroes.
#
#    In addition to that:
#
#    2. 0D parities are used for Stereobonds, Allenes, or Cumulenes if:
#
#    2a. A bond to the end-atom is shorter than MIN_BOND_LEN=0.000001
#    2b. A ratio of two bond lengths to the end-atom is smaller than MIN_SINE=0.03
#    2c. In case of a linear fragment X-A=B end-atom A is treated as satisfying 2a-b
#
#        0D parities are used if 2a or 2b or 2c applies to one or both end-atoms.
#
#    3. 0D parities are used for Tetrahedral Atoms if at least one of 3a-c is true:
#
#    3a. One of bonds to the central atom is shorter than MIN_BOND_LEN=0.000001
#    3b. A ratio of two bond lengths to the central atom is smaller than MIN_SINE=0.03
#    3c. The four neighbors are almost in one plane or the central atom and
#        its only 3 explicit neighbors are almost in one plane
#
#    Notes on 0D parities and 'undefined' stereogenic elements
#    =========================================================
#
#    If 0D parity is to be used according to 1-3 but    CH3     CH3
#    has not been provided then the corresponding         \    /
#    stereogenic element is considered 'undefined'.        C=CH
#                                                         /
#    For example, if in the structure (Fig. 3)           H
#    the explicit H has been moved so that it                Fig. 3
#    has same coordinates as atom >C= (that is,
#    the length of the bond H-C became zero)
#    then the double bond is assigned 'undefined'       CH3      CH3
#    parity which by default is omitted from the          \     /
#    Identifier.                                           CH=CH
#
#    However, the structure on Fig. 4 will have double        Fig. 4
#    bond parity 'o' and its parity in the Identifier is (-).
#
#    Notes on 0D parities in structures containing metals
#    ====================================================
#    Since InChI disconnects bonds to metals the 0D parities upon the
#    disconnection may change in several different ways:
#
#    1) previously non-stereogenic bond may become stereogenic:
#
#          \     /                            \     /
#           CH==CH          disconnection      CH==CH
#            \ /               ======>
#             M                                  M
#
#      before the disconnection:    after the disconnection:
#      atoms C have valence=5 and   the double bond may become
#      the double bond is not       stereogenic
#      recognized as stereogenic
#
#    2) previously stereogenic bond may become non-stereogenic:
#
#        M                           M(+)
#         \    /                             /
#          N==C      disconnection    (-)N==C
#              \        ======>              \
#
#    3) Oddball structures, usually resulting from projecting 3D
#       structures on the plane, may contain fragment like that
#       depicted on Fig. 5:
#
#               M   A                      M   A
#               |\ /   B                      /   B
#               | X   /     disconnection    /   /
#               |/ \ /         ======>      /   /
#               C===C                      C===C
#              Fig. 5
#      (X stands for bond intersection)
#
#      A-C=C-B parity is              A-C=C-B parity is
#      trans (e)                      cis (o) or undefined
#      because the bond               because C valence = 3,
#      orientation is same            not 4.
#      as on Fig, 6 below:
#
#           A       M
#            \     /     Removal of M from the structure
#             C===C      on Fig. 5 changes the geometry from trans
#            /     \     to cis.
#           M'      B    Removal of M and M' from the structure
#           Fig. 6       on Fig. 6 does not change the A-C=C-B
#                        geometry: it is trans.
#
#    To resolve the problem InChI API accepts the second parity
#    corresponding to the metal-disconnected structure.
#    To store both bond parities use left shift by 3 bits:
#
#    inchi_Stereo0D::parity = ParityOfConnected | (ParityOfDisconnected<<3)
#
#    In case when only disconnected structure parity exists set
#    ParityOfConnected = INCHI_PARITY_UNDEFINED.
#    This is the only case when INCHI_PARITY_UNDEFINED parity
#    may be fed to the InChI.
#
#    In cases when the bond parity in a disconnected structure exists and
#    differs from the parity in the connected structure the atoms A and B
#    should be non-metals.
#


lib = None

platform = get_platform()
if platform == 'win-amd64':
    opt_flag = '/'
    libname = 'libinchi.dll'
elif platform == 'linux-x86_64':
    opt_flag = '-'
    libname = 'libinchi.so'
elif platform.startswith('macosx') and platform.endswith('x86_64'):
    opt_flag = '-'
    libname = 'libinchi.dynlib'
elif platform.startswith('macosx') and platform.endswith('arm64'):
    opt_flag = '-'
    libname = 'libinchi_arm64.dylib'
else:
    warn('unsupported platform for libinchi', ImportWarning)
    libname = None

if libname:
    file = files(__package__).joinpath(libname)
    if file.is_file():
        with as_file(file) as f:
            try:
                lib = cdll.LoadLibrary(str(f))
            except OSError:
                warn('libinchi loading problem', ImportWarning)
    else:
        warn('broken package installation. libinchi not found', ImportWarning)


__all__ = ['inchi']
