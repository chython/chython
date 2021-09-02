# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from distutils.util import get_platform
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from itertools import count
from os import name
from pathlib import Path
from re import split
from sys import prefix, exec_prefix
from typing import List, Iterator
from warnings import warn
from ._mdl import Parser, common_isotopes, parse_error
from ..containers import MoleculeContainer
from ..containers.bonds import Bond
from ..exceptions import ValenceError, IsChiral, NotChiral
from ..periodictable import H


class INCHIRead(Parser):
    """
    INCHI separated per lines files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.
    line should be start with INCHI string and
    optionally continues with space/tab separated list of key:value [or key=value] data if header=None.
        example:
            InChI=1S/C2H5/c1-2/h1H2,2H3/q+1 id:123 key=value
    if header=True then first line of file should be space/tab separated list of keys including INCHI column key.
        example:
            ignored_inchi_key key1 key2
            InChI=1S/C2H5/c1-2/h1H2,2H3/q+1 1 2
    also possible to pass list of keys (without inchi_pseudo_key) for mapping space/tab separated list
    of INCHI and values: header=['key1', 'key2'] # order depended
    """
    def __init__(self, file, header=None, ignore_stereo=False, **kwargs):
        """
        :param ignore: Skip some checks of data or try to fix some errors.
        :param remap: Remap atom numbers started from one.
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        :param calc_cis_trans: Calculate cis/trans marks from 2d coordinates.
        :param ignore_stereo: Ignore stereo data.
        """
        if isinstance(file, str):
            self._file = open(file)
            self.__is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self.__is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO, FileInput)):
            self._file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)
        self.__file = iter(self._file.readline, '')

        if header is True:
            self.__header = next(self.__file).split()[1:]
        elif header:
            if not isinstance(header, (list, tuple)) or not all(isinstance(x, str) for x in header):
                raise TypeError('expected list (tuple) of strings')
            self.__header = header
        else:
            self.__header = None

        self.__ignore_stereo = ignore_stereo
        self._data = self.__data()

    def __data(self):
        file = self._file
        parse = self.parse
        try:
            seekable = file.seekable()
        except AttributeError:
            seekable = False
        pos = file.tell() if seekable else None
        for n, line in enumerate(self.__file):
            try:
                x = parse(line)
            except ValueError:
                yield parse_error(n, pos, self._format_log(), line)
            else:
                yield x
            if seekable:
                pos = file.tell()

    @classmethod
    def create_parser(cls, header=None, ignore_stereo=False, *args, **kwargs):
        """
        Create INCHI parser function configured same as INCHIRead object
        """
        obj = object.__new__(cls)
        obj._INCHIRead__header = header
        obj._INCHIRead__ignore_stereo = ignore_stereo
        super(INCHIRead, obj).__init__(*args, **kwargs)
        return obj.parse

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[MoleculeContainer]:
        """
        parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self) -> Iterator[MoleculeContainer]:
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self) -> MoleculeContainer:
        return next(iter(self))

    def parse(self, inchi: str) -> MoleculeContainer:
        """
        convert INCHI string into MoleculeContainer object. string should be start with INCHI and
        optionally continues with space/tab separated list of key:value [or key=value] data.
        """
        if not inchi:
            raise ValueError('Empty string')
        self._flush_log()

        inchi, *data = inchi.split()
        if self.__header is None:
            meta = {}
            for x in data:
                try:
                    k, v = split('[=:]', x, 1)
                    meta[k] = v
                except ValueError:
                    self._info(f'invalid metadata entry: {x}')
        else:
            meta = dict(zip(self.__header, data))

        record = self.__parse_inchi(inchi)
        record['meta'] = meta
        return self._convert_molecule(record)

    def _create_molecule(self, data, mapping):
        mol = super()._create_molecule(data, mapping, _skip_calc_implicit=True)
        atoms = mol._atoms
        bonds = mol._bonds
        hydrogens = mol._hydrogens

        # set hydrogen atoms. INCHI designed for hydrogens handling. hope correctly.
        free = count(max(mapping.values()) + 1)
        for n, atom in enumerate(data['atoms']):
            n = mapping[n]
            if atom['element'] != 'H':
                hydrogens[n] = atom['hydrogens']
            # in chython hydrogens never have implicit H.
            elif atom['hydrogens']:  # >[xH]-H case
                m = next(free)
                atoms[m] = a = H()
                a._attach_to_graph(mol, m)
                bonds[n][m] = b = Bond(1)
                bonds[m] = {n: b}
                hydrogens[n] = 0
                hydrogens[m] = 0
            else:  # H+, H* or >H-[xH] cases
                hydrogens[n] = 0
            # convert isotopic implicit hydrogens to explicit
            for i, k in enumerate(('p', 'd', 't'), 1):
                if atom[k]:
                    for _ in range(atom[k]):
                        m = next(free)
                        atoms[m] = a = H(i)
                        a._attach_to_graph(mol, m)
                        bonds[n][m] = b = Bond(1)
                        bonds[m] = {n: b}
                        hydrogens[m] = 0

        if self.__ignore_stereo or \
                not data['stereo_atoms'] and not data['stereo_cumulenes'] and not data['stereo_allenes']:
            return mol

        st = mol._stereo_tetrahedrons
        sa = mol._stereo_allenes
        ctt = mol._stereo_cis_trans_terminals

        stereo = []
        for n, ngb, s in data['stereo_atoms']:
            n = mapping[n]
            if n in st:
                stereo.append((mol.add_atom_stereo, n, [mapping[x] for x in ngb], s))
        for n, nn, mn, s in data['stereo_allenes']:
            n = mapping[n]
            if n in sa:
                stereo.append((mol.add_atom_stereo, n, mapping[nn], mapping[mn], s))
        for n, m, nn, nm, s in data['stereo_cumulenes']:
            n = mapping[n]
            if n in ctt:
                stereo.append((mol.add_cis_trans_stereo, n, mapping[m], mapping[nn], mapping[nm], s))

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
                    self._info('structure has errors, stereo data skipped')
                    mol.flush_cache()
                    break
            else:
                stereo = fail_stereo
                if len(stereo) == old_stereo:
                    break
                mol.flush_stereo_cache()
                continue
            break
        return mol

    @staticmethod
    def __parse_inchi(string):
        structure = INCHIStructure()
        if lib.GetStructFromINCHI(byref(InputINCHI(string)), byref(structure)):
            lib.FreeStructFromINCHI(byref(structure))
            raise ValueError('invalid INCHI')

        atoms, bonds = [], []
        seen = set()
        for n in range(structure.num_atoms):
            seen.add(n)
            atom = structure.atom[n]

            atoms.append({'element': atom.atomic_symbol, 'charge': atom.charge,
                          'mapping': 0, 'x': atom.x, 'y': atom.y, 'z': atom.z, 'isotope': atom.isotope,
                          'is_radical': atom.is_radical, 'hydrogens': atom.implicit_hydrogens,
                          'p': atom.implicit_protium, 'd': atom.implicit_deuterium, 't': atom.implicit_tritium})

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
        return {'atoms': atoms, 'bonds': bonds, 'stereo_atoms': stereo_atoms, 'stereo_allenes': stereo_allenes,
                'stereo_cumulenes': stereo_cumulenes}


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
        isotope = self.isotopic_mass
        if not isotope:
            isotope = None
        elif isotope > 9000:  # OVER NINE THOUSANDS!
            isotope += common_isotopes[self.atomic_symbol] - 10000
        return isotope

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


try:
    from site import getuserbase
except ImportError:
    prefixes = {prefix, exec_prefix}
else:
    user_prefix = getuserbase()
    if user_prefix:
        prefixes = {prefix, exec_prefix, user_prefix}
    else:
        prefixes = {prefix, exec_prefix}

sitepackages = []
for pr in prefixes:
    pr = Path(pr)
    if name == 'posix':
        sitepackages.append(pr / 'local' / 'lib')
    else:
        sitepackages.append(pr)
    sitepackages.append(pr / 'lib')

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
else:
    warn('unsupported platform for libinchi', ImportWarning)
    libname = None
    __all__ = []
    del INCHIRead

if libname:
    for site in sitepackages:
        lib_path = site / libname
        if lib_path.exists():
            try:
                lib = cdll.LoadLibrary(str(lib_path))
            except OSError:
                warn('libinchi loading problem', ImportWarning)
                __all__ = []
                del INCHIRead
                break
            __all__ = ['INCHIRead']
            break
    else:
        warn('broken package installation. libinchi not found', ImportWarning)
        __all__ = []
        del INCHIRead
