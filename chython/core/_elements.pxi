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
# Element tables: the symbol table, the MDL reference-mass table, and the isotope tables.
#
# NONE OF THE DATA IN THIS FILE IS EDITABLE HERE.  It is compiled from `elements.tsv` and
# `isotopes.tsv` -- hand-maintained data files, and the authority -- by `gen_element_tables.py`.
# The generated region is delimited below and a test fails if it drifts from the two files.  The six
# arrays have to agree with each other -- the offset array is the prefix sum of the count array,
# which is the run lengths of the three flat arrays, and every MDL reference mass number has to
# appear among them -- and as C initialisers those relationships were unstated and unchecked.  Two
# of them were broken.
#
# Standalone data with no arena and no query in it.  Keeping `_elements.pxi` ahead of its readers is
# a convention here, not a correctness requirement -- see RULES.md §7.3.  `sig_mask()` at the foot
# of this file reads `SIG_MASK` from `_features.pxi`, which is included *after* this file.  It sits
# here deliberately as the live evidence for §7.3's claim that `cdef extern from *` blocks are
# hoisted: a §1.1 sweep should not move it to `_features.pxi`.
#
# Layout of the isotope arrays: flat parallel arrays behind a 119-entry offset/count index, because
# the access pattern is "give me all isotopes of element Z" -- one range scan, no hashing.  Index 0
# of the per-element arrays is unused so that the index IS the atomic number.

# --- BEGIN GENERATED TABLES: python chython/core/test/gen_element_tables.py ---
# Compiled from chython/core/elements.tsv and chython/core/isotopes.tsv, which are the
# authority and are maintained by hand.  Do not edit here -- run the command above.
# Both files are in compiled order, so the k-th entry here is the k-th row there.
#
# 436 nuclides over 118 elements.  Nobody has a mass for 2 of them
# (Db-270, Ts-297), and those compile to 0.0.
#
# Every mass number in MDL_ISOTOPE has a row among them, and the compile step refuses a
# pair of tables where one does not: a file is entitled to state the isotope MDL itself
# hands out, and a missing row makes `element_mass` answer 0.0 for it -- an atom of
# bromine-80, the mass number every MDL bromine measures against, weighing nothing.
cdef tuple SYMBOLS = (
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og')


cdef extern from *:
    """
    /* mass number MDL measures its mass-difference field from; index 0 unused */
    static const unsigned short MDL_ISOTOPE[119] = {
    0, 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23,
    24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51,
    52, 55, 56, 59, 59, 64, 65, 70, 73, 75, 79, 80,
    84, 85, 88, 89, 91, 93, 96, 98, 101, 103, 106, 108,
    112, 115, 119, 122, 128, 127, 131, 133, 137, 139, 140, 141,
    144, 145, 150, 152, 157, 159, 163, 165, 167, 169, 173, 175,
    178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209,
    209, 210, 222, 223, 226, 227, 232, 231, 238, 237, 244, 243,
    247, 247, 251, 252, 257, 258, 259, 260, 261, 270, 269, 270,
    270, 278, 281, 281, 285, 278, 289, 289, 293, 297, 294
    };
    /* group number convention; 0 is the f block, which states none; index 0 unused */
    static const unsigned char VALENCE_ELECTRONS[119] = {
    0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8, 1, 2, 3,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8, 1, 2, 3, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 3, 4, 5, 6, 7, 8, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8
    };
    /* calculated atomic radius in angstroms; index 0 unused */
    static const double ATOMIC_RADIUS[119] = {
    0.0, 0.53, 0.31, 1.67, 1.12, 0.87, 0.67, 0.56,
    0.48, 0.42, 0.38, 1.9, 1.45, 1.18, 1.11, 0.98,
    0.87, 0.79, 0.71, 2.43, 1.94, 1.84, 1.76, 1.71,
    1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36,
    1.25, 1.14, 1.03, 0.94, 0.87, 2.65, 2.19, 2.12,
    2.06, 1.98, 1.9, 1.83, 1.78, 1.73, 1.69, 1.65,
    1.61, 1.56, 1.45, 1.33, 1.23, 1.15, 1.08, 2.98,
    2.53, 2.12, 2.12, 2.47, 2.06, 2.05, 2.38, 2.31,
    2.33, 2.25, 2.28, 2.26, 2.26, 2.22, 2.22, 2.17,
    2.08, 2.0, 1.93, 1.88, 1.85, 1.8, 1.77, 1.74,
    1.71, 1.56, 1.54, 1.43, 1.35, 1.27, 1.2, 2.98,
    2.53, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17,
    2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17,
    2.08, 2.0, 1.93, 1.88, 1.85, 1.8, 1.77, 1.74,
    1.71, 1.56, 1.54, 1.43, 1.35, 1.27, 1.2
    };
    /* first row of element Z in the flat arrays; prefix sum of ISOTOPE_COUNTS */
    static const unsigned short ISOTOPE_OFFSETS[119] = {
    0, 0, 3, 5, 7, 8, 10, 14, 17, 21, 24, 27,
    29, 32, 33, 36, 39, 44, 47, 50, 54, 62, 64, 69,
    71, 76, 78, 84, 89, 96, 100, 108, 113, 118, 121, 130,
    136, 143, 146, 152, 155, 161, 162, 170, 172, 180, 182, 190,
    195, 203, 206, 217, 220, 228, 235, 246, 248, 255, 257, 261,
    262, 269, 270, 279, 282, 290, 292, 299, 301, 307, 309, 317,
    320, 326, 328, 333, 337, 345, 348, 354, 357, 366, 369, 374,
    377, 379, 381, 382, 383, 387, 389, 391, 393, 396, 397, 400,
    402, 406, 408, 410, 411, 412, 413, 414, 416, 418, 420, 421,
    422, 423, 424, 425, 427, 428, 430, 431, 432, 433, 435
    };
    /* how many rows element Z has */
    static const unsigned char ISOTOPE_COUNTS[119] = {
    0, 3, 2, 2, 1, 2, 4, 3, 4, 3, 3, 2, 3, 1, 3, 3, 5, 3, 3, 4,
    8, 2, 5, 2, 5, 2, 6, 5, 7, 4, 8, 5, 5, 3, 9, 6, 7, 3, 6, 3,
    6, 1, 8, 2, 8, 2, 8, 5, 8, 3, 11, 3, 8, 7, 11, 2, 7, 2, 4, 1,
    7, 1, 9, 3, 8, 2, 7, 2, 6, 2, 8, 3, 6, 2, 5, 4, 8, 3, 6, 3,
    9, 3, 5, 3, 2, 2, 1, 1, 4, 2, 2, 2, 3, 1, 3, 2, 4, 2, 2, 1,
    1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2, 1
    };
    /* mass number of the k-th row */
    static const unsigned short ISOTOPE_NUMBERS[436] = {
    1, 2, 3, 3, 4, 6, 7, 9, 10, 11, 11, 12, 13, 14, 13,
    14, 15, 15, 16, 17, 18, 17, 18, 19, 20, 21, 22, 22, 23, 24,
    25, 26, 27, 28, 29, 30, 31, 32, 33, 32, 33, 34, 35, 36, 35,
    36, 37, 36, 38, 40, 39, 40, 41, 42, 40, 42, 43, 44, 45, 46,
    47, 48, 44, 45, 46, 47, 48, 49, 50, 50, 51, 50, 51, 52, 53,
    54, 52, 55, 54, 55, 56, 57, 58, 59, 55, 57, 58, 59, 60, 58,
    59, 60, 61, 62, 63, 64, 63, 64, 65, 67, 62, 64, 65, 66, 67,
    68, 69, 70, 67, 68, 69, 70, 71, 70, 72, 73, 74, 76, 75, 76,
    77, 73, 74, 75, 76, 77, 78, 79, 80, 82, 76, 77, 79, 80, 81,
    82, 78, 80, 81, 82, 83, 84, 86, 82, 85, 87, 84, 85, 86, 87,
    88, 89, 86, 89, 90, 89, 90, 91, 92, 94, 96, 93, 92, 94, 95,
    96, 97, 98, 99, 100, 98, 99, 96, 98, 99, 100, 101, 102, 104, 106,
    103, 105, 102, 103, 104, 105, 106, 108, 109, 110, 107, 108, 109, 110, 111,
    106, 108, 110, 111, 112, 113, 114, 116, 111, 113, 115, 112, 113, 114, 115,
    116, 117, 118, 119, 120, 122, 124, 121, 122, 123, 120, 122, 123, 124, 125,
    126, 128, 130, 123, 124, 125, 127, 129, 131, 135, 124, 126, 127, 128, 129,
    130, 131, 132, 133, 134, 136, 131, 133, 130, 132, 134, 135, 136, 137, 138,
    138, 139, 136, 138, 140, 142, 141, 142, 143, 144, 145, 146, 148, 150, 145,
    144, 145, 147, 148, 149, 150, 152, 153, 154, 151, 152, 153, 152, 153, 154,
    155, 156, 157, 158, 160, 159, 160, 156, 158, 160, 161, 162, 163, 164, 165,
    166, 162, 164, 166, 167, 168, 170, 169, 170, 168, 169, 170, 171, 172, 173,
    174, 176, 175, 176, 177, 174, 176, 177, 178, 179, 180, 180, 181, 180, 182,
    183, 184, 186, 185, 186, 187, 188, 184, 186, 187, 188, 189, 190, 191, 192,
    191, 192, 193, 190, 192, 194, 195, 196, 198, 195, 197, 198, 196, 197, 198,
    199, 200, 201, 202, 203, 204, 203, 204, 205, 204, 206, 207, 208, 210, 207,
    209, 210, 209, 210, 210, 211, 222, 223, 223, 226, 228, 233, 225, 227, 227,
    232, 231, 233, 234, 235, 238, 237, 239, 242, 244, 241, 243, 243, 244, 247,
    248, 247, 249, 249, 251, 252, 257, 258, 259, 260, 266, 261, 267, 268, 270,
    269, 270, 270, 278, 281, 281, 282, 285, 278, 286, 289, 289, 293, 293, 297,
    294
    };
    /* exact mass in daltons; 0.0 where none is known */
    static const double ISOTOPE_MASSES[436] = {
    1.007825, 2.014102, 3.016049, 3.016029, 4.002603, 6.015122,
    7.016004, 9.012182, 10.012937, 11.009305, 11.011432, 12.0,
    13.003355, 14.003242, 13.005738, 14.003074, 15.000109, 15.003065,
    15.994915, 16.999132, 17.99916, 17.002095, 18.000938, 18.998403,
    19.99244, 20.993847, 21.991386, 21.994437, 22.98977, 23.985042,
    24.985837, 25.982593, 26.981538, 27.976927, 28.976495, 29.97377,
    30.973762, 31.973908, 32.971726, 31.972071, 32.971458, 33.967867,
    34.969032, 35.967081, 34.968853, 35.968307, 36.965903, 35.967546,
    37.962732, 39.962383, 38.963707, 39.963999, 40.961826, 41.962402,
    39.962591, 41.958618, 42.958767, 43.955481, 44.956186, 45.953693,
    46.954541, 47.952534, 43.959403, 44.95591, 45.95263, 46.951764,
    47.947947, 48.947871, 49.944792, 49.947163, 50.943964, 49.94605,
    50.944767, 51.940512, 52.940654, 53.938885, 51.945566, 54.93805,
    53.939615, 54.938293, 55.934942, 56.935399, 57.933281, 58.934876,
    54.941999, 56.936291, 57.935753, 58.9332, 59.933817, 57.935348,
    58.9343467, 59.930791, 60.93106, 61.928349, 62.929669, 63.92797,
    62.929601, 63.929764, 64.927794, 66.92773, 61.93433, 63.929147,
    64.929241, 65.926037, 66.927131, 67.924848, 68.92655, 69.925325,
    66.928202, 67.92798, 68.925581, 69.926022, 70.924705, 69.92425,
    71.922076, 72.923459, 73.921178, 75.921403, 74.921596, 75.922394,
    76.920647, 72.926765, 73.922477, 74.922523, 75.919214, 76.919915,
    77.91731, 78.9184991, 79.916522, 81.9167, 75.924541, 76.921379,
    78.918338, 79.9185293, 80.916291, 81.916804, 77.920386, 79.916378,
    80.916592, 81.913485, 82.914136, 83.911507, 85.91061, 81.918209,
    84.911789, 86.909183, 83.913425, 84.912933, 85.909262, 86.908879,
    87.905614, 88.907451, 85.914886, 88.905848, 89.907152, 88.90889,
    89.904704, 90.905645, 91.90504, 93.906316, 95.908276, 92.906378,
    91.90681, 93.905088, 94.905841, 95.904679, 96.906021, 97.905408,
    98.907712, 99.907477, 97.907216, 98.906255, 95.907598, 97.905287,
    98.905939, 99.90422, 100.905582, 101.904349, 103.90543, 105.907329,
    102.905504, 104.905694, 101.905608, 102.906087, 103.904035, 104.905084,
    105.903483, 107.903894, 108.90595, 109.905152, 106.905093, 107.905956,
    108.904756, 109.906107, 110.905291, 105.906458, 107.904183, 109.903006,
    110.904182, 111.902757, 112.904401, 113.903358, 115.904755, 110.905103,
    112.904061, 114.903878, 111.904821, 112.905171, 113.902782, 114.903346,
    115.901744, 116.902954, 117.901606, 118.903309, 119.902197, 121.90344,
    123.905275, 120.903818, 121.9051737, 122.904216, 119.90402, 121.903047,
    122.904273, 123.90282, 124.904425, 125.903306, 127.904461, 129.906223,
    122.905589, 123.90621, 124.90463, 126.904468, 128.904988, 130.906125,
    134.910048, 123.905896, 125.904269, 126.905184, 127.90353, 128.904779,
    129.903508, 130.905082, 131.904155, 132.905911, 133.905394, 135.90722,
    130.905464, 132.905447, 129.90631, 131.905056, 133.904503, 134.905683,
    135.90457, 136.905821, 137.905241, 137.907107, 138.906348, 135.90714,
    137.905986, 139.905434, 141.90924, 140.907648, 141.907719, 142.90981,
    143.910083, 144.912569, 145.913112, 147.916889, 149.920887, 144.912749,
    143.911995, 144.91341, 146.914893, 147.914818, 148.91718, 149.917271,
    151.919728, 152.922097, 153.922205, 150.919846, 151.921744, 152.921226,
    151.919788, 152.92175, 153.920862, 154.922619, 155.92212, 156.923957,
    157.924101, 159.927051, 158.925343, 159.927168, 155.924278, 157.924405,
    159.925194, 160.92693, 161.926795, 162.928728, 163.929171, 164.930319,
    165.932284, 161.928775, 163.929197, 165.93029, 166.932045, 167.932368,
    169.93546, 168.934211, 169.935801, 167.933894, 168.93519, 169.934759,
    170.936322, 171.936378, 172.938207, 173.938858, 175.942568, 174.940768,
    175.942682, 176.943758, 173.94004, 175.941402, 176.94322, 177.943698,
    178.945815, 179.946549, 179.947466, 180.947996, 179.946706, 181.948206,
    182.950224, 183.950933, 185.954362, 184.952956, 185.954986, 186.955751,
    187.958114, 183.952491, 185.953838, 186.955748, 187.955836, 188.958145,
    189.958445, 190.96093, 191.961479, 190.960591, 191.962605, 192.962924,
    189.95993, 191.961035, 193.962664, 194.964774, 195.964935, 197.967876,
    194.965035, 196.966552, 197.968244, 195.965815, 196.967213, 197.966752,
    198.968262, 199.968309, 200.970285, 201.970626, 202.972873, 203.973476,
    202.972329, 203.9738635, 204.974412, 203.973029, 205.974449, 206.975881,
    207.976636, 209.984189, 206.978471, 208.980383, 209.98412, 208.9824304,
    209.982874, 209.987155, 210.987496, 222.017578, 223.019736, 223.018502,
    226.02541, 228.03107, 233.048065, 225.02323, 227.027752, 227.027704,
    232.03805, 231.035879, 233.040247, 234.040946, 235.043923, 238.050783,
    237.048173, 239.052163, 242.058743, 244.064204, 241.056829, 243.06138,
    243.061389, 244.062753, 247.070354, 248.072349, 247.070307, 249.074987,
    249.074854, 251.079587, 252.08298, 257.095106, 258.098431, 259.10103,
    260.1055, 266.11983, 261.10877, 267.12153, 268.125676, 0.0,
    269.128634, 270.133363, 270.134293, 278.15481, 281.164516, 281.16636,
    282.169127, 285.177444, 278.17058, 286.182555, 289.190444, 289.0,
    293.204555, 293.0, 0.0, 294.0
    };
    /* natural terrestrial fraction; 0.0 for a nuclide with none */
    static const double ISOTOPE_ABUNDANCES[436] = {
    0.999885, 0.000115, 0.0, 1e-06, 0.999999, 0.0759,
    0.9241, 1.0, 0.199, 0.801, 0.0, 0.9893,
    0.0107, 0.0, 0.0, 0.99632, 0.00368, 0.0,
    0.99757, 0.00038, 0.00205, 0.0, 0.0, 1.0,
    0.9048, 0.0027, 0.0925, 0.0, 1.0, 0.7899,
    0.1, 0.1101, 1.0, 0.922296, 0.046832, 0.030872,
    1.0, 0.0, 0.0, 0.9493, 0.0076, 0.0429,
    0.0, 0.0002, 0.7578, 0.0, 0.2422, 0.003365,
    0.000632, 0.996003, 0.932581, 0.000117, 0.067302, 0.0,
    0.96941, 0.00647, 0.00135, 0.02086, 0.0, 4e-05,
    0.0, 0.00187, 0.0, 1.0, 0.0825, 0.0744,
    0.7372, 0.0541, 0.0518, 0.0025, 0.9975, 0.04345,
    0.0, 0.83789, 0.09501, 0.02365, 0.0, 1.0,
    0.05845, 0.0, 0.91754, 0.02119, 0.00282, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.680769,
    0.0, 0.262231, 0.011399, 0.036345, 0.0, 0.009256,
    0.6917, 0.0, 0.3083, 0.0, 0.0, 0.4863,
    0.0, 0.279, 0.041, 0.1875, 0.0, 0.0062,
    0.0, 0.0, 0.60108, 0.0, 0.39892, 0.2084,
    0.2754, 0.0773, 0.3628, 0.0761, 1.0, 0.0,
    0.0, 0.0, 0.0089, 0.0, 0.0937, 0.0763,
    0.2377, 0.0, 0.4961, 0.0873, 0.0, 0.0,
    0.5069, 0.0, 0.4931, 0.0, 0.0035, 0.0228,
    0.0, 0.1158, 0.1149, 0.57, 0.173, 0.0,
    0.7217, 0.2783, 0.0056, 0.0, 0.0986, 0.07,
    0.8258, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.5145, 0.1122, 0.1715, 0.1738, 0.028, 1.0,
    0.1484, 0.0925, 0.1592, 0.1668, 0.0955, 0.2413,
    0.0, 0.0963, 0.0, 1.0, 0.0554, 0.0187,
    0.1276, 0.126, 0.1706, 0.3155, 0.1862, 0.0,
    1.0, 0.0, 0.0102, 0.0, 0.1114, 0.2233,
    0.2733, 0.2646, 0.0, 0.1172, 0.51839, 0.0,
    0.48161, 0.0, 0.0, 0.0125, 0.0089, 0.1249,
    0.128, 0.2413, 0.1222, 0.2873, 0.0749, 0.0,
    0.0429, 0.9571, 0.0097, 0.0, 0.0066, 0.0034,
    0.1454, 0.0768, 0.2422, 0.0859, 0.3258, 0.0463,
    0.0579, 0.5721, 0.0, 0.4279, 0.0009, 0.0255,
    0.0089, 0.0474, 0.0707, 0.1884, 0.3174, 0.3408,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0009, 0.0009, 0.0, 0.0192, 0.2644,
    0.0408, 0.2118, 0.2689, 0.0, 0.1044, 0.0887,
    0.0, 1.0, 0.00106, 0.00101, 0.02417, 0.06592,
    0.07854, 0.11232, 0.71698, 0.0009, 0.9991, 0.00185,
    0.00251, 0.8845, 0.11114, 1.0, 0.272, 0.122,
    0.238, 0.083, 0.172, 0.057, 0.056, 1.0,
    0.0307, 0.0, 0.1499, 0.1124, 0.1382, 0.0738,
    0.2675, 0.0, 0.2275, 0.4781, 0.0, 0.5219,
    0.002, 0.0, 0.0218, 0.148, 0.2047, 0.1565,
    0.2484, 0.2186, 1.0, 0.0, 0.0006, 0.001,
    0.0234, 0.1891, 0.2551, 0.249, 0.2818, 1.0,
    0.0, 0.0014, 0.0161, 0.3361, 0.2293, 0.2678,
    0.1493, 1.0, 0.0, 0.0013, 0.0, 0.0304,
    0.1428, 0.2183, 0.1613, 0.3183, 0.1276, 0.9741,
    0.0259, 0.0, 0.0016, 0.0526, 0.186, 0.2728,
    0.1362, 0.3508, 0.00012, 0.99988, 0.0012, 0.265,
    0.1431, 0.3064, 0.2843, 0.374, 0.0, 0.626,
    0.0, 0.0002, 0.0159, 0.0196, 0.1324, 0.1615,
    0.2626, 0.0, 0.4078, 0.373, 0.0, 0.627,
    0.00014, 0.00782, 0.32967, 0.33832, 0.25242, 0.07163,
    0.0, 1.0, 0.0, 0.0015, 0.0, 0.0997,
    0.1687, 0.231, 0.1318, 0.2986, 0.0, 0.0687,
    0.29524, 0.0, 0.70476, 0.014, 0.241, 0.221,
    0.524, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    1.0, 1.0, 0.0, 5.5e-05, 0.0072, 0.992745,
    1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 0.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
    1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
    1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 0.0, 1.0
    };
    """
    const uint16_t MDL_ISOTOPE[119]
    const uint8_t  VALENCE_ELECTRONS[119]
    const double   ATOMIC_RADIUS[119]
    const uint16_t ISOTOPE_OFFSETS[119]
    const uint8_t  ISOTOPE_COUNTS[119]
    const uint16_t ISOTOPE_NUMBERS[436]
    const double   ISOTOPE_MASSES[436]
    const double   ISOTOPE_ABUNDANCES[436]

DEF ISOTOPE_ROWS = 436
# --- END GENERATED TABLES ---


# 0 is the reserved "chython states no count", and it is reserved because no element has zero valence
# electrons -- the same argument H_UNKNOWN makes for 15 (RULES.md §6.1, and the domain is declared
# here and nowhere else).  It covers the 28 f-block rows and an out-of-range atomic number.
DEF VALENCE_ELECTRONS_UNKNOWN = 0


cdef inline uint32_t el_valence_electrons(uint32_t z) noexcept nogil:
    """Valence electrons of element `z` under the group number convention, or 0 for the f block.

    The convention, and its bounds, are stated once in the header of `elements.tsv`: group number for
    groups 1-12, group number minus 10 for groups 13-18, so zinc is 12 and chlorine is 7.  For every
    main-group element that is Kier and Hall's Zv, which is what the valence connectivity indices in
    `_descriptors.pxi` need.

    RETURNS THE SENTINEL, NEVER RAISES, and every caller has to test it: cerium through lutetium and
    thorium through lawrencium state no count, because the f electrons are neither reliably core nor
    reliably valence.  A descriptor built on this refuses on such an atom rather than treating it as
    an element with no valence electrons.
    """
    if z < 1 or z > 118:
        return VALENCE_ELECTRONS_UNKNOWN
    return VALENCE_ELECTRONS[z]


def valence_electrons_table():
    """Expose VALENCE_ELECTRONS as a plain tuple, index 0 unused (test helper)."""
    cdef uint32_t i
    cdef list out = []
    for i in range(119):
        out.append(int(VALENCE_ELECTRONS[i]))
    return tuple(out)


cdef inline double el_atomic_radius(uint32_t z) noexcept nogil:
    """The calculated atomic radius of element `z`, in angstroms.

    ONE RADIUS AND IT IS THE CALCULATED ONE -- an SCF orbital measure, neither covalent nor van der
    Waals.  The header of `elements.tsv` states the column, including the 32 rows the published set does
    not reach, which carry the group analogue one period up.  The Hall-Kier alpha table in
    `_descriptors.pxi` is `r_cov / 0.77 - 1` per hybridization and answers a different question.

    0.0 IS NOT A RADIUS, and it is what element 0 gets: the R marker carries no radius the way it
    carries no mass.  An out-of-range atomic number gets it too.  A renderer draws no sphere at 0.0,
    which is the answer a marker wants -- unlike the valence electron count, no caller here has to test
    for a sentinel before doing arithmetic.
    """
    if z < 1 or z > 118:
        return 0.0
    return ATOMIC_RADIUS[z]


def atomic_radius_table():
    """Expose ATOMIC_RADIUS as a plain tuple, index 0 unused (test helper)."""
    cdef uint32_t i
    cdef list out = []
    for i in range(119):
        out.append(ATOMIC_RADIUS[i])
    return tuple(out)


cdef inline uint32_t el_period(uint32_t z) noexcept nogil:
    """The principal quantum number of the valence shell: the element's period.

    NOT A COLUMN, and it must not become one.  The period boundaries ARE the periodic table -- 2, 10,
    18, 36, 54, 86 are where the shells close -- so this is arithmetic over a fact rather than a
    convention anybody could disagree with, and a column would be a seventh place for it to drift.
    Kier and Hall's electrotopological state reads it as N.

    ELEMENT 0, THE R MARKER, HAS NO PERIOD and gets 1 here, which is arithmetic and not an answer.  No
    caller reaches it: `el_period`'s one reader is the intrinsic state in `_descriptors.pxi`, which
    takes `el_valence_electrons` first and that refuses on z outside 1-118.
    """
    if z <= 2:
        return 1
    if z <= 10:
        return 2
    if z <= 18:
        return 3
    if z <= 36:
        return 4
    if z <= 54:
        return 5
    if z <= 86:
        return 6
    return 7


def element_period(uint32_t z):
    """Expose `el_period` (test helper)."""
    return el_period(z)


cdef dict _build_symbol_table():
    # an explicit loop, not a comprehension: warn.undeclared bans comprehensions in .pyx
    cdef dict out = {}
    cdef uint32_t i
    for i in range(len(SYMBOLS)):
        out[SYMBOLS[i]] = i + 1
    # R is the fragment marker, element 0.  It is a symbol the readers and writers spell, and NOT a
    # SMARTS primitive: `[R]` in a query goes on meaning ring count.
    out['R'] = 0
    return out


cdef dict SYMBOL_TO_NUMBER = _build_symbol_table()

# `SYMBOL_TO_NUMBER['R'] == 0` is a LEGAL answer, so 0 cannot be the not-found sentinel.  Every lookup
# that distinguishes "absent" from "the marker" uses this compile-time constant, matching
# `_to_atomic_number`'s `except?` value.
DEF NOT_AN_ELEMENT = 0xffffffff


cdef uint32_t _to_atomic_number(element) except? 0xffffffff:
    if isinstance(element, int):
        if element < 0 or element > 118:
            raise ValueError('element must be an atomic number in 0-118, 0 being R')
        return <uint32_t> element
    if isinstance(element, str):
        # An index past `R_INDEX_MAX` is not a spelling of an element, so it is an unknown symbol
        # rather than an out-of-range index: this answers "which element is this symbol".  The bound
        # is named in the message because a caller who wrote `R500` meant an index.
        if element.startswith('R') and len(element) > 1 and element[1:].isdigit():
            if int(element[1:]) > R_INDEX_MAX:
                raise ValueError(f'unknown element symbol {element!r}: an R index stops at '
                                 f'{int(R_INDEX_MAX)}')
            return 0
        if element == 'R':
            return 0
        if element not in SYMBOL_TO_NUMBER:
            raise ValueError(f'unknown element symbol {element!r}')
        return <uint32_t> SYMBOL_TO_NUMBER[element]
    raise NotImplementedError(f'element must be an atomic number or a symbol, got {type(element)}')


def mdl_isotope_table():
    """Expose MDL_ISOTOPE to the test suite as a plain tuple."""
    cdef uint32_t i
    cdef list out = []
    for i in range(119):
        out.append(MDL_ISOTOPE[i])
    return tuple(out)


def element_symbols():
    """``SYMBOLS`` as a plain tuple, index 0 is R (the fragment marker) so that index == atomic number.

    Not a test hook: file writers need the atomic-number-to-symbol direction, and `add_atom`
    only provides the inverse.  Exposing the core's table is what keeps there being one symbol
    table in the library (RULES.md §6) instead of a second copy in every writer.
    """
    cdef uint32_t i
    cdef list out = ['R']
    for i in range(len(SYMBOLS)):
        out.append(SYMBOLS[i])
    return tuple(out)


cdef inline str symbol_of(atom_t *a):
    """The symbol of a STORED atom: ``'R'``, ``'R12'``, or the element's own.

    Not ``element_symbols()[element]``: the index rides in ``reserved``, so the symbol of an R is
    not a function of its element number alone.  Every subscript of ``SYMBOLS`` that reads an atom
    is this call instead -- ``SYMBOLS[0 - 1]`` answers ``'Og'`` and says nothing about being wrong.
    """
    cdef uint8_t index
    if a.element == 0:
        index = at_r_index(a)
        return 'R%d' % index if index else 'R'
    return SYMBOLS[a.element - 1]


cdef uint8_t _parse_r_index(str symbol) except? 0xff:
    """`'R'` -> 0, `'R12'` -> 12.  The domain is 0..R_INDEX_MAX, two decimal digits.

    Called only when ``len(symbol) > 1``; a plain ``'R'`` never arrives here.
    """
    cdef object value = int(symbol[1:])
    if value > R_INDEX_MAX:
        raise ValueError('R index %d is past R_INDEX_MAX (%d)' % (value, R_INDEX_MAX))
    return <uint8_t> value


def isotope_data(uint32_t z):
    """Return isotope data for element z as a tuple of (mass_number, exact_mass, abundance) triples.

    Raises ValueError if z is not in 1-118.  Returns an empty tuple for elements with no data.
    A triple whose exact_mass is 0.0 is a nuclide `isotopes.tsv` lists but has no measured mass
    for; there are two, and they are listed in the generated block's header comment.
    """
    cdef uint32_t i, off, cnt
    cdef list out = []
    if z < 1 or z > 118:
        raise ValueError(f'atomic number {z} out of range 1-118')
    off = ISOTOPE_OFFSETS[z]
    cnt = ISOTOPE_COUNTS[z]
    for i in range(off, off + cnt):
        out.append((int(ISOTOPE_NUMBERS[i]), ISOTOPE_MASSES[i], ISOTOPE_ABUNDANCES[i]))
    return tuple(out)


cdef double element_mass(uint32_t z, uint32_t isotope) noexcept nogil:
    """One atom's mass in daltons: the exact mass of `isotope`, or the abundance-weighted average
    over the natural isotopes when `isotope` is 0.

    Zero when `z` is out of range, when `isotope` names a mass number `isotopes.tsv` has no row for,
    and for the two rows it has no measured mass for.  Deliberately not an exception: the only
    caller is `float(molecule)`, and a molecule the arena accepted may legitimately carry an element
    with no measured abundances or a synthetic isotope nobody has weighed.  A mass of zero for such
    an atom is a visibly wrong number in a sum of masses; raising would make `float()` unusable on
    the record instead, and the arena's whole contract is that a record it holds can be asked
    questions.  Matches chython 2's `atomic_mass`, which sums `isotopes_distribution` for an
    unspecified isotope and indexes `isotopes_masses` for a specified one -- except on an isotope
    with no measured mass, where chython 2 raises `KeyError` and this answers 0.

    THE ZERO CASE EXCLUDES THE MDL REFERENCE MASS NUMBERS ON PURPOSE.  `MDL_ISOTOPE` frequently
    names a rounded standard atomic weight rather than an abundant nuclide -- bromine's is 80, not
    79 -- so `M  ISO` 80 on a bromine is a thing files say, and while those rows were missing such
    an atom weighed nothing.  Every mass number in `MDL_ISOTOPE` has a row in `isotopes.tsv` now,
    and `assert_invariants` refuses to compile a pair of tables where one does not.
    """
    cdef uint32_t i, off, cnt
    cdef double acc = 0.0
    if z < 1 or z > 118:
        return 0.0
    off = ISOTOPE_OFFSETS[z]
    cnt = ISOTOPE_COUNTS[z]
    if isotope:
        for i in range(off, off + cnt):
            if ISOTOPE_NUMBERS[i] == isotope:
                return ISOTOPE_MASSES[i]
        return 0.0
    for i in range(off, off + cnt):
        acc += ISOTOPE_ABUNDANCES[i] * ISOTOPE_MASSES[i]
    return acc


def isotope_offsets_table():
    """Expose ISOTOPE_OFFSETS as a plain tuple (test helper)."""
    cdef uint32_t i
    cdef list out = []
    for i in range(119):
        out.append(int(ISOTOPE_OFFSETS[i]))
    return tuple(out)


def isotope_counts_table():
    """Expose ISOTOPE_COUNTS as a plain tuple (test helper)."""
    cdef uint32_t i
    cdef list out = []
    for i in range(119):
        out.append(int(ISOTOPE_COUNTS[i]))
    return tuple(out)


def sig_mask():
    """Expose SIG_MASK to the test suite as a plain tuple."""
    cdef uint32_t i
    cdef list out = []
    for i in range(4):
        out.append(SIG_MASK[i])
    return tuple(out)


cdef inline bint element_is_heteroatom(uint8_t z) noexcept nogil:
    """Is this element a heteroatom? Not carbon, not hydrogen, and not the R fragment marker.

    Element 0 is the R marker: it reads as carbon for a neighbour's derived features, so it is not
    a heteroatom of the atom it caps.
    """
    return z != 6 and z != 1 and z != 0
