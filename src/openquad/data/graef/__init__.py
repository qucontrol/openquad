# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

"""S2 and SO3 quadratures from numerical error optimization from M. Gräf.

The files in this directory have been converted from the data from Dr. Manuel
Gräf published at [1] and downloaded in January 2024.

[1]: https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/
"""

S2_GAUSS_METHODS = [
    # degree, size
    [1, 2],
    [2, 4],
    [3, 6],
    [4, 10],
    [5, 12],
    [6, 18],
    [7, 22],
    [8, 28],
    [9, 32],
    [10, 42],
    [11, 48],
    [12, 58],
    [13, 64],
    [14, 72],
    [15, 82],
    [16, 98],
    [17, 104],
    [18, 122],
    [19, 130],
    [20, 148],
    [21, 156],
    [22, 178],
    [23, 186],
    [24, 210],
    [25, 220],
    [26, 244],
    [27, 254],
    [28, 282],
    [29, 292],
    [30, 322],
    [32, 364],
    [34, 410],
    [35, 422],
    [36, 458],
    [37, 472],
    [38, 508],
    [39, 522],
    [44, 672],
]
"""list: Table of S2 Gauss-type quadratures from M. Gräf."""

S2_SPHERICAL_DESIGNS = [
    # degree, size
    [1, 2],
    [2, 4],
    [3, 6],
    [5, 12],
    [7, 24],
    [8, 36],
    [9, 48],
    [10, 60],
    [11, 70],
    [12, 84],
    [13, 94],
    [14, 108],
    [15, 120],
    [16, 144],
    [17, 156],
    [18, 180],
    [19, 192],
    [20, 216],
    [21, 234],
    [22, 264],
    [23, 278],
    [24, 312],
    [25, 328],
    [26, 360],
    [27, 380],
    [28, 420],
    [29, 432],
    [30, 480],
    [31, 498],
    [32, 540],
    [33, 564],
    [34, 612],
    [35, 632],
    [36, 684],
    [37, 706],
    [38, 756],
    [39, 782],
    [40, 840],
    [41, 864],
    [42, 924],
    [44, 1008],
    [46, 1104],
    [48, 1200],
    [50, 1296],
    [52, 1404],
    [54, 1512],
    [56, 1620],
    [58, 1740],
    [60, 1860],
    [62, 1980],
    [64, 2112],
    [66, 2244],
    [68, 2376],
    [70, 2520],
    [72, 2664],
    [74, 2808],
    [76, 2964],
    [78, 3120],
    [80, 3276],
    [82, 3444],
    [84, 3612],
    [86, 3780],
    [88, 3960],
    [90, 4140],
    [94, 4512],
    [98, 4896],
    [100, 5100],
    [114, 6612],
    [124, 7812],
    [200, 21000], # non-optimal
    [500, 130000], # non-optimal
    [1000, 520000], # non-optimal
    [1000, 1002000], # non-optimal
]
"""list: Table of S2 spherical designs from M. Gräf."""

SO3_GAUSS_METHODS = [
    # degree, size
    [1, 4],
    [2, 11],
    [3, 23],
    [4, 43],
    [5, 60], #TODO: weights are not normalized for this one
    [6, 116],
    [7, 168],
    [8, 240],
    [9, 300],
    [11, 504],
    [14, 960],
]
"""list: Table of SO3 Gauss-type quadratures from M. Gräf."""

SO3_CHEBYSHEV_METHODS = [
    # degree, size
    [1, 4],
    [2, 12],
    [3, 24],
    [4, 57],
    [5, 60],
    [6, 154],
    [7, 168],
    [8, 312],
    [9, 360],
    [11, 672],
    [13, 1176],
    [14, 1260],
    [15, 1776],
    [17, 2520],
    [19, 3300],
    [23, 5880],
]
"""list: Table of SO3 Chebyshev-type quadratures from M. Gräf."""
