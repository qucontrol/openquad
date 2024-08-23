# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

"""
Lebedev quadrature grids for integration on the S2 unit sphere.

The files in this directory have been converted from the data from the
`sphericalquadpy` package from Thomas Camminday [1], revision 0646547 from May
29, 2019, published with the MIT license.

Ref. [1] uses the quadrature points and weights originally provided by John
Burkhardt [2] under the GNU LGPL License.

[1]: https://github.com/thomascamminady/sphericalquadpy
[2]: https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
"""

LEBEDEV_LAIKOV_QUADRATURES = [
    [3, 6],   
    [5, 14],  
    [7, 26],  
    [9, 38],  
    [11, 50],  
    [13, 74],  
    [15, 86],  
    [17, 110], 
    [19, 146], 
    [21, 170], 
    [23, 194], 
    [25, 230], 
    [27, 266], 
    [29, 302], 
    [31, 350], 
    [35, 434], 
    [41, 590], 
    [47, 770], 
    [53, 974], 
    [59, 1202],
    [65, 1454],
    [71, 1730],
    [77, 2030],
    [83, 2354],
    [89, 2702],
    [95, 3074],
    [101, 3470],
    [107, 3890],
    [113, 4334],
    [119, 4802],
    [125, 5294],
    [131, 5810],
]
"""list: Table of Lebedev-Laikov quadrature rules."""
