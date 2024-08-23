# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

"""Nearly optimal SO3 coverings from C.F.F. Karney.

The files in this directory have been converted from the `orientation` package
from Charles F. F. Karney [1], revision db92f00 from October 12, 2015,
published with the MIT license.

These coverings have non-uniform weights and no degree of exactness.

Non-optimal coverings have been excluded from this directory.

[1]: https://github.com/cffk/orientation
"""

SO3_COVERINGS = [
    # size, original covering name
    [24, "c48u1"],
    [60, "c600v"],
    [216, "c48n9"],
    [360, "c600vc"],
    [648, "c48u27"],
    [1992, "c48u83"],
    [3768, "c48u157"],
    [4344, "c48u181"],
    [7416, "c48n309"],
    [12456, "c48u519"],
    [12648, "c48n527"],
    [19560, "c48u815"],
    [27672, "c48u1153"],
    [28824, "c48u1201"],
    [39384, "c48u1641"],
    [53256, "c48u2219"],
    [68808, "c48u2867"],
    [70728, "c48u2947"],
    [89592, "c48u3733"],
    [112824, "c48u4701"],
    [113976, "c48u4749"],
    [141096, "c48u5879"],
    [170664, "c48u7111"],
    [207576, "c48u8649"],
]
"""list: Table of SO3 coverings from C.F.F. Karney."""
