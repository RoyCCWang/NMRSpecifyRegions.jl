# SPDX-License-Identifier: MPL-2.0
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

module NMRSpecifyRegions

# Write your package code here.

# dependencies to NMRSpectraSimulator. Need to add explicitly since these packages are not on the Julia public registry.
#import GISSMOReader # https://github.com/AI4DBiological-Systems/GISSMOReader.jl

#import NMRSpectraSimulator # https://github.com/AI4DBiological-Systems/NMRSpectraSimulator.jl
import Graphs

include("../src/types.jl")
include("../src/utils.jl")

include("../src/IO/string_helpers.jl")
include("../src/IO/config_helpers.jl")


include("../src/regions/setup.jl")
include("../src/regions/overlap.jl")

end
