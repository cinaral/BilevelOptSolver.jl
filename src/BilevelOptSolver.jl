module BilevelOptSolver

using Symbolics
using SparseArrays
using LinearAlgebra
using Ipopt
using PATHSolver
using HiGHS
using Random
import Pkg

if haskey(ENV, "HSL_PATH")
    Pkg.develop(path=ENV["HSL_PATH"])
end
# uncomment this if you have obtained HSL:
import HSL_jll


include("solver_interfaces.jl")
include("construct_bop.jl")
include("solve_bop.jl")

export construct_bop, solve_bop
end # module BilevelOptSolver

# TODO 2025-07-04: 
# add parameters
# fritz john type
# verify stationary point for BOP_i