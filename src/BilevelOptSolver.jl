module BilevelOptSolver

using Symbolics
using SparseArrays
using LinearAlgebra
using Ipopt
using PATHSolver
using HiGHS
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
