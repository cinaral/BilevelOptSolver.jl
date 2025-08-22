module BilevelOptSolver

using Symbolics
using SparseArrays
using LinearAlgebra
using Ipopt
using PATHSolver
using Random
import Pkg

# OPTIONAL: Uncomment the following if you wish to use STFC LibHSL, which you must obtain separately 
#if haskey(ENV, "HSL_PATH")
#    Pkg.develop(path=ENV["HSL_PATH"])
#end
# import HSL_jll

include("solver_interfaces.jl")
include("construct_bop.jl")
include("solve_bop.jl")

export construct_bop, solve_bop
end # module BilevelOptSolver