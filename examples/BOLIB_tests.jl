# https://biopt.github.io/bolib/
module BOLIB

using BilevelOptSolver
using Test


bolib_dir = "dataset/converted"
bolib_files = filter(contains(r".jl$"), readdir(bolib_dir))
include("../" .* bolib_dir * "/problems_list.jl")
include.("../" .* filter(contains(r".jl$"), readdir(bolib_dir; join=true)))

for prob in problems
    @info prob
    p = getfield(Main.BOLIB, Symbol(prob))()

	#if "AiyoshiShimizu1984Ex2" == prob
		try
			bop = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g)
			sol, is_converged = solve_bop(bop; x_init=p.xy_init)
	
			if !is_converged
				@warn "failed to converge" is_converged
			end
	
			if (p.Ff_optimal[3] == 1)
				if !isapprox([p.F(sol); p.f(sol)], p.Ff_optimal[1:2])
					@warn "not optimal"
				end
			elseif (p.Ff_optimal[3] == 2)
				if [p.F(sol); p.f(sol)] .> p.Ff_optimal[1:2]
					@warn "worse than known"
				end
			end
		catch e
			@warn "failed to run"
			@info e
		end
	#end
   
end # module BOLIB_tests
end

#@testitem "simple_1" begin
#module MitsosBarton2006Ex328

#@test sol.obj_val ≈ 17.014017145179164 atol = 1e-5
#end