# https://biopt.github.io/bolib/
#module BOLIB

using BilevelOptSolver

bolib_dir = "dataset/converted"
bolib_files = filter(contains(r".jl$"), readdir(bolib_dir))
include("../" .* bolib_dir * "/problems_list.jl")
include.("../" .* filter(contains(r".jl$"), readdir(bolib_dir; join=true)))

for prob in problems
    if "DempeDutta2012Ex31" == prob
        p = getfield(Main.BOLIB, Symbol(prob))()

        bop = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g)
        sol, is_converged = solve_bop(bop; x_init=p.xy_init)

        if !is_converged
            @info "$prob: failed to converge"
        end

        if (p.Ff_optimal[3] == 1)
            if !isapprox([p.F(sol); p.f(sol)], p.Ff_optimal[1:2])
                @info "$prob: not optimal"
            end
        elseif (p.Ff_optimal[3] == 2)
            if any([p.F(sol); p.f(sol)] .> p.Ff_optimal[1:2])
                @info "$prob: worse than best known"
            end
        end
    end
end
