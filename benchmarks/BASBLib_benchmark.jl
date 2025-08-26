include("../examples/BASBLib_examples.jl")

function rate_BASBLib_result2(name, x, Ff; tol=1e-7)
    prob = getfield(Main.BASBLib, Symbol(name))()
    is_cost_optimal = isapprox(Ff, prob.Ff_optimal; atol=tol)
    is_xy_optimal = isapprox(x, prob.xy_optimal; atol=tol)
    
    if (is_cost_optimal && is_xy_optimal)
        return "optimal"
    else
        return "FUCKNOT optimal"
    end
end
(; info, Ff, is_sol_valid, x, λ, iter_count, status, elapsed_time, bop, syms) = solve_BASBLib_prob(name="as_2013_01", verbosity=5, rate_fun=rate_BASBLib_result2);
print("")