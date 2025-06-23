
if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) and set BOLIB_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB

function run_all_BOLIB_examples(; verbosity=0)
    prob_count = 0
    converged_count = 0
    optimalish_count = 0
    suboptimalish_count = 0
    bops = []
    sols = Vector{Float64}[]
    success_arr = Bool[]
    iter_counts = Int64[]
    elapsed_arr = Float64[]

    for prob in BOLIB.examples
        if "SinhaMaloDeb2014TP9" == prob || "SinhaMaloDeb2014TP10" == prob
            # these fail to compile for some reason so we skip
            continue
        end
        p = getfield(Main.BOLIB, Symbol(prob))()

        bop = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g, verbosity=0)
        # dry runs for @time...
        solve_bop(bop; max_iter=1, x_init=p.xy_init, verbosity=0, is_using_PATH=false)
        #solve_bop(bop; max_iter=1, x_init=p.xy_init, verbosity=0, is_using_PATH=true)
        elapsed_time = @elapsed begin
            x, is_success, iter_count = solve_bop(bop; max_iter=200, x_init=p.xy_init, verbosity, is_using_PATH=false)
        end
        push!(bops, bop)
        push!(sols, x)
        push!(iter_counts, iter_count)
        push!(success_arr, is_success)
        push!(elapsed_arr, elapsed_time)

        if is_success
            converged_count += 1
            print("$(prob_count+1)\t $prob\t $(iter_count) iterations:\t ")
         
            is_optimal, is_best = rate_BOLIB_success(p, bop, x, elapsed_time)
            if is_optimal || is_best
                optimalish_count += 1
            else
                suboptimalish_count += 1
            end
        else
            print("$(prob_count+1)\t $prob\tFailed to converge")
        end
		print("\t($elapsed_time s)\n")
        prob_count += 1
        #end
    end

    success_elapsed_sum = 0
    n_success = 0
    for (i, is_success) in enumerate(success_arr)
        if is_success
            success_elapsed_sum += elapsed_arr[i]
        end
        n_success += 1
    end

    print("Out of $prob_count problems, $converged_count ($(converged_count/prob_count*100)%) converged.\n")
    print("Out of converged solutions: $optimalish_count ($(optimalish_count/converged_count*100)%) were optimal or best known, while $suboptimalish_count ($(suboptimalish_count/converged_count*100)%) were suboptimal or worse than best known.\n")
    print("Elapsed min/max (s): $(minimum(elapsed_arr))/$(maximum(elapsed_arr)), success mean elapsed (s): $(success_elapsed_sum/n_success)\n")
end

function rate_BOLIB_success(BOLIB, bop, x, elapsed_time; tol=1e-6)
    is_optimal = false
    is_best = false
    Ff = [bop.F(x); bop.f(x)]
    Ff_star = BOLIB.Ff_optimal[1:2]

    if BOLIB.Ff_optimal[3] == 1 # star means optimal
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_optimal = true
            print("optimal")
        else
            print("SUBOPTIMAL:  Ff = $Ff \tFf*: $Ff_star")
        end
    elseif BOLIB.Ff_optimal[3] == 2 # star means best
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_best = true
            print("same as best known")
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_best = true
            print("better than best known: Ff = $Ff \tFf* = $Ff_star")
        else
            print("WORSE than best known: Ff = $Ff \tFf*: $Ff_star")
        end
    else
        print("no reference solution")
    end

    (; is_optimal, is_best)
end
