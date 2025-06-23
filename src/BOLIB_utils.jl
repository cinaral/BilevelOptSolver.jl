
if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) and set BOLIB_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB

function run_all_BOLIB_examples(; verbosity=0)
    prob_count = 0
    success_count = 0
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
        push!(elapsed_arr, elapsed_time)

        if prob_count % 20 == 0
            print("id name iterations: elapsed: rating, x, Ff (Ff*), result\n")
        end
        print("$(prob_count+1)\t $prob\t $(iter_count):\t ")

        is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(p, bop, x)
        print("$(round(elapsed_time, sigdigits=5)) s,\t" * rating * ",\t$(round.(x, sigdigits=5)), $(round.(Ff, sigdigits=5)) ($(round.(Ff_star, sigdigits=5)))\t")

        if is_optimal || is_best
            optimalish_count += 1
            if !is_success
                print("Success but didn't converge!\t")
                is_success = true
            end
        else
            if is_success # otherwise doesn't matter
                suboptimalish_count += 1
            end
        end

        if is_success
            print("SUCCESS\n")
            success_count += 1
        else
            print("FAIL\n")
        end

        push!(success_arr, is_success)
        prob_count += 1
    end

    success_elapsed_sum = 0
    n_success = 0
    for (i, is_success) in enumerate(success_arr)
        if is_success
            success_elapsed_sum += elapsed_arr[i]
        end
        n_success += 1
    end

    print("Out of $prob_count problems, $success_count ($(success_count/prob_count*100)%) converged.\n")
    print("Out of converged solutions: $optimalish_count ($(optimalish_count/success_count*100)%) were optimal or best known, while $suboptimalish_count ($(suboptimalish_count/success_count*100)%) were suboptimal or worse than best known.\n")
    print("Elapsed min/max (s): $(minimum(elapsed_arr))/$(maximum(elapsed_arr)), success mean elapsed (s): $(success_elapsed_sum/n_success)\n")
end

function rate_BOLIB_result(BOLIB, bop, x; tol=1e-2)
    is_optimal = false
    is_best = false
    Ff = [bop.F(x); bop.f(x)]
    Ff_star = BOLIB.Ff_optimal[1:2]
    rating = ""

    if BOLIB.Ff_optimal[3] == 1 # star means optimal
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_optimal = true
            rating = "optimal"
        elseif isapprox(Ff[1], Ff_star[1]; atol=2 * tol)
            is_optimal = true
            if Ff[2] < Ff_star[2] - tol
                rating = "optimal F and f is somehow BETTER"
            else
                rating = "optimal F but f is worse"
            end
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_optimal = true
            rating = "better than optimal?!"
        else
            rating = "SUBOPTIMAL"
        end
    elseif BOLIB.Ff_optimal[3] == 2 # star means best
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_best = true
            rating = "best known"
        elseif isapprox(Ff[1], Ff_star[1]; atol=2 * tol)
            is_best = true
            if Ff[2] < Ff_star[2] - tol
                rating = "best known F and f is somehow BETTER"
            else
                rating = "best known F but f is worse"
            end
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_best = true
            rating = "better than best known"
        else
            rating = "WORSE than best known"
        end
    else
        rating = "no reference solution"
    end
    is_optimal, is_best, Ff, Ff_star, rating
end
