
if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) and set BOLIB_PATH environment variable to run this example.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB

using DataFrames

function run_all_BOLIB_examples(; verbosity=0, max_iter=100, is_checking_x_agree=false, init_solver="PATH", solver="PATH", tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=false)
    prob_count = 0
    success_count = 0
    optimalish_count = 0
    suboptimalish_count = 0
    ps = []
    bops = []
    prob_names = []
    x_out_arr = Vector{Float64}[]
    status_arr = []
    success_arr = Bool[]
    iter_counts = Int64[]
    elapsed_arr = Float64[]
    ratings = []
    n1_arr = Int64[]
    n2_arr = Int64[]
    m1_arr = Int64[]
    m2_arr = Int64[]
    x_init_arr = Vector{Float64}[]
    Ff_star_arr = Vector{Float64}[]
    Ff_out_arr = Vector{Float64}[]

    for prob in BOLIB.examples
        if "SinhaMaloDeb2014TP9" == prob || "SinhaMaloDeb2014TP10" == prob
            # these fail to compile for some reason so we skip
            continue
        end
        p = getfield(Main.BOLIB, Symbol(prob))()

        bop, _ = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g, verbosity=0)

        # print iter info
        if prob_count % 20 == 0
            print("id name: success (status), iters (elapsed s): success, rating, x -> Ff (Ff*), result\n")
        end
        print("$(prob_count+1)\t $prob:\t ")

        # dry run for @elapsed...
        solve_bop(bop; max_iter=1, x_init=p.xy_init, verbosity=0)

        elapsed_time = @elapsed begin
            success, x, iter_count, status = solve_bop(bop; max_iter, x_init=p.xy_init, verbosity, is_checking_x_agree, tol, norm_dv_len, conv_dv_len, is_checking_min, init_solver, solver)
        end

        is_optimal, is_best, Ff, Ff_star, rating, is_probably_wrong = rate_BOLIB_result(p, bop, x)

        if success
            if is_optimal || is_best
                optimalish_count += 1
            else
                suboptimalish_count += 1
            end
        end

        if success
            print("success ")
            if is_probably_wrong
                print("(!) ")
            end
            success_count += 1
        else
            print("FAIL ")
        end
        print("($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s):\t" * rating * ",\t $(round.(x, sigdigits=5)) -> $(round.(Ff, sigdigits=5)) ($(round.(Ff_star, sigdigits=5)))\n")
        push!(ps, p)
        push!(bops, bop)
        push!(prob_names, prob)
        push!(n1_arr, bop.n1)
        push!(n2_arr, bop.n2)
        push!(m1_arr, bop.m1)
        push!(m2_arr, bop.m2)
        push!(x_init_arr, p.xy_init)
        push!(iter_counts, iter_count)
        push!(elapsed_arr, elapsed_time)
        push!(status_arr, status)
        push!(x_out_arr, x)
        push!(Ff_out_arr, Ff)
        push!(Ff_star_arr, Ff_star)
        push!(ratings, rating)
        push!(success_arr, success)

        prob_count += 1
    end

    df = DataFrame("name" => prob_names, "n1" => n1_arr, "n2" => n2_arr, "m1" => m1_arr, "m2" => m2_arr, "Status" => status_arr, "Success" => success_arr, "Rating" => ratings, "Iterations" => iter_counts, "Solve time (s)" => round.(elapsed_arr, sigdigits=3), "x" => map((x) -> round.(x, sigdigits=3), x_out_arr), "Ff" => map((Ff) -> round.(Ff, sigdigits=3), Ff_out_arr), "x_init" => map((x) -> round.(x, sigdigits=3), x_init_arr), "Ff*" => map((Ff) -> round.(Ff, sigdigits=3), Ff_star_arr))

    success_elapsed_arr = elapsed_arr[success_arr]

    print("Out of $(prob_count) problems, $(success_count) ($(round((success_count/prob_count*100),sigdigits=3))%) were successful.\n")
    print("Out of successful solutions: $(optimalish_count) ($(round(optimalish_count/success_count*100,sigdigits=3))%) were optimal or best known, while $(suboptimalish_count) ($(round(suboptimalish_count/success_count*100,sigdigits=3))%) were suboptimal or worse than best known.\n")
    print("Elapsed min-max: $(round(minimum(elapsed_arr),sigdigits=2))-$(round(maximum(elapsed_arr),sigdigits=2)) s, median: $(round(median(elapsed_arr),sigdigits=2)) s, mean: $(round(mean(elapsed_arr),sigdigits=2)) s\n")
    print("Elapsed successful min-max: $(round(minimum(success_elapsed_arr),sigdigits=2))-$(round(maximum(success_elapsed_arr),sigdigits=2)) s, median: $(round(median(success_elapsed_arr),sigdigits=2)) s, mean: $(round(mean(success_elapsed_arr),sigdigits=2))\n")

    (; prob_count, success_arr, success_count, optimalish_count, suboptimalish_count, ps, bops, x_out_arr, iter_counts, elapsed_arr, status_arr, df)
end

function rate_BOLIB_result(BOLIB, bop, x; tol=1e-2)
    is_probably_wrong = false
    is_optimal = false
    is_best = false
    Ff = [bop.F(x); bop.f(x)]
    Ff_star = BOLIB.Ff_optimal[1:2]
    rating = ""

    if BOLIB.Ff_optimal[3] == 1 # star means optimal
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_optimal = true
            rating = "optimal Ff"
        elseif isapprox(Ff[1], Ff_star[1]; atol=2 * tol)
            is_optimal = true
            if Ff[2] < Ff_star[2] - tol
                rating = "optimal F and f is somehow BETTER"
            else
                rating = "optimal F but f is worse"
            end
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_probably_wrong = true
            rating = "better than optimal?!"
        else
            rating = "SUBOPTIMAL Ff"
        end
    elseif BOLIB.Ff_optimal[3] == 2 # star means best
        if isapprox(Ff, Ff_star; atol=2 * tol)
            is_best = true
            rating = "best known Ff"
        elseif isapprox(Ff[1], Ff_star[1]; atol=2 * tol)
            is_best = true
            if Ff[2] < Ff_star[2] - tol
                rating = "best known F and f is somehow BETTER"
            else
                rating = "best known F but f is worse"
            end
        elseif Ff[1] < Ff_star[1] - tol || (isapprox(Ff[1], Ff_star[1]; atol=2 * tol) && Ff[2] < Ff_star[2] - tol)
            is_best = true
            is_probably_wrong = true
            rating = "better than best known Ff"
        else
            rating = "WORSE than best known Ff"
        end
    else
        rating = "no reference solution"
    end
    is_optimal, is_best, Ff, Ff_star, rating, is_probably_wrong
end
