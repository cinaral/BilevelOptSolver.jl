if !haskey(ENV, "BOLIB_PATH")
    error("You need to obtain [BOLIB.jl](https://github.com/cinaral/BOLIB.jl) (converted from [BOLIBver2](https://biopt.github.io/bolib/)) and set BOLIB_PATH environment variable to run this benchmark.\n")
end
Pkg.develop(path=ENV["BOLIB_PATH"])
import BOLIB
using BilevelOptSolver
using DataFrames
using Statistics
using Random
using Dates

"""
Usage:
```
julia> include("benchmarks/BOLIB_benchmark.jl")
julia> df = benchmark_BOLIB(example_ids=1:165);
```
"""
function benchmark_BOLIB(; example_ids=1:length(BOLIB.examples), verbosity=0, tol=1e-7, init_solver="IPOPT", solver="IPOPT", max_iter=50, conv_dv_len=3, do_force_hp_init=false, do_require_strict_min=false, do_check_x_agreem=true, max_rand_restart_ct=50, rng=MersenneTwister(), do_force_dry_run=false)
    dataframes = []
    success_arr = Bool[]
    elapsed_arr = Float64[]
    prob_count = 0

    for id in example_ids
        name = BOLIB.examples[id]
        print("$id $name\n")

        if "SinhaMaloDeb2014TP9" == name || "SinhaMaloDeb2014TP10" == name
            # these fail to compile for some reason so we skip
            continue
        end

        prob_count += 1
        prob = getfield(Main.BOLIB, Symbol(name))()
        bop, _ = construct_bop(prob.n1, prob.n2, prob.F, prob.G, prob.f, prob.g; verbosity=0)
        x_init = gen_x_init(name, rng)

        # dry run for @elapsed...
        if do_force_dry_run
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init=prob.xy_init, verbosity=0, tol, init_solver, solver, max_iter=2, conv_dv_len, do_force_hp_init, do_require_strict_min, do_check_x_agreem, max_rand_restart_ct)
        end

        elapsed_time = @elapsed begin
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; x_init, verbosity, tol, init_solver, solver, max_iter, conv_dv_len, do_force_hp_init, do_require_strict_min, do_check_x_agreem, max_rand_restart_ct)
        end

        Ff = [bop.F(x); bop.f(x)]
        success, rating, x_optimal = rate_BOLIB_result(name, x, Ff, is_sol_valid; tol)

        if is_sol_valid
            info_status = "success"
        else
            info_status = "FAIL"
        end
        info = "$info_status: $rating ($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s), x = $(round.(x, sigdigits=5)), Ff = $(round.(Ff, sigdigits=5)) (x* = $(round.(x_optimal, sigdigits=5)), Ff* = $(round.(prob.Ff_optimal, sigdigits=5)))\n"

        if verbosity > 0
            print(info)
        end

        dataframes = [dataframes; DataFrame("name" => name, "n1" => bop.n1, "n2" => bop.n2, "m1" => bop.m1, "m2" => bop.m2, "Iterations" => iter_count, "Status" => status, "Solve time (s)" => round.(elapsed_time, sigdigits=3), "Success" => success, "is_sol_valid" => is_sol_valid, "Rating" => rating, "Ff" => Ref(round.(Ff, sigdigits=3)), "Ff*" => Ref(round.(prob.Ff_optimal, sigdigits=3)), "x" => Ref(round.(x, sigdigits=3)), "x*" => Ref(round.(x_optimal, sigdigits=3)), "x_init" => Ref(round.(x_init, sigdigits=3)))]
        push!(success_arr, success)
        push!(elapsed_arr, elapsed_time)
    end

    success_elapsed_arr = elapsed_arr[success_arr]
    success_count = length(findall(success_arr))

    print("Out of $(prob_count) problems, $(success_count) ($(round((success_count/prob_count*100),sigdigits=3))%) were successful.\n")
    print("Elapsed min-max: $(round(minimum(elapsed_arr),sigdigits=2))-$(round(maximum(elapsed_arr),sigdigits=2)) s, median: $(round(median(elapsed_arr),sigdigits=2)) s, mean: $(round(mean(elapsed_arr),sigdigits=2)) s\n")
    if !isempty(success_elapsed_arr)
        print("Elapsed successful min-max: $(round(minimum(success_elapsed_arr),sigdigits=2))-$(round(maximum(success_elapsed_arr),sigdigits=2)) s, median: $(round(median(success_elapsed_arr),sigdigits=2)) s, mean: $(round(mean(success_elapsed_arr),sigdigits=2))\n")
    end

    df = vcat(dataframes...)
    date_time = Dates.format(now(), "yyyy-mm-dd_HHMM")
    df = hcat(df,
        DataFrame("Date Time:" => [
            date_time; fill("", length(df[!, 1]) - 1)]),
        DataFrame("Success Count:" => [
            success_count; fill("", length(df[!, 1]) - 1)]),
        DataFrame("All elapsed median (s):" => [median(elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("All elapsed mean (s):" => [mean(elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("All elapsed min (s):" => [minimum(elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("All elapsed max (s):" => [maximum(elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("Success elapsed median (s):" => [median(success_elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("Success elapsed mean (s):" => [mean(success_elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("Success elapsed min (s):" => [minimum(success_elapsed_arr); fill("", length(df[!, 1]) - 1)]),
        DataFrame("Success elapsed max (s):" => [maximum(success_elapsed_arr); fill("", length(df[!, 1]) - 1)])
    )
    df
end

function rate_BOLIB_result(name, x, Ff, is_sol_valid; tol)
    prob = getfield(Main.BOLIB, Symbol(name))()
    success = false
    rating = ""
    x_optimal = []

    if prob.Ff_optimal[3] == 0 # star means optimal
        rating = "no reference"
    else
        is_cost_optimal = isapprox(Ff, prob.Ff_optimal[1:2]; atol=tol) # looser cost tol
        if prob.Ff_optimal[3] == 1 # star means optimal
            if is_cost_optimal
                rating = "optimal"
            else
                rating = "NOT optimal"
            end
        elseif prob.Ff_optimal[3] == 2 # star means best
            if is_cost_optimal
                rating = "optimal (known)"
            else
                rating = "NOT optimal (known)"
            end
        end
    end

    # if optimal we consider a success 
    if rating == "optimal"
        success = true

        if !is_sol_valid
            @info "couldn't verify valid sol but the result is optimal"
        end
    end

    if name == "MitsosBarton2006Ex31"
        x_optimal = [1.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex32" # no sol
        if !is_sol_valid
            success = true
        end
    elseif name == "ClarkWesterberg1988"
        x_optimal = [19.0; 14.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "LiuHart1994"
        x_optimal = [4.0; 4.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Bard1984a"
        x_optimal = [0.8889; 2.222]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "AnandalinghamWhite1990"
        x_optimal = [16.0; 11.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Bard1991Ex2"
        x_optimal = [0.0; 0.0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "ClarkWesterberg1990b"
        x_optimal = [5.0; 4; 2]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "BardFalk1982Ex2"
        x_optimal = [2.0; 0; 1.5; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "CandlerTownsley1982"
        x_optimal = [0.0; 0.9; 0; 0.6; 0.4]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "LucchettiEtal1987"
        x_optimal = [1.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Yezza1996Ex31"
        x_optimal = [0.25; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Dempe1992b"
        x_optimal = [1.0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "ClarkWesterberg1990a" # multiple sols
        x_optimal = [1.0; 3]
        if is_sol_valid && (isapprox(x, x_optimal; atol=tol) || isapprox(x, [3.0; 5]; atol=tol) || isapprox(x, [4.4; 4.8]; atol=tol))
            success = true
        end
    elseif name == "TuyEtal2007" # 2025-08-27 TODO check
        x_optimal = [1.5; 1.5]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Bard1988Ex1"
        x_optimal = [1.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "ShimizuAiyoshi1981Ex1"
        x_optimal = [10.0; 10.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "DeSilva1978"
        x_optimal = [0.5; 0.5; 0.5; 0.5]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "FalkLiu1995"
        x_optimal = [0.866; 0.866; 0.866; 0.866]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "ShimizuAiyoshi1981Ex2" # multiple sols
        x_optimal = [0.0; 0; -10; -10]
        if is_sol_valid && (isapprox(x, x_optimal; atol=tol) || isapprox(x, [20.0; 5; 10; 5]; atol=tol))
            success = true
        end
    elseif name == "Bard1988Ex3"
        x_optimal = [0.0; 2; 1.875; 0.896]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "DempeDutta2012Ex31"
        x_optimal = [0.707; 0.707; 0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Bard1988Ex2"
        x_optimal = [7.0; 3; 12; 18; 0; 10; 30; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex35"
        x_optimal = [0.5]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex36"
        x_optimal = [-1.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex310" # multiple sols
        x_optimal = [0.1; 0.5]
        if is_sol_valid && 0.1 - tol < x[1] ≤ 1.0 + tol && isapprox(x[2], 0.5; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex313"
        x_optimal = [0.0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex315"
        x_optimal = [-1.0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex316" # multiple sols
        x_optimal = [0.0; 1]
        if is_sol_valid && (isapprox(x, x_optimal; atol=tol) || isapprox(x, [-0.5; -1]; atol=tol))
            success = true
        end
    elseif name == "KleniatiAdjiman2014Ex3"
        x_optimal = [0.0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex39"
        x_optimal = [-1.0; -1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "NieEtal2017Ex34"
        x_optimal = [2.0; 0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "GumusFloudas2001Ex5"
        x_optimal = [0.1936; 9.967; 10]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "GumusFloudas2001Ex3"
        x_optimal = [0.0; 0.9; 0; 0.6; 0.4]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex312"
        x_optimal = [0.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex314"
        x_optimal = [0.25; 0.5]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex317" # multiple sols
        x_optimal = [-0.25; 0.5]
        if is_sol_valid && (isapprox(x, x_optimal; atol=tol) || isapprox(x, [-0.25; -0.5]; atol=tol))
            success = true
        end
    elseif name == "MitsosBarton2006Ex318" # 2025-08-27 TODO: check this
        x_optimal = [1.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex319" # 2025-08-27 TODO: check this
        x_optimal = [0.189; 0.4343]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex320"
        x_optimal = [0.5; 0.5]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex321"
        x_optimal = [-0.5545; 0.4554]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex324"
        x_optimal = [0.2106; 1.799]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex322"
        x_optimal = [-0.5545; 0.4554]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "DempeDutta2012Ex24" # 2025-08-27 TODO: check this
        x_optimal = [1.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "ShimizuEtal1997b"
        x_optimal = [11.25; 5]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "YeZhu2010Ex42"
        x_optimal = [1.0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex38"
        x_optimal = [-0.567; 0.0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Colson2002BIPA2"
        x_optimal = [1.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Colson2002BIPA4"
        x_optimal = [0; 0.579]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Colson2002BIPA1"
        x_optimal = [6.082; 4.487]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Colson2002BIPA3"
        x_optimal = [4.0; 0]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "Colson2002BIPA5"
        x_optimal = [1.941; 0.0; 1.211]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "FloudasZlobec1998"
        x_optimal = [1.0; 0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "MitsosBarton2006Ex326" # multiple sols
        x_optimal = [-1.0; -1; 1; 1; -0.707]
        if is_sol_valid && (isapprox(x, x_optimal; atol=tol) || isapprox(x, [-1.0; -1; 1; -1; -0.707]; atol=tol))
            success = true
        end
    elseif name == "NieWangYe2017Ex52"
        x_optimal = [-1.0; -1.0; 1.097; 0.3143; -0.8284]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "NieWangYe2017Ex57"
        x_optimal = [1.0; 1; 0; 0; 1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "NieWangYe2017Ex54"
        x_optimal = [0.0; 0; -0.707; -0.707; 0.618; 0; -0.558; -0.554]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "NieWangYe2017Ex58"
        x_optimal = [0.544; 0.468; 0.490; 0.495; -0.783; -0.501; -0.288; -0.184]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    elseif name == "KleniatiAdjiman2014Ex4"
        x_optimal = [1.0; -1; -1; -1; -1; -1; -1; -1; -1; -1]
        if is_sol_valid && isapprox(x, x_optimal; atol=tol)
            success = true
        end
    else
        is_sol_valid # assuming it's local sol
        success = true
    end

    (; success, rating, x_optimal)
end

function gen_x_init(name, rng)
    prob = getfield(Main.BOLIB, Symbol(name))()
    x_init = prob.xy_init
    x_init
end