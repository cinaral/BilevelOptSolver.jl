# https://biopt.github.io/bolib/

using BilevelOptSolver

bolib_dir = "dataset/converted"
bolib_files = filter(contains(r".jl$"), readdir(bolib_dir))
include("../" .* bolib_dir * "/working_problems_list.jl")
include.("../" .* filter(contains(r".jl$"), readdir(bolib_dir; join=true)))

tol = 1e-3
prob_count = 1
converged_count = 0
optimalish_count = 0
suboptimalish_count = 0
bops = []
sols = Vector{Float64}[]
success_arr = Bool[]
iter_counts = Int64[]
elapsed_arr = Float64[]

for prob in problems
    #if "TuyEtal2007Ex3" == prob
    p = getfield(Main, Symbol(prob))()

    bop = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g, verbosity=0)
    # dry runs for @time...
    solve_bop(bop; max_iter=1, x_init=p.xy_init, verbosity=0, is_using_PATH=false)
    solve_bop(bop; max_iter=1, x_init=p.xy_init, verbosity=0, is_using_PATH=true)
    elapsed_time = @elapsed begin
        sol, is_success, iter_count = solve_bop(bop; max_iter=200, x_init=p.xy_init, verbosity=0, is_using_PATH=true)
    end
    push!(bops, bop)
    push!(sols, sol)
    push!(iter_counts, iter_count)
    push!(success_arr, is_success)
    push!(elapsed_arr, elapsed_time)

    if is_success
        global converged_count += 1
        print("$prob_count\t $prob\t $(iter_count) iterations:\t ")

        if p.Ff_optimal[3] == 1
            if isapprox([p.F(sol); p.f(sol)], p.Ff_optimal[1:2]; rtol=tol)
                print("optimal")
                global optimalish_count += 1
            else
                print("SUBOPTIMAL")
                global suboptimalish_count += 1
            end
        elseif p.Ff_optimal[3] == 2
            if isapprox([p.F(sol); p.f(sol)], p.Ff_optimal[1:2]; rtol=tol)
                print("same as best known")
                global optimalish_count += 1
            elseif p.F(sol) < p.Ff_optimal[1] - tol || (isapprox(p.F(sol), p.Ff_optimal[1]; rtol=tol) && p.f(sol) < p.Ff_optimal[2] - tol)
                print("better than best known: ours F=$(p.F(sol)), f=$(p.F(sol)), best known: F=$(p.Ff_optimal[1]), f=$(p.Ff_optimal[2])")
                global optimalish_count += 1
            else
                print("WORSE than best known")
                global suboptimalish_count += 1
            end
        else
            print("no reference solution")
        end
        #print("\n")
    else
        print("$prob_count\t $prob\tFailed to converge")
    end

    print("\tElapsed: $(elapsed_time) s\n")
    global prob_count += 1
    #end
end
prob_count -= 1

success_elapsed_sum = 0
n_success = 0
for (i, is_success) in enumerate(success_arr)
    if is_success
        global success_elapsed_sum += elapsed_arr[i]
    end
    global n_success += 1
end

print("Out of $prob_count problems, $converged_count ($(converged_count/prob_count*100)%) converged.\nOut of converged solutions: $optimalish_count ($(optimalish_count/converged_count*100)%) were optimal or best known, while $suboptimalish_count ($(suboptimalish_count/converged_count*100)%) were suboptimal or worse than best known.\nElapsed min/max (s): $(minimum(elapsed_arr))/$(maximum(elapsed_arr)), success mean elapsed (s): $(success_elapsed_sum/n_success)")


succ_elapsed_arr = elapsed_arr[success_arr .== true];
#import GLMakie
#fig = GLMakie.Figure()
#ax = GLMakie.Axis(fig[1, 1], xlabel = "t", xscale=log10)
#GLMakie.scatter!(ax, succ_elapsed_arr, zeros(length(succ_elapsed_arr)))
#scatter!(ax, ipopt_elapsed_arr[ipopt_success_arr.==true], ones(length(ipopt_elapsed_arr[ipopt_success_arr.==true])))