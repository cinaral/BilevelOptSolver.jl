using BilevelOptSolver

# AiyoshiShimizu1984Ex2
n₁::Int64 = 2
n₂::Int64 = 2

function F(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    2 * x₁[1] + 2 * x₁[2] - 3 * x₂[1] - 3 * x₂[2] - 60
end

function G(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    [
        2 * x₂[2] - x₁[2] - x₂[1] - x₁[1] + 40
        50 - x₁[1];
        50 - x₁[2]
        x₁[1]
        x₁[2]
    ]
end

function f(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    (x₂[1] - x₁[1] + 20)^2 + (x₂[2] - x₁[2] + 20)^2
end

function g(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    [
        x₁[1] - 2 * x₂[1] - 10
        x₁[2] - 2 * x₂[2] - 10
        x₂[1] + 10
        x₂[2] + 10
        20 - x₂[1]
        20 - x₂[2]
    ]
end

x_init = [10.0; 10; 20; 20]
Ff_optimal = [5.0; 0]

bop = construct_bop(n₁, n₂, F, G, f, g);
x, is_success, iter_count = solve_bop(bop; x_init, verbosity=6)

if is_success
    Ff = [bop.F(x); bop.f(x)]
    @info "success x = $x, Ff val = $Ff"
    @assert isapprox(Ff_optimal, Ff; rtol=1e-4) # currently not optimal 2025-06-20
end


