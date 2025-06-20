using BilevelOptSolver

# Bard1988Ex2
n₁::Int64 = 4
n₂::Int64 = 4

function F(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    (x₂[2] + x₂[4]) * (x₂[2] + x₂[4] - 160) + (x₂[1] + x₂[3]) * (x₂[1] + x₂[3] - 200)
end

function G(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    [
        40 - x₁[2] - x₁[3] - x₁[4] - x₁[1]
        10 - x₁[1]
        5 - x₁[2]
        15 - x₁[3]
        20 - x₁[4]
        x₁[1]
        x₁[2]
        x₁[3]
        x₁[4]
    ]
end

function f(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    (x₂[1] - 4)^2 + (x₂[4] - 2)^2 + (x₂[2] - 13)^2 + (x₂[3] - 35)^2
end

function g(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    [
        x₁[1] - (2 * x₂[1]) / 5 - (7 * x₂[2]) / 10; x₁[2] - (3 * x₂[1]) / 5 - (3 * x₂[2]) / 10;
        x₁[3] - (2 * x₂[3]) / 5 - (7 * x₂[4]) / 10; x₁[4] - (3 * x₂[3]) / 5 - (3 * x₂[4]) / 10;
        20 - x₂[1]
        20 - x₂[2]
        40 - x₂[3]
        40 - x₂[4]
        x₂[1]
        x₂[2]
        x₂[3]
        x₂[4]
    ]
end

x_init = [5.0; 5; 15; 15; 0; 0; 0; 0]
Ff_optimal = [-6600.0; 54]

bop = construct_bop(n₁, n₂, F, G, f, g);
x, is_success, iter_count = solve_bop(bop; x_init, verbosity=2)

if is_success
    Ff = [bop.F(x); bop.f(x)]
    @info "success x = $x, Ff val = $Ff"
    @assert isapprox(Ff_optimal, Ff; rtol=1e-4)
end