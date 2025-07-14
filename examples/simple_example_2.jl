using BilevelOptSolver

n₁::Int64 = 1
n₂::Int64 = 1

function F(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    x₁[1] + x₂[1]
end

function G(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    [
        -x₁[1] - 2 * x₂[1] + 10
        2 * x₁[1] - x₂[1]
        -x₁[1] + 2 * x₂[1]
        x₂[1]
    ]
end

function f(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    x₂[1]
end

function g(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁+n₂]
    [
        x₁[1] + x₂[1] - 3.0,
    ]
end

bop, _ = construct_bop(n₁, n₂, F, G, f, g);
x, status, iter_count = solve_bop(bop; x_init, is_checking_min=true, verbosity=5)

Ff = [bop.F(x); bop.f(x)]
print("status: [$status], $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5))\n");
