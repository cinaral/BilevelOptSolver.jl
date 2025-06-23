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

bop = construct_bop(n₁, n₂, F, G, f, g);
x, is_success, iter_count = solve_bop(bop; x_init, verbosity=2)

if is_success
    Ff = [bop.F(x); bop.f(x)]
    @info "success x = $x, Ff val = $Ff"
end

