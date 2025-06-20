using BilevelOptSolver

# Kleinert (2021)
n₁::Int64 = 1
n₂::Int64 = 1
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1] + 6y[1]
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        x[1] - 5y[1] + 12.5,
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    -y[1]
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        2x[1] - y[1]
        -x[1] - y[1] + 6
        -x[1] + 6y[1] + 3
        x[1] + 3y[1] - 3
    ]
end
Ff_ref = [3.0, 1.]
x_init = [2.4564338234981746; 0.9845259227566776]
bop = construct_bop(n₁, n₂, F, G, f, g);
sol, is_success, iter_count = solve_bop(bop; x_init, verbosity=0) # bilevel feasible init
if is_success
    @info "success" sol
end
