using BilevelOptSolver

# AllendeStill2013
n₁::Int64 = 2
n₂::Int64 = 2
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1]^2 - 2 * x[1] + x[2]^2 - 2*x[2] + y[1]^2 + y[2]^2
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
		x[1]
		x[2]
		y[1]
		y[2]
		-x[1] + 2
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
	y[1]^2 - 2*x[1]*y[1] + y[2]^2 - 2*x[2]*y[2]
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
		-(y[1] - 1)^2 + 0.25
		-(y[2] - 1)^2 + 0.25
    ]
end

x_optimal = [0.5 0.5 0.5 0.5]

bop = construct_bop(n₁, n₂, F, G, f, g; verbosity=2);
sol, is_success, iter_count = solve_bop(bop; x_init, verbosity=2)
if is_success
    @info "success" sol
end
