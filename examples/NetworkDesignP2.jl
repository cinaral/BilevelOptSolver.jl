using BilevelOptSolver

n1::Int64 = 5
n2::Int64 = 5

function F(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    x[1] + x[2] + x[3] + x[4] + x[5] + y[3] * (y[3] / (x[3] + 1) + 10) + y[1] * (y[1] / (x[1] + 1) + 50) + y[5] * (y[5] / (x[5] + 1) + 50) + (10 * y[2]^2) / (x[2] + 1) + (10 * y[4]^2) / (x[4] + 1)
end

function G(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [x[1] + 1; x[2] + 1; x[3] + 1; x[4] + 1; x[5] + 1]
end

function f(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    y[3] * (y[3] / (2 * (x[3] + 1)) + 10) + y[1] * (y[1] / (2 * (x[1] + 1)) + 50) + y[5] * (y[5] / (2 * (x[5] + 1)) + 50) + (5 * y[2]^2) / (x[2] + 1) + (5 * y[4]^2) / (x[4] + 1)
end

function g(xy)
    x = @view xy[1:n1]
    y = @view xy[n1+1:n1+n2]
    [6 - y[3] - y[5] - y[1]; y[3] - y[2] + y[5]; y[1] + y[3] - y[4]; y[1] + y[3] + y[5] - 6; y[2] - y[3] - y[5]; y[4] - y[3] - y[1]; y[1]; y[2]; y[3]; y[4]; y[5]]
end

xy_init = Float64[0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
Ff_optimal = Float64[142.9; 81.95; 2]

bop = construct_bop(n1, n2, F, G, f, g);
#sol, is_success, iter_count = @btime solve_bop(bop; verbosity=0)
sol, is_success, iter_count =  solve_bop(bop; x_init=xy_init, is_using_PATH=false, verbosity=2)
if is_success
    @info "success" sol
end
