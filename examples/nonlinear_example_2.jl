using BilevelOptSolver

# NieEtal2017Ex52
n₁::Int64 = 2
n₂::Int64 = 3
n::Int64 = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1] * y[1] + x[2] * y[2] + x[1] * x[2] * y[1] * y[2] * y[3]
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        x[1] + 1
        x[2] + 1
        -x[1] + 1
        -x[2] + 1
        -y[1] * y[2] + x[1]^2
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    x[1] * y[1]^2 + x[2]^2 * y[2] * y[3] - y[1] * y[3]^2
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n]
    [
        -1 + y[1]^2 + y[2]^2 + y[3]^2
        -y[1]^2 - y[2]^2 - y[3]^2 + 2
    ]
end

x_init =  [0; 0; 1.1097; 0.3143; −0.8184] # causes a bug
#global_opt_sol = [−1; −1; 1.1097; 0.3143; −0.8184]

bop = construct_bop(n₁, n₂, F, G, f, g);
sol, is_success, iter_count = solve_bop(bop; x_init, verbosity=0)
if is_success
    @info "success" sol
end
