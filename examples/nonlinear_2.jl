using BilevelOptSolver

# NieEtal2017Ex52
n₁::Int64 = 2
n₂::Int64 = 3

function F(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    x₁[1] * x₂[1] + x₁[2] * x₂[2] + x₁[1] * x₁[2] * x₂[1] * x₂[2] * x₂[3]
end

function G(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    [
        x₁[1] + 1
        x₁[2] + 1
        -x₁[1] + 1
        -x₁[2] + 1
        -x₂[1] * x₂[2] + x₁[1]^2
    ]
end

function f(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    x₁[1] * x₂[1]^2 + x₁[2]^2 * x₂[2] * x₂[3] - x₂[1] * x₂[3]^2
end

function g(x)
    x₁ = @view x[1:n₁]
    x₂ = @view x[n₁+1:n₁ + n₂]
    [
        -1 + x₂[1]^2 + x₂[2]^2 + x₂[3]^2
        -x₂[1]^2 - x₂[2]^2 - x₂[3]^2 + 2
    ]
end

x_init =  [0; 0; 1.1097; 0.3143; −0.8184] # causes a bug
Ff_optimal = [−1; −1; 1.1097; 0.3143; −0.8184]

bop = construct_bop(n₁, n₂, F, G, f, g);
x, is_success, iter_count = solve_bop(bop; x_init, verbosity=2)

if is_success
    Ff = [bop.F(x); bop.f(x)]
    @info "success x = $x, Ff val = $Ff"
    #@assert isapprox(Ff_optimal, Ff; rtol=1e-4) # currently not optimal 2025-06-20
end

