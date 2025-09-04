function RobustPortfolioP2(N; d=2)
    n1::Int64 = N + 1
    n2::Int64 = N

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        -x[end]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            -x[end] + y' * x[1:N];
            x[1:N];
            1.0 - sum(x[1:N]);
            sum(x[1:N]) - 1.0
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        y' * x[1:N] - x[end]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        i = (1:N)'
        si = ((0.05 / 3 / N) * sqrt.(2 * N * (N + 1) .* i))
        yi = 1.15 .+ (0.05 / N) .* i
        sx = 1.5 * (1 + sum((x[1:N] .- 1 / N) .^ 2))
        [
            sx^2 - sum((abs.(y .- yi)) .^ 2 ./ si);
            y
        ]
    end

    xy_init = ones(n1 + n2)
    Ff_optimal = Float64[-1.15; 0; 2]

    (; n1, n2, F, G, f, g, xy_init, Ff_optimal)
end

