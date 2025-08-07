""" 
[ka_2014_02](https://basblsolver.github.io/BASBLib/NLP-NLP/ka_2014_02)
"""
function ka_2014_02()
    n1::Int64 = 5
    n2::Int64 = 5

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        -x[1]^2 - x[2]^2 - x[3]^2 - x[4]^2 - x[5]^2 - y[1]^2 - y[2]^2 - y[3]^2 - y[4]^2 - y[5]^2
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [x[1] + 1; x[2] + 1; x[3] + 1; x[4] + 1; x[5] + 1; 1 - x[1]; 1 - x[2]; 1 - x[3]; 1 - x[4]; 1 - x[5]; x[1] - y[1] * y[2]; exp(x[2]) - y[3] - x[1]; -x[2] * y[1]^2]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        y[3] / 10 + y[2]^2 * (x[1] + x[2]) + y[1]^3 + x[3] * x[4] * x[5] * (y[4]^2 + y[5]^2)
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [y[1] + 1; y[2] + 1; y[3] + 1; y[4] + 1; y[5] + 1; 1 - y[1]; 1 - y[2]; 1 - y[3]; 1 - y[4]; 1 - y[5]; y[3]^2 - x[1] + 1 / 5]
    end

    xy_optimal = [1.0; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0]
    Ff_optimal = [-10.0; -3.1]
    name = "ka_2014_02"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end

