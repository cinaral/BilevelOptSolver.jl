""" 
[ka_2014_01](https://basblsolver.github.io/BASBLib/LP-NLP/ka_2014_01)

nonconvex f
note this is mb_2007_13 with swapped cubic term in f 
"""
function ka_2014_01()
    n1::Int64 = 1
    n2::Int64 = 1

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        x[1] - y[1]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            x[1] + 1.0
            1.0 - x[1]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        (x[1] * y[1]^2) / 2 - y[1]^3 * x[1]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] + 1.0
            1.0 - y[1]
        ]
    end

    xy_optimal = [0.0; 1.0]
    Ff_optimal = [-1.0; 0.0]
    name = "ka_2014_01"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end
