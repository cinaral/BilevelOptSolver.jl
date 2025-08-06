""" 
[as_1984_01](https://basblsolver.github.io/BASBLib/LP-QP/as_1984_01)

"""
function as_1984_01()
	n1::Int64 = 2
	n2::Int64 = 2

	function F(xy)
		x = @view xy[1:n1]
		y = @view xy[n1+1:n1+n2]
		2*x[1] + 2*x[2] - 3*y[1] - 3*y[2] - 60
	end

	function G(xy)
		x = @view xy[1:n1]
		y = @view xy[n1+1:n1+n2]
		[2*y[2] - x[2] - y[1] - x[1] + 40; 50 - x[1]; 50 - x[2]; x[1]; x[2]]
	end

	function f(xy)
		x = @view xy[1:n1]
		y = @view xy[n1+1:n1+n2]
		(y[1] - x[1] + 20)^2 + (y[2] - x[2] + 20)^2
	end

	function g(xy)
		x = @view xy[1:n1]
		y = @view xy[n1+1:n1+n2]
		[x[1] - 2*y[1] - 10; x[2] - 2*y[2] - 10; y[1] + 10; y[2] + 10; 20 - y[1]; 20 - y[2]]
	end

	xy_init = Float64[10;10;20;20]
	xy_optimal = Float64[50;50;-10;-10]
	Ff_optimal = Float64[5.;0.;1]
    name = "as_1984_01"

    (; n1, n2, F, G, f, g, xy_init, xy_optimal, Ff_optimal, name)
end
