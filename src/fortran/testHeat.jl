include("heat.jl")
N = 5
a = zeros(N, N)
a[[1, N],  :] .= 1.
a[:, [1, N]] .= 1.
b = copy(a)
err = heat(a, b)
show(stdout, "text/plain", b);
