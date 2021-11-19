# package up, and have a coordiante transform from shift,a,b to u,1,length(s).

include("../src/Thirian/fractional_shift.jl")

# include("../src/FractionalDelayFilters.jl")
# import .FractionalDelayFilters


using LinearAlgebra
using FFTW
import PyPlot

import Utilities # https://gitlab.com/RoyCCWang/utilities

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


# sine for debug.
f = sin
P = LinRange(-pi, pi, 3000)
s = f.(P)

## set up.
a = P[1]
b = P[end]

S = FiniteSignalType(s, a, b)
N_thirian = 15
b_thirian = Vector{Float64}(undef, N_thirian)


######## get shift coordinates from Δp, which is in units of P.
Δp = 2.0

resetsignal!(S)
shift = Δp*length(S.s)/(b-a)
applyshift!(S, shift, b_thirian)


# visualize.
evals_S = collect( evalsignal(S, j) for j = 1:length(S.s) )

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, s, label = "s")
PyPlot.plot(P, evals_S, label = "evals_S", "--")
#PyPlot.plot(P, evals_s_func, label = "evals_s_func", "--")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("fractional shift with Δp = $(Δp)")
