
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
U = P
f_U = f.(U)

## filter.
s = real.(f_U)

N_thirian = 15
b_thirian = Vector{Float64}(undef, N_thirian)

τ = 0.31
applydelaythirian!(b_thirian, s, τ)

## visualize.

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, s, label = "s", "--")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("reversible fractional delay filtering")


x_ranges = [ P ]
s_itp, d_s_itp, d2_s_itp = Utilities.setupcubicitp(real.(f_U), x_ranges, 1.0)
s_func = xx->s_itp([xx])

S = FiniteSignalType(real.(f_U), P[1], P[end])

# no shifts.
j = 14 # test this time index.
println("no shift: j ", j)
s0 = real.(f_U)
println("zero: evalsignal( S, j ) - s0[j] = ", evalsignal( S, j ) - s0[j])
println("zero: s0[j] - s_func(P[j] = ", s0[j] - s_func(P[j]))
println()

# integer shifts.
d = 19
applyintegershift!(S, d)

j = d+4 # test this time index.
println("shift: j ", j)
evalsignal( S, j ) # evaluates s0[j-d].
println("zero: s0[j-d] - evalsignal( S, j ) = ", s0[j-d] - evalsignal( S, j ))
println()

# fractional shift, τ ∈ [0, 0.5].
ΔP = P[2]-P[1]

d = 2 # this is the integer delay (positive).
applyintegershift!(S, d)

τ = 0.319 # fractional delay (positive).

resetsignal!(S)
applydelaythirian!(b_thirian, S.s, τ)
evals_S = collect( evalsignal(S, j) for j = 1:length(S.s) )

evals_s_func = collect( s_func( P[j]-ΔP*τ ) for j = 1:length(P) )

println("norm difference between cubic itp and fractional delay: ", norm(evals_s_func - evals_S))
#[evals_S evals_s_func]

# visualize.
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, evals_S, label = "evals_S", "--")
PyPlot.plot(P, evals_s_func, label = "evals_s_func", "--")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("cubic itp vs. fractional delay filtering")

### test reversible property of thirian delay filter. see https://lcondat.github.io/publis/condat_icip08_rotations.pdf
τ = 0.21
resetsignal!(S)
applydelaythirian!(b_thirian, S.s, τ)

tmp_s = copy(S.s)

S.s[:] = reverse(S.s)
applydelaythirian!(b_thirian, S.s, τ)
S.s[:] = reverse(S.s)


#reverse!(tmp_s)
applyadvancethirian!(b_thirian, tmp_s, τ) # bake reverse!() into this function.
#reverse!(tmp_s)
println("test applyadvancethirian!: ", norm(tmp_s-S.s))

## timing.
# import BenchmarkTools
# BenchmarkTools.@btime applydelaythirian!(b_thirian, tmp_s, τ)
# BenchmarkTools.@btime applyadvancethirian!(b_thirian, tmp_s, τ)

evals_S = collect( evalsignal(S, j) for j = 1:length(S.s) )
println("reversible test: ", norm(evals_S - s0))
println()

### arb shift, test reversible, frac part ∈ [0, 0.5]
resetsignal!(S)

shift = 1.23
applyshift!(S, shift, b_thirian)

shift = -1.23
applyshift!(S, shift, b_thirian)
evals_S = collect( evalsignal(S, j) for j = 1:length(S.s) )
println("reversible test: ", norm(evals_S - s0))
println()

### arb shift, test reversible, frac part ∈ [0.5, 1]
resetsignal!(S)

shift = 1.9
applyshift!(S, shift, b_thirian)

shift = -1.9
applyshift!(S, shift, b_thirian)
evals_S = collect( evalsignal(S, j) for j = 1:length(S.s) )
println("reversible test: ", norm(evals_S - s0))
println()


######## get shift coordinates from Δp, which is in units of P.
Δp = 2.0
a = P[1]
b = P[end]

resetsignal!(S)

shift = Δp*length(S.s)/(b-a)
applyshift!(S, shift, b_thirian)


# visualize.
evals_S = collect( evalsignal(S, j) for j = 1:length(S.s) )
evals_s_func = collect( s_func( P[j]-ΔP*τ ) for j = 1:length(P) )

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(f_U), label = "f")
PyPlot.plot(P, evals_S, label = "evals_S", "--")
#PyPlot.plot(P, evals_s_func, label = "evals_s_func", "--")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("fractional shift with Δp = $(Δp)")
