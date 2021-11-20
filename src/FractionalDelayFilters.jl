module FractionalDelayFilters

using LinearAlgebra, FFTW

include("../src/Thirian/fractional_shift.jl")

export FiniteSignalType, resetsignal!, applyshift!, evalsignal

end # module
