# FractionalDelayFilters

1-D fractional shift filters. Filters implemented: all-pass anti-causal Thirian from [1].

To install, start Julia REPL, press `]`, and run the command `add https://github.com/RoyCCWang/FractionalDelayFilters`

Example folder:
`shift_Thirian.jl` illustrates how to use this package to shift a sequence that starts at time `a` and ends at time `b` by `Δp` amount. Positive `Δp` is a delay (right shift), vice versa for negative `Δp`.

 ### References
 [1] Condat, Laurent, and Dimitri Van De Ville. "Fully reversible image rotation by 1-D filtering." 2008 15th IEEE International Conference on Image Processing. IEEE, 2008.
