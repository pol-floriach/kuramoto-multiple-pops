# Multiple populations of Kuramoto oscillators

Collection of Julia codes for the study of symmetry breaking states on the Kuramoto model oscillators, mainly in 2 and multiple populations, using the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) library.

## Usage

To run a file, either run:

```
julia <filename.jl>
```
on a terminal emulator or:
```
include("path/to/dir/filename.jl")
```
on an active Julia REPL.
## Contents
### One population
+ `kuramoto_1pop`: Simulation example of a Kuramoto-Sakaguchi system of oscillators, using both the microscopical equations and a mean field model (using the Ott-Antonsen ansatz).

+ `kuramoto_1pop_sim`: Integrates the microscopical system and outputs an animation of the N oscillators phase on the unit circle.

### Two populations
+ `kuramoto_2pops.jl`: Simulation of a system of 2 Kuramoto-Sakaguchi oscillator populations, with self coupling and cross coupling. Also using the microscopical equations and a reduced mean-field one.

+ `kuramoto_2pops_fixp.jl`: Finds fixed points of the reduced system, determines their stability and integrates the system to make sure results are correct. 

+ `kuramoto_2pops_symm.jl`: Simulates the system for a necessary amount of time for it to reach equilibrium. It does so with a given set of parameters such that both populations are symmetrical (same order parameter for both) and on an antiphase state (same order parameter modulus, but with a 180º shift).
It also computes the analytical equilibrium point for those parameters. 

+ `kuramoto_2pops_symm_phase_diagr`: Produces a contour plot 
### Multiple populations 

### Ring of oscillators


## Dependencies 
[Julia](https://github.com/JuliaLang/julia)

[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

[Plots.jl](https://github.com/JuliaPlots/Plots.jl)

[Distributions.jl](https://github.com/JuliaStats/Distributions.jl)

[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)

[Peaks.jl](https://github.com/halleysfifthinc/Peaks.jl)

[LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl)

## Photo gallery