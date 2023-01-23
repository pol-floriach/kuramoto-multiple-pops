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
The phase Kuramoto oscillator $\phi_j$ is given by: 
$$\dot\phi_j = \omega_j + K|R|\sin(\Phi - \phi_j + \alpha),$$
Where $\omega_j$ is the oscillator natural frequency, K is the coupling force, $R = |R|e^{i\Phi}$ is the complex order parameter, and $\alpha$ the Sakaguchi phase shift.

Given that $\omega_i$ follow a Lorentzian distribution with parameter $\Delta$, one can describe the N dimensional system with a 2 dimensional one with the evolution of the Kuramoto order parameter:
$$ \dot R = -\Delta R + \frac{K}{2} (R e^{i\alpha} - |R|^2R^*e^{-i\alpha})$$

+ `kuramoto_1pop`: Simulation example of a Kuramoto-Sakaguchi system of oscillators, using both the microscopical equations and a mean field model, defining an order parameter for the system (and using the Ott-Antonsen ansatz).

+ `kuramoto_1pop_sim`: Integrates the microscopical system and outputs an animation of the N oscillators phase on the unit circle.

### Two populations
Consider 2 populations, which we name A and, B of N oscillators each. The dynamics of an oscillator, described by their phase $\phi_j$ in population $\sigma$ are given by:
$$ \dot\phi_j = \omega_j^{\sigma} + K_sR_{\sigma}\sin(\Phi_{\sigma} - \phi_j^{(\sigma)} - \alpha) + K_cR^{(\kappa)}\sin(\Phi_{\kappa} - \phi_j^{(\sigma)} - \alpha)$$ 
Where $K_s$ is the coupling force among a population, and $K_c$ the one from the other population, $\kappa$ the other population and $Z_{\sigma} = R_{\sigma}e^{i\Phi_{\sigma}}$.

These set of $2N$ ODEs can be reduced to a single ODE for each:

$$ \dot Z_{\sigma} = -\Delta Z_{\sigma} + \frac{K_s}{2}( Z_{\sigma}e^{-i\alpha} - Z_\sigma^2Z_{\sigma}e^{i\alpha}) + \frac{K_c}{2}(Z_{\kappa}e^{-i\alpha} - Z_{\sigma}^2Z_{\kappa}^{*}e^{i\alpha})$$

Describing $K_s = Kp$ and $K_c = (1-p)K$, where p controls the amount of coupling within a population, and analyzing a symmetrical case given by $Z_{\sigma}$ is the same for both populations, one obtains the equilibrium parameter order modulus as:
$$R = \sqrt{1 - \frac{2\Delta}{K\cos(\alpha)}}$$ 

+ `kuramoto_2pops.jl`: Simulation of a system of 2 Kuramoto-Sakaguchi oscillator populations, with self coupling and cross coupling. Also using the microscopical equations and a reduced mean-field one.

+ `kuramoto_2pops_fixp.jl`: Finds fixed points of the reduced system, determines their stability and integrates the system to make sure results are correct. 

+ `kuramoto_2pops_symm.jl`: Simulates the system for a necessary amount of time for it to reach equilibrium. It does so with a given set of parameters such that both populations are symmetrical (same order parameter for both) and on an antiphase state (same order parameter modulus, but with a 180º shift).
It also computes the analytical equilibrium point for those parameters. 

### Multiple populations 
Consider $M$ populations composed of $N$ Kuramoto oscillatros each with an all-to-all topology, with a fraction of self-connectivity $p$ and $1-p$ for the mean-field. One can describe the Ott-Antonsen equations as:

$$ \dot Z_{\sigma} = -\Delta Z_{\sigma} + \frac{K}{2} (\overline Z e^{-i\alpha} - Z_{\sigma}^2 \overline Z^{*}e^{i\alpha})$$
where 
$$ \overline Z = \frac{1}{M}\sum_{\kappa = 1}^M Z_{\sigma}$$
 is the global Kuramoto order parameter for all the M populations

+ `kuramoto_Mpops_cartesian.jl`: Simulation of a system of M all-to-all coupled populations of N Kuramoto oscillators. System reduced from M·N equations to 2M equations, using the order parameter of each population in cartesian form.
+ `kuramoto_Mpops_polar.jl`: Essentially the same code as + `kuramoto_Mpops_cartesian.jl` but with the equations in polar form.

+ `kuramoto_Mpops_poincare_section.jl`: Produces a plot of a Poincaré section of the M population system, for a desired range of coupling constant K.

+ `kuramoto_Mpops_bifurcation_diagram.jl`: Produces a bifurcation diagram of the M population system. 

See photo gallery for images of Poincaré sections and bifurcation diagrams.
### Ring of oscillators
+ `kuramoto_ring.jl`: Simulation of a ring of N Kuramoto-Sakaguchi oscillators, coupled all to all. 

## Dependencies 
[Julia](https://github.com/JuliaLang/julia)

[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

[Plots.jl](https://github.com/JuliaPlots/Plots.jl)

[Distributions.jl](https://github.com/JuliaStats/Distributions.jl)

[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)

[Peaks.jl](https://github.com/halleysfifthinc/Peaks.jl)

[LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl)

## Photo gallery
Go to [gallery](gallery/) to see figures produced using the above codes.