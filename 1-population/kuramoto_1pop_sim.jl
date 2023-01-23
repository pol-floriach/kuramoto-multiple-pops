# Script to integrate a system N of Kuramoto-Sakaguchi Oscillators, using the microscopical equations 
# Saves an animation of the oscillator's phase on the unit circle, and the global order parameter 
using Plots, DifferentialEquations, Distributions

# Function to integrate 
function kuramoto_micro!(dϕ,ϕ,K,t)
    @. ϕ      = mod2pi(ϕ)
    cosvec    = @. cos(ϕ); 
    sinvec    = @. sin(ϕ); 
    real_part = 1/N*sum(cosvec); 
    imag_part = 1/N*sum(sinvec); 
    α         = atan(imag_part,real_part);
    R         = ((real_part)^2 + ( imag_part)^2)^(.5);
    cosag     = cos(α+γ); 
    sinag     = sin(α+γ);
    @. dϕ     = ωj + K*R*(-sinvec*cosag + sinag*cosvec);

end
#----------------------------------------------------------------

# Callbacks to save order parameter (z = R⋅texp(iα))
function save_R(u,t,integrator)
    R = ((1/N*sum(sin.(u)))^2 + (1/N*sum(cos.(u)))^2)^(.5);
    return R
end
function save_α(u,t,integrator)
    α = atan(sum(sin.(u)),sum(cos.(u))); α = α + π;
    return α
end

saved_R = SavedValues(Float64,Float64);
saved_α = SavedValues(Float64, Float64);
cb_R = SavingCallback(save_R,saved_R,saveat=0:.016667:100);
cb_α = SavingCallback(save_α,saved_α,saveat=0:.016667:100);

cb = CallbackSet(cb_R,cb_α);


# Number of oscillators
const N = 500;
#natural frequencies
ωj = rand(Cauchy(0,1),N);
#initial conditions and time span
ϕ₀ = rand(Cauchy(0,1),N);
tspan = (0.0,100.0);

# Sakaguchi phase shift
γ = 0;
# Coupling force
K = 7;

# Integration 
prob = ODEProblem(kuramoto_micro!,ϕ₀,tspan,K);
sol = solve(prob,saveat=0.016667,callback=cb);
Rvec = saved_R.saveval
αvec = saved_α.saveval
ϕ_j = sol.u
#=
R = sum(Rvec[50:end])/length(Rvec[50:end]);
α = sum(αvec[50:end])/length(αvec[50:end]);
=#

logocolors = Colors.JULIA_LOGO_COLORS;
thetaj = 0:.1:2π+0.1

anim = @animate for ii in eachindex(sol.t[2:end])
        
        plot(cos.(thetaj), sin.(thetaj), color =:black, legend = false, grid = false, axis = false, linewidth = 2.5)
        plot!(0.25*cos.(thetaj), 0.25*sin.(thetaj), color =:black, legend = false, grid = false, axis = false)
        plot!(0.5*cos.(thetaj), 0.5*sin.(thetaj), color =:black, legend = false, grid = false, axis = false)
        plot!(0.75*cos.(thetaj), 0.75*sin.(thetaj), color =:black, legend = false, grid = false, axis = false)


        scatter!((cos.(ϕ_j[ii]), sin.(ϕ_j[ii])),
        markersize = 5, 
        markershape = :circ, 
        markercolor = logocolors.blue,
        legend = false,
         )
        plot!([0,Rvec[ii-1]*-cos(αvec[ii])], [0,-Rvec[ii]*sin(αvec[ii])],arrow = true, color=:red, linewidth = 3, legend = false)

        annotate!(0.8,0.9,sol.t[ii])
        annotate!(0,0.19, "0.25", fontsize = 1)

end


@time gif(anim, "kuramoto-1pop-animation.gif", fps = 10)
