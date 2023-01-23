## Simulation of N Kuramoto-Sakaguchi oscillators

using DifferentialEquations, Plots, Distributions
#----------------------------------------------------------------
#ode system 
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
function kuramoto_macro!(du,u,K,t)
    # u = [R, α]
    du[1] = -Δ*u[1] + K*u[1]/2*cos(γ) - K*u[1]^3/2*cos(γ)
    du[2] =           K/2*sin(γ)      + K*u[1]^2/2*sin(γ)
    nothing
end
#----------------------------------------------------------------

#callbacks to save parameter order (z = R⋅texp(iα))
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
cb_R = SavingCallback(save_R,saved_R,saveat=0:.5:100);
cb_α = SavingCallback(save_α,saved_α, saveat= 0:.5:100);

cb = CallbackSet(cb_R,cb_α);
#----------------------------------------------------------------

# Number of oscillators
const N = 1000;
#natural frequencies
ωj = rand(Cauchy(0,1),N);
#initial conditions and time span
ϕ₀ = rand(Cauchy(0,1),N);
tspan = (0.0,100.0);
###-----------------------------------------------------------###
#                        MICROSCOPIC                            #
###-----------------------------------------------------------###

γ = 0;
Rvec1 = zeros(length(0:.5:10));
iter_number = 0;
for K = 0:.5:10
   global iter_number +=1;
    prob = ODEProblem(kuramoto_micro!,ϕ₀,tspan,K);
    sol = solve(prob,saveat=0.1,callback=cb);
    R_k = saved_R.saveval;
    α_k = saved_α.saveval;

    Rvec1[iter_number] = sum(R_k[50:end])/length(R_k[50:end]);
end 

#---------------------------------------------------------------------------

γ = π/4;
Rvec2 = zeros(length(0:.5:10));
iter_number = 0;
for K = 0:.5:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_micro!,ϕ₀,tspan,K);
    sol = solve(prob,saveat=0.1,callback=cb);

    R_k = saved_R.saveval;
    α_k = saved_α.saveval;

    Rvec2[iter_number] = sum(R_k[50:end])/length(R_k[50:end]);
end 

#---------------------------------------------------------------------------
# CANVIANT PARAMETRE D'ESCALA
ωj = rand(Cauchy(π,2),N);
ϕ₀ = rand(Cauchy(π,2),N);
tspan = (0.0,100.0);

γ = 0;
Rvec3 = zeros(length(0:.5:10));
iter_number = 0;
for K = 0:.5:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_micro!,ϕ₀,tspan,K);
    sol = solve(prob,saveat=0.1,callback=cb);

    R_k = saved_R.saveval;
    α_k = saved_α.saveval;

    Rvec3[iter_number] = sum(R_k[50:end])/length(R_k[50:end]);
end 

#---------------------------------------------------------------------------

γ = π/4;
Rvec4 = zeros(length(0:.5:10));
iter_number = 0;
for K = 0:.5:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_micro!,ϕ₀,tspan,K);
    sol = solve(prob,saveat=0.1,callback=cb_R);

    R_k = saved_R.saveval;
    α_k = saved_α.saveval;

    Rvec4[iter_number] = sum(R_k[50:end])/length(R_k[50:end]);
end 

#---------------------------------------------------------------------------

γ = π/2;
Rvec5 = zeros(length(0:.5:10));
iter_number = 0;
for K = 0:.5:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_micro!,ϕ₀,tspan,K);
    sol = solve(prob,saveat=0.1,callback=cb_R);

    R_k = saved_R.saveval;
    α_k = saved_α.saveval;
    Rvec5[iter_number] = sum(R_k[50:end])/length(R_k[50:end]);
end 

#---------------------------------------------------------------------------
##  MACROSCOPIC

Δ = 1;
ϕ₀ = rand(Cauchy(0,Δ),N);
tspan = (0.0,100.0);

R₀ = ((1/N*sum(sin.(ϕ₀)))^2 + (1/N*sum(cos.(ϕ₀)))^2)^(.5);
α₀ = atan(sum(sin.(ϕ₀)),sum(cos.(ϕ₀))); α₀ = mod2pi(α₀);
z₀ = [R₀;α₀];

γ = 0;

Rvec_macro1 = zeros(length(0:.1:10));
iter_number = 0;
for K = 0:.1:10
    global iter_number +=1;
    prob = ODEProblem(kuramoto_macro!,z₀,tspan,K,dtmax = 0.01);
    sol = solve(prob,saveat=0.1)

    Rvec_macro1[iter_number] = sum(sol[1,50:end])/length(sol[1,50:end])
end
#-------------------------------------
γ = π/4;
Rvec_macro2 = zeros(length(0:.1:10));
iter_number = 0;
for K = 0:.1:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_macro!,z₀,tspan,K,dtmax = 0.01);
    sol = solve(prob,saveat=0.1)

    Rvec_macro2[iter_number] = sum(sol[1,50:end])/length(sol[1,50:end])
end

#-------------------------------------

Δ = 2;
ϕ₀ = rand(Cauchy(0,Δ),N);
tspan = (0.0,100.0);

R₀ = ((1/N*sum(sin.(ϕ₀)))^2 + (1/N*sum(cos.(ϕ₀)))^2)^(.5);
α₀ = atan(sum(sin.(ϕ₀)),sum(cos.(ϕ₀))); α₀ = mod2pi(α₀);
z₀ = [R₀;α₀];

γ = 0;
Rvec_macro3 = zeros(length(0:.1:10));
iter_number = 0;

for K = 0:.1:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_macro!,z₀,tspan,K,dtmax = 0.01);
    sol = solve(prob,saveat=0.1)

    Rvec_macro3[iter_number] = sum(sol[1,50:end])/length(sol[1,50:end])
end
#-------------------------------------
γ = π/4;

Rvec_macro4 = zeros(length(0:.1:10));
iter_number = 0;
for K = 0:.1:10
    global iter_number +=1; 
    prob = ODEProblem(kuramoto_macro!,z₀,tspan,K,dtmax = 0.01);
    sol = solve(prob,saveat=0.1)

    Rvec_macro4[iter_number] = sum(sol[1,50:end])/length(sol[1,50:end])
end


#---------------------------------------------------------------------------
## Plots
gr();
scatter(0:.5:10,Rvec1,label="γ  = 0,   Δ = 1")
scatter!(0:.5:10,Rvec2,label="γ = π/4, Δ = 1")
scatter!(0:.5:10,Rvec3,label="γ = 0,   Δ = 2")
scatter!(0:.5:10,Rvec4,label="γ = π/4, Δ = 2")
scatter!(0:.5:10,Rvec5, label = "γ = π/2, Δ = 2",legend=:topleft)
plot!(0:.1:10,Rvec_macro1, ls = :dash,label="γ = 0,   Δ = 1")
plot!(0:.1:10,Rvec_macro2, ls = :dash,label="γ = π/4, Δ = 1")
plot!(0:.1:10,Rvec_macro3, ls = :dash,label="γ = 0,   Δ = 2")
plot!(0:.1:10,Rvec_macro4, ls = :dash,label="γ = π/2, Δ = 2")

title!("Order pareameter as a function of the coupling force")
title!("Paràmetre d'ordre en funció de la força d'acoblament"); 
xaxis!("K"); yaxis!("|R|")
#----------------------------------------------

