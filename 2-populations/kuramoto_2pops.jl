## Two populations of N Kuramoto oscillators each, with self-coupling and cross coupling
# Z_a = R_A⋅exp(iα)
# Z_b = R_B⋅exp(iβ) 

using DifferentialEquations, Plots, Distributions
#----------------------------------------------------------------
function kuramoto_2pop_micro(du,u,K_s,t)
    # system of two coupled populations of Kuramoto-Sakaguchi oscillators
    ϕ_A  = @view u[1:N];     @. ϕ_A = mod2pi(ϕ_A)
    ϕ_B  = @view u[N+1:end]; @. ϕ_B = mod2pi(ϕ_B)
    dϕ_A = @view du[1:N];
    dϕ_B = @view du[N+1:end];
    cosvec_A =  @. cos(ϕ_A); 
    sinvec_A =  @. sin(ϕ_A); 
    cosvec_B =  @. cos(ϕ_B); 
    sinvec_B =  @. sin(ϕ_B); 

    # computation of population A's R
    real_part_A = 1/N*sum(cosvec_A); 
    imag_part_A = 1/N*sum(sinvec_A); 
    α = atan(imag_part_A,real_part_A);
    R_A = ((real_part_A)^2 + (imag_part_A)^2)^(.5);
    
    # computation of population B's R
    real_part_B = 1/N*sum(cosvec_B); 
    imag_part_B = 1/N*sum(sinvec_B); 
    β = atan(imag_part_B,real_part_B);
    R_B = ((real_part_B)^2 + (imag_part_B)^2)^(.5);


    # auxiliary sin and cos for sin(α + γ - ϕⱼ) = sin(α + γ)cos(ϕⱼ) - sin(ϕⱼ)cos(α + γ), for γ_s and γ_c. 
    # α and γ_s
    cosαγ_s = cos(α+γ_s);
    sinαγ_s = sin(α+γ_s);
    # α and γ_c
    cosαγ_c = cos(α+γ_c);
    sinαγ_c = sin(α+γ_c);
    # β and γ_s 
    cosβγ_s = cos(β +γ_s);
    sinβγ_s = sin(β +γ_s);
    # α and γ_c
    cosβγ_c = cos(β+γ_c);
    sinβγ_c = sin(β+γ_c);
    
    # Dynamics

     @. dϕ_A = ω_A + K_s*R_A*(sinαγ_s*cosvec_A - cosαγ_s*sinvec_A) + K_c*R_B*(sinβγ_c*cosvec_A - cosβγ_c.*sinvec_A )

     @. dϕ_B = ω_B + K_s*R_B*(sinβγ_s*cosvec_B - cosβγ_s*sinvec_B) + K_c*R_A*(sinαγ_c*cosvec_B - cosαγ_c.*sinvec_B)

   # @. dϕ_A = ω_A + K_s*R_A*sin(α + γ_s -ϕ_A) + K_c*R_B*sin(β + γ_c - ϕ_A)
   # @. dϕ_B = ω_B + K_s*R_B*sin(β + γ_s -ϕ_B) + K_c*R_A*sin(α + γ_c - ϕ_B)
end

function kuramoto_2pop_macro(du,u,K_s,t)
    # mean field of a system of two coupled populations of Kuramoto-Sakaguchi oscillators
    # u = [R_a, R_b, α, β]
    R_A =   u[1]
    R_B =   u[2]
    α   =   u[3]
    β   =   u[4]
    
    du[1] = -Δa*R_A + 0.5*( K_s*R_A*cos(γ_s) + K_c*R_B*cos(γ_c+β-α)     - K_s*R_A^3*cos(γ_s) - K_c*R_A^2*R_B*cos(α-γ_c-β))    
    du[2] = -Δb*R_B + 0.5*( K_s*R_B*cos(γ_s) + K_c*R_A*cos(γ_c+α-β)     - K_s*R_B^3*cos(γ_s) - K_c*R_B^2*R_A*cos(β-γ_c-α)) 
    du[3] =          0.5*( -K_s*sin(γ_s)     + K_c*R_B/R_A*sin(γ_c+β-α) + K_s*R_A^2*sin(γ_s) - K_c*R_A*R_B*sin(α-γ_c-β))
    du[4] =          0.5*( -K_s*sin(γ_s)     + K_c*R_A/R_B*sin(γ_c+α-β) + K_s*R_B^2*sin(γ_s) - K_c*R_B*R_A*sin(β-γ_c-α))
    nothing     
end

function kuramoto_2pop_macro2(du,u,K_s,t)
    # 
    # u = [R_a, R_b, θ]
    R_A = u[1]
    R_B = u[2]
    θ   = u[3]

    du[1] = -Δa*R_A + 0.5*( K_s*R_A*cos(γ_s) + K_c*R_B*cos(γ_c - θ) - K_s*R_A^3*cos(γ_s) - K_c*R_A^2*R_B*cos(θ-γ_c))  
    du[2] = -Δb*R_B + 0.5*( K_s*R_B*cos(γ_s) + K_c*R_A*cos(γ_c + θ) - K_s*R_B^3*cos(γ_s) - K_c*R_B^2*R_A*cos(γ_c+θ)) 
    du[3] = 0.5*( K_s*sin(γ_s)*(R_A^2-R_B^2) + K_c*(R_B/R_A*sin(γ_c-θ) - R_A/R_B*sin(γ_c+θ) - R_A*R_B*(sin(θ - γ_c) + sin(θ+γ_c))))
    nothing
end

#---------------------------------------------------------------------------------------------------------------------
# callback function to save parameter order of both 
function save_R_A(u,t,integrator)
    R_A = ((1/N*sum(sin.(u[1:N])))^2 + (1/N*sum(cos.(u[1:N])))^2)^(.5)
    return R_A
end

function save_R_B(u,t,integrator)
    R_B = ((1/N*sum(sin.(u[N+1:end])))^2 + (1/N*sum(cos.(u[N+1:end])))^2)^(.5)
    return R_B
end

 function save_α(u,t,integrator)
    α = atan(sum(sin.(u[1:N])),sum(cos.(u[1:N])));
    return α
 end

 function save_β(u,t,integrator)
    β = atan(sum(sin.(u[N+1:end])),sum(cos.(u[N+1:end])));    
    return β
end

saved_R_A = SavedValues(Float64,Float64);
saved_R_B = SavedValues(Float64,Float64);
saved_α   = SavedValues(Float64, Float64);
saved_β   = SavedValues(Float64, Float64);

cb_R_A = SavingCallback(save_R_A,saved_R_A,saveat=0:.5:100);
cb_R_B = SavingCallback(save_R_B,saved_R_B,saveat=0:.5:100);
cb_α   = SavingCallback(save_α,saved_α, saveat= 0:.5:100);
cb_β   = SavingCallback(save_β,saved_β, saveat= 0:.5:100);

cb = CallbackSet(cb_R_A,cb_R_B,cb_α,cb_β);

#---------------------------------------------------------------------------------------------------------------------


# of oscillators
const N = 1000;


# parameters of Lorentzian distribution
Δa = 1; Δb = 1;
#natural frequencies
ω_A = rand(Cauchy(0,Δa),N) 
ω_B = rand(Cauchy(0,Δb),N)

#initial conditions
ϕ₀a = rand(Cauchy(0,Δa),N)
ϕ₀b = rand(Cauchy(0,Δb),N)
ϕ₀  = [ϕ₀a; ϕ₀b]


tspan = (0.0,100.0);


# parameters of Kuramoto-Sakaguchi model
γ_s = γ_c = 1.2;
K_c = 0.8;
#-------------------------------------------------------------------------------
# MICROSCOPIC

RAvec = zeros(length(0:.5:10));
RBvec = copy(RAvec);
αvec  = copy(RAvec);  
βvec  = copy(RAvec);
ii = 0;
for K_s = 0:.5:10
    global ii +=1;
    prob      = ODEProblem(kuramoto_2pop_micro,ϕ₀,tspan,K_s);
    sol       = solve(prob,saveat=0.1, callback = cb);
    RA_k      = saved_R_A.saveval;
    RB_k      = saved_R_B.saveval;
    α_k       = saved_α.saveval;
    β_k       = saved_β.saveval;
    RAvec[ii] = sum(RA_k[50:end])/length(RA_k[50:end]);
    RBvec[ii] = sum(RB_k[50:end])/length(RA_k[50:end]);
    αvec[ii]  = sum(α_k[50:end])/length(α_k[50:end]);
    βvec[ii]  = sum(β_k[50:end])/length(β_k[50:end]);
end 
plotly()
α_k = mod2pi.(α_k)
β_k = mod2pi.(β_k)
plot(α_k); plot!(β_k)
#-------------------------------------------------------------------------------
# MACROSCOPIC

R_A₀ = ((1/N*sum(sin.(ϕ₀a)))^2 + (1/N*sum(cos.(ϕ₀a)))^2)^(.5)
α₀ = atan(sum(sin.(ϕ₀a)),sum(cos.(ϕ₀a))) 
α₀ = mod2pi(α₀)

R_B₀ = ((1/N*sum(sin.(ϕ₀b)))^2 + (1/N*sum(cos.(ϕ₀b)))^2)^(.5)
β₀ = atan(sum(sin.(ϕ₀b)),sum(cos.(ϕ₀b))); β₀ = mod2pi(β₀)

u₀ = [R_A₀; R_B₀; α₀; β₀]


RAvec_macro = zeros(length(0:.1:10));
RBvec_macro = copy(RAvec_macro);
αvec_macro  = copy(RAvec_macro);  
βvec_macro  = copy(RAvec_macro);
ii = 0;
for K_s = 0:.1:10
    global ii +=1;
    prob = ODEProblem(kuramoto_2pop_macro,u₀,tspan,K_s,dtmax = 0.01);
    sol = solve(prob,saveat=0.1)

    RAvec_macro[ii] = sum(sol[1,50:end])/length(sol[1,50:end]); 
    RBvec_macro[ii] = sum(sol[2,50:end])/length(sol[2,50:end])
    αvec_macro[ii]  = sum(sol[3,50:end])/length(sol[3,50:end])
    βvec_macro[ii]  = sum(sol[4,50:end])/length(sol[4,50:end])
end
#-------------------------------------------------------------------------------
# MACROSCOPIC, REDUCED SYSTEM

θ₀ = atan(sum(sin.(ϕ₀a)),sum(cos.(ϕ₀a))) - atan(sum(sin.(ϕ₀b)),sum(cos.(ϕ₀b)))


u₀ = [R_A₀; R_B₀; θ₀]


RAvec_macro_2 = zeros(length(0:.1:10));
RBvec_macro_2 = copy(RAvec_macro_2);
θvec_macro_2  = copy(RAvec_macro_2);  
ii = 0;
for K_s = 0:.1:10
    global ii +=1;
    prob = ODEProblem(kuramoto_2pop_macro2,u₀,tspan,K_s,dtmax = 0.01);
    sol = solve(prob,saveat=0.1)

    RAvec_macro_2[ii] = sum(sol[1,50:end])/length(sol[1,50:end]); 
    RBvec_macro_2[ii] = sum(sol[2,50:end])/length(sol[2,50:end])
    θvec_macro_2[ii]  = sum(sol[3,50:end])/length(sol[3,50:end])
end
#-------------------------------------------------------------------------------
gr()
p1 = scatter(0:.5:10,RAvec, label = "R_A");
     plot!(0:.1:10,RAvec_macro, label = "R_A", ls = :dash)
     scatter!(0:.5:10,RBvec, label = "R_B", title = "R",legend=:bottomright)
     plot!(0:.1:10,RBvec_macro, label = "R_B", title = "R", ls = :dash)

    xlabel!("K_s");
p2 = scatter(0:.5:10, αvec, label = "α");
     plot!(0:.1:10, αvec_macro, label = "α", ls = :dash); 
     scatter!(0:.5:10,βvec, label = "β", title = "α, β");
     plot!(0:.1:10,βvec_macro, label = "β", ls = :dash, legend = :bottomleft)
     xlabel!("K_s")

ptot = plot(p1,p2, layout = (2,1))

plot(0:.1:10, mod2pi.(θvec_macro_2), label = "R_B, 3 eqs")

