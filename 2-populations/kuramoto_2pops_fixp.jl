## FIXED POINTS AND STABILITY OF SYSTEM OF 2 POPULATIONS OF KURAMOTO-SAKAGUCHI OSCILLATORS,
## HAVING REDUCED THE ODE SYSTEM TO A 3 DIMENSOPNAL ONE USING THE OTT-ANTONSEN ANSATZ AND 
## EXPLOITED A SYMMETRY OF THE PROBLEM


#-------------------------------------------------------------------------
using NonlinearSolve, StaticArrays, Plots, LinearAlgebra, DifferentialEquations, Distributions, LaTeXStrings

# F(R_A, R_B, θ) = 0 to find fixed points
function kuramoto_2pops_fixed_points(u,_)
    R_A = u[1]
    R_B = u[2]
    θ   = u[3]

    F1 = -Δa*R_A + 0.5*( K_s*R_A*cos(γ_s) + K_c*R_B*cos(γ_c - θ) - K_s*R_A^3*cos(γ_s) - K_c*R_A^2*R_B*cos(θ-γ_c))  
    F2 = -Δb*R_B + 0.5*( K_s*R_B*cos(γ_s) + K_c*R_A*cos(γ_c + θ) - K_s*R_B^3*cos(γ_s) - K_c*R_B^2*R_A*cos(γ_c+θ)) 
    F3 = 0.5*( K_s*sin(γ_s)*(R_A^2-R_B^2) + K_c*(R_B/R_A*sin(γ_c-θ) - R_A/R_B*sin(γ_c+θ) - R_A*R_B*(sin(θ - γ_c) + sin(θ+γ_c))))

    @SVector[F1, F2, F3]
    #probN = NonlinearProblem{false}(f, u0)
    #solver = solve(probN, NewtonRaphson(), tol = 1e-9)
end

# jacobian of the ODE system to determine the stability of the fixed points
function kuramoto_2pops_jacobian(u)

    R_A = u[1]
    R_B = u[2]
    θ   = u[3]
    J = zeros(3,3)
    J[1,1] = -Δ + 0.5*K_s*cos(γ_s) - 3/2*K_s*R_A^2*cos(γ_s) - K_c*R_A*R_B*cos(θ-γ_c)
    J[1,2] = 0.5*(K_c*cos(γ_c-θ) - K_c*R_A^2*cos(θ - γ_c ))
    J[1,3] = 0.5*(K_c*R_B*sin(γ_c-θ) + K_c*R_B*R_A^2*sin(θ-γ_c))

    J[2,1] = 0.5*(K_c*cos(γ_c+θ) -K_c*R_B^2*cos(θ+γ_c))
    J[2,2] = -Δ + 0.5*K_s*cos(γ_s) - 3/2*K_s*R_B^2*cos(γ_s) - K_c*R_A*R_B*cos(θ+γ_c)
    J[2,3] = 0.5*(-K_c*R_A*sin(θ+γ_c) + K_c*R_A*R_B^2*sin(θ+γ_c))
    
    J[3,1] =  K_s*sin(γ_s)*R_A +K_c*(-R_B/R_A^2*sin(γ_c-θ) - sin(θ+γ_c)/(2*R_B) - R_B/2*(sin(θ-γ_c) + sin(θ+γ_c)))
    J[3,2] = -K_s*sin(γ_s)*R_B +K_c*( R_A/R_B^2*sin(γ_c+θ) + sin(γ_c-θ)/(2*R_A) - R_A/2*(sin(θ-γ_c) + sin(θ+γ_c)))
    J[3,3] = -0.5*K_c*(+R_B/R_A*cos(γ_c-θ) + R_A/R_B*cos(γ_c+θ) + R_A*R_B*(cos(θ-γ_c) + cos(θ+γ_c)))
    return J
end

function kuramoto_2pop_macro_reduced(du,u,_,t)
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
 #-------------------------------------------------------------------------
const Δ = 1; Δa = Δb = Δ
#-------------------------------------------------------------------------
# ROOT FINDING

# Parameters
K_s = 2.0;
K_c = 0.2;
γ_s = 0;
γ_c = 0;

u₀ = @SVector[0.68,0.678,2.5];
# fixed point
prob = NonlinearProblem{false}(kuramoto_2pops_fixed_points,u₀);
fixed_point = solve(prob,NewtonRaphson(), tol = 1e-12)

# eigenvalues of jacobian at fixed point
J = kuramoto_2pops_jacobian(fixed_point);
vals = eigvals(J)
plotly();
#=
scatter(real(vals), imag(vals), markersize = 10)
plot!(zeros(length(-5:.1:5)), -5:.1:5, color = :black)
plot!(-1:.1:1,zeros(length(-1:.1:1)), color = :black)
=#

#-----------------------------------------------------------
# SIMULACIÓ PER COMPROVAR COMPORTAMENT
#=
R_A₀ = 0.2
R_B₀ = 0.1
θ₀   = 0.2
=#
N = 1000;

distr(ii)=tan(π*0.5*(2*ii-N-1)/(N+1));              
ω_A = Δa*distr.(1:N); ω_B = Δb*distr.(1:N);


#ω_A = rand(Cauchy(0,Δ),N) 
#ω_B = rand(Cauchy(0,Δ),N)
ϕ₀a = rand(Cauchy(0,Δ),N)
ϕ₀b = rand(Cauchy(0,Δ),N)

R_A₀ = ((1/N*sum(sin.(ϕ₀a)))^2 + (1/N*sum(cos.(ϕ₀a)))^2)^(.5)
R_B₀ = ((1/N*sum(sin.(ϕ₀b)))^2 + (1/N*sum(cos.(ϕ₀b)))^2)^(.5)
α₀ = atan(sum(sin.(ϕ₀a)),sum(cos.(ϕ₀a))); α₀ = mod2pi.(α₀)
β₀ = atan(sum(sin.(ϕ₀b)),sum(cos.(ϕ₀b))); β₀ = mod2pi.(β₀)
θ₀ = α₀ - β₀; θ₀ = mod2pi.(θ₀)

u₀ = [R_A₀; R_B₀; θ₀]
tspan = (0.0,1000.0)


K_s = 7;
K_c = 0.8;
γ_s = 1.2;
γ_c = 1.2;

# MACRO, 3 EQS
u₀ = [R_A₀; R_B₀; θ₀]
prob = ODEProblem(kuramoto_2pop_macro_reduced,u₀,tspan,K_s,dtmax = 0.01);
sol_macro = solve(prob,saveat=0.1)
# MACRO, 4 EQS
u₀ = [R_A₀; R_B₀;α₀; β₀]
prob = ODEProblem(kuramoto_2pop_macro,u₀, tspan, K_s, dtmax = 0.01);
sol_2 = solve(prob, saveat = 0.1)
# MICRO
ϕ₀  = [ϕ₀a; ϕ₀b];
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

cb_R_A = SavingCallback(save_R_A,saved_R_A,saveat=0:1:1000);
cb_R_B = SavingCallback(save_R_B,saved_R_B,saveat=0:1:1000);
cb_α   = SavingCallback(save_α,saved_α, saveat= 0:1:1000);
cb_β   = SavingCallback(save_β,saved_β, saveat= 0:1:1000);

cb = CallbackSet(cb_R_A,cb_R_B,cb_α,cb_β);

prob      = ODEProblem(kuramoto_2pop_micro,ϕ₀,tspan,K_s);
sol       = solve(prob,saveat=0.1, callback = cb);
RA_k      = saved_R_A.saveval;
RB_k      = saved_R_B.saveval;
α_k       = saved_α.saveval; β_k = saved_β.saveval;
θ_k       = @. mod2pi(α_k - β_k);


p1 = plot(sol_macro.t,sol_macro[1,:], label = "R_A", legend = :bottomright)
     plot!(sol_macro.t,sol_macro[2,:], label = "R_B")
     scatter!(0:1:1000,RA_k, label = "R_A, micro")
     scatter!(0:1:1000,RB_k, label = "R_B, micro")
p2 = plot(0:1:1000,RA_k, label = "R_A, micro")
     plot!(0:1:1000,RB_k, label = "R_B, micro")

     yaxis!("R")
     title!("B")



p2 = plot(sol_macro.t,mod2pi.(sol[3,:]), label = "θ")
xaxis!("t"); yaxis!("θ");
scatter!(0:5:500,θ_k, label = "θ, micro")


plot(p1,p2, layout = (2,1))


#------------------------------------------------------------
#------------------------------------------------------------
# ESCOMBRAT MÈTODE NEWOTN
u₀ = @SVector[0.42, 0.41, 0.0]
K_s = 2.0;
K_c = 1
γ_s = π/4;
γ_c = 0; 

prob = NonlinearProblem{false}(kuramoto_2pops_fixed_points,u₀);
fixed_point = solve(prob,NewtonRaphson(), tol = 1e-12);
fixed_point = fixed_point.u;
fixed_point_vec = [fixed_point_vec, fixed_point];
u₀ = fixed_point; 
K_c += 0.5

fixed_point_vec = [@SVector[0.0, 0.0, 0.0]];
u₀ = @SVector[0.42, 0.41, 0.0];
for ii = 1:.5:5
    global K_c = ii
    prob = NonlinearProblem{false}(kuramoto_2pops_fixed_points,u₀);
    fixed_point = solve(prob,NewtonRaphson(), tol = 1e-12)
    fixed_point = fixed_point.u; println(fixed_point)
   # fixed_point_vec = [fixed_point_vec, fixed_point];
   push!(fixed_point_vec, fixed_point)
    u₀ = fixed_point
end


#------------------------------------------------------------------
# ESCOMBRAT AMB SCIML.JL STEADY STATE PROBLEM 

K_s = 8 
γ_s = 1.2
γ_c = 1.1
u₀ = [0.5, 0.5, π]

fixed_point_vec = [[0.0+0*im; 0.0+0*im; 0.0+0*im]]
vaps = copy(fixed_point_vec)
maxs_a = [];
maxs_b = [];
mins_a = [];
mins_b = [];
for ii = 0.5:.1:1.5
    global K_c = ii; println(K_c)
    steadyp = SteadyStateProblem(kuramoto_2pop_macro_reduced,u₀);
    sol_newt = solve(steadyp);
    fixed_point = sol_newt.u;
    if ii == 1.1
    u₀ = [0.6, 0.6, 0];
    else
    u₀ = fixed_point; 
    end
    J = kuramoto_2pops_jacobian(fixed_point);
    vals = eigvals(J);

    prob = ODEProblem(kuramoto_2pop_macro_reduced, u₀, tspan, K_s, dtmax = 0.01)
    sol_odeprob = solve(prob, saveat = 0.1)

    max_a = max(sol_odeprob[1,:]); 
    max_b = max(sol_odeprob[2,:]); 
    min_a = 
    

    push!(fixed_point_vec,fixed_point);
    push!(vaps,vals);
    push!(maxs_a, maxs_a);
    push!(mins_a, min_a);
    push!(maxs_b, maxs_b);
    push!(mins_b, min_b);
end

fixed_point_vec = fixed_point_vec[2:end]
vaps = vaps[2:end]
 
paleta1 = palette([:yellow, :green], 12);
paleta2 = palette(:Blues_6); paleta2 = paleta2[2:end];
p1 = plot(zeros(length(-2:.1:2)), -2:.1:2, color = :black, label = false);
paleta1 = palette(:GnBu);
for jj = 1:7
    scatter!(vaps[jj], palette = paleta1, label = Kcvec[jj], legend = :topleft, legend_title = L"K_c") 
end

for jj = 8:11
    scatter!(vaps[jj], palette = paleta1, label = Kcvec[jj])
end

ylabel!("Im(λ)")
xlabel!("Re(λ)")

p1


