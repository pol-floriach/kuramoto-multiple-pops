using DifferentialEquations, Plots, LinearAlgebra, Alert, LaTeXStrings, Statistics



function kuramoto_Mpop!(du,u,p,t)
    K_s, K_c, α = p        #parameters

    R = @view u[1:M]
    φs = @view u[M+1:end] 
    dR = @view du[1:M]
    dφ    = @view du[M+1:end];

    sinφ = @. sin(φs)
    cosφ = @. cos(φs)
    #R    = abs.(R)
    # meanfield
    a=mean(R.*cosφ)
    b=mean(R.*sinφ);

    Rmean=sqrt(a^2+b^2);
    Φ=atan(b,a);
    R = abs.(R)
       cosAux=@. cos(α+φs-Φ);
       sinAux=@. sin(α+φs-Φ);
    cosα=cos(α);
    sinα=sin(α);

    #ODE equation system
    @. dR = -Δ*R + 0.5*(K_s*cosα*(R-R^3) + K_c*(1-R^2)*Rmean*cosAux )
    @. dφ = -0.5*(K_s*sinα*(1+R^2) + K_c*(1/R + R)*Rmean*sinAux) 

    nothing
end


# Computation parameters
M = 500;                 # number of populations
Δ = 1.0;                 # disorder parameter
p = 0.9;                 # coupling ratio
K = 7.3;
                         # coupling strength
K_s = p*K;               # self  coupling  force 
K_c = (1-p)*K;           # cross coupling  force
α = 1.2;                 # sakaguchi phase shift 

# Initial conditions
R₀ = rand(M)*0.3;  
φ₀ = rand(M)*2π; φ₀ = mod2pi.(φ₀)
u₀ = [R₀; φ₀]

tspan = (0.0,10000.0)    # time span
par = K_s, K_c, α;

prob = ODEProblem(kuramoto_Mpop!, u₀1, tspan,par, dtmax = 0.1, dtmin = 1e-5);
sol = @time solve(prob, saveat = 0.1);

x=sol[1:M,:].*cos.(sol[M+1:2*M,:]);
y=sol[1:M,:].*sin.(sol[M+1:2*M,:]);

Rmean=sqrt.(mean(x,dims=1).^2+mean(y,dims=1).^2)
Rmean = vec(Rmean)

# Time series plot
p = plot(  sol.t,sol[1,:], 
            title = "Time series \n"*"K  = "*string(K)*", p = "*string(p), 
            label = "osc. 1, polar")
plot!(  sol.t,sol[320,:], 
        label = "osc. 2, polar")
plot!(  sol.t,Rmean,
        xlim=(1e3,4e3),
        lw=2, 
        label = "Mean, polar")