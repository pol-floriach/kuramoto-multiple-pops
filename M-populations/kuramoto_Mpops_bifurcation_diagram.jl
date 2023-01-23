# BIFURCATION DIAGRAM OF A SYSTEM OF M POPULATIONS OF KURAMOTO-SAKAGUCHI  OSCILLATORS
# WITH SAKAGUCHI PHASE SHIFT = 0.2, 

# Package dependencies
using Plots, DifferentialEquations, Distributions, Interpolations, Peaks, Statistics

# Dynamical system in cartesian coordinates
function kuramoto_Mpop_carts!(du,u,p,t)
    # mean field of a system of P coupled populations of N Kuramoto-Sakaguchi oscillators, 2P equations using the Ott-Antonsen ansatz
    K_s, K_c, α = p        #parameters

    x = @view u[1:M]
    y = @view u[M+1:end] 
    dx = @view du[1:M]
    dy    = @view du[M+1:end];

    cosα=cos(α);
    sinα=sin(α);

    mult1 = K_s*(x*cosα - y*sinα)
    mult2 = K_s*(x*sinα + y*cosα)
    mult3 = K_c*(mean(x)*cosα - mean(y)*sinα)
    mult4 = K_c*(mean(x)*sinα + mean(y)*cosα)

    #ODE equation system
    @. dx = -Δ*x + 0.5*mult1 + 0.5*mult3 - (x^2 - y^2)*(mult1 + mult3) - 2*x*y*(mult2+mult4)
    @. dy = -Δ*y + 0.5*mult2 + 0.5*mult4 + (x^2 - y^2)*(mult2 + mult4) - 2*x*y*(mult1+mult3)

    nothing
end
# Bifurcation Diagram initialization
p1 = plot(0,0,xlim = (9.0,13), ylim = (0,1), legend = false)

# PARAMETERS
global M = 500;                 # number of populations
global Δ = 1.0;                 # disorder parameter
p = 0.7;                        # coupling ratio
# p = 0.9            
if p == 0.7
    global Kvec = 9.4:.01:15      # Desired range of K
elseif p== 0.9
    global Kvec = 6.1:.01:7.3
end
α = 1.2;                        # sakaguchi phase shift 
Tend = 50000.0                   # elapsed time

# Variables for informative purposes
total_its = string(length(Kvec));
itnumber = 0;

for K in Kvec
    # Parameters
    K_s = p * K;                 # self  coupling  force 
    K_c = (1 - p) * K;           # cross coupling  force

    # IVP
    R₀ = rand(P) * 0.3;
    φ₀ = mod2pi.(rand(P) * 2π);
    x₀ = @. R₀*cos(φ₀)
    y₀ = @. R₀*sin(φ₀)
    u₀ = [x₀; y₀];
    tspan = (0.0, Tend);    # time span
    global par = K_s, K_c, α;

    prob = ODEProblem(kuramoto_Ppop_carts!, u₀, tspan, par, dtmax=0.1, dtmin = 1e-5);
    sol = @time solve(prob, saveat=0.1);

    # Rmean computation
    x = sol[1:P,:]
    y = sol[P+1:end,:]
    R = @. sqrt(x^2 + y^2)
    Rmean = sqrt.(mean(x, dims=1) .^ 2 + mean(y, dims=1) .^ 2)
    global Rmean = vec(Rmean)


    # Interpolation
    itp = interpolate!(Rmean, BSpline(Quadratic(InPlace(OnCell()))))
    global Rmean_itp = [itp[ii] for ii in 2e4:.01:Int(tspan[2]/0.1)]

    # Local maxima computation
    pks, vals =  findmaxima(Rmean_itp, strict = true)
    
    # Filtering equal valued peaks
    vals_trunc = trunc.(vals, digits = 7)
    unique!(vals_trunc)
    plot!(  [K for ii in eachindex(vals_trunc)],
            vals_trunc, 
            lt = :scatter,
            markersize = 1,
            markerstrokewidth = 0,
            color = :steelblue2,
            label = ""
        )
    
# Terminal information of progress 
    itnumber +=1
    println("Iterations completed: "*string(itnumber)*"/"*total_its) 
end

#bifurcation diagram
title!("Bifurcation diagram, "*latexstring(P)*" populations\np = "*latexstring(p)*", "*L"\alpha"*" = "*latexstring(α), )
xlabel!(L"K")
ylabel!(L"R_{mean}")
display(p2)

name = "kuramoto_Mpops_bifurcation_"*string(p)*".png"
savefig(p2,name)
println("Bifurcation diagram saved in "*pwd()*" as "*name)