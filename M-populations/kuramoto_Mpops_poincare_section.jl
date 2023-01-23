# POINCARE SECTION OF A SYSYEM OF M POPULATIONS OF KURAMOTO-SAKAGUCHI OSCILLATORS 
using Plots, DifferentialEquations, Distributions, Interpolations, Peaks, Statistics, LaTeXStrings

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
# Specify p 
p = 0.7;                 # coupling ratio

# Specify range of K for which you want to plot the Poincaré section:
K_first = 17.1;
K_last  = 18.0; 
steplength = 0.1;
Kvec = collect(K_first:steplength:K_last);

peaks_vec = [];
iter = 0;
for K in Kvec
    global M = 500;                 # number of populations
    global Δ = 1.0;                 # disorder parameter
    K_s = p*K;                      # self  coupling  force 
    K_c = (1-p)*K;                  # cross coupling  force
    α = 1.2;                        # sakaguchi phase shift 
    Tend = 15000.0

    R₀ = rand(M)*0.3;        # initial conditions for IVP
    φ₀ = rand(M)*2π; φ₀ = mod2pi.(φ₀)
    x₀ = @. R₀*cos(φ₀)
    y₀ = @. R₀*sin(φ₀)
    u₀ = [x₀; y₀];

    tspan = (0.0,Tend)    # time span
    global par = K_s, K_c, α;

    prob = ODEProblem(kuramoto_Ppop_carts!, u₀, tspan,par, dtmax = 0.1, dtmin = 1e-5);
    sol = @time solve(prob, saveat = 0.1);

    x = sol[1:M,:];             
    y = sol[M+1:end,:];

    global Rmean = vec(sqrt.(mean(x,dims=1).^2+mean(y,dims=1).^2));   # mean modulus of order parameter

    itp = interpolate!(Rmean, BSpline(Quadratic(InPlace(OnCell())))); # interpolation of Rmean
    global Rmean_itp = [itp[ii] for ii in 2e4:.01:Int(tspan[2]/0.1)]; 

    pks, vals =  findmaxima(Rmean_itp, strict = true);                # local maxima of Rmean
    

    push!(peaks_vec, vals)
    iter +=1

    p1 = plot( peaks_vec[iter][2:end], peaks_vec[iter][1:end-1],
        lt = :scatter,
        markersize = 1,
        markerstrokewidth = 0,
        color = :steelblue2,
        xlabel = L"Z_{n+1}",
        ylabel = L"Z_{n}",
        alpha = 0.7,
        label = "K = "*string(Kvec[iter]),
        background_color_legend = nothing,
        foreground_color_legend = nothing
    ) 
    savefig(p1,"kuramoto_poinc_"*string(Kvec[iter])*".png")
    println("Iterations completed = "*string(iter)*"/"*string(length(Kvec)))
end


using Alert
alert("DONE!")



# iter = 0
# for ii in 13.1:.1:15
#     iter +=1
#     p1 = plot( peaks_vec[iter][2:end], peaks_vec[iter][1:end-1],
#     lt = :scatter,
#     markersize = 1,
#     markerstrokewidth = 0,
#     color = :steelblue2,
#     xlabel = L"Z_{n+1}",
#     ylabel = L"Z_{n}",
#     alpha = 0.6,
#     label = "K = "*string(Kvec[iter]),
#     background_color_legend = nothing,
#     foreground_color_legend = nothing
#     ) 
#     savefig(p1,"kuramoto_poinc_"*string(Kvec[iter])*".png")
# end
