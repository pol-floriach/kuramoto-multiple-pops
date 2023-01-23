using DifferentialEquations, Plots, LinearAlgebra, Alert, LaTeXStrings, Statistics


# Mean field of a system of M coupled populations of N 
# Kuramoto-Sakaguchi oscillators, 2M equations using the Ott-Antonsen ansatz
# Note: the equations are on cartesian form
function kuramoto_Mpop_carts!(du,u,p,t)
       K_s, K_c, α = p        
   
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

# Computation parameters
M = 500;                 # number of populations
Δ = 1.0;                 # disorder parameter
p = 0.9;                 # coupling ratio
K = 7.3;
                         # coupling strength
K_s = p*K;               # self  coupling  force 
K_c = (1-p)*K;           # cross coupling  force
α = 1.2;                 # sakaguchi phase shift 


# initial conditions
R₀ = rand(M)*0.3;  
φ₀ = rand(M)*2π; φ₀ = mod2pi.(φ₀)

x₀ = @. R₀*cos(φ₀)
y₀ = @. R₀*sin(φ₀)
u₀ = [x₀; y₀];
tspan = (0.0,10000.0)    # time span
par = K_s, K_c, α;

prob = ODEProblem(kuramoto_Mpop_carts!, u₀2, tspan,par, dtmax = 0.1, dtmin = 1e-5);
sol = @time solve(prob, saveat = 0.1);


x = sol2[1:M,:];
y = sol2[M+1:end,:];
R = @. sqrt(x^2 + y^2)
Rmean=sqrt.(mean(x,dims=1).^2+mean(y,dims=1).^2)
Rmean = vec(Rmean)

# If you want to plot the mean phase
#Φ = atan.(mean(y,dims=1),mean(x,dims=1))'

# Time series plot
p = plot(  sol2.t,R[1,:], 
        title = "Time series \n"*"K  = "*string(K)*", p = "*string(p), 
        label = "osc. 1, cartesian"
        )
plot!(  sol2.t,R[38,:], 
    label = "osc. 2, cartesian")
plot!(  sol2.t,Rmean2,
    xlim=(1e3,4e3),
    lw=2, 
    label = "Mean")


  
# Use this part of the code to obtain a heatmap, with time on the x axis, 
# oscillator number on the y axis and the order parameter on the 
# normal direction of the screen
# p3 = heatmap(sol.t,collect(1:M),R, 
#               colorbar_title = " R", 
#               c = :Accent_3,
#               title = "Local R \n"*"K = "*string(K)*", p = "*string(p)*", α = "*string(α) 
#          )



# Use this part of the code if you want to create a gif animation of the time series. 
# anim = @animate for ii ∈ 9000:10000
#     ii = 10000
#   plot!(  x1[:,ii],y1[:,ii],
#             lt=:scatter,
#             title = "Local order parameter \n t = "*string(round(0.1*ii)),
#             xlim = (-1,1), 
#             ylim = (-1,1), 
#             label = false)   
#       plot!([cos(θ) for θ in 0:.1:2π],[sin(θ) for θ in 0:.1:2π],
#              ls = :dash, 
#              color = :black, 
#              label = false)
#       plot!([0],[0], 
#              lt = :scatter, 
#              color = :black, 
#              label = false)
             
# end
# gif(anim)




# Use this part of the code to obtain the lyapunov exponents of the system
# Note: computing more than the maximum lyapunov exponent takes a really long time, 
# I would recommend using another code 
#using DynamicalSystems
#dynsys = ContinuousDynamicalSystem(prob)
# Maximum lyapunov exponent
#λ = @time lyapunov(dynsys, 10000)
# Lyapunov spectrum
#λ_vec = @time lyapunovspectrum(dynsys,10000,1; show_progress = true)



