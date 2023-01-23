using DifferentialEquations, LinearAlgebra, Distributions, LaTeXStrings, Plots, Alert

# ODE system
function kuramoto_ring!(dϕ,ϕ,_,t)
    @. ϕ = mod2pi(ϕ)
    sinvec = @. sin(ϕ)
    cosvec = @. cos(ϕ)
    # més allocations 
    # Gdotcosvec = G'*cosvec
    # Gdotsinvec = G'*sinvec

    # sumatorivec = @. (sinvec*cosα-sinα*cosvec)*Gdotcosvec - (cosvec*cosα + sinvec*sinα)*Gdotsinvec
    # #sumatorivec = @. (cosvec*cosα+sinα*sinvec)*Gdotcosvec + (sinvec*cosα - cosvec*sinα)*Gdotsinvec
    # @. dϕ = ωi - K*2π/N*sumatorivec
    mul!(Gdotcosvec, G',cosvec)
    mul!(Gdotsinvec, G,sinvec)
    @. dϕ = ωi - K*2π/N* ((sinvec*cosα-sinα*cosvec)*Gdotcosvec - (cosvec*cosα + sinvec*sinα)*Gdotsinvec)
end

# callback to save local order parameter
function save_R(u,t,integrator)
    @. u = mod2pi(u)
    
    R_k = ((1/N*sum(sin.(u)))^2 + (1/N*sum(cos.(u)))^2)^(.5);
    return R_k
end

function save_R_local(u,t,integrator)
    @.u = mod2pi(u)
     R_k = 2π/N*G'*cos.(u) .+ 2π/N*G'*sin.(u)*im 
     R_k = real.(abs.(R_k))
     return R_k
end
saved_Rl = SavedValues(Float64,Array{Float64}) # local
saved_R = SavedValues(Float64,Float64)         # global

Tend = 1000
cb_Rl = SavingCallback(save_R_local,saved_Rl,saveat=0:1:Tend);
cb_R = SavingCallback(save_R,saved_R,saveat=0:.1:Tend);

N = 1000;

α = 1.5
cosα = cos(α);
sinα = sin(α);

A = 0.95

K = 4.0

G = zeros(N,N); # coupling matrix
for ii = 1:N, jj = 1:N
    x = 2π*abs(ii-jj)/N
    G[ii,jj] = (1 + A*cos(x))/2π    
end

global  G = Symmetric(G);

ωi = mod2pi.(rand(Cauchy(0,0.015), N)); 
ϕ₀ = mod2pi.(rand(Cauchy(0,1), N));
    #= per estat quimera
        ωi = zeros(N)
        ϕ₀ = zeros(N)
        for ii = 1:N
        ϕ₀[ii] = 6*exp(-0.76*(2π/N*ii)^2)*(rand(1)[1]-1)/2
        end =#

Gdotsinvec =  zeros(N)
Gdotcosvec =  zeros(N)
tspan = (0.0, Tend);
prob = ODEProblem(kuramoto_ring!, ϕ₀, tspan);

sol = @time solve(prob, saveat=0.1, callback = cb_Rl, dtmin = 1e-5); println("~>Simulation completed") # local


R_loc = reduce(hcat, saved_Rl.saveval); println("~>Order Parameter extracted");
# local
#R_k = saved_R.saveval;                 # global


#name = "K_"*string(K)*"-A_"*string(A)*"-alph_"*string(α)*".png"
#fig = plot(sol.t, R_k, title = "Global order parameter \n A ="*string(A)*", K = "*string(K)*", α = "*string(α), label = "R")

## HEATMAP  
x = 2π/N*collect(1:N);  

p1 = heatmap(sol.t[1:10:end],x,R_loc, colorbar_title = "Rₖ", c = :Accent_3);
xaxis!("time");
yaxis!("x");
title!("Time evolution of local order parameter Rₖ, \n "*"A = "*string(A)*", K = "*string(K)*", α = "*string(α)*"\n"); println("~>Heatmap done"); 
p1
#savefig(p1, "K_"*string(K)*"_alph_"*string(α)*".png"); alert("Heatmap!!!");


R_loc_mean = zeros(N);
for ii in 1:N
    R_loc_mean[ii] = sum(R_loc[ii,:])/length(R_loc[ii,:])
end
p2 = plot(x,R_loc_mean, xlim = (0,2π), ylim = (0.0, 1.0))
fig = plot(p1,p2, layout = (2,1))
alert("DONE!")


p1
N = 2000;

G = zeros(N,N); # coupling matrix
for ii = 1:N, jj = 1:N
    x = 2π*abs(ii-jj)/N
    G[ii,jj] = (1 + A*cos(x))/2π    
end

global  G = Symmetric(G);

ωi = mod2pi.(rand(Cauchy(0,0.015), N)); 
ϕ₀ = mod2pi.(rand(Cauchy(0,1), N));
    #= per estat quimera
        ωi = zeros(N)
        ϕ₀ = zeros(N)
        for ii = 1:N
        ϕ₀[ii] = 6*exp(-0.76*(2π/N*ii)^2)*(rand(1)[1]-1)/2
        end =#


tspan = (0.0, Tend);
prob = ODEProblem(kuramoto_ring!, ϕ₀, tspan);

sol = @time solve(prob, saveat=0.1, callback = cb_Rl); println("~>Simulation completed") # local
#sol = @time solve(prob, saveat=0.1, callback = cb_R);   # global


R_loc = reduce(hcat, saved_Rl.saveval); println("~>Order Parameter extracted");
# local
#R_k = saved_R.saveval;                 # global


#name = "K_"*string(K)*"-A_"*string(A)*"-alph_"*string(α)*".png"
#fig = plot(sol.t, R_k, title = "Global order parameter \n A ="*string(A)*", K = "*string(K)*", α = "*string(α), label = "R")

## HEATMAP  
x = 2π/N*collect(1:N);  

p3 = heatmap(sol.t[1:10:end],x,R_loc, colorbar_title = "Rₖ", c = :Accent_3);
xaxis!("time");
yaxis!("x");
title!("Time evolution of local order parameter Rₖ, \n "*"A = "*string(A)*", K = "*string(K)*", α = "*string(α)*"\n"); println("~>Heatmap done"); 
p3
#savefig(p1, "K_"*string(K)*"_alph_"*string(α)*".png"); alert("Heatmap!!!");




fig = plot(p1, p3, layout = (2,1))
alert("ALl DONE")
# a = 1
# b = 1000
# p3= heatmap(sol.t[end-1000:end],collect(a:b),reduce(hcat,mod2pi.(sol).u[end-1000:end])[a:b,:], colorbar_title = "ϕ_i", c = :Blues)
#     xlabel!("time")
#     ylabel!("Oscillator index")
# title!("Time series of a ring of "*string(N)*" oscillators. \n"*"A = "*string(A)*", K = "*string(K)*", α = "*string(α)*"\n");

# fig = plot(p1,p3, layout = (2,1))

# savefig(fig,name)

