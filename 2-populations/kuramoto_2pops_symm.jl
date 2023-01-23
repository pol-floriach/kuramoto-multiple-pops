using DifferentialEquations, Plots


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


# SYMMETRIC

K_s = 8
K_c = 1.5
γ_s = 1.2
γ_c = 1.1


u₀ = [0.3; 0.4; 3]
tspan = (0,1000)
prob = ODEProblem(kuramoto_2pop_macro_reduced,u₀, tspan, dtmax = 0.01)
sol = solve(prob, saveat = 0.1)

p1 = plot(sol.t, sol[1,:])
    plot!(sol.t, sol[2,:])
p2 = plot(sol.t, mod2pi.(sol[3,:]))

R_analytic = sqrt((-2*Δ + K_s*cos(γ_s) + K_c*cos(γ_c))/(K_s*cos(γ_s) + K_c*cos(γ_c)))
R_numeric = sum(sol[1,300:end]) / length(sol[1,300:end])


# ANTIPHAASE


K_s = 8
K_c = 0.4
γ_s = 1.2
γ_c = 1.1


u₀ = [0.3; 0.4; 3]
tspan = (0,1000)
prob = ODEProblem(kuramoto_2pop_macro_reduced,u₀, tspan, dtmax = 0.01)
sol = solve(prob, saveat = 0.1)

p1 = plot(sol.t, sol[1,:])
    plot!(sol.t, sol[2,:])
p2 = plot(sol.t, mod2pi.(sol[3,:]))

R_analytic = sqrt((2*Δ - K_s*cos(γ_s) + K_c*cos(γ_c))/(-K_s*cos(γ_s) + K_c*cos(γ_c)))
R_numeric = sum(sol[1,300:end]) / length(sol[1,300:end])

p1
xlabel!("t")
ylabel!("R")

