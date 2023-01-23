using Plots, LinearAlgebra, Contour
Contour
import Contour: contours, levels, level, lines, coordinates

function kuramoto_eigenvalues_symmetric_J1(K_c, K_s)
    R = sqrt((-2*Δ + K_s*cos(γ_s) + K_c*cos(γ_c))/(K_s*cos(γ_s) + K_c*cos(γ_c)))

    a = -Δ + 0.5*(K_s*cos(γ_s) + K_c*cos(γ_c))*(1-3*R^2)
    b = -K_c*R^3*sin(γ_c)
    c = R*(K_s*sin(γ_s) + K_c*sin(γ_c))
    d = 0

    J₁ = [a  b; 
      c  d]

    tr_J₁ = a 
    det_J₁ = - b*c

    vaps_J₁ = 0.5 * (tr_J₁ + sqrt(Complex(tr_J₁^2 - 4*det_J₁))); 
    re_vaps_J₁ = real(vaps_J₁);
    return re_vaps_J₁
end

function kuramoto_eigenvalues_symmetric_J2(K_c, K_s)
    R = sqrt((-2*Δ + K_s*cos(γ_s) + K_c*cos(γ_c))/(K_s*cos(γ_s) + K_c*cos(γ_c)))

   # e = -Δ + 0.5*(K_s*cos(γ_s)*(1-3*R^2) - K_c*cos(γ_c)*(1+R^2))
    #f = K_c*R*sin(γ_c)
    #g = K_s*R*sin(γ_s) - K_c/R*sin(γ_c)
    #h = - K_c*cos(γ_c)*(1+R^2)

    e = -Δ + 0.5*K_s*cos(γ_s) - 3/2*K_s*R^2*cos(γ_s) - K_c*R^2*cos(γ_c) - 0.5*K_c*cos(γ_c) + 0.5*K_c*R^2*cos(γ_c)
    f = K_c*R*sin(γ_c) - K_c*R^3*sin(γ_c)
    g = K_s*R*sin(γ_s) - K_c/R*sin(γ_c)
    h = - K_c*cos(γ_c)*(1+R^2)
    
    J₂ = [e  f; 
      g  h]

    tr_J₂ = e + h
    det_J₂ = e*h - g*f

    vaps_J₂ = 0.5 * (tr_J₂ + sqrt(Complex(tr_J₂^2 - 4*det_J₂)))
    re_vaps_J₂ = real(vaps_J₂)
    return re_vaps_J₂
end


Δ = 1
γ_s = γ_c = 1.2


K_s = K_c = range(5,25, length = 1000)

#contour(K_s, K_c,kuramoto_eigenvalues_symmetric_J1, levels = [0], ratio = 1)

R(K_s::Float64) = sqrt((-2*Δ + K_s*cos(γ_s) + K_c*cos(γ_c))/(K_s*cos(γ_s) + K_c*cos(γ_c)))

Rvec = [0.0]
λvec_sym = [0.0]
λvec_anti = [0.0]
for K_s = 4.8:.01:25
    R_ks = R(K_s)
    λ_ks_sym =  kuramoto_eigenvalues_symmetric_J1(K_s, 0.8)
    λ_ks_anti=  kuramoto_eigenvalues_symmetric_J2(K_s, 0.8)
    push!(Rvec, R_ks)
    push!(λvec_sym, λ_ks_sym)
    push!(λvec_anti, λ_ks_anti)
end

p1 = plot(4.8:.01:25, Rvec[2:end], xlabel= K_s, ylabel = "R_a")

p2 = plot(4.8:.01:25, λvec_sym[2:end])
     plot!(4.8:0.01:25, λvec_anti[2:end])

plot(p1,p2, layout = (2,1))

# K_s = 6.164

## CONTOUR PLOT
K_s = 5:.5:8; K_c = 0.1:1.5


K_s = range(5.6,8, length = 100)
K_c = range(0.1,1.5, length = 100)

cont = contour(K_c,K_s,kuramoto_eigenvalues_symmetric_J2, 0.0)
