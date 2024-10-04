using DifferentialEquations
using LinearAlgebra
using SpecialFunctions
using Plots
using JLD






println("Loading RTA Data.")

y = load("Data/MDRTA_FCK_T_1_Lam_1_n_0_3_L_150.jld")["y"]



t     = y[1]
ρ     = y[2]
nₐᵣ   = y[3]


nₙ = findall(x->x ==1,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==2,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

ϵᵈ  = ρ[:,nₑ,1]
nᵈ  = ρ[:,nₙ,1]

T  = (1/3)*(ϵᵈ ./nᵈ )



L = 150


plot(t,T, xaxis=:log, xlabel="τ", ylabel="T",label="Λ = 1",dpi=300)




println("Loading RTA Data.")

y = load("Data/MDRTA_FCK_T_1_Lam_0.5_n_0_3_L_150.jld")["y"]



t     = y[1]
ρ     = y[2]
nₐᵣ   = y[3]


nₙ = findall(x->x ==1,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==2,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

ϵᵈ  = ρ[:,nₑ,1]
nᵈ  = ρ[:,nₙ,1]

T  = (1/3)*(ϵᵈ ./nᵈ )



L = 150


plot!(t,T, xaxis=:log, xlabel="τ", ylabel="T",label="Λ = 0.5",dpi=300)


println("Loading RTA Data.")

y = load("Data/MDRTA_FCK_T_1_Lam_0.25_n_0_3_L_150.jld")["y"]



t     = y[1]
ρ     = y[2]
nₐᵣ   = y[3]


nₙ = findall(x->x ==1,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==2,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

ϵᵈ  = ρ[:,nₑ,1]
nᵈ  = ρ[:,nₙ,1]

T  = (1/3)*(ϵᵈ ./nᵈ )



L = 150


plot!(t,T, xaxis=:log, xlabel="τ", ylabel="T",label="Λ = 0.25",dpi=300)


println("Loading RTA Data.")

y = load("Data/MDRTA_FCK_T_1_Lam_0_n_0_3_L_150.jld")["y"]



t     = y[1]
ρ     = y[2]
nₐᵣ   = y[3]


nₙ = findall(x->x ==1,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==2,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

ϵᵈ  = ρ[:,nₑ,1]
nᵈ  = ρ[:,nₙ,1]

T  = (1/3)*(ϵᵈ ./nᵈ )



L = 150


plot!(t,T, xaxis=:log, xlabel="τ", ylabel="T",label="Λ = 0",dpi=300)



savefig("Plots/The_Problem/Temp_MDRTA.pdf")


