using SpecialFunctions
using Plots
using JLD


#------------------------------------RTA_RK4--------------------------------------------------------------
# Data/RTA_T_1.jld              --------- RTA Usual Collision kernel -- n ϵ {1,2,3}
nₙ  = 1
nₑ  = 2

T₀  = 1
m   = 0
μ₀  = 0

tₛ = 0.1
tₑ = 100
ξ  = 0.01

L = 10

η₀ = (tₛ*T₀)/ξ

println("Loading RTA Data.")

y = load("Data/RTA_T_1.jld")["y"]
t_RTA = y[1]
ρ_RTA = y[2]

ϵᵈ_RTA  = ρ_RTA[:,nₑ,1]
nᵈ_RTA  = ρ_RTA[:,nₙ,1]

T_RTA  = (1/3)*(ϵᵈ_RTA ./nᵈ_RTA )
τ_RTA  = t_RTA #((T_RTA .^(1+0.0)).*(t_RTA)./η₀)



#------------------------------------RTA_Dcol_RK4------------------------------------------------
# RTA_Dcol_T_1.jld              --------- RTA Denicol's Collision kernel, Relaxing Lowest moment.
nₐᵣ = [0,1,2,3]


nₙ  = 2
nₑ  = 3

T₀  = 1
Λ   = 1
m   = 0
μ₀  = 0

tₛ = 0.1
tₑ = 100
ξ  = 0.01

L = 10

Λ = 0.0
#--------------------------------
η₀ = (tₛ/ξ)*(T₀)
#--------------------------------


println("Loading Denicol RTA Data.")

y = load("Data/RTA_Dcol_T_1.jld")["y"]
t_DRTA = y[1]
ρ_DRTA = y[2]

ϵᵈ_DRTA  = ρ_DRTA[:,nₑ,1]
nᵈ_DRTA  = ρ_DRTA[:,nₙ,1]

T_DRTA  = (1/3)*(ϵᵈ_DRTA ./nᵈ_DRTA )
τ_DRTA  = t_DRTA #((T_DRTA.^(1+0.0)).*(t_DRTA)./η₀)



#----------------------------------------MDRTA_Simplified_RK4_Z_1----------------------------------------------------------
# Data/MDRTA_T_1_Λ_1_Z_1.jld      --------- MDRTA Λ = 1.0, T = 1, nₘᵢₙ =  1.0, Free streaming lowest moment n ϵ {1,2}
nₐᵣ = [1,2]


nₙ  = 1
nₑ  = 2

T₀  = 1
m   = 0
μ₀  = 0

L = 10


tₛ = 0.1
tₑ = 100
ξ  = 0.01

Λ = 1
#--------------------------------
η₀ = 10 #(tₛ/ξ)*((T₀)^(1+1))
#--------------------------------



println("Loading MDRTA Data. Λ = 1.0, T = 1, nₘᵢₙ =  1.0, Free streaming lowest moment n ϵ {1,2}")

y = load("Data/MDRTA_T_1_t0_$(tₛ)_Λ_1_ZZ.jld")["y"]
t_MDRTA_Z_1 = y[1]
ρ_MDRTA_Z_1 = y[2]

ϵᵈ_MDRTA_Z_1  = ρ_MDRTA_Z_1[:,nₑ,1]
nᵈ_MDRTA_Z_1 = ρ_MDRTA_Z_1[:,nₙ,1]

T_MDRTA_Z_1  = (1/3)*(ϵᵈ_MDRTA_Z_1 ./nᵈ_MDRTA_Z_1 )
τ_MDRTA_Z_1  = t_MDRTA_Z_1 # ((T_MDRTA_Z_1.^(1+1)).*(t_MDRTA_Z_1)./η₀)

#-----MDRTA_Simplified_RK4_ZZ_0.5------------
# Data/MDRTA_T_1_Λ_0.5_Z.jld    --------- MDRTA Λ = 0.5, T = 1, nₘᵢₙ =  0.5, Free streaming lowest moment n ϵ {0.5,1,1.5,2}

nₐᵣ = [1,1.5,2]


nₙ  = 1
nₑ  = 3

T₀  = 1
Λ   = 1
m   = 0
μ₀  = 0

L = 10


tₛ = 0.1
tₑ = 100
ξ  = 0.01


Λ = 0.5
#--------------------------------
η₀ = 10 # (tₛ/ξ)*((T₀)^(1+0.5))
#--------------------------------

println("Loading MDRTA Data. Λ = 0.5, T = 1, nₘᵢₙ =  0.5, Free streaming lowest moment n ϵ {0.5,1,1.5,2}")

y = load("Data/MDRTA_T_1_t0_$(tₛ)_Λ_0.5_ZZ.jld")["y"]
t_MDRTA_Z_5 = y[1]
ρ_MDRTA_Z_5 = y[2]


ϵᵈ_MDRTA_Z_5  = ρ_MDRTA_Z_5[:,nₑ,1]
nᵈ_MDRTA_Z_5 = ρ_MDRTA_Z_5[:,nₙ,1]

T_MDRTA_Z_5  = (1/3)*(ϵᵈ_MDRTA_Z_5 ./nᵈ_MDRTA_Z_5 )
τ_MDRTA_Z_5  = t_MDRTA_Z_5 #((T_MDRTA_Z_5.^(1+0.5)).*(t_MDRTA_Z_5)./η₀)

println("\n----------------------\n Loading Finished \n----------------------\n")


#-----MDRTA_Simplified_RK4_ZZ_0.25------------
# Data/MDRTA_T_1_Λ_0.5_Z.jld    --------- MDRTA Λ = 0.5, T = 1, nₘᵢₙ =  0.5, Free streaming lowest moment n ϵ {1,1.25,1.5,1.75,2}

nₐᵣ = [1,1.25,1.5,1.75,2]


nₙ  = 1
nₑ  = 5

T₀  = 1
Λ   = 1
m   = 0
μ₀  = 0

L = 10


tₛ = 0.1
tₑ = 100
ξ  = 0.01


Λ = 0.25
#--------------------------------
η₀ = 10 # (tₛ/ξ)*((T₀)^(1+0.5))
#--------------------------------

println("Loading MDRTA Data. Λ = 0.5, T = 1, nₘᵢₙ =  0.5, Free streaming lowest moment n ϵ {0.5,1,1.5,2}")

y = load("Data/MDRTA_T_1_t0_$(tₛ)_Λ_0.25_ZZ.jld")["y"]
t_MDRTA_Z_25 = y[1]
ρ_MDRTA_Z_25 = y[2]


ϵᵈ_MDRTA_Z_25  = ρ_MDRTA_Z_25[:,nₑ,1]
nᵈ_MDRTA_Z_25 = ρ_MDRTA_Z_25[:,nₙ,1]

T_MDRTA_Z_25  = (1/3)*(ϵᵈ_MDRTA_Z_25 ./nᵈ_MDRTA_Z_25 )
τ_MDRTA_Z_25  = t_MDRTA_Z_25 #((T_MDRTA_Z_5.^(1+0.5)).*(t_MDRTA_Z_5)./η₀)

println("\n----------------------\n Loading Finished \n----------------------\n")



#-----------------Plotting

plot(τ_RTA,T_RTA, xaxis=:log, xlabel="τ", ylabel="T",label="RTA",dpi=300)
plot!(τ_DRTA,T_DRTA, xaxis=:log, label="Λ = 0")



plot!(τ_MDRTA_Z_1,T_MDRTA_Z_1, xaxis=:log, label="Λ = 1")
plot!(τ_MDRTA_Z_5,T_MDRTA_Z_5, xaxis=:log, label="Λ = 0.5")
plot!(τ_MDRTA_Z_25,T_MDRTA_Z_25, xaxis=:log, label="Λ = 0.25")

savefig("Plots/Temperature_Plots.png")

