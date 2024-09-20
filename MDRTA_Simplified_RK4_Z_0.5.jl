using LinearAlgebra
using SpecialFunctions
using Plots
using JLD


function trange(tspan,N,interpolation)
    if interpolation == "log"
        tₛ       = log(tspan[1])
        tₑ      = log(tspan[2])
        step    = (tₑ -  tₛ)/N
        return [exp(tₛ + i*step) for i in 1:N ]

    elseif interpolation == "lin"
        tₛ       = tspan[1]
        tₑ      = tspan[2]
        step    = (tₑ -  tₛ)/N
        return range(tsapn[2],tspan[1],step = step)
    else
        println("Invalid Interpolation input.") 
    end    

    
end


function RK4(y₀::Matrix,ts::Vector,fun,p = nothing)
    print("Array RK4 Solver Called. \n")
    
    ysize = (size(ts)...,size(y₀)...)  # Get the size of the time series amtrix
    y  = Array{Float64}(undef,ysize )  # Create the time series matrix
    dy = similar(y₀)                   # Create matrix to store dy
    view(y,1,:,:) .= y₀                # Initialise t = tᵢ value of y to y₀
    
    # Checks if the function takes any extra parameters.
    if isnothing(p) 
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:,:)              # Creates a view of the Matrix of variables

            # Computes the RK4 parameters
            k₁ = fun(dy, yᵢ,t)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt)
            k₄ = fun(dy, yᵢ .+     dt*k₃, t +     dt)

            view(y,i+1,:,:) .= yᵢ .+ (dt*( k₁ .+ 2k₂ .+ 2k₃ .+ k₄ )/6)

        end

    else
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:,:)             # Creates a view of the Matrix of variables

            # Computes the RK4 parameters
            k₁ = fun(dy, yᵢ,t,p)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt,p)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt,p)
            k₄ = fun(dy, yᵢ .+     dt*k₃, t +     dt,p)

            view(y,i+1,:,:) .= yᵢ .+ (dt*( k₁ .+ 2k₂ .+ 2k₃ .+ k₄ )/6)

        end
    end
    return y
    
end


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
function P(n,l) 
    return (n+2l)*( (2l*(2l-1))/((4l+1)*(4l-1)) )   
end

function Q(n,l) 
    return ((2/3 ) + n*( (8*l^2 + 4l -1)/((4l-1)*(4l+3)) ) + ( (2l*(2l+1))/(3*(4l-1)*(4l+3)))) 
end

function R(n,l) 
    return (n-2l-1)*( ((2l+1)*(2l+2))/((4l+1)*(4l+3)) )     
end


function ρeq(T,μ,n,l)
    if l == 0
        return (gamma(n+2)*((T^(n+2))/(2*π^2))*exp(μ/T))/(2*l+1)        
    else
        return 0 
    end
end


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

function MDRTA_S(dρ::Matrix,ρ,t::Float64,p)
    N   = 3
    L   = p[2]
    η₀   = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]

    nᵈ = ρ[nₙ,1]
    ϵᵈ = ρ[nₑ,1]

    # Computing Temperature and Chemical potential
    T = (1/3)*(ϵᵈ/nᵈ)    
    μ = T*log( 27*(π^2)*nᵈ*((nᵈ/ϵᵈ)^3))

    # The reciprocal of relaxation time
    Λ = 0.5
    ωᵣ = (T^(1+0.5))/η₀


    

    #Free Streaming of n = 0.
    # n = 0, l = 0
    dρ[1,1]     = -(1/t)*(                      Q(0.5,0)*ρ[1,1]   + R(0.5,0)*ρ[1,1+1]   ) 
    for j in 2:L-1
        # n = 0, l = j -1
        dρ[1,j] = -(1/t)*(  P(0.5,j-1)*ρ[1,j-1] + Q(0.5,j-1)*ρ[1,j] + R(0.5,j-1)*ρ[1,j+1]   )
    end
    # n = 0, l = L -1
    dρ[1,L]     = -(1/t)*(  P(0.5,L-1)*ρ[1,L-1] +  Q(0.5,L-1)*ρ[1,L]                    )

    for (i,nv) in enumerate(nₐᵣ[begin+1:end])
        n = i+1
        # Generalised Collision kernal

        Cₙₗ = -ωᵣ*( (ρ[n-1,1]*ρeq(T,μ,1 - Λ,0) -ρ[nₙ-1,1]* ρeq(T,μ,nv- Λ,0) )*(  -ρeq(T,μ,1 - Λ,0)*(ρeq(T,μ,2 - Λ,0)^2) +  ρeq(T,μ,3 - Λ,0)*(ρeq(T,μ,1 - Λ,0)^2) ) - 
        ρeq(T,μ,1 - Λ,0)*( ρ[nₙ-1,1]*ρeq(T,μ,2 - Λ,0) -  ρ[nₑ-1,1]*ρeq(T,μ,1 - Λ,0)  )*(  ρeq(T,μ,nv - Λ,0)*ρeq(T,μ,2 - Λ,0)   - ρeq(T,μ,1 - Λ,0)*ρeq(T,μ,nv+1 - Λ,0) )  )/( ρeq(T,μ,1 - Λ,0)*( -ρeq(T,μ,1 - Λ,0)*(ρeq(T,μ,2 - Λ,0)^2) +  ρeq(T,μ,3 - Λ,0)*(ρeq(T,μ,1 - Λ,0)^2)) )

        """if nv == 2
            println("Cₙₗ E : ",Cₙₗ)
        elseif nv == 1
            println("Cₙₗ N : ",Cₙₗ)
        end"""
        # Collision kernal zero for number and energy desnity, n = 1,2 (Array index 2,3) and l = 0 (Array index 1).
        # n, l = 0
        dρ[n,1] = -(1/t)*(                      Q(nv,0)*ρ[n,1] + R(nv,0)*ρ[n,2]    ) + Cₙₗ 
        # Running through all l > 0 moments.
        for j in 2:L-1 
            # n = i-1, l = j-1                                                                     # Collision dependance.
            dρ[n,j] = -(1/t)*( P(nv,j-1)*ρ[n,j-1] +  Q(nv,j-1)*ρ[n,j] + R(nv,j-1)*ρ[n,j+1]    )  - ωᵣ*ρ[n-1,j] 
        end
        # Truncation for l = L moment.
        # n = i-1, l = L-1
        dρ[n,L]     = -(1/t)*( P(nv,L-1)*ρ[n,L-1] +  Q(nv,L-1)*ρ[n,L]                      )  - ωᵣ*ρ[n-1,L]
    end

    return dρ
end

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------





nₐᵣ = [0.5,1,1.5,2]


nₙ  = 2
nₑ  = 4

T₀  = 1
Λ   = 1
m   = 0
μ₀  = 0


N = size(nₐᵣ)[1]
L = 10


tₛ = 0.1
tₑ = 200
ξ  = 0.01


Λ = 0.5
#--------------------------------
η₀ = (tₛ/ξ)*((T₀)^(1+0.5))
#--------------------------------

tspan = (tₛ,tₑ)

ρ₀ = Matrix{Float64}(undef, N, L)

for (n,nv) in enumerate(nₐᵣ)    
    for l in 0:L-1        
        ρ₀[n,l+1] = ρeq(T₀,μ₀,nv,l)
    end
end

ϵᵈ₀ = ρ₀[nₑ,1]
nᵈ₀ = ρ₀[nₙ,1]

println(" MDRTA Λ = 0.5 \n Free streaming lowest moment n ϵ {0.5,1,1.5,2} \n")

η = η₀*T₀^(0.5)
println("m     : ", m)
println("T₀    : ",(1/3)*(ϵᵈ₀/nᵈ₀))
println("t₀    : ",tₛ)
println("η₀/s₀ : ",η)


p = (N,L,η₀,nₐᵣ,nₙ,nₑ)

tspan = trange((tₛ,tₑ),100,"log")


#------------------------------------------------

ρₜ = RK4(ρ₀,tspan,MDRTA_S,p)

#------------------------------------------------

ϵᵈ = ρₜ[:,nₑ,1]
nᵈ = ρₜ[:,nₙ,1]

T = (1/3)*(ϵᵈ./nᵈ)
#τ = ((T.^(1+0.5)).*(tspan)./η₀)
τ = tspan
println("Saving data.")

y = (tspan,ρₜ)
save("MDRTA_T_1_Λ_0.5_Z.jld","y",y)


plot(τ,T, xaxis=:log,xlabel="τ", ylabel="T",label="η = $(η), Λ = 0.5",dpi=300)
#plot(tspan,log.(nᵈ), xaxis=:log)
#plot(tspan,T, xaxis=:log)