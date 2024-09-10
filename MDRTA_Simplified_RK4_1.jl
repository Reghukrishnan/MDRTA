using LinearAlgebra
using SpecialFunctions
using Plots


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






function RTA(dρ::Matrix,ρ,t::Float64,p)
    N   = p[1]
    L   = p[2]
    η₀   = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]

    Nd = ρ[nₙ,1]
    Ed = ρ[nₑ,1]
    T = (1/3)*(Ed/Nd)
    #print(T)
    μ = T*log( 27*(π^2)*Nd*((Nd/Ed)^3))
    

    for (n,nv) in enumerate(nₐᵣ)
        for l =0:L-1
            if l == 0
                dρ[n,l+1] = - (0  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2] + (T/η₀)*t*(ρ[n,l+1] - ρeq(T,μ,nv,l) ) )/t
                #if nv ==1 || nv == 2
                #    print((ρ[n,l+1] - ρeq(T,μ,nv,l) ),"\n")
                #end
            elseif l == L-1
                dρ[n,l+1] = - (P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] + 0 + (T/η₀)*t*(ρ[n,l+1] - ρeq(T,μ,nv,l)  ) )/t
            else
                dρ[n,l+1] = - (P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2] + (T/η₀)*t*(ρ[n,l+1] - ρeq(T,μ,nv,l)  ) )/t
                
            end
        end
    end
    return dρ
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

    nᵈ = ρ[2,1]
    ϵᵈ = ρ[3,1]

    # Computing Temperature and Chemical potential
    T = (1/3)*(ϵᵈ/nᵈ)    
    μ = T*log( 27*(π^2)*nᵈ*((nᵈ/ϵᵈ)^3))

    # The reciprocal of relaxation time
    ωᵣ = (T^2)/η₀

    #Free Streaming of n = 0.
    # n = 0, l = 0
    dρ[1,1]     = -(1/t)*(                      Q(0,0)*ρ[1,1]   + R(0,0)*ρ[1,1+1]   ) 
    for j in 2:L-1
        # n = 0, l = j -1
        dρ[1,j] = -(1/t)*(  P(0,j-1)*ρ[1,j-1] + Q(0,j-1)*ρ[1,j] + R(0,j-1)*ρ[1,j+1]   )
    end
    # n = 0, l = L -1
    dρ[1,L]     = -(1/t)*(  P(0,L-1)*ρ[1,L-1] +  Q(0,L-1)*ρ[1,L]                    )

    for i in 2:3
            # Collision kernal zero for number and energy desnity, n = 1,2 (Array index 2,3) and l = 0 (Array index 1).
            # n = 1,2, l = 0
            dρ[i,1] = -(1/t)*(                      Q(i-1,0)*ρ[i,1] + R(i-1,0)*ρ[i,2]    )
        # Running through all l > 0 moments.
        for j in 2:L-1 
            # n = i-1, l = j-1                                                                     # Collision dependance.
            dρ[i,j] = -(1/t)*( P(i-1,j-1)*ρ[i,j-1] +  Q(i-1,j-1)*ρ[i,j] + R(i-1,j-1)*ρ[i,j+1]    )  - ωᵣ*ρ[i-1,j] 
        end
        # Truncation for l = L moment.
        # n = i-1, l = L-1
        dρ[i,L]     = -(1/t)*( P(i-1,L-1)*ρ[i,L-1] +  Q(i-1,L-1)*ρ[i,L]                      )  - ωᵣ*ρ[i-1,L]
    end

    return dρ
end

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------



function MDRTA_ST(dρ::Matrix,ρ,t::Float64,p)
    N   = 3
    L   = p[2]
    η₀   = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]

    nᵈ = ρ[2,1]
    ϵᵈ = ρ[3,1]

    # Computing Temperature and Chemical potential
    T = (1/3)*(ϵᵈ/nᵈ)    
    μ = T*log( 27*(π^2)*nᵈ*((nᵈ/ϵᵈ)^3))

    # The reciprocal of relaxation time
    ωᵣ = (T^(1+1))/η₀

    #Free Streaming of n = 0.
    # n = 0, l = 0
    dρ[1,1]     = -(1/t)*(                      Q(0,0)*ρ[1,1]   + R(0,0)*ρ[1,1+1]   ) 
    for j in 2:L-1
        # n = 0, l = j -1
        dρ[1,j] = -(1/t)*(  P(0,j-1)*ρ[1,j-1] + Q(0,j-1)*ρ[1,j] + R(0,j-1)*ρ[1,j+1]   )
    end
    # n = 0, l = L -1
    dρ[1,L]     = -(1/t)*(  P(0,L-1)*ρ[1,L-1] +  Q(0,L-1)*ρ[1,L]                    )

    for i in 2:3
            # Collision kernal zero for number and energy desnity, n = 1,2 (Array index 2,3) and l = 0 (Array index 1).
            # n = 1,2, l = 0
            dρ[i,1] = -(1/t)*(                      Q(i-1,0)*ρ[i,1] + R(i-1,0)*ρ[i,2]    )
        # Running through all l > 0 moments.
        for j in 2:L-1 
            # n = i-1, l = j-1                                                                     # Collision dependance.
            dρ[i,j] = -(1/t)*( P(i-1,j-1)*ρ[i,j-1] +  Q(i-1,j-1)*ρ[i,j] + R(i-1,j-1)*ρ[i,j+1]    )  - ωᵣ*ρ[i,j] 
        end
        # Truncation for l = L moment.
        # n = i-1, l = L-1
        dρ[i,L]     = -(1/t)*( P(i-1,L-1)*ρ[i,L-1] +  Q(i-1,L-1)*ρ[i,L]                      )  - ωᵣ*ρ[i,L]
    end

    return dρ
end



#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


nₐᵣ = [0,1,2]


nₙ  = 2
nₑ  = 3

T₀  = 1
Λ   = 1
m   = 0
μ₀  = 0


N = size(nₐᵣ)[1]
L = 10


tₛ = 0.1
tₑ = 50
ξ  = 0.01

#--------------------------------
η₀ = (tₛ/ξ)*((T₀)^(1+1))
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

println("m     : ", m)
println("T₀    : ",(1/3)*(ϵᵈ₀/nᵈ₀))
println("t₀    : ",tₛ)
println("η₀/s₀ : ",η₀)


p = (N,L,η₀,nₐᵣ,nₙ,nₑ)

tspan = trange((tₛ,tₑ),100,"log")


#------------------------------------------------

ρₜ = RK4(ρ₀,tspan,MDRTA_S,p)

#------------------------------------------------

ϵᵈ = ρₜ[:,nₑ,1]
nᵈ = ρₜ[:,nₙ,1]

T = (1/3)*(ϵᵈ./nᵈ)
τ = ((T.^(1+1)).*(tspan)./η₀)

plot(τ,T, xaxis=:log)
#plot(tspan,log.(nᵈ), xaxis=:log)
#plot(tspan,T, xaxis=:log)