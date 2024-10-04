using DifferentialEquations
using LinearAlgebra
using SpecialFunctions
using Plots
using JLD

Γ = gamma


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
        return [tₛ + i*step for i in 1:N ]
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








function RTA(dρ::Matrix,ρ,t::Float64,p)
    N   = p[1]
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
        return (Γ(n+2)*((T^(n+2))/(2*(π^2)))*exp(μ/T))/(2*l+1)        
    else
        return 0 
    end
end



function MDRTA(dρ::Matrix,ρ,t::Float64,p)

    #print("Time : " , t,"\n")

    N   = p[1]     # Number of n values
    L   = p[2]     # Number of l values
    η   = p[3]     # Initial viscosity
    nₐᵣ = p[4]     # Array with n values
    nₙ  = p[5]     # position of number density moment
    nₑ  = p[6]     # position of energy density moment
    Λ   = p[7]     # The exponent of momentum dependance 

    if Λ != 0
        s = 1
    else
        s = 0
    end

    nₘ = nₐᵣ[1]    # minimum n value

    nᵈ = ρ[nₙ,1]
    ϵᵈ = ρ[nₑ,1]

    # Computing Temperature and Chemical potential
    T = (1/3)*(ϵᵈ/nᵈ)    
    μ = T*log( 27*(π^2)*nᵈ*((nᵈ/ϵᵈ)^3))

    #print(" Comp : ",ρeq(T,μ,2,0),"\n")
    #print(" Comp : ",ρeq(T,μ,1,0),"\n")
    

    a₁ = ρ[nₙ-s,1]/ρeq(T,μ,1-Λ,0) 
    a₂ = ρeq(T,μ,1-Λ,0)/ρeq(T,μ,2-Λ,0)
    a₃ = (- ρeq(T,μ,1-Λ,0)*ρeq(T,μ,2-Λ,0)^2 + (ρeq(T,μ,1-Λ,0)^2)*ρeq(T,μ,3-Λ,0))/(ρeq(T,μ,2-Λ,0)^2)
    B  = (ρ[nₙ-1,1]*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρ[nₑ-1,1])/(a₃*ρeq(T,μ,2-Λ,0))

    ωᵣ = (T^(1+Λ))/η₀

    # Truncating coupling to lowest n value.
    #println("ρ[0,0] :", ρ[1,1])
    dρ[1,1] = - (0  
                + Q(nₘ,0)*ρ[1,1] 
                    + R(nₘ,0)*ρ[1,2] 
                         )/t  #- B*ωᵣ*( a₂*ρeq(T,μ,nₘ +1 -Λ,0))
    for l =1:L-2
        dρ[1,l+1] = - (P(nₘ,l)*ρ[1,l]  
                    + Q(nₘ,l)*ρ[1,l+1] 
                     + R(nₘ,l)*ρ[1,l+2] )/t
    end
    dρ[1,L] = - (P(nₘ,L-1)*ρ[1,L-1]  
                    + Q(nₘ,L-1)*ρ[1,L] 
                    + 0 )/t


    for (i,nv) in enumerate(nₐᵣ[begin+1:end])
        n = i+1
        
        Cnl = ωᵣ*t*(   (ρ[n-s,1]* ρeq(T,μ,1-Λ,0)   - ρ[nₙ-s,1]*ρeq(T,μ,nv -Λ,0) )*ρeq(T,μ,2-Λ,0)* (a₃*ρeq(T,μ,2-Λ,0))
                    -(ρ[nₙ-s,1]*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρ[nₑ-s,1])*ρeq(T,μ,1-Λ,0)*( ρeq(T,μ,nv -Λ,0)*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρeq(T,μ,nv +1-Λ,0) )   )/(ρeq(T,μ,2-Λ,0)*ρeq(T,μ,1-Λ,0)*(a₃*ρeq(T,μ,2-Λ,0)))
        """if nv==1 
            print("Col_N : ",Cnl," ","\n")
        elseif nv ==2
            print("Col_E : ",Cnl," ","\n")
        end"""

        #dρ[n,l+1] = - (0  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2] + Cnl )/t

        dρ[n,1] = - (0  + Q(nv,0)*ρ[n,1] + R(nv,0)*ρ[n,2] + Cnl )/t  # For  l = 0 we have the complicated collisio term
        
        """if nv ==1 || nv ==2
            print("dρ : ",dρ[n,1],"\n")
        end"""
        
        for l = 1:L-1                
                
            if l == L-1
                dρ[n,l+1] = - ( P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] +         0        + ωᵣ*t*(ρ[n-s,l+1]) )/t
            else
                dρ[n,l+1] = - ( P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2] + ωᵣ*t*(ρ[n-s,l+1]) )/t
                
            end
        end
    end
    return dρ
end







Λ = 0
nₘₐₓ = 3

# Number of n values to evaluate
if Λ != 0
     N= nₘₐₓ/Λ
else
    N= nₘₐₓ
end

# Constructing the n values
if Λ != 0
    nₐᵣ = [Λ*i for i in 0:N]
else
    nₐᵣ = [i for i in 0:N]
end

        

nₙ = findall(x->x ==1,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==2,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density


println(N)
println(nₐᵣ)



T₀  = 1
m   = 0
μ₀  = 0


N = size(nₐᵣ)[1]
L = 150


tₛ = 0.1
tₑ = 100
ξ  = 0.01


#------------------------------------------

η₀ = (tₛ/ξ)*((T₀)^(1+Λ))

#------------------------------------------



ρ₀ = Matrix{Float64}(undef, N, L)

for (n,nv) in enumerate(nₐᵣ)    
    for l in 0:L-1        
        ρ₀[n,l+1] = ρeq(T₀,μ₀,nv,l)
    end
end


ϵᵈ₀ = ρ₀[nₑ,1]
nᵈ₀ = ρ₀[nₙ,1]

println("RTA  \n Relaxing lowest moment. n ϵ {0,0.5,1,1.5,2,2.5,3,3.5} \n")


η₀ = (tₛ/ξ)*T₀

println("m     : ", m)
println("T₀    : ",(1/3)*(ϵᵈ₀/nᵈ₀))
println("t₀    : ",tₛ)
println("η₀/s₀ : ",η₀)


p = (N,L,η₀,nₐᵣ,nₙ,nₑ,Λ)

Nₚ = 1000
tspan = trange((tₛ,tₑ),Nₚ,"log")

#------------------------------------------------
println("\nSolving Denicol RTA.")


ρₜ = RK4(ρ₀,tspan,MDRTA,p)

ϵᵈ = ρₜ[:,nₑ,1]
nᵈ = ρₜ[:,nₙ,1]

T = (1/3)*(ϵᵈ./nᵈ)
τ = ((T.^(1+Λ)).*(tspan)./η₀)


y = (tspan,ρₜ,nₐᵣ,p)
println("Saving data.")
save("Data/MDRTA_FCK_T_$(T₀)_Lam_$(Λ)_n_0_3_L_$(L).jld","y",y)


s = Int(1)
e = Nₚ
println("Plotting Temperature proper time graph.")
plot(tspan[s:e],T[s:e], xaxis=:log ,xlabel="τ", ylabel="T",label="η = $(η₀), Λ = $(Λ)",dpi=300)




