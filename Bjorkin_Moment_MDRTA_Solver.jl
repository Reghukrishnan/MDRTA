using DifferentialEquations
using LinearAlgebra
using SpecialFunctions
using Plots

Γ = gamma

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





function RTA(dρ,ρ,p,t)

    print("Time : " , t,"\n")

    N   = p[1]     # Number of n values
    L   = p[2]     # Number of l values
    η   = p[3]     # Initial viscosity
    nₐᵣ = p[4]     # Array with n values
    nₙ  = p[5]     # position of number density moment
    nₑ  = p[6]     # position of energy density moment
    Λ   = p[7]     # The exponent of momentum dependance 

    nₘ = nₐᵣ[1]    # minimum n value

    Nᵈ = ρ[nₙ,1]             # Number density
    Eᵈ = ρ[nₑ,1]            # Energy density

    print("E : ", Eᵈ,"\n")
    print("N : ", Nᵈ,"\n")

    T = (1/3)*(Eᵈ/Nᵈ)       # Temperature
    

    μ = T*log( 27*(π^2)*Nᵈ*((Nᵈ/Eᵈ)^3)) # Chemical potential

    #print(" Comp : ",ρeq(T,μ,2,0),"\n")
    print(" Comp : ",ρeq(T,μ,1,0),"\n")
    

    a₁ = ρ[nₙ-1,1]/ρeq(T,μ,1-Λ,0) 
    a₂ = ρeq(T,μ,1-Λ,0)/ρeq(T,μ,2-Λ,0)
    a₃ = (- ρeq(T,μ,1-Λ,0)*ρeq(T,μ,2-Λ,0)^2 + (ρeq(T,μ,1-Λ,0)^2)*ρeq(T,μ,3-Λ,0))/(ρeq(T,μ,2-Λ,0)^2)
    B  = (ρ[nₙ-1,1]*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρ[nₑ-1,1])/(a₃*ρeq(T,μ,2-Λ,0))

    ωᵣ = (T^(1+Λ))/η₀

    # Truncating coupling to lowest n value.
    ρ[1,1] = - (0  
                + Q(nₘ,0)*ρ[1,1] 
                    + R(nₘ,0)*ρ[1,2] 
                        - B*ωᵣ*t*( a₂*ρeq(T,μ,nₘ +1 -Λ,0)) )/t  
    for l =1:L-2
        ρ[1,l+1] = - (P(nₘ,l)*ρ[1,l]  
                    + Q(nₘ,l)*ρ[1,l+1] 
                     + R(nₘ,l)*ρ[1,l+2] )/t
    end
    ρ[1,L] = - (P(nₘ,L-1)*ρ[1,L-1]  
                    + Q(nₘ,L-1)*ρ[1,L] 
                    + 0 )/t


    for (i,nv) in enumerate(nₐᵣ[begin+1:end])
        n = i+1
        "Cnl = ωᵣ*t*(    (ρ[n-1,1]* ρeq(T,μ,1-Λ,0)   - ρ[nₙ-1,1]*ρeq(T,μ,nv -Λ,0) )*ρeq(T,μ,2-Λ,0)
                    -B*ρeq(T,μ,1-Λ,0)*( ρeq(T,μ,nv -Λ,0)*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρeq(T,μ,nv +1-Λ,0) )   )/(ρeq(T,μ,2-Λ,0)*ρeq(T,μ,1-Λ,0))" # The Collision kernel
        Cnl = ωᵣ*t*(   (ρ[n-1,1]* ρeq(T,μ,1-Λ,0)   - ρ[nₙ-1,1]*ρeq(T,μ,nv -Λ,0) )*ρeq(T,μ,2-Λ,0)* (a₃*ρeq(T,μ,2-Λ,0))
                    -(ρ[nₙ-1,1]*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρ[nₑ-1,1])*ρeq(T,μ,1-Λ,0)*( ρeq(T,μ,nv -Λ,0)*ρeq(T,μ,2-Λ,0) - ρeq(T,μ,1-Λ,0)*ρeq(T,μ,nv +1-Λ,0) )   )/(ρeq(T,μ,2-Λ,0)*ρeq(T,μ,1-Λ,0)*(a₃*ρeq(T,μ,2-Λ,0)))
        if nv==1 
            print("Col_N : ",Cnl," ","\n")
        elseif nv ==2
            #print("Col_E : ",Cnl," ","\n")
        end

        #dρ[n,l+1] = - (0  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2] + Cnl )/t

        dρ[n,1] = - (0  + Q(nv,0)*ρ[n,1] + R(nv,0)*ρ[n,2] + Cnl )/t  # For  l = 0 we have the complicated collisio term
        
        if nv ==1 || nv ==2
            print("dρ : ",dρ[n,1],"\n")
        end
        
        for l = 1:L-1                
                
            if l == L-1
                dρ[n,l+1] = - ( P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] +         0        + ωᵣ*t*(ρ[n-1,l+1]) )/t
            else
                dρ[n,l+1] = - ( P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2] + ωᵣ*t*(ρ[n-1,l+1]) )/t
                
            end
        end
    end
    return dρ
end


function RT4(y₀,ts,fun,p)

    tᵢ = size(ts)[1]
    n,l = size(y₀)

    y  = Array{Float64}(undef,tᵢ,n,l )
    dy = similar(y₀)

    for (i,t) in enumerate(ts[begin:end-1])
        dt = ts[i+1] - ts[i]
        
        k₁ = fun(dy, y[i],t)
        k₂ = fun(dy, y[i] + 0.5*dt*k₁, p, t + 0.5*dt)
        k₃ = fun(dy, y[i] + 0.5*dt*k₂, p, t + 0.5*dt)
        k₅ = fun(dy, y[i] +     dt*k₂, p, t +     dt)

        y[i+1] = y[i] + (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

    end
    return y
    
end



Λ = 0.5
nₐᵣ = [0,0.5,1,1.5,2,2.5,3,3.5]


nₙ  = 3
nₑ  = 5

T₀  = 1
m   = 0
μ₀  = 0

N = size(nₐᵣ)[1]
L = 10


ti = 0.1
tf = 20
ξ  = 0.1

η = (ti*T₀)/ξ

tspan = (ti,tf)

ρ₀ = Matrix{Float64}(undef, N, L)

for (n,nv) in enumerate(nₐᵣ)    
    for l in 0:L-1        
        ρ₀[n,l+1] = ρeq(T₀,μ₀,nv,l)
    end
end



p = [N,L,η,nₐᵣ,nₙ,nₑ,Λ]

#prob = ODEProblem(RTA,ρ₀,tspan,p)
#sol = solve(prob,RK4(),dt =0.1)

#T = (1/3)*(sol[nₑ,1,:]./sol[nₙ,1,:])
#τ = (T.*(sol.t)./η)
print("Zero")
#plot(τ,T, xaxis=:log)
#plot(sol.t,sol[nₙ,1,:], xaxis=:log)
