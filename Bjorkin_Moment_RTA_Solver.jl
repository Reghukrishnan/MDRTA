using DifferentialEquations
using LinearAlgebra
using SpecialFunctions
using Plots


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





function RTA(dρ,ρ,p,t)
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




nₐᵣ = [1,2,3]


nₙ  = 1
nₑ  = 2

T₀  = 1
m   = 0
μ₀  = 0

N = size(nₐᵣ)[1]
L = 10


ti = 0.1
tf = 500
ξ  = 0.01

η₀ = (ti*T₀)/ξ

tspan = (ti,tf)

ρ₀ = Matrix{Float64}(undef, N, L)

for (n,nv) in enumerate(nₐᵣ)    
    for l in 0:L-1        
        ρ₀[n,l+1] = ρeq(T₀,μ₀,nv,l)
    end
end



p = [N,L,η₀,nₐᵣ,nₙ,nₑ]

prob = ODEProblem(RTA,ρ₀,tspan,p)
sol = solve(prob)

T = (1/3)*(sol[nₑ,1,:]./sol[nₙ,1,:])
τ = (T.*(sol.t)./η)

plot(τ,T, xaxis=:log)
#plot(sol.t,sol[nₙ,1,:], xaxis=:log)
