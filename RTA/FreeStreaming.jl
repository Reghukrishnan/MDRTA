using DifferentialEquations
using LinearAlgebra
using SpecialFunctions
using Plots

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


function FreeStream(dρ::Matrix,ρ,t::Float64)


    for (i,nv) in enumerate(nₐᵣ)
        n = i

        dρ[n,1] = - (0  + Q(nv,0)*ρ[n,1] + R(nv,0)*ρ[n,2]  )/t  

        for l = 1:L-1 
            if l == L-1
                dρ[n,l+1] = - ( P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] +         0         )/t
            else
                dρ[n,l+1] = - ( P(nv,l)*ρ[n,l]  + Q(nv,l)*ρ[n,l+1] + R(nv,l)*ρ[n,l+2]  )/t
                
            end
        end
    end
    return dρ
    
end




nₐᵣ = [0,0.5,1,1.5,2,2.5,3,3.5]


nₙ  = 3
nₑ  = 5

T₀  = 1
m   = 0
μ₀  = 0




N = size(nₐᵣ)[1]
L = 100



tₛ = 0.1
tₑ = 1000
ξ  = 0.01

η₀ = (tₛ/ξ)*T₀


ρ₀ = Matrix{Float64}(undef, N, L)

for (n,nv) in enumerate(nₐᵣ)    
    for l in 0:L-1        
        ρ₀[n,l+1] = ρeq(T₀,μ₀,nv,l)
    end
end

ϵᵈ₀ = ρ₀[nₑ,1]
nᵈ₀ = ρ₀[nₙ,1]

println("RTA  \n Relaxing lowest moment. n ϵ {0,0.5,1,1.5,2,2.5,3,3.5} \n")

println("m     : ", m)
println("T₀    : ",(1/3)*(ϵᵈ₀/nᵈ₀))
println("t₀    : ",tₛ)
println("η₀/s₀ : ",η₀)




tspan = trange((tₛ,tₑ),1000,"log")

#------------------------------------------------
println("\nSolving FreeStreaming System.")


ρₜ = RK4(ρ₀,tspan,FreeStream)

ϵᵈ = ρₜ[:,nₑ,1]
nᵈ = ρₜ[:,nₙ,1]

T = (1/3)*(ϵᵈ./nᵈ)
τ = (T.*(tspan)./η₀)

function pltmom(n,l)

    ρ_half0 = ρₜ[:,n+1,l+1]
    plot!(tspan,ρ_half0, xaxis=:log)    
end



println("Plotting Temperature proper time graph.")
n = 1
s = 1
e = 1000
#plot(tspan,T, xaxis=:log, xlabel="τ", ylabel="T",label="η = $(η₀), Λ = 0",dpi=300)
plot(tspan[s:e],ρₜ[s:e,n,3], xaxis=:log,label="ρ[$(nₐᵣ[n]),3]",xlabel="t", ylabel="ρ", title="Moment Evolution $(nₐᵣ[n]), max l = $(L)", fmt = :pdf,dpi = 300)

#plot!(tspan[s:e],[0 for i in tspan[s:e]],colour = "black",label="") 

#plot(tspan,ρₜ[:,3,2])

#plot!(tspan,ρₜ[100:,3,4])
plot!(tspan[s:e],ρₜ[s:e,n,5],label="ρ[$(nₐᵣ[n]),5]")
plot!(tspan[s:e],ρₜ[s:e,n,7],label="ρ[$(nₐᵣ[n]),7]") 
plot!(tspan[s:e],ρₜ[s:e,n,9],label="ρ[$(nₐᵣ[n]),9]")
plot!(tspan[s:e],ρₜ[s:e,n,17],label="ρ[$(nₐᵣ[n]),17]") 
plot!(tspan[s:e],ρₜ[s:e,n,37],label="ρ[$(nₐᵣ[n]),37]")
#plot!(tspan[s:e],ρₜ[s:e,n,57],label="ρ[$(nₐᵣ[n]),17]") 

savefig("Plots/The_Problem/Rn$(nₐᵣ[n])l_$(L)_E.pdf")