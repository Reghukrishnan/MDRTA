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
        return [tₛ + i*step for i in 1:N ]
    else
        println("Invalid Interpolation input.") 
    end    

    
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

lr = range(0,10,100)

n = 4
Pr = [P(n,i) for i in lr]
Qr = [Q(n,i) for i in lr]
Rr = [R(n,i) for i in lr]

plot(lr,Pr,xlabel="l", ylabel="Cf",label="P",dpi=300)
plot!(lr,Qr,xlabel="l", ylabel="Cf",label="Q",dpi=300)
plot!(lr,Rr,xlabel="l", ylabel="Cf",label="R",dpi=300)
plot!(lr,Pr.-Rr,xlabel="l", ylabel="Cf",label="P-R",dpi=300)