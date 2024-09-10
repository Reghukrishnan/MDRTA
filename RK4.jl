using LinearAlgebra
using Plots





function dExp!(dy,y,t,γ)
    dy = γ*y

    return dy
end



#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
"""
    Runge–Kutta implimentation for first order differential equations of one variable.
    Returns an array with [y(t)]
"""
function RK4!(y₀::Number,ts,fun,p = nothing)
    print("Number RK4 Solver Called. \n")

    ysize = size(ts)  # Get the size of the time series amtrix
    y  = Array{Float64}(undef,ysize )   # Create the time series matrix
    
    y[1] = y₀               # Initialise t = tᵢ value of y to y₀
    dy = 0.0
    if isnothing(p)
        for (i,t) in enumerate(ts[begin:end-1])

            dt = (ts[i+1] - ts[i])[1]
            
            

            k₁ = fun( dy, y[i]            , t          )            
            k₂ = fun( dy, y[i] + 0.5*dt*k₁, t + 0.5*dt )
            k₃ = fun( dy, y[i] + 0.5*dt*k₂, t + 0.5*dt )
            k₄ = fun( dy, y[i] +     dt*k₃, t +     dt )

            y[i+1] = y[i] + (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

        end

    else
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]            
            
            
            k₁ = fun(dy, y[i], t, p)  

            k₂ = fun(dy, y[i] + 0.5*dt*k₁, t + 0.5*dt, p)            
            k₃ = fun(dy, y[i] + 0.5*dt*k₂, t + 0.5*dt, p)
            k₄ = fun(dy, y[i] +     dt*k₃, t +     dt, p)
            
            y[i+1] = y[i] + (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

        end
    end
    return y
    
end

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
"""
    Runge–Kutta implimentation for first order coupled differential equations organised in n x 1 matrix form.
    Returns an array with [y(t)] 
"""
function dLinear(dy::Vector,y,t::Number,A::Matrix)
    dy = A*y 
    return dy

end

function RK4!(y₀::Vector,ts,fun,p = nothing)
    print("Vector RK4 Solver Called. \n")

    ysize = (size(ts)...,size(y₀)...)  # Get the size of the time series amtrix
    y  = zeros(ysize )   # Create the time series matrix

    dy = similar(y₀)                  # Create matrix to store dy
    view(y,1,:) .= y₀               # Initialise t = tᵢ value of y to y₀


    
    if isnothing(p)
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:)

            k₁ = fun(dy, yᵢ,t)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt)
            k₄ = fun(dy, yᵢ .+     dt*k₃, t +     dt)

            view(y,i+1,:) .= yᵢ .+ (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

        end

    else
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:)

            k₁ = fun(dy, yᵢ,t,p)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt,p)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt,p)
            k₄ = fun(dy, yᵢ .+     dt*k₃, t +     dt,p)

            view(y,i + 1,:) .= yᵢ .+ (dt*( k₁ .+ 2k₂ .+ 2k₃ .+ k₄ )/6)

        end
    end
    return y
    
end


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

"""
    Runge–Kutta implimentation for first order coupled differential equations organised in n x l matrix form.
    Returns an array with [y(t)] 
"""
function RK4!(y₀::Matrix,ts::Vector,fun,p = nothing)
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






# Eigen Values = -1,3
# Eigen Vectors = (1/√2)(±1,1)


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

ts = range(0.0,10,step = 0.01)

"""A = 1
γ = 0.1
y = RK4( A,ts,dExp,γ)

plot(ts,[y,A*exp.(γ*ts)],label=["comp" "exact"])"""


"""
    Solution to the harmonic oscillator Equation.

    The second order Equation is rewritten as two coupled
    first order differential equations.
"""
function ExampleVector()

    ω = 2*π/5
    A = [ 0.0 1.0
        -ω^2 0.0] 

    y₀ = [0
        2]
    y = RK4(y₀,ts,dLinear,A)

    U = (1/sqrt(2))*[1 -1
        1 1]

    Uᵀ = (1/sqrt(2))*[1 1
        -1 1]

    plot(ts,y[:,1])
end

