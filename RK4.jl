function dExp(dy,y,t,γ)
    dy = γ*y
    return dy
end


"function RK4(y₀,ts,fun,p = nothing)
    
    ysize = (size(ts)...,size(y₀)...)  # Get the size of the time series amtrix
    y  = Array{Float64}(undef,ysize )   # Create the time series matrix
    dy = similar(y₀)                   # Create matrix to store dy
    view(y,1,:,:) .= y₀               # Initialise t = tᵢ value of y to y₀
    
    if isnothing(p)
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:,:)

            k₁ = fun(dy, yᵢ,t)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt)
            k₅ = fun(dy, yᵢ .+     dt*k₂, t +     dt)

            view(y,i+1,:,:) .= yᵢ .+ (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

        end

    else
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:,:)

            k₁ = fun(dy, yᵢ,t)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt,p)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt,p)
            k₅ = fun(dy, yᵢ .+     dt*k₂, t +     dt,p)

            view(y,i+1,:,:) = yᵢ .+ (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

        end
    end
    return y
    
end"


function RK4(y₀::Number,ts,fun,p = nothing)
    
    ysize = size(ts)  # Get the size of the time series amtrix
    y  = Array{Float64}(undef,ysize )   # Create the time series matrix
    
5gy[1] = y₀               # Initialise t = tᵢ value of y to y₀
    dy = 0
    if isnothing(p)
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            

            k₁ = fun(dy, y[i],t)
            
            k₂ = fun(dy, y[i] + 0.5*dt*k₁, t + 0.5*dt)
            k₃ = fun(dy, y[i] + 0.5*dt*k₂, t + 0.5*dt)
            k₅ = fun(dy, y[i] +     dt*k₂, t +     dt)

            y[i+1] .= y[i] + (dt*( k₁ + 2k₂ + 2k₃ + k₄ )/6)

        end

    else
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]            
            

            k₁ = fun(dy, y[i],t)
            k₂ = fun(dy, y[i] .+ 0.5*dt*k₁, t + 0.5*dt,p)
            k₃ = fun(dy, y[i] .+ 0.5*dt*k₂, t + 0.5*dt,p)
            k₄ = fun(dy, y[i] .+     dt*k₂, t +     dt,p)

            y[i+1] = y[i] .+ (dt*( k₁ .+ 2k₂ .+ 2k₃ .+ k₄ )/6)

        end
    end
    return y
    
end


ts = range(0.0,10.0,step = 0.01)
y = RK4(1.0,ts,dExp,0.10)

plot(ts,y)
