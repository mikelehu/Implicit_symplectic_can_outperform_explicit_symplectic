function InitialSchwarzschildv2(T = Float64)
    
    parms = BigFloat[0.995,4.6,8.9e-5]

    r = 11
    θ =  pi/2
    pr = 0

    u0 = BigFloat[r, θ, pr, 0]

    pθ = r * sqrt(-1 -2*Ham_Schwarzschild(u0,parms))
    parms = convert.(T, parms)
 
    x=r*cos(θ)
    y=r*sin(θ)
    #u0 = convert.(T, [r, θ, pr, pθ])
    u0 = convert.(T, [x, y, pr, pθ])

    return u0, parms

end

