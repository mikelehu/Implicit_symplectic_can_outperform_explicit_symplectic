
#
#  File: utils.jl
#      - run_many_Splitting
#      - plots_Splitting

function run_many_Splitting(r,s,ddt0,tspan)
    
    nruns=100
    alg=Splitting_alg(r=r)
    
    cpus=similar(ddt0)
    nflowsH3=similar(ddt0)
    retcodes=[true for k in ddt0]
    MaxΔH=[0. for i in ddt0]

    sols=Array{ODESolution}(undef,length(ddt0))
    flows=[flowH1Schwarzschild!, flowH2Schwarzschild!, flowH3Schwarzschild!]
    prob=SplittingProblem(flows, u0, tspan, parms)
    H0=Ham_Schwarzschild(u0_B,parms_B)

    for i in 1:length(ddt0)

        print(",",ddt0[i])
        dt0=ddt0[i]  

        # save_everystep=true
        m0=1
        sols[i]=solver_Splitting(prob, alg, dt=dt0, msteps=m0)
        
        if sols[i].success
           nflowsH3[i]=sols[i].stats.naccept*s           
        else
           retcodes[i]=false
           nflowsH3[i]=Inf
        end

        m0 = max(1,div(Int64(ceil((tF-t0)/ddt0[i])),1000))
        ΔH0 = map(x->Ham_Schwarzschild(BigFloat.(x),parms_B), sols[i].u)./H0.-1 
        MaxΔH[i]=maximum(abs.(ΔH0[2:end]-ΔH0[1:end-1]))
        #MaxΔH[i]=maximum(abs.(ΔH0))

        #save_everystep=false
        solx=solver_Splitting(prob, alg, dt=dt0, save_everystep=false)

        if solx.success
           cpus[i]=0.
           for k in 1:nruns
               cpus[i]+=@elapsed solver_Splitting(prob, alg, dt=dt0, save_everystep=false)
           end
           cpus[i]=cpus[i]/nruns
        else
           cpus[i]=Inf
        end

    end
    
    return sols, retcodes, nflowsH3, cpus, MaxΔH

end
    

function plots_Splitting(title,r, ddt0,sols,nflowsH3,cpus,MaxΔH)
    
    H0=Ham_Schwarzschild(u0_B,parms_B)
    MaxΔH=[0. for i in ddt0]

    pl1=plot(ddt0,cpus, seriestype=:scatter,label="", 
             title="CPU-time", xlabel="dt", ylabel="CPU");
    pl2=plot(ddt0,nflowsH3, seriestype=:scatter, label="",
             title="nflowsH3", xlabel="dt", ylabel="nflowsH3");

    pl3=plot(title="Error in Ham r=$r",xlabel="t ", ylabel="log10(H/H0)", 
             yscale=:log10, label="")

    yrange=(1e-30,1e-1)
    for  i in 1:length(ddt0)
       ΔH0 = map(x->Ham_Schwarzschild(BigFloat.(x),parms_B), sols[i].u)./H0.-1
       MaxΔH[i]=maximum(abs.(ΔH0))
       pl3=plot!(sols[i].t[2:end],abs.(ΔH0[2:end])/ddt0[i]^r,  ylims=yrange, labels="")
    end

    fig=plot(pl1,pl2,pl3, layout=(1,3), size=(950,300), plot_title=title, plot_titlevspan=0.2)
    
    return fig
    
end

