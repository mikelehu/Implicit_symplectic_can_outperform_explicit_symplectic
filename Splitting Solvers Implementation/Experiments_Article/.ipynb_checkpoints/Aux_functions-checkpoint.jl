
struct tsolutions_IRKGL16
    dts
    sols
    retcodes
    iters
    cpus
    MaxΔHglobal
    MaxΔHlocal
end

struct tsolutions_Splitting
   s
   alg
   dts
   sols
   retcodes
   nflowsH3
   cpus
   MaxΔHglobal
   MaxΔHlocal
end

function run_many_IRKGL16(alg, prob, ddt0, HAM; nruns=1)


    u0=prob.u0
    parms=prob.p
    
    u0_B=BigFloat.(u0)
    parms_B=BigFloat.(parms)
    
    cpus=similar(ddt0)
    iters=similar(ddt0)
    retcodes=[true for k in ddt0]
    MaxΔHglobal=[0. for i in ddt0]
    MaxΔHlocal=[0. for i in ddt0]
    
    sols=Array{IRKGaussLegendre.ODESolution}(undef,length(ddt0))
    H0=HAM(u0_B,parms_B)
    
    for i in 1:length(ddt0)
        
        print(",",ddt0[i])    
        dt0=ddt0[i]
        
        sols[i]=solve(prob,alg,dt=dt0,adaptive=false)
        if sols[i].retcode==ReturnCode.Success
           iters[i]=sols[i].stats.nfpiter
        else
           retcodes[i]=false
           iters[i]=Inf
        end

        H = [HAM(BigFloat.(u),parms_B) for u in sols[i].u]
        ΔH0 = @. Float64(abs(H/H0-1))
        H_lerr = @. Float64(abs((H[2:end] / H[1:end-1]) - 1))
        MaxΔHglobal[i]=maximum(ΔH0)
        MaxΔHlocal[i]=maximum(H_lerr)
 
        # save_everystep=false
        solx=solve(prob,alg,dt=dt0,adaptive=false,save_everystep = false)
        if solx.retcode==ReturnCode.Success
           cpus[i]=0.
           for k in 1:nruns
               cpus[i]+=@elapsed solve(prob,alg,dt=dt0,adaptive=false, save_everystep = false)
           end
           cpus[i]=cpus[i]/nruns
        else
           cpus[i]=Inf
        end
        
    end

    return tsolutions_IRKGL16(ddt0, sols,retcodes,iters,cpus,MaxΔHglobal,MaxΔHlocal)

end    


function plots_IRKGL16(title, prob, HAM, solutions)

    dts=solutions.dts
    sols=solutions.sols
    iters=solutions.iters
    cpus=solutions.cpus
    MaxΔH=solutions.MaxΔHlocal

    u0=prob.u0
    parms=prob.p
    
    u0_B=BigFloat.(u0)
    parms_B=BigFloat.(parms)
    
    H0=HAM(u0_B,parms_B)
    MaxΔH=[0. for i in dts]

    pl1=plot(dts,cpus, seriestype=:scatter,label="", 
             title="CPU-time", xlabel="dt", ylabel="CPU");
    pl2=plot(dts,iters, seriestype=:scatter, label="",
         title="Iterations", xlabel="dt", ylabel="iter");

    pl3=plot(title="Error in Ham",xlabel="t ", ylabel="log10(H/H0)", 
              yscale=:log10, label="")

     for  i in 1:length(dts)
          m0 = max(1,div(Int64(ceil((tF-t0)/dts[i])),1000))
          ΔH0 = map(x->HAM(BigFloat.(x),parms_B), sols[i].u)./H0.-1 
          pl3=plot!(sols[i].t[2:m0:end],abs.(ΔH0[2:m0:end]), labels="")
     end

     fig=plot(pl1,pl2,pl3, layout=(1,3), size=(950,300),plot_title=title, plot_titlevspan=0.2)
     return fig
    
end


function run_many_Splitting(s, alg, prob, ddt0, HAM; nruns=1)
    

    u0=prob.u0
    parms=prob.p
    
    u0_B=BigFloat.(u0)
    parms_B=BigFloat.(parms)
    
    cpus=similar(ddt0)
    nflowsH3=similar(ddt0)
    retcodes=[true for k in ddt0]
    MaxΔHglobal=[0. for i in ddt0]
    MaxΔHlocal=[0. for i in ddt0]

    sols=Array{SplittingMethods.ODESolution}(undef,length(ddt0))
    H0=HAM(u0_B, parms_B)

    for i in 1:length(ddt0)

        print(",",ddt0[i])
        dt0=ddt0[i]  

        # save_everystep=true
        m0=1
        sols[i]=solver_Splitting(prob, alg, dt=dt0)

        if sols[i].success
           nflowsH3[i]=sols[i].stats.naccept*s           
        else
           retcodes[i]=false
           nflowsH3[i]=Inf
        end

        H = [HAM(BigFloat.(u),parms_B) for u in sols[i].u]
        ΔH0 = @. Float64(abs(H/H0-1))
        H_lerr = @. Float64(abs((H[2:end] / H[1:end-1]) - 1))
        MaxΔHglobal[i]=maximum(ΔH0)
        MaxΔHlocal[i]=maximum(H_lerr)

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
    
    return tsolutions_Splitting(s, alg, ddt0, sols,retcodes,nflowsH3,cpus,MaxΔHglobal,MaxΔHlocal)
 
end


function plots_Splitting(title, alg, prob, HAM, solutions)
    
    dts=solutions.dts
    sols=solutions.sols
    nflowsH3=solutions.nflowsH3
    cpus=solutions.cpus
    MaxΔH=solutions.MaxΔHlocal
    
    u0=prob.u0
    parms=prob.p

    r=alg.r

    u0_B=BigFloat.(u0)
    parms_B=BigFloat.(parms)

    H0=HAM(u0_B, parms_B)
    MaxΔH=[0. for i in dts]

    pl1=plot(dts,cpus, seriestype=:scatter,label="", 
             title="CPU-time", xlabel="dt", ylabel="CPU");
    pl2=plot(dts,nflowsH3, seriestype=:scatter, label="",
             title="nflowsH3", xlabel="dt", ylabel="nflowsH3");

    pl3=plot(title="Error in Ham r=$r",xlabel="t ", ylabel="log10(H/H0)", 
             yscale=:log10, label="")

    yrange=(1e-8,1)
    for  i in 1:length(dts)
       m0 = max(1,div(Int64(ceil((tF-t0)/dts[i])),1000))
       ΔH0 = map(x->HAM(BigFloat.(x),parms_B), sols[i].u)./H0.-1
       pl3=plot!(sols[i].t[2:m0:end],abs.(ΔH0[1:m0:end]),  labels="")  
    end

    fig=plot(pl1,pl2,pl3, layout=(1,3), size=(950,300), plot_title=title, plot_titlevspan=0.2)
    return fig
    
end


function find_index_ge_cpu(lookup_value, v)

   idx=0
   for k in length(v):-1:1

       if v[k]>=lookup_value
          idx=k
          break
       end

   end

   return idx

end

function compute_index_ge_cpu!(indexes, lookup_value, manysols)

   for i  in 1:length(manysols)

      idx=find_index_ge_cpu(lookup_value, manysols[i].cpus)
      push!(indexes,idx)

   end

   return nothing

end


function lcm_floats_vector(numbers::Vector{Float64})
   # Convert the numbers to rationals
   rationals = rationalize.(numbers)
   # Compute the LCM iteratively
   lcm_rationals = reduce(lcm, rationals)
   return lcm_rationals
end