struct tcoeffs{tType}
    s::Int64
    a::Array{tType,1}
    b::Array{tType,1}
    c::Array{tType,1}
end

struct tcache{UT,tT,fT1,fT2,fT3,pT}
    flowH1::fT1
    flowH2::fT2
    flowH3::fT3
    parms::pT
    du::UT
    step_number::Array{Int64,0}
    length_u::Int64
    tf::tT
end

struct ODESolution{uType,tType}
    t::Array{tType,1}
    u::Array{uType,1}
    stats::SciMLBase.DEStats
    success::Bool
 end


struct SplittingProblem{fT,uT,tspanT,pT}
    flows::Array{fT,1}
    u0::uT
    tspan::tspanT
    p::pT
end


struct Splitting_alg{T, U}
    r::T
    rkn::U
end


function Splitting_alg(; r::Int=10, rkn::Bool=false)
    Splitting_alg{Int64, Bool}(r, rkn)
end

function solver_Splitting(
                  prob::SplittingProblem{fT,uType,tspanType,pType},
                  alg::Splitting_alg{r,rkn},
                  args...;
                  dt = zero(eltype(tspanType)),
                  msteps=1,
                  save_everystep=true
                )where{fT,uType,tspanType,pType,r,rkn}
#

    tType = eltype(tspanType)

    stats = SciMLBase.DEStats(0)
    stats.nf=0
    stats.nf2=0
    stats.naccept=0

    u0 = prob.u0
    tspan = prob.tspan
    parms = prob.p

    t0=tspan[1]
    tf= tspan[2]
    step_number = Array{Int64,0}(undef)
    step_number[] = 0
    length_u = length(u0)

    dt = min(abs(dt), abs(tf - t0))
    signdt = sign(tspan[2] - tspan[1])
    dts = Array{tType}(undef, 1)
    dts = [dt, dt, signdt]

    step_fun::Function=empty
    
    checks=true

    du=zero(u0)

    if length(prob.flows)==3
        step_fun=Base_3flows_mstep!
        cache=tcache(prob.flows[1], prob.flows[2], prob.flows[3], parms, du, step_number, length_u, tf)
    elseif length(prob.flows)==2
        step_fun=Base_2flows_mstep!
        cache=tcache(prob.flows[1], prob.flows[2], empty, parms, du, step_number, length_u, tf)
    else  
        @warn("There is not defined Splitting method based on $(length(flows)) flows")   
        checks=false     
    end

    if checks

        if !(alg.rkn)
            s, a, b, c = Splitting_Coefficients(alg.r,dt)
        else
            s, a, b, c = Splitting_rkn_Coefficients(alg.r,dt)
        end

        coeffs=tcoeffs(s,a,b,c)

        uu = uType[]
        tt = tType[]

        push!(uu, copy(u0))
        push!(tt, t0)

        ttj = [t0, zero(t0)]
        uj = copy(u0)
        ej = zero(u0)

    
        cont = true
        error_warn=0
        step_retcode=true


        while cont

            #println("step=", step_number[], ",msteps=", msteps, ",dt=",dt, ",tj=", ttj[1])
            
            step_number[]+= msteps
            
            step_retcode=step_fun(ttj,uj,ej,dts,msteps,stats,coeffs,cache)

            if !step_retcode
                error_warn=1
                cont=false
                break
            end

            if save_everystep 
                push!(tt, ttj[1])
                push!(uu, copy(uj))
            end

            if ttj[1] == tf 
                cont=false
                break
            end
    

        end

        stats.naccept = step_number[]

        if error_warn!=0

            @warn("Error during the integration warn=$error_warn")
            sol=ODESolution(tt,uu,stats,false)
        
        else

            if !save_everystep
                push!(uu,copy(uj))
                push!(tt,ttj[1])
            end

            sol=ODESolution(tt,uu,stats,true)

        end
    
    else
        sol=ODESolution(tt,uu,stats,false)
    end


    return (sol)


end

