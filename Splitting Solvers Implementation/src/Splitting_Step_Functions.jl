

function Base_3flows_mstep!(ttj::Array{tType,1},
                            uj::uType,
                            ej::uType,
                            dts::Array{tType,1},
                            msteps::Int64,
                            stats::SciMLBase.DEStats,
                            coeffs::tcoeffs{tType},
                            cache::tcache{uType, tType, fT1,fT2,fT3, pT}) where  {uType, tType,fT1,fT2,fT3, pT}

    flowH1=cache.flowH1
    flowH2=cache.flowH2
    flowH3=cache.flowH3
    parms=cache.parms
    step_number=cache.step_number
    len = cache.length_u
    lenq=div(len,2)
    tf=cache.tf

    s=coeffs.s
    a=coeffs.a
    b=coeffs.b
    c=coeffs.c

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = signdt*dt

    step_retcode=true

    dt = dts[1]
    sdt = signdt*dt


    for i in 1:msteps

        for is in 1:s

            flowH1(uj,ej,sdt*a[is],parms)
            flowH2(uj,ej,sdt*b[is],parms)
            flowH3(uj,ej,sdt*c[is],parms)
            flowH2(uj,ej,sdt*b[is],parms)
  
        end

        flowH1(uj,ej,sdt*a[s+1],parms)


        dtmax = abs((tf-ttj[1])-ttj[2])
        if abs(sdt) >= dtmax
            ttj[1] = tf
            ttj[2] = 0
        else
            res = Base.TwicePrecision(ttj[1], ttj[2]) + sdt
            ttj[1] = res.hi
            ttj[2] = res.lo
        end
    
        dts[1]=min(abs(sdt),dtmax)
        dts[2]=dt

    end

    return(step_retcode)

end




function Base_2flows_mstep!(ttj::Array{tType,1},
                            uj::uType,
                            ej::uType,
                            dts::Array{tType,1},
                            msteps::Int64,
                            stats::SciMLBase.DEStats,
                            coeffs::tcoeffs{tType},
                            cache::tcache{uType, tType, fT1,fT2,fT3, pT}) where  {uType, tType,fT1,fT2,fT3, pT}

    flowH1=cache.flowH1
    flowH2=cache.flowH2
    flowH3=cache.flowH3
    parms=cache.parms
    step_number=cache.step_number
    len = cache.length_u
    lenq=div(len,2)
    tf=cache.tf

    s=coeffs.s
    a=coeffs.a
    b=coeffs.b
    c=coeffs.c

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = signdt*dt

    step_retcode=true

    dt = dts[1]
    sdt = signdt*dt


    for i in 1:msteps

        for is in 1:s

            flowH1(uj,ej,sdt*a[is],parms)
            flowH2(uj,ej,sdt*c[is],parms)
  
        end

        flowH1(uj,ej,sdt*a[s+1],parms)


        dtmax = abs((tf-ttj[1])-ttj[2])
        if abs(sdt) >= dtmax
            ttj[1] = tf
            ttj[2] = 0
        else
            res = Base.TwicePrecision(ttj[1], ttj[2]) + sdt
            ttj[1] = res.hi
            ttj[2] = res.lo
        end
    
        dts[1]=min(abs(sdt),dtmax)
        dts[2]=dt

    end

    return(step_retcode)

end