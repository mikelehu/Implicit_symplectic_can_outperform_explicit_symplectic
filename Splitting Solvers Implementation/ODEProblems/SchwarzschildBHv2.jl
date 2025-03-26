
#
#  Schwarzschild problem
#  (cartesian coordinates)
#  


function Ham_Schwarzschildv2(u,parms)
     
   E=parms[1]
   L=parms[2]
   β=parms[3]
   x=u[1]
   y=u[2]
   pr=u[3]
   pθ=u[4]

   r2=x^2+y^2
   E2=E*E
   r=sqrt(r2)
   pr2=pr*pr
   pθ2=pθ*pθ
   sinθ=y/r
   sinθ2=sinθ*sinθ
   invsinθ2=1/sinθ2
   invr=1/r
   invr2=invr*invr


   H1=1/2*invr2*invsinθ2*(L-β/2*r2*sinθ2)^2-1/2*(r/(r-2))*E2
   H2= 1/2*(pr2+pθ2*invr2)
   H3= -invr*pr2

   return H1+H2+H3

end

function H1v2(u,parms)
     
   E=parms[1]
   L=parms[2]
   β=parms[3]
   x=u[1]
   y=u[2]
   pr=u[3]
   pθ=u[4]

   r2=x^2+y^2
   E2=E*E
   r=sqrt(r2)
   pr2=pr*pr
   pθ2=pθ*pθ
   sinθ=y/r
   sinθ2=sinθ*sinθ
   invsinθ2=1/sinθ2
   invr=1/r
   invr2=invr*invr


   H1=1/2*invr2*invsinθ2*(L-β/2*r2*sinθ2)^2-1/2*(r/(r-2))*E2
   H2= 1/2*(pr2+pθ2*invr2)
   H3= -invr*pr2

   return H1

end

function H2v2(u,parms)
     
   E=parms[1]
   L=parms[2]
   β=parms[3]
   x=u[1]
   y=u[2]
   pr=u[3]
   pθ=u[4]

   r2=x^2+y^2
   E2=E*E
   r=sqrt(r2)
   pr2=pr*pr
   pθ2=pθ*pθ
   sinθ=y/r
   sinθ2=sinθ*sinθ
   invsinθ2=1/sinθ2
   invr=1/r
   invr2=invr*invr


   H1=1/2*invr2*invsinθ2*(L-β/2*r2*sinθ2)^2-1/2*(r/(r-2))*E2
   H2= 1/2*(pr2+pθ2*invr2)
   H3= -invr*pr2

   return H2

end

function H3v2(u,parms)
     
   E=parms[1]
   L=parms[2]
   β=parms[3]
   x=u[1]
   y=u[2]
   pr=u[3]
   pθ=u[4]

   r2=x^2+y^2
   E2=E*E
   r=sqrt(r2)
   pr2=pr*pr
   pθ2=pθ*pθ
   sinθ=y/r
   sinθ2=sinθ*sinθ
   invsinθ2=1/sinθ2
   invr=1/r
   invr2=invr*invr


   H1=1/2*invr2*invsinθ2*(L-β/2*r2*sinθ2)^2-1/2*(r/(r-2))*E2
   H2= 1/2*(pr2+pθ2*invr2)
   H3= -invr*pr2

   return H3

end


function flowH1Schwarzschildv2!(uj,ej,h,parms)


   E=parms[1]
   L=parms[2]
   β=parms[3]

   x=uj[1]
   y=uj[2]
   pr=uj[3]
   pθ=uj[4]

   r2=x^2+y^2
   E2=E*E
   r=sqrt(r2)

   gradpr=-1/4*β^2*y^2/r-E2/(r-2)^2+L^2/(y^2*r)
   gradpθ=L^2*x/y^3-1/4*β^2*x*y

   uj[3]=pr+h*gradpr
   uj[4]=pθ+h*gradpθ

   return nothing


end


function flowH2Schwarzschildv2!(uj,ej,h,parms)

   x=uj[1]
   y=uj[2]
   pr=uj[3]
   pθ=uj[4]

   r2=x^2+y^2
   r=sqrt(r2)
   cosθ=x/r
   sinθ=y/r
    
   aux = pθ/r
   x = aux*h
   y = r + pr*h
   r1 = sqrt(x^2+y^2)  
   pr = (aux*x + pr*y)/r1  
   uj[1]=(cosθ*x+sinθ*y)
   uj[2]=(sinθ*x-cosθ*y)
   uj[3]=pr

   return(nothing)

end


function flowH3Schwarzschildv2!(uj,ej,h,parms)

   x=uj[1]
   y=uj[2]
   pr=uj[3]
   pθ=uj[4]

   r2=x^2+y^2
   r=sqrt(r2)

   k=-pr^2/r
   lag=pr*r+3*k*h
   pr1=cbrt(-k*lag)
   r1=lag/pr1
   lag2=r1/r

   uj[1]=lag2*x
   uj[2]=lag2*y
   uj[3]=pr1

   return nothing

end

