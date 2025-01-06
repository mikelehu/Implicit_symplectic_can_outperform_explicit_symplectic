


function Ham_Schwarzschild(u,parms)
     
   E=parms[1]
   L=parms[2]
   β=parms[3]
   r=u[1]
   θ=u[2]
   pr=u[3]
   pθ=u[4]

   E2=E*E
   r2=r*r
   pr2=pr*pr
   pθ2=pθ*pθ
   sinθ=sin(θ)
   sinθ2=sinθ*sinθ
   invsinθ2=1/sinθ2
   invr=1/r
   invr2=invr*invr


   H1=1/2*invr2*invsinθ2*(L-β/2*r2*sinθ2)^2-1/2*(r/(r-2))*E2
   H2= 1/2*(pr2+pθ2*invr2)
   H3= -invr*pr2

   return H1+H2+H3

end

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
   #sinθ=x/r
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

function H1(u,parms)
     
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
   #sinθ=x/r
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

function H2(u,parms)
     
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
   #sinθ=x/r
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

function H3(u,parms)
     
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
   #sinθ=x/r
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



function flowH1Schwarzschild_!(uj,ej,h,parms)

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
   #sinθ=x/r
   #cosθ=y/r
   cosθ=x/r
   sinθ=y/r

   invr=1/r
   invr2=invr*invr
   invr3=invr2*invr
   sinθ2=sinθ*sinθ
   invsinθ=1/sinθ
   invsinθ2=invsinθ*invsinθ
   invsinθ3=invsinθ2*invsinθ

   gradpr=-1/4*β^2*r*sinθ2-E2/(r-2)^2+L^2*invsinθ2*invr3
   gradpθ=L^2*cosθ*invsinθ3*invr2-1/4*β^2*r2*sinθ*cosθ


   uj[3]=pr+h*gradpr
   uj[4]=pθ+h*gradpθ

   return nothing


end


function flowH1Schwarzschild!(uj,ej,h,parms)

# 2024-12-16

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
   #cosθ=x/r
   #sinθ=y/r

   #invr=1/r
   #invr2=invr*invr
   #invr3=invr2*invr
   #sinθ2=sinθ*sinθ
   #invsinθ=1/sinθ
   #invsinθ2=invsinθ*invsinθ
   #invsinθ3=invsinθ2*invsinθ

   #gradpr=-1/4*β^2*r*sinθ2-E2/(r-2)^2+L^2*invsinθ2*invr3
   #gradpθ=L^2*cosθ*invsinθ3*invr2-1/4*β^2*r2*sinθ*cosθ

   gradpr=-1/4*β^2*y^2/r-E2/(r-2)^2+L^2/(y^2*r)
   gradpθ=L^2*x/y^3-1/4*β^2*x*y

   uj[3]=pr+h*gradpr
   uj[4]=pθ+h*gradpθ

   return nothing


end


function flowH2Schwarzschild!(uj,ej,h,parms)

   x=uj[1]
   y=uj[2]
   pr=uj[3]
   pθ=uj[4]

   r2=x^2+y^2
   r=sqrt(r2)
   #sinθ=x/r
   #cosθ=y/r
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


function flowH3Schwarzschild!(uj,ej,h,parms)

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

