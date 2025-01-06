function Potential_U(q1,q2)


    return 1/2*(q1*q1+q2*q2)+q1*q1*q2-1/3*q2*q2*q2

end

function HenonHeliesHam(u,parms)

    q1=u[1]
    q2=u[2]
    p1=u[3]
    p2=u[4]
    
    U=Potential_U(q1,q2)
    
    return 1/2*(p1*p1+p2*p2)+U

end


function H1(u)

    p1=u[3]
    p2=u[4]
      
    return 1/2*(p1*p1+p2*p2)

end

function H2(u)

    q1=u[1]
    q2=u[2]

    U=Potential_U(q1,q2)

    return U

end


function HenonHeilesODE!(du,u,parms,t)

  q1=u[1]
  q2=u[2]
  p1=u[3]
  p2=u[4]

  du[1]=p1
  du[2]=p2
  
  du[3]=-q1-2*q1*q2
  du[4]=-q2-q1*q1+q2*q2
  
  return nothing
  
end

function flowH1HenonHeiles!(uj,ej,h,parms)

   q1=uj[1]
   q2=uj[2]
   p1=uj[3]
   p2=uj[4]

   uj[1]=q1+h*p1
   uj[2]=q2+h*p2

   return nothing

end


function flowH2HenonHeiles!(uj,ej,h,parms)

    q1=uj[1]
    q2=uj[2]
    p1=uj[3]
    p2=uj[4]
 
    uj[3]=p1+h*(-q1-2*q1*q2)
    uj[4]=p2+h*(-q2-q1*q1+q2*q2)

    return nothing
 
 end


