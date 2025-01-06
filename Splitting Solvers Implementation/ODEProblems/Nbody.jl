function NbodyEnergy(u,Gm)
     N = length(Gm)
     zerouel = zero(eltype(u))
     T = zerouel
     U = zerouel
     for i in 1:N
        qi = u[:,i,1]
        vi = u[:,i,2]
        Gmi = Gm[i]
        T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
        for j in (i+1):N
           qj = u[:,j,1]
           Gmj = Gm[j]
           qij = qi - qj
           U -= Gmi*Gmj/norm(qij)
        end
     end
    1/2*T + U
end


function H1(u,Gm)

   N = length(Gm)
   zerouel = zero(eltype(u))
   T = zerouel
   for i in 1:N
      vi = u[:,i,2]
      Gmi = Gm[i]
      T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
   end

   return 1/2*T 

end

function H2(u,Gm)

   N = length(Gm)
   zerouel = zero(eltype(u))
   U = zerouel
   for i in 1:N
      qi = u[:,i,1]
      vi = u[:,i,2]
      Gmi = Gm[i]
      for j in (i+1):N
         qj = u[:,j,1]
         Gmj = Gm[j]
         qij = qi - qj
         U -= Gmi*Gmj/norm(qij)
      end
   end

   return  U

end


function NbodyBarycenter(u,Gm)
     N = length(Gm)
     dim = size(u,1)
     A = zeros(dim)
     B = zeros(dim)
     for i in 1:N
        qi = u[:,i,1]
        vi = u[:,i,2]
        Gmi = Gm[i]
        A += Gmi*qi
        B += Gmi*vi
     end
     return A, B
end




function flowH1Nbody!(uj,ej,h,Gm)

     N = length(Gm)

     for i in 1:3, j in 1:N
        uj[i,j,1] = uj[i,j,1] + h*uj[i,j,2]
     end

    return nothing
end



function flowH2Nbody!(uj,ej,h,Gm)

     N = length(Gm)
     for i in 1:N
        xi = uj[1,i,1]
        yi = uj[2,i,1]
        zi = uj[3,i,1]
        Gmi = h*Gm[i]
        for j in i+1:N
            xij = xi - uj[1,j,1]
            yij = yi - uj[2,j,1]
            zij = zi - uj[3,j,1]
            Gmj = h*Gm[j]
            dotij = (xij*xij+yij*yij+zij*zij)
            auxij = 1/(sqrt(dotij)*dotij)
            Gmjauxij = Gmj*auxij
            uj[1,i,2] -= Gmjauxij*xij
            uj[2,i,2] -= Gmjauxij*yij
            uj[3,i,2] -= Gmjauxij*zij
            Gmiauxij = Gmi*auxij
            uj[1,j,2] += Gmiauxij*xij
            uj[2,j,2] += Gmiauxij*yij
            uj[3,j,2] += Gmiauxij*zij
        end
     end

    return nothing

end


