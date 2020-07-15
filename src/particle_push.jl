function push_v!( p )

exext = 0.0
eyext = 0.0

for ipart = 1:p.nbpart

   #*** Changement de variable u = gamma*vit

   if( relativ ) then

      u2    =   p.vit[ipart,1]*p.vit[ipart,1]	+ p.vit[ipart,2]*p.vit[ipart,2]

      if ( u2 >= csq )
         println("Erreur : u2 >= c2 dans le calcul de la vitesse")
      else
         gamma = 1/sqrt( 1 - u2/csq )
      end

      p.vit[ipart,1] = gamma*p.vit[ipart,1]
      p.vit[ipart,2] = gamma*p.vit[ipart,2]

   else

       gamma = 1

   end


   #*** Separation des effets electriques et magnetiques

   #*** On ajoute la moitie de l'effet champ p.ctrique E

   dum = 0.5 * dt * p.q_over_m

   p.vit[ipart,1] = p.vit[ipart,1] + dum*(p.epx[ipart]+exext)
   p.vit[ipart,2] = p.vit[ipart,2] + dum*(p.epy[ipart]+eyext)

   #*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (p.bpz[ipart]+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   p.vit[ipart,1] = p.vit[ipart,1] + p.vit[ipart,2] * tantheta
   p.vit[ipart,2] = p.vit[ipart,2] - p.vit[ipart,1] * sintheta
   p.vit[ipart,1] = p.vit[ipart,1] + p.vit[ipart,2] * tantheta

   #*** Autre moitie de l'effet du champ electrique E

   p.vit[ipart,1] = p.vit[ipart,1] + dum*(p.epx[ipart]+exext)
   p.vit[ipart,2] = p.vit[ipart,2] + dum*(p.epy[ipart]+eyext)

   #*** On repasse a la vitesse (changement de variable inverse)

   if ( relativ )

      u2 = p.vit[ipart,1]*p.vit[ipart,1] + p.vit[ipart,2]*p.vit[ipart,2]

      gamma = sqrt( 1. + u2/csq )
 
      p.vit[ipart,1] = p.vit[ipart,1] / gamma
      p.vit[ipart,2] = p.vit[ipart,2] / gamma

   end

end

end 


function push_x!( p, mesh, dt )  # Avancee de coef * dt

    for ipart=1:p.nbpart

       p.pos[ipart,1] = mod(p.pos[ipart,1] + p.vit[ipart,1] * dt, mesh.dimy)
       p.pos[ipart,2] = mod(p.pos[ipart,2] + p.vit[ipart,2] * dt, mesh.dimy)

       p.case[ipart,1] = trunc(Int, p.pos[ipart,1] / dimx * nx) + 1
       p.case[ipart,2] = trunc(Int, p.pos[ipart,2] / dimy * ny) + 1

    end 

end 


function compute_current!( f, p, m, jx, jy )

   fill!(jx, 0)
   fill!(jy, 0)

   one_over_s = m.dx * m.dy
   
   for ipart=1:nbpart

      i = p.case[ipart,1]
      j = p.case[ipart,2]

      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]

      dum = p.p[ipart] * one_over_s

      a1 = (x[i+1]-xp) * (y[j+1]-yp) * dum
      a2 = (xp-x[i]) * (y[j+1]-yp) * dum
      a3 = (xp-x[i]) * (yp-y[j]) * dum
      a4 = (x[i+1]-xp) * (yp-y[j]) * dum

      dum = p.vit[ipart,1] * one_over_s

      jx[i,j]     = jx[i,j]     + a1*dum  
      jx[i+1,j]   = jx[i+1,j]   + a2*dum 
      jx[i+1,j+1] = jx[i+1,j+1] + a3*dum 
      jx[i,j+1]   = jx[i,j+1]   + a4*dum 

      dum = p.vit[ipart,2] * one_over_s

      jy[i,j]     = jy[i,j]     + a1*dum  
      jy[i+1,j]   = jy[i+1,j]   + a2*dum 
      jy[i+1,j+1] = jy[i+1,j+1] + a3*dum 
      jy[i,j+1]   = jy[i,j+1]   + a4*dum 

   end

   for i=1:nx
      jx[i,1]  = jx[i,1] + jx[i,ny]
      jx[i,ny] = jx[i,1]
      jy[i,1]  = jy[i,1] + jy[i,ny]
      jy[i,ny] = jy[i,1]
   end
   for j=1:ny
      jx[1,j]  = jx[1,j] + jx[nx,j]
      jx[nx,j] = jx[1,j]
      jy[1,j]  = jy[1,j] + jy[nx,j]
      jy[nx,j] = jy[1,j]
   end


    for j=1:ny+1, i=1:nx
       f.jx[i,j] = 0.5 * (jx[i,j]+jx[i+1,j])
    end
    
    for j=1:ny, i=1:nx+1
       f.jy[i,j] = 0.5 * (jy[i,j]+jy[i,j+1])
    end

end 
