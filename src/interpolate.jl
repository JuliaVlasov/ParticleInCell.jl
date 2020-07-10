
function interpol_eb( m, f, p )

    for ipart=1:p.nbpart

        i = p.case(ipart,1)
        j = p.case(ipart,2)
        xp = p.pos(ipart,1)
        yp = p.pos(ipart,2)
    
        dum = 1/(m.hx[i]*m.hy[j])
        a1 = (x[i+1]-xp) * (y(j+1)-yp) * dum
        a2 = (xp-x[i]) * (y(j+1)-yp) * dum
        a3 = (xp-x[i]) * (yp-y[j]) * dum
        a4 = (x[i+1]-xp) * (yp-y[j]) * dum
    
        p.epx[ipart] = a1 * f.ex[i,j] + a2 * f.ex[i+1,j] + a3 * f.ex[i+1,j+1] + a4 * f.ex[i,j+1] 
        p.epy[ipart] = a1 * f.ey[i,j] + a2 * f.ey[i+1,j] + a3 * f.ey[i+1,j+1] + a4 * f.ey[i,j+1] 
        p.bpz[ipart] = a1 * f.bz[i,j] + a2 * f.bz[i+1,j] + a3 * f.bz[i+1,j+1] + a4 * f.bz[i,j+1] 

    end

end 


function avancee_vitesse( p )

for ipart = 1, nbpart

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

   end if


   #*** Separation des effets p.ctriques et magnetiques

   #*** On ajoute la moitie de l'effet champ p.ctrique E

   dum = 0.5 * dt * q_sur_m
   p.vit[ipart,1] = p.vit[ipart,1] + dum*(p.epx[ipart]+exext)
   p.vit[ipart,2] = p.vit[ipart,2] + dum*(p.epy[ipart]+eyext)

   #*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (p.bpz[ipart]+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   p.vit[ipart,1] = p.vit[ipart,1] + p.vit[ipart,2] * tantheta
   p.vit[ipart,2] = p.vit[ipart,2] - p.vit[ipart,1] * sintheta
   p.vit[ipart,1] = p.vit[ipart,1] + p.vit[ipart,2] * tantheta

   #*** Autre moitie de l'effet du champ p.ctrique E

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


function avancee_part( p, coef )  # Avancee de coef * dt

    for ipart=1:p.nbpart
       p.pos[ipart,1] = p.pos[ipart,1] + p.vit[ipart,1] * dt * coef
       p.pos[ipart,2] = p.pos[ipart,2] + p.vit[ipart,2] * dt * coef
    end
    
    #*** Mise a jour des "cases"
    
    for ipart=1:p.nbpart
       i = 0
       while (p.pos(ipart,1) >= x[i] .and. p.pos(ipart,1)<dimx) 
          i=i+1
       end
       if ( p.pos(ipart,1) >= dimx ) i=nx+1
       p.case(ipart,1) = i-1
       j = 0
       while (p.pos(ipart,2) >= y[j] .and. p.pos(ipart,2)<dimy) 
          j=j+1 
       end
       if ( p.pos(ipart,2) >= dimy ) j=ny+1
       p.case(ipart,2) = j-1
    end 

end 


function sortie_part( p, m )

   for ipart=1:nbpart
      if( p.pos[ipart,1] >= dimx ) p.pos[ipart,1] = p.pos[ipart,1] - m.dimx
      if( p.pos[ipart,2] >= dimy ) p.pos[ipart,2] = p.pos[ipart,2] - m.dimy
      if( p.pos[ipart,1] < 0 )  p.pos[ipart,1] = p.pos[ipart,1] + m.dimx
      if( p.pos[ipart,2] < 0 )  p.pos[ipart,2] = p.pos[ipart,2] + m.dimy
   end

   for ipart=1:nbpart
       i = 0
       while (p.pos[ipart,1] >= x[i] && p.pos[ipart,1] <= dimx) 
          i=i+1
       end
       p.case[ipart,1] = i-1
       j = 0
       while (p.pos[ipart,2] >= y[j] && p.pos[ipart,2] <= dimy) 
          j=j+1 
       end
       p.case[ipart,2] = j-1
   end

end 

function calcul_j_cic( m, p, f, jx, jy )

   fill!(jx, 0)
   fill!(jy, 0)
   
   for ipart=1:nbpart
      i = p.case[ipart,1]
      j = p.case[ipart,2]
      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]
      dum = p.p[ipart] / (m.hx[i]*m.hy[j])
      a1 = (x[i+1]-xp) * (y(j+1)-yp) * dum
      a2 = (xp-x[i]) * (y(j+1)-yp) * dum
      a3 = (xp-x[i]) * (yp-y[j]) * dum
      a4 = (x[i+1]-xp) * (yp-y[j]) * dum
      dum = p.vit[ipart,1] / (m.hx[i]*m.hy[j]) 
      jx[i,j]     = jx[i,j]     + a1*dum  
      jx[i+1,j]   = jx[i+1,j]   + a2*dum 
      jx[i+1,j+1] = jx[i+1,j+1] + a3*dum 
      jx[i,j+1]   = jx[i,j+1]   + a4*dum 
      dum = p.vit[ipart,2] / (m.hx[i]*m.hy[j]) 
      jy[i,j]     = jy[i,j]     + a1*dum  
      jy[i+1,j]   = jy[i+1,j]   + a2*dum 
      jy[i+1,j+1] = jy[i+1,j+1] + a3*dum 
      jy[i,j+1]   = jy[i,j+1]   + a4*dum 
   end

   for i=0:nx
      jx[i,0]  = jx[i,0] + jx[i,ny]
      jx[i,ny] = jx[i,0]
      jy[i,0]  = jy[i,0] + jy[i,ny]
      jy[i,ny] = jy[i,0]
   end
   for j=0:ny
      jx[0,j]  = jx[0,j] + jx[nx,j]
      jx[nx,j] = jx[0,j]
      jy[0,j]  = jy[0,j] + jy[nx,j]
      jy[nx,j] = jy[0,j]
   end


    for i=0:nx-1
    for j=0:ny
       f.jx[i,j] = 0.5 * (jx[i,j]+jx[i+1,j])
    end 
    end
    
    for i=0:nx
    for j=0:ny-1
       f.jy[i,j] = 0.5 * (jy[i,j]+jy[i,j+1])
    end
    end

end 
