function push_v!( p, dt )

    for ipart = 1:p.nbpart
    
        v1 = p.vit[ipart,1]
        v2 = p.vit[ipart,2]
        e1 = p.epx[ipart]
        e2 = p.epy[ipart]

        # u = gamma*vit
    
        if( relativ ) then
    
           u2 = v1*v1+v2*v2
    
           if ( u2 >= csq )
              println("Erreur : u2 >= c2 dans le calcul de la vitesse")
           else
              gamma = 1/sqrt( 1 - u2/csq )
           end
    
           v1 = gamma*v1
           v2 = gamma*v2
    
        else
    
            gamma = 1
    
        end
    
        #*** Separation des effets electriques et magnetiques
        #*** On ajoute la moitie de l'effet champ p.ctrique E
    
        dum = 0.5 * dt * p.q_over_m
    
        v1 = v1 + dum * e1
        v2 = v2 + dum * e2
    
        #*** Algorithme de Buneman pour les effets magnetiques
     
        tantheta = dum * p.bpz[ipart] / gamma 
        sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)
    
        v1 = v1 + v2 * tantheta
        v2 = v2 - v1 * sintheta
        v1 = v1 + v2 * tantheta
    
        #*** Autre moitie de l'effet du champ electrique E
    
        v1 = v1 + dum * e1
        v2 = v2 + dum * e2
    
        #*** On repasse a la vitesse (changement de variable inverse)
    
        if ( relativ )
    
           u2 = v1*v1 + v2*v2
    
           gamma = sqrt( 1. + u2/csq )
     
           v1 = v1 / gamma
           v2 = v2 / gamma
    
        end

        p.vit[ipart, 1] = v1
        p.vit[ipart, 2] = v2
    
    end

end 


function push_x!( p, dt )

   p.pos .+= p.vit .* dt

end 


