function interpol_eb!( ex, ey, bz, p :: Particles, m :: Mesh )

    for ipart=1:p.nbpart

        i  = p.case(ipart,1)
        j  = p.case(ipart,2)
        xp = p.pos(ipart,1)
        yp = p.pos(ipart,2)
    
        dum = 1/(m.hx[i]*m.hy[j])

        a1 = (x[i+1]-xp) * (y(j+1)-yp) * dum
        a2 = (xp-x[i]) * (y(j+1)-yp) * dum
        a3 = (xp-x[i]) * (yp-y[j]) * dum
        a4 = (x[i+1]-xp) * (yp-y[j]) * dum
    
        ex[ipart] = a1 * ex[i,j] + a2 * ex[i+1,j] + a3 * ex[i+1,j+1] + a4 * ex[i,j+1] 
        ey[ipart] = a1 * ey[i,j] + a2 * ey[i+1,j] + a3 * ey[i+1,j+1] + a4 * ey[i,j+1] 
        bz[ipart] = a1 * bz[i,j] + a2 * bz[i+1,j] + a3 * bz[i+1,j+1] + a4 * bz[i,j+1] 

    end

end 


