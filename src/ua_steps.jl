export preparation!

function preparation!( ua         :: UA, 
                       dt         :: Float64, 
                       particles  :: Particles, 
                       xt         :: Array{ComplexF64,3}, 
                       yt         :: Array{ComplexF64,3}) 

    ε      = ua.ε
    ntau   = ua.ntau
    nbpart = particles.nbpart
    tau    = ua.tau
    ltau   = ua.ltau
    
    for m = 1:nbpart

        x1 = particles.x[1,m]
        x2 = particles.x[2,m]

        particles.b[m] = 1 + 0.5 * sin(x1) * sin(x2)
        particles.t[m] = dt * particles.b[m]

        ua.pl[1,m] = particles.t[m]
        ua.ql[1,m] = particles.t[m]^2 / 2

        t = particles.t[m]
        b = particles.b[m]

        for n=2:ntau
            elt = exp(-1im*ua.ltau[n]*t/ε) 
            ua.pl[n,m] = ε * 1im*(elt-1)/ua.ltau[n]
            ua.ql[n,m] = ε * (ε*(1-elt) -1im*ua.ltau[n]*t)/ua.ltau[n]^2
        end

        ex,  ey  = particles.e[1:2,m]
        vx,  vy  = particles.v[1:2,m]
        vxb, vyb = vx/b, vy/b

        for n = 1:ntau

            τ = tau[n]

            h1 = ε * (sin(τ) * vxb - cos(τ) * vyb)
            h2 = ε * (sin(τ) * vyb + cos(τ) * vxb)

            xt1 = x1 + h1 + ε * vyb
            xt2 = x2 + h2 - ε * vxb

            xt[n,1,m] = xt1
            xt[n,2,m] = xt2

            interv=(1+0.5*sin(xt1)*sin(xt2)-b)/ε

            exb = ((  cos(τ)*vy - sin(τ)*vx) * interv + ex)/b
            eyb = (( -cos(τ)*vx - sin(τ)*vy) * interv + ey)/b

            ua.r[n,1] = cos(τ)* exb - sin(τ) * eyb
            ua.r[n,2] = sin(τ)* exb + cos(τ) * eyb

        end

        mul!(ua.r̃, ua.rtau, ua.r)

        for n = 2:ntau
            ua.r̃[n,1] = -1im * ua.r̃[n,1]/ltau[n]
            ua.r̃[n,2] = -1im * ua.r̃[n,2]/ltau[n]
        end

        ldiv!(ua.r, ua.rtau, ua.r̃)

        for n = 1:ntau
            yt[n,1,m] = vx + (ua.r[n,1] - ua.r[1,1]) * ε
            yt[n,2,m] = vy + (ua.r[n,2] - ua.r[1,2]) * ε
        end

    end

end

export update_particles_e!

function update_particles_e!( particles :: Particles, 
                              et        :: Array{Float64,3}, 
                              fields    :: MeshFields, 
                              ua        :: UA, 
                              xt        :: Array{ComplexF64,3}) 

    interpol_eb_m6!( et, fields, xt, particles.nbpart, ua.ntau) 

end

export update_particles_x!

function update_particles_x!( particles :: Particles, 
                              fields    :: MeshFields, 
                              ua        :: UA, 
                              xt        :: Array{ComplexF64,3}) 

    compute_rho_m6!( fields, particles, xt, ua )

end

export compute_f!

function compute_f!( fx        :: Array{ComplexF64,3}, 
                     fy        :: Array{ComplexF64,3}, 
                     ua        :: UA, 
                     particles :: Particles, 
                     xt        :: Array{ComplexF64,3}, 
                     yt        :: Array{ComplexF64,3},
                     et        :: Array{Float64,3} )

    for m=1:particles.nbpart
    
        b = particles.b[m]

        for n=1:ua.ntau
    
            xt1 = real(xt[n,1,m])
            xt2 = real(xt[n,2,m])
    
            yt1 = yt[n,1,m] 
            yt2 = yt[n,2,m] 
    
            τ = ua.tau[n]
    
            fx[n,1,m] = ( cos(τ) * yt1 + sin(τ) * yt2)/b
            fx[n,2,m] = (-sin(τ) * yt1 + cos(τ) * yt2)/b
    
            interv = (1 + 0.5*sin(xt1)*sin(xt2)-b)/ua.ε
    
            tmp1 = et[n,1,m]+(  cos(τ)*yt2 - sin(τ)*yt1)*interv
            tmp2 = et[n,2,m]+(- cos(τ)*yt1 - sin(τ)*yt2)*interv
    
            fy[n,1,m] = (cos(τ)*tmp1-sin(τ)*tmp2)/b
            fy[n,2,m] = (sin(τ)*tmp1+cos(τ)*tmp2)/b
    
        end
    
    end

    fft!(fx,1)
    fft!(fy,1)

end

export ua_step!

function ua_step!( xt        :: Array{ComplexF64, 3}, 
                   x̃t        :: Array{ComplexF64, 3}, 
                   ua        :: UA, 
                   particles :: Particles, 
                   fx        :: Array{ComplexF64, 3} )


    for m=1:particles.nbpart

        t = particles.t[m]

        for n=1:ua.ntau

            elt = exp(-1im*ua.ltau[n]*t/ua.ε) 
            xt[n,1,m] = elt * x̃t[n,1,m] + ua.pl[n,m] * fx[n,1,m]
            xt[n,2,m] = elt * x̃t[n,2,m] + ua.pl[n,m] * fx[n,2,m]

        end

    end

end

function ua_step!( xt        :: Array{ComplexF64, 3}, 
                   x̃t        :: Array{ComplexF64, 3}, 
                   ua        :: UA, 
                   particles :: Particles, 
                   fx        :: Array{ComplexF64, 3},
                   gx        :: Array{ComplexF64, 3} )

     for m=1:particles.nbpart

         t = particles.t[m]

         for n=1:ua.ntau

             elt = exp(-1im*ua.ltau[n]*t/ua.ε) 
             fx1 = fx[n,1,m]
             fx2 = fx[n,2,m]
             xt1 = elt * x̃t[n,1,m] + ua.pl[n,m] * fx1
             xt2 = elt * x̃t[n,2,m] + ua.pl[n,m] * fx2
             xt1 += ua.ql[n,m] * (gx[n,1,m] - fx1) / t
             xt2 += ua.ql[n,m] * (gx[n,2,m] - fx2) / t

             xt[n,1,m] = xt1
             xt[n,2,m] = xt2

         end

    end

end

export compute_v!

function compute_v!( yt       :: Array{ComplexF64, 3},
                    particles :: Particles, 
                    ua        :: UA )

        for m=1:particles.nbpart

            t = particles.t[m]

            px, py = 0.0im, 0.0im
            for n = 1:ua.ntau
                elt = exp(1im*ua.ltau[n]*t/ua.ε) 
                px += yt[n,1,m]/ua.ntau * elt
                py += yt[n,2,m]/ua.ntau * elt
            end

            particles.v[1,m] = real(cos(t/ua.ε)*px+sin(t/ua.ε)*py)
            particles.v[2,m] = real(cos(t/ua.ε)*py-sin(t/ua.ε)*px)

        end

end
