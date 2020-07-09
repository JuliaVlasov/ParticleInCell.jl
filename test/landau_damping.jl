using OffsetArrays
using Sobol
using Plots

@testset "Landau Damping" begin

struct MeshFields

   ex
   ey
   bz
   jx
   jy
   r0
   r1

   function MeshFields(nx, ny)

       ex = OffsetArray{Float64}(undef, 0:nx-1,0:ny)
       ey = OffsetArray{Float64}(undef, 0:nx,0:ny-1)
       bz = OffsetArray{Float64}(undef, 0:nx-1,0:ny-1)
       jx = OffsetArray{Float64}(undef, 0:nx-1,0:ny)
       jy = OffsetArray{Float64}(undef, 0:nx,0:ny-1)
       r0 = OffsetArray{Float64}(undef, 0:nx,0:ny) 
       r1 = OffsetArray{Float64}(undef, 0:nx,0:ny)

       new( ex, ey, bz, jx, jy, r0, r1)

   end

end


nx        = 128		    # nombre de pts suivant x
ny        = 16   		# nombre de pts suivant y
cfl       = 0.9 		# nombre de Courant-Friedrich-Levy
tfinal    = 50      	# temps final
nstepmax  = 20000000	# nbre d'iterations maxi
nomcas    = "plasma"	# nom du cas traite (6 lettres)
jname     = "jcico1"	# methode de calcul des courants
icrea     = 1			# frequence d'emission des particules
idiag     = 50    	    # frequence des sorties graphiques
bcname    = "period"	# type de conditions limites
exext     = 0			# champ electrique exterieur suivant x
eyext     = 0	 		# champ electriaue exterieur suivant y
bzext     = 0			# champ magnetique exterieur
charge    = 1           # charge d'une macro particule
masse     = 1           # masse d'une macro particule
c         = 8	 	    # vitesse de la lumiere
e0        = 1           # permittivite du vide

csq = c*c
q_sur_m = charge / masse
poids = charge

alpha = 0.1
kx = 0.5
ky = 0.
dimx = 2*pi/kx
dimy = 1  
poids = dimx * dimy 

x   = OffsetArray{Float64}(undef, -1:nx+1) 
y   = OffsetArray{Float64}(undef, -1:ny+1)
hx  = OffsetArray{Float64}(undef, -1:nx)
hy  = OffsetArray{Float64}(undef, -1:ny)
hhx = OffsetArray{Float64}(undef,  0:nx)
hhy = OffsetArray{Float64}(undef,  0:ny)

dx = dimx / nx
dy = dimy / ny

x[0] = 0.
y[0] = 0.

for i=1:nx
    x[i] = (i*dx) *(i*dx+1)/(1+dimx)
end

for j=1:ny
    y[j] = (j*dy) *(j*dy+1)/(1+dimy)
end

for i=0:nx-1
    hx[i] = x[i+1]-x[i]
end

for j=0:ny-1
    hy[j] = y[j+1]-y[j]
end

hx[nx] = hx[0]  
hx[-1] = hx[nx-1]
hy[ny] = hy[0]
hy[-1] = hy[ny-1]

x[-1]   = x[0] - hx[nx-1]  # points utiles pour le cas period
x[nx+1] = x[nx] + hx[0]
y[-1]   = y[0] - hy[ny-1]
y[ny+1] = y[ny] + hy[0]

hhx[0]  =  0.5 * ( hx[0] + hx[nx-1] ) 
hhx[nx] =  0.5 * ( hx[0] + hx[nx-1] ) 
for i=1:nx-1
   hhx[i] = 0.5 * ( hx[i] + hx[i-1] )
end

hhy[0]  = 0.5 * ( hy[0] + hy[ny-1] ) 
hhy[ny] = 0.5 * ( hy[0] + hy[ny-1] ) 

for j=1:ny-1
   hhy[j] = 0.5 * ( hy[j] + hy[j-1] )
end

dx = hx[0]
for i=1:nx-1
   (hx[i]<dx)  && (dx = hx[i])
end

dy = hy[0]
for j=1:ny-1
   (hy[j]<dy)  && ( y = hy[j] )
end

dt    = cfl  / sqrt(1/(dx*dx)+1/(dy*dy)) / c

nstep = trunc(Int, tfinal ÷ dt)

println(" cfl = $cfl ")
println(" dx = $dx dy = $dy dt = $dt")
println(" Nombre d'iteration nstep $nstep ")

f = MeshFields(nx, ny)

ex = OffsetArray{Float64}(undef, 0:nx,0:ny)
ey = OffsetArray{Float64}(undef, 0:nx,0:ny)
bz = OffsetArray{Float64}(undef, 0:nx,0:ny)


time  = 0
iplot = 0

nstep > nstepmax && ( nstep = nstepmax )

istep = 1

fill!(f.ex, 0)
fill!(f.ey, 0)
fill!(f.bz, 0)
fill!(f.jx, 0)
fill!(f.jy, 0)
fill!(f.r0, 0)
fill!(f.r1, 0)

for i=0:nx-1
    aux1 = alpha/kx * sin(kx*x[i])
    aux2 = alpha * cos(kx*x[i])
    for j=0:ny
        f.ex[i,j] = aux1
        f.r1[i,j] = aux2
    end
end
      
struct Particles

    nbpart
    pos
    case
    vit
    epx
    epy
    bpz
    p

    function Particles( nbpart )

        pos = zeros(nbpart,2)
        case = zeros(nbpart,2)
        vit = zeros(nbpart,2)
        epx = zeros(nbpart)
        epy = zeros(nbpart)
        bpz = zeros(nbpart)
        p = zeros(nbpart)

        new( nbpart, pos, case, vit, epx, epy, bpz, p )

    end 

end

nbpart = 100*nx*ny

pg = Particles(nbpart)


function landau_sampling!( pg, alpha, kx )

    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1-x0) > 1e-12)
            p = x0 + alpha * sin( kx * x0) / kx
            f = 1 + alpha * cos( kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end

    s = Sobol.SobolSeq(3)
    nbpart = pg.nbpart

    for i=1:nbpart

        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2, r3 = Sobol.next!(s)
        θ = r1 * 2π
        pg.pos[i,1] = newton(r2)
        pg.pos[i,2] = r3 * dimy
        pg.vit[i,1] = v * cos(θ)
        pg.vit[i,2] = v * sin(θ)
        pg.p[i] = 1 / nbpart
    end

end

landau_sampling!( pg, alpha, kx )

xp = view(pg.pos, :, 1)
vp = pg.vit

pp = plot(layout=(3,1))
histogram!(pp[1,1], xp, normalize=true, bins = 100, lab="x")
plot!(pp[1,1], x -> (1+alpha*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
histogram!(pp[2,1], vp[:,1], normalize=true, bins = 100, lab="vx")
plot!(pp[2,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
histogram!(pp[3,1], vp[:,2], normalize=true, bins = 100, lab="vy")
plot!(pp[3,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
savefig("particles.png")

#=

for istep = 1:nstep

   if (istep > 1) then
      call faraday( f0 ) 	!Calcul de B(n-1/2) --> B(n)			
   end if

   call decalage( f0, f1 )
   call interpol_eb( f1, p )

   write(*,*) p%pos(1,1:2), p%vit(1,1:2)
   call avancee_vitesse( p )

   if (jname == 'jcico1') then
      call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
      call sortie_part( p )
      call calcul_j_cic( p, f0 )
      call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
      call sortie_part( p )
   else if (jname == 'jcoco1') then
      call avancee_part( p, 1.d0 )
      call calcul_j_villa( p, f0 )
      call sortie_part( p )
   else
      call avancee_part( p, 1.d0 )
      call sortie_part( p )
   end if
        
   !call calcul_rho( p, f0 )

   call faraday( f0 )   !Calcul de B(n) --> B(n+1/2)
   call ampere( f0 )    !Calcul de E(n) --> E(n+1)
   call conditions_limites( f0, time )

   time = time + dt
   print*,'time = ',time, ' nbpart = ', nbpart

   if ( istep==1 .or. mod(istep,idiag) == 0 .or. istep==nstep ) then
      iplot = iplot + 1
     ! call diag_coc( f0, p, time, iplot )
     ! call diag_champ_part( p, time, iplot )
     ! call plot_champ( f0, iplot, time )
     ! call plot_phases( p, iplot, time )
     ! call distribution_v( p, iplot, time )  
     ! call distribution_x( p, iplot, time )
      call modeE( f0, iplot, time )
   endif

end 

=#


end
