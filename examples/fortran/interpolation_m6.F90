module m6_interpolation_m

use m6_compute_rho_m
use particles_m
use mesh_fields_m

implicit none

integer, private :: im3
integer, private :: im2
integer, private :: im1
integer, private :: ip1
integer, private :: ip2
integer, private :: ip3
integer, private :: jm3
integer, private :: jm2
integer, private :: jm1
integer, private :: jp1
integer, private :: jp2
integer, private :: jp3

real(8), private :: cm3x
real(8), private :: cp3x
real(8), private :: cm2x
real(8), private :: cp2x
real(8), private :: cm1x
real(8), private :: cp1x
real(8), private :: cx  
real(8), private :: cy  
real(8), private :: cm3y
real(8), private :: cp3y
real(8), private :: cm2y
real(8), private :: cp2y
real(8), private :: cm1y
real(8), private :: cp1y
    
     
contains

subroutine interpolate_eb_m6_complex( e, fields, x, nbpart, ntau )

    real(8),           intent(out) :: e(:,:,:)
    type(fields_2d_t), intent(in)  :: fields
    complex(8),        intent(in)  :: x(:,:,:)
    integer(8),        intent(in)  :: nbpart
    integer,           intent(in)  :: ntau

    real(8)             :: xmin
    real(8)             :: xmax
    real(8)             :: ymin
    real(8)             :: ymax
    real(8)             :: dimx
    real(8)             :: dimy
    integer             :: i
    integer             :: j
    integer             :: k
    integer             :: l
    integer             :: nx
    integer             :: ny
    real(8)             :: dx
    real(8)             :: dy
    real(8)             :: px
    real(8)             :: py
    real(8)             :: dpx
    real(8)             :: dpy
    real(8)             :: s
    integer             :: n

    nx = fields%mesh%nx
    ny = fields%mesh%ny

    dx = fields%mesh%dx
    dy = fields%mesh%dy

    xmin = fields%mesh%xmin
    ymin = fields%mesh%ymin
    xmax = fields%mesh%xmax 
    ymax = fields%mesh%ymax

    dimx = xmax - xmin
    dimy = ymax - ymin

    do k=1,nbpart

        do n=1,ntau
    
            px = real(x(n,1,k))/dx
            py = real(x(n,2,k))/dy

            px = modulo(px, real(nx,8))
            py = modulo(py, real(ny,8))

            i   = floor(px)
            dpx = px - real(i, kind=8)
            j   = floor(py)
            dpy = py - real(j, kind=8)
    
            im3 = modulo(i-3,nx) + 1
            im2 = modulo(i-2,nx) + 1
            im1 = modulo(i-1,nx) + 1
            ip1 = modulo(i+1,nx) + 1
            ip2 = modulo(i+2,nx) + 1
            ip3 = modulo(i+3,nx) + 1
            jm3 = modulo(j-3,ny) + 1
            jm2 = modulo(j-2,ny) + 1
            jm1 = modulo(j-1,ny) + 1
            jp1 = modulo(j+1,ny) + 1
            jp2 = modulo(j+2,ny) + 1
            jp3 = modulo(j+3,ny) + 1

            i = i + 1
            j = j + 1

            cm3x = f_m6(3+dpx)
            cp3x = f_m6(3-dpx)
            cm2x = f_m6(2+dpx)
            cp2x = f_m6(2-dpx)
            cm1x = f_m6(1+dpx)
            cp1x = f_m6(1-dpx)
            cx   = f_m6(  dpx)
            cy   = f_m6(  dpy)
            cm3y = f_m6(3+dpy)
            cp3y = f_m6(3-dpy)
            cm2y = f_m6(2+dpy)
            cp2y = f_m6(2-dpy)
            cm1y = f_m6(1+dpy)
            cp1y = f_m6(1-dpy)
    
     
            do l = 1,2

                s = 0d0
                s = s + cm3x * cm3y * fields%e(l,im3,jm3)   
                s = s + cm3x * cm2y * fields%e(l,im3,jm2)   
                s = s + cm3x * cm1y * fields%e(l,im3,jm1)   
                s = s + cm3x * cy   * fields%e(l,im3,j  )   
                s = s + cm3x * cp1y * fields%e(l,im3,jp1)   
                s = s + cm3x * cp2y * fields%e(l,im3,jp2)   
                s = s + cm3x * cp3y * fields%e(l,im3,jp3)   
                s = s + cm2x * cm3y * fields%e(l,im2,jm3)   
                s = s + cm2x * cm2y * fields%e(l,im2,jm2)   
                s = s + cm2x * cm1y * fields%e(l,im2,jm1)   
                s = s + cm2x * cy   * fields%e(l,im2,j  )   
                s = s + cm2x * cp1y * fields%e(l,im2,jp1)   
                s = s + cm2x * cp2y * fields%e(l,im2,jp2)   
                s = s + cm2x * cp3y * fields%e(l,im2,jp3)   
                s = s + cm1x * cm3y * fields%e(l,im1,jm3)   
                s = s + cm1x * cm2y * fields%e(l,im1,jm2)   
                s = s + cm1x * cm1y * fields%e(l,im1,jm1)   
                s = s + cm1x * cy   * fields%e(l,im1,j  )   
                s = s + cm1x * cp1y * fields%e(l,im1,jp1)   
                s = s + cm1x * cp2y * fields%e(l,im1,jp2)   
                s = s + cm1x * cp3y * fields%e(l,im1,jp3)   
                s = s + cx   * cm3y * fields%e(l,i  ,jm3)   
                s = s + cx   * cm2y * fields%e(l,i  ,jm2)   
                s = s + cx   * cm1y * fields%e(l,i  ,jm1)   
                s = s + cx   * cy   * fields%e(l,i  ,j  )   
                s = s + cx   * cp1y * fields%e(l,i  ,jp1)   
                s = s + cx   * cp2y * fields%e(l,i  ,jp2)   
                s = s + cx   * cp3y * fields%e(l,i  ,jp3)   
                s = s + cp1x * cm3y * fields%e(l,ip1,jm3)   
                s = s + cp1x * cm2y * fields%e(l,ip1,jm2)   
                s = s + cp1x * cm1y * fields%e(l,ip1,jm1)   
                s = s + cp1x * cy   * fields%e(l,ip1,j  )   
                s = s + cp1x * cp1y * fields%e(l,ip1,jp1)   
                s = s + cp1x * cp2y * fields%e(l,ip1,jp2)   
                s = s + cp1x * cp3y * fields%e(l,ip1,jp3)   
                s = s + cp2x * cm3y * fields%e(l,ip2,jm3)   
                s = s + cp2x * cm2y * fields%e(l,ip2,jm2)   
                s = s + cp2x * cm1y * fields%e(l,ip2,jm1)   
                s = s + cp2x * cy   * fields%e(l,ip2,j  )   
                s = s + cp2x * cp1y * fields%e(l,ip2,jp1)   
                s = s + cp2x * cp2y * fields%e(l,ip2,jp2)   
                s = s + cp2x * cp3y * fields%e(l,ip2,jp3)   
                s = s + cp3x * cm3y * fields%e(l,ip3,jm3)   
                s = s + cp3x * cm2y * fields%e(l,ip3,jm2)   
                s = s + cp3x * cm1y * fields%e(l,ip3,jm1)   
                s = s + cp3x * cy   * fields%e(l,ip3,j  )   
                s = s + cp3x * cp1y * fields%e(l,ip3,jp1)   
                s = s + cp3x * cp2y * fields%e(l,ip3,jp2)   
                s = s + cp3x * cp3y * fields%e(l,ip3,jp3)

                e(n,l,k) = s

            end do
    
        end do

    end do

end subroutine interpolate_eb_m6_complex 

subroutine interpolate_eb_m6_real( particles, fields )

    type(particles_t) :: particles
    type(fields_2d_t) :: fields

    integer :: nx
    integer :: ny
    real(8) :: dx
    real(8) :: dy
    real(8) :: dimx
    real(8) :: dimy
    real(8) :: dpx
    real(8) :: dpy
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    real(8) :: px
    real(8) :: py
    real(8) :: e

    nx = fields%mesh%nx
    ny = fields%mesh%ny

    dx = fields%mesh%dx
    dy = fields%mesh%dy

    dimx = fields%mesh%xmax - fields%mesh%xmin
    dimy = fields%mesh%ymax - fields%mesh%ymin

    do k=1,particles%nbpart
    
       px = particles%x(1,k)/dx
       py = particles%x(2,k)/dy

       px = modulo(px, real(nx,8))
       py = modulo(py, real(ny,8))

       i   = floor(px)
       dpx = px - real(i, kind=8)
       j   = floor(py)
       dpy = py - real(j, kind=8)
    
       im3 = modulo(i-3,nx) + 1
       im2 = modulo(i-2,nx) + 1
       im1 = modulo(i-1,nx) + 1
       ip1 = modulo(i+1,nx) + 1
       ip2 = modulo(i+2,nx) + 1
       ip3 = modulo(i+3,nx) + 1
       jm3 = modulo(j-3,ny) + 1
       jm2 = modulo(j-2,ny) + 1
       jm1 = modulo(j-1,ny) + 1
       jp1 = modulo(j+1,ny) + 1
       jp2 = modulo(j+2,ny) + 1
       jp3 = modulo(j+3,ny) + 1

       i = i + 1
       j = j + 1

       cm3x = f_m6(3+dpx)
       cp3x = f_m6(3-dpx)
       cm2x = f_m6(2+dpx)
       cp2x = f_m6(2-dpx)
       cm1x = f_m6(1+dpx)
       cp1x = f_m6(1-dpx)
       cx   = f_m6(  dpx)
       cy   = f_m6(  dpy)
       cm3y = f_m6(3+dpy)
       cp3y = f_m6(3-dpy)
       cm2y = f_m6(2+dpy)
       cp2y = f_m6(2-dpy)
       cm1y = f_m6(1+dpy)
       cp1y = f_m6(1-dpy)
    
     
       do l = 1,2

           e = 0d0
           e = e + cm3x * cm3y * fields%e(l,im3,jm3)   
           e = e + cm3x * cm2y * fields%e(l,im3,jm2)   
           e = e + cm3x * cm1y * fields%e(l,im3,jm1)   
           e = e + cm3x * cy   * fields%e(l,im3,j  )   
           e = e + cm3x * cp1y * fields%e(l,im3,jp1)   
           e = e + cm3x * cp2y * fields%e(l,im3,jp2)   
           e = e + cm3x * cp3y * fields%e(l,im3,jp3)   
           e = e + cm2x * cm3y * fields%e(l,im2,jm3)   
           e = e + cm2x * cm2y * fields%e(l,im2,jm2)   
           e = e + cm2x * cm1y * fields%e(l,im2,jm1)   
           e = e + cm2x * cy   * fields%e(l,im2,j  )   
           e = e + cm2x * cp1y * fields%e(l,im2,jp1)   
           e = e + cm2x * cp2y * fields%e(l,im2,jp2)   
           e = e + cm2x * cp3y * fields%e(l,im2,jp3)   
           e = e + cm1x * cm3y * fields%e(l,im1,jm3)   
           e = e + cm1x * cm2y * fields%e(l,im1,jm2)   
           e = e + cm1x * cm1y * fields%e(l,im1,jm1)   
           e = e + cm1x * cy   * fields%e(l,im1,j  )   
           e = e + cm1x * cp1y * fields%e(l,im1,jp1)   
           e = e + cm1x * cp2y * fields%e(l,im1,jp2)   
           e = e + cm1x * cp3y * fields%e(l,im1,jp3)   
           e = e + cx   * cm3y * fields%e(l,i  ,jm3)   
           e = e + cx   * cm2y * fields%e(l,i  ,jm2)   
           e = e + cx   * cm1y * fields%e(l,i  ,jm1)   
           e = e + cx   * cy   * fields%e(l,i  ,j  )   
           e = e + cx   * cp1y * fields%e(l,i  ,jp1)   
           e = e + cx   * cp2y * fields%e(l,i  ,jp2)   
           e = e + cx   * cp3y * fields%e(l,i  ,jp3)   
           e = e + cp1x * cm3y * fields%e(l,ip1,jm3)   
           e = e + cp1x * cm2y * fields%e(l,ip1,jm2)   
           e = e + cp1x * cm1y * fields%e(l,ip1,jm1)   
           e = e + cp1x * cy   * fields%e(l,ip1,j  )   
           e = e + cp1x * cp1y * fields%e(l,ip1,jp1)   
           e = e + cp1x * cp2y * fields%e(l,ip1,jp2)   
           e = e + cp1x * cp3y * fields%e(l,ip1,jp3)   
           e = e + cp2x * cm3y * fields%e(l,ip2,jm3)   
           e = e + cp2x * cm2y * fields%e(l,ip2,jm2)   
           e = e + cp2x * cm1y * fields%e(l,ip2,jm1)   
           e = e + cp2x * cy   * fields%e(l,ip2,j  )   
           e = e + cp2x * cp1y * fields%e(l,ip2,jp1)   
           e = e + cp2x * cp2y * fields%e(l,ip2,jp2)   
           e = e + cp2x * cp3y * fields%e(l,ip2,jp3)   
           e = e + cp3x * cm3y * fields%e(l,ip3,jm3)   
           e = e + cp3x * cm2y * fields%e(l,ip3,jm2)   
           e = e + cp3x * cm1y * fields%e(l,ip3,jm1)   
           e = e + cp3x * cy   * fields%e(l,ip3,j  )   
           e = e + cp3x * cp1y * fields%e(l,ip3,jp1)   
           e = e + cp3x * cp2y * fields%e(l,ip3,jp2)   
           e = e + cp3x * cp3y * fields%e(l,ip3,jp3)

           particles%e(l,k) = e
    
        end do

    end do
    
end subroutine interpolate_eb_m6_real

end module m6_interpolation_m
