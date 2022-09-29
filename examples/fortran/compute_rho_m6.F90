module m6_compute_rho_m

use ua_type_m
use particles_m
use mesh_fields_m

implicit none

integer, private :: im3, im2, im1, ip1, ip2, ip3
integer, private :: jm3, jm2, jm1, jp1, jp2, jp3
real(8), private :: cm3x, cp3x, cm2x, cp2x, cm1x, cp1x, cx
real(8), private :: cy, cp1y, cm1y, cp2y, cm2y, cp3y, cm3y
real(8), private :: dpx, dpy
real(8), private :: px
real(8), private :: py
real(8), private :: dx
real(8), private :: dy
integer, private :: nx
integer, private :: ny
real(8), private :: weight
real(8), private :: rho_total
integer, private :: i, j, k
real(8), private :: dimx
real(8), private :: dimy

contains

pure function f_m6( q )

    real(8), intent(in)  :: q
    real(8) :: f_m6

    if ( q < 1d0 )  then
        f_m6 = (3d0-q)**5-6d0*(2d0-q)**5+15d0*(1d0-q)**5
    else if ( q >= 1d0 .and. q < 2d0 ) then
        f_m6 = (3d0-q)**5-6*(2d0-q)**5
    else if ( q >= 2d0 .and. q < 3d0 ) then
        f_m6 = (3d0-q)**5
    else
    f_m6 = 0d0
    end if
    
    f_m6 = f_m6 / 120d0

end function f_m6 

subroutine compute_rho_m6_complex( fields, particles, xt, ua ) 

    type(fields_2d_t)  :: fields
    type(particles_t)  :: particles
    complex(8)         :: xt(:,:,:)
    type(ua_t)         :: ua

    complex(8)         :: t
    real(8)            :: xt1
    real(8)            :: xt2

    fields%rho = 0d0

    nx = fields%mesh%nx
    ny = fields%mesh%ny
    nx = fields%mesh%nx
    ny = fields%mesh%ny
    dx = fields%mesh%dx
    dy = fields%mesh%dy

    dimx = fields%mesh%xmax - fields%mesh%xmin
    dimy = fields%mesh%ymax - fields%mesh%ymin
    
    do k = 1,particles%nbpart
    
        t = cmplx(particles%t(k), 0d0, kind=8)

        call fftw_execute_dft(ua%fw, xt(:,1,k), ua%ft) 

        ua%ft = ua%ft * exp(cmplx(0d0,1d0,kind=8)*ua%ltau*t/ua%eps)/cmplx(ua%ntau,0d0,kind=8)

        xt1 = real(sum(ua%ft))

        call fftw_execute_dft(ua%fw, xt(:,2,k), ua%ft) 

        ua%ft = ua%ft * exp(cmplx(0d0,1d0,kind=8)*ua%ltau*t/ua%eps)/cmplx(ua%ntau,0d0,kind=8)

        xt2 = real(sum(ua%ft))

        particles%x(1,k) = xt1
        particles%x(2,k) = xt2

        px = xt1/dx
        py = xt2/dy

        px = modulo(px, real(nx,8))
        py = modulo(py, real(ny,8))

        i   = floor(px)
        dpx = px - real(i ,kind=8)
        j   = floor(py)
        dpy = py - real(j ,kind=8)

        weight = particles%w
      
        im3 = modulo(i -3,nx)+1
        im2 = modulo(i -2,nx)+1
        im1 = modulo(i -1,nx)+1
        ip1 = modulo(i +1,nx)+1
        ip2 = modulo(i +2,nx)+1
        ip3 = modulo(i +3,nx)+1
        jm3 = modulo(j -3,ny)+1
        jm2 = modulo(j -2,ny)+1
        jm1 = modulo(j -1,ny)+1
        jp1 = modulo(j +1,ny)+1
        jp2 = modulo(j +2,ny)+1
        jp3 = modulo(j +3,ny)+1

        i  = i +1
        j  = j +1
      
        cm3x = f_m6(3+dpx)
        cp3x = f_m6(3-dpx)
        cm2x = f_m6(2+dpx)
        cp2x = f_m6(2-dpx)
        cm1x = f_m6(1+dpx)
        cp1x = f_m6(1-dpx)
        cx   = f_m6(dpx)
        cy   = f_m6(dpy)
        cp1y = f_m6(1-dpy)
        cm1y = f_m6(1+dpy)
        cp2y = f_m6(2-dpy)
        cm2y = f_m6(2+dpy)
        cp3y = f_m6(3-dpy)
        cm3y = f_m6(3+dpy)
      
        fields%rho(im3,jm3) = fields%rho(im3,jm3) + cm3x * cm3y * weight
        fields%rho(im3,jm2) = fields%rho(im3,jm2) + cm3x * cm2y * weight
        fields%rho(im3,jm1) = fields%rho(im3,jm1) + cm3x * cm1y * weight
        fields%rho(im3,j  ) = fields%rho(im3,j  ) + cm3x * cy   * weight
        fields%rho(im3,jp1) = fields%rho(im3,jp1) + cm3x * cp1y * weight
        fields%rho(im3,jp2) = fields%rho(im3,jp2) + cm3x * cp2y * weight
        fields%rho(im3,jp3) = fields%rho(im3,jp3) + cm3x * cp3y * weight
                                                   
        fields%rho(im2,jm3) = fields%rho(im2,jm3) + cm2x * cm3y * weight
        fields%rho(im2,jm2) = fields%rho(im2,jm2) + cm2x * cm2y * weight
        fields%rho(im2,jm1) = fields%rho(im2,jm1) + cm2x * cm1y * weight
        fields%rho(im2,j  ) = fields%rho(im2,j  ) + cm2x * cy   * weight
        fields%rho(im2,jp1) = fields%rho(im2,jp1) + cm2x * cp1y * weight
        fields%rho(im2,jp2) = fields%rho(im2,jp2) + cm2x * cp2y * weight
        fields%rho(im2,jp3) = fields%rho(im2,jp3) + cm2x * cp3y * weight
                                                   
        fields%rho(im1,jm3) = fields%rho(im1,jm3) + cm1x * cm3y * weight
        fields%rho(im1,jm2) = fields%rho(im1,jm2) + cm1x * cm2y * weight
        fields%rho(im1,jm1) = fields%rho(im1,jm1) + cm1x * cm1y * weight
        fields%rho(im1,j  ) = fields%rho(im1,j  ) + cm1x * cy   * weight
        fields%rho(im1,jp1) = fields%rho(im1,jp1) + cm1x * cp1y * weight
        fields%rho(im1,jp2) = fields%rho(im1,jp2) + cm1x * cp2y * weight
        fields%rho(im1,jp3) = fields%rho(im1,jp3) + cm1x * cp3y * weight
                                                   
        fields%rho(i  ,jm3) = fields%rho(i  ,jm3) + cx   * cm3y * weight
        fields%rho(i  ,jm2) = fields%rho(i  ,jm2) + cx   * cm2y * weight
        fields%rho(i  ,jm1) = fields%rho(i  ,jm1) + cx   * cm1y * weight
        fields%rho(i  ,j  ) = fields%rho(i  ,j  ) + cx   * cy   * weight
        fields%rho(i  ,jp1) = fields%rho(i  ,jp1) + cx   * cp1y * weight
        fields%rho(i  ,jp2) = fields%rho(i  ,jp2) + cx   * cp2y * weight
        fields%rho(i  ,jp3) = fields%rho(i  ,jp3) + cx   * cp3y * weight
                                                   
        fields%rho(ip1,jm3) = fields%rho(ip1,jm3) + cp1x * cm3y * weight
        fields%rho(ip1,jm2) = fields%rho(ip1,jm2) + cp1x * cm2y * weight
        fields%rho(ip1,jm1) = fields%rho(ip1,jm1) + cp1x * cm1y * weight
        fields%rho(ip1,j  ) = fields%rho(ip1,j  ) + cp1x * cy   * weight
        fields%rho(ip1,jp1) = fields%rho(ip1,jp1) + cp1x * cp1y * weight
        fields%rho(ip1,jp2) = fields%rho(ip1,jp2) + cp1x * cp2y * weight
        fields%rho(ip1,jp3) = fields%rho(ip1,jp3) + cp1x * cp3y * weight
                                                   
        fields%rho(ip2,jm3) = fields%rho(ip2,jm3) + cp2x * cm3y * weight
        fields%rho(ip2,jm2) = fields%rho(ip2,jm2) + cp2x * cm2y * weight
        fields%rho(ip2,jm1) = fields%rho(ip2,jm1) + cp2x * cm1y * weight
        fields%rho(ip2,j  ) = fields%rho(ip2,j  ) + cp2x * cy   * weight
        fields%rho(ip2,jp1) = fields%rho(ip2,jp1) + cp2x * cp1y * weight
        fields%rho(ip2,jp2) = fields%rho(ip2,jp2) + cp2x * cp2y * weight
        fields%rho(ip2,jp3) = fields%rho(ip2,jp3) + cp2x * cp3y * weight
                                                   
        fields%rho(ip3,jm3) = fields%rho(ip3,jm3) + cp3x * cm3y * weight
        fields%rho(ip3,jm2) = fields%rho(ip3,jm2) + cp3x * cm2y * weight
        fields%rho(ip3,jm1) = fields%rho(ip3,jm1) + cp3x * cm1y * weight
        fields%rho(ip3,j  ) = fields%rho(ip3,j  ) + cp3x * cy   * weight
        fields%rho(ip3,jp1) = fields%rho(ip3,jp1) + cp3x * cp1y * weight
        fields%rho(ip3,jp2) = fields%rho(ip3,jp2) + cp3x * cp2y * weight
        fields%rho(ip3,jp3) = fields%rho(ip3,jp3) + cp3x * cp3y * weight

    end do
    
    fields%rho(1:nx,ny+1) = fields%rho(1:nx,1)
    fields%rho(nx+1,1:ny) = fields%rho(1,1:ny)
    fields%rho(nx+1,ny+1) = fields%rho(1,1)
    
    fields%rho = fields%rho / (dx*dy)
    
    rho_total = sum(fields%rho(1:nx,1:ny)) * dx * dy
    print*, " rho total : ", rho_total

    fields%rho = fields%rho - rho_total/dimx/dimy


end subroutine compute_rho_m6_complex

subroutine compute_rho_m6_real( fields, particles) 

    type(fields_2d_t) :: fields
    type(particles_t) :: particles

    fields%rho = 0d0

    nx = fields%mesh%nx
    ny = fields%mesh%ny

    dx = fields%mesh%dx
    dy = fields%mesh%dy

    dimx = fields%mesh%xmax - fields%mesh%xmin
    dimy = fields%mesh%ymax - fields%mesh%ymin

    do k = 1,particles%nbpart
    
        px = particles%x(1,k)/dx
        py = particles%x(2,k)/dy

        px = modulo(px, real(nx,8))
        py = modulo(py, real(ny,8))

        i   = floor(px)
        dpx = px - real(i , kind=8)
        j   = floor(py)
        dpy = py - real(j , kind=8)

        weight = particles%w
      
        im3 = modulo(i -3,nx)+1
        im2 = modulo(i -2,nx)+1
        im1 = modulo(i -1,nx)+1
        ip1 = modulo(i +1,nx)+1
        ip2 = modulo(i +2,nx)+1
        ip3 = modulo(i +3,nx)+1
        jm3 = modulo(j -3,ny)+1
        jm2 = modulo(j -2,ny)+1
        jm1 = modulo(j -1,ny)+1
        jp1 = modulo(j +1,ny)+1
        jp2 = modulo(j +2,ny)+1
        jp3 = modulo(j +3,ny)+1

        i  = i+1
        j  = j+1
      
        cm3x = f_m6(3+dpx)
        cp3x = f_m6(3-dpx)
        cm2x = f_m6(2+dpx)
        cp2x = f_m6(2-dpx)
        cm1x = f_m6(1+dpx)
        cp1x = f_m6(1-dpx)
        cx   = f_m6(dpx)
        cy   = f_m6(dpy)
        cp1y = f_m6(1-dpy)
        cm1y = f_m6(1+dpy)
        cp2y = f_m6(2-dpy)
        cm2y = f_m6(2+dpy)
        cp3y = f_m6(3-dpy)
        cm3y = f_m6(3+dpy)
      
        fields%rho(im3,jm3) = fields%rho(im3,jm3)+cm3x * cm3y * weight
        fields%rho(im3,jm2) = fields%rho(im3,jm2)+cm3x * cm2y * weight
        fields%rho(im3,jm1) = fields%rho(im3,jm1)+cm3x * cm1y * weight
        fields%rho(im3,j  ) = fields%rho(im3,j  )+cm3x * cy   * weight
        fields%rho(im3,jp1) = fields%rho(im3,jp1)+cm3x * cp1y * weight
        fields%rho(im3,jp2) = fields%rho(im3,jp2)+cm3x * cp2y * weight
        fields%rho(im3,jp3) = fields%rho(im3,jp3)+cm3x * cp3y * weight

        fields%rho(im2,jm3) = fields%rho(im2,jm3)+cm2x * cm3y * weight
        fields%rho(im2,jm2) = fields%rho(im2,jm2)+cm2x * cm2y * weight
        fields%rho(im2,jm1) = fields%rho(im2,jm1)+cm2x * cm1y * weight
        fields%rho(im2,j  ) = fields%rho(im2,j  )+cm2x * cy   * weight
        fields%rho(im2,jp1) = fields%rho(im2,jp1)+cm2x * cp1y * weight
        fields%rho(im2,jp2) = fields%rho(im2,jp2)+cm2x * cp2y * weight
        fields%rho(im2,jp3) = fields%rho(im2,jp3)+cm2x * cp3y * weight

        fields%rho(im1,jm3) = fields%rho(im1,jm3)+cm1x * cm3y * weight
        fields%rho(im1,jm2) = fields%rho(im1,jm2)+cm1x * cm2y * weight
        fields%rho(im1,jm1) = fields%rho(im1,jm1)+cm1x * cm1y * weight
        fields%rho(im1,j  ) = fields%rho(im1,j  )+cm1x * cy   * weight
        fields%rho(im1,jp1) = fields%rho(im1,jp1)+cm1x * cp1y * weight
        fields%rho(im1,jp2) = fields%rho(im1,jp2)+cm1x * cp2y * weight
        fields%rho(im1,jp3) = fields%rho(im1,jp3)+cm1x * cp3y * weight

        fields%rho(i  ,jm3) = fields%rho(i  ,jm3)+cx   * cm3y * weight
        fields%rho(i  ,jm2) = fields%rho(i  ,jm2)+cx   * cm2y * weight
        fields%rho(i  ,jm1) = fields%rho(i  ,jm1)+cx   * cm1y * weight
        fields%rho(i  ,j  ) = fields%rho(i  ,j  )+cx   * cy   * weight
        fields%rho(i  ,jp1) = fields%rho(i  ,jp1)+cx   * cp1y * weight
        fields%rho(i  ,jp2) = fields%rho(i  ,jp2)+cx   * cp2y * weight
        fields%rho(i  ,jp3) = fields%rho(i  ,jp3)+cx   * cp3y * weight

        fields%rho(ip1,jm3) = fields%rho(ip1,jm3)+cp1x * cm3y * weight
        fields%rho(ip1,jm2) = fields%rho(ip1,jm2)+cp1x * cm2y * weight
        fields%rho(ip1,jm1) = fields%rho(ip1,jm1)+cp1x * cm1y * weight
        fields%rho(ip1,j  ) = fields%rho(ip1,j  )+cp1x * cy   * weight
        fields%rho(ip1,jp1) = fields%rho(ip1,jp1)+cp1x * cp1y * weight
        fields%rho(ip1,jp2) = fields%rho(ip1,jp2)+cp1x * cp2y * weight
        fields%rho(ip1,jp3) = fields%rho(ip1,jp3)+cp1x * cp3y * weight

        fields%rho(ip2,jm3) = fields%rho(ip2,jm3)+cp2x * cm3y * weight
        fields%rho(ip2,jm2) = fields%rho(ip2,jm2)+cp2x * cm2y * weight
        fields%rho(ip2,jm1) = fields%rho(ip2,jm1)+cp2x * cm1y * weight
        fields%rho(ip2,j  ) = fields%rho(ip2,j  )+cp2x * cy   * weight
        fields%rho(ip2,jp1) = fields%rho(ip2,jp1)+cp2x * cp1y * weight
        fields%rho(ip2,jp2) = fields%rho(ip2,jp2)+cp2x * cp2y * weight
        fields%rho(ip2,jp3) = fields%rho(ip2,jp3)+cp2x * cp3y * weight

        fields%rho(ip3,jm3) = fields%rho(ip3,jm3)+cp3x * cm3y * weight
        fields%rho(ip3,jm2) = fields%rho(ip3,jm2)+cp3x * cm2y * weight
        fields%rho(ip3,jm1) = fields%rho(ip3,jm1)+cp3x * cm1y * weight
        fields%rho(ip3,j  ) = fields%rho(ip3,j  )+cp3x * cy   * weight
        fields%rho(ip3,jp1) = fields%rho(ip3,jp1)+cp3x * cp1y * weight
        fields%rho(ip3,jp2) = fields%rho(ip3,jp2)+cp3x * cp2y * weight
        fields%rho(ip3,jp3) = fields%rho(ip3,jp3)+cp3x * cp3y * weight

    end do
    
    fields%rho(1:nx,ny+1) = fields%rho(1:nx,1)
    fields%rho(nx+1,1:ny) = fields%rho(1,1:ny)
    fields%rho(nx+1,ny+1) = fields%rho(1,1)
    
    fields%rho = fields%rho / (dx*dy)
    
    rho_total = sum(fields%rho(1:nx,1:ny)) * dx * dy

    fields%rho = fields%rho - rho_total/dimx/dimy

end subroutine compute_rho_m6_real

end module m6_compute_rho_m
