module interpolation_cic_m

use mesh_fields_m
use particles_m

implicit none

contains

subroutine interpolate_eb_cic( p, f )

type(particles_t), intent(inout) :: p
type(fields_3d_t), intent(in)    :: f

real(8) :: a1, a2, a3, a4, a5, a6, a7, a8
real(8) :: xp, yp, zp
integer :: i, j, k
real(8) :: dx,  dy,  dz
real(8) :: dxp, dyp, dzp
integer :: d

integer(8) :: m

dx = f%mesh%dx
dy = f%mesh%dy
dz = f%mesh%dz

do m = 1_8, p%nbpart

   xp = p%x(1,m)/dx
   yp = p%x(2,m)/dy
   zp = p%x(3,m)/dz

   i = floor(xp)
   j = floor(yp)
   k = floor(zp)

   dxp = xp - real(i, kind=8)
   dyp = yp - real(j, kind=8)
   dzp = zp - real(k, kind=8)

   a1 = (1d0 - dxp) * (1d0 - dyp) * (1d0 - dzp) 
   a2 = dxp         * (1d0 - dyp) * (1d0 - dzp) 
   a3 = (1d0 - dxp) * dyp         * (1d0 - dzp) 
   a4 = dxp         * dyp         * (1d0 - dzp)
   a5 = (1d0 - dxp) * (1d0 - dyp) * dzp
   a6 = dxp         * (1d0 - dyp) * dzp
   a7 = (1d0 - dxp) * dyp         * dzp
   a8 = dxp         * dyp         * dzp

   i = i + 1
   j = j + 1
   k = k + 1

   do d = 1, 3

       p%e(d,m) = a1 * f%e(d, i,  j,   k  ) + a2 * f%e(d, i+1, j,   k  ) &
                + a3 * f%e(d, i,  j+1, k  ) + a4 * f%e(d, i+1, j+1, k  ) &
                + a5 * f%e(d, i,  j,   k+1) + a6 * f%e(d, i+1, j,   k+1) &
                + a7 * f%e(d, i,  j+1, k+1) + a8 * f%e(d, i+1, j+1, k+1) 

   end do

end do

end subroutine interpolate_eb_cic

end module interpolation_cic_m
