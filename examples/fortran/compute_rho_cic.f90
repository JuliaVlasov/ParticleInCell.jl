module compute_rho_cic_m

use mesh_fields_m
use particles_m

implicit none

contains


subroutine compute_rho_cic( f, p)

type(particles_t), intent(in)    :: p
type(fields_3d_t), intent(inout) :: f

real(8)    :: a1, a2, a3, a4, a5, a6, a7, a8
real(8)    :: xp, yp, zp, vol
real(8)    :: dxp, dyp, dzp
integer    :: ip, jp, kp
real(8)    :: dx, dy, dz
integer    :: nx, ny, nz
integer(8) :: m

dx = f%mesh%dx
dy = f%mesh%dy
dz = f%mesh%dz

f%rho = 0.0d0

vol   = p%w / (dx*dy*dz)

do m = 1_8, p%nbpart

   xp = p%x(1, m)/dx
   yp = p%x(2, m)/dy
   zp = p%x(3, m)/dz

   ip = floor(xp)
   jp = floor(yp)
   kp = floor(zp)

   dxp = xp - real(ip, kind=8)
   dyp = yp - real(jp, kind=8)
   dzp = zp - real(kp, kind=8)

   a1 = (1d0 - dxp) * (1d0 - dyp) * (1d0 - dzp) 
   a2 = dxp         * (1d0 - dyp) * (1d0 - dzp) 
   a3 = (1d0 - dxp) * dyp         * (1d0 - dzp) 
   a4 = dxp         * dyp         * (1d0 - dzp)
   a5 = (1d0 - dxp) * (1d0 - dyp) * dzp
   a6 = dxp         * (1d0 - dyp) * dzp
   a7 = (1d0 - dxp) * dyp         * dzp
   a8 = dxp         * dyp         * dzp

   ip = ip+1
   jp = jp+1
   kp = kp+1

   f%rho(ip  ,jp  ,kp  ) = f%rho(ip  ,jp  ,kp  ) + a1 * vol
   f%rho(ip+1,jp  ,kp  ) = f%rho(ip+1,jp  ,kp  ) + a2 * vol
   f%rho(ip  ,jp+1,kp  ) = f%rho(ip  ,jp+1,kp  ) + a3 * vol
   f%rho(ip+1,jp+1,kp  ) = f%rho(ip+1,jp+1,kp  ) + a4 * vol
   f%rho(ip  ,jp  ,kp+1) = f%rho(ip  ,jp  ,kp+1) + a5 * vol
   f%rho(ip+1,jp  ,kp+1) = f%rho(ip+1,jp  ,kp+1) + a6 * vol
   f%rho(ip  ,jp+1,kp+1) = f%rho(ip  ,jp+1,kp+1) + a7 * vol
   f%rho(ip+1,jp+1,kp+1) = f%rho(ip+1,jp+1,kp+1) + a8 * vol

end do

nx = f%mesh%nx
ny = f%mesh%ny
nz = f%mesh%nz

f%rho(nx+1,:,:) = f%rho(1,:,:)
f%rho(:,ny+1,:) = f%rho(:,1,:)
f%rho(:,:,nz+1) = f%rho(:,:,1)


end subroutine compute_rho_cic


end module compute_rho_cic_m
