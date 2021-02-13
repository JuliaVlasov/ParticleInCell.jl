module Maxwell

use iso_c_binding
use zone

implicit none

integer, private :: i, j

real(c_double), private :: dex_dy, dey_dx
real(c_double), private :: dbz_dx, dbz_dy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( f, dt )

   real(8), intent(inout) :: f(:,:,:)
   real(8), intent(in) :: dt

   do j=1,ny
   do i=1,nx
      dex_dy  = (f(1,i,mod1(j+1,ny))-f(1,i,j)) / dy
      dey_dx  = (f(2,mod1(i+1,nx),j)-f(2,i,j)) / dx
      f(3,i,j) = f(3,i,j) + dt * (dex_dy - dey_dx)
   end do
   end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere( f, j0, dt )

   real(c_double), intent(inout) :: f(:,:,:), j0(:,:,:)
   real(c_double), intent(in) :: dt

   do i=1,nx
   do j=1,ny
      dbz_dy = (f(3,i,j)-f(3,i,mod1(j-1,ny))) / dy
      f(1,i,j) = f(1,i,j) + dt * dbz_dy - dt * j0(1,i,j)
   end do
   end do

   do i=1,nx
   do j=1,ny
      dbz_dx = (f(3,i,j)-f(3,mod1(i-1,nx),j)) / dx
      f(2,i,j) = f(2,i,j) - dt * dbz_dx - dt * j0(2,i,j)
   end do
   end do

end subroutine ampere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decalage( f0, f1 )

real(c_double), intent(in) :: f0(:,:,:)
real(c_double), intent(out) :: f1(:,:,:)

do i=1,nx
   do j=1,ny
      f1(1,i,j) = 0.50 * ( f0(1,mod1(i-1,nx),j) + f0(1,i,j) )
      f1(2,i,j) = 0.50 * ( f0(2,i,mod1(j-1,ny)) + f0(2,i,j) )
      f1(3,i,j) = 0.25 * ( f0(3,mod1(i-1,nx),mod1(j-1,ny)) &
                         + f0(3,i,mod1(j-1,ny)) &
                         + f0(3,mod1(i-1,nx),j) &
                         + f0(3,i,j) )
   end do
end do

end subroutine decalage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Maxwell
