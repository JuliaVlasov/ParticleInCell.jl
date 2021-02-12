module Maxwell

use zone
implicit none

integer, private :: i, j

real(kind=prec), private :: dex_dy, dey_dx
real(kind=prec), private :: dbz_dx, dbz_dy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( ex, ey, bz, dt )

   real(8), intent(in) :: ex(:,:), ey(:,:)
   real(8), intent(in) :: dt
   real(8), intent(inout) :: bz(:,:)

   !*** On utilise l'equation de Faraday sur un demi pas
   !*** de temps pour le calcul du champ magnetique  Bz 
   !*** a l'instant n puis n+1/2 apres deplacement des
   !*** particules

   do j=1,ny
   do i=1,nx
      dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
      dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
      bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
   end do
   end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere( ex, ey, bz, jx, jy, dt )

   real(8), intent(inout) :: ex(:,:), ey(:,:)
   real(8), intent(in) :: dt
   real(8), intent(in) :: bz(:,:), jx(:,:), jy(:,:)
   !*** Calcul du champ electrique E au temps n+1
   !*** sur les points internes du maillage
   !*** Ex aux points (i+1/2,j)
   !*** Ey aux points (i,j+1/2)

   do i=1,nx
   do j=2,ny
      dbz_dy = (bz(i,j)-bz(i,j-1)) / dy
      ex(i,j) = ex(i,j) + dt * dbz_dy - dt * jx(i,j)
   end do
   end do

   do i=2,nx
   do j=1,ny
      dbz_dx = (bz(i,j)-bz(i-1,j)) / dx
      ey(i,j) = ey(i,j) - dt * dbz_dx - dt * jy(i,j)
   end do
   end do

   do i=1,nx
      dbz_dy = (bz(i,1)-bz(i,ny)) / dy
      ex(i,1)  = ex(i,1) + dt * dbz_dy - dt * jx(i,1)
      ex(i,ny+1) = ex(i,1)
   end do
   
   do j=1,ny
      dbz_dx = (bz(1,j)-bz(nx,j)) / dx
      ey(1,j)  = ey(1,j) - dt * dbz_dx - dt * jy(1,j)
      ey(nx+1,j) = ey(1,j)
   end do

end subroutine ampere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decalage( tm, tm1 )

type(mesh_fields) :: tm, tm1

!*** Calcul des composantes des champs 
!*** sur les noeuds du maillage de rho
!*** par interpolation lineaire

do i=1,nx-1
   do j=1,ny-1
      tm1%ex(i,j) = 0.5*( tm%ex(i-1,j) + tm%ex(i,j) )
      tm1%ey(i,j) = 0.5*( tm%ey(i,j-1) + tm%ey(i,j) )
      tm1%bz(i,j) = 0.25 * ( tm%bz(i-1,j-1) + tm%bz(i,j-1) + tm%bz(i-1,j) + tm%bz(i,j) )
   end do
end do

do i=1,nx-1 !Sud et Nord
   tm1%ex(i,0) = 0.5 * ( tm%ex(i-1,0) + tm%ex(i,0) ) 
   tm1%ey(i,0) = 0.5 * ( tm%ey(i,ny-1) + tm%ey(i,0) ) 
   tm1%bz(i,0) = 0.25 * ( tm%bz(i-1,ny-1) + tm%bz(i,ny-1) + tm%bz(i-1,0) + tm%bz(i,0) )
   tm1%ex(i,ny) = tm1%ex(i,0) 
   tm1%ey(i,ny) = tm1%ey(i,0) 
   tm1%bz(i,ny) = tm1%bz(i,0) 
end do

do j=1,ny-1 !Ouest et Est
   tm1%ex(0,j) = 0.5 * ( tm%ex(nx-1,j) + tm%ex(0,j) )
   tm1%ey(0,j) = 0.5 * ( tm%ey(0,j-1) + tm%ey(0,j) )
   tm1%bz(0,j) = 0.25 * ( tm%bz(nx-1,j-1) + tm%bz(0,j-1) + tm%bz(nx-1,j) + tm%bz(0,j) ) 
   tm1%ex(nx,j) = tm1%ex(0,j) 
   tm1%ey(nx,j) = tm1%ey(0,j) 
   tm1%bz(nx,j) = tm1%bz(0,j) 
end do

!Coins
tm1%ex(0,0) = 0.5 * ( tm%ex(nx-1,0) + tm%ex(0,0) )
tm1%ey(0,0) = 0.5 * ( tm%ey(0,ny-1) + tm%ey(0,0) )
tm1%bz(0,0) = 0.25 * ( tm%bz(nx-1,ny-1) + tm%bz(0,ny-1) + tm%bz(nx-1,0) + tm%bz(0,0) ) 

tm1%ex(nx,0)  = tm1%ex(0,0) 
tm1%ex(nx,ny) = tm1%ex(0,0) 
tm1%ex(0,ny)  = tm1%ex(0,0) 

tm1%ey(nx,0)  = tm1%ey(0,0) 
tm1%ey(nx,ny) = tm1%ey(0,0) 
tm1%ey(0,ny)  = tm1%ey(0,0) 

tm1%bz(nx,0)  = tm1%bz(0,0) 
tm1%bz(nx,ny) = tm1%bz(0,0) 
tm1%bz(0,ny)  = tm1%bz(0,0) 

end subroutine decalage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Maxwell
