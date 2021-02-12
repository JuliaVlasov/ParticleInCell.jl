module Maxwell

use zone
implicit none

integer, private :: i, j

real(kind=prec), private :: dex_dy, dey_dx
real(kind=prec), private :: dbz_dx, dbz_dy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( tm )

type( mesh_fields ) :: tm

   !*** On utilise l'equation de Faraday sur un demi pas
   !*** de temps pour le calcul du champ magnetique  Bz 
   !*** a l'instant n puis n+1/2 apres deplacement des
   !*** particules

   do i=0,nx-1
   do j=0,ny-1
      dex_dy     = (tm%ex(i,j+1)-tm%ex(i,j)) / dy
      dey_dx     = (tm%ey(i+1,j)-tm%ey(i,j)) / dx
      tm%bz(i,j) = tm%bz(i,j) + 0.5 * dt * (dex_dy - dey_dx)
   end do
   end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere( tm )

type( mesh_fields ) :: tm

   !*** Calcul du champ electrique E au temps n+1
   !*** sur les points internes du maillage
   !*** Ex aux points (i+1/2,j)
   !*** Ey aux points (i,j+1/2)

   do i=0,nx-1
   do j=1,ny-1
      dbz_dy = (tm%bz(i,j)-tm%bz(i,j-1)) / dy
      tm%ex(i,j) = tm%ex(i,j) + dt * dbz_dy - dt * tm%jx(i,j)
   end do
   end do

   do i=1,nx-1
   do j=0,ny-1
      dbz_dx = (tm%bz(i,j)-tm%bz(i-1,j)) / dx
      tm%ey(i,j) = tm%ey(i,j) - dt * dbz_dx - dt * tm%jy(i,j)
   end do
   end do

   do i=0,nx-1
      dbz_dy = (tm%bz(i,0)-tm%bz(i,ny-1)) / dy
      tm%ex(i,0)  = tm%ex(i,0) + dt * dbz_dy - dt * tm%jx(i,0)
      tm%ex(i,ny) = tm%ex(i,0)
   end do
   
   do j=0,ny-1
      dbz_dx = (tm%bz(0,j)-tm%bz(nx-1,j)) / dx
      tm%ey(0,j)  = tm%ey(0,j) - dt * dbz_dx - dt * tm%jy(0,j)
      tm%ey(nx,j) = tm%ey(0,j)
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
