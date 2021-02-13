module quietstart

use iso_c_binding

implicit none

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure real(c_double) function bit_reversing( n )

integer(c_int32_t), intent(in) :: n
integer(c_int32_t) :: deci, k, div
real(c_double) :: miroir

k = 25
miroir = 0.d0
deci = n

if (deci > 33554432) then 
   !print*,'too big > 33554432=2**25'
   stop
else 
   do while (k >= 0)
      div = 2**k
      if (deci/div == 1) then
        ! a = a + 10**k
         miroir = miroir + 2.**(-k-1)
         deci = deci - div
      endif
      k = k-1
   enddo
endif

bit_reversing = miroir

end function bit_reversing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure real(c_double) function trinary_reversing( n )

integer(c_int32_t), intent(in) :: n
integer(c_int32_t) :: deci, k, div
real(c_double) :: miroir

!a = 0
k = 16
miroir = 0.d0
deci = n

if (deci > 43046721) then 
   !print*, ' too big > 43046721=3**16 '
   stop
else 
   do while (k >= 0)
      div = 3**k
      if (deci/div == 1) then
        ! a = a + 10**k
         miroir = miroir + 3.**(-k-1)
         deci = deci - div
      else if (deci/div == 2) then
        ! a = a + 2 * 10**k
         miroir = miroir + 2 * 3.**(-k-1)
         deci = deci - 2 * div
      endif
      k = k-1
   enddo
endif

trinary_reversing = miroir

end function trinary_reversing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure real(c_double) function penta_reversing( n )

integer, intent(in) :: n
integer :: deci, k, div
real(c_double) :: miroir

!a = 0
k = 11
miroir = 0.d0
deci = n

if (deci > 48828125) then 
   !print*,'too big > 48828125=5**11'
   stop
else 
   do while (k >= 0)
      div = 5**k
      if (deci/div == 1) then
        ! a = a + 10**k
         miroir = miroir + 5.**(-k-1)
         deci = deci - div
      else if (deci/div == 2) then
        ! a = a + 2 * 10**k
         miroir = miroir + 2 * 5.**(-k-1)
         deci = deci - 2 * div
      else if (deci/div == 3) then
        ! a = a + 3 * 10**k
         miroir = miroir + 3 * 5.**(-k-1)
         deci = deci - 3 * div
      else if (deci/div == 4) then
        ! a = a + 4 * 10**k
         miroir = miroir + 4 * 5.**(-k-1)
         deci = deci - 4 * div
      endif
      k = k-1
   enddo
endif

penta_reversing = miroir

end function penta_reversing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dichotomie_x( alpha, kx, a, b, R, eps ) 

! il faut D(a)<R<D(b), on cherche x tq R=D(x), resu dans a 

real(c_double), intent(in) :: alpha, kx, R, eps
real(c_double) :: a, b
real(c_double) :: x, D, pi

pi = 4d0 * atan(1d0)

D = ( kx*a + alpha * sin(kx*a) ) / (2*pi)
do while ( D<R-eps .or. D>R+eps )
   x = (a+b)/2
   D = ( kx*x + alpha * sin(kx*x) ) / (2*pi)
   if ( D<R-eps ) then
      a = x
   else if ( D>R+eps ) then 
      b = x
   else
      a = x
   endif
end do

end subroutine dichotomie_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module quietstart
