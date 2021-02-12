module heat
    use iso_c_binding, only: c_int32_t, c_double
    public :: kernel
contains
    subroutine kernel(m, n, u_in,  u_out, error) bind( C, name="heatKernel" )

        implicit none
        integer(c_int32_t), intent(in) :: m, n
        real(c_double), dimension( 1:m, 1:n ), intent(in) :: u_in
        real(c_double), dimension( 1:m, 1:n ), intent(out) :: u_out
        real(c_double), intent(out) :: error

        integer(c_int32_t) :: i, j

        error = 0.d0
        u_out(2:m-1,2:n-1) = 4.d0 * u_in(2:m-1, 2:n-1) &
                                  - u_in(1:m-2, 2:n-1) &
                                  - u_in(3:m, 2:n-1)   &
                                  - u_in(2:m-1,1:n-2)  &
                                  - u_in(2:m-1,3:n)
        error  =  sum((u_out(:,:) - u_in(:,:))**2)

    end subroutine kernel
end module heat
