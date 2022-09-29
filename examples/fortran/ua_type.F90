module ua_type_m

use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f03'

type :: ua_t
  
    integer                                :: ntau 
    real(8)                                :: eps 
    real(8), allocatable                   :: tau(:)
    real(8), allocatable                   :: ltau(:)

    complex(8), allocatable                :: pl(:,:)
    complex(8), allocatable                :: ql(:,:)

    type(C_PTR)                            :: fw
    type(C_PTR)                            :: bw
    complex(C_DOUBLE_COMPLEX), allocatable :: ft(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: rt(:,:)
    complex(C_DOUBLE_COMPLEX), allocatable :: rf(:,:)

end type ua_t

real(8)               :: pi
integer, private      :: n

contains

subroutine init_ua( ua, ntau, eps, nbpart )

    type(ua_t)             :: ua
    integer,    intent(in) :: ntau
    real(8),    intent(in) :: eps
    integer(8), intent(in) :: nbpart

    real(8) :: dtau


    pi = 4d0 * atan(1d0)

    ua%ntau = ntau
    ua%eps  = cmplx(eps,0d0,kind=8)

    dtau = 2d0 * pi / real(ntau,kind=8)
    
    allocate(ua%ltau(ntau))

    do n = 1, ntau/2
        ua%ltau(n) = cmplx(n-1,0d0,kind=8)
    end do
    do n = ntau/2+1, ntau
        ua%ltau(n) = cmplx(n-1-ntau,0d0,kind=8)
    end do
    
    allocate(ua%tau(ntau))

    do n = 1, ntau
        ua%tau(n)  = real(n-1,8) * dtau
    end do

    allocate(ua%pl(ntau, nbpart))
    allocate(ua%ql(ntau, nbpart))

    allocate(ua%ft(ntau))
    allocate(ua%rt(ntau,2))
    allocate(ua%rf(ntau,2))

    ua%fw = fftw_plan_dft_1d(ntau, ua%rt(:,1), ua%rf(:,1), &
    &         FFTW_FORWARD, FFTW_ESTIMATE)
    ua%bw = fftw_plan_dft_1d(ntau, ua%rf(:,1), ua%rt(:,1), &
    &         FFTW_BACKWARD,FFTW_ESTIMATE)

end subroutine init_ua

end module ua_type_m
