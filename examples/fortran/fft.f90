module fft_m

use, intrinsic :: iso_c_binding

implicit none
include "fftw3.f03" 

integer     :: nfft
type(C_PTR) :: plan_fw
type(C_PTR) :: plan_bw

interface fft
    module procedure fft_1d
    module procedure fft_2d
end interface fft

interface ifft
    module procedure ifft_1d
    module procedure ifft_2d
end interface ifft

contains

subroutine init_fft( n, source, destination )

    integer,    intent(in)     :: n
    complex(8), intent(inout)  :: source(:)
    complex(8), intent(inout)  :: destination(:)

    nfft = n

    plan_fw = fftw_plan_dft_1d(n, source, destination, FFTW_FORWARD, FFTW_ESTIMATE)
    plan_bw = fftw_plan_dft_1d(n, destination, source, FFTW_BACKWARD,FFTW_ESTIMATE)

end subroutine init_fft

subroutine fft_2d( source, destination)

    complex(8), intent(inout)  :: source(:,:)
    complex(8), intent(out)    :: destination(:,:)

    call fftw_execute_dft(plan_fw, source(:,1), destination(:,1))
    call fftw_execute_dft(plan_fw, source(:,2), destination(:,2))

    destination = destination / nfft

end subroutine fft_2d

subroutine ifft_2d( source, destination)

    complex(8), intent(inout) :: source(:,:)
    complex(8), intent(out)   :: destination(:,:)

    call fftw_execute_dft(plan_bw, source(:,1), destination(:,1))
    call fftw_execute_dft(plan_bw, source(:,2), destination(:,2))

end subroutine ifft_2d


subroutine fft_1d( source, destination)

    complex(8), intent(inout)  :: source(:)
    complex(8), intent(out)    :: destination(:)

    call fftw_execute_dft(plan_fw, source, destination)

    destination = destination / nfft

end subroutine fft_1d

subroutine ifft_1d( source, destination)

    complex(8), intent(inout) :: source(:)
    complex(8), intent(out)   :: destination(:)

    call fftw_execute_dft(plan_bw, source, destination)

end subroutine ifft_1d


end module fft_m
