module ua_steps_m

use ua_type_m
use mesh_fields_m
use particles_m
use m6_compute_rho_m
use m6_interpolation_m

implicit none

integer, private :: m, n

contains

subroutine preparation( ua, dt, particles, xt, yt) 

    type(ua_t)         , intent(inout) :: ua 
    real(8)            , intent(in)    :: dt
    type(particles_t)  , intent(inout) :: particles 
    complex(8)         , intent(out)   :: xt(:,:,:)         
    complex(8)         , intent(out)   :: yt(:,:,:)         

    real(8)    :: x1
    real(8)    :: x2
    complex(8) :: elt
    real(8)    :: t
    real(8)    :: b
    real(8)    :: h1
    real(8)    :: h2
    real(8)    :: ex
    real(8)    :: ey
    real(8)    :: vx
    real(8)    :: vy
    real(8)    :: vxb
    real(8)    :: vyb
    real(8)    :: xt1
    real(8)    :: xt2
    integer    :: ntau
    integer    :: nbpart
    real(8)    :: eps
    real(8)    :: interv
    real(8)    :: exb
    real(8)    :: eyb

    eps    = ua%eps
    ntau   = ua%ntau
    nbpart = particles%nbpart

    do m = 1,nbpart

        x1 = particles%x(1,m)
        x2 = particles%x(2,m)

        b = 1d0 + 0.5d0 * sin(x1) * sin(x2)
        t = dt * b

        particles%b(m) = b
        particles%t(m) = t

        ua%pl(1,m) = cmplx(t, 0d0, kind=8)
        ua%ql(1,m) = cmplx(t**2 / 2d0, 0d0, kind=8)

        do n=2,ua%ntau
            elt = exp((0d0,-1d0) * ua%ltau(n)*t/eps)
            ua%pl(n,m) = eps * cmplx(0d0,1d0,kind=8) * (elt-1d0)/ua%ltau(n)
            ua%ql(n,m) = eps * (eps*(1d0-elt) -cmplx(0d0,1d0,kind=8)*ua%ltau(n)*t)/ua%ltau(n)**2
        end do

        ex  = particles%e(1,m)
        ey  = particles%e(2,m)
        vx  = particles%v(1,m)
        vy  = particles%v(2,m)
        vxb = vx/b
        vyb = vy/b

        do n = 1,ntau

            h1 = eps * (sin(ua%tau(n)) * vxb - cos(ua%tau(n)) * vyb)
            h2 = eps * (sin(ua%tau(n)) * vyb + cos(ua%tau(n)) * vxb)

            xt1 = x1 + h1 + eps * vyb
            xt2 = x2 + h2 - eps * vxb

            xt(n,1,m) = cmplx(xt1,0d0,kind=8)
            xt(n,2,m) = cmplx(xt2,0d0,kind=8)

            interv=(1d0+0.5d0*sin(xt1)*sin(xt2)-b)/eps

            exb = ((  cos(ua%tau(n))*vy - sin(ua%tau(n))*vx) * interv + ex)/b
            eyb = (( -cos(ua%tau(n))*vx - sin(ua%tau(n))*vy) * interv + ey)/b

            ua%rt(n,1) = cmplx(cos(ua%tau(n))*exb-sin(ua%tau(n))*eyb,0d0,kind=8)
            ua%rt(n,2) = cmplx(sin(ua%tau(n))*exb+cos(ua%tau(n))*eyb,0d0,kind=8)

        end do

        call fftw_execute_dft(ua%fw, ua%rt(:,1), ua%rf(:,1)) 
        call fftw_execute_dft(ua%fw, ua%rt(:,2), ua%rf(:,2)) 

        do n = 2,ntau
            ua%rf(n,1) = -(0d0, 1d0) / ua%ltau(n) * ua%rf(n,1) / real(ntau, kind=8)
            ua%rf(n,2) = -(0d0, 1d0) / ua%ltau(n) * ua%rf(n,2) / real(ntau, kind=8)
        end do

        call fftw_execute_dft(ua%bw, ua%rf(:,1), ua%rt(:,1)) 
        call fftw_execute_dft(ua%bw, ua%rf(:,2), ua%rt(:,2)) 

        do n = 1,ntau
            yt(n,1,m) = vx + (ua%rt(n,1) - ua%rt(1,1)) * eps
            yt(n,2,m) = vy + (ua%rt(n,2) - ua%rt(1,2)) * eps
        end do

    end do

end subroutine preparation

subroutine interpolation( particles, et, fields, ua, xt ) 

    type(particles_t), intent(inout)  :: particles
    real(8),           intent(out)    :: et(:,:,:)
    type(fields_2d_t), intent(in)     :: fields
    type(ua_t),        intent(in)     :: ua
    complex(8),        intent(in)     :: xt(:,:,:)

    call interpolate_eb_m6_complex( et, fields, xt, particles%nbpart, ua%ntau) 

end subroutine interpolation

subroutine deposition( particles, fields, ua, xt) 

    type(particles_t), intent(inout) :: particles
    type(fields_2d_t), intent(inout) :: fields
    type(ua_t),        intent(inout) :: ua
    complex(8),        intent(in)    :: xt(:,:,:)

    call compute_rho_m6_complex( fields, particles, xt, ua )

end subroutine deposition

subroutine compute_f( fx, fy, ua, particles, xt, yt, et )

    complex(8)        :: fx(:,:,:)
    complex(8)        :: fy(:,:,:)
    type(ua_t)        :: ua 
    type(particles_t) :: particles
    complex(8)        :: xt(:,:,:)
    complex(8)        :: yt(:,:,:)
    real(8)           :: et(:,:,:)

    real(8)    :: interv
    real(8)    :: tau
    complex(8) :: tmp1
    complex(8) :: tmp2
    real(8)    :: xt1
    real(8)    :: xt2
    complex(8) :: yt1
    complex(8) :: yt2
    real(8)    :: b

    do m=1,particles%nbpart
    
        b = particles%b(m)

        do n=1,ua%ntau
    
            xt1 = real(xt(n,1,m))
            xt2 = real(xt(n,2,m))
    
            yt1 = yt(n,1,m) 
            yt2 = yt(n,2,m) 
    
            tau = ua%tau(n)
    
            fx(n,1,m) = ( cos(tau) * yt1 + sin(tau) * yt2)/b
            fx(n,2,m) = (-sin(tau) * yt1 + cos(tau) * yt2)/b
    
            interv = (1d0 + 0.5d0*sin(xt1)*sin(xt2)-b)/ua%eps
    
            tmp1 = et(n,1,m)+(  cos(tau)*yt2 - sin(tau)*yt1)*interv
            tmp2 = et(n,2,m)+(- cos(tau)*yt1 - sin(tau)*yt2)*interv
    
            fy(n,1,m) = (cos(tau)*tmp1-sin(tau)*tmp2)/b
            fy(n,2,m) = (sin(tau)*tmp1+cos(tau)*tmp2)/b
    
        end do

        call fftw_execute_dft(ua%fw, fx(:,1,m), fx(:,1,m)) 
        call fftw_execute_dft(ua%fw, fx(:,2,m), fx(:,2,m)) 
        call fftw_execute_dft(ua%fw, fy(:,1,m), fy(:,1,m)) 
        call fftw_execute_dft(ua%fw, fy(:,2,m), fy(:,2,m)) 

    end do

    fx = fx / real(ua%ntau, kind=8)
    fy = fy / real(ua%ntau, kind=8)


end subroutine compute_f

subroutine ua_step1( xt, xf, ua, particles, fx )

    complex(8)        , intent(inout)  :: xt(:,:,:)
    complex(8)        , intent(inout)  :: xf(:,:,:)
    type(ua_t)        , intent(inout)  :: ua
    type(particles_t) , intent(in)     :: particles 
    complex(8)        , intent(in)     :: fx(:,:,:)

    real(8)           :: t
    integer           :: m
    integer           :: n
    complex(8)        :: elt

    complex(8), allocatable :: tmp(:,:)

    do m = 1, particles%nbpart

        call dfftw_execute_dft( ua%fw, xt(:,1,m), xf(:,1,m))
        call dfftw_execute_dft( ua%fw, xt(:,2,m), xf(:,2,m))

        t = particles%t(m)

        do n=1,ua%ntau

            elt = exp(-cmplx(0d0,1d0,kind=8)*ua%ltau(n)*t/ua%eps)
            elt = elt / real(ua%ntau, kind=8)
            ua%rf(n,1) = elt * xf(n,1,m) + ua%pl(n,m) * fx(n,1,m)
            ua%rf(n,2) = elt * xf(n,2,m) + ua%pl(n,m) * fx(n,2,m)

        end do

        call fftw_execute_dft(ua%bw, ua%rf(:,1), xt(:,1,m)) 
        call fftw_execute_dft(ua%bw, ua%rf(:,2), xt(:,2,m)) 

    end do

end subroutine ua_step1

subroutine ua_step2( xt, xf, ua, particles, fx, gx )

    complex(8)        , intent(inout)  :: xt(:,:,:)
    complex(8)        , intent(in)     :: xf(:,:,:)
    type(ua_t)        , intent(inout)  :: ua
    type(particles_t) , intent(in)     :: particles 
    complex(8)        , intent(in)     :: fx(:,:,:)
    complex(8)        , intent(in)     :: gx(:,:,:)

    real(8)           :: t
    integer           :: m
    integer           :: n
    complex(8)        :: elt

    do m = 1, particles%nbpart

        t = particles%t(m)

        do n=1,ua%ntau

            elt = exp(-cmplx(0d0,1d0,kind=8)*ua%ltau(n)*t/ua%eps)
            elt = elt / real(ua%ntau, kind=8)
            ua%rf(n,1) = elt * xf(n,1,m) + ua%pl(n,m) * fx(n,1,m) &
                     + ua%ql(n,m) * (gx(n,1,m) - fx(n,1,m)) / t
            ua%rf(n,2) = elt * xf(n,2,m) + ua%pl(n,m) * fx(n,2,m) &
                     + ua%ql(n,m) * (gx(n,2,m) - fx(n,2,m)) / t

        end do

        call fftw_execute_dft(ua%bw, ua%rf(:,1), xt(:,1,m)) 
        call fftw_execute_dft(ua%bw, ua%rf(:,2), xt(:,2,m)) 

    end do

end subroutine ua_step2

subroutine compute_v( ua, particles, yt, yf )

    type(ua_t)        , intent(inout)  :: ua
    type(particles_t) , intent(inout)  :: particles 
    complex(8)        , intent(inout)  :: yt(:,:,:)
    complex(8)        , intent(inout)  :: yf(:,:,:)

    real(8)    :: t
    complex(8) :: elt
    complex(8) :: px
    complex(8) :: py

    do m=1, particles%nbpart

        t = particles%t(m)

        call fftw_execute_dft(ua%fw, yt(:,1,m), yf(:,1,m)) 
        call fftw_execute_dft(ua%fw, yt(:,2,m), yf(:,2,m)) 

        px = (0d0, 0d0)
        py = (0d0, 0d0)

        do n = 1, ua%ntau
            elt = exp(cmplx(0d0,1d0,kind=8)*ua%ltau(n)*t/ua%eps) 
            px = px + yf(n,1,m)/real(ua%ntau,kind=8) * elt
            py = py + yf(n,2,m)/real(ua%ntau,kind=8) * elt
        end do

        particles%v(1,m) = real(cos(t/ua%eps)*px+sin(t/ua%eps)*py)
        particles%v(2,m) = real(cos(t/ua%eps)*py-sin(t/ua%eps)*px)

    end do

end subroutine compute_v

end module ua_steps_m
