module particles_m

use mesh_fields_m

implicit none

type :: particles_t

   integer(8)           :: nbpart
   real(8), allocatable :: x(:,:)
   real(8), allocatable :: v(:,:)
   real(8), allocatable :: e(:,:)
   real(8), allocatable :: b(:)
   real(8), allocatable :: t(:)
   real(8)              :: w

end type particles_t

contains

subroutine init_particles_2d( self, nbpart, mesh, alpha, kx )

    type(particles_t), intent(out)  :: self
    integer(8),        intent(in)   :: nbpart
    type(mesh_t),      intent(in)   :: mesh
    real(8),           intent(in)   :: alpha
    real(8),           intent(in)   :: kx

    real(8)              :: dimx
    real(8)              :: dimy
    integer              :: m
    integer              :: i
    integer              :: j
    integer              :: k
    real(8), parameter   :: eps = 1.d-12
    real(8)              :: xi, yi, zi
    real(8)              :: temm
    integer              :: nseed
    integer, allocatable :: seed(:)

    self%nbpart = nbpart

    allocate(self%x(2,nbpart))
    allocate(self%v(2,nbpart))
    allocate(self%e(2,nbpart))
    allocate(self%b(nbpart))
    allocate(self%t(nbpart))

    dimx = mesh%xmax - mesh%xmin
    dimy = mesh%ymax - mesh%ymin

    self%w = dimx * dimy / nbpart

    nseed = 33
    allocate(seed(nseed))
    
    seed = [-1584649339, -1457681104, 1579121008,   -819547200, &
              249798090, -517237887,   177452147,   -981503238, &
             1418301473,  1989625004, 2065424384,   -296364178, &
             1658790794, -435188152, -1643185032,   1461389312, &
             1869073641,  1321930686,  483734018,   1269936416, &
            -1999561453,  906251506,   782514880,    428753705, &
            -2031262823,  263953581,  1026600222,  -1118515860, &
             1633712916, -464192498, -1860714528,   1436611533, 0]
    
    call random_seed(put  = seed)

    k = 1
    do while (k<=nbpart)

        call random_number(xi)
        xi   = xi * dimx
        call random_number(yi)
        yi   = yi * dimy
        call random_number(zi)
        zi   = (2d0+alpha)*zi
        temm = 1d0+sin(yi)+alpha*cos(kx*xi)
        if (temm>=zi) then
            self%x(1,k) = xi
            self%x(2,k) = yi
            k = k + 1
        end if
    end do
    
    k = 1
    do while (k<=nbpart)

        call random_number(xi)
        xi = (xi-0.5d0) * 10d0
        call random_number(yi)
        yi = (yi-0.5d0) * 10d0
        call random_number(zi)

        temm = (exp(-((xi-2d0)**2+yi**2)/2d0) &
        &      +exp(-((xi+2d0)**2+yi**2)/2d0))/2d0

        if (temm>=zi) then
            self%v(1,k) = xi
            self%v(2,k) = yi
            k = k + 1
        end if

    end do

end subroutine init_particles_2d


subroutine init_particles_3d( self, nbpart, mesh )

    type(particles_t), intent(out)  :: self
    integer(8),        intent(in)   :: nbpart
    type(mesh_t),      intent(in)   :: mesh

    real(8)              :: dimx, dimy, dimz
    integer              :: m
    integer              :: i, j, k
    real(8), parameter   :: eps = 1.d-12
    real(8)              :: xi, yi, wi, zi 
    real(8)              :: temm
    integer              :: nseed
    integer, allocatable :: seed(:)
    real(8)              :: pi

    pi = 4d0 * atan(1d0)

    self%nbpart = nbpart

    allocate(self%x(3,nbpart))
    allocate(self%v(3,nbpart))
    allocate(self%e(3,nbpart))
    allocate(self%b(nbpart))
    allocate(self%t(nbpart))

    dimx = mesh%xmax - mesh%xmin
    dimy = mesh%ymax - mesh%ymin
    dimz = mesh%zmax - mesh%zmin

    self%w = dimx * dimy * dimz / nbpart

    nseed = 33
    allocate(seed(nseed))
    
    seed = [-1584649339, -1457681104, 1579121008,   -819547200, &
              249798090, -517237887,   177452147,   -981503238, &
             1418301473,  1989625004, 2065424384,   -296364178, &
             1658790794, -435188152, -1643185032,   1461389312, &
             1869073641,  1321930686,  483734018,   1269936416, &
            -1999561453,  906251506,   782514880,    428753705, &
            -2031262823,  263953581,  1026600222,  -1118515860, &
             1633712916, -464192498, -1860714528,   1436611533, 0]
    
    call random_seed(put  = seed)

    call random_number(self%x)
    self%x(3,:) = mesh%zmin + dimz * self%x(3,:)

    m = 1
    do while (m <= nbpart)
        call random_number(xi)
        xi = 9.d0 * xi
        call random_number(yi)
        yi = 2.d0 * pi * yi
        call random_number(zi)
        zi = (1.0d0+0.02d0) * zi
        temm = (1d0+0.02d0*cos(4d0*yi))*exp(-5d0*(xi-4.8d0)**2)
        if (temm>=zi) then
            self%x(1,m) = cos(yi)*xi + 9d0
            self%x(2,m) = sin(yi)*xi + 9d0
            m = m + 1
        end if
    end do

    m = 1
    do while (m <= nbpart)
        call random_number(xi)
        xi = (xi-0.5d0) * 8d0
        call random_number(yi)
        yi = (yi-0.5d0) * 8d0
        call random_number(wi)
        wi = (wi-0.5d0) * 8d0
        call random_number(zi)
        temm = exp(-2d0*(xi**2+yi**2+wi**2))
        if (temm >= zi) then
            self%v(1, m) = xi
            self%v(2, m) = yi
            self%v(3, m) = wi
            m = m+1
        end if
    end do

end subroutine init_particles_3d

end module particles_m
