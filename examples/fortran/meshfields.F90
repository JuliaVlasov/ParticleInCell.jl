module mesh_fields_m

implicit none

type :: mesh_t

    real(8) :: xmin 
    real(8) :: xmax 
    integer :: nx   
    real(8) :: dx   
    real(8) :: ymin 
    real(8) :: ymax 
    integer :: ny   
    real(8) :: dy   
    real(8) :: zmin 
    real(8) :: zmax 
    integer :: nz   
    real(8) :: dz   

end type mesh_t

type :: fields_2d_t

    type(mesh_t)         :: mesh 
    real(8), allocatable :: e(:,:,:)
    real(8), allocatable :: rho(:,:)

end type fields_2d_t

type :: fields_3d_t

    type(mesh_t)         :: mesh 
    real(8), allocatable :: e(:,:,:,:)
    real(8), allocatable :: rho(:,:,:)

end type fields_3d_t

interface init_mesh

     module procedure init_mesh_2d
     module procedure init_mesh_3d

end interface init_mesh

interface init_fields

     module procedure init_fields_2d
     module procedure init_fields_3d

end interface init_fields

contains

    subroutine init_mesh_2d( self, xmin, xmax, nx, ymin, ymax, ny )

        type(mesh_t) :: self

        real(8), intent(in) :: xmin 
        real(8), intent(in) :: xmax 
        integer, intent(in) :: nx   
        real(8), intent(in) :: ymin 
        real(8), intent(in) :: ymax 
        integer, intent(in) :: ny   

        self%xmin = xmin
        self%xmax = xmax
        self%nx   = nx
        self%ymin = ymin
        self%ymax = ymax
        self%ny   = ny
        self%dx   = (xmax - xmin) / real(nx,8)
        self%dy   = (ymax - ymin) / real(ny,8)

        self%zmin = 0d0
        self%zmax = 0d0
        self%nz   = 1
        self%dz   = 0d0

    end

    subroutine init_mesh_3d( self, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz )

        type(mesh_t) :: self

        real(8), intent(in) :: xmin 
        real(8), intent(in) :: xmax 
        integer, intent(in) :: nx   
        real(8), intent(in) :: ymin 
        real(8), intent(in) :: ymax 
        integer, intent(in) :: ny   
        real(8), intent(in) :: zmin 
        real(8), intent(in) :: zmax 
        integer, intent(in) :: nz   

        self%xmin = xmin
        self%ymin = ymin
        self%zmin = zmin

        self%xmax = xmax
        self%ymax = ymax
        self%zmax = zmax

        self%nx   = nx
        self%ny   = ny
        self%nz   = nz

        self%dx   = (xmax - xmin) / real(nx,8)
        self%dy   = (ymax - ymin) / real(ny,8)
        self%dz   = (zmax - zmin) / real(nz,8)

    end

    subroutine init_fields_2d( self, mesh )

         type(fields_2d_t) :: self
         type(mesh_t)      :: mesh
	    
         self%mesh = mesh

         allocate(self%e(2, mesh%nx+1, mesh%ny+1))
         allocate(self%rho(mesh%nx+1, mesh%ny+1))

    end

    subroutine init_fields_3d( self, mesh )

         type(fields_3d_t) :: self
         type(mesh_t)      :: mesh
	    
         self%mesh = mesh

         allocate(self%e(3, mesh%nx+1, mesh%ny+1, mesh%nz+1))
         allocate(self%rho(mesh%nx+1, mesh%ny+1, mesh%nz+1))

    end


end module mesh_fields_m
