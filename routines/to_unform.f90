

SUBROUTINE WRITE2FILE(dir, nz, z_scale, u_star, dt_dim, nt, tini, tend, T_scale, outdata)

IMPLICIT NONE

! Parameters
REAL(KIND=8),PARAMETER     :: kappa=0.4D0         ! Von Karman constant

INTEGER         :: nz            ! Number of z-nodes
REAL(KIND=8)    :: z_scale       ! Length Scale
REAL(KIND=8)    :: u_star        ! u* scale used in program
REAL(KIND=8)    :: T_scale       ! theta scale
REAL(KIND=8)    :: dt_dim ! dimensional dt
REAL(KIND=8)    :: dt     ! non-dimensional timestep

INTEGER      :: nt               ! Total number of time outputs
INTEGER      :: tend             ! Final timestep for averaging (output/2.5)
INTEGER      :: tini             ! Initial timestep for averaging (output/2.5)

! Main variables
CHARACTER(len=200)           :: file                ! File name
CHARACTER(len=100)         :: dir                 ! Folder
REAL(KIND=8),ALLOCATABLE   :: z(:,:)              ! Vertical coordinates
REAL(KIND=8),ALLOCATABLE   :: avgtx(:,:)          ! Averages

! Variables to determine ustar
INTEGER                    :: num                 ! Number of points for time averaging
REAL(KIND=8),ALLOCATABLE   :: t(:), tt(:)         ! Time
REAL(KIND=8),ALLOCATABLE   :: dat(:,:)          ! Data
REAL(KIND=8),ALLOCATABLE   :: ustar(:)            ! Data averaged in x
REAL(KIND=8),ALLOCATABLE   :: tstar(:)            ! Data averaged in x
REAL(KIND=8)               :: tstar_avg_curr            ! Data averaged in x and t
REAL(KIND=8),ALLOCATABLE   :: ustarr(:)           ! Data averaged in x

! Output variable
REAL(KIND=8),INTENT(OUT),DIMENSION (23,nz)   :: outdata       ! Output array

! Auxiliar variables
INTEGER                    :: i             ! Counters
REAL                       :: aux           ! Auxiliar
logical :: exist ! file exist flag

dt=dt_dim*u_star/z_scale            ! Timestep

!print *, 'tinit = ', tini
!print *, 'tend = ', tend

!
! BEGINNING CODE   
!

! Memory allocation
ALLOCATE(z(3,nz),avgtx(25,nz))
ALLOCATE(t(nt),tt(nt),dat(nt,3),ustar(nt),tstar(nt),ustarr(nt))

real, dimension(nx,ny,nz) :: array
  do i=1,nx
    do j=1,ny
       do k=1,nz
          array(i,j,k) = (i+j+k)
       end do
    end do
 end do
 inquire (IOLENGTH=rec_len) array  
 open(1,file='myarray.raw',status = 'unknown',
form='unformatted', access='direct',recl=rec_len)
write(1,rec=1) array
close(1)

END SUBROUTINE


