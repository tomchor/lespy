SUBROUTINE cond_density4D(Vars, nxc, nyc, nt, nv, vCorr)
! Correlates 4D array assuming that x, y are dims 1 and 2
! Returns a 4D array with shape 2*nxc+1, 2*nyc+1, nt, nv
IMPLICIT NONE
INTEGER, PARAMETER  ::  dp = SELECTED_REAL_KIND (13)
real(kind=8), intent(in), dimension(:,:,:,:) :: Vars
real(kind=8) :: dummysize
integer, intent(in) :: nxc, nyc
integer :: ii, i, iix, iiy, iv, it, dims(4), nv, nt, nx, ny
integer, dimension(2*nxc+1) :: xdel
integer, dimension(2*nyc+1) :: ydel
real(kind=8), intent(out) :: vCorr(2*nxc+1, 2*nyc+1, nt, nv)
real(kind=8), dimension(:,:,:,:), allocatable :: rolled, prerolled
real(kind=8), dimension(:,:), allocatable :: Mean
real(kind=8), dimension(:,:,:), allocatable :: Mean3d

dims = shape(Vars)
nx=dims(1)
ny=dims(1)
dummysize=nx*ny
allocate(rolled(nx, ny, nt, nv))
allocate(prerolled(nx, ny, nt, nv))
allocate(Mean3d(ny, nt, nv))
allocate(Mean(nt, nv))

!------
! Calculate mean, var and correlation matrix for calculations
Mean3d = sum(Vars, dim=1)/size(Vars, dim=1)
Mean = sum(Mean3d, dim=1)/size(Mean3d, dim=1)
!------

!------
! Here we set the size of the correlation matrix
ii=1
do i=-nxc,+nxc
    xdel(ii)=i
    ii=ii+1
enddo
ii=1
do i=-nyc,+nyc
    ydel(ii)=i
    ii=ii+1
enddo
!------

!---------
! Calculate the correlation
do iiy=1,size(ydel)
    print*,'fortran loop:',iiy,' of', size(ydel)
    prerolled = cshift(Vars, ydel(iiy), dim=2)
    do iix=1,size(xdel)
        rolled = cshift(prerolled, xdel(iix), dim=1)
        forall (iv=1:nv)
            forall (it=1:nt)
                vCorr(iix,iiy,it,iv) = sum(Vars(:,:,it,iv) * rolled(:,:,it,iv))/(dummysize * Mean(it,iv)**2)
            endforall
        endforall
    enddo
enddo
!---------

END SUBROUTINE



