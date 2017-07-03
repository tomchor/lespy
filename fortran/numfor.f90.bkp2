SUBROUTINE correlate4D(Vars, nxc, nyc, nv, nt, vCorr)
! Correlates 4D array assuming that x, y are dims 3 and 4
! Returns a 4D array with shape nv, nt, 2*nxc+1, 2*nyc+1
IMPLICIT NONE
real(kind=8), intent(in), dimension(:,:,:,:) :: Vars
real(kind=8) :: dummysize
integer, intent(in) :: nxc, nyc
integer :: ii, i, iix, iiy, iv, it, dims(4), nv, nt, nx, ny
integer, dimension(2*nxc+1) :: xdel
integer, dimension(2*nyc+1) :: ydel
real(kind=8), intent(out) :: vCorr(nv, nt, 2*nxc+1, 2*nyc+1)
real(kind=8), dimension(:,:,:,:), allocatable :: rolled, prerolled
real(kind=8), dimension(:,:), allocatable :: Mean
real(kind=8), dimension(:,:,:), allocatable :: Mean3d

dims = shape(Vars)
nx=dims(3)
ny=dims(4)
dummysize=nx*ny
allocate(rolled(nv, nt, nx, ny))
allocate(prerolled(nv, nt, nx, ny))
allocate(Mean3d(nv, nt, nx))
allocate(Mean(nv, nt))

!------
! Calculate mean, var and correlation matrix for calculations
! Mean = Vars.mean(axis=(2,3))
Mean3d = sum(Vars, dim=4)/size(Vars, dim=4)
Mean = sum(Mean3d, dim=3)/size(Mean3d, dim=3)
!------

!------
! Here we set the size of the correlation matrix
! xdel = _np.arange(-nxc,nxc+1)
ii=1
do i=-nxc,+nxc
    xdel(ii)=i
    ii=ii+1
enddo
! ydel = _np.arange(-nyc,nyc+1)
ii=1
do i=-nyc,+nyc
    ydel(ii)=i
    ii=ii+1
enddo
!------

!---------
! Calculate the correlation
do iiy=1,size(ydel)
    print*,iiy,' of', size(ydel)
    prerolled = cshift(Vars, ydel(iiy), dim=4)
    do iix=1,size(xdel)
        rolled = cshift(prerolled, xdel(iix), dim=3)
        forall (it=1:nt)
            forall (iv=1:nv)
                vCorr(iv,it,iix,iiy) = (sum(Vars(iv,it,:,:) * rolled(iv,it,:,:))/dummysize) / (Mean(iv,it)**2)
            endforall
        endforall
    enddo
enddo
!---------

END SUBROUTINE

SUBROUTINE correlate5D(Vars, nxc, nyc, nv, nt, ns, vCorr)
! Correlates 4D array assuming that x, y are dims 3 and 4
! Returns a 4D array with shape nv, nt, 2*nxc+1, 2*nyc+1
IMPLICIT NONE
real(kind=8), intent(in), dimension(:,:,:,:) :: Vars
real(kind=8) :: dummysize
integer, intent(in) :: nxc, nyc
integer :: ii, i, iix, iiy, iv, it, dims(4), nv, nt, nx, ny, ns
integer, dimension(2*nxc+1) :: xdel
integer, dimension(2*nyc+1) :: ydel
real(kind=8), intent(out) :: vCorr(nv, nt, 2*nxc+1, 2*nyc+1)
real(kind=8), dimension(:,:,:,:), allocatable :: rolled, prerolled
real(kind=8), dimension(:,:), allocatable :: Mean
real(kind=8), dimension(:,:,:), allocatable :: Mean3d

dims = shape(Vars)
nx=dims(3)
ny=dims(4)
dummysize=nx*ny
allocate(rolled(nv, nt, nx, ny))
allocate(prerolled(nv, nt, nx, ny))
allocate(Mean3d(nv, nt, nx))
allocate(Mean(nv, nt))

!------
! Calculate mean, var and correlation matrix for calculations
! Mean = Vars.mean(axis=(2,3))
Mean3d = sum(Vars, dim=4)/size(Vars, dim=4)
Mean = sum(Mean3d, dim=3)/size(Mean3d, dim=3)
!------

!------
! Here we set the size of the correlation matrix
! xdel = _np.arange(-nxc,nxc+1)
ii=1
do i=-nxc,+nxc
    xdel(ii)=i
    ii=ii+1
enddo
! ydel = _np.arange(-nyc,nyc+1)
ii=1
do i=-nyc,+nyc
    ydel(ii)=i
    ii=ii+1
enddo
!------

!---------
! Calculate the correlation
do iiy=1,size(ydel)
    print*,iiy,' of', size(ydel)
    prerolled = cshift(Vars, ydel(iiy), dim=4)
    do iix=1,size(xdel)
        rolled = cshift(prerolled, xdel(iix), dim=3)
        forall (it=1:nt)
            forall (iv=1:nv)
                vCorr(iv,it,iix,iiy) = (sum(Vars(iv,it,:,:) * rolled(iv,it,:,:))/dummysize) / (Mean(iv,it)**2)
            endforall
        endforall
    enddo
enddo
!---------

END SUBROUTINE


