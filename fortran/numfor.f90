
SUBROUTINE correlate4D(Vars, nxc, nyc, vCorr)
! Correlates 4D array assuming that x, y are dims 3 and 4
! Returns a 4D array with shape nv, nt, 2*nxc+1, 2*nyc+1

IMPLICIT NONE
real(kind=8), intent(in), dimension(:,:,:,:) :: Vars
integer, intent(in) :: nxc, nyc
integer :: ii, i, iix, iiy, iv, it, ix, iy, dims(4), nv, nt, nx, ny
integer, dimension(2*nxc+1) :: xdel
integer, dimension(2*nyc+1) :: ydel
real(kind=8), intent(out), allocatable :: vCorr(:,:,:,:)
real(kind=8), dimension(:,:,:,:), allocatable :: rolled
real(kind=8), dimension(:,:), allocatable :: dummy, Mean
real(kind=8), dimension(:,:,:), allocatable :: Mean3d


dims = shape(Vars)
nv=dims(1)
nt=dims(2)
nx=dims(3)
ny=dims(4)
allocate(vCorr(nv, nt, 2*nxc+1, 2*nyc+1))
allocate(rolled(nv, nt, nx, ny))
allocate(Mean3d(nv, nt, nx))
allocate(Mean(nv, nt))
allocate(dummy(nx, ny))

!------
! Calculate mean, var and correlation matrix for calculations
! Mean = Vars.mean(axis=(2,3))
Mean3d = sum(Vars, dim=4)/size(Vars, dim=4)
Mean = sum(Mean3d, dim=3)/size(Mean3d, dim=3)
!------
print*,Vars
print*,size(Vars)

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
do iix=1,size(xdel)
    ix=xdel(iix)
    do iiy=1,size(ydel)
!        print*,iix,iiy
        rolled = cshift(Vars, ix, dim=3)
        rolled = cshift(rolled, iy, dim=4)
!        print*
!        print '(3i3)', rolled(1,:,:,:)
!        print '(3i3)', rolled(2,:,:,:)
!        print '(3i3)', rolled(3,:,:,:)
        do iv=1,nv
            do it=1,nt
!                print *,'got here'
                dummy = Vars(iv,it,:,:) * rolled(iv,it,:,:)
                !vCorr[iv,it,iix,iiy] = mean(Vars[iv,it] * rolled[iv,it]) / (Mean[iv,it]**2.)
                vCorr(iv,it,iix,iiy) = (sum(dummy)/size(dummy)) / (Mean(iv,it)**2.)
            enddo
        enddo
    enddo
enddo


END SUBROUTINE


