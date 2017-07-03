

SUBROUTINE correlate4D(Vars, nxc, nyc)
! Correlates 4D array assuming that x, y are dims 3 and 4
! Returns a 4D array with shape nv, nt, 2*nxc+1, 2*nyc+1

IMPLICIT NONE
real(kind=8), intent(in), dimension(:,:,:,:) :: Vars
integer, intent(in) :: nxc, nyc
integer :: ii, i, iix, iiy, ix, iy, dims(4), nv, nt, nx, ny
integer, dimension(2*nxc+1) :: xdel
integer, dimension(2*nyc+1) :: ydel
real(kind=8), dimension(:,:,:,:), allocatable :: vCorr, rolled
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
!Mean = Vars.mean(axis=(2,3))
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
do iix=1,size(xdel)
    ix=xdel(iix)
    do iiy=1,size(ydel)
        print*,iix,iiy
        rolled = roll(Vars, ix, axis=2)
        rolled = roll(rolled, iy, axis=3)
        do iv=1,size(nv)
            do it=1,size(nt)
                print *,'got here'
                dummy = mean(Vars(iv,it) * rolled(iv,it))
                !vCorr[iv,it,iix,iiy] = mean(Vars[iv,it] * rolled[iv,it]) / (Mean[iv,it]**2.)
                vCorr(iv,it,iix,iiy) = (sum(dummy)/size(dummy)) / (Mean[iv,it]**2.)
            enddo
        enddo
    enddo
enddo


END SUBROUTINE


!function mean3D(array, axis=0)
!implicit none
!real(kind=8), dimension(:,:,:,:) :: mean3D, array
!integer :: axis

!if (dim==0) then
!    mean3D = sum(array)/size(array)
!else
!    mean3D = sum(array, dim=axis)/size(array,dim=axis)
!endif
!endfunction

