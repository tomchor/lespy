
program test
implicit none
interface 
       subroutine correlate4d(Vars, nxc, nyc, vCorr)
          implicit none
          real, intent(in), dimension(:,:) :: Vars
          real, intent(out),dimension(:,:), allocatable :: vCorr
          integer, intent(in) :: nxc, nyc
       end subroutine
end interface
real, dimension(3,3) :: a
real, dimension(:,:), allocatable :: corr

a = reshape( (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), (/ 3, 3 /))
print '(3f4.1)', a(1,:)
print '(3f4.1)', a(2,:)
print '(3f4.1)', a(3,:)
print *
call correlate4D(a, 2, 2, corr)
end program test
 


SUBROUTINE correlate4D(Vars, nxc, nyc, vCorr)

IMPLICIT NONE
real(kind=8), intent(in), dimension(:,:) :: Vars
integer, intent(in) :: nxc, nyc
integer :: ii, i, iix, iiy, iv, it, ix, iy
integer, dimension(2*nxc+1) :: xdel
integer, dimension(2*nyc+1) :: ydel
real(kind=8), dimension(:,:), allocatable :: vCorr

allocate(vCorr(2*nxc+1, 2*nyc+1))

print*,Vars
print*,size(Vars)

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

