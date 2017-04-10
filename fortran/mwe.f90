SUBROUTINE mwe(Vars, nxc, nyc, dims, vCorr)

IMPLICIT NONE
real(kind=8), dimension(:,:,:,:) :: Vars
integer :: nxc, nyc
integer, intent(in) :: dims(4)
real(kind=8), intent(out) :: vCorr(dims(1), dims(2), 2*nxc+1, 2*nyc+1)


print*,size(vCorr)
print*,size(Vars)

END SUBROUTINE


