SUBROUTINE read_vel_sc(nx, ny, nz_tot, u_tot, v_tot, w_tot, theta_tot)
! Reads binary vel_sc files from the model
implicit none

INTEGER, INTENT(IN):: nx
INTEGER, INTENT(IN):: ny
INTEGER, INTENT(IN):: nz_tot
!INTEGER, INTENT(IN), DIMENSION(::

!LOGICAL, INTENT(IN), OPTIONAL:: S_FLAG, PCON_FLAG
LOGICAL:: SFLAG, PCONFLAG

REAL(kind=8), intent(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: u_tot
REAL(kind=8), intent(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: v_tot
REAL(kind=8), intent(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: w_tot
REAL(kind=8), intent(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: theta_tot
REAL(kind=8), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: PCon_tot
REAL(kind=8), DIMENSION(nx,ny)   :: ustar_avg
REAL(kind=8), DIMENSION(nx,ny)   :: P_surf_flux   ! Surface pollen flux     
REAL(kind=8), DIMENSION(nx,ny)   :: deposition    ! Surface pollen deposition (this is actually the net flux=deposition-source)
REAL(kind=8), DIMENSION(nx,ny)   :: P_surf_flux_dep! Surface pollen flux for Cr=0 everywhere
REAL(kind=8), DIMENSION(nx,ny)   :: Real_dep      ! This is the real deposition (using Cr=0 everywhere)


SFLAG=.TRUE.
PCONFLAG=.FALSE.
!print *,PCONFLAG, SFLAG
!print*,PRESENT(s_flag)
!IF (PRESENT(S_FLAG)) SFLAG=S_FLAG
!IF (PRESENT(PCON_FLAG)) PCONFLAG=PCON_FLAG
!print *,PCONFLAG, SFLAG


OPEN(UNIT = 12, FORM = 'unformatted', FILE = 'vel_sc00400000.out')

IF (SFLAG .AND. PCONFLAG) THEN
    read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), &
             theta_tot(:,:,1:nz_tot), PCon_tot(:,:,1:nz_tot), &
             deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
             P_surf_flux(:,:), P_surf_flux_dep(:,:)
ELSEIF (SFLAG) THEN
    read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot)
ELSE
    read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot)
ENDIF
close(12)

END
