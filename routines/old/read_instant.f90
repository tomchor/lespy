SUBROUTINE read_vel_sc(fname, nx, ny, nz_tot, u_tot, v_tot, w_tot, theta_tot)
! Reads binary vel_sc files from the model
implicit none

INTEGER, INTENT(IN):: nx
INTEGER, INTENT(IN):: ny
INTEGER, INTENT(IN):: nz_tot
!INTEGER, INTENT(IN), DIMENSION(::

!LOGICAL, INTENT(IN), OPTIONAL:: S_FLAG, PCON_FLAG
LOGICAL:: SFLAG, PCONFLAG
CHARACTER(len=200)::fname

REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: u_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: v_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: w_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: theta_tot
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


OPEN(UNIT = 12, FORM = 'unformatted', FILE = fname)
!OPEN(UNIT = 12, FORM = 'unformatted', FILE = 'vel_sc00400000.out')

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



SUBROUTINE read_vel_pcon(fname, nx, ny, nz_tot, n_con, u_tot, v_tot, w_tot, theta_tot,&
        PCon_tot, ustar_avg, P_surf_flux, deposition, P_surf_flux_dep, Real_dep)
! Reads binary vel_sc files from the model
implicit none

INTEGER, INTENT(IN):: nx
INTEGER, INTENT(IN):: ny
INTEGER, INTENT(IN):: nz_tot, n_con
!INTEGER, INTENT(IN), DIMENSION(::

!LOGICAL, INTENT(IN), OPTIONAL:: S_FLAG, PCON_FLAG
LOGICAL:: SFLAG, PCONFLAG
CHARACTER(len=200)::fname

REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: u_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: v_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: w_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: theta_tot
!REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: PCon_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot, n_con) :: PCon_tot
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: ustar_avg
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: P_surf_flux   ! Surface pollen flux     
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: deposition    ! Surface pollen deposition (this is actually the net flux=deposition-source)
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: P_surf_flux_dep! Surface pollen flux for Cr=0 everywhere
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: Real_dep      ! This is the real deposition (using Cr=0 everywhere)


SFLAG=.TRUE.
PCONFLAG=.TRUE.
!print *,PCONFLAG, SFLAG
!print*,PRESENT(s_flag)
!IF (PRESENT(S_FLAG)) SFLAG=S_FLAG
!IF (PRESENT(PCON_FLAG)) PCONFLAG=PCON_FLAG
!print *,PCONFLAG, SFLAG


!OPEN(UNIT = 12, FORM = 'unformatted', FILE = '/data/1/tomaschor/LES/LV3_test/output/vel_sc01516000.out')
OPEN(UNIT = 12, FORM = 'unformatted', FILE = fname)

IF (SFLAG .AND. PCONFLAG) THEN
    read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), &
             theta_tot(:,:,1:nz_tot), PCon_tot(:,:,1:nz_tot,:), &
!             theta_tot(:,:,1:nz_tot), PCon_tot(:,:,1:nz_tot), &
             deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
             P_surf_flux(:,:), P_surf_flux_dep(:,:)
ELSEIF (SFLAG) THEN
    read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot)
ELSE
    read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot)
ENDIF
close(12)

END
