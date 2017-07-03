SUBROUTINE read_binary(fname, nx, ny, nz_tot, n_con, s_flag, pcon_flag, &
        flag_endless, source, u_tot, v_tot, w_tot, theta_tot, PCon_tot)

IMPLICIT NONE

INTEGER, INTENT(IN):: nx
INTEGER, INTENT(IN):: ny
INTEGER, INTENT(IN):: nz_tot, n_con

LOGICAL, INTENT(IN) :: s_flag, pcon_flag, flag_endless
CHARACTER(len=200), INTENT(IN)     ::fname
CHARACTER(len=10), INTENT(IN)      ::source

REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: u_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: v_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: w_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: theta_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot, n_con) :: PCon_tot
REAL(kind=8), DIMENSION(nx,ny)   :: ustar_avg
REAL(kind=8), DIMENSION(nx,ny)   :: P_surf_flux   ! Surface pollen flux     
REAL(kind=8), DIMENSION(nx,ny)   :: deposition    ! Surface pollen deposition (this is actually the net flux=deposition-source)
REAL(kind=8), DIMENSION(nx,ny)   :: P_surf_flux_dep! Surface pollen flux for Cr=0 everywhere
REAL(kind=8), DIMENSION(nx,ny)   :: Real_dep      ! This is the real deposition (using Cr=0 everywhere)


OPEN(UNIT = 12, FORM = 'unformatted', FILE = fname)

IF (flag_endless) THEN
    !print*,'Reading endless output file'
    IF (source=='vel') THEN
        !print*,'reading vel_t'
        read(12) u_tot, v_tot(:,:,:), w_tot(:,:,:)
    ELSE IF (source=='temp') THEN
        !print*,'reading temp_t'
        read(12) theta_tot
    ELSE IF (source=='con_tt') THEN
        !print*,'reading con_tt'
        READ(12)  PCon_tot(:,:,:,:)
        !print*, 'maxvals', maxval(PCon_tot), minval(pcon_tot)
    ENDIF

ELSE
    print*,'Reading vel_sc file'
    IF (s_flag .AND. pcon_flag) THEN
        !print*,'reading long one with s_flag and pcon_flag'
        print*,nx,ny,nz_tot, n_con
        read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:,:,:), &
             theta_tot(:,:,1:nz_tot), PCon_tot(:,:,1:nz_tot,:), &
             deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
             P_surf_flux(:,:), P_surf_flux_dep(:,:)
    ELSEIF (s_flag) THEN
        read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), theta_tot(:,:,1:nz_tot)
    ELSE
        read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot)
    ENDIF
ENDIF
close(12)

END


