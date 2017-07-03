SUBROUTINE read_binary(fname, nx, ny, nz_tot, n_con, s_flag, pcon_flag, &
        flag_endless, source, &
        u_tot, v_tot, w_tot, theta_tot,&
        PCon_tot, ustar_avg, P_surf_flux, deposition, P_surf_flux_dep, Real_dep)

IMPLICIT NONE

INTEGER, INTENT(IN):: nx
INTEGER, INTENT(IN):: ny
INTEGER, INTENT(IN):: nz_tot, n_con

LOGICAL, INTENT(IN) :: s_flag, pcon_flag, flag_endless
CHARACTER(len=200), INTENT(IN)     ::fname
CHARACTER(len=4), INTENT(IN)       ::source

REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: u_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: v_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: w_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: theta_tot
REAL(kind=8), INTENT(out), DIMENSION(2*(nx/2+1), ny, 1:nz_tot, n_con) :: PCon_tot
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: ustar_avg
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: P_surf_flux   ! Surface pollen flux     
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: deposition    ! Surface pollen deposition (this is actually the net flux=deposition-source)
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: P_surf_flux_dep! Surface pollen flux for Cr=0 everywhere
REAL(kind=8), INTENT(out), DIMENSION(nx,ny)   :: Real_dep      ! This is the real deposition (using Cr=0 everywhere)

REAL(kind=8), DIMENSION(nx,ny)   :: sgs_t3
REAL(kind=8), DIMENSION(2*(nx/2+1), ny, 1:nz_tot) :: RHS_T_tot!, RHSx_tot, RHSy_tot, RHSz_tot
!REAL(kind=8), DIMENSION(2*(nx/2+1), ny, nz_tot, n_con) :: RHS_PCon_tot
REAL(kind=8) :: psi_m!, Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot

OPEN(UNIT = 12, FORM = 'unformatted', FILE = fname)

IF (flag_endless) THEN
    print*,'Reading endless output file'
    IF (source=='vel') THEN
        print*,'reading vel_t'
        read(12) u_tot, v_tot(:,:,:), w_tot(:,:,:)
    ELSE IF (source=='temp') THEN
        print*,'reading temp_t'
        read(12) theta_tot!(:,:,:)!, RHS_T_tot(:,:,1:nz_tot), sgs_t3(:,:), psi_m
    ELSE IF (source=='con') THEN
        print*,'reading con_tt'
        READ(12)  PCon_tot(:,:,:,1)!, RHS_PCon_tot(:,:,1:nz_tot,1)
    ENDIF

ELSE
    print*,'Reading vel_sc file'
    IF (s_flag .AND. pcon_flag) THEN
        read(12) u_tot(:, :, 1:nz_tot), v_tot(:, :, 1:nz_tot), w_tot(:, :, 1:nz_tot), &
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










SUBROUTINE read_vel_sc(fname, nx, ny, nz_tot, u_tot, v_tot, w_tot, theta_tot)
! Reads binary vel_sc files from the model
implicit none

INTEGER, INTENT(IN):: nx
INTEGER, INTENT(IN):: ny
INTEGER, INTENT(IN):: nz_tot
!INTEGER, INTENT(IN), DIMENSION(::

!LOGICAL, INTENT(IN), OPTIONAL:: s_flag, pcon_flag
LOGICAL:: SFLAG, pcon_flag
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
pcon_flag=.FALSE.
!print *,pcon_flag, SFLAG
!print*,PRESENT(s_flag)
!IF (PRESENT(s_flag)) SFLAG=s_flag
!IF (PRESENT(pcon_flag)) pcon_flag=pcon_flag
!print *,pcon_flag, SFLAG


OPEN(UNIT = 12, FORM = 'unformatted', FILE = fname)
!OPEN(UNIT = 12, FORM = 'unformatted', FILE = 'vel_sc00400000.out')

IF (SFLAG .AND. pcon_flag) THEN
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

!LOGICAL, INTENT(IN), OPTIONAL:: s_flag, pcon_flag
LOGICAL:: SFLAG, pcon_flag
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
pcon_flag=.TRUE.
!print *,PCONFLAG, SFLAG
!print*,PRESENT(s_flag)
!IF (PRESENT(s_flag)) SFLAG=s_flag
!IF (PRESENT(pcon_flag)) PCONFLAG=pcon_flag
!print *,PCONFLAG, SFLAG


!OPEN(UNIT = 12, FORM = 'unformatted', FILE = '/data/1/tomaschor/LES/LV3_test/output/vel_sc01516000.out')
OPEN(UNIT = 12, FORM = 'unformatted', FILE = fname)

IF (SFLAG .AND. PCON_FLAG) THEN
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
