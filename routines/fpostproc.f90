

SUBROUTINE POSTPROC(dir, nz, z_scale, u_star, dt_dim, nt, tini, tend, T_scale, outdata)

IMPLICIT NONE

! Parameters
REAL(KIND=8),PARAMETER     :: kappa=0.4D0         ! Von Karman constant

INTEGER         :: nz            ! Number of z-nodes
REAL(KIND=8)    :: z_scale       ! Length Scale
REAL(KIND=8)    :: u_star        ! u* scale used in program
REAL(KIND=8)    :: T_scale       ! theta scale
REAL(KIND=8)    :: dt_dim ! dimensional dt
REAL(KIND=8)    :: dt     ! non-dimensional timestep

INTEGER      :: nt               ! Total number of time outputs
INTEGER      :: tend             ! Final timestep for averaging (output/2.5)
INTEGER      :: tini             ! Initial timestep for averaging (output/2.5)

! Main variables
CHARACTER(len=200)           :: file                ! File name
CHARACTER(len=100)         :: dir                 ! Folder
REAL(KIND=8),ALLOCATABLE   :: z(:,:)              ! Vertical coordinates
REAL(KIND=8),ALLOCATABLE   :: avgtx(:,:)          ! Averages

! Variables to determine ustar
INTEGER                    :: num                 ! Number of points for time averaging
REAL(KIND=8),ALLOCATABLE   :: t(:), tt(:)         ! Time
REAL(KIND=8),ALLOCATABLE   :: dat(:,:)          ! Data
REAL(KIND=8),ALLOCATABLE   :: ustar(:)            ! Data averaged in x
REAL(KIND=8),ALLOCATABLE   :: tstar(:)            ! Data averaged in x
REAL(KIND=8)               :: tstar_avg_curr            ! Data averaged in x and t
REAL(KIND=8),ALLOCATABLE   :: ustarr(:)           ! Data averaged in x

! Output variable
REAL(KIND=8),INTENT(OUT),DIMENSION (23,nz)   :: outdata       ! Output array

! Auxiliar variables
INTEGER                    :: i             ! Counters
REAL                       :: aux           ! Auxiliar
logical :: exist ! file exist flag

dt=dt_dim*u_star/z_scale            ! Timestep

!print *, 'tinit = ', tini
!print *, 'tend = ', tend

!
! BEGINNING CODE   
!

! Memory allocation
ALLOCATE(z(3,nz),avgtx(25,nz))
ALLOCATE(t(nt),tt(nt),dat(nt,3),ustar(nt),tstar(nt),ustarr(nt))


! Read data from file
file=TRIM(dir)//'/aver_txz.out'
OPEN(UNIT=10,FILE=file,STATUS='OLD',ACTION='READ') 
file=TRIM(dir)//'/aver_tyz.out'
OPEN(UNIT=11,FILE=file,STATUS='OLD',ACTION='READ')    
! Read only first vertical level
DO i=1,nt
    READ(10,*)t(i),dat(i,1)
    READ(11,*)t(i),dat(i,2)
    tt(i)=INT(t(i)/dt)
END DO


CLOSE(10)
CLOSE(11)

ustar(:)=SQRT(dat(:,1)**2+dat(:,2)**2) !u*^2
ustar=SQRT(ustar) !u*
print*,'chaging USTAR from', sum(ustar)/nt,' to', u_star
ustar=u_star

num=COUNT(tt(:)>=tini .AND. tt(:)<=tend)
aux = SUM(ustar(:),MASK=(tt(:)>=tini .AND. tt(:)<=tend))/num!*u_star

PRINT*,'u_star_med = ',SUM(ustar)/nt!*u_star
PRINT*,'u_star = ',aux


! Read data from file
file=TRIM(dir)//'/aver_sgs_t3.out'
inquire(file=file, exist=exist)
if (exist) then
    OPEN(UNIT=10,FILE=file,STATUS='OLD',ACTION='READ') 
    ! Read only first vertical level
    DO i=1,nt
      READ(10,*)t(i),dat(i,1)
    END DO
    CLOSE(10)

    tstar(:)=dat(:,1)/ustar(:)
    tstar_avg_curr=SUM(tstar(:),MASK=(tt(:)>=tini .AND. tt(:)<=tend))/num*T_scale
    !tstar for average period
    
    PRINT*,'theta_star = ',SUM(tstar)/nt*T_scale
    PRINT*,'wtheta_0 = ',SUM(tstar)/nt*T_scale*SUM(ustar)/nt*u_star
else
    tstar_avg_curr=1.d0
end if

ustarr = ustar

! Average data from all files
file=TRIM(dir)//'/aver_u.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(1,:))

file=TRIM(dir)//'/aver_v.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(2,:))

file=TRIM(dir)//'/aver_w.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(3,:))

file=TRIM(dir)//'/aver_dudz.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(4,:))

file=TRIM(dir)//'/aver_dvdz.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(5,:))
avgtx(5,1) = avgtx(5,2)


! Normalize by u_star^2
ustar=ustarr*ustarr
file=TRIM(dir)//'/aver_u2.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(6,:))

file=TRIM(dir)//'/aver_v2.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(7,:))

file=TRIM(dir)//'/aver_w2.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(8,:))

file=TRIM(dir)//'/aver_uw.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(9,:))

file=TRIM(dir)//'/aver_vw.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(10,:))

file=TRIM(dir)//'/aver_txx.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(11,:))

file=TRIM(dir)//'/aver_tyy.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(12,:))

file=TRIM(dir)//'/aver_tzz.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(13,:))

file=TRIM(dir)//'/aver_txz.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(14,:))

file=TRIM(dir)//'/aver_tyz.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(15,:))

! These are not normalized
ustar=1.D0
file=TRIM(dir)//'/aver_Cs.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(16,:))

file=TRIM(dir)//'/aver_beta_sgs.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(17,:))

file=TRIM(dir)//'/aver_betaclip_sgs.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(18,:))

file=TRIM(dir)//'/aver_Cs_Ssim.out'
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(19,:))



! Average <u> - dimensional
file=TRIM(dir)//'/aver_u.out'
ustar=1.0/u_star
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(20,:))

! Average <v> - dimensional
file=TRIM(dir)//'/aver_v.out'
ustar=1.0/u_star
CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(21,:))



! Some calculations for output

! Determine z-coordinates
DO i=1,nz
  
  ! uv-nodes
    z(1,i)=(i-0.5D0)/nz
    
    ! w-nodes
    z(2,i)=(1.D0*i)/nz
    
END DO

! Average temperature terms - theta_star

! Average <T> - dimensional
file=TRIM(dir)//'/aver_theta.out'
ustar=1.0
inquire(file=file, exist=exist)
if (exist) then
  CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(22,:))
end if

file=TRIM(dir)//'/aver_dTdz.out'
ustar=1.0/u_star
if(SUM(tstar)/nt.eq.0.0) ustar=1.0/T_scale !nao normaliza se tstar eh igual a zero
inquire(file=file, exist=exist)
if (exist) then
  CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(23,:))
end if


! Average <wT> - nondimensional
ustar=ustarr*tstar_avg_curr/T_scale
if(SUM(tstar)/nt.eq.0.0) ustar=1.0/T_scale !nao normaliza se tstar eh igual a zero


file=TRIM(dir)//'/aver_wt.out'
inquire(file=file, exist=exist)
if (exist) then
  CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(24,:))
end if

file=TRIM(dir)//'/aver_sgs_t3.out'
inquire(file=file, exist=exist)
if (exist) then
  CALL AVG_T(nz,nt,dt,tini,tend,file,ustar,avgtx(25,:))
end if

 

! cs-nodes
z(3,1)=z(1,1)
z(3,2:nz)=z(2,1:nz-1)

! Output results
DO i=1,nz
  if (i == 1) then

        outdata(1,i) = z(1,i)*z_scale   ! z
        outdata(2,i) = z(2,i)*z_scale   ! z
        outdata(3,i) = z(3,i)*z_scale   ! z
        outdata(4,i) = avgtx(1,i)   ! <U>/u*
        outdata(5,i) =        avgtx(2,i)   ! <V>/u*
        outdata(6,i) =        avgtx(3,i)   ! <W>/u*
        outdata(7,i) =        avgtx(4,i)   ! <dUdz>/u**dz*
        outdata(8,i) =        avgtx(5,i)   ! <dVdz>/u**dz*
        outdata(9,i) =        avgtx(20,i)  ! <U>
        outdata(10,i) =        avgtx(21,i)  ! <V>
        outdata(11,i) =        avgtx(22,i)  ! <Theta>
        outdata(12,i) =        avgtx(6,i)-avgtx(1,i)*avgtx(1,i)    !<u^2>/u*^2
        outdata(13,i) =        avgtx(7,i)-avgtx(2,i)*avgtx(2,i)    !<v^2>/u*^2
        outdata(14,i) =        avgtx(8,i)-avgtx(3,i)*avgtx(3,i)    !<w^2>/u*^2
        outdata(15,i) =        avgtx(9,i)-avgtx(1,i)*avgtx(3,i)+avgtx(14,i) !<uw>/u*^2
        outdata(16,i) =        avgtx(10,i)-avgtx(2,i)*avgtx(3,i)+avgtx(15,i) !<vw>/u*^2
        outdata(17,i) =        avgtx(24,i)-avgtx(3,i)*avgtx(22,i)/tstar_avg_curr+avgtx(25,i) !<wT>/u*T*
        outdata(18,i) =        avgtx(16,i)   ! cs
        outdata(19,i) =        avgtx(17,i)   ! beta
        outdata(20,i) =        avgtx(18,i)   ! beta_clip
        outdata(21,i) =        avgtx(19,i)   ! cs_rns
        outdata(22,i) =        avgtx(14,i)  !<txz>/u*^2
        outdata(23,i) =        avgtx(15,i) !<tyz>/u*^2

  else

        outdata(1,i) = z(1,i)*z_scale   ! z
        outdata(2,i) = z(2,i)*z_scale   ! z
        outdata(3,i) = z(3,i)*z_scale   ! z
        outdata(4,i) = avgtx(1,i)   ! <U>/u*
        outdata(5,i) =        avgtx(2,i)   ! <V>/u*
        outdata(6,i) =        avgtx(3,i)   ! <W>/u*
        outdata(7,i) =        avgtx(4,i)   ! <dUdz>/u**dz*
        outdata(8,i) =        avgtx(5,i)   ! <dVdz>/u**dz*
        outdata(9,i) =        avgtx(20,i)  ! <U>
        outdata(10,i) =        avgtx(21,i)  ! <V>
        outdata(11,i) =        avgtx(22,i)  ! <Theta>
        outdata(12,i) =        avgtx(6,i)-avgtx(1,i)*avgtx(1,i)    !<u^2>/u*^2
        outdata(13,i) =        avgtx(7,i)-avgtx(2,i)*avgtx(2,i)    !<v^2>/u*^2
        outdata(14,i) =        avgtx(8,i)-avgtx(3,i)*avgtx(3,i)    !<w^2>/u*^2
        outdata(15,i) =        avgtx(9,i)-0.5*(avgtx(1,i)+avgtx(1,i-1))*avgtx(3,i)+avgtx(14,i)  !<uw>/u*^2
        outdata(16,i) =        avgtx(10,i)-0.5*(avgtx(2,i)+avgtx(2,i-1))*avgtx(3,i)+avgtx(15,i) !<vw>/u*^2
        outdata(17,i) =        avgtx(24,i)-avgtx(3,i)*avgtx(22,i)/tstar_avg_curr+avgtx(25,i)    !<wT>/u*T*
        outdata(18,i) =        avgtx(16,i)   ! cs
        outdata(19,i) =        avgtx(17,i)   ! beta
        outdata(20,i) =        avgtx(18,i)   ! beta_clip
        outdata(21,i) =        avgtx(19,i)   ! cs_rns
        outdata(22,i) =        avgtx(14,i)  !<txz>/u*^2
        outdata(23,i) =        avgtx(15,i) !<tyz>/u*^2




  endif
END DO
CLOSE(10) 

END SUBROUTINE




SUBROUTINE AVG_T(nz,nt,dt,tini,tend,file,scale_factor,avg)
!    ################################################################################
!    ##                                  AVG_X_T                                   ##
!    ##                                                                            ##
!    ##                               Developed by                                 ##
!    ##                             Marcelo Chamecki                               ##
!    ##                                03/30/2006                                  ##
!    ##                                                                            ##
!    ################################################################################
!    PURPOSE:  This routine average LES output in x and time.
!    ################################################################################

! No implicit variables
IMPLICIT none

! Main variables
INTEGER                    :: nz              ! Number of z-nodes
INTEGER                    :: nt              ! Number of time outputs
INTEGER                    :: num             ! Number of points for time averaging
INTEGER                    :: tini            ! Initial time for averaging
INTEGER                    :: tend            ! Final time for averaging
REAL(KIND=8)               :: dt              ! Time step
REAL(KIND=8)               :: time            ! Time
CHARACTER(LEN=*)           :: file            ! File name
INTEGER,DIMENSION(nt)            :: t         ! Time
REAL(KIND=8),DIMENSION(nt,nz) :: dat       ! Data
REAL(KIND=8),DIMENSION(nz)       :: avg       ! Data averaged in t
REAL(KIND=8),DIMENSION(nt)       :: scale_factor     ! u*(t)

! Auxiliar variables
INTEGER                    :: i         ! Counters

! Read data from file
OPEN(UNIT=10,FILE=trim(file),STATUS='OLD',ACTION='READ') 

DO i=1,nt
    READ(10,*)time,dat(i,1:nz)

    ! Time to iteration number
    t(i)=INT(time/dt)
    
    ! Proper normalization
    dat(i,:)=dat(i,:)/scale_factor(i)
    
END DO
CLOSE(10)

! Average in time
num=COUNT(t(:)>=tini .AND. t(:)<=tend)
DO i=1,nz
    avg(i)=SUM(dat(:,i),MASK=(t(:)>=tini .AND. t(:)<=tend))/num
END DO

RETURN
END SUBROUTINE AVG_T
