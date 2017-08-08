!---------------------------------------------------------------------------------------------------
! interpolate initial conditions--in physical space
! this is not intended to be a subroutine...
!
! 3.28.13, STS: Modified for output files w/ theta as scalar, and for simulations being run w/ MPI.
!---------------------------------------------------------------------------------------------------
Program interp
Implicit None

!Variable Declarations -----------------------------------------------------------------------------
    
integer                        :: ix,iy,iz
integer,parameter                :: nx_s = 160, ny_s =160, nz_s = 75+1        !# nodes on small grid (w/ MPI)
integer,parameter                :: nx_b = 320, ny_b =320, nz_b = 150+1        !# nodes on big grid (w/ MPI)

integer                        :: S_s, S_b, S_s2, S_b2, i

character(len=36)                :: path = '/data/1/tomaschor/LES/convB/restart/'
character(len=30)                 :: input_file = 'vel_tt00400000.out.coarse'        !Coarse grid
character(len=30)                 :: input_file_temp = 'temp_tt00400000.out.coarse'        !Coarse grid
character(len=23)                 :: output_file ='vel_tt00400000.out'        !Restart file for fine grid
character(len=23)                 :: output_file_temp ='temp_tt00400000.out'        !Restart file for fine grid

real(kind=8)                    :: dummy

!Variables on coarse grid
real(kind=8),dimension(nx_s+2,ny_s,1:nz_s)    :: u_s, v_s, w_s, theta_s, RHSx_s, RHSy_s, RHSz_s, RHSt_s
real(kind=8),dimension(nx_s+2,ny_s,1)    :: sgs_t3_s                !
real(kind=8),dimension(nx_s,ny_s)        :: psi_m_s
real(kind=8),dimension(nx_s+2,ny_s)        :: psi_m_s_aux
real(kind=8),dimension(nx_s+2,ny_s,nz_s)    :: cs_s, FLM_s,FMM_s,FQN_s,FNN_s

!Variables on fine grid
real(kind=8),dimension(nx_b+2,ny_b,1:nz_b)    :: u_b,v_b,w_b,theta_b,RHSx_b,RHSy_b,RHSz_b,RHSt_b
real(kind=8),dimension(nx_b+2,ny_b,1)    :: sgs_t3_b    
real(kind=8),dimension(nx_b,ny_b)        :: psi_m_b
real(kind=8),dimension(nx_b+2,ny_b)        :: psi_m_b_aux
real(kind=8),dimension(nx_b+2,ny_b,nz_b)    :: cs_b, FLM_b,FMM_b,FQN_b,FNN_b

!Vertical profiles (for debugging purposes)
real(kind=8),dimension(nz_s)            :: u_s_avg,v_s_avg,w_s_avg,theta_s_avg,RHSx_s_avg,RHSy_s_avg,RHSz_s_avg,RHSt_s_avg
real(kind=8),dimension(nz_s)            :: cs_s_avg,FLM_s_avg,FMM_s_avg,FQN_s_avg,FNN_s_avg
real(kind=8),dimension(3,nz_s)            :: z_s

real(kind=8),dimension(nz_b)            :: u_b_avg,v_b_avg,w_b_avg,theta_b_avg,RHSx_b_avg,RHSy_b_avg,RHSz_b_avg,RHSt_b_avg
real(kind=8),dimension(nz_b)            :: cs_b_avg,FLM_b_avg,FMM_b_avg,FQN_b_avg,FNN_b_avg
real(kind=8),dimension(3,nz_b)            :: z_b

!Begin Code ----------------------------------------------------------------------------------------
!Initialize
sgs_t3_s(:,:,:) = 0.D0
sgs_t3_b(:,:,:) = 0.D0
psi_m_b = 0.D0

!NOTE: Will probably have to do 2D interpolation for psi_m if sfc heat flux is nonzero. Here not an issue, since
!we're imposing heat flux @ canopy top and sfc. Obukhov length is +infty. 

write(*,*)'Going to interpolate ',input_file,'and ',input_file_temp
write(*,*)'Coarse grid resolution: ',nx_s,' x ',ny_s,' x ',nz_s
write(*,*)'Fine grid resolution: ',nx_b,' x ',ny_b,' x ',nz_b


write(*,*)'Input file'
open(unit=10,file=path//input_file,form='unformatted',action='read')
!read(10)u_s(:, :, 1:nz_s), v_s(:, :, 1:nz_s), w_s(:, :, 1:nz_s), theta_s(:,:,1:nz_s), &
!          RHSx_s(:, :, 1:nz_s), RHSy_s(:, :, 1:nz_s), RHSz_s(:, :, 1:nz_s),&
!          RHSt_s(:,:,1:nz_s), sgs_t3_s(:,:,1), psi_m_s,Cs_s, FLM_s, FMM_s,&
read(10) u_s(:, :, 1:nz_s), v_s(:, :, 1:nz_s), w_s(:, :, 1:nz_s), &
          RHSx_s(:, :, 1:nz_s), RHSy_s(:, :, 1:nz_s), RHSz_s(:, :, 1:nz_s),&
          Cs_s, FLM_s, FMM_s, FQN_s, FNN_s
close(10)

write(*,*)'Input file'
open(unit=12,file=path//input_file_temp,form='unformatted',action='read')
read(12) theta_s(:,:,1:nz_s), RHSt_s(:,:,1:nz_s), sgs_t3_s(:,:,1), psi_m_s
close(12)


   

write(*,*)'u'
call interpolate(u_s,nx_s,ny_s,nz_s,u_b,nx_b,ny_b,nz_b)
write(*,*)'v'
call interpolate(v_s,nx_s,ny_s,nz_s,v_b,nx_b,ny_b,nz_b)
write(*,*)'w'
call interpolate(w_s,nx_s,ny_s,nz_s,w_b,nx_b,ny_b,nz_b)
write(*,*)'Theta'
call interpolate(theta_s,nx_s,ny_s,nz_s,theta_b,nx_b,ny_b,nz_b)
write(*,*)'RHSx'
call interpolate(RHSx_s,nx_s,ny_s,nz_s,RHSx_b,nx_b,ny_b,nz_b)
write(*,*)'RHSy'
call interpolate(RHSy_s,nx_s,ny_s,nz_s,RHSy_b,nx_b,ny_b,nz_b)
write(*,*)'RHSz'
call interpolate(RHSz_s,nx_s,ny_s,nz_s,RHSz_b,nx_b,ny_b,nz_b)
write(*,*)'RHSt'
call interpolate(RHSt_s,nx_s,ny_s,nz_s,RHSt_b,nx_b,ny_b,nz_b)
write(*,*)'Cs'
call interpolate(Cs_s,nx_s,ny_s,nz_s,Cs_b,nx_b,ny_b,nz_b)
write(*,*)'psi_m'
!call interpolate(psi_m_s,nx_s,ny_s,1,psi_m_b,nx_b,ny_b,1)
! TOMAS CHOR: had to change this for the interpolation to work as it is
psi_m_s_aux(1:nx_s,:)=psi_m_s(1:nx_s,:)
psi_m_b_aux=-9999
call interpolate(psi_m_s_aux,nx_s,ny_s,1,psi_m_b_aux,nx_b,ny_b,1)
psi_m_b(1:nx_b,:)=psi_m_b_aux(1:nx_b,:)
write(*,*)'sgs_t3'
call interpolate(sgs_t3_s,nx_s,ny_s,1,sgs_t3_b,nx_b,ny_b,1)
write(*,*)'FLM'
call interpolate(FLM_s,nx_s,ny_s,nz_s,FLM_b,nx_b,ny_b,nz_b)
write(*,*)'FMM'
call interpolate(FMM_s,nx_s,ny_s,nz_s,FMM_b,nx_b,ny_b,nz_b)
write(*,*)'FQN'
call interpolate(FQN_s,nx_s,ny_s,nz_s,FQN_b,nx_b,ny_b,nz_b)
write(*,*)'FNN'
call interpolate(FNN_s,nx_s,ny_s,nz_s,FNN_b,nx_b,ny_b,nz_b)

u_b(nx_b+1:nx_b+2,:,:) = 0.D0
v_b(nx_b+1:nx_b+2,:,:) = 0.D0
w_b(nx_b+1:nx_b+2,:,:) = 0.D0
theta_b(nx_b+1:nx_b+2,:,:) = 0.D0
RHSx_b(nx_b+1:nx_b+2,:,:) = 0.D0
RHSy_b(nx_b+1:nx_b+2,:,:) = 0.D0
RHSz_b(nx_b+1:nx_b+2,:,:) = 0.D0
RHSt_b(nx_b+1:nx_b+2,:,:) = 0.D0
Cs_b(nx_b+1:nx_b+2,:,:) = 1E-32
sgs_t3_b(nx_b+1:nx_b+2,:,:) = 0.D0
FLM_b(nx_b+1:nx_b+2,:,:) = 1E-32
FMM_b(nx_b+1:nx_b+2,:,:) = 1E-32
FQN_b(nx_b+1:nx_b+2,:,:) = 1E-32
FNN_b(nx_b+1:nx_b+2,:,:) = 1E-32

S_s = nx_s*ny_s*nz_s
S_b = nx_b*ny_b*nz_b
S_s2 = nx_s*ny_s
S_b2 = nx_b*ny_b

print*, 'u', SUM(u_s)/S_s,SUM(u_b)/S_b
print*, 'v', SUM(v_s)/S_s,SUM(v_b)/S_b
print*, 'w', SUM(w_s)/S_s,SUM(w_b)/S_b
print*, 'RHSx', SUM(RHSx_s)/S_s,SUM(RHSx_b)/S_b
print*, 'RHSy', SUM(RHSy_s)/S_s,SUM(RHSy_b)/S_b
print*, 'RHSz', SUM(RHSz_s)/S_s,SUM(RHSz_b)/S_b
print*, 'Cs', SUM(Cs_s)/S_s,SUM(Cs_b)/S_b
print*, 'FLM', SUM(FLM_s)/S_s,SUM(FLM_b)/S_b
print*, 'FMM', SUM(FMM_s)/S_s,SUM(FMM_b)/S_b
print*, 'FQN', SUM(FQN_s)/S_s,SUM(FQN_b)/S_b
print*, 'FNN', SUM(FNN_s)/S_s,SUM(FNN_b)/S_b
print*, 'theta', SUM(theta_s)/S_s,SUM(theta_b)/S_b
print*, 'RHSt', SUM(RHSt_s)/S_s,SUM(RHSt_b)/S_b
print*, 'sgs_t3', SUM(sgs_t3_s)/S_s2,SUM(sgs_t3_b)/S_b2
print*, 'psi_m', SUM(psi_m_s)/S_s2,SUM(psi_m_b)/S_b2

write(*,*)'Writing output file'
open(unit=11,file=output_file,form='unformatted',action='write')
write(11) u_b(:, :, 1:nz_b), v_b(:, :, 1:nz_b), w_b(:, :, 1:nz_b), &
          rhsx_b(:, :, 1:nz_b), rhsy_b(:, :, 1:nz_b), rhsz_b(:, :, 1:nz_b),&
          Cs_b, flm_b, fmm_b,&
          fqn_b, fnn_b
close(11)


open(unit=13,file=output_file_temp,form='unformatted',action='write')
write(13) theta_b(:,:,1:nz_b), RHSt_b(:,:,1:nz_b), sgs_t3_b(:,:,1), psi_m_b
close(13)


!STOP    !Stop after interpolation... no debugging for now

!--------------------------------------------------------------
!DEBUG BLOCK --------------------------------------------------
!--------------------------------------------------------------
!Calculate average vertical profiles on coarse and fine grids

write(*,*)'Averaging on coarse grid'
!Coarse grid
Call avg_xy(u_s,nx_s,ny_s,nz_s,u_s_avg)
Call avg_xy(v_s,nx_s,ny_s,nz_s,v_s_avg)
Call avg_xy(w_s,nx_s,ny_s,nz_s,w_s_avg)
Call avg_xy(theta_s,nx_s,ny_s,nz_s,theta_s_avg)
Call avg_xy(RHSx_s,nx_s,ny_s,nz_s,RHSx_s_avg)
Call avg_xy(RHSy_s,nx_s,ny_s,nz_s,RHSy_s_avg)
Call avg_xy(RHSz_s,nx_s,ny_s,nz_s,RHSz_s_avg)
Call avg_xy(RHSt_s,nx_s,ny_s,nz_s,RHSt_s_avg)
Call avg_xy(cs_s,nx_s,ny_s,nz_s,cs_s_avg)
Call avg_xy(FLM_s,nx_s,ny_s,nz_s,FLM_s_avg)
Call avg_xy(FMM_s,nx_s,ny_s,nz_s,FMM_s_avg)
Call avg_xy(FQN_s,nx_s,ny_s,nz_s,FQN_s_avg)
Call avg_xy(FNN_s,nx_s,ny_s,nz_s,FNN_s_avg)

DO i=1,nz_s !Calculate z-coordinates
    ! uv-nodes
    z_s(1,i)=(i-0.5D0)/nz_s
    ! w-nodes
    z_s(2,i)=(1.D0*i)/nz_s
END DO
! cs-nodes
z_s(3,1)=z_s(1,1)
z_s(3,2:nz_s)=z_s(2,1:nz_s-1)  


write(*,*)'Averaging fine grid'
!Fine grid
Call avg_xy(u_b,nx_b,ny_b,nz_b,u_b_avg)
Call avg_xy(v_b,nx_b,ny_b,nz_b,v_b_avg)
Call avg_xy(w_b,nx_b,ny_b,nz_b,w_b_avg)
Call avg_xy(theta_b,nx_b,ny_b,nz_b,theta_b_avg)
Call avg_xy(RHSx_b,nx_b,ny_b,nz_b,RHSx_b_avg)
Call avg_xy(RHSy_b,nx_b,ny_b,nz_b,RHSy_b_avg)
Call avg_xy(RHSz_b,nx_b,ny_b,nz_b,RHSz_b_avg)
Call avg_xy(RHSt_b,nx_b,ny_b,nz_b,RHSt_b_avg)
Call avg_xy(cs_b,nx_b,ny_b,nz_b,cs_b_avg)
Call avg_xy(FLM_b,nx_b,ny_b,nz_b,FLM_b_avg)
Call avg_xy(FMM_b,nx_b,ny_b,nz_b,FMM_b_avg)
Call avg_xy(FQN_b,nx_b,ny_b,nz_b,FQN_b_avg)
Call avg_xy(FNN_b,nx_b,ny_b,nz_b,FNN_b_avg)

DO i=1,nz_b !Calculate z-coordinates
    ! uv-nodes
    z_b(1,i)=(i-0.5D0)/nz_b
    ! w-nodes
    z_b(2,i)=(1.D0*i)/nz_b
END DO
! cs-nodes
z_b(3,1)=z_b(1,1)
z_b(3,2:nz_b)=z_b(2,1:nz_b-1)  

write(*,*)'Writing mean profiles'
!Write mean vertical profiles to file
open(unit=20,file='coarse_profiles.txt',action='write',form='formatted')
do i=1,nz_s
    write(20,'(16(F,X))')z_s(1,i),z_s(2,i),z_s(3,i),u_s_avg(i),v_s_avg(i),w_s_avg(i),theta_s_avg(i),&
        RHSx_s_avg(i),RHSy_s_avg(i),RHSz_s_avg(i),RHSt_s_avg(i),cs_s_avg(i),&
        FLM_s_avg(i),FMM_s_avg(i),FQN_s_avg(i),FNN_s_avg(i)
end do
close(20)

open(unit=21,file='fine_profiles.txt',action='write',form='formatted')
do i=1,nz_b
    write(21,'(16(F,X))')z_b(1,i),z_b(2,i),z_b(3,i),u_b_avg(i),v_b_avg(i),w_b_avg(i),theta_b_avg(i),&
        RHSx_b_avg(i),RHSy_b_avg(i),RHSz_b_avg(i),RHSt_b_avg(i),cs_b_avg(i),&
        FLM_b_avg(i),FMM_b_avg(i),FQN_b_avg(i),FNN_b_avg(i)
end do
close(21)






!--------------------------------------------------------------
!END DEBUG BLOCK ----------------------------------------------
!--------------------------------------------------------------
print*,'Dont forget to issue "rm *avg_stats.dat"'

end program interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------     
subroutine avg_xy(dat,nx,ny,nz,avg_dat)
implicit none

integer :: nx,ny,nz,i,j,k
real(kind=8) :: dat(nx+2,ny,nz)
real(kind=8) :: avg_dat(nz)

avg_dat(:) = 0.D0    

do k=1,nz
do i=1,nx
do j=1,ny
    avg_dat(k) = avg_dat(k) + dat(i,j,k)/(real(nx,8)*real(ny,8))
end do
end do
end do

return
end subroutine avg_xy
!-------------

subroutine interpolate(u_s,nx_s,ny_s,nz_s,u_b,nx_b,ny_b,nz_b)

implicit none
integer :: nx_s, ny_s, nz_s, nx_b, ny_b, nz_b
double precision, dimension(nx_b+2,ny_b,nz_b) :: u_b
double precision, dimension(nx_s+2,ny_s,nz_s) :: u_s
double precision :: rx, ry, rz, P1, P2, P3, P4, P5, P6
double precision :: dx, dy, dz
integer :: jx, jy, jz, i, j, k, i_wrap, j_wrap                !!????????

rx = real(nx_s)/real(nx_b)
ry = real(ny_s)/real(ny_b)
rz = real(nz_s-1)/real(nz_b-1)

do jz=1,nz_b-1  ! we already know the last point!
    k = floor((jz-1)*rz) + 1
    dz = (jz-1)*rz-(k-1)
    do jy=1,ny_b
        j = floor((jy-1)*ry) + 1
        j_wrap = modulo(j+1-1,ny_s) + 1
        dy = (jy-1)*ry-(j-1)
        do jx=1,nx_b
            i = floor((jx-1)*rx) + 1
            i_wrap = modulo(i+1-1,nx_s) + 1
            dx = (jx-1)*rx-(i-1)

            P1 = (1.-dx)*u_s(i,j,k) + dx*u_s(i_wrap,j,k)
            P2 = (1.-dx)*u_s(i,j,k+1) + dx*u_s(i_wrap,j,k+1)
            P3 = (1.-dx)*u_s(i,j_wrap,k) + dx*u_s(i_wrap,j_wrap,k)
            P4 = (1.-dx)*u_s(i,j_wrap,k+1) + dx*u_s(i_wrap,j_wrap,k+1)
            P5 = (1.-dy)*P1 + dy*P3 
            P6 = (1.-dy)*P2 + dy*P4

            u_b(jx,jy,jz) = (1.-dz)*P5 + dz*P6

        end do
    end do
end do

! still have to do nz_b level
do jy=1,ny_b
    j = floor((jy-1)*ry) + 1
    j_wrap = modulo(j+1-1,ny_s) + 1
    dy = (jy-1)*ry-(j-1)
    do jx=1,nx_b
        i = 1 + floor((jx-1)*rx)
        i_wrap = modulo(i+1-1,nx_s) + 1
        dx = (jx-1)*rx-(i-1)

        P1 = (1.-dx)*u_s(i,j,nz_s) + dx*u_s(i_wrap,j,nz_s)
        P2 = (1.-dx)*u_s(i,j_wrap,nz_s) + dx*u_s(i_wrap,j_wrap,nz_s)

        u_b(jx,jy,nz_b) = (1.-dy)*P1 + dy*P2

    end do

end do

end subroutine interpolate


