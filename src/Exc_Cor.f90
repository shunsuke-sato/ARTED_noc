!
!  Copyright 2016 ARTED developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
subroutine Exc_Cor
  use global_variables
  implicit none
  integer :: i
  real(8) :: trho,rs,rssq,rsln,V_xc,E_xc

  select case(cEex_Cor)
  case('PZ')
    call Exc_Cor_PZ
  case('TBmBJ')
    call Exc_Cor_TBmBJ_GS
  case default
    err_message='Err. in Eex_Cor' ; call err_finalize
  end select

end subroutine Exc_Cor
!=========================================================================================
subroutine Exc_Cor_PZ
  use global_variables
  implicit none
  integer :: i
  real(8) :: trho,rs,rssq,rsln,V_xc,E_xc

  do i=1,NL
    trho=rho_e(i)+1d-10
    rs=(3d0/(4*Pi*trho))**(1d0/3d0)
    V_xc=-4d0/3d0*0.4582d0/rs
    E_xc=-.4582d0/rs
    if (rs>1d0) then
      rssq=sqrt(rs)
      V_xc=V_xc+gammaU*(1d0+7d0/6d0*beta1U*rssq + 4d0/3d0*beta2U*rs)/(1d0+beta1U*rssq+beta2U*rs)**2
      E_xc=E_xc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
    else
      rsln=log(rs)
      V_xc=V_xc+AU*rsln+(BU-AU/3d0) + 2d0/3d0*CU*rs*rsln+(2d0*DU-CU)/3d0*rs
      E_xc=E_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
    end if
    Vxc(i)=V_xc
    Exc(i)=E_xc
  enddo

  return
end subroutine Exc_Cor_PZ
!=========================================================================================
Subroutine Exc_Cor_TBmBJ_GS !(GS_RT)
  use Global_Variables
  implicit none
  real(8),parameter :: alpha=-0.012d0,beta=1.023d0,gamma=0.80d0
  character(2) :: GS_RT
  real(8) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  real(8) :: c,tau_s_jrho,D_s_jrho,Q_s,rhs,x_s,b_s,Vx_BR,Vx_MBJ
  real(8) :: trho,rs,rhos,ec,dec_drhoa,dec_drhob
  integer :: i


  call rho_j_tau('GS',rho_s,tau_s,j_s,grho_s,lrho_s)
  if(cVal_mBJ > 1d-10) then
    c=cVal_mBJ ! use c-value given by input file
  else
!    c=sum(sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)/rho_s(:))*Hxyz/aLxyz
!    c=alpha+beta*sqrt(c)
    stop
  endif

  do i=1,NL
    tau_s_jrho=tau_s(i)-(j_s(i,1)**2+j_s(i,2)**2+j_s(i,3)**2)/rho_s(i)/2
    D_s_jrho=2*tau_s_jrho-0.25d0*(grho_s(i,1)**2+grho_s(i,2)**2+grho_s(i,3)**2)/rho_s(i)
    Q_s=(lrho_s(i)-2*gamma*D_s_jrho)/6d0
    rhs=2d0/3d0*pi**(2d0/3d0)*rho_s(i)**(5d0/3d0)/Q_s
    call BR_Newton(rhs,x_s)
    b_s=(x_s**3*exp(-x_s)/(8*pi*rho_s(i)))**(1d0/3d0)
    Vx_BR=-(1-exp(-x_s)-0.5d0*x_s*exp(-x_s))/b_s
    Vx_MBJ=c*Vx_BR+(3*c-2)/pi*sqrt(5d0/12d0)*sqrt(2*tau_s_jrho/rho_s(i))

    trho=rho_e(i)+1d-10
    rs=(3d0/(4*Pi*trho))**(1d0/3d0)
    Exc(i)=-.4582d0/rs
    Vxc(i)=Vx_MBJ
    rhos=rho_s(i)
    call PWc(rhos,rhos,ec,dec_drhoa,dec_drhob)
    Vxc(i)=Vxc(i)+ec+rho_e(i)*dec_drhoa
    Exc(i)=Exc(i)+ec
    Exc(i)=Exc(i)*trho
  enddo

  return
End Subroutine Exc_Cor_TBmBJ_GS
!=========================================================================================
Subroutine rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  use Global_Variables
  use PSE_variables
  implicit none
  character(2) :: GS_RT
  real(8) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  integer :: ikb,ik,ib,i
  real(8) :: tau_s_l(NL),j_s_l(NL,3),ss(3)
  complex(8) :: zs(3),zgk0psi(NL,3)
  integer :: thr_id,omp_get_thread_num
  real(8) :: fact

!  allocate(tau_s_l_omp(NL,0:NUMBER_THREADS-1),j_s_l_omp(NL,3,0:NUMBER_THREADS-1)) ! sato
  rho_s=rho_e*0.5d0 + rho_nlcc*0.5d0

  tau_s_l=0d0
  j_s_l=0d0

  if(GS_RT == 'GS')then

    K_point : do ik=NK_s,NK_e
      Band : do ib=1,NB
        tpsi(:)=zu_GS(:,ib,ik)
        call PSE_gk0psi_DFT(ik,zgk0psi)

        fact = (occ(ib,ik)*0.5d0)*0.5d0
        tau_s_l(:)=tau_s_l(:) + fact*( &
          abs(zgk0psi(:,1))**2 + abs(zgk0psi(:,2))**2 + abs(zgk0psi(:,3))**2)

        fact = (occ(ib,ik)*0.5d0)
        j_s_l(:,1)=j_s_l(:,1)+fact*real(conjg(tpsi(:))*(-zgk0psi(:,1)))
        j_s_l(:,2)=j_s_l(:,2)+fact*real(conjg(tpsi(:))*(-zgk0psi(:,2)))
        j_s_l(:,3)=j_s_l(:,3)+fact*real(conjg(tpsi(:))*(-zgk0psi(:,3)))

      end do Band
    end do K_point

!$omp end parallel

  else  if(GS_RT == 'RT')then
    stop
  else
    stop
  end if

  call MPI_ALLREDUCE(tau_s_l,tau_s,NL,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(j_s_l,j_s,NL*3,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)

  tau_s = tau_s + 0.5d0*tau_nlcc

  tpsi = rho_s
  call PSE_grho_lrho_DFT(zgk0psi)
  grho_s(:,1) = real(zgk0psi(:,1))
  grho_s(:,2) = real(zgk0psi(:,2))
  grho_s(:,3) = real(zgk0psi(:,3))
  lrho_s(:)= real(tpsi)


  return
End Subroutine rho_j_tau


  SUBROUTINE PWc(rhoa,rhob,ec,dec_drhoa,dec_drhob)
  implicit none
  real(8),parameter :: Pi=3.141592653589793d0
  real(8),intent(IN) :: rhoa,rhob
  real(8),intent(OUT) :: ec,dec_drhoa,dec_drhob
  real(8),parameter :: pec0=1d0,pec1=1d0,par=1d0
  real(8),parameter :: Aec0=0.031091d0,Aec1=0.015545d0,Aar=0.016887d0
  real(8),parameter :: a1ec0=0.21370d0,a1ec1=0.20548d0,a1ar=0.11125d0
  real(8),parameter :: b1ec0=7.5957d0,b1ec1=14.1189d0,b1ar=10.357d0
  real(8),parameter :: b2ec0=3.5876d0,b2ec1=6.1977d0,b2ar=3.6231d0
  real(8),parameter :: b3ec0=1.6382d0,b3ec1=3.3662d0,b3ar=0.88026d0
  real(8),parameter :: b4ec0=0.49294d0,b4ec1=0.62517d0,b4ar=0.49671d0
  real(8),parameter :: f2d=1.709921d0,eps=1d-15
  real(8) rs,sqrs,zeta,fzeta,dfzeta,Q0ec0,Q0ec1,Q0ar,Q1ec0,Q1ec1,Q1ar
  real(8) Qdec0,Qdec1,Qdar,Gec0,Gec1,Gar,dGec0,dGec1,dGar,decdrs,decdzeta

  rs=(3d0/(4d0*Pi*(rhoa+rhob)))**(1d0/3d0)
  sqrs=sqrt(rs)
  zeta=(rhoa-rhob)/(rhoa+rhob)
  if(abs(zeta-1d0) < eps) zeta=1d0-eps
  if(abs(zeta+1d0) < eps) zeta=-1d0+eps
  fzeta=((1+zeta)**(4d0/3d0)+(1-zeta)**(4d0/3d0)-2)/(2d0**(4d0/3d0)-2d0)
  dfzeta=4d0/3d0*((1+zeta)**(1d0/3d0)-(1-zeta)**(1d0/3d0))/(2d0**(4d0/3d0)-2d0)

  Q0ec0=-2*Aec0*(1+a1ec0*rs)
  Q0ec1=-2*Aec1*(1+a1ec1*rs)
  Q0ar =-2*Aar *(1+a1ar *rs)
  Q1ec0=2*Aec0*(b1ec0*sqrs+b2ec0*rs+b3ec0*rs*sqrs+b4ec0*rs**(pec0+1))
  Q1ec1=2*Aec1*(b1ec1*sqrs+b2ec1*rs+b3ec1*rs*sqrs+b4ec1*rs**(pec1+1))
  Q1ar =2*Aar *(b1ar *sqrs+b2ar *rs+b3ar *rs*sqrs+b4ar *rs**(par +1))
  Qdec0=Aec0*(b1ec0/sqrs+2*b2ec0+3*b3ec0*sqrs+2*(pec0+1)*b4ec0*rs**pec0)
  Qdec1=Aec1*(b1ec1/sqrs+2*b2ec1+3*b3ec1*sqrs+2*(pec1+1)*b4ec1*rs**pec1)
  Qdar =Aar *(b1ar /sqrs+2*b2ar +3*b3ar *sqrs+2*(par +1)*b4ar *rs**par )
  Gec0=Q0ec0*log(1+1d0/Q1ec0)
  Gec1=Q0ec1*log(1+1d0/Q1ec1)
  Gar =Q0ar *log(1+1d0/Q1ar )
  dGec0=-2*Aec0*a1ec0*log(1+1d0/Q1ec0)-Q0ec0*Qdec0/(Q1ec0**2+Q1ec0)
  dGec1=-2*Aec1*a1ec1*log(1+1d0/Q1ec1)-Q0ec1*Qdec1/(Q1ec1**2+Q1ec1)
  dGar =-2*Aar *a1ar *log(1+1d0/Q1ar )-Q0ar *Qdar /(Q1ar **2+Q1ar )

  ec=Gec0-Gar*fzeta/f2d*(1-zeta**4)+(Gec1-Gec0)*fzeta*zeta**4
  decdrs=dGec0*(1-fzeta*zeta**4)+dGec1*fzeta*zeta**4-dGar*fzeta/f2d*(1-zeta**4)
  decdzeta=dfzeta*(-Gar/f2d*(1-zeta**4)+(Gec1-Gec0)*zeta**4) &
&         +4*zeta**3*(fzeta*(Gec1-Gec0)+Gar*fzeta/f2d)

  dec_drhoa=(-rs/3*decdrs-(zeta-1)*decdzeta)/(rhoa+rhob)
  dec_drhob=(-rs/3*decdrs-(zeta+1)*decdzeta)/(rhoa+rhob)

  return
end SUBROUTINE PWc
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine PSE_gk0psi_DFT(ik,zgk0psi)
  use global_variables
  use PSE_variables
  implicit none
  complex(8) :: zgk0psi(NL,3),ztmp(0:NL1-1,0:NL2-1,0:NL3-1)
  integer :: ik
  integer :: i,ix,iy,iz
  real(8) :: kAc2_2
  integer :: ilma,ia,j
  complex(8) :: uVpsi


  do i=1,NL
    zft2(iLx(1,i),iLx(2,i),iLx(3,i))=tpsi(i)
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(cexp_x3(:,iz)*zft2(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(cexp_x2(:,iy)*zft1(ix,:,iz))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(cexp_x1(:,ix)*zft2(:,iy,iz))
  end do
  end do
  end do

  ztmp = zft1
! for x
  zft1(:,:,:)=( -Grad_x_zI(:,:,:) + kAc0_Cvec(1,ik) )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    zgk0psi(i,1)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do


! for y
  zft1(:,:,:)=( -Grad_y_zI(:,:,:) + kAc0_Cvec(2,ik) )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    zgk0psi(i,2)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do

! for z
  zft1(:,:,:)=( -Grad_z_zI(:,:,:) + kAc0_Cvec(3,ik) )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    zgk0psi(i,3)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do

  return
end subroutine PSE_gk0psi_DFT
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine PSE_grho_lrho_DFT(zgk0psi)
  use global_variables
  use PSE_variables
  implicit none
  complex(8) :: zgk0psi(NL,3),ztmp(0:NL1-1,0:NL2-1,0:NL3-1)
  integer :: i,ix,iy,iz
  real(8) :: kAc2_2
  integer :: ilma,ia,j
  complex(8) :: uVpsi


  do i=1,NL
    zft2(iLx(1,i),iLx(2,i),iLx(3,i))=tpsi(i)
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(cexp_x3(:,iz)*zft2(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(cexp_x2(:,iy)*zft1(ix,:,iz))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(cexp_x1(:,ix)*zft2(:,iy,iz))
  end do
  end do
  end do

  ztmp = zft1
! for x
  zft1(:,:,:)=( Grad_x_zI(:,:,:)/zI  )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    zgk0psi(i,1)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do


! for y
  zft1(:,:,:)=( Grad_y_zI(:,:,:) /zI )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    zgk0psi(i,2)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do

! for z
  zft1(:,:,:)=( Grad_z_zI(:,:,:)/zI )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    zgk0psi(i,3)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do

! for Lap
  zft1(:,:,:)=- ( Grad_x_zI(:,:,:)**2 + Grad_y_zI(:,:,:)**2 + Grad_z_zI(:,:,:)**2 &
    )*ztmp(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    tpsi(i)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do


  return
end subroutine PSE_grho_lrho_DFT
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
SUBROUTINE BR_Newton(rhs,x_s)
  implicit none
  real(8),intent(IN) :: rhs
  real(8),intent(OUT) :: x_s
  integer iter
  real(8) :: xmin,xmax,x,fx,dfx

! find xmax
  xmin=0d0
  x=1d0
  do
    fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
    if(fx < 0) then
      xmax=x
      exit
    endif
    x=x*2
  enddo
! bi-section
  do
    x=0.5d0*(xmin+xmax)
    fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
    if(fx < 0) then
      xmax=x
    else
      xmin=x
    endif
    if(xmax-xmin < 1d-4) exit
  enddo
! Newton-Raphson
  do iter=1,7
    fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
    dfx=(1-2d0/3d0*x)*exp(-2d0/3d0*x)/rhs-1d0
    x=x-fx/dfx
  enddo
  x_s=x
  return
end subroutine BR_Newton

