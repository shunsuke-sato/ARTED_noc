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
! Cvec: cartesian 
! Lvec: Lattice
! Rvec: Reciprocal lattice
module global_variables
  implicit none
!ARTED version
  character(50),parameter :: CODE_ver='crab.2014.09.13.0'

! constants
  real(8),parameter :: Pi=3.141592653589793d0
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
  real(8),parameter :: umass=1822.9d0

! DFT parameters
  character(10) :: cEex_Cor
  real(8) :: cVal_mBJ
  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
  real(8),parameter :: CU=0.002d0,DU=-0.0116d0

! Control parameter
  character(3) :: method ! 'PAW'(projector augmented wave) or 'PSE'(pseudo-potential) or 'CHK'(check)

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr
  integer :: NEW_COMM_WORLD,NEWPROCS,NEWRANK ! sato
  integer :: NK_ave,NG_ave,NK_s,NK_e,NG_s,NG_e
  integer :: NK_remainder,NG_remainder


! grid
  integer :: NL1,NL2,NL3,NL
  integer :: NK1,NK2,NK3,NK
  real(8) :: a_Cvec(3,3),a_Cvec_d(3,3),b_Cvec(3,3)
  real(8) :: norm_a_Cvec(3),norm_b_Cvec(3)
  real(8) :: mat_vv_a_Cvec(3,3)
  real(8) :: A_matrix(3,3),B_matrix(3,3),B_t_matrix(3,3)
  real(8) :: aL,aL1,aL2,aL3
  real(8),allocatable :: kAc_Rvec(:,:),kAc0_Rvec(:,:)
  real(8),allocatable :: kAc_Cvec(:,:),kAc0_Cvec(:,:)
  real(8) :: H1,H2,H3,H123,Vcell
  real(8) :: dk1,dk2,dk3
  real(8),allocatable :: Lx(:,:) !,Rxyz_Cvec(:,:)
  integer,allocatable :: iLx(:,:),iLx123(:,:,:)
  integer,allocatable :: ifdx1(:,:),ifdx2(:,:),ifdx3(:,:)  
  real(8),allocatable :: lap(:),nab(:)
  real(8),allocatable :: lap1(:),lap2(:),lap3(:)
  real(8),allocatable :: nab1(:),nab2(:),nab3(:)
  character(2) :: type_spatial_difference
  

! wave function, work
  integer :: NB,NB_TD
  real(8),allocatable :: occ_TD(:,:),occ_TD_l(:,:)
  complex(8),allocatable :: zu(:,:,:),zu_GS(:,:,:),zu_GS0(:,:,:)
  complex(8),allocatable :: tpsi(:),htpsi(:),tpsi_g1(:),tpsi_g2(:)

! density, potential
  real(8),allocatable :: rho_c(:),rho_p(:),rho_e(:),rho_e_l(:)
  real(8),allocatable :: rho_c_3D(:,:,:)
  real(8),allocatable :: Vloc(:),Vh(:),Vxc(:),Exc(:)

! physical quantity
  real(8),allocatable :: esp(:,:),esp_l(:,:)
  real(8) :: Etot_GS,Ekin_GS,Etot_RT

! material
  integer :: NI,NE,Nelec
  integer,allocatable :: Zatom(:),Kion(:)
  real(8),allocatable :: Rion_Lvec(:,:)
  real(8),allocatable :: occ(:,:)

! GS parameter
  integer :: Ncg,Nscf

! RT parameter
  character(1) :: option_IP='N'
  integer :: Nt
  real(8) :: dt


! laser
  character(2) :: Tr_Lo
  character(7) :: laser_type 
  real(8) :: IWcm2_1,tpulsefs_1,omegaev_1,phi_CEP_1,IWcm2_2,tpulsefs_2,omegaev_2,phi_CEP_2,Epdir_1(3),Epdir_2(3)
  real(8) :: dAc
  real(8) :: T1_T2fs
  real(8),allocatable :: Acext_Cvec(:,:),Acind_Cvec(:,:),Actot_Cvec(:,:),javt_Cvec(:,:)

! reentrance
  character(1) :: entrance_option
  real(8) :: Time_shoutdown

! control or/and I/O
  character(2) :: calc_mode
  character(200) :: err_message
  character(50) :: SYSname
  character(1) :: option_TD_band_dos
  integer :: Nstep_TD_band_dos

! Discrete Fourier Transformation
  complex(8),allocatable :: exp_x1(:,:),exp_x2(:,:),exp_x3(:,:),cexp_x1(:,:),cexp_x2(:,:),cexp_x3(:,:)
  complex(8),allocatable :: zft1(:,:,:),zft2(:,:,:),zft3(:,:,:)
  real(8),allocatable :: Lap_k(:,:,:),InLap_k(:,:,:),Grad_x_zI(:,:,:),Grad_y_zI(:,:,:),Grad_z_zI(:,:,:)

! Conjugate Gradient calculation
  complex(8),allocatable :: xk(:),hxk(:),pk(:),gk(:),pko(:)
  real(8),allocatable :: esp_var(:,:)

! subsupase diag
  complex(8),allocatable :: zutmp_diag(:,:),za_diag(:,:)

! linear mixing
  real(8),allocatable :: rho_e_old(:)

! predicotr-corrector
  integer :: Npred_corr
  complex(8),allocatable :: zu_t(:,:,:)
  real(8),allocatable :: Vloc_t(:)

! Modified Broyden's method
  integer,parameter :: MaxMem_MB = 8
  real(8),parameter :: alpha_MB = 0.35d0
  real(8),allocatable :: rho_MB_in(:,:),rho_MB_out(:,:)

! Basis expansion method
  integer :: NB_basis, NK_shift, NB_basis_main,NB_basis_shift
  real(8),allocatable :: kshift(:,:)
  complex(8),allocatable :: zu_basis(:,:,:)
  complex(8),allocatable :: zH_loc(:,:,:),zPi_loc(:,:,:)
  complex(8),allocatable :: zV_NL(:,:,:,:),zPi_NL(:,:,:,:)
  complex(8),allocatable :: zH_tot(:,:,:),zPi_tot(:,:,:)
  complex(8),allocatable :: zH0_tot(:,:,:),zdH_tot(:,:,:)
  real(8),allocatable :: H0_eigval(:,:)
  integer,parameter :: NAmax = 20
  real(8) :: Amax,dAmax
  real(8),allocatable :: Actot_BE(:),javt_BE(:)

  complex(8),allocatable :: zCt(:,:,:),zACt_tmp(:),ztCt_tmp(:)
  complex(8),allocatable :: zC_eig(:,:,:)
  integer,parameter :: NLanczos = 16
  complex(8),allocatable :: zLanCt(:,:,:),zACt_Lan(:,:),ztCt_Lan(:,:)
  real(8),parameter :: epsilon_Lan = 1d-12

! Houson basis decomposition
  logical,parameter :: switch_Houston_probe_decomposition = .false. !! .true.
  real(8),allocatable :: Ac_pump_BE(:),Ac_probe_BE(:)
end Module Global_Variables
