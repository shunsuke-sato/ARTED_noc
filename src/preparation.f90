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
subroutine preparation
  use global_variables
  use PSE_variables
  implicit none
  integer :: i,j
  character(10) :: cEex_Cor_tmp
! == MPI initialize
  call MPI_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,Myrank,ierr)
  NEW_COMM_WORLD=MPI_COMM_WORLD

  if(Myrank == 0)then
    call today
    write(*,'(2A)')'CODE ver. =',CODE_ver  
  end if
! == MPI initialize

!== input start
  if(myrank == 0)then
    write(*,'(A)') '!-----------------------------------------------------------------------------------------------------'
    write(*,'(A)') '<< Input date >>'

    read(*,*)calc_mode
    write(*,'(A,2x,A)')'calc_mode = ',calc_mode

    read(*,*)entrance_option
    write(*,'(A,2x,A)')'entrance_option = ',entrance_option
  end if

  call MPI_BCAST(calc_mode,64,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(entrance_option,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  if(entrance_option == 'N')then
!  else if(entrance_option == 'R')then
!    call Read_reentrance_data
!    return
!    write(*,'(A)') '!-----------------------------------------------------------------------------------------------------'
  else
    err_message='invalid entrance_option'
    call err_finalize
  end if
  
  if(myrank == 0)then
    read(*,*) Time_shoutdown
    read(*,*) SYSname
    read(*,*) cEex_Cor,cVal_mBJ
    read(*,*) ps_format 
    read(*,*) method, type_spatial_difference
    read(*,*) a_Cvec_d(1,1),a_Cvec_d(2,1),a_Cvec_d(3,1)
    read(*,*) a_Cvec_d(1,2),a_Cvec_d(2,2),a_Cvec_d(3,2)
    read(*,*) a_Cvec_d(1,3),a_Cvec_d(2,3),a_Cvec_d(3,3)
    read(*,*) aL,aL1,aL2,aL3
    read(*,*) NL1,NL2,NL3
    read(*,*) NK1,NK2,NK3
    read(*,*) NB,Nelec
    read(*,*) Ncg,Nscf
    read(*,*) Nt,dt,Npred_corr
    read(*,*) option_TD_band_dos,Nstep_TD_band_dos
    read(*,*) option_IP
    read(*,*) dAc
    read(*,*) laser_type,Tr_Lo
    read(*,*) IWcm2_1,tpulsefs_1,omegaev_1,phi_CEP_1
    read(*,*) Epdir_1(1),Epdir_1(2),Epdir_1(3)
    read(*,*) IWcm2_2,tpulsefs_2,omegaev_2,phi_CEP_2
    read(*,*) Epdir_2(1),Epdir_2(2),Epdir_2(3)
    read(*,*) T1_T2fs
    read(*,*) NI,NE
    write(*,'(A,2x,e16.6e3)') 'Time_shoutdown = ',Time_shoutdown
    write(*,'(A,2x,A)') 'SYSname = ',SYSname
    write(*,'(A,2x,A,2x,e16.6e3)') 'cEex_Cor, ,cVal_mBJ = ',cEex_Cor,cVal_mBJ
    write(*,*) 'ps_format =',ps_format 
    write(*,'(A,2x,A,2x,A)') 'method, type_spatial_difference',method, type_spatial_difference
    write(*,'(A,2x,3e16.6E3)') 'a_Cvec_d(1,1),a_Cvec_d(2,1),a_Cvec_d(3,1) = ',a_Cvec_d(1,1),a_Cvec_d(2,1),a_Cvec_d(3,1)
    write(*,'(A,2x,3e16.6E3)') 'a_Cvec_d(1,2),a_Cvec_d(2,2),a_Cvec_d(3,2) = ',a_Cvec_d(1,2),a_Cvec_d(2,2),a_Cvec_d(3,2)
    write(*,'(A,2x,3e16.6E3)') 'a_Cvec_d(1,3),a_Cvec_d(2,3),a_Cvec_d(3,3) = ',a_Cvec_d(1,3),a_Cvec_d(2,3),a_Cvec_d(3,3)
    write(*,'(A,2x,4e16.6E3)') 'aL,aL1,aL2,aL3 = ',aL,aL1,aL2,aL3
    write(*,'(A,4(2x,I0))') 'NL1,NL2,NL3 = ',NL1,NL2,NL3
    write(*,'(A,4(2x,I0))') 'NK1,NK2,NK3 = ',NK1,NK2,NK3
    write(*,'(A,2(2x,I0))') 'NB,Nelec = ',NB,Nelec
    write(*,'(A,2(2x,I0))') 'Ncg,Nscf = ',Ncg,Nscf
    write(*,'(A,2x,I0,2x,e16.6e3)') 'Nt,dt = ',Nt,dt
    write(*,'(A,2x,A,2x,I6)')'option_TD_band_dos,Nstep_TD_band_dos = ',option_TD_band_dos,Nstep_TD_band_dos
    write(*,'(A,2x,A)')'option_IP = ',option_IP
    write(*,'(A,2x,e16.6e3)') 'dAc = ',dAc
    write(*,'(A,2x,A,2x,A)') 'laser_type, Tr_Lo = ',laser_type,Tr_Lo
    write(*,'(A,2x,4(2x,e16.6e3))') 'IWcm2_1,tpulsefs_1,omegaev_1,phi_CEP_1 = ',IWcm2_1,tpulsefs_1,omegaev_1,phi_CEP_1
    write(*,'(A,2x,3(2x,e16.6e3))') 'Epdir_1(1),Epdir_1(2),Epdir_1(3) = ',Epdir_1(1),Epdir_1(2),Epdir_1(3)
    write(*,'(A,2x,4(2x,e16.6e3))') 'IWcm2_2,tpulsefs_2,omegaev_2,phi_CEP_2 = ',IWcm2_2,tpulsefs_2,omegaev_2,phi_CEP_2
    write(*,'(A,2x,3(2x,e16.6e3))') 'Epdir_1(1),Epdir_1(2),Epdir_1(3) = ',Epdir_2(1),Epdir_2(2),Epdir_2(3)
    write(*,'(A,2x,e16.6e3)') 'T1_T2fs = ',T1_T2fs
    write(*,'(A,2(2x,I0))') 'NI,NE = ',NI,NE

  end if

  call MPI_BCAST(Time_shoutdown,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(SYSname,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(cEex_Cor,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(cVal_mBJ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ps_format,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(method,3,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(type_spatial_difference,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(a_Cvec_d,3*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(aL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(aL1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(aL2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(aL3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NL1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NL2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NL3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NK1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NK2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NK3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nelec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Ncg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nscf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(option_TD_band_dos,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Nstep_TD_band_dos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(option_IP,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Npred_corr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dAc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Tr_Lo,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(laser_type,7,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(IWcm2_1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tpulsefs_1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(omegaev_1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(phi_CEP_1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Epdir_1,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(IWcm2_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tpulsefs_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(omegaev_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(phi_CEP_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Epdir_2,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(T1_T2fs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)    
  call MPI_BCAST(NE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)    
  
  allocate(Zatom(NE),Kion(NI),Rion_Lvec(3,NI))

  if(myrank == 0)then
    read(*,*)(Zatom(j),j=1,NE)
    do i=1,NI
      read(*,*)j,Rion_Lvec(1,i),Rion_Lvec(2,i),Rion_Lvec(3,i),Kion(i)
    end do

    write(*,'(A,999(2x,I0))')'Zatom(:)',(Zatom(j),j=1,NE)
    do i=1,NI
      write(*,'(A,2x,i0,2x,3(e16.6e3,2x),i0)') &
        &'i,Rion_Lvec(1,i),Rion_Lvec(2,i),Rion_Lvec(3,i),Kion(i)',i,Rion_Lvec(1,i),Rion_Lvec(2,i),Rion_Lvec(3,i),Kion(i)
    end do


    write(*,'(A)') '!-----------------------------------------------------------------------------------------------------'
  end if

  call MPI_BCAST(Zatom,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Kion,NI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(Rion_Lvec,3*NI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!== input end


  call init_grid
  call prep_finite_difference
  call NK_split
  call pre_allocation
  call prep_Discrete_Fourier_Transformation
  call prep_Conjugate_Gradient
  call prep_subspace_diag

  select case(method)
  case('PAW')
    err_message='In preparation PAW method' ; call err_finalize
    call err_finalize
  case('PSE')
    call PSE_input_pseudopotential_YS
    call PSE_prep_ps_periodic
  case('CHK')
  case default
    err_message='Invalid option "mehtod"' ; call err_finalize
  end select

  call init_rho_p

  call init_wf
  call Gram_Schmidt

  call psi_rho('GS')
  allocate(rho_e_old(NL)) ! linear mixing
  rho_e_old=rho_e  ! linear mixing
  rho_c=rho_e-rho_p

  call Hartree
  cEex_Cor_tmp = cEex_Cor; cEex_Cor = 'PZ'
  call Exc_Cor;  cEex_Cor = cEex_Cor_tmp
  call local_potential

  if(Myrank == 0)write(*,'(A)')'preparation is completed' 
  return
end subroutine preparation
