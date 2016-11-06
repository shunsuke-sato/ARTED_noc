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
subroutine PSE_ground_state_calculation
  use global_variables
  implicit none
  integer :: iter_scf,ib
  real(8) :: rmix=0.2d0
  real(8) :: jav(3)
  real(8) :: Etot,Ekin
  real(8),allocatable :: Vloc_old(:)
  character(10) :: cEex_Cor_tmp

  if(Myrank == 0)write(*,'(A)')'Ground state calculation start'
  if(myrank == 0)write(*,*)'charge =',sum(rho_c)*H123,sum(rho_e)*H123,sum(rho_p)*H123  

  if(Myrank == 0)open(21,file='density_conv_GS.out')

  allocate(rho_MB_in(NL,MaxMem_MB),rho_MB_out(NL,MaxMem_MB)) ! For Broiden's method
  allocate(Vloc_old(NL)); Vloc_old = Vloc

  do iter_scf=1,Nscf
    if(myrank == 0)write(*,*)kAc_Cvec(:,1)
    if(Myrank == 0)write(*,'(A,2X,I5)')'iter_scf = ',iter_scf

    rho_MB_in(:,min(iter_scf,MaxMem_MB)) = rho_e(:) ! For Broiden's method
    
    call Gram_Schmidt
    call PSE_Conjugate_Gradient(Ncg)
    call Gram_Schmidt
    call PSE_subspace_diag
    call Gram_Schmidt

    call psi_rho('GS')
!    rho_e=rho_e*rmix+rho_e_old*(1d0-rmix) ! linear mixing 
!    rho_e_old=rho_e ! linear mixing
    rho_c=rho_e-rho_p
!    rho_MB_out(:,min(iter_scf,MaxMem_MB)) = rho_e(:) ! For Broiden's method
!    if(myrank == 0)write(21,'(I7,2x,999e26.16e3)')iter_scf, &
!      & sum((rho_MB_out(:,min(iter_scf,MaxMem_MB)) - rho_MB_in(:,min(iter_scf,MaxMem_MB)))**2)*H123
!    call modified_Broyden_mixing(iter_scf) ! For Broiden's method

    call Hartree
    cEex_Cor_tmp = cEex_Cor
    if(cEex_Cor_tmp == 'TBmBJ'.and. iter_scf < 20)cEex_Cor = 'PZ'
    call Exc_Cor
    cEex_Cor = cEex_Cor_tmp
    call local_potential
    Vloc = rmix*Vloc + (1d0-rmix)*Vloc_old
    Vloc_old = Vloc
    call PSE_current_GS(jav)

    call PSE_energy(Etot,Ekin,'GS')
    if(myrank == 0)write(*,'(A,e26.16e3)')'Etot =',Etot
    if(myrank == 0)write(*,'(A,e26.16e3)')'Ekin =',Ekin
    if(myrank == 0)write(*,'(3(A6,e16.6e3))')'jav(1)= ',jav(1),'jav(2)= ',jav(2),'jav(3)= ',jav(3)

!    Vloc=Vxc ! debug
    if(myrank == 0)write(*,*)'charge =',sum(rho_c)*H123,sum(rho_e)*H123,sum(rho_p)*H123
    if(myrank == 0)write(*,*)sum(Vh**2)*H123

    if(Myrank == 0)then
      write(*,'(4(I3,e16.6e3))')(ib,esp(ib,NK_s),ib=1,NB)
    end if

    
  end do

  if(Myrank == 0)close(21)

  deallocate(rho_MB_in,rho_MB_out) ! For Broiden's method

  if(myrank == 0)write(*,'(A)')'GS final'
  call psi_rho('GS')
  rho_e=rho_e
  rho_c=rho_e

  call Hartree
  call Exc_Cor
  call local_potential
  call PSE_current_GS(jav)

  call PSE_energy(Etot,Ekin,'GS')
  if(myrank == 0)write(*,'(A,e26.16e3)')'Etot =',Etot
  if(myrank == 0)write(*,'(A,e26.16e3)')'Ekin =',Ekin
  if(myrank == 0)write(*,'(3(A6,e16.6e3))')'jav(1)= ',jav(1),'jav(2)= ',jav(2),'jav(3)= ',jav(3)
  
  Etot_GS=Etot
  Ekin_GS=Ekin
  if(Myrank == 0)write(*,'(A)')'Ground state calculation is completed'
  call Band_DoS

  return
end subroutine PSE_ground_state_calculation
