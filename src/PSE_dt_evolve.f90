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
subroutine PSE_dt_evolve(iter)
  use global_variables
  use PSE_variables
  integer :: iter,ipred_corr
  integer :: iexp


  kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Actot_Cvec(1,iter)
  kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Actot_Cvec(2,iter)
  kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Actot_Cvec(3,iter)
  call PSE_pre_hpsi
  zu_t(:,:,:)=zu(:,:,:)
  Vloc_t(:)=Vloc(:)

  do iprep_corr=1,Npred_corr

!    call PSE_split_operator_4th
    call PSE_Taylor
    call psi_rho('RT') 

    rho_c=rho_e-rho_p
    call Hartree
    call Exc_Cor
    call local_potential

    Vloc(:)=0.5d0*(Vloc(:)+Vloc_t(:))
    zu(:,:,:)=zu_t(:,:,:)
    if(iprep_corr == 1)then
      kAc_Cvec(1,:)=kAc0_Cvec(1,:)+0.5d0*(Actot_Cvec(1,iter)+Actot_Cvec(1,iter+1))
      kAc_Cvec(2,:)=kAc0_Cvec(2,:)+0.5d0*(Actot_Cvec(2,iter)+Actot_Cvec(2,iter+1))
      kAc_Cvec(3,:)=kAc0_Cvec(3,:)+0.5d0*(Actot_Cvec(3,iter)+Actot_Cvec(3,iter+1))
      call PSE_pre_hpsi
    end if


  end do

!  call PSE_split_operator_4th
  call PSE_Taylor
  call psi_rho('RT') 
  rho_c=rho_e-rho_p

  if(myrank == 0)write(*,*)'charge =',sum(rho_c)*H123,sum(rho_e)*H123,sum(rho_p)*H123  
    call Hartree
    call Exc_Cor
  call local_potential


  return
end subroutine PSE_dt_evolve
  
