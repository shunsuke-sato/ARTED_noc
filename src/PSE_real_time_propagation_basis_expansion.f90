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
subroutine PSE_real_time_propagation_basis_expansion
  use global_variables
  implicit none
  integer :: iter,iter_t
  real(8) :: jav,Act_t

  call init_wf_basis_expansion
  call init_Ac_basis_expansion

!== Start current
  Act_t = Actot_BE(0)
  kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Act_t*Epdir_1(1)
  kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Act_t*Epdir_1(2)
  kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Act_t*Epdir_1(3)
!  call BE_current(jav)
  javt_BE(0)=jav
!== End current

  do iter=0,Nt
!== Start dt_evolve
    Act_t = 0.5d0*( Actot_BE(iter+1) + Actot_BE(iter) )
    kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Act_t*Epdir_1(1)
    kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Act_t*Epdir_1(2)
    kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Act_t*Epdir_1(3)
!  call BE_dt_evolve
!== End dt_evolve

!== Start current
  Act_t = Actot_BE(iter+1)
  kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Act_t*Epdir_1(1)
  kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Act_t*Epdir_1(2)
  kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Act_t*Epdir_1(3)
!  call BE_current(jav)
  javt_BE(iter+1)=jav
!== End current

  end do

  return
end subroutine PSE_real_time_propagation_basis_expansion
