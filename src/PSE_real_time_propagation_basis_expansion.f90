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

  if(myrank == 0)write(*,"(A)")"== Start real-time propagation with basis expansion."

  call init_wf_basis_expansion
  call init_Ac_basis_expansion

!== Start current
  Act_t = Actot_BE(0)
  call BE_current(jav)
  javt_BE(0)=jav
!== End current

  do iter=0,Nt
!== Start dt_evolve
    Act_t = 0.5d0*( Actot_BE(iter+1) + Actot_BE(iter) )
    call BE_dt_evolve(Act_t)
!== End dt_evolve

!== Start current
  Act_t = Actot_BE(iter+1)
  call BE_current(jav)
  javt_BE(iter+1)=jav
!== End current

  end do

  if(myrank == 0)write(*,"(A)")"== End real-time propagation with basis expansion."

  return
end subroutine PSE_real_time_propagation_basis_expansion
