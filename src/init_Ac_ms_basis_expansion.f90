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
subroutine init_Ac_ms_basis_expansion
  use global_variables
  use ms_maxwell_ks_variables
  implicit none
  integer :: iter
  real(8) :: tt
  real(8) :: f0_1,f0_2,omega_1,omega_2,tpulse_1,tpulse_2,T1_T2

  if(myrank == 0)write(*,"(A)")"== Start: Initialization of vector potential."
  allocate(Ac_m(nx_s:nx_e), jt_m(Mx))
  allocate(Ac_m_n(nx_s:nx_e),Ac_m_o(nx_s:nx_e))

  if(myrank == 0)write(*,"(A)")"== End: Initialization of vector potential."

  return
end subroutine init_Ac_ms_basis_expansion
