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
  use PSE_variables
  implicit none
  integer :: iter,iter_t
  real(8) :: jav

  call init_Ac_basis_expansion

  return
end subroutine PSE_real_time_propagation_basis_expansion
