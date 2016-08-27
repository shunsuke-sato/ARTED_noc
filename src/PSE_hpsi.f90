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
subroutine PSE_hpsi(ik)
  use global_variables
  use PSE_variables
  implicit none
  integer :: ik


  select case(type_spatial_difference)
  case('FD')
    call PSE_hpsi_grid(ik)
  case('FT')
    call PSE_hpsi_DFT(ik)
  case default
    err_message = 'invalid parameter in type_spatial_difference'
    call err_finalize
  end select


  return
end subroutine PSE_hpsi
