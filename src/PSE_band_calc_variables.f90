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
module PSE_band_calc_variables

  integer :: Nsym_point, NB_band_calc
  real(8),allocatable :: b0_band_vec(:), band_Cvec_d(:,:)
  character(10),allocatable :: Name_sym_point(:)
  integer,allocatable :: Ndist_band(:)

end module PSE_band_calc_variables
