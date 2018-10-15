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
subroutine MS_preparation
  use global_variables
  use ms_maxwell_ks_variables
  implicit none
  integer :: ix

  write(*,*)myrank
! Parameter for maxwell-kohn-sham
  Nx_s = -2000
  Nx_e =  2000
  Mx = 4
  dx_m = 40d-9/0.529177d-10

  if(mod(Nprocs,Mx) /=0)stop 'Error Mx is not dividable by Nprocs'

  nprocs_per_Mpoint = Nprocs/Mx
  do ix = 1, Mx
    if(nprocs_per_Mpoint*ix > myrank)then
      macro_point_id = ix
      exit
    end if
  end do

  write(*,*)"myrak, macro_point_id",myrank,macro_point_id

end subroutine MS_preparation
