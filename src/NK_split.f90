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
subroutine NK_split
  use global_variables
  implicit none

  NK_ave=NK/Nprocs; NK_remainder=mod(NK,Nprocs)
  if(Myrank < NK_remainder)then
    NK_s=(NK_ave+1)*Myrank+1
    NK_e=(NK_ave+1)*Myrank+(NK_ave+1)
  else if(Myrank >= NK_remainder)then
    NK_s=(Myrank-NK_remainder)*NK_ave+(NK_ave+1)*NK_remainder+1
    NK_e=(Myrank-NK_remainder)*NK_ave+(NK_ave+1)*NK_remainder+NK_ave
  end if

  if(Myrank == 0)write(*,'(A)')'NK-split is completed'
  if(Myrank == 0)write(*,'(2(A,2x,I0,2x))')'NK_ave =',NK_ave,'NK_remainder =',NK_remainder

  return
end subroutine NK_split
