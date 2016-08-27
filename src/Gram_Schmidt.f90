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
subroutine Gram_Schmidt
  use global_variables
  implicit none
  complex(8) :: zs
  real(8) :: s
  integer :: ik,ib,ib2
  do ik=NK_s,NK_e
    do ib=1,NB

      do ib2=1,ib-1
        zs=sum(conjg(zu_GS(:,ib2,ik))*zu_GS(:,ib,ik))*H123
        zu_GS(:,ib,ik)=zu_GS(:,ib,ik)-zs*zu_GS(:,ib2,ik)
      end do

      s=1d0/sqrt(sum(abs(zu_GS(:,ib,ik))**2)*H123)
      zu_GS(:,ib,ik)=zu_GS(:,ib,ik)*s
    end do
  end do


  return
end Subroutine Gram_Schmidt
