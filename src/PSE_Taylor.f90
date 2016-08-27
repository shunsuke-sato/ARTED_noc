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
subroutine PSE_Taylor
  use global_variables
  use PSE_variables
  integer :: ik,ib,iexp
  complex(8) :: zcoef


! exp(-zI*Dt/2*\hat{V})  
  K_points1 :do ik=NK_s,NK_e
  Band1 :    do ib=1,NB_TD
    tpsi(:)=zu(:,ib,ik)
    zcoef=1d0
    do iexp=1,4
      zcoef=zcoef*(-zI*Dt)/dble(iexp)
      call PSE_hpsi(ik)
      zu(:,ib,ik)=zu(:,ib,ik)+zcoef*htpsi(:)
      tpsi(:)=htpsi(:)
    end do
  end do Band1
  end do K_points1

  return
end subroutine PSE_Taylor
