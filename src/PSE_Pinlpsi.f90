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
subroutine PSE_Pinlpsi(ik)
  use global_variables
  use PSE_variables
  implicit none
  integer :: ik
  integer :: i,ix,iy,iz
  integer :: ilma,ia,j
  complex(8) :: uVpsi

  htpsi= 0d0 !htpsi+(Vloc+kAc2_2)*tpsi

!  return
!Calculating nonlocal part
  do ilma=1,Nlma
    ia=a_tbl(ilma)
    uVpsi=0.d0
    do j=1,Mps(ia)
      i=Jxyz(j,ia)
      uVpsi=uVpsi+uV(j,ilma)*ekr(j,ia,ik)*tpsi(i)
    enddo
    uVpsi=uVpsi*H123*iuV(ilma)
    do j=1,Mps(ia)
      i=Jxyz(j,ia)
      htpsi(i)=htpsi(i)+conjg(ekr(j,ia,ik))*uVpsi*uV(j,ilma)
    enddo
  enddo

  return
end subroutine PSE_Pinlpsi
