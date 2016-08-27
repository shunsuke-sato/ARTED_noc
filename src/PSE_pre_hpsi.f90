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
subroutine PSE_pre_hpsi
  use global_variables
  use PSE_variables
  implicit none
  integer :: ik
  integer :: ia,j,i,ix,iy,iz
  real(8) :: kr,x,y,z,r
  real(8) :: dist_Lvec(3)
  
  do ik=NK_s,NK_e
    do ia=1,NI
      do j=1,Mps(ia)
        i=Jxyz(j,ia)
        dist_Lvec(1)=Lx(1,i)-dble(Jxx(j,ia))
        dist_Lvec(2)=Lx(2,i)-dble(Jyy(j,ia))
        dist_Lvec(3)=Lx(3,i)-dble(Jzz(j,ia))

        x=sum(A_matrix(1,:)*dist_Lvec(:))
        y=sum(A_matrix(2,:)*dist_Lvec(:))
        z=sum(A_matrix(3,:)*dist_Lvec(:))
        r=sqrt(x**2+y**2+z**2)+1d-50

        kr=kAc_Cvec(1,ik)*x+kAc_Cvec(2,ik)*y+kAc_Cvec(3,ik)*z
        ekr(j,ia,ik)=exp(zI*kr)
      enddo
    enddo
  end do


  return
end subroutine PSE_pre_hpsi
