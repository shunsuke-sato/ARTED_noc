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
subroutine PSE_ppsi_DFT(ik)
  use global_variables
  use PSE_variables
  implicit none
  integer :: ik
  integer :: i,ix,iy,iz
  real(8) :: kAc2_2
  integer :: ilma,ia,j
  complex(8) :: uVpsi

  kAc2_2=sum(kAc_Cvec(:,ik)**2)*0.5

  do i=1,NL
    zft2(iLx(1,i),iLx(2,i),iLx(3,i))=tpsi(i)
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(cexp_x3(:,iz)*zft2(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(cexp_x2(:,iy)*zft1(ix,:,iz))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(cexp_x1(:,ix)*zft2(:,iy,iz))
  end do
  end do
  end do

  zft1(:,:,:)=( &
    -Epdir_1(1)*Grad_x_zI(:,:,:) &
    -Epdir_1(2)*Grad_y_zI(:,:,:) &
    -Epdir_1(3)*Grad_z_zI(:,:,:) &
    )*zft1(:,:,:)/dble(NL1*NL2*NL3)

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*zft1(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    htpsi(i)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do

  htpsi=htpsi+( & 
    +Epdir_1(1)*kAc_Cvec(1,ik) &
    +Epdir_1(2)*kAc_Cvec(2,ik) &
    +Epdir_1(3)*kAc_Cvec(3,ik) &
    )*tpsi
  return

  return
end subroutine PSE_ppsi_DFT
