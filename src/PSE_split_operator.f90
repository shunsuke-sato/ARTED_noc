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
subroutine PSE_split_operator
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
      zcoef=zcoef*(-zI*Dt*0.5d0)/dble(iexp)
      call PSE_Vpsi(ik)
      zu(:,ib,ik)=zu(:,ib,ik)+zcoef*htpsi(:)
      tpsi(:)=htpsi(:)
    end do
  end do Band1
  end do K_points1

! exp(-zI*Dt*\hat{T})  
  K_points2 :do ik=NK_s,NK_e

    zft2(:,:,:)=-0.5d0*Lap_k(:,:,:) &
      -kAc_Cvec(1,ik)*Grad_x_zI(:,:,:) &
      -kAc_Cvec(2,ik)*Grad_y_zI(:,:,:) &
      -kAc_Cvec(3,ik)*Grad_z_zI(:,:,:) &
      +sum(kAc_Cvec(:,ik)**2)*0.5
    
    zft3(:,:,:)=exp(-zI*Dt*zft2(:,:,:))/dble(NL1*NL2*NL3)

  Band2 :    do ib=1,NB_TD
    tpsi(:)=zu(:,ib,ik)

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
    
    zft1(:,:,:)=zft3(:,:,:)*zft1(:,:,:)
  
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
      zu(i,ib,ik)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
    end do



   end do Band2
   end do K_points2


! exp(-zI*Dt/2*\hat{V})  
  K_points3 :do ik=NK_s,NK_e
  Band3 :    do ib=1,NB_TD
    tpsi(:)=zu(:,ib,ik)
    zcoef=1d0
    do iexp=1,4
      zcoef=zcoef*(-zI*Dt*0.5d0)/dble(iexp)
      call PSE_Vpsi(ik)
      zu(:,ib,ik)=zu(:,ib,ik)+zcoef*htpsi(:)
      tpsi(:)=htpsi(:)
    end do
   end do Band3
  end do K_points3

  return
end subroutine PSE_split_operator
