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
subroutine prep_Discrete_Fourier_Transformation
  use global_variables
  implicit none
  integer :: i,j,k
  real(8) :: an1,an2,an3
  real(8) :: k1(0:NL1-1),k2(0:NL2-1),k3(0:NL3-1),Bt(3,3)
  real(8) :: kx,ky,kz ! debug
  real(8) :: k1_t,k2_t,k3_t

  an1=norm_a_Cvec(1)
  an2=norm_a_Cvec(2)
  an3=norm_a_Cvec(3)

  allocate(exp_x1(0:NL1-1,0:NL1-1),exp_x2(0:NL2-1,0:NL2-1),exp_x3(0:NL3-1,0:NL3-1))
  allocate(cexp_x1(0:NL1-1,0:NL1-1),cexp_x2(0:NL2-1,0:NL2-1),cexp_x3(0:NL3-1,0:NL3-1))
  allocate(zft1(0:NL1-1,0:NL2-1,0:NL3-1),zft2(0:NL1-1,0:NL2-1,0:NL3-1),zft3(0:NL1-1,0:NL2-1,0:NL3-1))
  allocate(Lap_k(0:NL1-1,0:NL2-1,0:NL3-1),InLap_k(0:NL1-1,0:NL2-1,0:NL3-1))
  allocate(Grad_x_zI(0:NL1-1,0:NL2-1,0:NL3-1),Grad_y_zI(0:NL1-1,0:NL2-1,0:NL3-1),Grad_z_zI(0:NL1-1,0:NL2-1,0:NL3-1))




  do i=0,NL1-1
    do j=0,NL1-1
      exp_x1(i,j)=exp(zI*2d0*pi*dble(i*j)/dble(NL1))
      cexp_x1(i,j)=exp(-zI*2d0*pi*dble(i*j)/dble(NL1))
    end do
  end do
  
  do i=0,NL2-1
    do j=0,NL2-1
      exp_x2(i,j)=exp(zI*2d0*pi*dble(i*j)/dble(NL2))
      cexp_x2(i,j)=exp(-zI*2d0*pi*dble(i*j)/dble(NL2))
    end do
  end do
  
  do i=0,NL3-1
    do j=0,NL3-1
      exp_x3(i,j)=exp(zI*2d0*pi*dble(i*j)/dble(NL3))
      cexp_x3(i,j)=exp(-zI*2d0*pi*dble(i*j)/dble(NL3))
    end do
  end do


  do i=0,NL1-1
    k1(i)=2d0*pi*dble(i)
    if(i >= NL1/2) k1(i)=2d0*pi*dble(i-NL1)
  end do
  do i=0,NL2-1
    k2(i)=2d0*pi*dble(i)
    if(i >= NL2/2) k2(i)=2d0*pi*dble(i-NL2)
  end do
  do i=0,NL3-1
    k3(i)=2d0*pi*dble(i)
    if(i >= NL3/2) k3(i)=2d0*pi*dble(i-NL3)
  end do
  
  
  do i=1,3
  do j=1,3
    Bt(i,j)=sum(B_matrix(i,:)*B_matrix(j,:))
  end do
  end do

  if(myrank == 0)then
    write(*,*)A_matrix(1,1),A_matrix(1,2),A_matrix(1,3)
    write(*,*)A_matrix(2,1),A_matrix(2,2),A_matrix(2,3)
    write(*,*)A_matrix(3,1),A_matrix(3,2),A_matrix(3,3)
    write(*,*)B_matrix(1,1),B_matrix(1,2),B_matrix(1,3)
    write(*,*)B_matrix(2,1),B_matrix(2,2),B_matrix(2,3)
    write(*,*)B_matrix(3,1),B_matrix(3,2),B_matrix(3,3)
    write(*,*)Bt(1,1),Bt(1,2),Bt(1,3)
    write(*,*)Bt(2,1),Bt(2,2),Bt(2,3)
    write(*,*)Bt(3,1),Bt(3,2),Bt(3,3)
  end if

  an1=0d0 ! debug
  InLap_k=0d0
  do i=0,NL1-1
  do j=0,NL2-1
  do k=0,NL3-1

    k1_t = k1(i)
    k2_t = k2(j)
    k3_t = k3(k)
!    if(i==NL1/2)k1_t = 0d0
!    if(j==NL2/2)k2_t = 0d0
!    if(k==NL3/2)k3_t = 0d0

    Lap_k(i,j,k)=-Bt(1,1)*k1(i)*k1(i)-Bt(1,2)*k1_t*k2_t-Bt(1,3)*k1_t*k3_t &
                 -Bt(2,1)*k2_t*k1_t-Bt(2,2)*k2(j)*k2(j)-Bt(2,3)*k2_t*k3_t &
                 -Bt(3,1)*k3_t*k1_t-Bt(3,2)*k3_t*k2_t-Bt(3,3)*k3(k)*k3(k) 

    Grad_x_zI(i,j,k)=-B_matrix(1,1)*k1_t-B_matrix(2,1)*k2_t-B_matrix(3,1)*k3_t
    Grad_y_zI(i,j,k)=-B_matrix(1,2)*k1_t-B_matrix(2,2)*k2_t-B_matrix(3,2)*k3_t
    Grad_z_zI(i,j,k)=-B_matrix(1,3)*k1_t-B_matrix(2,3)*k2_t-B_matrix(3,3)*k3_t

!== debug
    if(myrank == 0)then
      kx=k1(i)*B_matrix(1,1)+k2(j)*B_matrix(2,1)+k3(k)*B_matrix(3,1)
      ky=k1(i)*B_matrix(1,2)+k2(j)*B_matrix(2,2)+k3(k)*B_matrix(3,2)
      kz=k1(i)*B_matrix(1,3)+k2(j)*B_matrix(2,3)+k3(k)*B_matrix(3,3)
!      write(*,'(3(I4,2x),3e16.6e3)')i,j,k,kx+Grad_x_zI(i,j,k),ky+Grad_y_zI(i,j,k),kz+Grad_z_zI(i,j,k)
      an1=an1+abs(kx+Grad_x_zI(i,j,k))+abs(ky+Grad_y_zI(i,j,k))+abs(kz+Grad_z_zI(i,j,k))
    end if
!== debug
    
!    if(i+j+k == 0)cycle
    if(abs(Lap_k(i,j,k)) < 1d-10)then
      if(myrank == 0)write(*,*)'i,j,k,Lap',i,j,k,Lap_k(i,j,k)
      InLap_k(i,j,k)=0d0
    else
      InLap_k(i,j,k)=1d0/Lap_k(i,j,k)
    end if
  end do
  end do
  end do

  if(myrank == 0)write(*,*)'an1=',an1

  return
end subroutine prep_Discrete_Fourier_Transformation
