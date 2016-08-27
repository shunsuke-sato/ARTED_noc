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
subroutine init_grid
  use global_variables
  implicit none
  integer :: i,j
  real(8) :: tvec(3)
  integer :: i1,i2,i3

  NL=NL1*NL2*NL3
  NK=NK1*NK2*NK3

  a_Cvec(:,1)=aL*aL1*a_Cvec_d(:,1)
  a_Cvec(:,2)=aL*aL2*a_Cvec_d(:,2)
  a_Cvec(:,3)=aL*aL3*a_Cvec_d(:,3)

  tvec(1)=a_Cvec(2,2)*a_Cvec(3,3)-a_Cvec(3,2)*a_Cvec(2,3)
  tvec(2)=a_Cvec(3,2)*a_Cvec(1,3)-a_Cvec(1,2)*a_Cvec(3,3)
  tvec(3)=a_Cvec(1,2)*a_Cvec(2,3)-a_Cvec(2,2)*a_Cvec(1,3)

  norm_a_Cvec(1)=sqrt(sum(a_Cvec(:,1)**2))
  norm_a_Cvec(2)=sqrt(sum(a_Cvec(:,2)**2))
  norm_a_Cvec(3)=sqrt(sum(a_Cvec(:,3)**2))

  do i=1,3
    do j=1,3
      mat_vv_a_Cvec(i,j)=sum(a_Cvec(:,i)*a_Cvec(:,j))
    end do
  end do


  call AB_matrix


  H1=1d0/dble(NL1) 
  H2=1d0/dble(NL2) 
  H3=1d0/dble(NL3) 

  Vcell=abs(sum(a_Cvec(:,1)*tvec(:)))
  H123=Vcell/dble(NL)

  b_Cvec(1,1)=2d0*pi*(a_Cvec(2,2)*a_Cvec(3,3)-a_Cvec(3,2)*a_Cvec(2,3))
  b_Cvec(2,1)=2d0*pi*(a_Cvec(3,2)*a_Cvec(1,3)-a_Cvec(1,2)*a_Cvec(3,3))
  b_Cvec(3,1)=2d0*pi*(a_Cvec(1,2)*a_Cvec(2,3)-a_Cvec(2,2)*a_Cvec(1,3))

  b_Cvec(1,2)=2d0*pi*(a_Cvec(2,3)*a_Cvec(3,1)-a_Cvec(3,3)*a_Cvec(2,1))
  b_Cvec(2,2)=2d0*pi*(a_Cvec(3,3)*a_Cvec(1,1)-a_Cvec(1,3)*a_Cvec(3,1))
  b_Cvec(3,2)=2d0*pi*(a_Cvec(1,3)*a_Cvec(2,1)-a_Cvec(2,3)*a_Cvec(1,1))

  b_Cvec(1,3)=2d0*pi*(a_Cvec(2,1)*a_Cvec(3,2)-a_Cvec(3,1)*a_Cvec(2,2))
  b_Cvec(2,3)=2d0*pi*(a_Cvec(3,1)*a_Cvec(1,2)-a_Cvec(1,1)*a_Cvec(3,2))
  b_Cvec(3,3)=2d0*pi*(a_Cvec(1,1)*a_Cvec(2,2)-a_Cvec(2,1)*a_Cvec(1,2))

  b_Cvec=b_Cvec/sum(a_Cvec(:,1)*tvec(:))

  if(myrank == 0)then
    write(*,*)'check axb'
    write(*,*)sum(a_Cvec(:,1)*b_Cvec(:,1))/(2d0*pi),sum(a_Cvec(:,1)*b_Cvec(:,2))/(2d0*pi),sum(a_Cvec(:,1)*b_Cvec(:,3))/(2d0*pi)
    write(*,*)sum(a_Cvec(:,2)*b_Cvec(:,1))/(2d0*pi),sum(a_Cvec(:,2)*b_Cvec(:,2))/(2d0*pi),sum(a_Cvec(:,2)*b_Cvec(:,3))/(2d0*pi)
    write(*,*)sum(a_Cvec(:,3)*b_Cvec(:,1))/(2d0*pi),sum(a_Cvec(:,3)*b_Cvec(:,2))/(2d0*pi),sum(a_Cvec(:,3)*b_Cvec(:,3))/(2d0*pi)
  end if
    

  dk1=1d0/dble(NK1) 
  dk2=1d0/dble(NK2) 
  dk3=1d0/dble(NK3) 

  norm_b_Cvec(1)=sqrt(sum(b_Cvec(:,1)**2))
  norm_b_Cvec(2)=sqrt(sum(b_Cvec(:,2)**2))
  norm_b_Cvec(3)=sqrt(sum(b_Cvec(:,3)**2))

  allocate(kAc_Rvec(3,NK),kAc0_Rvec(3,NK),Lx(3,NL),iLx(3,NL),iLx123(0:NL1-1,0:NL2-1,0:NL3-1))
  allocate(kAc_Cvec(3,NK),kAc0_Cvec(3,NK))

  i=0
  do i1=0,NL1-1
    do i2=0,NL2-1
      do i3=0,NL3-1
        i=i+1
        iLx(1,i)=i1; iLx(2,i)=i2; iLx(3,i)=i3
        iLx123(i1,i2,i3)=i
      end do
    end do
  end do
  Lx(1,:)=dble(iLx(1,:))*H1; Lx(2,:)=dble(iLx(2,:))*H2; Lx(3,:)=dble(iLx(3,:))*H3


!  do i=1,NL
!    Rxyz_Cvec(1,i)=Lx(1,i)*a_Cvec(1,1) &
!            &+Lx(2,i)*a_Cvec(1,2) &
!            &+Lx(3,i)*a_Cvec(1,3)
!    Rxyz_Cvec(2,i)=Lx(1,i)*a_Cvec(2,1) &
!            &+Lx(2,i)*a_Cvec(2,2) &
!            &+Lx(3,i)*a_Cvec(2,3) 
!    Rxyz_Cvec(3,i)=Lx(1,i)*a_Cvec(3,1) &
!            &+Lx(2,i)*a_Cvec(3,2) &
!            &+Lx(3,i)*a_Cvec(3,3)
!  end do

  
  i=0
  do i1=0,NK1-1
    do i2=0,NK2-1
      do i3=0,NK3-1
        i=i+1
        kAc0_Rvec(1,i)=dble(i1)*dk1-dk1*(dble(NK1/2)-0.5d0)
        kAc0_Rvec(2,i)=dble(i2)*dk2-dk2*(dble(NK2/2)-0.5d0)
        kAc0_Rvec(3,i)=dble(i3)*dk3-dk3*(dble(NK3/2)-0.5d0)
      end do
    end do
  end do

  kAc0_Cvec(1,:)=b_Cvec(1,1)*kAc0_Rvec(1,:)+b_Cvec(1,2)*kAc0_Rvec(2,:)+b_Cvec(1,3)*kAc0_Rvec(3,:)
  kAc0_Cvec(2,:)=b_Cvec(2,1)*kAc0_Rvec(1,:)+b_Cvec(2,2)*kAc0_Rvec(2,:)+b_Cvec(2,3)*kAc0_Rvec(3,:)
  kAc0_Cvec(3,:)=b_Cvec(3,1)*kAc0_Rvec(1,:)+b_Cvec(3,2)*kAc0_Rvec(2,:)+b_Cvec(3,3)*kAc0_Rvec(3,:)
  if(myrank == 0)write(*,*)'kAc0',kAc0_Cvec(:,1)
  if(myrank == 0)write(*,*)'b',b_Cvec(1,:)
  if(myrank == 0)write(*,*)'b',b_Cvec(2,:)
  if(myrank == 0)write(*,*)'b',b_Cvec(3,:)
  kAc_Rvec=kAc0_Rvec
  kAc_Cvec=kAc0_Cvec

  if(myrank == 0)then
    open(10,file='kac.out')
    do i=1,NK
      write(10,'(100e26.16e3)')kAc_Cvec(1,i),kAc_Cvec(2,i),kAc_Cvec(3,i)
    end do
    close(10)
  end if
  

  if(Myrank == 0)write(*,'(A)')'Real and k-space grid initialization is completed'

  return
end subroutine init_grid
