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
subroutine PSE_hpsi_grid(ik)
  use global_variables
  use PSE_variables
  implicit none
  integer :: ik
  integer :: i,ix,iy,iz
  real(8) :: kAc2_2
  integer :: ilma,ia,j
  complex(8) :: uVpsi
  complex(8) :: zs(3),zt(3)
  real(8) :: b11,b22,b33,b12,b23,b31
  real(8) :: bb1,bb2,bb3

  kAc2_2=sum(kAc_Cvec(:,ik)**2)*0.5d0
  b11=B_t_matrix(1,1)
  b22=B_t_matrix(1,1)
  b33=B_t_matrix(1,1)
  b12=B_t_matrix(1,2)
  b23=B_t_matrix(2,3)
  b31=B_t_matrix(3,1)
  bb1=kAc_Cvec(1,ik)*B_matrix(1,1)+kAc_Cvec(2,ik)*B_matrix(1,2)+kAc_Cvec(3,ik)*B_matrix(1,3)
  bb2=kAc_Cvec(1,ik)*B_matrix(2,1)+kAc_Cvec(2,ik)*B_matrix(2,2)+kAc_Cvec(3,ik)*B_matrix(2,3)
  bb3=kAc_Cvec(1,ik)*B_matrix(3,1)+kAc_Cvec(2,ik)*B_matrix(3,2)+kAc_Cvec(3,ik)*B_matrix(3,3)


    do i=1,NL
      zs(1)=lap1(0)*tpsi(i)&
          &+lap1(1)*(tpsi(ifdx1(1,i))+tpsi(ifdx1(-1,i)))&
          &+lap1(2)*(tpsi(ifdx1(2,i))+tpsi(ifdx1(-2,i)))&
          &+lap1(3)*(tpsi(ifdx1(3,i))+tpsi(ifdx1(-3,i)))&
          &+lap1(4)*(tpsi(ifdx1(4,i))+tpsi(ifdx1(-4,i)))

      zt(1)= &!nab1(0)*tpsi(i)&
          &+nab1(1)*(tpsi(ifdx1(1,i))-tpsi(ifdx1(-1,i)))&
          &+nab1(2)*(tpsi(ifdx1(2,i))-tpsi(ifdx1(-2,i)))&
          &+nab1(3)*(tpsi(ifdx1(3,i))-tpsi(ifdx1(-3,i)))&
          &+nab1(4)*(tpsi(ifdx1(4,i))-tpsi(ifdx1(-4,i)))

      zs(2)=lap2(0)*tpsi(i)&
          &+lap2(1)*(tpsi(ifdx2(1,i))+tpsi(ifdx2(-1,i)))&
          &+lap2(2)*(tpsi(ifdx2(2,i))+tpsi(ifdx2(-2,i)))&
          &+lap2(3)*(tpsi(ifdx2(3,i))+tpsi(ifdx2(-3,i)))&
          &+lap2(4)*(tpsi(ifdx2(4,i))+tpsi(ifdx2(-4,i)))

      zt(2)= &!nab2(0)*tpsi(i)&
          &+nab2(1)*(tpsi(ifdx2(1,i))-tpsi(ifdx2(-1,i)))&
          &+nab2(2)*(tpsi(ifdx2(2,i))-tpsi(ifdx2(-2,i)))&
          &+nab2(3)*(tpsi(ifdx2(3,i))-tpsi(ifdx2(-3,i)))&
          &+nab2(4)*(tpsi(ifdx2(4,i))-tpsi(ifdx2(-4,i)))

      zs(3)=lap3(0)*tpsi(i)&
          &+lap3(1)*(tpsi(ifdx3(1,i))+tpsi(ifdx3(-1,i)))&
          &+lap3(2)*(tpsi(ifdx3(2,i))+tpsi(ifdx3(-2,i)))&
          &+lap3(3)*(tpsi(ifdx3(3,i))+tpsi(ifdx3(-3,i)))&
          &+lap3(4)*(tpsi(ifdx3(4,i))+tpsi(ifdx3(-4,i)))

      zt(3)= &!nab3(0)*tpsi(i)&
          &+nab3(1)*(tpsi(ifdx3(1,i))-tpsi(ifdx3(-1,i)))&
          &+nab3(2)*(tpsi(ifdx3(2,i))-tpsi(ifdx3(-2,i)))&
          &+nab3(3)*(tpsi(ifdx3(3,i))-tpsi(ifdx3(-3,i)))&
          &+nab3(4)*(tpsi(ifdx3(4,i))-tpsi(ifdx3(-4,i)))

      htpsi(i)=-0.5d0*(b11*zs(1)+b22*zs(2)+b33*zs(3))-zI*(bb1*zt(1)+bb2*zt(2)+bb3*zt(3))
      tpsi_g1(i)=zt(1)
      tpsi_g2(i)=zt(2)
    enddo

    do i=1,NL
!(1,2) >> d/dx2 * d/dx1
      zt(1)= &!nab1(0)*tpsi(i)&
          &+nab2(1)*(tpsi_g1(ifdx2(1,i))-tpsi_g1(ifdx2(-1,i)))&
          &+nab2(2)*(tpsi_g1(ifdx2(2,i))-tpsi_g1(ifdx2(-2,i)))&
          &+nab2(3)*(tpsi_g1(ifdx2(3,i))-tpsi_g1(ifdx2(-3,i)))&
          &+nab2(4)*(tpsi_g1(ifdx2(4,i))-tpsi_g1(ifdx2(-4,i)))

!(2,3) >> d/dx3 * d/dx2

      zt(2)= &!nab3(0)*tpsi(i)&
          &+nab3(1)*(tpsi_g2(ifdx3(1,i))-tpsi_g2(ifdx3(-1,i)))&
          &+nab3(2)*(tpsi_g2(ifdx3(2,i))-tpsi_g2(ifdx3(-2,i)))&
          &+nab3(3)*(tpsi_g2(ifdx3(3,i))-tpsi_g2(ifdx3(-3,i)))&
          &+nab3(4)*(tpsi_g2(ifdx3(4,i))-tpsi_g2(ifdx3(-4,i)))

!(3,1) >> d/dx3 * d/dx1
      zt(3)= &!nab1(0)*tpsi(i)&
          &+nab3(1)*(tpsi_g1(ifdx3(1,i))-tpsi_g1(ifdx3(-1,i)))&
          &+nab3(2)*(tpsi_g1(ifdx3(2,i))-tpsi_g1(ifdx3(-2,i)))&
          &+nab3(3)*(tpsi_g1(ifdx3(3,i))-tpsi_g1(ifdx3(-3,i)))&
          &+nab3(4)*(tpsi_g1(ifdx3(4,i))-tpsi_g1(ifdx3(-4,i)))

      htpsi(i)=htpsi(i)-0.5d0*2d0*(b12*zt(1)+b23*zt(2)+b31*zt(3))

    end do


    htpsi=htpsi+(Vloc+kAc2_2)*tpsi

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
end subroutine PSE_hpsi_grid
