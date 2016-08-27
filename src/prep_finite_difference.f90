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
subroutine prep_finite_difference
  use global_variables
  implicit none
  integer,parameter :: Nd=4
  integer :: i,j,k,n
  real(8) :: t,s,s0

  allocate(lap(-Nd:Nd),nab(-Nd:Nd))
  allocate(lap1(-Nd:Nd),lap2(-Nd:Nd),lap3(-Nd:Nd))
  allocate(nab1(-Nd:Nd),nab2(-Nd:Nd),nab3(-Nd:Nd))
  allocate(ifdx1(-Nd:Nd,1:NL),ifdx2(-Nd:Nd,1:NL),ifdx3(-Nd:Nd,1:NL))

  lap=0.d0
  nab=0.d0

  n=Nd
  s=0.d0
  do i=-n,n
    if ( i==0 ) cycle
    do j=-n,n
      if ( j==i .or. j==0 ) cycle
      s=s+1.d0/dble(i*j)
    end do
  end do
  lap(0)=s

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=0.d0
    do k=-n,n
      if ( k==j .or. k==0 ) cycle
      s0=1.d0
      do i=-n,n
        if ( i==k .or. i==j .or. i==0 ) cycle
        s0=s0*(-i)
      end do
      s=s+s0
    end do
    lap( j)=2.d0*s/t
    lap(-j)=lap(j)
  end do

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=1.d0
    do i=-n,n
      if ( i==j .or. i==0 ) cycle
      s=s*(-i)
    end do
    nab( j)=s/t
    nab(-j)=-nab(j)
  end do

  lap1(-Nd:Nd)=lap(-Nd:Nd)/H1**2
  lap2(-Nd:Nd)=lap(-Nd:Nd)/H2**2
  lap3(-Nd:Nd)=lap(-Nd:Nd)/H3**2
  nab1(-Nd:Nd)=nab(-Nd:Nd)/H1
  nab2(-Nd:Nd)=nab(-Nd:Nd)/H2
  nab3(-Nd:Nd)=nab(-Nd:Nd)/H3

  do i=1,NL
    do n=-Nd,Nd
      ifdx1(n,i)=iLx123(mod(iLx(1,i)+n+NL1,NL1),iLx(2,i),iLx(3,i))
      ifdx2(n,i)=iLx123(iLx(1,i),mod(iLx(2,i)+n+NL2,NL2),iLx(3,i))
      ifdx3(n,i)=iLx123(iLx(1,i),iLx(2,i),mod(iLx(3,i)+n+NL3,NL3))
    end do
  end do



  return
end subroutine prep_finite_difference
