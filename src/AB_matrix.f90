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
subroutine AB_matrix
  use global_variables
  implicit none
  integer :: i,j
  real(8) :: detA
  real(8) :: a(3,3),b(3,3)

  A_matrix(:,:)=a_Cvec(:,:)
  a=A_matrix

  detA=a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3) &
    -a(1,3)*a(2,2)*a(3,1)-a(2,3)*a(3,2)*a(1,1)-a(3,3)*a(1,2)*a(2,1)

  b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
  b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
  b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)

  b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
  b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
  b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)

  b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)

  B_matrix=b/detA

  if(myrank == 0)then
    write(*,*)'AB check'
    write(*,*)sum(A_matrix(1,:)*B_matrix(:,1)),sum(A_matrix(1,:)*B_matrix(:,2)),sum(A_matrix(1,:)*B_matrix(:,3))
    write(*,*)sum(A_matrix(2,:)*B_matrix(:,1)),sum(A_matrix(2,:)*B_matrix(:,2)),sum(A_matrix(2,:)*B_matrix(:,3))
    write(*,*)sum(A_matrix(3,:)*B_matrix(:,1)),sum(A_matrix(3,:)*B_matrix(:,2)),sum(A_matrix(3,:)*B_matrix(:,3))
  end if

  do i=1,3
  do j=1,3
    B_t_matrix(i,j)=sum(B_matrix(i,:)*B_matrix(j,:))
  end do
  end do

  return
end subroutine AB_matrix
  
