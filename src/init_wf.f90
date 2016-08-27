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
Subroutine init_wf
  use Global_Variables
  implicit none
  integer :: iseed,ib,ik,i
  real(8) :: r2,x1,y1,z1,rnd
    real(8) :: x0,y0,z0
  real(8) :: r22,x2,y2,z2

  zu_GS=0.d0
  iseed=123
  do ik=1,NK
    do ib=1,NB
      call quickrnd(iseed,rnd)
      x0=rnd
      call quickrnd(iseed,rnd)
      y0=rnd
      call quickrnd(iseed,rnd)
      z0=rnd
      if(ik >= NK_s .and. ik<=NK_e)then
        do i=1,NL
          x1=Lx(1,i)-x0
          y1=Lx(2,i)-y0
          z1=Lx(3,i)-z0

          x2=A_matrix(1,1)*x1+A_matrix(1,2)*y1+A_matrix(1,3)*z1
          y2=A_matrix(2,1)*x1+A_matrix(2,2)*y1+A_matrix(2,3)*z1
          z2=A_matrix(3,1)*x1+A_matrix(3,2)*y1+A_matrix(3,3)*z1
          r2=x2**2+y2**2+z2**2

          zu_GS(i,ib,ik)=exp(-0.5d0*r2)
        enddo
      end if
    enddo
  enddo

  if(Myrank == 0)write(*,'(A)')'wave-functions initialization is completed'

  return
End Subroutine init_wf
