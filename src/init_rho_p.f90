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
subroutine init_rho_p
  use global_variables
  implicit none
  real(8),parameter :: sigma_g=2.0d0
  integer :: ia1,ia2,ia3,i,j
  real(8) :: an1,an2,an3
  integer :: iNI
  real(8) :: r2
  real(8) :: dist_Lvec(3)

  rho_p=0d0
  select case(method)
  case('PAW')
    err_message='In preparation PAW method' ; call err_finalize
    call err_finalize
  case('PSE')
    return
  case('CHK')
    do iNI=1,NI
      do ia1=-8,8
        do ia2=-8,8
          do ia3=-8,8
            do i=1,NL
              dist_Lvec(1)=Lx(1,i)-(Rion_Lvec(1,iNI)+dble(ia1))
              dist_Lvec(2)=Lx(2,i)-(Rion_Lvec(2,iNI)+dble(ia2))
              dist_Lvec(3)=Lx(3,i)-(Rion_Lvec(3,iNI)+dble(ia3))
              
              r2=0d0
              do j=1,3
                r2=r2+dist_Lvec(j)*sum(dist_Lvec(:)*mat_vv_a_Cvec(j,:))
              end do
              
              rho_p(i)=rho_p(i)+dble(Zatom(Kion(iNi)))*(1d0/sqrt(2d0*pi*sigma_g**2))**3*exp(-0.5d0*r2/(sigma_g**2))
            end do
          end do
        end do
      end do
    end do
  case default
    err_message='Invalid option "mehtod"' ; call err_finalize
  end select

  return
end subroutine init_rho_p
