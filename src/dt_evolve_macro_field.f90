!
!  Copyright 2018 S.A. Sato
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
subroutine dt_evolve_macro_field
  use global_variables
  use ms_maxwell_ks_variables
  implicit none
  real(8) :: lap0, lap1
  real(8) :: Lap_Ac(nx_s:nx_e)
  integer :: ix

  lap0 = -2d0/dx_m**2
  lap1 = 1d0/dx_m**2

  Lap_Ac(nx_s) = lap0*Ac_m(nx_s) + lap1Ac_m(nx_s+1)
  do ix = nx_s+1,nx_e-1
    Lap_Ac(ix) = lap0*Ac_m(ix) + lap1*(Ac_m(ix+1)+Ac_m(ix-1))
  end do
  Lap_Ac(nx_e) = lap0*Ac_m(nx_e) + lap1*Ac_m(nx_e-1)

  Ac_m_n = 2d0*Ac_m - Ac_m_o + (dt*clight)**2*Lap_Ac
  Ac_m_n(1:Mx) = Ac_m_n(1:Mx) -4d0*pi*dt**2*jt_m(1:Mx)


end subroutine dt_evolve_macro_field

