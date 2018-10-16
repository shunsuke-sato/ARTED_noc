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
subroutine init_Ac_ms_basis_expansion
  use global_variables
  use ms_maxwell_ks_variables
  implicit none
  integer :: iter, ix
  real(8) :: tt
  real(8) :: f0_1,f0_2,omega_1,omega_2,tpulse_1,tpulse_2,T1_T2

  if(myrank == 0)write(*,"(A)")"== Start: Initialization of vector potential."
  allocate(Ac_m(nx_s:nx_e), jt_m(Mx))
  allocate(Ac_m_n(nx_s:nx_e), Ac_m_o(nx_s:nx_e))
  allocate(x_m(nx_s:nx_e))

  do ix = nx_s, nx_e
    x_m(ix) = ix*dx_m -0.5d0*dx_m
  end do


  Ac_m = 0d0

  f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
  omega_1=omegaev_1/(2d0*Ry)  ! frequency in a.u.
  tpulse_1=tpulsefs_1/0.024189d0 ! pulse duration in a.u.
  f0_2=5.338d-9*sqrt(IWcm2_2)      ! electric field in a.u.
  omega_2=omegaev_2/(2d0*Ry)  ! frequency in a.u.
  tpulse_2=tpulsefs_2/0.024189d0 ! pulse duration in a.u.
  T1_T2=T1_T2fs/0.024189d0 ! pulse duration in a.u.

  
  select case(laser_type)
  case('cos2cos')
! pump
    do ix = nx_s, nx_e

!t=0
      tt = 0d0 -x_m(ix)/clight-0.5d0*tpulse_1
      if(abs(tt)<0.5d0*tpulse_1)then
        Ac_m(ix) = -f0_1/omega_1*cos(pi*tt/tpulse_1)**2*sin(omega_1*tt)
      end if

!t=-dt
      tt = -dt -x_m(ix)/clight-0.5d0*tpulse_1
      if(abs(tt)<0.5d0*tpulse_1)then
        Ac_m_o(ix) = -f0_1/omega_1*cos(pi*tt/tpulse_1)**2*sin(omega_1*tt)
      end if

    end do
! pump
    do ix = nx_s, nx_e

!t=0
      tt = 0d0 -x_m(ix)/clight-0.5d0*tpulse_1 - T1_T2
      if(abs(tt)<0.5d0*tpulse_2)then
        Ac_m(ix) = Ac_m(ix) -f0_2/omega_2*cos(pi*tt/tpulse_2)**2*sin(omega_2*tt)
      end if

      tt = -dt -x_m(ix)/clight-0.5d0*tpulse_1 - T1_T2
      if(abs(tt)<0.5d0*tpulse_2)then
        Ac_m_o(ix) = Ac_m_o(ix) -f0_2/omega_2*cos(pi*tt/tpulse_2)**2*sin(omega_2*tt)
      end if


    end do


  case('cos4cos')
! pump
    do ix = nx_s, nx_e

!t=0
      tt = 0d0 -x_m(ix)/clight-0.5d0*tpulse_1
      if(abs(tt)<0.5d0*tpulse_1)then
        Ac_m(ix) = -f0_1/omega_1*cos(pi*tt/tpulse_1)**4*sin(omega_1*tt)
      end if

!t=-dt
      tt = -dt -x_m(ix)/clight-0.5d0*tpulse_1
      if(abs(tt)<0.5d0*tpulse_1)then
        Ac_m_o(ix) = -f0_1/omega_1*cos(pi*tt/tpulse_1)**4*sin(omega_1*tt)
      end if

    end do
! pump
    do ix = nx_s, nx_e

!t=0
      tt = 0d0 -x_m(ix)/clight-0.5d0*tpulse_1 - T1_T2
      if(abs(tt)<0.5d0*tpulse_2)then
        Ac_m(ix) = Ac_m(ix) -f0_2/omega_2*cos(pi*tt/tpulse_2)**4*sin(omega_2*tt)
      end if

      tt = -dt -x_m(ix)/clight-0.5d0*tpulse_1 - T1_T2
      if(abs(tt)<0.5d0*tpulse_2)then
        Ac_m_o(ix) = Ac_m_o(ix) -f0_2/omega_2*cos(pi*tt/tpulse_2)**4*sin(omega_2*tt)
      end if


    end do
  case default
    err_message='error in init_Ac'
    call err_finalize
  end select



  if(myrank == 0)write(*,"(A)")"== End: Initialization of vector potential."

  return
end subroutine init_Ac_ms_basis_expansion
