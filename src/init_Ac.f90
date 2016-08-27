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
subroutine init_Ac
  use global_variables
  implicit none
  integer :: iter
  real(8) :: tt
  real(8) :: f0_1,f0_2,omega_1,omega_2,tpulse_1,tpulse_2,T1_T2

  allocate(Acext_Cvec(3,0:Nt+2),Acind_Cvec(3,0:Nt+2) &
    ,Actot_Cvec(3,0:Nt+2),javt_Cvec(3,0:Nt+2))

  f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
  omega_1=omegaev_1/(2d0*Ry)  ! frequency in a.u.
  tpulse_1=tpulsefs_1/0.02418d0 ! pulse duration in a.u.
  f0_2=5.338d-9*sqrt(IWcm2_2)      ! electric field in a.u.
  omega_2=omegaev_2/(2d0*Ry)  ! frequency in a.u.
  tpulse_2=tpulsefs_2/0.02418d0 ! pulse duration in a.u.
  T1_T2=T1_T2fs/0.02418d0 ! pulse duration in a.u.
  javt_Cvec=0.d0
  Acext_Cvec=0.d0; Acind_Cvec=0.d0; Actot_Cvec=0.d0

! 'cos2cos'                          ! laser_type 'cos2cos', 'cos4cos' or 'impulse'

  select case(laser_type)
  case('impulse')
    Acext_Cvec(1,:)=Epdir_1(1)*dAc
    Acext_Cvec(2,:)=Epdir_1(2)*dAc
    Acext_Cvec(3,:)=Epdir_1(3)*dAc
  case('cos2cos')
! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
    do iter=0,Nt+2
      tt=iter*dt
      if (tt<tpulse_1) then
        Acext_Cvec(:,iter)=-Epdir_1(:)*f0_1/omega_1*(sin(pi*tt/tpulse_1))**2*cos(omega_1*tt+phi_CEP_1*2d0*pi)
      end if
    enddo
! probe laser
    do iter=0,Nt+2
      tt=iter*dt
      if ( (tt-T1_T2>0d0) .and. (tt-T1_T2<tpulse_2) ) then
        Acext_Cvec(:,iter)=Acext_Cvec(:,iter) &
          &-Epdir_2(:)*f0_2/omega_2*(sin(pi*(tt-T1_T2)/tpulse_2))**2*cos(omega_2*(tt-T1_T2)+phi_CEP_2*2d0*pi)
      endif
    enddo
  case('cos4cos')
! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
    do iter=0,Nt+2
      tt=iter*dt
      if (tt<tpulse_1) then
        Acext_Cvec(:,iter)=-Epdir_1(:)*f0_1/omega_1*(cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1))**4*sin(omega_1*(tt-0.5d0*tpulse_1)+phi_CEP_1*2d0*pi)
      end if
    enddo
! probe laser
    do iter=0,Nt+2
      tt=iter*dt
      if ( (tt>0.5d0*(tpulse_1-tpulse_2)+T1_T2) .and. (tt>0.5d0*(tpulse_1-tpulse_2)+T1_T2+tpulse_2) ) then
        Acext_Cvec(:,iter)=Acext_Cvec(:,iter) &
          &-Epdir_2(:)*f0_2/omega_2*(cos(pi*(tt-(0.5d0*tpulse_1+T1_T2))/tpulse_2))**4*sin(omega_2*(tt-(0.5d0*tpulse_1+T1_T2))+phi_CEP_2*2d0*pi)
      endif
    enddo
  case default
    err_message='error in init_Ac'
    call err_finalize
  end select




  return
end subroutine init_Ac
