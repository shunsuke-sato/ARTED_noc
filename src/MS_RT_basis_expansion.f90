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
subroutine MS_RT_basis_expansion
  use global_variables
  use ms_maxwell_ks_variables
  implicit none
  integer :: iter,iter_t,ix
  real(8) :: jav,Act_t
  character(64) :: cit, cfilename

  if(myrank == 0)write(*,"(A)")"== Start real-time propagation with basis expansion."

  call init_wf_basis_expansion
  call init_Ac_ms_basis_expansion

  call MS_current(jt_m,Ac_m)
  if(myrank == 0)open(30,file='Ac_vac_rear.out')
  if(myrank == 0)write(30,"(999e26.16e3)")dt*0,Ac_m(0),Ac_m(Mx+1)

  do iter=0,Nt
    if(myrank == 0 .and. (mod(iter,200)==0 .or. iter == Nt))then
      write(cit,"(I9.9)")iter
      cfilename = 'Act_'//trim(cit)//'.out'
      open(41,file=cfilename)
      do ix = nx_s, nx_e
        write(41,"(999e26.16e3)")x_m(ix),Ac_m(ix)
      end do
      close(41)

    end if

! Compute Ac_m_n from Ac_m and Ac_m_o
    call dt_evolve_macro_field

    Act_t = 0.5d0*(Ac_m_n(macro_point_id) + Ac_m(macro_point_id))
    call BE_dt_evolve_for_MS(Act_t)

    call MS_current(jt_m,Ac_m_n)

    Ac_m_o = Ac_m
    Ac_m   = Ac_m_n
    if(myrank == 0)write(30,"(999e26.16e3)")dt*(iter+1),Ac_m(0),Ac_m(Mx+1)



  end do

  if(myrank == 0)close(30)
  if(myrank == 0)write(*,"(A)")"== End real-time propagation with basis expansion."

  return
end subroutine MS_RT_basis_expansion
