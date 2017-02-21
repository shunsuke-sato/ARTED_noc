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
subroutine PSE_real_time_propagation_basis_expansion
  use global_variables
  implicit none
  integer :: iter,iter_t
  real(8) :: jav,Act_t

  if(myrank == 0)write(*,"(A)")"== Start real-time propagation with basis expansion."

  call init_wf_basis_expansion
  if(switch_Houston_probe_decomposition)then
    call init_Ac_basis_expansion_Houston_probe_decomp
  else
    call init_Ac_basis_expansion
  end if

  if(myrank == 0)open(104,file=trim(SYSname)//'_eex.out')
  Act_t = 0d0
  call BE_energy(Act_t); Etot_GS = Etot_RT
  if(myrank==0)write(104,"(999e26.16e3)")0d0,Etot_RT,Etot_RT-Etot_GS

!== Start current
  Act_t = Actot_BE(0)
  call BE_current(jav,Act_t)
  javt_BE(0)=jav
!== End current

  do iter=0,Nt
!== Start dt_evolve
    if(switch_Houston_probe_decomposition)then
      Act_t = 0.5d0*( Actot_BE(iter+1) + Actot_BE(iter) )
      call BE_dt_evolve_Houston_probe_decomp(iter,Act)
    else
      Act_t = 0.5d0*( Actot_BE(iter+1) + Actot_BE(iter) )
      call BE_dt_evolve(Act_t)
    end if
!== End dt_evolve

!== Start current
  Act_t = Actot_BE(iter+1)
  call BE_current(jav,Act_t)
  javt_BE(iter+1)=jav
!== End current

!== Start energy
  if(mod(iter,10) == 0 .or. iter==Nt)then
    call BE_energy(Act_t)
    if(myrank==0)write(104,"(999e26.16e3)")dt*dble(iter+1),Etot_RT,Etot_RT-Etot_GS
  end if

!== Start writing section
  if(mod(iter,400) == 0 .or. iter == Nt)then
    if(myrank == 0)then
      open(102,file=trim(SYSname)//'_jac.out')
      do iter_t=0,Nt+1
        write(102,'(100e26.16e3)')Dt*dble(iter_t),Actot_BE(iter_t),javt_BE(iter_t)
      end do
      close(102)
    end if
  end if
!== End writing section

  end do

  call BE_excited_electrons
  call BE_matrix_element_for_dielectric_function
  if(myrank == 0)close(104)
  if(myrank == 0)write(*,"(A)")"== End real-time propagation with basis expansion."

  return
end subroutine PSE_real_time_propagation_basis_expansion
