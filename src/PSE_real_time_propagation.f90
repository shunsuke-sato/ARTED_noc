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
subroutine PSE_real_time_propagation
  use global_variables
  use PSE_variables
  implicit none
  integer :: iter,iter_t
  real(8) :: jav(3)
  real(8) :: Etot,Ekin

  zu(:,1:NB_TD,:)=zu_GS(:,1:NB_TD,:)
  zu_GS0(:,:,:)=zu_GS(:,:,:)

!== START open files
  if(myrank == 0)then
    open(103,file=trim(SYSname)//'_nex.out')
    open(104,file=trim(SYSname)//'_Eex.out')
  end if
!== END open files

  call init_Ac
  Actot_Cvec(:,0) = Acext_Cvec(:,0) + Acind_Cvec(:,0)

!== Start current
  kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Actot_Cvec(1,0)
  kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Actot_Cvec(2,0)
  kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Actot_Cvec(3,0)
  call PSE_current_RT(jav)
  javt_Cvec(:,0)=jav(:)
!== End current

  if(Tr_Lo == 'Tr')then
    Acind_Cvec(:,1) = 0d0
  else if(Tr_Lo == 'Lo')then
    Acind_Cvec(:,1) = 2d0*Acind_Cvec(:,0) -Acind_Cvec(:,0) - 4d0*pi*javt_Cvec(:,0)*dt**2
  end if
  Actot_Cvec(:,1) = Acext_Cvec(:,1) + Acind_Cvec(:,1)

  do iter=0,Nt

    call PSE_dt_evolve(iter)

!== Start current
    kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Actot_Cvec(1,iter+1)
    kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Actot_Cvec(2,iter+1)
    kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Actot_Cvec(3,iter+1)

    call PSE_current_RT(jav)
    javt_Cvec(:,iter+1)=jav(:)
!== End current

    if(Tr_Lo == 'Tr')then
      Acind_Cvec(:,iter+2) = 0d0
    else if(Tr_Lo == 'Lo')then
      Acind_Cvec(:,iter+2) = 2d0*Acind_Cvec(:,iter+1) - Acind_Cvec(:,iter)&
        - 4d0*pi*javt_Cvec(:,iter+1)*dt**2
    end if
    Actot_Cvec(:,iter+2) = Acext_Cvec(:,iter+2) + Acind_Cvec(:,iter+2)

!== Start write section
!== current
    if(mod(iter,400) == 0 .or. iter == Nt)then
      if(myrank == 0)then
        open(102,file=trim(SYSname)//'_jac.out')
        do iter_t=0,Nt+1
          write(102,'(100e26.16e3)')Dt*dble(iter_t),javt_Cvec(1,iter_t),javt_Cvec(2,iter_t),javt_Cvec(3,iter_t) &
            ,Acext_Cvec(1,iter_t),Acext_Cvec(2,iter_t),Acext_Cvec(3,iter_t)&
            ,Acind_Cvec(1,iter_t),Acind_Cvec(2,iter_t),Acind_Cvec(3,iter_t)&
            ,Actot_Cvec(1,iter_t),Actot_Cvec(2,iter_t),Actot_Cvec(3,iter_t)
        end do
        close(102)
      end if
    end if

!== band
    if( ( mod(iter,Nstep_TD_band_dos) == 0 .and. option_TD_band_dos == 'Y' ) .or. iter == Nt) call TD_Band_DOS(iter)

!== Eex
    if( mod(iter,50) == 0 )then
      call PSE_energy(Etot,Ekin,'RT')
      if( myrank == 0 )write(104,'(100e26.16e3)')Dt*dble(iter),Etot,Ekin,Etot-Etot_GS,Ekin-Ekin_GS
    end if


!== Start write section

  end do

!== START close files
  if(myrank == 0)then
    close(103)
    close(104)
  end if
!== END close files

  return
end subroutine PSE_real_time_propagation
