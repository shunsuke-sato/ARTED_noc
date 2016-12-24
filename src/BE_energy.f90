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
subroutine BE_energy(Act_t)
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: Act_t
  real(8) :: Etot_RT_l
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib

!  kAc_Cvec(1,:)=kAc0_Cvec(1,:)+Act_t*Epdir_1(1)
!  kAc_Cvec(2,:)=kAc0_Cvec(2,:)+Act_t*Epdir_1(2)
!  kAc_Cvec(3,:)=kAc0_Cvec(3,:)+Act_t*Epdir_1(3)

!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_t-dble(iav)*dAmax) < diff)then
      diff = abs(Act_t-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_t) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  xx = (Act_t-dble(iav_t)*dAmax)/dAmax
  zH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot = zH_tot + zH_loc + zPi_loc*Act_t
  do ib = 1,NB_basis
    zH_tot(ib,ib,:) = zH_tot(ib,ib,:) + 0.5d0*Act_t**2
  end do
!== construct current matrix end

  Etot_RT_l = 0d0
  K_point : do ik=NK_s,NK_e
  Band : do ib=1,NB_TD

    zACt_tmp(:) = matmul(zH_tot(:,:,ik),zCt(:,ib,ik))
    Etot_RT_l = Etot_RT_l + occ(ib,ik)*sum(conjg(zCt(:,ib,ik))*zACt_tmp(:))

  end do Band
  end do K_point

  call MPI_ALLREDUCE(Etot_RT_l,Etot_RT,1,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)  

  return
end subroutine BE_energy
