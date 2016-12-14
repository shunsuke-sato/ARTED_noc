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
subroutine BE_dt_evolve(Act_t)
  use global_variables
  use PSE_variables
  real(8) :: Act_t
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact

!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_t-dble(iav)*dAmax) < diff)then
      diff = abs(Act_t-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(abs(iav_t) == NAmax)then
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

  K_point : do ik=NK_s,NK_e
  Band : do ib=1,NB_TD

    zfact = 1d0
    ztCt_tmp(:)=zCt(:,ib,ik)
    do iexp = 1,4
      zfact = zfact*(-zI*dt)/dble(iexp)
      zACt_tmp(:) = matmul(zH_tot(:,:,ik),ztCt_tmp(:))
      zCt(:,ib,ik) = zCt(:,ib,ik) + zfact*zACt_tmp(:)
      ztCt_tmp(:) = zACt_tmp(:)
    end do

  end do Band
  end do K_point

  return
end subroutine BE_dt_evolve
  
