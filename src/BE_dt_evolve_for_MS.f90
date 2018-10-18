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
subroutine BE_dt_evolve_for_MS(Act_t)
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: Act_t
  real(8),parameter :: eps_Act = 1d-6
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

  call BE_dt_full_evolve_Taylor


  return
  contains

    subroutine BE_dt_full_evolve_Taylor
      implicit none
      integer,parameter :: nexp_taylor = 4
      complex(8) :: zvec0(NB_basis, NB_TD)
      complex(8) :: zhvec0(NB_basis, NB_TD)
      

      K_point : do ik=NK_s,NK_e

        zfact = 1d0
        zvec0 = zCt(1:NB_basis, 1:NB_TD, ik)
        do iexp = 1, nexp_taylor
          zfact = zfact*(-zI*dt)/iexp
!          zhvec0(1:NB_basis, 1:NB_TD) = matmul(&
!            zH_tot(1:NB_basis,1:NB_basis,ik),&
!            zvec0(1:NB_TD))
          call zhemm('L', 'U', NB_basis, NB_TD, 1d0, &
            zH_tot(1:NB_basis,1:NB_basis,ik), &
            NB_basis,&
            zvec0(1:NB_basis,1:NB_TD), &
            NB_basis, 0d0, &
            zhvec0(1:NB_basis,1:NB_TD), &
            NB_basis)
            

          zCt(1:NB_basis, 1:NB_TD, ik) = zCt(1:NB_basis, 1:NB_TD, ik) &
            +zfact*zhvec0(1:NB_basis, 1:NB_TD)
          if(iexp /= nexp_taylor)zvec0 = zhvec0

        end do
      end do K_point

    end subroutine BE_dt_full_evolve_Taylor


  end subroutine BE_dt_evolve_for_MS
