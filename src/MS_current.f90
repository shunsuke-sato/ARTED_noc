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
subroutine MS_current(jt_m_out,Ac_m_in)
  use global_variables
  use PSE_variables
  use ms_maxwell_ks_variables
  implicit none
  real(8),intent(out) :: jt_m_out(Mx)
  real(8),intent(in) :: Ac_m_in(Nx_s:Nx_e)
  real(8) :: jt_m_l(Mx)
  real(8) :: jav,Act_t
  real(8) :: jav_l
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib
  complex(8) :: zvec0(NB_basis, NB_TD)
  integer :: ix_m


  jt_m_l = 0d0
  M_point : do ix_m=Mx_s, Mx_e
  

  Act_t = Ac_m_in(ix_m)
  if(Act_t /= Act_t)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_t) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if
  
!== construct current matrix start
  diff = 1d10
  iav_t = -NAmax
  do iav = -NAmax,NAmax
    if( abs(Act_t-dble(iav)*dAmax) < diff)then
      diff = abs(Act_t-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1

  xx = (Act_t-dble(iav_t)*dAmax)/dAmax
  zPi_tot(:,:,:) = 0.5d0*zPi_NL(:,:,:,iav_t+1)*(xx**2+xx) &
                  +0.5d0*zPi_NL(:,:,:,iav_t-1)*(xx**2-xx) &
                  +      zPi_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zPi_tot = zPi_tot + zPi_loc
  do ib = 1,NB_basis
    zPi_tot(ib,ib,:) = zPi_tot(ib,ib,:) + Act_t
  end do
!== construct current matrix end

  jav_l = 0d0
!$omp parallel do private(ik, ib, zvec0) reduction(+:jav_l)
  K_point : do ik=NK_s,NK_e

    zvec0(1:NB_basis, 1:NB_TD) = matmul(zPi_tot(:,:,ik),zCt_Mpoint(:,1:NB_TD,ik, ix_m))
    Band : do ib=1,NB_TD
      jav_l = jav_l + occ(ib,ik)*sum(conjg(zCt_Mpoint(:,ib,ik,ix_m))*zvec0(:, ib))
    end do Band
  end do K_point

  jav_l=jav_l/Vcell

  jt_m_l(ix_m) = jt_m_l(ix_m) + jav_l
  end do M_point
  

  call MPI_ALLREDUCE(jt_m_l,jt_m_out,Mx,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)  

  return
end subroutine MS_current
