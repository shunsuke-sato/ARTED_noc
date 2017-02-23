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
subroutine BE_dt_evolve_Houston_probe_decomp(iter,Act_t)
  use global_variables
  use PSE_variables
  implicit none
  integer,intent(in) :: iter
  real(8),intent(in) :: Act_t
  real(8),parameter :: eps_Act = 1d-6
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact
  real(8) :: Act_tmp,Act_probe_tmp
  complex(8) :: zvec_t(NB_basis,NB_TD)
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  complex(8),allocatable :: zMat_diag(:,:)

! Normal propagation
  if(Ac_probe_BE(iter) == 0d0 .and. Ac_probe_BE(iter+1) == 0d0)then
    call BE_dt_evolve(Act_t)
    return
  end if


  lwork=6*NB_basis
  allocate(work_lp(lwork),rwork(3*NB_basis-2),w(NB_basis))
  allocate(zMat_diag(NB_basis,NB_basis))

! Propagation with Houston decomposition for probe Hamiltonian
! Here, we employ etrs propagation scheme

!!=== start the propagation (t => t + dt/2) ===

  Act_tmp = 0.5d0*(Ac_pump_BE(iter) + Ac_pump_BE(iter+1) )
  Act_probe_tmp = 0.5d0*(Ac_probe_BE(iter) + Ac_probe_BE(iter+1) )
!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if


  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot = zH_tot + zH_loc + zPi_loc*Act_tmp
  do ib = 1,NB_basis
    zH_tot(ib,ib,:) = zH_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do


!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zPi_tot(:,:,:) = 0.5d0*zPi_NL(:,:,:,iav_t+1)*(xx**2+xx) &
                  +0.5d0*zPi_NL(:,:,:,iav_t-1)*(xx**2-xx) &
                  +      zPi_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zPi_tot = zPi_tot + zPi_loc
  do ib = 1,NB_basis
    zPi_tot(ib,ib,:) = zPi_tot(ib,ib,:) + Act_tmp
  end do
!== construct current matrix end

  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_tot(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zdH_tot(:,:,ik) = matmul(zPi_tot(:,:,ik),zMat_diag(:,:))
    zPi_tot(:,:,ik) = matmul(transpose(conjg(zMat_diag(:,:))),zdH_tot(:,:,ik))
    zdH_tot(:,:,ik) = zPi_tot(:,:,ik)*Act_probe_tmp*Mask_probe(:,:)

    zvec_t(:,:) = matmul(transpose(conjg(zMat_diag(:,:))),zCt(:,:,ik))

    ztCt_Lan = zvec_t
    zfact = 1d0
    do iexp = 1,4
      zfact = zfact*(-zI*0.5d0*dt)/iexp
      zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
      zvec_t = zvec_t + zfact * zACt_Lan
      ztCt_Lan = zACt_Lan
    end do

    do ib = 1,NB_TD
      zvec_t(:,ib) = exp(-zI*w(:)*dt)*zvec_t(:,ib)
    end do

    ztCt_Lan = zvec_t
    zfact = 1d0
    do iexp = 1,4
      zfact = zfact*(-zI*0.5d0*dt)/iexp
      zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
      zvec_t = zvec_t + zfact * zACt_Lan
      ztCt_Lan = zACt_Lan
    end do

    zCt(:,:,ik) = matmul(zMat_diag(:,:),zvec_t(:,:))

  end do
!!=== end the propagation (t => t + dt/2) ===


  return
end subroutine BE_dt_evolve_Houston_probe_decomp
  

subroutine BE_dt_evolve_Houston_probe_decomp_org(iter,Act_t)
  use global_variables
  use PSE_variables
  implicit none
  integer,intent(in) :: iter
  real(8),intent(in) :: Act_t
  real(8),parameter :: eps_Act = 1d-6
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact
  real(8) :: Act_tmp
  complex(8) :: zvec_t(NB_basis,NB_TD)
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  complex(8),allocatable :: zMat_diag(:,:)

! Normal propagation
  if(Ac_probe_BE(iter) == 0d0 .and. Ac_probe_BE(iter+1) == 0d0)then
    call BE_dt_evolve(Act_t)
    return
  end if


  lwork=6*NB_basis
  allocate(work_lp(lwork),rwork(3*NB_basis-2),w(NB_basis))
  allocate(zMat_diag(NB_basis,NB_basis))

! Propagation with Houston decomposition for probe Hamiltonian
! Here, we employ etrs propagation scheme

!!=== start the propagation (t => t + dt/2) ===

  Act_tmp = Ac_pump_BE(iter)
!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if


  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot = zH_tot + zH_loc + zPi_loc*Act_tmp
  do ib = 1,NB_basis
    zH_tot(ib,ib,:) = zH_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do


!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zPi_tot(:,:,:) = 0.5d0*zPi_NL(:,:,:,iav_t+1)*(xx**2+xx) &
                  +0.5d0*zPi_NL(:,:,:,iav_t-1)*(xx**2-xx) &
                  +      zPi_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zPi_tot = zPi_tot + zPi_loc
  do ib = 1,NB_basis
    zPi_tot(ib,ib,:) = zPi_tot(ib,ib,:) + Act_tmp
  end do
!== construct current matrix end

  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_tot(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zdH_tot(:,:,ik) = matmul(zPi_tot(:,:,ik),zMat_diag(:,:))
    zPi_tot(:,:,ik) = matmul(transpose(conjg(zMat_diag(:,:))),zdH_tot(:,:,ik))
    zdH_tot(:,:,ik) = zPi_tot(:,:,ik)*Ac_probe_BE(iter)*Mask_probe(:,:)

    zvec_t(:,:) = matmul(transpose(conjg(zMat_diag(:,:))),zCt(:,:,ik))

    ztCt_Lan = zvec_t
    zfact = 1d0
    do iexp = 1,4
      zfact = zfact*(-zI*0.5d0*0.5d0*dt)/iexp
      zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
      zvec_t = zvec_t + zfact * zACt_Lan
      ztCt_Lan = zACt_Lan
    end do

    do ib = 1,NB_TD
      zvec_t(:,ib) = exp(-zI*w(:)*0.5d0*dt)*zvec_t(:,ib)
    end do

    ztCt_Lan = zvec_t
    zfact = 1d0
    do iexp = 1,4
      zfact = zfact*(-zI*0.5d0*0.5d0*dt)/iexp
      zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
      zvec_t = zvec_t + zfact * zACt_Lan
      ztCt_Lan = zACt_Lan
    end do

    zCt(:,:,ik) = matmul(zMat_diag(:,:),zvec_t(:,:))

  end do
!!=== end the propagation (t => t + dt/2) ===


!!=== start the propagation (t +dt/2 => t + dt) ===

  Act_tmp = Ac_pump_BE(iter+1)
!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if


  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot = zH_tot + zH_loc + zPi_loc*Act_tmp
  do ib = 1,NB_basis
    zH_tot(ib,ib,:) = zH_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do


!== construct current matrix start
  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1
  if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zPi_tot(:,:,:) = 0.5d0*zPi_NL(:,:,:,iav_t+1)*(xx**2+xx) &
                  +0.5d0*zPi_NL(:,:,:,iav_t-1)*(xx**2-xx) &
                  +      zPi_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zPi_tot = zPi_tot + zPi_loc
  do ib = 1,NB_basis
    zPi_tot(ib,ib,:) = zPi_tot(ib,ib,:) + Act_tmp
  end do
!== construct current matrix end

  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_tot(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zdH_tot(:,:,ik) = matmul(zPi_tot(:,:,ik),zMat_diag(:,:))
    zPi_tot(:,:,ik) = matmul(transpose(conjg(zMat_diag(:,:))),zdH_tot(:,:,ik))
    zdH_tot(:,:,ik) = zPi_tot(:,:,ik)*Ac_probe_BE(iter+1)*Mask_probe(:,:)

    zvec_t(:,:) = matmul(transpose(conjg(zMat_diag(:,:))),zCt(:,:,ik))

    ztCt_Lan = zvec_t
    zfact = 1d0
    do iexp = 1,4
      zfact = zfact*(-zI*0.5d0*0.5d0*dt)/iexp
      zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
      zvec_t = zvec_t + zfact * zACt_Lan
      ztCt_Lan = zACt_Lan
    end do

    do ib = 1,NB_TD
      zvec_t(:,ib) = exp(-zI*w(:)*0.5d0*dt)*zvec_t(:,ib)
    end do

    ztCt_Lan = zvec_t
    zfact = 1d0
    do iexp = 1,4
      zfact = zfact*(-zI*0.5d0*0.5d0*dt)/iexp
      zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
      zvec_t = zvec_t + zfact * zACt_Lan
      ztCt_Lan = zACt_Lan
    end do

    zCt(:,:,ik) = matmul(zMat_diag(:,:),zvec_t(:,:))

  end do
!!=== end the propagation (t + dt/2 => t + dt) ===



  return
end subroutine BE_dt_evolve_Houston_probe_decomp_org
  
