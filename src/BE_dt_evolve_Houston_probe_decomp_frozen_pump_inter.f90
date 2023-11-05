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
subroutine BE_dt_evolve_Houston_probe_decomp_frozen_pump_inter(iter)
  use global_variables
  use PSE_variables
  implicit none
  integer,intent(in) :: iter
  real(8),parameter :: eps_Act = 1d-6
  integer,parameter :: NTaylor = 4
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp,ib1,ib2
  complex(8) :: zfact
  real(8) :: Act_tmp,Act_probe_tmp
  complex(8) :: zvec_t(NB_basis,NB_TD),zvec_t2(NB_basis,NB_TD)
  integer :: ierr_lan
  complex(8),allocatable :: zMat_tmp(:,:)
  complex(8),allocatable :: zMat_tmp2(:,:)
  complex(8),allocatable :: zUm(:,:)
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:),w2(:)
  integer :: info
  complex(8),allocatable :: zMat_diag(:,:)
  complex(8),allocatable :: zMat_diag2(:,:)


  lwork=6*NB_basis
  allocate(work_lp(lwork),rwork(3*NB_basis-2),w(NB_basis),w2(NB_basis))
  allocate(zMat_diag(NB_basis,NB_basis))
  allocate(zMat_diag2(NB_basis,NB_basis))
  allocate(zMat_tmp(NB_basis,NB_basis))
  allocate(zMat_tmp2(NB_basis,NB_basis))
  allocate(zUm(NB_basis,NB_basis))

! Propagation with Houston decomposition for probe Hamiltonian
! Here, we employ etrs propagation scheme

!==construct Pump hamiltonian at t==

  Act_tmp = Ac_pump_BE(iter)
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot = zH_tot + zH_loc + Act_tmp*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zH_tot(ib,ib,:) = zH_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do

!==construct probe hamiltonian at t==

  Act_tmp = Ac_pump_BE(iter) + Ac_probe_BE(iter)
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zdH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zdH_tot = zdH_tot + zH_loc + Act_tmp*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zdH_tot(ib,ib,:) = zdH_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do

  zdH_tot = zdH_tot - zH_tot


!==construct Pump hamiltonian at t+dt ==

  Act_tmp = Ac_pump_BE(iter+1)
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot2(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot2 = zH_tot2 + zH_loc + Act_tmp*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zH_tot2(ib,ib,:) = zH_tot2(ib,ib,:) + 0.5d0*Act_tmp**2
  end do


!==construct probe hamiltonian at t + dt ==

  Act_tmp = Ac_pump_BE(iter+1) + Ac_probe_BE(iter+1)
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

  diff = 1d10
  do iav = -NAmax,NAmax
    if( abs(Act_tmp-dble(iav)*dAmax) < diff)then
      diff = abs(Act_tmp-dble(iav)*dAmax)
      iav_t = iav
    end if
  end do
  if(iav_t == NAmax)iav_t=NAmax -1
  if(iav_t == -NAmax)iav_t=-NAmax +1

  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zdH_tot2(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zdH_tot2 = zdH_tot2 + zH_loc + Act_tmp*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zdH_tot2(ib,ib,:) = zdH_tot2(ib,ib,:) + 0.5d0*Act_tmp**2
  end do

  zdH_tot2 = zdH_tot2 - zH_tot2



  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_tot(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zMat_diag2(:,:)=zH_tot2(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag2, NB_basis, w2, work_lp, lwork, rwork, info)

! probe masking at t
    zMat_tmp(:,:) = matmul(zdH_tot(:,:,ik),zMat_diag(:,:))
    zdH_tot(:,:,ik) = matmul(transpose(conjg(zMat_diag(:,:))),zMat_tmp(:,:))
    zdH_tot(:,:,ik) = zdH_tot(:,:,ik)*Mask_probe(:,:)

    zMat_tmp(:,:) = matmul(zdH_tot(:,:,ik),transpose(conjg(zMat_diag(:,:))))
    zdH_tot(:,:,ik) = matmul(zMat_diag(:,:),zMat_tmp(:,:))

! probe masking at t + dt
    zMat_tmp2(:,:) = matmul(zdH_tot2(:,:,ik),zMat_diag2(:,:))
    zdH_tot2(:,:,ik) = matmul(transpose(conjg(zMat_diag2(:,:))),zMat_tmp2(:,:))
    zdH_tot2(:,:,ik) = zdH_tot2(:,:,ik)*Mask_probe(:,:)

    zMat_tmp2(:,:) = matmul(zdH_tot2(:,:,ik),transpose(conjg(zMat_diag2(:,:))))
    zdH_tot2(:,:,ik) = matmul(zMat_diag2(:,:),zMat_tmp2(:,:))

! propagation with probe hamiltonian from t to t+dt/2
    call BE_dt_probe_ham_Taylor(zdH_tot(:,:,ik), zCt(:,:,ik), 0.5d0*dt)

!-- Start: Pump propagation

! construct the propagation matrix
    do ib1=1,nb_basis
      do ib2=1,nb_basis

        zUm(ib1,ib2) = sum(conjg(zMat_diag2(:,ib1))*zMat_diag(:,ib2)) &
            *exp(-zi*0.5d0*dt*(w2(ib1)+w(ib2)))

      end do
    end do

    zUm=zUm*Mask_pump(:,:)
    call gram_schmidt(nb_basis, zUm)

    zvec_t = matmul(transpose(conjg(zMat_diag(:,:))),zCt(:,:,ik))
    zvec_t2 = matmul(zUm,zvec_t)
    zCt(:,:,ik) = matmul(zMat_diag2(:,:),zvec_t2)


!-- End:   Pump propagation





! propagation with probe hamiltonian from t+dt/2 to t+dt
    call BE_dt_probe_ham_Taylor(zdH_tot2(:,:,ik), zCt(:,:,ik), 0.5d0*dt)

  end do


  return
contains

  subroutine BE_dt_probe_ham_Taylor(zH_ham, zpsi, dt_t)
    implicit none
    complex(8),intent(in) :: zH_ham(NB_basis, NB_basis)
    complex(8),intent(inout) :: zpsi(NB_basis, NB_TD)
    real(8),intent(in) :: dt_t
    complex(8) :: zvec(NB_basis, NB_TD)
    complex(8) :: zhvec(NB_basis, NB_TD)
    integer,parameter :: ntaylor = 4
    integer :: iexp
    complex(8) :: zfact

    zfact = 1d0


    zvec = zpsi
    do iexp = 1, ntaylor
      zfact = zfact*(-zi*dt_t)/iexp
      call zhemm('L', 'U', NB_basis, NB_TD, (1d0,0d0), &
          zH_ham(1:NB_basis,1:NB_basis), &
          NB_basis,&
          zvec(1:NB_basis,1:NB_TD), &
          NB_basis, (0d0,0d0), &
          zhvec(1:NB_basis,1:NB_TD), &
          NB_basis)

      zpsi = zpsi + zfact*zhvec
      zvec = zhvec

    end do
  end subroutine BE_dt_probe_ham_Taylor

  subroutine gram_schmidt(n, zmat)
    implicit none
    integer,intent(in) :: n
    complex(8),intent(inout) :: zmat(n,n)
    real(8) :: ss
    complex(8) :: zs
    integer :: i,j


    do i = 1, n

      ss = sqrt( sum(abs(zmat(:,i))**2) )
      zmat(:,i) = zmat(:,i)/ss

      do j = 1,i-1

        zs = sum(conjg(zmat(:,j))*zmat(:,i))
        zmat(:,i) = zmat(:,i)-zs*zmat(:,j)

      end do

      ss = sqrt( sum(abs(zmat(:,i))**2) )
      zmat(:,i) = zmat(:,i)/ss

    end do

  end subroutine gram_schmidt
end subroutine BE_dt_evolve_Houston_probe_decomp_frozen_pump_inter

