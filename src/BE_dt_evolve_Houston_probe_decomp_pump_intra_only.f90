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
subroutine BE_dt_evolve_Houston_probe_decomp_pump_intra_only(iter,Act_t)
  use global_variables
  use PSE_variables
  implicit none
  integer,intent(in) :: iter
  real(8),intent(in) :: Act_t
  real(8),parameter :: eps_Act = 1d-6
  integer,parameter :: NTaylor = 4
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact
  real(8) :: Act_tmp,Act_probe_tmp
  complex(8) :: zvec_t(NB_basis,NB_TD)
  integer :: ierr_lan
  complex(8),allocatable :: zMat_tmp(:,:)
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  complex(8),allocatable :: zMat_diag(:,:)

  lwork=6*NB_basis
  allocate(work_lp(lwork),rwork(3*NB_basis-2),w(NB_basis))
  allocate(zMat_diag(NB_basis,NB_basis))
  allocate(zMat_tmp(NB_basis,NB_basis))

! Propagation with Houston decomposition for probe Hamiltonian
! Here, we employ etrs propagation scheme

!==construct Pump hamiltonian ==

  Act_tmp = 0.5d0*(Ac_pump_BE(iter) + Ac_pump_BE(iter+1) )
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

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



!==construct full hamiltonian ==

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
  zdH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zdH_tot = zdH_tot + zH_loc + Act_t*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zdH_tot(ib,ib,:) = zdH_tot(ib,ib,:) + 0.5d0*Act_t**2
  end do

! construct probe Hamiltonian
  zdH_tot = zdH_tot - zH_tot


!== construct current matrix end

  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_tot(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zMat_tmp(:,:) = matmul(zdH_tot(:,:,ik),zMat_diag(:,:))
    zdH_tot(:,:,ik) = matmul(transpose(conjg(zMat_diag(:,:))),zMat_tmp(:,:))
    zdH_tot(:,:,ik) = zdH_tot(:,:,ik)*Mask_probe(:,:)

    zMat_tmp(:,:) = matmul(zdH_tot(:,:,ik),transpose(conjg(zMat_diag(:,:))))
    zdH_tot(:,:,ik) = matmul(zMat_diag(:,:),zMat_tmp(:,:))

!    zH_tot(:,:,ik) = zH_tot(:,:,ik) + zdH_tot(:,:,ik)
    zH_tot(:,:,ik) = zdH_tot(:,:,ik) ! probe only
  end do



  call BE_half_dt_evolution_only_with_intraband(Ac_pump_BE(iter) &
      ,0.5d0*(Ac_pump_BE(iter)+Ac_pump_BE(iter+1)))
 
  call BE_dt_full_evolve_Krylov_exact_diag

  call BE_half_dt_evolution_only_with_intraband( &
      0.5d0*(Ac_pump_BE(iter)+Ac_pump_BE(iter+1)) &
      ,Ac_pump_BE(iter+1))

  return
contains
  subroutine BE_dt_full_evolve_Krylov_exact_diag
    implicit none
    integer,parameter :: nvec = 16
    complex(8) :: zvec(NB_basis, NB_TD, nvec)
    complex(8) :: zhvec(NB_basis, NB_TD, nvec)
    complex(8) :: zham_m(nvec, nvec),zUprop_m(nvec, nvec)
    integer :: ib, ivec, jvec, ik
    real(8) :: ss
    complex(8) :: zs
!LAPACK
    integer :: lwork
    complex(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info

    lwork=6*nvec+128
    allocate(work_lp(lwork),rwork(3*nvec-2),w(nvec))
      
    K_point : do ik=NK_s,NK_e

      zvec(1:NB_basis,1:NB_TD,1) = zCt(1:NB_basis, 1:NB_TD, ik)
!normalize
      do ib = 1, nb_td
        ss = sum(abs(zvec(:,ib,1))**2); ss = 1d0/sqrt(ss)
        zvec(:,ib,1)=zvec(:,ib,1)*ss
      end do

!Construction of Krylov subspace
      do ivec = 1, nvec

        call zhemm('L', 'U', NB_basis, NB_TD, (1d0,0d0), &
            zH_tot(1:NB_basis,1:NB_basis,ik), &
            NB_basis,&
            zvec(1:NB_basis,1:NB_TD, ivec), &
            NB_basis, (0d0,0d0), &
            zhvec(1:NB_basis,1:NB_TD,ivec), &
            NB_basis)
        
        if(ivec /= nvec)then
          do ib = 1, nb_td
            ss = sum(conjg(zvec(:,ib,ivec))*zhvec(:,ib,ivec))
            zvec(:,ib,ivec+1) = zhvec(:,ib,ivec)-ss*zvec(:,ib,ivec)
            ss = sum(abs(zvec(:,ib,ivec+1))**2)
            if(ss == 0d0)then
              write(*,"(A)")'Warning: (a) linear dependency in BE_dt_full_evolve_Krylov_exact_diag'
              zvec(:,ib,ivec+1)=1d0/sqrt(dble(NB_basis))
            else
              ss = 1d0/sqrt(ss)
              zvec(:,ib,ivec+1)=zvec(:,ib,ivec+1)*ss
            end if
          end do

!Gram-Schmidt orthonormalization
          do ib = 1, nb_td
            do jvec = 1, ivec
              zs = sum(conjg(zvec(:,ib,jvec))*zvec(:,ib,ivec+1))
              zvec(:,ib,ivec+1) = zvec(:,ib,ivec+1) -zs*zvec(:,ib,jvec)
            end do
            ss = sum(abs(zvec(:,ib,ivec+1))**2)
            if(ss == 0d0)then
              write(*,"(A)")'Warning: (b) linear dependency in BE_dt_full_evolve_Krylov_exact_diag'
              stop
            end if
            ss = 1d0/sqrt(ss)
            zvec(:,ib,ivec+1)=zvec(:,ib,ivec+1)*ss
          end do
          
        end if

      end do

      
      do ib = 1, nb_td
        do ivec = 1, nvec
          zham_m(ivec,ivec) = sum( conjg(zhvec(:,ib,ivec))*zvec(:,ib,ivec))
          do jvec = ivec+1, nvec
            zham_m(ivec,jvec) = sum( conjg(zhvec(:,ib,ivec))*zvec(:,ib,jvec))
            zham_m(jvec,ivec) = conjg(zham_m(ivec,jvec))
          end do
        end do
        
!diag
        call zheev('V', 'U', nvec, zham_m, nvec, w, work_lp, lwork, rwork, info)

        zUprop_m = 0d0
        do ivec = 1, nvec
          zUprop_m(ivec, ivec) = exp(-zi*dt*w(ivec))
        end do
        zUprop_m = matmul(matmul(zham_m,zUprop_m),conjg(transpose(zham_m)))
        zCt(1:NB_basis, ib, ik) = zvec(:,ib,1)*zUprop_m(1,1)
        do ivec = 2, nvec
          zCt(1:NB_basis, ib, ik) = &
              zCt(1:NB_basis, ib, ik) + zvec(:,ib,ivec)*zUprop_m(ivec,1)
        end do
        
      end do


    end do K_point

  end subroutine BE_dt_full_evolve_Krylov_exact_diag

end subroutine BE_dt_evolve_Houston_probe_decomp

subroutine BE_half_dt_evolution_only_with_intraband(act_1, act_2)
  use global_variables
  implicit none
  real(8),intent(in) :: act_1, act_2
  real(8) :: Act_tmp
  complex(8),allocatable :: zH_tot_1(:,:,:), zH_tot_2(:,:,:)
  complex(8),allocatable :: zvec(NB_basis, NB_TD)
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  real(8),allocatable :: w1(:),w2(:)
  integer :: info
  complex(8),allocatable :: zMat_diag_1(:,:),zMat_diag_2(:,:)
  
  lwork=6*NB_basis
  allocate(work_lp(lwork),rwork(3*NB_basis-2),w(NB_basis))
  allocate(w1(NB_basis),w2(NB_basis))
  allocate(zMat_diag_1(NB_basis,NB_basis))
  allocate(zMat_diag_2(NB_basis,NB_basis))
  allocate(zMat_tmp(NB_basis,NB_basis))
    
  allocate(zH_tot_1(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zH_tot_2(NB_basis,NB_basis,NK_s:NK_e))


!== act_1 start ==
  Act_tmp = act_1
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

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



  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot_1(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot_1 = zH_tot_1 + zH_loc + Act_tmp*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zH_tot_1(ib,ib,:) = zH_1_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do
!== act_1 end ==
!== act_2 start ==
  Act_tmp = act_2
  if(Act_tmp /= Act_tmp)then
    err_message='Act_t is NaN.'
    call err_finalize
  else if(abs(Act_tmp) > Amax)then
    err_message='Amax is too small.'
    call err_finalize
  end if

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



  xx = (Act_tmp-dble(iav_t)*dAmax)/dAmax
  zH_tot_2(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
    +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
    +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
  zH_tot_2 = zH_tot_2 + zH_loc + Act_tmp*(zPi_loc(:,:,:,1)*Epdir_1(1) &
                                     +zPi_loc(:,:,:,2)*Epdir_1(2) &
                                     +zPi_loc(:,:,:,3)*Epdir_1(3) )
  do ib = 1,NB_basis
    zH_tot_2(ib,ib,:) = zH_2_tot(ib,ib,:) + 0.5d0*Act_tmp**2
  end do
!== act_2 end ==


  do ik = nk_s, nk_e
    zMat_diag_1(:,:)=zH_tot_1(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag_1, NB_basis, w1, work_lp, lwork, rwork, info)

    zMat_diag_2(:,:)=zH_tot_2(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag_1, NB_basis, w2, work_lp, lwork, rwork, info)

    w = 0.5d0*(w1+w2)

    do ib = 1, nb_basis
      zs = sum(conjg(zMat_diag_2(:,ib))*zMat_diag_1(:,ib))
      phi = aimag(log(zs))
      zMat_diag_2(:,ib) = zMat_diag_2(:,ib)*exp(zi*phi)
    end do

! projection onto the basis 1
    do ib = 1, nb_td
      do ib2 = 1, nb_basis
        zvec(ib2,ib) = sum(conjg(zMat_diag_1(:,ib2))*zCt(:,ib,ik))
      end do
    end do

! propagation from t-> t+dt/2
    do ib2 = 1, nb_basis
      zvec(ib2,:) = zvec(ib2,:)*exp(-zi*dt*0.5d0*w(ib2))
    end do

! projection back from the basis 2
    zCt(:,:,ik) = matmul(zMat_diag_2(:,:),zvec(:,:))

  end do


end subroutine BE_half_dt_evolution_only_with_intraband
