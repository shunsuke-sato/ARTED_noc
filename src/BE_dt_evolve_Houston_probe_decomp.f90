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
  integer,parameter :: NTaylor = 4
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact
  real(8) :: Act_tmp,Act_probe_tmp
  complex(8) :: zvec_t(NB_basis,NB_TD)
  integer :: ierr_lan
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

!==construct Pump hamiltonian ==

  Act_tmp = 0.5d0*(Ac_pump_BE(iter) + Ac_pump_BE(iter+1) )
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
  zdH_tot = zdH_tot + zH_loc + zPi_loc*Act_t
  do ib = 1,NB_basis
    zdH_tot(ib,ib,:) = zdH_tot(ib,ib,:) + 0.5d0*Act_t**2
  end do

! construct probe Hamiltonian
  zdH_tot = zdH_tot - zH_tot


!== construct current matrix end

  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_tot(:,:,ik)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zPi_tot(:,:,ik) = matmul(zdH_tot(:,:,ik),zMat_diag(:,:))
    zdH_tot(:,:,ik) = matmul(transpose(conjg(zMat_diag(:,:))),zPi_tot(:,:,ik))
    zdH_tot(:,:,ik) = zdH_tot(:,:,ik)*Mask_probe(:,:)

    zPi_tot(:,:,ik) = matmul(zdH_tot(:,:,ik),transpose(conjg(zMat_diag(:,:))))
    zdH_tot(:,:,ik) = matmul(zMat_diag(:,:),zPi_tot(:,:,ik))

    zH_tot(:,:,ik) = zH_tot(:,:,ik) + zdH_tot(:,:,ik)
  end do


  if(abs(Act_t) < eps_Act)then
    zdH_tot(:,:,:) = zH_tot(:,:,:) - (zH_loc(:,:,:) + zV_NL(:,:,:,0))
    call BE_dt_half_evolve_Taylor
!  call BE_dt_half_evolve_Lanczos
    call BE_dt_evolve_Free
    call BE_dt_half_evolve_Taylor
!  call BE_dt_half_evolve_Lanczos
  else
    zCt_tmp = zCt
    call BE_dt_full_evolve_Lanczos(ierr_lan)
    if(ierr_lan /= 0)then
      zCt = zCt_tmp
      zdH_tot(:,:,:) = zH_tot(:,:,:) - (zH_loc(:,:,:) + zV_NL(:,:,:,0))
      call BE_dt_half_evolve_Taylor
      call BE_dt_evolve_Free
      call BE_dt_half_evolve_Taylor

    end if
  end if

  return
end subroutine BE_dt_evolve_Houston_probe_decomp

subroutine BE_dt_evolve_Free
  use global_variables
  implicit none
  complex(8) :: zvec0(NB_basis)
  integer ik, ib

  K_point : do ik=NK_s,NK_e
    Band : do ib=1,NB_TD

      zvec0 = matmul(transpose(conjg(zC_eig(:,:,ik))),zCt(:,ib,ik))
      zvec0(:) = exp(-zI*H0_eigval(:,ik)*dt)*zvec0(:)
      zCt(:,ib,ik) = matmul(zC_eig(:,:,ik),zvec0(:))

    end do Band
  end do K_point
  
  return
end subroutine BE_dt_evolve_Free


subroutine BE_dt_half_evolve_Taylor
  use global_variables
  implicit none
  integer :: ik,ib,iexp
  complex(8) :: zfact

  K_point : do ik=NK_s,NK_e
    Band : do ib=1,NB_TD

      zfact = 1d0
      ztCt_tmp(:)=zCt(:,ib,ik)
      do iexp = 1,8
        zfact = zfact*((-zI*dt)*0.5d0)/dble(iexp)
        zACt_tmp(:) = matmul(zdH_tot(:,:,ik),ztCt_tmp(:))
        zCt(:,ib,ik) = zCt(:,ib,ik) + zfact*zACt_tmp(:)
        ztCt_tmp(:) = zACt_tmp(:)
      end do

    end do Band
  end do K_point
end subroutine BE_dt_half_evolve_Taylor


subroutine BE_dt_full_evolve_Lanczos(ierr_lan)
  use global_variables
  implicit none
  integer :: iLan
  real(8) :: alpha(NLanczos,NB_TD),beta(NLanczos,NB_TD)
  real(8) :: Hmat_Lan(NLanczos,NLanczos,NB_TD)
  complex(8),allocatable :: zvec(:)
  real(8) :: ss
  integer :: ik,ib
  integer :: ierr_lan,ierr_lan_l
!LAPACK ==
  integer :: lwork,Nmat
  real(8),allocatable :: work_lp(:),Amat(:,:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
!LAPACK ==

  ierr_lan_l = 0

  K_point : do ik=NK_s,NK_e

    Nmat = 0
    do ib=1,NB_TD        
      ss = sum(abs(zCt(:,ib,ik))**2)
      zLanCt(:,ib,1) = zCt(:,ib,ik)/sqrt(ss)
    end do
    beta(1,:) = 0d0

    do iLan=1,NLanczos-1
      ztCt_Lan(:,:) = zLanCt(:,:,iLan)
      zACt_Lan(:,:) = matmul(zH_tot(:,:,ik),ztCt_Lan(:,:))

      do ib=1,NB_TD
        alpha(iLan,ib)=sum(conjg(ztCt_Lan(:,ib))*zACt_Lan(:,ib) )
      end do
      
      if(iLan /= 1)then
        do ib=1,NB_TD
          zLanCt(:,ib,iLan+1) = zAct_Lan(:,ib) &
            -alpha(iLan,ib)*zLanCt(:,ib,iLan) &
            -beta(iLan,ib)*zLanCt(:,ib,iLan-1)
        end do
      else
        do ib=1,NB_TD
          zLanCt(:,ib,iLan+1) = zAct_Lan(:,ib) &
            -alpha(iLan,ib)*zLanCt(:,ib,iLan) 
        end do
      end if
      
      do ib=1,NB_TD
        beta(iLan+1,ib)=sqrt(sum(abs(zLanCt(:,ib,iLan+1))**2)  )
      end do
      if(minval(beta(iLan+1,:)) < epsilon_Lan)then
        Nmat = iLan
        write(*,"(A,2x,3I7)")"Lanzcos break, Nmat,ik,myrank",Nmat,ik,myrank
        ierr_lan_l = 1
        exit
      end if
      do ib=1,NB_TD
        zLanCt(:,ib,iLan+1) = zLanCt(:,ib,iLan+1)/beta(iLan+1,ib)
      end do
      
    end do
    
    if(Nmat == 0)then
      Nmat = Nlanczos
      ztCt_Lan(:,:) = zLanCt(:,:,Nmat)
      zACt_Lan(:,:) = matmul(zH_tot(:,:,ik),ztCt_Lan(:,:))
      do ib=1,NB_TD
        alpha(Nmat,ib)=sum(conjg(ztCt_Lan(:,ib))*zACt_Lan(:,ib) )
      end do
    end if

    lwork=6*(Nmat)**2
    allocate(work_lp(lwork),rwork(3*Nmat-2),w(Nmat),Amat(Nmat,Nmat),zvec(Nmat))
    
    do ib=1,NB_TD
      Amat = 0d0
      do iLan=1,Nmat
        Amat(iLan,iLan)=alpha(iLan,ib)
      end do
      do iLan = 2,Nmat
        Amat(iLan-1,iLan) = beta(iLan,ib)
        Amat(iLan,iLan-1) = beta(iLan,ib)
      end do
      
      Call dsyev('V', 'U', Nmat, Amat, Nmat, w, work_lp, lwork, info) 
      zvec = 0d0; zvec(1) = 1d0
      zvec = matmul(transpose(Amat),zvec)
      zvec(:) = exp(-zI*dt*w(:))*zvec(:)
      zvec = matmul(Amat,zvec)
      zCt(:,ib,ik)=matmul(zLanCt(:,ib,1:Nmat),zvec(1:Nmat))
    end do
    
    deallocate(work_lp,rwork,w,Amat,zvec)

    call MPI_ALLREDUCE(ierr_lan_l,ierr_lan,1,MPI_INTEGER,MPI_SUM,NEW_COMM_WORLD,ierr)      

    
  end do K_point
end subroutine BE_dt_full_evolve_Lanczos

  

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
  
