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
        call zhemm('L', 'U', NB_basis, NB_TD, (1d0,0d0), &
          zH_tot(1:NB_basis,1:NB_basis,ik), &
          NB_basis,&
          zvec0(1:NB_basis,1:NB_TD), &
          NB_basis, (0d0,0d0), &
          zhvec0(1:NB_basis,1:NB_TD), &
          NB_basis)
            

        zCt(1:NB_basis, 1:NB_TD, ik) = zCt(1:NB_basis, 1:NB_TD, ik) &
          +zfact*zhvec0(1:NB_basis, 1:NB_TD)
        if(iexp /= nexp_taylor)zvec0 = zhvec0

      end do
    end do K_point

  end subroutine BE_dt_full_evolve_Taylor

end subroutine BE_dt_evolve
!---------------------------------------------------------------------------------
subroutine BE_dt_evolve_old(Act_t)
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: Act_t
  real(8),parameter :: eps_Act = 1d-6
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact

  if(Act_t == 0d0)then
    ! free propagation
  call BE_dt_evolve_Free
    return
  end if

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

  if(abs(Act_t) < eps_Act)then
    xx = (Act_t-dble(iav_t)*dAmax)/dAmax
    if(iav_t == 0)then
      zdH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx+1d0) &
        +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx-1d0) &
        +      zV_NL(:,:,:,iav_t)*( - xx)
      zdH_tot(:,:,:) = zdH_tot(:,:,:)/dAmax

    else

      zdH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
        +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
        +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
      zdH_tot(:,:,:) = (zdH_tot(:,:,:) - zV_NL(:,:,:,0))/Act_t
    end if
    zdH_tot = zdH_tot + zPi_loc
    do ib = 1,NB_basis
      zdH_tot(ib,ib,:) = zdH_tot(ib,ib,:) + 0.5d0*Act_t
    end do
!== construct current matrix end

    call BE_dt_half_evolve_Taylor
!  call BE_dt_half_evolve_Lanczos
    call BE_dt_evolve_Free
    call BE_dt_half_evolve_Taylor
!  call BE_dt_half_evolve_Lanczos
  else

    xx = (Act_t-dble(iav_t)*dAmax)/dAmax
    zH_tot(:,:,:) = 0.5d0*zV_NL(:,:,:,iav_t+1)*(xx**2+xx) &
      +0.5d0*zV_NL(:,:,:,iav_t-1)*(xx**2-xx) &
      +      zV_NL(:,:,:,iav_t)*(1d0 - xx**2)
    zH_tot = zH_tot + zH_loc + zPi_loc*Act_t
    do ib = 1,NB_basis
      zH_tot(ib,ib,:) = zH_tot(ib,ib,:) + 0.5d0*Act_t**2
    end do

    call BE_dt_full_evolve_Lanczos
  end if

  return
  contains

    subroutine BE_dt_evolve_Free
      implicit none
      complex(8) :: zvec0(NB_basis)

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
      implicit none

      K_point : do ik=NK_s,NK_e
        Band : do ib=1,NB_TD

          zfact = 1d0
          ztCt_tmp(:)=zCt(:,ib,ik)
          do iexp = 1,8
            zfact = zfact*(-zI*(dt*Act_t)*0.5d0)/dble(iexp)
            zACt_tmp(:) = matmul(zdH_tot(:,:,ik),ztCt_tmp(:))
            zCt(:,ib,ik) = zCt(:,ib,ik) + zfact*zACt_tmp(:)
            ztCt_tmp(:) = zACt_tmp(:)
          end do

        end do Band
      end do K_point
    end subroutine BE_dt_half_evolve_Taylor

    subroutine BE_dt_half_evolve_Lanczos
      implicit none
      integer :: iLan
      real(8) :: alpha(NLanczos,NB_TD),beta(NLanczos,NB_TD)
      real(8) :: Hmat_Lan(NLanczos,NLanczos,NB_TD)
      complex(8),allocatable :: zvec(:)
      real(8) :: ss
!LAPACK ==
      integer :: lwork,Nmat
      real(8),allocatable :: work_lp(:),Amat(:,:)
      real(8),allocatable :: rwork(:),w(:)
      integer :: info
!LAPACK ==



      K_point : do ik=NK_s,NK_e

        Nmat = 0
        do ib=1,NB_TD        
          ss = sum(abs(zCt(:,ib,ik))**2)
          zLanCt(:,ib,1) = zCt(:,ib,ik)/sqrt(ss)
        end do
        beta(1,:) = 0d0

        do iLan=1,NLanczos-1
          ztCt_Lan(:,:) = zLanCt(:,:,iLan)
          zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))

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
            exit
          end if
          do ib=1,NB_TD
            zLanCt(:,ib,iLan+1) = zLanCt(:,ib,iLan+1)/beta(iLan+1,ib)
          end do

        end do

        if(Nmat == 0)then
          Nmat = Nlanczos
          ztCt_Lan(:,:) = zLanCt(:,:,Nmat)
          zACt_Lan(:,:) = matmul(zdH_tot(:,:,ik),ztCt_Lan(:,:))
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
          zvec(:) = exp(-zI*(dt*Act_t)*0.5d0*w(:))*zvec(:)
          zvec = matmul(Amat,zvec)
          zCt(:,ib,ik)=matmul(zLanCt(:,ib,1:Nmat),zvec(1:Nmat))
        end do

        deallocate(work_lp,rwork,w,Amat,zvec)

      end do K_point
    end subroutine BE_dt_half_evolve_Lanczos


    subroutine BE_dt_full_evolve_Lanczos
      implicit none
      integer :: iLan
      real(8) :: alpha(NLanczos,NB_TD),beta(NLanczos,NB_TD)
      real(8) :: Hmat_Lan(NLanczos,NLanczos,NB_TD)
      complex(8),allocatable :: zvec(:)
      real(8) :: ss
!LAPACK ==
      integer :: lwork,Nmat
      real(8),allocatable :: work_lp(:),Amat(:,:)
      real(8),allocatable :: rwork(:),w(:)
      integer :: info
!LAPACK ==



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

      end do K_point
    end subroutine BE_dt_full_evolve_Lanczos

  end subroutine BE_dt_evolve_old
  
