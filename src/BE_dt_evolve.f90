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

!  call BE_dt_evolve_Taylor
  call BE_dt_evolve_Lanczos

  return
  contains
    subroutine BE_dt_evolve_Taylor
      implicit none

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
    end subroutine BE_dt_evolve_Taylor

    subroutine BE_dt_evolve_Lanczos
      implicit none
      integer :: iLan
      real(8) :: alpha(NLanczos,NB_TD),beta(NLanczos,NB_TD)
      real(8) :: Hmat_Lan(NLanczos,NLanczos,NB_TD)
      complex(8) :: zvec(Nlanczos)
      real(8) :: ss
!LAPACK ==
      integer :: lwork,Nmat
      real(8),allocatable :: work_lp(:),Amat(:,:)
      real(8),allocatable :: rwork(:),w(:)
      integer :: info
      Nmat = Nlanczos
      lwork=6*(Nmat)**2
      allocate(work_lp(lwork),rwork(3*Nmat-2),w(Nmat),Amat(Nmat,Nmat))
!LAPACK ==

      K_point : do ik=NK_s,NK_e

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

          do ib=1,NB_TD
            zLanCt(:,ib,iLan+1) = zAct_Lan(:,ib) &
              -alpha(iLan,ib)*zLanCt(:,ib,iLan) &
              -beta(iLan,ib)*zLanCt(:,ib,iLan-1)
          end do

          do ib=1,NB_TD
            beta(iLan+1,ib)=sqrt(sum(abs(zLanCt(:,ib,iLan+1))**2)  )
            zLanCt(:,ib,iLan+1) = zLanCt(:,ib,iLan+1)/beta(iLan+1,ib)
          end do

        end do

        ztCt_Lan(:,:) = zLanCt(:,:,NLanczos)
        zACt_Lan(:,:) = matmul(zH_tot(:,:,ik),ztCt_Lan(:,:))
        do ib=1,NB_TD
          alpha(NLanczos,ib)=sum(conjg(ztCt_Lan(:,ib))*zACt_Lan(:,ib) )
        end do

        do ib=1,NB_TD
          Amat = 0d0
          do iLan=1,NLanczos
            Amat(iLan,iLan)=alpha(iLan,ib)
          end do
          do iLan = 2,Nlanczos
            Amat(iLan-1,iLan) = beta(iLan,ib)
            Amat(iLan,iLan-1) = beta(iLan,ib)
          end do

          Call dsyev('V', 'U', Nmat, Amat, Nmat, w, work_lp, lwork, info) 
          zvec = 0d0; zvec(1) = 1d0
          zvec = matmul(transpose(Amat),zvec)
          zvec(:) = exp(-zI*dt*w(:))*zvec(:)
          zvec = matmul(Amat,zvec)
          zCt(:,ib,ik)=matmul(zLanCt(:,ib,1:NLanczos),zvec(1:NLanczos))
        end do


      end do K_point
    end subroutine BE_dt_evolve_Lanczos

end subroutine BE_dt_evolve
  
