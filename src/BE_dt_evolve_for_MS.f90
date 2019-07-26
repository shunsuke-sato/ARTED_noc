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
subroutine BE_dt_evolve_for_MS(Act_m_t)
  use global_variables
  use ms_maxwell_ks_variables
  use PSE_variables
  implicit none
  real(8),intent(in) :: Act_m_t(nx_s:nx_e)
  real(8) :: Act_t
  real(8),parameter :: eps_Act = 1d-6
  integer :: iav, iav_t
  real(8) :: diff,xx
  integer :: ik,ib,iexp
  complex(8) :: zfact
  integer :: ix_m

  M_point : do ix_m=Mx_s, Mx_e
    Act_t = Act_m_t(ix_m)

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

  call BE_dt_full_evolve_Taylor(ix_m)

  end do M_point


  return
  contains

    subroutine BE_dt_full_evolve_Taylor(ix_m)
      implicit none
      integer,intent(in) :: ix_m
      integer,parameter :: nexp_taylor = 4
      complex(8) :: zvec0(NB_basis, NB_TD)
      complex(8) :: zhvec0(NB_basis, NB_TD)
      
!$omp parallel do private(ik, iexp, zfact, zvec0, zhvec0)
      K_point : do ik=NK_s,NK_e

        zfact = 1d0
        zvec0 = zCt_Mpoint(1:NB_basis, 1:NB_TD, ik, ix_m)
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
            

          zCt_Mpoint(1:NB_basis, 1:NB_TD, ik, ix_m) = &
            zCt_Mpoint(1:NB_basis, 1:NB_TD, ik, ix_m) &
            +zfact*zhvec0(1:NB_basis, 1:NB_TD)
          if(iexp /= nexp_taylor)zvec0 = zhvec0

        end do
      end do K_point

    end subroutine BE_dt_full_evolve_Taylor

    subroutine BE_dt_full_evolve_Krylov_exact_diag(ix_m)
      implicit none
      integer,intent(in) :: ix_m
      integer,parameter :: nvec = 16
      complex(8) :: zvec(NB_basis, NB_TD, nvec)
      complex(8) :: zhvec(NB_basis, NB_TD, nvec)
      complex(8) :: zham_m(nvec, nvec),zUprop_m(nvec, nvec)
      integer :: ib, ivec, jvec
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

        zvec(1:NB_basis,1:NB_TD,1) = zCt_Mpoint(1:NB_basis, 1:NB_TD, ik, ix_m)
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
            zvec(:,:,ivec+1) = zhvec(:,:,ivec)
!Gram-Schmidt orthonormalization
            do ib = 1, nb_td
              do jvec = 1, ivec
                zs = sum(conjg(zvec(:,ib,jvec))*zvec(:,ib,ivec+1))
                zvec(:,ib,ivec+1) = zvec(:,ib,ivec+1) -zs*zvec(:,ib,jvec)
              end do
              ss = sum(abs(zvec(:,ib,ivec+1))**2); ss = 1d0/sqrt(ss)
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
          zCt_Mpoint(1:NB_basis, ib, ik, ix_m) = zvec(:,ib,1)*zUprop_m(1,1)
          do ivec = 2, nvec
            zCt_Mpoint(1:NB_basis, ib, ik, ix_m) = zvec(:,ib,ivec)*zUprop_m(ivec,1)
          end do

        end do


      end do K_point

    end subroutine BE_dt_full_evolve_Krylov_exact_diag


  end subroutine BE_dt_evolve_for_MS
