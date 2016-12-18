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
subroutine init_wf_basis_expansion
  use global_variables
  implicit none
  complex(8),allocatable :: zMat_diag(:,:)
  integer :: ik
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info

  if(myrank == 0)write(*,"(A)")"== Start: Initialization of wavefunctions."

  lwork=6*NB_basis
  allocate(work_lp(lwork),rwork(3*NB_basis-2),w(NB_basis))
  allocate(zMat_diag(NB_basis,NB_basis))

  allocate(zCt(NB_basis,NB_TD,NK_s:NK_e))
  allocate(ztCt_tmp(NB_basis),zACt_tmp(NB_basis))
  allocate(zLanCt(NB_basis,NB_TD,NLanczos))
  allocate(ztCt_Lan(NB_basis,NB_TD),zACt_Lan(NB_basis,NB_TD))

  do ik = NK_s,NK_e
    zMat_diag(:,:)=zH_loc(:,:,ik)+zV_NL(:,:,ik,0)
    call zheev('V', 'U', NB_basis, zMat_diag, NB_basis, w, work_lp, lwork, rwork, info)
    zCt(1:NB_basis,1:NB_TD,ik)=zMat_diag(1:NB_basis,1:NB_TD)
  end do

  if(myrank == 0)write(*,"(A)")"== End: Initialization of wavefunctions."

  return
end subroutine init_wf_basis_expansion
