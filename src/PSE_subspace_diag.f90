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
subroutine PSE_subspace_diag
  use global_variables
  use PSE_variables
  implicit none
  integer:: ik,ib1,ib2
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info


  lwork=6*NB
  allocate(work_lp(lwork),rwork(3*NB-2),w(NB))

  esp_l=0d0

  call PSE_pre_hpsi

  do ik=NK_s,NK_e
    do ib1=1,NB
      tpsi(1:NL)=zu_GS(1:NL,ib1,ik)
      call PSE_hpsi(ik)
      do ib2=ib1+1,NB
        za_diag(ib2,ib1)=sum(conjg(zu_GS(:,ib2,ik))*htpsi(:))*H123
        za_diag(ib1,ib2)=conjg(za_diag(ib2,ib1))
      enddo
      za_diag(ib1,ib1)=real(sum(conjg(zu_GS(:,ib1,ik))*htpsi(:))*H123)
    enddo
    call zheev('V', 'U', NB, za_diag, NB, w, work_lp, lwork, rwork, info)

    zutmp_diag=0.d0
    do ib1=1,NB
      do ib2=1,NB
        zutmp_diag(:,ib1)=zutmp_diag(:,ib1)+zu_GS(:,ib2,ik)*za_diag(ib2,ib1)
      enddo
    enddo
    zu_GS(:,:,ik)=zutmp_diag(:,:)
    esp_l(:,ik)=w(:)
  enddo


  call MPI_ALLREDUCE(esp_l,esp,NB*NK,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  return
end subroutine PSE_subspace_diag
  
