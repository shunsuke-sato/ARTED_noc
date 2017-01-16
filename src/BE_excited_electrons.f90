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
subroutine BE_excited_electrons
  use global_variables
  implicit none
  real(8),allocatable :: occ_t(:,:,:),occ_t_l(:,:,:)
  integer :: ik,ib1,ib2

  if(myrank == 0)write(*,"(A)")"== Start: Distribution of electron-hole."

  allocate(occ_t(NB_basis,NB_TD,NK),occ_t_l(NB_basis,NB_TD,NK))
  occ_t_l = 0d0

  do ik = NK_s,NK_e
    do ib1 = 1,NB_TD
      do ib2 = 1,NB_basis

        occ_t_l(ib2,ib1,ik) = abs(sum(conjg(zC_eig(:,ib2,ik))*zCt(:,ib1,ik)))**2/dble(NK)

      end do
    end do
  end do

  call MPI_ALLREDUCE(occ_t_l,occ_t,NB_basis*NB_TD*NK,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)  

  if(myrank == 0)then
    open(701,file="dist_nex.out")
    write(701,"(A,2x,3I7)")"#NB_basis,NB_TD,NK=",NB_basis,NB_TD,NK
    do ik = NK_s,NK_e
      do ib1 = 1,NB_TD
        do ib2 = 1,NB_basis
          write(701,"(999e26.16e3)")occ_t(ib2,ib1,ik)
        end do
      end do
    end do
    close(701)
  end if

  if(myrank == 0)write(*,"(A)")"== Start: Distribution of electron-hole."

  return
end subroutine BE_excited_electrons
