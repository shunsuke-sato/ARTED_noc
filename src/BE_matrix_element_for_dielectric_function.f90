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
subroutine BE_matrix_element_for_dielectric_function
! This routine should be merged to the initial preparation of the matrix elements.
  use global_variables
  implicit none
  real(8),allocatable :: Pz2(:,:)
  integer :: ik,ib1,ib2
  character(99) :: cik,filename

  if(myrank == 0)write(*,"(A)")"== Start: Matrix-element for dielectric function."

  allocate(Pz2(NB_basis,NB_basis))
  zPi_tot(:,:,:) =  zPi_loc(:,:,:) + zPi_NL(:,:,:,0)

  do ik = NK_s,NK_e

    do ib1=1,NB_basis
      do ib2=ib1,NB_basis
        zACt_tmp(:) = matmul(zPi_tot(:,:,ik),zC_eig(:,ib2,ik))
        Pz2(ib1,ib2) = abs(sum(conjg(zC_eig(:,ib2,ik))*zACt_tmp(:)))**2
        Pz2(ib2,ib1) = Pz2(ib1,ib2)
      end do
    end do


    write(cik,"(I9.9)")ik
    filename=trim(cik)//"_matrix_elements_for_eps.out"
    open(201,file=filename,form='unformatted')
    write(201)Pz2
    close(201)
  end do

  if(myrank == 0)write(*,"(A)")"== End: Matrix-element for dielectric function."

  return
end subroutine BE_matrix_element_for_dielectric_function
