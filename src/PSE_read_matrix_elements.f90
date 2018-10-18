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
subroutine PSE_read_matrix_elements
  use global_variables
  implicit none
  integer :: ik
  character(50) :: cik, filename

  if(myrank == 0)write(*,"(A)")"== Start reading matrix elements."
  if(myrank == 0)then
    open(200,file="matrix_element/basis_exp_basic.out",form='unformatted')
    read(200)NB_basis
    read(200)Amax,dAmax
    read(200)Epdir_1
    close(200)
  end if

  call MPI_BCAST(NB_basis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Amax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dAmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Epdir_1,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  if(abs(dble(NAmax) -Amax/dAmax) > 0.01d0)then
    err_message='NAmax is not consistent.'
    call err_finalize
  end if

  allocate(zH_loc(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zPi_loc(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zV_NL(NB_basis,NB_basis,NK_s:NK_e,-NAmax:NAmax))
  allocate(zPi_NL(NB_basis,NB_basis,NK_s:NK_e,-NAmax:NAmax))
  allocate(zH_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zPi_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zH0_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zdH_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(H0_eigval(NB_basis,NK_s:NK_e))
  do ik = NK_s,NK_e
    write(cik,"(I9.9)")ik
    filename="matrix_element/"//trim(cik)//"_matrix_elements.out"
    open(201,file=filename,form='unformatted')
    read(201)zH_loc(:,:,ik)
    read(201)zPi_loc(:,:,ik)
    read(201)zV_NL(:,:,ik,:)
    read(201)zPi_NL(:,:,ik,:)
    close(201)
  end do

  if(myrank == 0)write(*,"(A)")"== End reading matrix elements."
  return
end subroutine PSE_read_matrix_elements
