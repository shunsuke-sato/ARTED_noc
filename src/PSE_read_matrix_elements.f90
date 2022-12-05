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
  integer :: NK_s_comm, NK_e_comm, iproc
  complex(8),allocatable :: zH_loc_comm(:,:,:),zPi_loc_comm(:,:,:,:)
  complex(8),allocatable :: zV_NL_comm(:,:,:,:),zPi_NL_comm(:,:,:,:,:)
  integer :: ista(MPI_STATUS_SIZE)
  

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
  allocate(zPi_loc(NB_basis,NB_basis,NK_s:NK_e,3))
  allocate(zV_NL(NB_basis,NB_basis,NK_s:NK_e,-NAmax:NAmax))
  allocate(zPi_NL(NB_basis,NB_basis,NK_s:NK_e,3,-NAmax:NAmax))
  allocate(zH_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zPi_tot(NB_basis,NB_basis,NK_s:NK_e,3))
  allocate(zH0_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zdH_tot(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zH_tot2(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zdH_tot2(NB_basis,NB_basis,NK_s:NK_e))
  allocate(H0_eigval(NB_basis,NK_s:NK_e))

  if(myrank == 0)then
    filename="matrix_element/matrix_elements.out"
    open(201,file=filename,form='unformatted')
    do ik = NK_s,NK_e
      read(201)zH_loc(:,:,ik)
      read(201)zPi_loc(:,:,ik,:)
      read(201)zV_NL(:,:,ik,:)
      read(201)zPi_NL(:,:,ik,:,:)
    end do
  end if

  do iproc = 1, Nprocs-1
     if(myrank == iproc)then
        call MPI_Send(NK_s, 1, MPI_INTEGER, 0, iproc, MPI_COMM_WORLD, ierr)
        call MPI_Send(NK_e, 1, MPI_INTEGER, 0, iproc, MPI_COMM_WORLD, ierr)

        call MPI_Recv(zH_loc,size(zH_loc), MPI_DOUBLE_COMPLEX, 0, iproc, MPI_COMM_WORLD, ista, ierr)
        call MPI_Recv(zPi_loc,size(zPi_loc), MPI_DOUBLE_COMPLEX, 0, iproc, MPI_COMM_WORLD, ista, ierr)
        call MPI_Recv(zV_NL,size(zV_NL), MPI_DOUBLE_COMPLEX, 0, iproc, MPI_COMM_WORLD, ista, ierr)
        call MPI_Recv(zPi_NL,size(zPi_NL), MPI_DOUBLE_COMPLEX, 0, iproc, MPI_COMM_WORLD, ista, ierr)

     else if(myrank == 0)then
        call MPI_Recv(NK_s_comm, 1, MPI_INTEGER, iproc, iproc, MPI_COMM_WORLD, ista, ierr)
        call MPI_Recv(NK_e_comm, 1, MPI_INTEGER, iproc, iproc, MPI_COMM_WORLD, ista, ierr)

        allocate(zH_loc_comm(NB_basis,NB_basis,NK_s_comm:NK_e_comm))
        allocate(zPi_loc_comm(NB_basis,NB_basis,NK_s_comm:NK_e_comm,3))
        allocate(zV_NL_comm(NB_basis,NB_basis,NK_s_comm:NK_e_comm,-NAmax:NAmax))
        allocate(zPi_NL_comm(NB_basis,NB_basis,NK_s_comm:NK_e_comm,3,-NAmax:NAmax))        
        
        do ik = NK_s_comm,NK_e_comm
           read(201)zH_loc_comm(:,:,ik)
           read(201)zPi_loc_comm(:,:,ik,:)
           read(201)zV_NL_comm(:,:,ik,:)
           read(201)zPi_NL_comm(:,:,ik,:,:)
        end do
        
        call MPI_Send(zH_loc_comm,size(zH_loc_comm), MPI_DOUBLE_COMPLEX, iproc, iproc, MPI_COMM_WORLD, ierr)
        call MPI_Send(zPi_loc_comm,size(zPi_loc_comm), MPI_DOUBLE_COMPLEX, iproc, iproc, MPI_COMM_WORLD, ierr)   
        call MPI_Send(zV_NL_comm,size(zV_NL_comm), MPI_DOUBLE_COMPLEX, iproc, iproc, MPI_COMM_WORLD, ierr)
        call MPI_Send(zPi_NL_comm,size(zPi_NL_comm), MPI_DOUBLE_COMPLEX, iproc, iproc, MPI_COMM_WORLD, ierr)

        deallocate(zH_loc_comm,zPi_loc_comm,zV_NL_comm,zPi_NL_comm)
        
     end if
     
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end do
  
  if(myrank == 0)close(201)
  if(myrank == 0)write(*,"(A)")"== End reading matrix elements."
  return
end subroutine PSE_read_matrix_elements
