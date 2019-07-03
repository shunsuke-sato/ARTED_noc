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
subroutine MS_preparation
  use global_variables
  use ms_maxwell_ks_variables
  implicit none
  integer :: ix
  integer :: ik, irank_m, i
  character(50) :: cik, filename

  write(*,*)myrank

  call init_grid
  deallocate(Lx, iLx, iLx123)

! Parameter for maxwell-kohn-sham
  Mx = 224
  dx_m = 50d-9/0.529177d-10/Mx
  Nx_s = -(15d-6/0.529177d-10)/dx_m
  Nx_e = +(15d-6/0.529177d-10)/dx_m


  if(mod(Nprocs,Mx) /=0)stop 'Error Mx is not dividable by Nprocs'

  nprocs_per_Mpoint = Nprocs/Mx
  do ix = 1, Mx
    if(nprocs_per_Mpoint*ix > myrank)then
      macro_point_id = ix
      irank_m = myrank - nprocs_per_Mpoint*(ix-1)
      exit
    end if
  end do

  write(*,*)"myrak, macro_point_id,irank",myrank,macro_point_id,irank_m


  NK_ave=NK/Nprocs_per_Mpoint; NK_remainder=mod(NK,Nprocs_per_Mpoint)
  if(irank_m < NK_remainder)then
    NK_s=(NK_ave+1)*irank_m+1
    NK_e=(NK_ave+1)*irank_m+(NK_ave+1)
  else if(irank_m >= NK_remainder)then
    NK_s=(irank_m-NK_remainder)*NK_ave+(NK_ave+1)*NK_remainder+1
    NK_e=(irank_m-NK_remainder)*NK_ave+(NK_ave+1)*NK_remainder+NK_ave
  end if

  if(Myrank == 0)write(*,'(A)')'NK-split is completed'
  if(Myrank == 0)write(*,'(2(A,2x,I0,2x))')'NK_ave =',NK_ave,'NK_remainder =',NK_remainder


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


end subroutine MS_preparation
