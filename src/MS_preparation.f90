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
  integer :: ik
  character(50) :: cik, filename

  write(*,*)myrank
! Parameter for maxwell-kohn-sham
  Nx_s = -2000
  Nx_e =  2000
  Mx = 160
  dx_m = 40d-9/0.529177d-10/Mx

  if(mod(Nprocs,Mx) /=0)stop 'Error Mx is not dividable by Nprocs'

  nprocs_per_Mpoint = Nprocs/Mx
  do ix = 1, Mx
    if(nprocs_per_Mpoint*ix > myrank)then
      macro_point_id = ix
      exit
    end if
  end do

  write(*,*)"myrak, macro_point_id",myrank,macro_point_id


  NK_ave=NK/Nprocs; NK_remainder=mod(NK,Nprocs)
  if(Myrank < NK_remainder)then
    NK_s=(NK_ave+1)*Myrank+1
    NK_e=(NK_ave+1)*Myrank+(NK_ave+1)
  else if(Myrank >= NK_remainder)then
    NK_s=(Myrank-NK_remainder)*NK_ave+(NK_ave+1)*NK_remainder+1
    NK_e=(Myrank-NK_remainder)*NK_ave+(NK_ave+1)*NK_remainder+NK_ave
  end if

  if(Myrank == 0)write(*,'(A)')'NK-split is completed'
  if(Myrank == 0)write(*,'(2(A,2x,I0,2x))')'NK_ave =',NK_ave,'NK_remainder =',NK_remainder


  if(myrank == 0)then
    open(200,file="basis_exp_basic.out",form='unformatted')
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
    filename=trim(cik)//"_matrix_elements.out"
    open(201,file=filename,form='unformatted')
    read(201)zH_loc(:,:,ik)
    read(201)zPi_loc(:,:,ik)
    read(201)zV_NL(:,:,ik,:)
    read(201)zPi_NL(:,:,ik,:)
    close(201)
  end do

  if(myrank == 0)write(*,"(A)")"== End reading matrix elements."


end subroutine MS_preparation
