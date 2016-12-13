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

! Copyright (C) 2013 Shunsuke A. Sato.
! Device number 10-99: temporary files, 100: permanent files
! 102:trim(SYSname)//'_jac.out'
! 103:trim(SYSname)//'_nex.out'

program main
  use global_variables
  implicit none

  call preparation

  select case(calc_mode)
  case('RT') ! real-time propagation
    call PSE_ground_state_calculation
    call PSE_real_time_propagation
  case('BD') ! band structure
    call PSE_ground_state_calculation
    call PSE_band_calculation
  case('GS') ! ground state
    call PSE_ground_state_calculation
    if(myrank == 0)then
      open(99,file="Vloc_gs.out",form='unformatted')
      write(99)Vloc
      close(99)
    end if
  case('MS') ! Matrix element preparation
    if(myrank == 0)then
      open(99,file="Vloc_gs.out",form='unformatted')
      read(99)Vloc
      close(99)
    end if
    call MPI_BCAST(Vloc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call PSE_preparation_basis_set
    call PSE_preparation_matrix

  case('MT') ! Time-propagation with basis expansion

    call PSE_read_matrix_elements

  case default
    err_message='invalid calc_mode'
    call err_finalize
  end select

  call MPI_FINALIZE(ierr)
end program main

