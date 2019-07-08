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
! Cvec: cartesian 
! Lvec: Lattice
! Rvec: Reciprocal lattice
module ms_maxwell_ks_variables
  implicit none
! grid
  integer :: Nx_s, Nx_e
  integer :: Mx
  real(8) :: dx_m
  real(8),allocatable :: x_m(:)

! parallelization
  integer :: Mpoints_per_procs
  integer :: nprocs_per_Mpoint_group
  integer :: id_Mpoint_group
  integer :: Mx_s, Mx_e

! orbital
  integer :: NK_s_m, NK_e_m

! electromagnetic fields
  real(8),allocatable :: Ac_m(:), jt_m(:)
  real(8),allocatable :: Ac_m_n(:),Ac_m_o(:)




end module ms_maxwell_ks_variables
