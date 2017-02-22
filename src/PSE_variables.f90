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
module PSE_variables
  implicit none
! pseudopotential
  integer,parameter :: Nrmax=3000,Lmax=4
  character(2) :: ps_type
  character(10) :: ps_format 
  integer :: Nps,Nlma
  integer,allocatable :: Mps(:),Jxyz(:,:),Jxx(:,:),Jyy(:,:),Jzz(:,:)
  integer,allocatable :: Mlps(:),Lref(:),Zps(:),NRloc(:)
  integer,allocatable :: NRps(:),inorm(:,:),iuV(:),a_tbl(:)
  real(8),allocatable :: rad(:,:),Rps(:),vloctbl(:,:),udVtbl(:,:,:)
  real(8),allocatable :: radnl(:,:)
  real(8),allocatable :: Rloc(:),uV(:,:),duV(:,:,:),anorm(:,:)
  real(8),allocatable :: dvloctbl(:,:),dudVtbl(:,:,:)
  real(8),allocatable :: rho_nlcc_tbl(:,:),tau_nlcc_tbl(:,:)
  real(8),allocatable :: rho_nlcc(:),tau_nlcc(:)
  complex(8),allocatable :: ekr(:,:,:)
  real(8),allocatable :: Vpsl(:)
  logical,allocatable :: flag_nlcc(:)

end Module PSE_variables



