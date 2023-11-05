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
subroutine blocking_matrix_for_Houston_probe_decomp
  use global_variables
  implicit none
  integer :: ib1,ib2
  integer,allocatable :: nflag_include(:)

  if(myrank == 0)write(*,"(A)")&
    "== Start: constructing Blocking Matrix for Houston probe decomposition."

  allocate(Mask_probe(NB_basis,NB_basis))
  Mask_probe = 1d0


!  do ib1 = 1,NB_basis
!    do ib2 = 1,NB_basis
!
!      if(ib1 > 5 .and. ib2 > 5 .and. ib1 /= ib2)then
!        Mask_probe(ib1,ib2)=0d0
!      end if
!
!    end do
!  end do



  if(myrank == 0)write(*,"(A)")&
    "== End: constructing Blocking Matrix for Houston probe decomposition."


  if(switch_pump_induced_interband_transition)then
    allocate(Mask_pump(NB_basis,NB_basis))
    allocate(nflag_include(nb_basis))

! 1:include, 2:exclude
    nflag_include(:) = 1


    mask_pump = 1d0
    do ib1 = 1, nb_basis
      do ib2 = 1, nb_basis
        if(ib1 /= ib2)then
          mask_pump(ib1,ib2)=nflag_include(ib1)*nflag_include(ib2)
        end if
      end do
    end do

  end if

  return
end subroutine blocking_matrix_for_Houston_probe_decomp
