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
subroutine psi_rho(GS_RT) 
  use global_variables
  implicit none
  character(2) :: GS_RT
  integer :: ib,ik

  rho_e_l=0d0

  select case (GS_RT)
  case('GS')
    do ik=NK_s,NK_e
      do ib=1,NB
        rho_e_l(:)=rho_e_l(:)+occ(ib,ik)*abs(zu_GS(:,ib,ik))**2
      end do
    end do
  case('RT')
    do ik=NK_s,NK_e
      do ib=1,NB_TD
        rho_e_l(:)=rho_e_l(:)+occ(ib,ik)*abs(zu(:,ib,ik))**2
      end do
    end do
  case default
    err_message='Err in psi_rho';call err_finalize
  end select

  call MPI_ALLREDUCE(rho_e_l,rho_e,NL,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)


  return
end subroutine psi_rho
