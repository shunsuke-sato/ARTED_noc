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
subroutine TD_Band_DOS(iter)
  use global_variables
  implicit none
  integer :: iter
  integer :: ik,ib,ibt
  character(100) :: file_name
  integer :: iter_gs
  real(8) :: s,nocc
  character(50) :: citer

  zu_GS(:,:,:)=zu_GS0(:,:,:)

  select case(method)
  case('PAW')
    err_message='In preparation PAW method' ; call err_finalize
    call err_finalize
  case('PSE')
    do iter_gs=1,3
      call PSE_Conjugate_Gradient(Ncg)
      call Gram_Schmidt
      call PSE_subspace_diag
    end do
    occ_TD_l=0d0
    do ik=NK_s,NK_e
      do ib=1,NB
        s=0d0
        do ibt=1,NB_TD
          s=s+occ(ibt,ik)*abs(sum(conjg(zu(:,ibt,ik))*zu_GS(:,ib,ik))*H123)**2
        end do
        occ_TD_l(ib,ik)=s
      end do
    end do
    call MPI_ALLREDUCE(occ_TD_l,occ_TD,NB*NK,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  case('CHK')
  case default
    err_message='Invalid option "mehtod"' ; call err_finalize
  end select



  if(myrank == 0)then
    write(citer,'(I7.7)')iter
    file_name=trim(SYSname)//'_TD_Band_DoS_raw_'//trim(citer)//'.out'
    open(10,file=trim(file_name))
    write(10,'(A)')'# Time-resolved Band and Density of states calculation'
    write(10,'(A)')"# write(10,'10(I7,2x)')NK,NB,NK1,NK2,NK3 "
    write(10,'(A)')"# write(10,'3e26.16e3')Acext(1,iter+1),Acext(2,iter+1),Acext(3,iter+1)"
    write(10,'(A)')'# do ik=1,NK'
    write(10,'(A)')'# do ib=1,NB'
    write(10,'(A)')"# write(10,'(2e26.16e3)')esp(ib,ik),occ(ib,ik)"
    write(10,'(A)')'#end do'
    write(10,'(A)')'#end do'
    write(10,'(5(I7,2x))')NK,NB,NK1,NK2,NK3
    write(10,'(3e26.16e3)')Acext_Cvec(1,iter+1),Acext_Cvec(2,iter+1),Acext_Cvec(3,iter+1)
    do ik=1,NK
      do ib=1,NB
        write(10,'(2e26.16e3)')esp(ib,ik),occ_TD(ib,ik)
      end do
    end do
    close(10)
  end if


  if(myrank == 0)then
    if(NB_TD < NB)then
      write(103,'(3e26.16e3)')Dt*dble(iter+1),sum(occ(:,:))-sum(occ_TD(1:NB_TD,:)),sum(occ_TD(NB_TD+1:NB,:))
    else
      write(103,'(2e26.16e3)')Dt*dble(iter+1),sum(occ(:,:))-sum(occ_TD(1:NB_TD,:))
    end if
  end if

  return
end subroutine TD_Band_DOS
