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
subroutine Band_DoS
  use global_variables
  implicit none
  integer :: ik,ib
  character(100) :: file_name

  if(myrank == 0)then
    file_name=trim(SYSname)//'_Band_DoS_raw.out'
    open(10,file=trim(file_name))
    write(10,'(A)')'# Density of states'
    write(10,'(A)')"# write(10,'2(I7,2x)')NK,NB,NK1,NK2,NK3 "
    write(10,'(A)')'# do ik=1,NK'
    write(10,'(A)')'# do ib=1,NB'
    write(10,'(A)')"# write(10,'(2e26.16e3)')esp(ib,ik),occ(ib,ik)"
    write(10,'(A)')'end do'
    write(10,'(A)')'end do'
    write(10,'(5(I7,2x))')NK,NB,NK1,NK2,NK3
    do ik=1,NK
      do ib=1,NB
        write(10,'(2e26.16e3)')esp(ib,ik),occ(ib,ik)
      end do
    end do

    close(10)
  end if


  return
end subroutine Band_DoS
