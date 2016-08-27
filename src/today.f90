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
subroutine today
  implicit none
  character(8) :: date
  character(10) :: time
  character(5) :: zone 
  integer :: values(8)

  call date_and_time( date , time, zone, values )
  write(*,'(4(2x,A))')'day =',date,'area =',zone

  return
end subroutine today
