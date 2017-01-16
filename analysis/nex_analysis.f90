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
program nex_analysis
  implicit none
  real(8),parameter :: dw=0.1d0/27.2d0 ! Energy resolution [eV/27.2 =a.u.]
  real(8),parameter :: Pi=3.141592653589793d0
  real(8),allocatable :: esp(:,:),occ(:,:)
  real(8),allocatable :: DoS(:),DoS_occ(:)
  real(8) :: esp_min,esp_max,esp_d
  real(8) :: wi,wf,w
  integer :: Nw,iw
  integer :: ik,ib,NK,NB

  read(*,*)
  read(*,*)NB,NK
  allocate(esp(NB,NK),occ(NB,NK))

  do ik=1,NK
    do ib=1,NB
      read(*,*)esp(ib,ik),occ(ib,ik)
    end do
  end do

  esp_min=minval(esp)
  esp_max=maxval(esp)

  esp_d=(esp_max-esp_min)*1.0d0
  wi=esp_min-esp_d*0.1
  wf=esp_max+esp_d*0.1

  Nw=aint((wf-wi)/(dw/10d0))+1
  allocate(DoS(0:Nw),DoS_occ(0:Nw))
  DoS=0d0
  DoS_occ=0d0

  do ik=1,NK
    do ib=1,NB
      do iw=0,Nw
        w=wi+dble(iw)*(dw/10d0)
        DoS(iw)=DoS(iw)+(2d0/dble(NK))*(1d0/sqrt(2d0*pi*dw**2))*exp(-0.5d0*((w-esp(ib,ik))/dw)**2)
        DoS_occ(iw)=DoS_occ(iw)+occ(ib,ik)*(1d0/sqrt(2d0*pi*dw**2))*exp(-0.5d0*((w-esp(ib,ik))/dw)**2)
      end do
    end do
  end do

  
  do iw=0,Nw
    w=wi+dble(iw)*(dw/10d0)
    write(*,'(100e26.16e3)')w,DoS(iw),DoS_occ(iw)
  end do

end program nex_analysis
  

