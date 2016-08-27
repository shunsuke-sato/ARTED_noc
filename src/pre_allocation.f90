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
subroutine pre_allocation
  use global_variables
  implicit none

  NB_TD=Nelec/2

  allocate(zu(NL,NB_TD,NK_s:NK_e),zu_GS(NL,NB,NK_s:NK_e),zu_GS0(NL,NB,NK_s:NK_e))
  allocate(tpsi(NL),htpsi(NL),tpsi_g1(NL),tpsi_g2(NL))
  allocate(rho_c(NL),rho_p(NL),rho_e(NL),rho_e_l(NL))
  allocate(rho_c_3D(0:NL1-1,0:NL2-1,0:NL3-1))
  allocate(Vloc(NL),Vh(NL),Vxc(NL),Exc(NL))
  allocate(occ(NB,NK),esp(NB,NK),esp_l(NB,NK))
  allocate(zu_t(NL,NB_TD,NK_s:NK_e),Vloc_t(NL))
  allocate(occ_TD(NB,NK),occ_TD_l(NB,NK))


  occ=0d0; occ(1:Nelec/2,:)=2d0/dble(NK)

  if(Myrank == 0)write(*,'(A)')'pre allocations is completed'

  return
end subroutine pre_allocation
