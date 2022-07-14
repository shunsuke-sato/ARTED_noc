!
!  Copyright 2022 Shunsuke A. Sato
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
subroutine occupation_Fermi_Dirac_dist
  use global_variables
  implicit none
  real(8) :: beta, mu, mu_min, mu_max
  real(8) :: num_elec_t

  if(Telec < 0d0)return

  mu_max = maxval(esp)
  mu_min = minval(esp)

  beta = 1d0/(Telec/(2d0*Ry))


  if(myrank == 0)then

    do
      mu = 0.5d0*(mu_max + mu_min)
      occ = (2d0/NK)/(exp(beta*(esp-mu))+1d0)
      num_elec_t = sum(occ)
      if(num_elec_t > real(Nelec))then
        mu_max = mu
      else
        mu_min = mu
      end if

      if(mu_max - mu_min < 1d-10)exit

    end do

    mu = 0.5d0*(mu_max + mu_min)
    occ = (2d0/NK)/(exp(beta*(esp-mu))+1d0)


  end if

  call MPI_BCAST(occ,NB*NK,MPI_REAL8,0,MPI_COMM_WORLD,ierr)  

  return
end subroutine occupation_Fermi_Dirac_dist
