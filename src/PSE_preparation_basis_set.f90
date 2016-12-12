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
subroutine PSE_preparation_basis_set
  use global_variables
  implicit none
  integer :: iter_scf,ib
  real(8) :: rmix=0.2d0
  real(8) :: jav(3)
  real(8) :: Etot,Ekin
  real(8),allocatable :: Vloc_old(:)
  character(10) :: cEex_Cor_tmp
  integer :: NBs,NBe,ib1,ib2,ik,ikshift
  real(8) :: ss
  complex(8) :: zs

  if(Myrank == 0)write(*,'(A)')'Preparation matrix-element start'
  if(myrank == 0)write(*,*)'charge =',sum(rho_c)*H123,sum(rho_e)*H123,sum(rho_p)*H123  

  NK_shift = 2
  NB_basis_main = 14
  NB_basis_shift = 14
  NB_basis = NB_basis_main + NK_shift*NB_basis_shift
  allocate(kshift(3,NK_shift))
  kshift(1,1)=0.05d0/sqrt(2d0); kshift(2,1)=0.05d0/sqrt(2d0); kshift(3,1)=0d0
  kshift(:,2)=-kshift(:,1)
  if(NB < max(NB_basis_main,NB_basis_shift))stop "NB is too small."
  allocate(zu_basis(NL,NB_basis,NK_s:NK_e))


  do iter_scf=1,Nscf
    if(myrank == 0)write(*,*)kAc_Cvec(:,1)
    if(Myrank == 0)write(*,'(A,2X,I5)')'iter_scf = ',iter_scf

    rho_MB_in(:,min(iter_scf,MaxMem_MB)) = rho_e(:) ! For Broiden's method
    
    call Gram_Schmidt
    call PSE_Conjugate_Gradient(Ncg)
    call Gram_Schmidt
    call PSE_subspace_diag
    call Gram_Schmidt


    call PSE_current_GS(jav)
    call PSE_energy(Etot,Ekin,'GS')
    if(myrank == 0)write(*,'(A,e26.16e3)')'Etot =',Etot
    if(myrank == 0)write(*,'(A,e26.16e3)')'Ekin =',Ekin
    if(myrank == 0)write(*,'(3(A6,e16.6e3))')'jav(1)= ',jav(1),'jav(2)= ',jav(2),'jav(3)= ',jav(3)

!    Vloc=Vxc ! debug
    if(myrank == 0)write(*,*)'charge =',sum(rho_c)*H123,sum(rho_e)*H123,sum(rho_p)*H123
    if(myrank == 0)write(*,*)sum(Vh**2)*H123

    if(Myrank == 0)then
      write(*,'(4(I3,e16.6e3))')(ib,esp(ib,NK_s),ib=1,NB)
    end if
  end do
  zu_basis(:,1:NB_basis_main,NK_s:NK_e) = zu_GS(:,1:NB_basis_main,NK_s:NK_e)

  do ikshift = 1,NK_shift
    if(myrank == 0)write(*,"(A,2x,I7)")"ikshift =",ikshift
    kAc_Cvec(1,:)=kAc0_Cvec(1,:)+kshift(1,ikshift)
    kAc_Cvec(2,:)=kAc0_Cvec(2,:)+kshift(2,ikshift)
    kAc_Cvec(3,:)=kAc0_Cvec(3,:)+kshift(3,ikshift)

    do iter_scf=1,Nscf
      if(myrank == 0)write(*,*)kAc_Cvec(:,1)
      if(Myrank == 0)write(*,'(A,2X,I5)')'iter_scf = ',iter_scf

      call Gram_Schmidt
      call PSE_Conjugate_Gradient(Ncg)
      call Gram_Schmidt
      call PSE_subspace_diag
      call Gram_Schmidt


      call PSE_current_GS(jav)
      call PSE_energy(Etot,Ekin,'GS')
      if(myrank == 0)write(*,'(A,e26.16e3)')'Etot =',Etot
      if(myrank == 0)write(*,'(A,e26.16e3)')'Ekin =',Ekin
      if(myrank == 0)write(*,'(3(A6,e16.6e3))')'jav(1)= ',jav(1),'jav(2)= ',jav(2),'jav(3)= ',jav(3)

!    Vloc=Vxc ! debug
      if(myrank == 0)write(*,*)'charge =',sum(rho_c)*H123,sum(rho_e)*H123,sum(rho_p)*H123
      if(myrank == 0)write(*,*)sum(Vh**2)*H123

      if(Myrank == 0)then
        write(*,'(4(I3,e16.6e3))')(ib,esp(ib,NK_s),ib=1,NB)
      end if
    end do
    NBs = NB_basis_main + NB_basis_shift*(ikshift-1) + 1
    NBe = NB_basis_main + NB_basis_shift*ikshift
    zu_basis(:,NBs:NBe,NK_s:NK_e) = zu_GS(:,NBs:NBe,NK_s:NK_e)
  end do

  do ik=NK_s,NK_e
    do ib1 = 1,NB_basis

      do ib2 = 1,ib1-1
        zs = sum(conjg(zu_basis(:,ib2,ik))*zu_basis(:,ib1,ik))*H123
        zu_basis(:,ib1,ik) = zu_basis(:,ib1,ik) - zs*zu_basis(:,ib2,ik)
      end do

      ss = sum(abs(zu_basis(:,ib1,ik))**2)*H123
      zu_basis(:,ib1,ik)=zu_basis(:,ib1,ik)/sqrt(ss)

    end do
  end do

  return
end subroutine PSE_preparation_basis_set
