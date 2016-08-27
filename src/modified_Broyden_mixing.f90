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
subroutine modified_Broyden_mixing(iter_scf)
! Refference
! 1. D.D. Johnson, Phys. Rev. B 38, 12807 (1988)
! 2. Rok Zitko, Phys. Rev. B 80, 125125 (2009)
  use global_variables
  implicit none
  real(8),parameter :: omega0_MB = 0.01d0
  integer :: iter_scf
  integer :: Nmax
  real(8),allocatable :: delta_rho_MB(:,:),F_MB(:,:),delta_F_MB(:,:)
  real(8),allocatable :: Amat(:,:),beta_MB(:,:),ck_MB(:)
  integer :: k,n,i
  real(8) :: ss
! == LAPACK ==
  integer :: lwork,info
  integer,allocatable :: ipiv(:)
  real(8),allocatable :: work(:)
! == LAPACK ==

  if(iter_scf <= MaxMem_MB )then
    rho_e(:) = rho_MB_in(:,iter_scf) + alpha_MB*(rho_MB_out(:,iter_scf) - rho_MB_in(:,iter_scf))
  else


  Nmax = min(iter_scf,MaxMem_MB)-1
  allocate(delta_rho_MB(NL,Nmax+1),F_MB(NL,Nmax+1),delta_F_MB(NL,Nmax+1))
  allocate(Amat(Nmax,Nmax),beta_MB(Nmax,Nmax),ck_MB(Nmax))

! == LAPACK ==
  lwork = 64*Nmax
  allocate(work(lwork),ipiv(Nmax))
! == LAPACK ==

  do n = 1,Nmax+1
    F_MB(:,n) = rho_MB_out(:,n) - rho_MB_in(:,n)
  end do

  do n = 1,Nmax
    ss = sum( (F_MB(:,n+1) - F_MB(:,n))**2 ); ss = 1d0/sqrt(ss)
    delta_rho_MB(:,n) = (rho_MB_in(:,n+1)-rho_MB_in(:,n))*ss
    delta_F_MB(:,n) = (F_MB(:,n+1) - F_MB(:,n))*ss
  end do



  do k = 1,Nmax
    do n = 1,Nmax
      Amat(k,n) = 1d0*1d0*sum( delta_F_MB(:,n)*delta_F_MB(:,k) )
      if(k == n) Amat(k,n) = Amat(k,n) + omega0_MB**2
    end do
  end do

! == LAPACK ==
  call dgetrf(Nmax,Nmax,Amat,Nmax,ipiv,info)
  call dgetri(Nmax,Amat,Nmax,ipiv,work,lwork,info)
! == LAPACK ==
  do k = 1,Nmax
    do n = 1,Nmax
      beta_MB(k,n) = Amat(k,n)
    end do
  end do

  do k = 1,Nmax
    ck_MB(k) = sum(delta_F_MB(:,k)*F_MB(:,Nmax+1))
  end do

  rho_e(:) = rho_MB_in(:,Nmax+1) + alpha_MB*F_MB(:,Nmax+1)

  do n = 1,Nmax
    ss = sum(ck_MB(:)*beta_MB(:,n))
    rho_e(:) = rho_e(:) - ss*(alpha_MB*delta_F_MB(:,n) + delta_rho_MB(:,n))
  end do
  
  do i = 1,NL
    if(rho_e(i) <= 0d0)then
      if(myrank == 0)write(*,'(A)')'Negative density in modified Broyden'
      rho_e(:) = rho_MB_in(:,Nmax+1) + alpha_MB*F_MB(:,Nmax+1)
      exit
    end if
  end do

  end if

  if(iter_scf >= MaxMem_MB)then
    do k = 2,MaxMem_MB
      rho_MB_in(:,k-1) = rho_MB_in(:,k)
      rho_MB_out(:,k-1) = rho_MB_out(:,k)
    end do
  end if

  return
end subroutine modified_Broyden_mixing
  
