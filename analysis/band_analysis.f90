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
program band_analysis
  implicit none
  real(8),parameter :: Pi=3.141592653589793d0
  complex(8),parameter :: zI=(0d0,1d0)
  real(8),allocatable :: esp(:,:),occ(:,:),esp_NB(:)
  real(8),allocatable :: esp_k(:,:,:,:),occ_k(:,:,:,:)
  complex(8),allocatable :: esp_x(:,:,:,:),occ_x(:,:,:,:)
  complex(8),allocatable :: exp_ikx(:,:),exp_iky(:,:),exp_ikz(:,:)
  complex(8),allocatable :: cexp_ikx(:,:),cexp_iky(:,:),cexp_ikz(:,:)
  complex(8),allocatable :: zft1(:,:,:),zft2(:,:,:)
  integer :: NK,NB,NK1,NK2,NK3,Np
  integer :: ik,ib,ik1,ik2,ik3,i,ik1t,ik2t,ik3t
  real(8) :: Amat(3,3)
  real(8) :: Kxyz,kx,ky,kz,kx0,ky0,kz0,kxf,kyf,kzf,k1,k2,k3,klength
  real(8) :: Ef
  

  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*)
  read(*,*)NK,NB !,NK1,NK2,NK3
  NK1=24
  NK2=24
  NK3=24
  allocate(esp(NB,NK),occ(NB,NK))

  do ik=1,NK
    do ib=1,NB
      read(*,*)esp(ib,ik),occ(ib,ik)
    end do
  end do

  allocate(esp_k(0:NK1-1,0:NK2-1,0:NK3-1,NB),occ_k(0:NK1-1,0:NK2-1,0:NK3-1,NB))
  allocate(esp_x(0:NK1-1,0:NK2-1,0:NK3-1,NB),occ_x(0:NK1-1,0:NK2-1,0:NK3-1,NB))

  Ef=maxval(esp(1:4,:))
  esp=esp-Ef
  ik=0
  do ik1=0,NK1-1
  do ik2=0,NK2-1
  do ik3=0,NK3-1
    ik=ik+1
    do ib=1,NB
      esp_k(ik1,ik2,ik3,ib)=esp(ib,ik)
      occ_k(ik1,ik2,ik3,ib)=occ(ib,ik)
    end do
  end do
  end do
  end do


!  ik3=NK3/2
!  do ik1=0,NK1-1
!  do ik2=0,NK2-1
!    write(*,'(100e26.16e3)')dble(ik1)/dble(NK1)-(0.5d0+0.5d0/dble(NK1)) &
!      ,dble(ik2)/dble(NK2)-(0.5d0+0.5d0/dble(NK2)) &
!      ,esp_k(ik1,ik2,ik3,1),esp_k(ik1,ik2,ik3,2)
!  do ik3=0,NK3-1
!    write(*,'(100e16.6e3)')dble(ik1)/dble(NK1)-(0.5d0+0.5d0/dble(NK1)) &
!      ,dble(ik2)/dble(NK2)-(0.5d0+0.5d0/dble(NK2)) &
!      ,dble(ik3)/dble(NK3)-(0.5d0+0.5d0/dble(NK3)) &
!      ,esp_k(ik1,ik2,ik3,1)
 ! end do
!  end do
!  end do
!  stop

  allocate(esp_NB(NB))
  allocate(exp_ikx(0:NK1-1,0:NK1-1),exp_iky(0:NK2-1,0:NK2-1),exp_ikz(0:NK3-1,0:NK3-1))
  allocate(cexp_ikx(0:NK1-1,0:NK1-1),cexp_iky(0:NK2-1,0:NK2-1),cexp_ikz(0:NK3-1,0:NK3-1))
  allocate(zft1(0:NK1-1,0:NK2-1,0:NK3-1),zft2(0:NK1-1,0:NK2-1,0:NK3-1))

  do ik1=0,NK1-1
    do ik=0,NK1-1
      exp_ikx(ik,ik1)=exp(zI*2d0*pi*dble(ik*ik1)/dble(NK1))
      cexp_ikx(ik,ik1)=exp(-zI*2d0*pi*dble(ik*ik1)/dble(NK1))
    end do
  end do

  do ik2=0,NK2-1
    do ik=0,NK2-1
      exp_iky(ik,ik2)=exp(zI*2d0*pi*dble(ik*ik2)/dble(NK2))
      cexp_iky(ik,ik2)=exp(-zI*2d0*pi*dble(ik*ik2)/dble(NK2))
    end do
  end do

  do ik3=0,NK3-1
    do ik=0,NK3-1
      exp_ikz(ik,ik3)=exp(zI*2d0*pi*dble(ik*ik3)/dble(NK3))
      cexp_ikz(ik,ik3)=exp(-zI*2d0*pi*dble(ik*ik3)/dble(NK3))
    end do
  end do



  do ib=1,NB

  do ik1=0,NK1-1
  do ik2=0,NK2-1
  do ik3=0,NK3-1
    zft1(ik1,ik2,ik3)=sum(cexp_ikz(:,ik3)*esp_k(ik1,ik2,:,ib))
  end do
  end do
  end do
  do ik1=0,NK1-1
  do ik2=0,NK2-1
  do ik3=0,NK3-1
    zft2(ik1,ik2,ik3)=sum(cexp_iky(:,ik2)*zft1(ik1,:,ik3))
  end do
  end do
  end do
  do ik1=0,NK1-1
  do ik2=0,NK2-1
  do ik3=0,NK3-1
    zft1(ik1,ik2,ik3)=sum(cexp_ikx(:,ik1)*zft2(:,ik2,ik3))
  end do
  end do
  end do

  esp_x(:,:,:,ib)=zft1(:,:,:)/dble(NK1*NK2*NK3)
  
  end do

!== diamond structure ==!
  Np=30
  Amat(1,1)=0d0 ; Amat(1,2)=1d0 ; Amat(1,3)=1d0
  Amat(2,1)=1d0 ; Amat(2,2)=0d0 ; Amat(2,3)=1d0
  Amat(3,1)=1d0 ; Amat(3,2)=1d0 ; Amat(3,3)=0d0
  Amat=Amat/2d0
! from L-point(0.5,0.5,0.5) to \Gamma point(0.0,0.0,0.0)
  kx0=0.5d0 ; ky0=0.5d0 ; kz0=0.5d0
  kxf=0.0d0 ; kyf=0.0d0 ; kzf=0.0d0
  klength=sqrt((kxf-kx0)**2+(kyf-ky0)**2+(kzf-kz0)**2)
  kxyz=-klength/dble(Np)
  write(*,'(A)')'# L-point'
  do i=0,Np
    kxyz=kxyz+klength/dble(Np)
    kx=kx0+(kxf-kx0)*dble(i)/dble(Np)
    ky=ky0+(kyf-ky0)*dble(i)/dble(Np)
    kz=kz0+(kzf-kz0)*dble(i)/dble(Np)

    k1=Amat(1,1)*kx+Amat(1,2)*ky+Amat(1,3)*kz+(0.5d0-0.5d0/dble(NK1))
    k2=Amat(2,1)*kx+Amat(2,2)*ky+Amat(2,3)*kz+(0.5d0-0.5d0/dble(NK2))
    k3=Amat(3,1)*kx+Amat(3,2)*ky+Amat(3,3)*kz+(0.5d0-0.5d0/dble(NK3))
!    k1=0.52d0
!    k2=0.52d0
!    k3=0.52d0
!    write(*,'(A,100e16.6e3)')'#',k1,k2,k3,kx,ky,kz
    esp_NB(:)=0d0
    do ik1=0,NK1-1
      ik1t=ik1
      if(ik1 >= NK1/2)ik1t=ik1-NK1
      do ik2=0,NK2-1
        ik2t=ik2
        if(ik2 >= NK2/2)ik2t=ik2-NK2
        do ik3=0,NK3-1
          ik3t=ik3
          if(ik3 >= NK3/2)ik3t=ik3-NK3
          esp_NB(:)=esp_NB(:)+exp(zI*2d0*pi*(k1*dble(ik1t)+k2*dble(ik2t)+k3*dble(ik3t)))*esp_x(ik1,ik2,ik3,:)
        end do
      end do
    end do
      

    write(*,'(999e26.16e3)')kxyz,esp_NB(:)
  end do

! from \Gamma-point(0.0,0.0,0.0) to X point(1.0,0.0,0.0)
  kx0=0.0d0 ; ky0=0.0d0 ; kz0=0.0d0
  kxf=1.0d0 ; kyf=0.0d0 ; kzf=0.0d0
  klength=sqrt((kxf-kx0)**2+(kyf-ky0)**2+(kzf-kz0)**2)
  write(*,'(A)')'# \Gamma-point'
  do i=1,Np
    kxyz=kxyz+klength/dble(Np)
    kx=kx0+(kxf-kx0)*dble(i)/dble(Np)
    ky=ky0+(kyf-ky0)*dble(i)/dble(Np)
    kz=kz0+(kzf-kz0)*dble(i)/dble(Np)

    k1=Amat(1,1)*kx+Amat(1,2)*ky+Amat(1,3)*kz+(0.5d0-0.5d0/dble(NK1))
    k2=Amat(2,1)*kx+Amat(2,2)*ky+Amat(2,3)*kz+(0.5d0-0.5d0/dble(NK2))
    k3=Amat(3,1)*kx+Amat(3,2)*ky+Amat(3,3)*kz+(0.5d0-0.5d0/dble(NK3))
!    k1=0.52d0
!    k2=0.52d0
!    k3=0.52d0
!    write(*,'(A,100e16.6e3)')'#',k1,k2,k3,kx,ky,kz
    esp_NB(:)=0d0
    do ik1=0,NK1-1
      ik1t=ik1
      if(ik1 >= NK1/2)ik1t=ik1-NK1
      do ik2=0,NK2-1
        ik2t=ik2
        if(ik2 >= NK2/2)ik2t=ik2-NK2
        do ik3=0,NK3-1
          ik3t=ik3
          if(ik3 >= NK3/2)ik3t=ik3-NK3
          esp_NB(:)=esp_NB(:)+exp(zI*2d0*pi*(k1*dble(ik1t)+k2*dble(ik2t)+k3*dble(ik3t)))*esp_x(ik1,ik2,ik3,:)
        end do
      end do
    end do
      

    write(*,'(999e26.16e3)')kxyz,esp_NB(:)
  end do

! from X-point(1.0,0.0,0.0) to U point(1,1/4,1/4)
  kx0=1.0d0 ; ky0=0.0d0 ; kz0=0.0d0
  kxf=1d0 ; kyf=1d0/4d0 ; kzf=1d0/4d0
  klength=sqrt((kxf-kx0)**2+(kyf-ky0)**2+(kzf-kz0)**2)
  write(*,'(A)')'# X-point'
  do i=1,Np
    kxyz=kxyz+klength/dble(Np)
    kx=kx0+(kxf-kx0)*dble(i)/dble(Np)
    ky=ky0+(kyf-ky0)*dble(i)/dble(Np)
    kz=kz0+(kzf-kz0)*dble(i)/dble(Np)

    k1=Amat(1,1)*kx+Amat(1,2)*ky+Amat(1,3)*kz+(0.5d0-0.5d0/dble(NK1))
    k2=Amat(2,1)*kx+Amat(2,2)*ky+Amat(2,3)*kz+(0.5d0-0.5d0/dble(NK2))
    k3=Amat(3,1)*kx+Amat(3,2)*ky+Amat(3,3)*kz+(0.5d0-0.5d0/dble(NK3))
!    k1=0.52d0
!    k2=0.52d0
!    k3=0.52d0
!    write(*,'(A,100e16.6e3)')'#',k1,k2,k3,kx,ky,kz
    esp_NB(:)=0d0
    do ik1=0,NK1-1
      ik1t=ik1
      if(ik1 >= NK1/2)ik1t=ik1-NK1
      do ik2=0,NK2-1
        ik2t=ik2
        if(ik2 >= NK2/2)ik2t=ik2-NK2
        do ik3=0,NK3-1
          ik3t=ik3
          if(ik3 >= NK3/2)ik3t=ik3-NK3
          esp_NB(:)=esp_NB(:)+exp(zI*2d0*pi*(k1*dble(ik1t)+k2*dble(ik2t)+k3*dble(ik3t)))*esp_x(ik1,ik2,ik3,:)
        end do
      end do
    end do
      

    write(*,'(999e26.16e3)')kxyz,esp_NB(:)
  end do

! from K-point(3/4,3/4,0) to \Gamma-point(0.0,0.0,0.0)
  kx0=3d0/4d0 ; ky0=3d0/4d0 ; kz0=0d0
  kxf=0.0d0 ; kyf=0.0d0 ; kzf=0.0d0

  klength=sqrt((kxf-kx0)**2+(kyf-ky0)**2+(kzf-kz0)**2)
  write(*,'(A)')'# U,K-point'
  do i=1,Np
    kxyz=kxyz+klength/dble(Np)
    kx=kx0+(kxf-kx0)*dble(i)/dble(Np)
    ky=ky0+(kyf-ky0)*dble(i)/dble(Np)
    kz=kz0+(kzf-kz0)*dble(i)/dble(Np)

    k1=Amat(1,1)*kx+Amat(1,2)*ky+Amat(1,3)*kz+(0.5d0-0.5d0/dble(NK1))
    k2=Amat(2,1)*kx+Amat(2,2)*ky+Amat(2,3)*kz+(0.5d0-0.5d0/dble(NK2))
    k3=Amat(3,1)*kx+Amat(3,2)*ky+Amat(3,3)*kz+(0.5d0-0.5d0/dble(NK3))
!    k1=0.52d0
!    k2=0.52d0
!    k3=0.52d0
!    write(*,'(A,100e16.6e3)')'#',k1,k2,k3,kx,ky,kz
    esp_NB(:)=0d0
    do ik1=0,NK1-1
      ik1t=ik1
      if(ik1 >= NK1/2)ik1t=ik1-NK1
      do ik2=0,NK2-1
        ik2t=ik2
        if(ik2 >= NK2/2)ik2t=ik2-NK2
        do ik3=0,NK3-1
          ik3t=ik3
          if(ik3 >= NK3/2)ik3t=ik3-NK3
          esp_NB(:)=esp_NB(:)+exp(zI*2d0*pi*(k1*dble(ik1t)+k2*dble(ik2t)+k3*dble(ik3t)))*esp_x(ik1,ik2,ik3,:)
        end do
      end do
    end do
      

    write(*,'(999e26.16e3)')kxyz,esp_NB(:)
  end do
  write(*,'(A)')'# \Gamma-point'
!== diamond structure ==!

end program band_analysis
  
