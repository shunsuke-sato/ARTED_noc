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
subroutine PSE_preparation_matrix
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: f0_1,f0_2,omega_1,omega_2,tpulse_1,tpulse_2,T1_T2
  integer :: ik,ib1,ib2,iav
  complex(8) :: zs
  integer :: ilma,ia,i,j,ix,iy,iz
  complex(8) :: zjx1t,zjx2t,zjx3t,zjxt,zjyt,zjzt
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz
  real(8) :: x1,x2,x3

  f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
  omega_1=omegaev_1/(2d0*Ry)  ! frequency in a.u.
  Amax = f0_1/omega_1
  dAmax = Amax/dble(NAmax)

  allocate(zH_loc(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zPi_loc(NB_basis,NB_basis,NK_s:NK_e))
  allocate(zV_NL(NB_basis,NB_basis,NK_s:NK_e,-NAmax:NAmax))
  allocate(zPi_NL(NB_basis,NB_basis,NK_s:NK_e,-NAmax:NAmax))

  kAc_Cvec(:,:)=kAc0_Cvec(:,:)
  call PSE_pre_hpsi
! compute zH_loc
  zH_loc = 0d0
  do ik=NK_s,NK_e
    do ib2=1,NB_basis
      tpsi(:)=zu_basis(:,ib2,ik)
      call PSE_hpsi(ik)
      do ib1=ib2,NB_basis
        zs=sum(conjg(zu_basis(:,ib1,ik))*htpsi(:))*H123
        zH_loc(ib1,ib2,ik)=zs
        zH_loc(ib2,ib1,ik)=conjg(zs)
        if(ib1 == ib2) zH_loc(ib1,ib2,ik)=real(zs)
      end do
    end do
  end do

! compute zPi_loc
  zPi_loc = 0d0
  do ik=NK_s,NK_e
    do ib2=1,NB_basis
      tpsi(:)=zu_basis(:,ib2,ik)
      call PSE_ppsi(ik)
      do ib1=ib2,NB_basis
        zs=sum(conjg(zu_basis(:,ib1,ik))*htpsi(:))*H123
        zPi_loc(ib1,ib2,ik)=zs
        zPi_loc(ib2,ib1,ik)=conjg(zs)
        if(ib1 == ib2) zPi_loc(ib1,ib2,ik)=real(zs)
      end do
    end do
  end do

! compute zV_NL
  do iav = -NAmax,NAmax
    kAc_Cvec(1,:)=kAc0_Cvec(1,:) + dAmax*dble(iav)*Epdir_1(1)
    kAc_Cvec(2,:)=kAc0_Cvec(2,:) + dAmax*dble(iav)*Epdir_1(2)
    kAc_Cvec(3,:)=kAc0_Cvec(3,:) + dAmax*dble(iav)*Epdir_1(3)
    call PSE_pre_hpsi
    do ik=NK_s,NK_e
      do ib2=1,NB_basis
        tpsi(:)=zu_basis(:,ib2,ik)
        call PSE_Vnlpsi(ik)
        do ib1=ib2,NB_basis
          zs=sum(conjg(zu_basis(:,ib1,ik))*htpsi(:))*H123
          zV_NL(ib1,ib2,ik,iav)=zs
          zV_NL(ib2,ib1,ik,iav)=conjg(zs)
          if(ib1 == ib2) zV_NL(ib1,ib2,ik,iav)=real(zs)
        end do
      end do
    end do
  end do


! compute zPi_NL
  do iav = -NAmax,NAmax
    kAc_Cvec(1,:)=kAc0_Cvec(1,:) + dAmax*dble(iav)*Epdir_1(1)
    kAc_Cvec(2,:)=kAc0_Cvec(2,:) + dAmax*dble(iav)*Epdir_1(2)
    kAc_Cvec(3,:)=kAc0_Cvec(3,:) + dAmax*dble(iav)*Epdir_1(3)
    call PSE_pre_hpsi

    do ik=NK_s,NK_e
      do ib2=1,NB_basis
        do ib1=ib2,NB_basis
! <b1| r*VNL | b2>
          zjx1t=0d0;zjx2t=0d0;zjx3t=0d0
          do ilma=1,Nlma
            ia=a_tbl(ilma)
            uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
            do j=1,Mps(ia)
              i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
              x1=Lx(1,i)-dble(Jxx(j,ia))
              x2=Lx(2,i)-dble(Jyy(j,ia))
              x3=Lx(3,i)-dble(Jzz(j,ia))

              uVpsi =uVpsi +uV(j,ilma)*ekr(j,ia,ik)  *zu_basis(i,ib2,ik)
              uVpsix=uVpsix+uV(j,ilma)*ekr(j,ia,ik)*x1*zu_basis(i,ib1,ik)
              uVpsiy=uVpsiy+uV(j,ilma)*ekr(j,ia,ik)*x2*zu_basis(i,ib1,ik)
              uVpsiz=uVpsiz+uV(j,ilma)*ekr(j,ia,ik)*x3*zu_basis(i,ib1,ik)
            end do
            uVpsi =uVpsi *H123*iuV(ilma)
            uVpsix=uVpsix*H123
            uVpsiy=uVpsiy*H123
            uVpsiz=uVpsiz*H123
            
            zjx1t=zjx1t+conjg(uVpsix)*uVpsi/zI
            zjx2t=zjx2t+conjg(uVpsiy)*uVpsi/zI
            zjx3t=zjx3t+conjg(uVpsiz)*uVpsi/zI
          enddo

          zjxt=A_matrix(1,1)*zjx1t+A_matrix(1,2)*zjx2t+A_matrix(1,3)*zjx3t
          zjyt=A_matrix(2,1)*zjx1t+A_matrix(2,2)*zjx2t+A_matrix(2,3)*zjx3t
          zjzt=A_matrix(3,1)*zjx1t+A_matrix(3,2)*zjx2t+A_matrix(3,3)*zjx3t

          zs = Epdir_1(1)*zjxt + Epdir_1(2)*zjyt+ Epdir_1(3)*zjzt

! <b1| VNL*r | b2>
          zjx1t=0d0;zjx2t=0d0;zjx3t=0d0
          do ilma=1,Nlma
            ia=a_tbl(ilma)
            uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
            do j=1,Mps(ia)
              i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
              x1=Lx(1,i)-dble(Jxx(j,ia))
              x2=Lx(2,i)-dble(Jyy(j,ia))
              x3=Lx(3,i)-dble(Jzz(j,ia))

              uVpsi =uVpsi +uV(j,ilma)*ekr(j,ia,ik)  *zu_basis(i,ib1,ik)
              uVpsix=uVpsix+uV(j,ilma)*ekr(j,ia,ik)*x1*zu_basis(i,ib2,ik)
              uVpsiy=uVpsiy+uV(j,ilma)*ekr(j,ia,ik)*x2*zu_basis(i,ib2,ik)
              uVpsiz=uVpsiz+uV(j,ilma)*ekr(j,ia,ik)*x3*zu_basis(i,ib2,ik)
            end do
            uVpsi =uVpsi *H123*iuV(ilma)
            uVpsix=uVpsix*H123
            uVpsiy=uVpsiy*H123
            uVpsiz=uVpsiz*H123
            
            zjx1t=zjx1t+conjg(uVpsi)*uVpsix/zI
            zjx2t=zjx2t+conjg(uVpsi)*uVpsiy/zI
            zjx3t=zjx3t+conjg(uVpsi)*uVpsiz/zI
          enddo

          zjxt=A_matrix(1,1)*zjx1t+A_matrix(1,2)*zjx2t+A_matrix(1,3)*zjx3t
          zjyt=A_matrix(2,1)*zjx1t+A_matrix(2,2)*zjx2t+A_matrix(2,3)*zjx3t
          zjzt=A_matrix(3,1)*zjx1t+A_matrix(3,2)*zjx2t+A_matrix(3,3)*zjx3t

          zs = zs - (Epdir_1(1)*zjxt + Epdir_1(2)*zjyt+ Epdir_1(3)*zjzt)
          zPi_NL(ib1,ib2,ik,iav)=zs
          zPi_NL(ib2,ib1,ik,iav)=conjg(zs)
          if(ib1 == ib2) zPi_NL(ib1,ib2,ik,iav)=real(zs)

        end do
      end do
    end do
  end do

  return
end subroutine PSE_preparation_matrix
