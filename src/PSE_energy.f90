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
subroutine PSE_energy(Etot,Ekin,GS_RT)
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: Etot,Ekin,E_h,Eexc,Eloc,Enl
  real(8) :: Etot_l,Ekin_l,Enl_l
  character(2) :: GS_RT
  integer :: ik,ib,NBt
  integer :: ix,iy,iz
  integer :: ilma,ia,i,j
  complex(8) :: uVpsi

  if(GS_RT == 'GS') then
    NBt=NB
  else if(GS_RT == 'RT') then
    NBt=NB_TD
  end if

  Ekin_l=0d0
!== Start Ekin
  K_points2 :do ik=NK_s,NK_e

    zft2(:,:,:)=-0.5d0*Lap_k(:,:,:) &
      -kAc_Cvec(1,ik)*Grad_x_zI(:,:,:) &
      -kAc_Cvec(2,ik)*Grad_y_zI(:,:,:) &
      -kAc_Cvec(3,ik)*Grad_z_zI(:,:,:) &
      +sum(kAc_Cvec(:,ik)**2)*0.5
    
    zft3(:,:,:)=zft2(:,:,:)/(dble(NL1*NL2*NL3))

  Band2 :    do ib=1,NBt
    if(GS_RT == 'GS') then
      tpsi(:)=zu_GS(:,ib,ik)
    else if(GS_RT == 'RT') then
      tpsi(:)=zu(:,ib,ik)
    end if


    do i=1,NL
      zft2(iLx(1,i),iLx(2,i),iLx(3,i))=tpsi(i)
    end do

    do ix=0,NL1-1
    do iy=0,NL2-1
    do iz=0,NL3-1
      zft1(ix,iy,iz)=sum(cexp_x3(:,iz)*zft2(ix,iy,:))
    end do
    end do
    end do

    do ix=0,NL1-1
    do iy=0,NL2-1
    do iz=0,NL3-1
      zft2(ix,iy,iz)=sum(cexp_x2(:,iy)*zft1(ix,:,iz))
    end do
    end do
    end do

    do ix=0,NL1-1
    do iy=0,NL2-1
    do iz=0,NL3-1
      zft1(ix,iy,iz)=sum(cexp_x1(:,ix)*zft2(:,iy,iz))
    end do
    end do
    end do


    Ekin_l=Ekin_l+occ(ib,ik)*sum( zft3(:,:,:)*abs(zft1)**2 )*H123

    end do Band2
    end do K_points2
    call MPI_ALLREDUCE(Ekin_l,Ekin,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!== End Ekin

!== Start Hartree
    E_h=0.5d0*sum(rho_e(:)*Vh(:))*H123
!== End Hartree
!== Start Exc
    Eexc=sum(rho_e(:)*Exc(:))*H123
!== End  Exc
!== Start Exc
    Eloc=sum(rho_e(:)*Vpsl(:))*H123
!== End  Exc


!== Start Nonlocal
    Enl_l=0d0
    call PSE_pre_hpsi
   
    do ik=NK_s,NK_e
      do ib=1,NBt
        if(GS_RT == 'GS') then
          tpsi(:)=zu_GS(:,ib,ik)
        else if(GS_RT == 'RT') then
          tpsi(:)=zu(:,ib,ik)
        end if

        do ilma=1,Nlma
          ia=a_tbl(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            i=Jxyz(j,ia)
            uVpsi=uVpsi+uV(j,ilma)*ekr(j,ia,ik)*tpsi(i)
          enddo
          uVpsi=uVpsi*H123
          Enl_l=Enl_l+occ(ib,ik)*abs(uVpsi)**2*iuV(ilma)
        enddo
      end do
    end do
    call MPI_ALLREDUCE(Enl_l,Enl,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
!== End Nonlocal
    Etot=Ekin+E_h+Eexc+Eloc+Enl


  return
end subroutine PSE_energy
  
