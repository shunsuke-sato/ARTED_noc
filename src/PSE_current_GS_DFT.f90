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
subroutine PSE_current_GS_DFT(jav)
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: jav(3)
  integer :: ik,ib,i,ix,iy,iz,ilma,ia,j
  real(8) :: x1,x2,x3,jxt,jyt,jzt,jx1t,jx2t,jx3t
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz
  real(8) :: curr(3),curr_l(3)


  curr_l=0d0

  K_point : do ik=NK_s,NK_e
  Band : do ib=1,NB
    tpsi(:)=zu_GS(:,ib,ik)

    do i=1,NL
      zft3(iLx(1,i),iLx(2,i),iLx(3,i))=tpsi(i)
    end do
    
    do ix=0,NL1-1
      do iy=0,NL2-1
        do iz=0,NL3-1
          zft1(ix,iy,iz)=sum(cexp_x3(:,iz)*zft3(ix,iy,:))
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

    zft1(:,:,:)=zft1(:,:,:)/dble(NL1*NL2*NL3)


   curr_l(1)=curr_l(1)+occ(ib,ik)*sum(abs(zft1(:,:,:))**2*(-Grad_x_zI(:,:,:)+kAc_Cvec(1,ik)))
   curr_l(2)=curr_l(2)+occ(ib,ik)*sum(abs(zft1(:,:,:))**2*(-Grad_y_zI(:,:,:)+kAc_Cvec(2,ik)))
   curr_l(3)=curr_l(3)+occ(ib,ik)*sum(abs(zft1(:,:,:))**2*(-Grad_z_zI(:,:,:)+kAc_Cvec(3,ik)))
    
    jx1t=0d0;jx2t=0d0;jx3t=0d0
    
    do ilma=1,Nlma
      ia=a_tbl(ilma)
      uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
      do j=1,Mps(ia)
        i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
        x1=Lx(1,i)-Jxx(j,ia)
        x2=Lx(2,i)-Jyy(j,ia)
        x3=Lx(3,i)-Jzz(j,ia)

        uVpsi =uVpsi +uV(j,ilma)*ekr(j,ia,ik)  *zu_GS(i,ib,ik)
        uVpsix=uVpsix+uV(j,ilma)*ekr(j,ia,ik)*x1*zu_GS(i,ib,ik)
        uVpsiy=uVpsiy+uV(j,ilma)*ekr(j,ia,ik)*x2*zu_GS(i,ib,ik)
        uVpsiz=uVpsiz+uV(j,ilma)*ekr(j,ia,ik)*x3*zu_GS(i,ib,ik)
      end do
      uVpsi =uVpsi *H123*iuV(ilma)
      uVpsix=uVpsix*H123
      uVpsiy=uVpsiy*H123
      uVpsiz=uVpsiz*H123
      jx1t=jx1t+occ(ib,ik)*2d0*imag(conjg(uVpsix)*uVpsi)
      jx2t=jx2t+occ(ib,ik)*2d0*imag(conjg(uVpsiy)*uVpsi)
      jx3t=jx3t+occ(ib,ik)*2d0*imag(conjg(uVpsiz)*uVpsi)
    enddo

    jxt=A_matrix(1,1)*jx1t+A_matrix(1,2)*jx2t+A_matrix(1,3)*jx3t
    jyt=A_matrix(2,1)*jx1t+A_matrix(2,2)*jx2t+A_matrix(2,3)*jx3t
    jzt=A_matrix(3,1)*jx1t+A_matrix(3,2)*jx2t+A_matrix(3,3)*jx3t

    curr_l(1)=curr_l(1)+jxt/Vcell
    curr_l(2)=curr_l(2)+jyt/Vcell
    curr_l(3)=curr_l(3)+jzt/Vcell


  end do Band
  end do K_point


  call MPI_ALLREDUCE(curr_l,curr,3,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)  

  jav(:)=curr(:)

  return
end subroutine PSE_current_GS_DFT
