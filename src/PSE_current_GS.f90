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
subroutine PSE_current_GS(jav)
  use global_variables
  use PSE_variables
  implicit none
  real(8) :: jav(3)
  integer :: ik,ib,i,ix,iy,iz,ilma,ia,j
  real(8) :: x1,x2,x3,jxt,jyt,jzt,jx1t,jx2t,jx3t,x,y,z
  complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz
  real(8) :: curr(3),curr_l(3)
  real(8) :: curr_Lvec(3),curr_l_Lvec(3)
  complex(8) :: zt(3)
  real(8) :: bb1,bb2,bb3
  real(8) :: jx1,jx2,jx3

  call PSE_pre_hpsi  
  curr_l=0d0

  jx1=0d0
  jx2=0d0
  jx3=0d0

  jxt=0d0
  jyt=0d0
  jzt=0d0

  K_point : do ik=NK_s,NK_e
  Band : do ib=1,NB
    jxt=jxt+kAc_Cvec(1,ik)*occ(ib,ik)
    jyt=jyt+kAc_Cvec(2,ik)*occ(ib,ik)
    jzt=jzt+kAc_Cvec(3,ik)*occ(ib,ik)

    jx1t=0d0
    jx2t=0d0
    jx3t=0d0
    do i=1,NL
      zt(1)= &!nab1(0)*tpsi(i)&
        &+nab1(1)*(zu_GS(ifdx1(1,i),ib,ik)-zu_GS(ifdx1(-1,i),ib,ik))&
        &+nab1(2)*(zu_GS(ifdx1(2,i),ib,ik)-zu_GS(ifdx1(-2,i),ib,ik))&
        &+nab1(3)*(zu_GS(ifdx1(3,i),ib,ik)-zu_GS(ifdx1(-3,i),ib,ik))&
        &+nab1(4)*(zu_GS(ifdx1(4,i),ib,ik)-zu_GS(ifdx1(-4,i),ib,ik))

      zt(2)= &!nab1(0)*tpsi(i)&
        &+nab2(1)*(zu_GS(ifdx2(1,i),ib,ik)-zu_GS(ifdx2(-1,i),ib,ik))&
        &+nab2(2)*(zu_GS(ifdx2(2,i),ib,ik)-zu_GS(ifdx2(-2,i),ib,ik))&
        &+nab2(3)*(zu_GS(ifdx2(3,i),ib,ik)-zu_GS(ifdx2(-3,i),ib,ik))&
        &+nab2(4)*(zu_GS(ifdx2(4,i),ib,ik)-zu_GS(ifdx2(-4,i),ib,ik))
      
      zt(3)= &!nab1(0)*tpsi(i)&
        &+nab3(1)*(zu_GS(ifdx3(1,i),ib,ik)-zu_GS(ifdx3(-1,i),ib,ik))&
        &+nab3(2)*(zu_GS(ifdx3(2,i),ib,ik)-zu_GS(ifdx3(-2,i),ib,ik))&
        &+nab3(3)*(zu_GS(ifdx3(3,i),ib,ik)-zu_GS(ifdx3(-3,i),ib,ik))&
        &+nab3(4)*(zu_GS(ifdx3(4,i),ib,ik)-zu_GS(ifdx3(-4,i),ib,ik))
      
      jx1t=jx1t+aimag(conjg(zu_GS(i,ib,ik))*zt(1))
      jx2t=jx2t+aimag(conjg(zu_GS(i,ib,ik))*zt(2))
      jx3t=jx3t+aimag(conjg(zu_GS(i,ib,ik))*zt(3))
    end do
    
    jx1=jx1+occ(ib,ik)*jx1t*H123
    jx2=jx2+occ(ib,ik)*jx2t*H123
    jx3=jx3+occ(ib,ik)*jx3t*H123

  end do Band
  end do K_point

  jxt=jxt+B_matrix(1,1)*jx1+B_matrix(2,1)*jx2+B_matrix(3,1)*jx3
  jyt=jyt+B_matrix(1,2)*jx1+B_matrix(2,2)*jx2+B_matrix(3,2)*jx3
  jzt=jzt+B_matrix(1,3)*jx1+B_matrix(2,3)*jx2+B_matrix(3,3)*jx3


  jx1t=0d0;jx2t=0d0;jx3t=0d0

  K_point2 : do ik=NK_s,NK_e
  Band2 : do ib=1,NB_TD
    
    do ilma=1,Nlma
      ia=a_tbl(ilma)
      uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
      do j=1,Mps(ia)
        i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
        x1=Lx(1,i)-dble(Jxx(j,ia))
        x2=Lx(2,i)-dble(Jyy(j,ia))
        x3=Lx(3,i)-dble(Jzz(j,ia))

        uVpsi =uVpsi +uV(j,ilma)*ekr(j,ia,ik)  *zu_GS(i,ib,ik)
        uVpsix=uVpsix+uV(j,ilma)*ekr(j,ia,ik)*x1*zu_GS(i,ib,ik)
        uVpsiy=uVpsiy+uV(j,ilma)*ekr(j,ia,ik)*x2*zu_GS(i,ib,ik)
        uVpsiz=uVpsiz+uV(j,ilma)*ekr(j,ia,ik)*x3*zu_GS(i,ib,ik)
      end do
      uVpsi =uVpsi *H123*iuV(ilma)
      uVpsix=uVpsix*H123
      uVpsiy=uVpsiy*H123
      uVpsiz=uVpsiz*H123
      jx1t=jx1t+occ(ib,ik)*2d0*aimag(conjg(uVpsix)*uVpsi)
      jx2t=jx2t+occ(ib,ik)*2d0*aimag(conjg(uVpsiy)*uVpsi)
      jx3t=jx3t+occ(ib,ik)*2d0*aimag(conjg(uVpsiz)*uVpsi)
    enddo


  end do Band2
  end do K_point2


  jxt=jxt+A_matrix(1,1)*jx1t+A_matrix(1,2)*jx2t+A_matrix(1,3)*jx3t
  jyt=jyt+A_matrix(2,1)*jx1t+A_matrix(2,2)*jx2t+A_matrix(2,3)*jx3t
  jzt=jzt+A_matrix(3,1)*jx1t+A_matrix(3,2)*jx2t+A_matrix(3,3)*jx3t

  curr_l(1)=jxt/Vcell
  curr_l(2)=jyt/Vcell
  curr_l(3)=jzt/Vcell


  call MPI_ALLREDUCE(curr_l,curr,3,MPI_REAL8,MPI_SUM,NEW_COMM_WORLD,ierr)  

  jav(:)=curr(:)

  return
end subroutine PSE_current_GS
