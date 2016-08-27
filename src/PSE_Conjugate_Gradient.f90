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
subroutine PSE_Conjugate_Gradient(iter_cg_max)
  use global_variables
  implicit none
  real(8),parameter :: delta_cg=1d-15
  integer :: iter_cg,iter_cg_max
  real(8) :: xkHxk,gkgk,pkHpk
  real(8) :: s,uk,ev
  complex(8) :: cx,cp,xkHpk
  complex(8) :: zs
  integer :: ik,ib,ibt
  real(8) :: esp_var_l(NB,NK)

  esp_var_l=0d0

  call PSE_pre_hpsi

  K_points :do ik=NK_s,NK_e
  Band :    do ib=1,NB


    do ibt=1,ib-1
      zs=sum(conjg(zu_GS(:,ibt,ik))*zu_GS(:,ib,ik))*H123
      zu_GS(1:NL,ib,ik)=zu_GS(1:NL,ib,ik)-zu_GS(1:NL,ibt,ik)*zs
    end do
    s=sum(abs(zu_GS(:,ib,ik))**2)*H123
    xk(1:NL)=zu_GS(1:NL,ib,ik)/sqrt(s)
    tpsi(:)=xk(:)
    call PSE_hpsi(ik)
    hxk(:)=htpsi(:)
    xkHxk=sum(conjg(xk(:))*hxk(:))*H123

    do iter_cg=1,iter_cg_max
      gk(:)=hxk(:)-xkHxk*xk(:)
      do ibt=1,ib-1
        zs=sum(conjg(zu_GS(:,ibt,ik))*gk(:))*H123
        gk(1:NL)=gk(1:NL)-zu_GS(1:NL,ibt,ik)*zs
      end do
      s=sum(abs(gk(:))**2)*H123
      select case (iter_cg)
      case(1)
        pk(1:NL)=gk(1:NL)
      case default
        uk=s/gkgk
        pk(1:NL)=gk(1:NL)+uk*pk(1:NL)
      end select
      gkgk=s
      zs=sum(conjg(xk(:))*pk(:))*H123
      pko(1:NL)=pk(1:NL)-xk(1:NL)*zs
      s=sum(abs(pko(:))**2)*H123
      pko(1:NL)=pko(1:NL)/sqrt(s)
      tpsi(1:NL)=pko(1:NL); call PSE_hpsi(ik)
      xkHpk=sum(conjg(xk(:))*htpsi(:))*H123
      pkHpk=sum(conjg(pko(:))*htpsi(:))*H123
      ev=0.5d0*((xkHxk+pkHpk)-sqrt((xkHxk-pkHpk)**2+4*abs(xkHpk)**2))
      cx=xkHpk/(ev-xkHxk)
      cp=1.d0/sqrt(1.d0+abs(cx)**2)
      cx=cx*cp
      if(abs(ev-xkHxk)<delta_cg) exit
      xk(1:NL)=cx*xk(1:NL)+cp*pko(1:NL)
      hxk(1:NL)=cx*hxk(1:NL)+cp*htpsi(1:NL)
      xkHxk=sum(conjg(xk(:))*hxk(:))*H123
    end do
      s=sum(abs(xk(:))**2)*H123
      zu_GS(1:NL,ib,ik)=xk(1:NL)/sqrt(s)
      tpsi(1:NL)=zu_GS(1:NL,ib,ik); call PSE_hpsi(ik)
      xkHxk=sum(conjg(tpsi(:))*htpsi(:))*H123
      esp_var_l(ib,ik)=sqrt(sum(abs(htpsi(:)-xkHxk*tpsi(:))**2)*H123)*occ(ib,ik)
  end do Band
  end do K_points

  call MPI_ALLREDUCE(esp_var_l,esp_var,NB*NK,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  return
end subroutine PSE_Conjugate_Gradient
  
