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
!This file conatain one soubroutine.
!SUBROUTINE prep_ps_periodic(property)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine PSE_prep_ps_periodic
  use Global_Variables
  use PSE_Variables
  implicit none
  real(8) :: k1(0:NL1-1),k2(0:NL2-1),k3(0:NL3-1)
  real(8) :: kx,ky,kz
  
  character(11) :: property
  integer :: ik,n,i,a,j,ix,iy,iz,lma,l,m,lm,ir,intr,k
  integer :: lma_tbl((Lmax+1)**2,NI)
  real(8) :: G2sq,s,Vpsl_l(NL),G2,Gd,Gr,x,y,z,r,ratio1,ratio2,r2
  real(8) :: Ylm,dYlm,uVr(0:Lmax),duVr(0:Lmax)
  complex(8),allocatable :: Vion_G(:,:,:),dVloc_G(:,:,:,:)
  real(8) :: dist_Lvec(3)  
  real(8) :: detA
  integer :: ilma,ia

  detA=A_matrix(1,1)*A_matrix(2,2)*A_matrix(3,3)+A_matrix(2,1)*A_matrix(3,2)*A_matrix(1,3)+A_matrix(3,1)*A_matrix(1,2)*A_matrix(2,3) &
    -A_matrix(1,3)*A_matrix(2,2)*A_matrix(3,1)-A_matrix(2,3)*A_matrix(3,2)*A_matrix(1,1)-A_matrix(3,3)*A_matrix(1,2)*A_matrix(2,1)
  if(myrank == 0)then
    write(*,*)'detA',detA
    write(*,*)'Vcell',Vcell
    write(*,*)'Vcell/detA',Vcell/detA
  end if


  allocate(Vion_G(0:NL1-1,0:NL2-1,0:NL3-1),dVloc_G(0:NL1-1,0:NL2-1,0:NL3-1,NE))

  do i=0,NL1-1
    k1(i)=2d0*pi*dble(i)
    if(i >= NL1/2) k1(i)=2d0*pi*dble(i-NL1)
  end do
  do i=0,NL2-1
    k2(i)=2d0*pi*dble(i)
    if(i >= NL2/2) k2(i)=2d0*pi*dble(i-NL2)
  end do
  do i=0,NL3-1
    k3(i)=2d0*pi*dble(i)
    if(i >= NL3/2) k3(i)=2d0*pi*dble(i-NL3)
  end do

! local potential
  do ik=1,NE
    do i=0,NL1-1
    do j=0,NL2-1
    do k=0,NL3-1
      kx=k1(i)*B_matrix(1,1)+k2(j)*B_matrix(2,1)+k3(k)*B_matrix(3,1)
      ky=k1(i)*B_matrix(1,2)+k2(j)*B_matrix(2,2)+k3(k)*B_matrix(3,2)
      kz=k1(i)*B_matrix(1,3)+k2(j)*B_matrix(2,3)+k3(k)*B_matrix(3,3)
      G2sq=sqrt(kx**2+ky**2+kz**2)
      s=0.d0
      if (i+j+k == 0) then
        do ir=2,NRloc(ik)
          r=rad(ir,ik)
          s=s+4*Pi*r**2*(vloctbl(ir,ik)+Zps(ik)/r)*(rad(ir,ik)-rad(ir-1,ik))
        end do
      else
        do ir=2,NRloc(ik)
          r=rad(ir,ik)
          s=s+4*Pi*r**2*sin(G2sq*r)/(G2sq*r)*(vloctbl(ir,ik)+Zps(ik)/r)*(rad(ir,ik)-rad(ir-1,ik))
        end do
      endif
      dVloc_G(i,j,k,ik)=s
    end do
    end do
    end do
  end do
  Vion_G=0.d0
  do a=1,NI
    ik=Kion(a)
    do i=0,NL1-1
    do j=0,NL2-1
    do k=0,NL3-1
      kx=k1(i)*B_matrix(1,1)+k2(j)*B_matrix(2,1)+k3(k)*B_matrix(3,1)
      ky=k1(i)*B_matrix(1,2)+k2(j)*B_matrix(2,2)+k3(k)*B_matrix(3,2)
      kz=k1(i)*B_matrix(1,3)+k2(j)*B_matrix(2,3)+k3(k)*B_matrix(3,3)
      G2=kx**2+ky**2+kz**2
      Gd=k1(i)*Rion_Lvec(1,a)+k2(j)*Rion_Lvec(2,a)+k3(k)*Rion_Lvec(3,a)
      Vion_G(i,j,k)=Vion_G(i,j,k)+dVloc_G(i,j,k,ik)*exp(-zI*Gd)
      if(i+j+k == 0) cycle
      Vion_G(i,j,k)=Vion_G(i,j,k)-4*Pi/G2*Zps(ik)*exp(-zI*Gd)
!      Vion_G(i,j,k)=Vion_G(i,j,k)+4*Pi*InLap_k(i,j,k)*Zps(ik)*exp(-zI*Gd)
    end do
    end do
    end do
  end do

  Vion_G=Vion_G/abs(detA)

  Vpsl=0.d0
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x3(:,iz)*Vion_G(ix,iy,:))
  end do
  end do
  end do
  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft1(ix,iy,iz)=sum(exp_x2(:,iy)*zft2(ix,:,iz))
  end do
  end do
  end do

  do ix=0,NL1-1
  do iy=0,NL2-1
  do iz=0,NL3-1
    zft2(ix,iy,iz)=sum(exp_x1(:,ix)*zft1(:,iy,iz))
  end do
  end do
  end do

  do i=1,NL
    Vpsl(i)=zft2(iLx(1,i),iLx(2,i),iLx(3,i))
  end do


! nonlocal potential
  if (Myrank == 0 ) then
    write(*,*) ''
    write(*,*) '============nonlocal grid data=============='
  endif
  do a=1,NI
    ik=Kion(a)
    j=0
    do ix=-2,2
    do iy=-2,2
    do iz=-2,2
      do i=1,NL
        dist_Lvec(1)=Lx(1,i)-(Rion_Lvec(1,a)+dble(ix))
        dist_Lvec(2)=Lx(2,i)-(Rion_Lvec(2,a)+dble(iy))
        dist_Lvec(3)=Lx(3,i)-(Rion_Lvec(3,a)+dble(iz))
              
        r2=0d0
        do n=1,3
          r2=r2+dist_Lvec(n)*sum(dist_Lvec(:)*mat_vv_a_Cvec(n,:))
        end do
        r=sqrt(r2)

        if (r<Rps(ik)) then
          j=j+1
        endif
      enddo
    enddo
    enddo
    enddo
    Mps(a)=j
    if (Myrank == 0) then
      write(*,*) 'a =',a,'Mps(a) =',Mps(a)
    endif
  end do
  Nps=maxval(Mps(:))
  allocate(Jxyz(Nps,NI),Jxx(Nps,NI),Jyy(Nps,NI),Jzz(Nps,NI))
  allocate(ekr(Nps,NI,NK_s:NK_e)) ! sato

  do a=1,NI
    ik=Kion(a)
    j=0
    do ix=-4,4
    do iy=-4,4
    do iz=-4,4
      do i=1,NL
        dist_Lvec(1)=Lx(1,i)-(Rion_Lvec(1,a)+dble(ix))
        dist_Lvec(2)=Lx(2,i)-(Rion_Lvec(2,a)+dble(iy))
        dist_Lvec(3)=Lx(3,i)-(Rion_Lvec(3,a)+dble(iz))
              
        r2=0d0
        do n=1,3
          r2=r2+dist_Lvec(n)*sum(dist_Lvec(:)*mat_vv_a_Cvec(n,:))
        end do
        r=sqrt(r2)

        if (r<Rps(ik)) then
          j=j+1
          if (j<=Nps) then
            Jxyz(j,a)=i
            Jxx(j,a)=ix
            Jyy(j,a)=iy
            Jzz(j,a)=iz
          endif
        endif
      enddo
    enddo
    enddo
    enddo
  end do

  lma=0
  do a=1,NI
    ik=Kion(a)
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  Nlma=lma

  allocate(a_tbl(Nlma),uV(Nps,Nlma),iuV(Nlma),duV(Nps,Nlma,3))

  lma=0
  do a=1,NI
    ik=Kion(a)
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        a_tbl(lma)=a
        lma_tbl(lm,a)=lma
      enddo
    enddo
  enddo

  do a=1,NI
    ik=Kion(a)
    do j=1,Mps(a)
      i=Jxyz(j,a)
      dist_Lvec(1)=Lx(1,i)-(Rion_Lvec(1,a)+dble(Jxx(j,a)))
      dist_Lvec(2)=Lx(2,i)-(Rion_Lvec(2,a)+dble(Jyy(j,a)))
      dist_Lvec(3)=Lx(3,i)-(Rion_Lvec(3,a)+dble(Jzz(j,a)))
      x=sum(A_matrix(1,:)*dist_Lvec(:))
      y=sum(A_matrix(2,:)*dist_Lvec(:))
      z=sum(A_matrix(3,:)*dist_Lvec(:))
      r=sqrt(x**2+y**2+z**2)+1d-50

      do ir=1,NRps(ik)
        if(radnl(ir,ik).gt.r) exit
      enddo
      intr=ir-1
      if (intr.lt.0.or.intr.ge.NRps(ik))stop 'bad intr at prep_ps'
      ratio1=(r-radnl(intr,ik))/(radnl(intr+1,ik)-radnl(intr,ik))
      ratio2=1-ratio1
      do l=0,Mlps(ik)
        uVr(l)=ratio1*udVtbl(intr+1,l,ik)+ratio2*udVtbl(intr,l,ik)
!        duVr(l)=ratio1*dudVtbl(intr+1,l,ik)+ratio2*dudVtbl(intr,l,ik)
      enddo
      lm=0
      do l=0,Mlps(ik)
        if(inorm(l,ik)==0) cycle
        do m=-l,l
          lm=lm+1
          uV(j,lma_tbl(lm,a))=uVr(l)*Ylm(x,y,z,l,m)
!          duV(j,lma_tbl(lm,a),1)=duVr(l)*(x/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,1)
!          duV(j,lma_tbl(lm,a),2)=duVr(l)*(y/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,2)
!          duV(j,lma_tbl(lm,a),3)=duVr(l)*(z/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,3)
        enddo
      enddo
    enddo


    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        iuV(lma_tbl(lm,a))=inorm(l,ik)
      enddo
    enddo
  enddo


  if(myrank == 0)write(*,*)'PSE_prep_ps_periodic is complete'
  return
End Subroutine PSE_prep_ps_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!This file is "Ylm_dYlm.f"
!This file contain two functions.
!Function Ylm
!Function dYlm
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
! This function Ylm is related to the real spherical harmonics Ylm0 by
!c Ylm=sqrt(4*pi/(2l+1))*r**l*Ylm0 and is a monomial of x,y,z
  Function Ylm(x,y,z,il,im)
    implicit none
    real(8),intent(IN) :: x,y,z
    integer,intent(IN) :: il,im
    integer :: ilm
    real(8) :: Ylm,r2
    
    ilm=il*il+il+1+im
    r2=x*x+y*y+z*z
    select case( ilm )
    case(1)  ; Ylm=1.d0
    case(2)  ; Ylm=-y
    case(3)  ; Ylm=z
    case(4)  ; Ylm=-x
    case(5)  ; Ylm=sqrt(3.d0)*x*y                        ! lm=5  (2 -2)
    case(6)  ; Ylm=-sqrt(3.d0)*y*z                       ! lm=6  (2 -1)
    case(7)  ; Ylm=(3*z*z-r2)/2.d0                       ! lm=7  (2 0)
    case(8)  ; Ylm=-sqrt(3.d0)*x*z                       ! lm=8  (2 1)
    case(9)  ; Ylm=sqrt(3.d0/4.d0)*(x*x-y*y)             ! lm=9  (2 2)
    case(10) ; Ylm=-sqrt(5.d0/8.d0)*y*(3*x*x-y*y)        ! lm=10 (3 -3)
    case(11) ; Ylm=sqrt(15.d0)*x*y*z                     ! lm=11 (3 -2)
    case(12) ; Ylm=-sqrt(3.d0/8.d0)*y*(5*z*z-r2)         ! lm=12 (3 -1)
    case(13) ; Ylm=z*(5*z*z-3*r2)/2.d0                   ! lm=13 (3 0)
    case(14) ; Ylm=-sqrt(3.d0/8.d0)*x*(5*z*z-r2)         ! lm=14 (3 1)
    case(15) ; Ylm=sqrt(15.d0/4.d0)*z*(x*x-y*y)          ! lm=15 (3 2)
    case(16) ; Ylm=-sqrt(5.d0/8.d0)*x*(x*x-3*y*y)        ! lm=16 (3 3)
    case(17) ; Ylm=sqrt(35.d0)/2.d0*x*y*(x*x-y*y)        ! lm=17 (4 -4)
    case(18) ; Ylm=-sqrt(35.d0/8.d0)*y*z*(3*x*x-y*y)     ! lm=18 (4 -3)
    case(19) ; Ylm=sqrt(5.d0)/2.d0*x*y*(7*z*z-r2)        ! lm=19 (4 -2)
    case(20) ; Ylm=-sqrt(5.d0/8.d0)*y*z*(7*z*z-3*r2)     ! lm=20 (4 -1)
    case(21) ; Ylm=(35*z**4-30*z*z*r2+3.d0*r2*r2)/8.d0   ! lm=21 (4 0)
    case(22) ; Ylm=-sqrt(5.d0/8.d0)*x*z*(7*z**2-3*r2)    ! lm=22 (4 1)
    case(23) ; Ylm=sqrt(5.d0)/4.d0*(7*z*z-r2)*(x*x-y*y)  !lm=23 (4 2)
    case(24) ; Ylm=-sqrt(35.d0/8.d0)*x*z*(x*x-3*y*y)     ! lm=24 (4 3)
    case(25) ; Ylm=sqrt(35.d0)/8.d0*(x**4+y**4-6*x*x*y*y)!lm=25 (4 4)
    end select
    
    return
  End Function Ylm
  !--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
  Function dYlm(x,y,z,il,im,idir)
    implicit none
    real(8),intent(IN) :: x,y,z
    integer,intent(IN) :: il,im,idir
    integer :: ilm
    real(8) :: dYlm
    
    ilm=il*il+il+1+im
    if ((ilm > 9) .or. (il >= 3)) then
      write(*,*) 'dYlm routine not prepared for il>=3'
      stop
    endif
    dYlm=0.d0
    if(idir.eq.1.and.ilm.eq.4) dYlm=-1.d0
    if(idir.eq.2.and.ilm.eq.2) dYlm=-1.d0
    if(idir.eq.3.and.ilm.eq.3) dYlm=1.d0
    if(idir.eq.1.and.ilm.eq.5) dYlm=sqrt(3.d0)*y  ! lm=5 (2 -2)
    if(idir.eq.1.and.ilm.eq.7) dYlm=-x            ! lm=7 (2  0)
    if(idir.eq.1.and.ilm.eq.8) dYlm=-sqrt(3.d0)*z ! lm=8 (2  1)
    if(idir.eq.1.and.ilm.eq.9) dYlm=sqrt(3.d0)*x  ! lm=9 (2  2)
    if(idir.eq.2.and.ilm.eq.5) dYlm=sqrt(3.d0)*x  ! lm=5 (2 -2)
    if(idir.eq.2.and.ilm.eq.6) dYlm=-sqrt(3.d0)*z ! lm=6 (2 -1)
    if(idir.eq.2.and.ilm.eq.7) dYlm=-y            ! lm=7 (2  1)
    if(idir.eq.2.and.ilm.eq.9) dYlm=-sqrt(3.d0)*y ! lm=9 (2  2)
    if(idir.eq.3.and.ilm.eq.6) dYlm=-sqrt(3.d0)*y ! lm=6 (2 -1)
    if(idir.eq.3.and.ilm.eq.7) dYlm=2.d0*z        ! lm=7 (2  0)
    if(idir.eq.3.and.ilm.eq.8) dYlm=-sqrt(3.d0)*x ! lm=8 (2  1)
    
    return
  End Function dYlm
  !--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
  



