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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine PSE_input_pseudopotential_YS
  use Global_Variables
  use PSE_variables
  implicit none
  integer,parameter :: Lmax0=4,Nrmax0=50000
  real(8),parameter :: Eps0=1d-10
  integer :: ik,Mr,l,i,j
  real(8) :: rRC(0:Lmax0)
  real(8) :: r1,r2,r3,r4
  real(8) :: vpp(0:Nrmax0,0:Lmax0),upp(0:Nrmax0,0:Lmax0)   !zero in radial index for taking derivative
  real(8) :: dvpp(0:Nrmax0,0:Lmax0),dupp(0:Nrmax0,0:Lmax0) !zero in radial index for taking derivative
  real(8) :: rhor_nlcc(0:Nrmax0,0:2)   !zero in radial index for taking derivative
  character(2) :: atom_symbol
  character(50) :: ps_file
  character(10) :: ps_postfix

  if(myrank == 0)write(*,*)'start input_pseudopotential_KY'

  allocate(Rps(NE),NRps(NE))
  allocate(Zps(NE),NRloc(NE),Rloc(NE))
  allocate(Mps(NI),Lref(NE),Mlps(NE))
  allocate(anorm(0:Lmax,NE),inorm(0:Lmax,NE))
  allocate(rad(Nrmax,NE),vloctbl(Nrmax,NE))
  allocate(rho_nlcc_tbl(Nrmax,NE),tau_nlcc_tbl(Nrmax,NE))
  allocate(rho_nlcc(NL),tau_nlcc(NL))
  allocate(radnl(Nrmax,NE))
  allocate(udVtbl(Nrmax,0:Lmax,NE))
  allocate(Vpsl(NL))
  allocate(flag_nlcc(NE))
  flag_nlcc(:) = .false.

  if(Myrank == 0)read(*,*) (Lref(j),j=1,NE)
  call MPI_BCAST(Lref,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (Myrank == 0) then
! --- Making prefix ---
    select case (ps_format)
    case('KY')        ; ps_postfix = '_rps.dat'
    case('ABINIT')    ; ps_postfix = '.pspnc'
    case('ABINITFHI') ; ps_postfix = '.fhi'
    case('FHI')       ; ps_postfix = '.cpi'
!    case('ATOM')      ; ps_postfix = '.psf' !Not implemented yet
    case default ; stop 'Unprepared ps_format is required input_pseudopotential_YS'
    end select

! --- input pseudopotential and wave function ---
    do ik=1,NE
      select case (Zatom(ik))
      case (1) ; atom_symbol = 'H ' !; Mass(ik)=1.d0
      case (2) ; atom_symbol = 'He' !; Mass(ik)=4.d0
      case (3) ; atom_symbol = 'Li' !; Mass(ik)=7.d0
      case (4) ; atom_symbol = 'Be' !; Mass(ik)=9.d0
      case (5) ; atom_symbol = 'B ' !; Mass(ik)=11.d0
      case (6) ; atom_symbol = 'C ' !; Mass(ik)=12.d0
      case (7) ; atom_symbol = 'N ' !; Mass(ik)=14.d0
      case (8) ; atom_symbol = 'O ' !; Mass(ik)=16.d0
      case (9) ; atom_symbol = 'F ' !; Mass(ik)=19.d0
      case(10) ; atom_symbol = 'Ne' !; Mass(ik)=20.d0
      case(11) ; atom_symbol = 'Na' !; Mass(ik)=23.d0
      case(12) ; atom_symbol = 'Mg' !; Mass(ik)=24.d0
      case(13) ; atom_symbol = 'Al' !; Mass(ik)=27.d0
      case(14) ; atom_symbol = 'Si' !; Mass(ik)=28.d0
      case(15) ; atom_symbol = 'P ' !; Mass(ik)=31.d0
      case(16) ; atom_symbol = 'S ' !; Mass(ik)=32.d0
      case(17) ; atom_symbol = 'Cl' !; Mass(ik)=35.d0
      case(18) ; atom_symbol = 'Al' !; Mass(ik)=40.d0
      case(19) ; atom_symbol = 'K ' !; Mass(ik)=39.d0
      case(20) ; atom_symbol = 'Ca' !; Mass(ik)=40.d0
      case(21) ; atom_symbol = 'Sc' !; Mass(ik)=45.d0
      case(22) ; atom_symbol = 'Ti' !; Mass(ik)=48.d0
      case(23) ; atom_symbol = 'V ' !; Mass(ik)=51.d0
      case(24) ; atom_symbol = 'Cr' !; Mass(ik)=52.d0
      case(25) ; atom_symbol = 'Mn' !; Mass(ik)=55.d0
      case(26) ; atom_symbol = 'Fe' !; Mass(ik)=56.d0
      case(27) ; atom_symbol = 'Co' !; Mass(ik)=59.d0
      case(28) ; atom_symbol = 'Ni' !; Mass(ik)=59.d0
      case(29) ; atom_symbol = 'Cu' !; Mass(ik)=63.d0
      case(30) ; atom_symbol = 'Zn' !; Mass(ik)=65.d0
      case(31) ; atom_symbol = 'Ga' !; Mass(ik)=69.d0
      case(32) ; atom_symbol = 'Ge' !; Mass(ik)=73.d0
      case(33) ; atom_symbol = 'As' !; Mass(ik)=75.d0
      case(34) ; atom_symbol = 'Se' !; Mass(ik)=79.d0
      case(35) ; atom_symbol = 'Br' !; Mass(ik)=80.d0
      case(36) ; atom_symbol = 'Kr' !; Mass(ik)=84.d0
      case(37) ; atom_symbol = 'Rb' !; Mass(ik)=85.d0
      case(38) ; atom_symbol = 'Sr' !; Mass(ik)=88.d0
      case(39) ; atom_symbol = 'Y ' !; Mass(ik)=89.d0
      case(40) ; atom_symbol = 'Zr' !; Mass(ik)=91.d0
      case(41) ; atom_symbol = 'Nb' !; Mass(ik)=93.d0
      case(42) ; atom_symbol = 'Mo' !; Mass(ik)=96.d0
      case(43) ; atom_symbol = 'Tc' !; Mass(ik)=98.d0
      case(44) ; atom_symbol = 'Ru' !; Mass(ik)=101.d0
      case(45) ; atom_symbol = 'Rh' !; Mass(ik)=103.d0
      case(46) ; atom_symbol = 'Pd' !; Mass(ik)=106.d0
      case(47) ; atom_symbol = 'Ag' !; Mass(ik)=108.d0
      case(48) ; atom_symbol = 'Cd' !; Mass(ik)=112.d0
      case(49) ; atom_symbol = 'In' !; Mass(ik)=115.d0
      case(50) ; atom_symbol = 'Sn' !; Mass(ik)=119.d0
      case(51) ; atom_symbol = 'Sb' !; Mass(ik)=122.d0
      case(52) ; atom_symbol = 'Te' !; Mass(ik)=128.d0
      case(53) ; atom_symbol = 'I ' !; Mass(ik)=127.d0
      case(54) ; atom_symbol = 'Xe' !; Mass(ik)=131.d0
      case(55) ; atom_symbol = 'Cs' !; Mass(ik)=133.d0
      case(56) ; atom_symbol = 'Ba' !; Mass(ik)=137.d0
      case(57) ; atom_symbol = 'La' !; Mass(ik)=139.d0
      case(58) ; atom_symbol = 'Ce' !; Mass(ik)=140.d0
      case(59) ; atom_symbol = 'Pr' !; Mass(ik)=141.d0
      case(60) ; atom_symbol = 'Nd' !; Mass(ik)=144.d0
      case(61) ; atom_symbol = 'Pm' !; Mass(ik)=145.d0
      case(62) ; atom_symbol = 'Sm' !; Mass(ik)=150.d0
      case(63) ; atom_symbol = 'Eu' !; Mass(ik)=152.d0
      case(64) ; atom_symbol = 'Gd' !; Mass(ik)=157.d0
      case(65) ; atom_symbol = 'Tb' !; Mass(ik)=159.d0
      case(66) ; atom_symbol = 'Dy' !; Mass(ik)=164.d0
      case(67) ; atom_symbol = 'Ho' !; Mass(ik)=165.d0
      case(68) ; atom_symbol = 'Er' !; Mass(ik)=167.d0
      case(69) ; atom_symbol = 'Tm' !; Mass(ik)=169.d0
      case(70) ; atom_symbol = 'Yb' !; Mass(ik)=173.d0
      case(71) ; atom_symbol = 'Lu' !; Mass(ik)=175.d0
      case(72) ; atom_symbol = 'Hf' !; Mass(ik)=178.d0
      case(73) ; atom_symbol = 'Ta' !; Mass(ik)=181.d0
      case(74) ; atom_symbol = 'W ' !; Mass(ik)=184.d0
      case(75) ; atom_symbol = 'Re' !; Mass(ik)=186.d0
      case(76) ; atom_symbol = 'Os' !; Mass(ik)=190.d0
      case(77) ; atom_symbol = 'Ir' !; Mass(ik)=192.d0
      case(78) ; atom_symbol = 'Pt' !; Mass(ik)=195.d0
      case(79) ; atom_symbol = 'Au' !; Mass(ik)=197.d0
      case(80) ; atom_symbol = 'Hg' !; Mass(ik)=201.d0
      case(81) ; atom_symbol = 'Tl' !; Mass(ik)=204.d0
      case(82) ; atom_symbol = 'Pb' !; Mass(ik)=207.d0
      case(83) ; atom_symbol = 'Bi' !; Mass(ik)=209.d0
      case default ; stop 'Unprepared atomic data is called input_pseudopotential_YS'
      end select

      ps_file=trim(atom_symbol)//trim(ps_postfix)

      write(*,*) '===================pseudopotential data==================='
      write(*,*) 'ik ,atom_symbol=',ik, atom_symbol
      write(*,*) 'ps_format =',ps_format
      write(*,*) 'ps_file =',ps_file
      select case (ps_format)
      case('KY')        ; call Read_PS_KY(Lmax0,Nrmax0,Mr,rRC,upp,vpp,ik,ps_file)
      case('ABINIT')    ; call Read_PS_ABINIT(Lmax0,Nrmax0,Mr,rRC,upp,vpp,ik,ps_file)
      case('ABINITFHI') ; call Read_PS_ABINITFHI(Lmax0,Nrmax0,Mr,rRC,upp,vpp,rhor_nlcc,ik,ps_file)
      case('FHI')       ; call Read_PS_FHI(Lmax0,Nrmax0,Mr,rRC,upp,vpp,ik,ps_file)
!      case('ATOM')      ; call Read_PS_ATOM
      case default ; stop 'Unprepared ps_format is required input_pseudopotential_YS'
      end select
! Set meaning domain in the arrays 
      Rps(ik)=maxval(rRC(0:Mlps(ik)))
      do i=1,Nrmax
        if(rad(i,ik).gt.Rps(ik)) exit
      enddo
      NRps(ik)=i
      if(NRps(ik).ge.Nrmax) stop 'NRps>Nrmax at input_pseudopotential_YS'
      NRloc(ik)=NRps(ik)
      Rloc(ik)=Rps(ik)
      radnl(:,ik)=rad(:,ik)
      do l=0,Mlps(ik)
        anorm(l,ik) = 0.d0
        do i=1,Mr-1
          r1 = rad(i+1,ik)-rad(i,ik)
          anorm(l,ik) = anorm(l,ik)  &
               + (upp(i,l)**2*(vpp(i,l)-vpp(i,Lref(ik)))+upp(i-1,l)**2*(vpp(i-1,l)-vpp(i-1,Lref(ik))))*r1
        end do
        anorm(l,ik) = 0.5d0*anorm(l,ik)
        inorm(l,ik)=+1
        if(anorm(l,ik).lt.0.d0) then
          anorm(l,ik)=-anorm(l,ik)
          inorm(l,ik)=-1
        endif
        if(abs(anorm(l,ik)).lt.Eps0) inorm(l,ik)=0
        anorm(l,ik)=sqrt(anorm(l,ik))
      enddo

      write(*,*) 'Zps(ik), Mlps(ik) =',Zps(ik), Mlps(ik)
      write(*,*) 'Rps(ik), NRps(ik) =',Rps(ik), NRps(ik)
      write(*,*) 'Lref(ik) =',Lref(ik)
      write(*,*) 'anorm(ik,l) =',(anorm(l,ik),l=0,Mlps(ik))
      write(*,*) 'inorm(ik,l) =',(inorm(l,ik),l=0,Mlps(ik))
!      write(*,*) 'Mass(ik) =',Mass(ik)
      write(*,*) '=========================================================='

      call Making_PS_without_masking 

    enddo

  endif


  CALL MPI_BCAST(Zps,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Mlps,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Rps,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NRps,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NRloc,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Rloc,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(anorm,(Lmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inorm,(Lmax+1)*NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rad,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(radnl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(vloctbl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rho_nlcc_tbl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(tau_nlcc_tbl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(dvloctbl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(udVtbl,Nrmax*(Lmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(dudVtbl,Nrmax*(Lmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(Mass,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  return
  contains
!====
    Subroutine Making_PS_without_masking
      implicit none

! multiply sqrt((2l+1)/4pi)/r**(l+1) for radial w.f.
      do l=0,Mlps(ik)
        do i=1,Mr
          upp(i,l)=upp(i,l)*sqrt((2*l+1.d0)/(4*pi))/(rad(i+1,ik))**(l+1)
        enddo
        upp(0,l)=upp(1,l)
!        upp(0,l)=2.d0*upp(1,l)-upp(2,l)
      enddo

      do l=0,Mlps(ik)
        do i=1,Mr-1
          r1 = rad(i+1,ik)-rad(i,ik)
          r2 = rad(i+1,ik)-rad(i+2,ik)
          r3 = rad(i+2,ik)-rad(i,ik)
          r4 = r1/r2
          dvpp(i,l)=(r4+1.d0)*(vpp(i,l)-vpp(i-1,l))/r1-(vpp(i+1,l)-vpp(i-1,l))/r3*r4
          dupp(i,l)=(r4+1.d0)*(upp(i,l)-upp(i-1,l))/r1-(upp(i+1,l)-upp(i-1,l))/r3*r4
        end do
        dvpp(0,l)=dvpp(1,l)
        dvpp(Mr,l)=dvpp(Mr-1,l)
        dupp(0,l)=dupp(1,l)
        dupp(Mr,l)=dupp(Mr-1,l)
!        dvpp(0,l)=2.d0*dvpp(1,l)-dvpp(2,l)
!        dvpp(Mr,l)=2.d0*dvpp(Mr-1,l)-dvpp(Mr-2,l)
!        dupp(0,l)=2.d0*dupp(1,l)-dupp(2,l)
!        dupp(Mr,l)=2.d0*dupp(Mr-1,l)-dupp(Mr-2,l)
      end do

      do l=0,Mlps(ik)
        do i=1,NRps(ik)
          vloctbl(i,ik)=vpp(i-1,Lref(ik))
!          dvloctbl(i,ik)=dvpp(i-1,Lref(ik))
          udVtbl(i,l,ik)=(vpp(i-1,l)-vpp(i-1,Lref(ik)))*upp(i-1,l)
!          dudVtbl(i,l,ik)=(dvpp(i-1,l)-dvpp(i-1,Lref(ik)))*upp(i-1,l) + (vpp(i-1,l)-vpp(i-1,Lref(ik)))*dupp(i-1,l)
        enddo
        if (inorm(l,ik) == 0) cycle
        udVtbl(1:NRps(ik),l,ik)=udVtbl(1:NRps(ik),l,ik)/anorm(l,ik)
!        dudVtbl(1:NRps(ik),l,ik)=dudVtbl(1:NRps(ik),l,ik)/anorm(l,ik)
      enddo

      rho_nlcc_tbl(:,ik)=0d0; tau_nlcc_tbl(:,ik)=0d0
      if(.not.flag_nlcc(ik))return
      do i=1,NRps(ik)
        if(rhor_nlcc(i-1,0)/rhor_nlcc(0,0) < 1d-7)exit
        rho_nlcc_tbl(i,ik)=rhor_nlcc(i-1,0)
        tau_nlcc_tbl(i,ik)=0.25d0*rhor_nlcc(i-1,1)**2/rhor_nlcc(i-1,0)
      end do

      
      return
    End Subroutine Making_PS_without_masking
  End Subroutine PSE_input_pseudopotential_YS
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Read_PS_KY(Lmax0,Nrmax0,Mr,rRC,upp,vpp,ik,ps_file)
  use Global_Variables
  use PSE_variables
  implicit none
!argument
  integer,intent(in) :: Lmax0,Nrmax0,ik
  integer,intent(out) :: Mr
  real(8),intent(out) :: rRC(0:Lmax0)
  real(8),intent(out) :: vpp(0:Nrmax0,0:Lmax0),upp(0:Nrmax0,0:Lmax0)
  character(50),intent(in) :: ps_file
!local variable
  integer :: l,i,irPC
  real(8) :: step,rPC,r,rhopp(0:Nrmax0),rZps

  open(4,file=ps_file,status='old')
  read(4,*) Mr,step,Mlps(ik),rZps
  Zps(ik)=int(rZps+1d-10)
  if(Mr.gt.Nrmax0) stop 'Mr>Nrmax0 at Read_PS_KY'
  if(Mlps(ik).gt.Lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_KY'
  if(Mlps(ik).gt.Lmax) stop 'Mlps(ik)>Lmax at Read_PS_KY'
  read(4,*) irPC,(rRC(l),l=0,Mlps(ik))
  rPC=real(irPC) !Radius for partial core correction: not working in this version
  do i=0,Mr
    read(4,*) r,rhopp(i),(vpp(i,l),l=0,Mlps(ik))
  end do
  do i=0,Mr
    read(4,*) r,(upp(i,l),l=0,Mlps(ik))
  end do
  close(4)
! change to atomic unit
  step=step/a_B
  rRC(0:Mlps(ik))=rRC(0:Mlps(ik))/a_B
  vpp(0:Mr,0:Mlps(ik))=vpp(0:Mr,0:Mlps(ik))/(2*Ry)
  upp(0:Mr,0:Mlps(ik))=upp(0:Mr,0:Mlps(ik))*sqrt(a_B)

  do i=1,Nrmax
    rad(i,ik)=(i-1)*step
  enddo
  return
End Subroutine Read_PS_KY
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Read_PS_ABINIT(Lmax0,Nrmax0,Mr,rRC,upp,vpp,ik,ps_file)
  use Global_Variables
  use PSE_variables
!See http://www.abinit.org/downloads/psp-links/psp-links/lda_tm
  implicit none
!argument
  integer,intent(in) :: Lmax0,Nrmax0,ik
  integer,intent(out) :: Mr
  real(8),intent(out) :: rRC(0:Lmax0)
  real(8),intent(out) :: vpp(0:Nrmax0,0:Lmax0),upp(0:Nrmax0,0:Lmax0)
  character(50),intent(in) :: ps_file
!local variable
  integer :: i
  real(8) :: rZps
  integer :: ll
  real(8) :: zatom0, zion, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well,l
  real(8) :: e99_0,e99_9,nproj,rcpsp,rms,ekb1,ekb2,epsatm,rchrg,fchrg,qchrg
  character(1) :: dummy_text

  open(4,file=ps_file,status='old')
  read(4,*) dummy_text
  read(4,*) zatom0, zion, pspdat
  rZps = zion
  Zps(ik)=int(rZps+1d-10)
  read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
  Mlps(ik)=lmaxabinit
  if(lloc .ne. Lref(ik)) write(*,*) "Warning! Lref(ik=",ik,") is different from intended one in ",ps_file
  Mr = mmax - 1
  if(Mr.gt.Nrmax0) stop 'Mr>Nrmax0 at Read_PS_ABINIT'
  if(Mlps(ik).gt.Lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_ABINIT'
  if(Mlps(ik).gt.Lmax) stop 'Mlps(ik)>Lmax at Read_PS_ABINIT'
  do ll=0,Mlps(ik)
    read(4,*) l,e99_0,e99_9,nproj,rcpsp
    read(4,*) rms,ekb1,ekb2,epsatm
    rRC(ll) = rcpsp
  end do
  read(4,*) rchrg,fchrg,qchrg
  do ll=0,Mlps(ik)
    read(4,*) dummy_text
    do i=1,(Mr+1)/3
      read(4,*) vpp(3*(i-1),ll),vpp(3*(i-1)+1,ll),vpp(3*(i-1)+2,ll)
    end do
  end do
  do ll=0,Mlps(ik)
    read(4,*) dummy_text
    do i=1,(Mr+1)/3
      read(4,*) upp(3*(i-1),ll),upp(3*(i-1)+1,ll),upp(3*(i-1)+2,ll)
    end do
  end do
  close(4)

  do i=0,Nrmax-1
    rad(i+1,ik) = 1.0d2*(dble(i)/dble(mmax-1)+1.0d-2)**5 - 1.0d-8
  end do

  return
End Subroutine Read_PS_ABINIT
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Read_PS_ABINITFHI(Lmax0,Nrmax0,Mr,rRC,upp,vpp,rhor_nlcc,ik,ps_file)
  use Global_Variables
  use PSE_variables
!This is for  FHI pseudopotential listed in abinit web page and not for original FHI98PP.
!See http://www.abinit.org/downloads/psp-links/lda_fhi
  implicit none
!argument
  integer,intent(in) :: Lmax0,Nrmax0,ik
  integer,intent(out) :: Mr
  real(8),intent(out) :: rRC(0:Lmax0)
  real(8),intent(out) :: vpp(0:Nrmax0,0:Lmax0),upp(0:Nrmax0,0:Lmax0)
  real(8),intent(out) :: rhor_nlcc(0:Nrmax0,0:2)
  character(50),intent(in) :: ps_file
!local variable
  character(50) :: temptext
  integer :: i,j
  real(8) :: step,rZps,dummy
  integer :: Mr_l(0:Lmax0),l,ll
  real(8) :: step_l(0:Lmax),rRC_mat(0:Lmax,0:Lmax)=-1.d0

  open(4,file=ps_file,status='old')
  write(*,*) '===================Header of ABINITFHI pseudo potential==================='
  do i=1,7
    read(4,'(a)') temptext
    write(*,*) temptext
  end do
  write(*,*) '===================Header of ABINITFHI pseudo potential==================='
  read(4,*) rZps,Mlps(ik)
  Zps(ik)=int(rZps+1d-10)
  Mlps(ik) = Mlps(ik)-1
  if(Mlps(ik).gt.Lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_FHI'
  if(Mlps(ik).gt.Lmax) stop 'Mlps(ik)>Lmax at Read_PS_FHI'
  do i=1,10
     read(4,*) dummy
  end do
  do l=0,Mlps(ik)
    read(4,*) Mr_l(l),step_l(l)
    if(Mr_l(l).gt.Nrmax0) stop 'Mr>Nrmax0 at Read_PS_FHI'
    do i=1,Mr_l(l)
      read(4,*) j,rad(i+1,ik),upp(i,l),vpp(i,l) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
    end do
    rad(1,ik)=0.d0
    upp(0,l)=0.d0
    vpp(0,l)=vpp(1,l)-(vpp(2,l)-vpp(1,l))/(rad(3,ik)-rad(2,ik))*(rad(2,ik))
  end do
!Nonlinear core-correction
  do i=1,Mr_l(0)
    read(4,*,err=940) rad(i+1,ik),rhor_nlcc(i,0),rhor_nlcc(i,1),rhor_nlcc(i,2)
  end do
  rhor_nlcc(0,:)=rhor_nlcc(1,:)-(rhor_nlcc(2,:)-rhor_nlcc(1,:)) &
    /(rad(3,ik)-rad(2,ik))*(rad(2,ik))
  rhor_nlcc = rhor_nlcc/(4d0*pi)
  flag_nlcc(ik) = .true.
!Nonlinear core-correction
940  close(4)

  if(minval(Mr_l(0:Mlps(ik))).ne.maxval(Mr_l(0:Mlps(ik)))) then
    stop 'Mr are diffrent at Read_PS_FHI'
  else 
    Mr = minval(Mr_l(0:Mlps(ik)))
  end if
  if((maxval(step_l(0:Mlps(ik)))-minval(step_l(0:Mlps(ik)))).ge.1.d-14) then
    stop 'step are different at Read_PS_FHI'
  else 
    step = minval(step_l(0:Mlps(ik)))
  end if

  do i=Mr+1,Nrmax
    rad(i+1,ik) = rad(i,ik)*step
  end do

  do l=0,Mlps(ik)
    do ll=0,Mlps(ik)
      do i=Mr,1,-1
        if(abs(vpp(i,l)-vpp(i,ll)).gt.1.d-10) then
          rRC_mat(l,ll) = rad(i+1+1,ik)
          exit
        end if
      end do
    end do
    rRC(l)=maxval(rRC_mat(l,:))
  end do

  return
End Subroutine Read_PS_ABINITFHI
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Read_PS_FHI(Lmax0,Nrmax0,Mr,rRC,upp,vpp,ik,ps_file)
  use Global_Variables
  use PSE_variables
!This is for original FHI98PP and not for FHI pseudopotential listed in abinit web page
!See http://th.fhi-berlin.mpg.de/th/fhi98md/fhi98PP/
  implicit none
!argument
  integer,intent(in) :: Lmax0,Nrmax0,ik
  integer,intent(out) :: Mr
  real(8),intent(out) :: rRC(0:Lmax0)
  real(8),intent(out) :: vpp(0:Nrmax0,0:Lmax0),upp(0:Nrmax0,0:Lmax0)
  character(50),intent(in) :: ps_file
!local variable
  integer :: i,j
  real(8) :: step,rZps,dummy
  integer :: Mr_l(0:Lmax0),l,ll
  real(8) :: step_l(0:Lmax),rRC_mat(0:Lmax,0:Lmax)=-1.d0

  open(4,file=ps_file,status='old')
  read(4,*) rZps,Mlps(ik)
  Zps(ik)=int(rZps+1d-10)
  Mlps(ik) = Mlps(ik)-1
  if(Mlps(ik).gt.Lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_FHI'
  if(Mlps(ik).gt.Lmax) stop 'Mlps(ik)>Lmax at Read_PS_FHI'
  do i=1,10
     read(4,*) dummy
  end do
  do l=0,Mlps(ik)
    read(4,*) Mr_l(l),step_l(l)
    if(Mr_l(l).gt.Nrmax0) stop 'Mr>Nrmax0 at Read_PS_FHI'
    do i=1,Mr_l(l)
      read(4,*) j,rad(i+1,ik),upp(i,l),vpp(i,l) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
    end do
    rad(1,ik)=0.d0
    upp(0,l)=0.d0
    vpp(0,l)=vpp(1,l)-(vpp(2,l)-vpp(1,l))/(rad(3,ik)-rad(2,ik))*(rad(2,ik))
  end do
  close(4)

  if(minval(Mr_l(0:Mlps(ik))).ne.maxval(Mr_l(0:Mlps(ik)))) then
    stop 'Mr are diffrent at Read_PS_FHI'
  else 
    Mr = minval(Mr_l(0:Mlps(ik)))
  end if
  if((maxval(step_l(0:Mlps(ik)))-minval(step_l(0:Mlps(ik)))).ge.1.d-14) then
    stop 'step are different at Read_PS_FHI'
  else 
    step = minval(step_l(0:Mlps(ik)))
  end if

  do i=Mr+1,Nrmax
    rad(i+1,ik) = rad(i,ik)*step
  end do

  do l=0,Mlps(ik)
    do ll=0,Mlps(ik)
      do i=Mr,1,-1
        if(abs(vpp(i,l)-vpp(i,ll)).gt.1.d-10) then
          rRC_mat(l,ll) = rad(i+1+1,ik)
          exit
        end if
      end do
    end do
    rRC(l)=maxval(rRC_mat(l,:))
  end do

  return
End Subroutine Read_PS_FHI
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
!    Subroutine Read_PS_ATOM !.psf format created by ATOM for SIESTA
!      implicit none
!      return
!    End Subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130



Subroutine PSE_input_pseudopotential_KY
  use Global_Variables
  use PSE_variables
  implicit none
  integer,parameter :: Lmax0=4,Nrmax0=50000
  real(8),parameter :: Eps0=1d-10
  integer :: ik,Mr,l,i,j,irPC
  real(8) :: step,rPC,rRC(0:Lmax0),r,rhopp(0:Nrmax0),rZps
  real(8) :: vpp(0:Nrmax0,0:Lmax0),upp(0:Nrmax0,0:Lmax0)
  real(8) :: dvpp(0:Nrmax0,0:Lmax0),dupp(0:Nrmax0,0:Lmax0)
  character(50) :: ps_file

  if(myrank == 0)write(*,*)'start input_pseudopotential_KY'

  allocate(Rps(NE),NRps(NE))
  allocate(Zps(NE),NRloc(NE),Rloc(NE))
  allocate(Mps(NI),Lref(NE),Mlps(NE))
  allocate(anorm(0:Lmax,NE),inorm(0:Lmax,NE))
  allocate(rad(Nrmax,NE),vloctbl(Nrmax,NE))
  allocate(radnl(Nrmax,NE))
  allocate(udVtbl(Nrmax,0:Lmax,NE))
  allocate(Vpsl(NL))

  if(Myrank == 0)read(*,*) (Lref(j),j=1,NE)
  call MPI_BCAST(Lref,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! --- input pseudopotential and wave function ---
  if(Myrank == 0) then
    do ik=1,NE
      select case (Zatom(ik))
      case (1)
        ps_file='H_rps.dat'
!        Mass(ik)=1.d0
      case (3)
        ps_file='Li_rps.dat'
!        Mass(ik)=7.d0
      case (6)
        ps_file='C_rps.dat'
!        Mass(ik)=12.d0
      case (7)
        ps_file='N_rps.dat'
!        Mass(ik)=14.d0
      case (8)
        ps_file='O_rps.dat'
!        Mass(ik)=16.d0
      case (9)
        ps_file='F_crps.dat'
      case(11)
        ps_file='Na_rps.dat'
!        Mass(ik)=23.d0
      case(13)
        ps_file='Al_rps.dat'
!        Mass(ik)=27.d0
      case(14)
        ps_file='Si_rps.dat'
!        Mass(ik)=28.d0
      case (29)
        ps_file='Cu_rps.dat'
!        Mass(ik)=63.d0
      case (31)
        ps_file='Ga_rps.dat'
      case (33)
        ps_file='As_rps.dat'
      case(51)
        ps_file='Sb_rps.dat'
!        Mass(ik)=122.d0
      case(83)
        ps_file='Bi_rps.dat'
!        Mass(ik)=209.d0
      case default
        err_message='No such pseudo-potential'
        call err_finalize
      end select

      open(4,file=ps_file,status='old')
      read(4,*) Mr,step,Mlps(ik),rZps
      Zps(ik)=int(rZps+1d-10)
      if(Mr.gt.Nrmax0) stop 'Mr>Nrmax0 at input_pseudopotential_KY'
      if(Mlps(ik).gt.Lmax0) stop 'Mlps(ik)>Lmax0 at input_pseudopotential_KY'
      if(Mlps(ik).gt.Lmax) stop 'Mlps(ik)>Lmax at input_pseudopotential_KY'
      read(4,*) irPC,(rRC(l),l=0,Mlps(ik))
      rPC=real(irPC)
      do i=0,Mr
        read(4,*) r,rhopp(i),(vpp(i,l),l=0,Mlps(ik))
      end do
      do i=0,Mr
        read(4,*) r,(upp(i,l),l=0,Mlps(ik))
      end do
      close(4)

! change to atomic unit
      step=step/a_B
      rRC(0:Mlps(ik))=rRC(0:Mlps(ik))/a_B
      vpp(0:Mr,0:Mlps(ik))=vpp(0:Mr,0:Mlps(ik))/(2*Ry)
      upp(0:Mr,0:Mlps(ik))=upp(0:Mr,0:Mlps(ik))*sqrt(a_B)

      do l=0,Mlps(ik)
        anorm(l,ik)=sum(upp(0:Mr,l)**2*(vpp(0:Mr,l)-vpp(0:Mr,Lref(ik)))*step)
        inorm(l,ik)=+1
        if(abs(anorm(l,ik)).lt.Eps0) then
          inorm(l,ik)=0
        else 
          if(anorm(l,ik).lt.0.d0) then
            anorm(l,ik)=-anorm(l,ik)
            inorm(l,ik)=-1
          endif
        endif
        anorm(l,ik)=sqrt(anorm(l,ik))
      enddo
! multiply sqrt((2l+1)/4pi)/r**(l+1) for radial w.f.
      do l=0,Mlps(ik)
        do i=1,Mr
          upp(i,l)=upp(i,l)*sqrt((2*l+1.d0)/(4*pi))/(i*step)**(l+1)
        enddo
        upp(0,l)=upp(1,l)
      enddo

!      do l=0,Mlps(ik)
!        do i=1,Mr-1
!          dvpp(i,l)=(vpp(i+1,l)-vpp(i-1,l))/(2d0*step)
!          dupp(i,l)=(upp(i+1,l)-upp(i-1,l))/(2d0*step)
!        end do
!        dvpp(0,l)=dvpp(1,l)
!        dvpp(Mr,l)=dvpp(Mr-1,l)
!        dupp(0,l)=dupp(1,l)
!        dupp(Mr,l)=dupp(Mr-1,l)
!      end do

      Rps(ik)=maxval(rRC(0:Mlps(ik)))
      do i=1,Nrmax
        rad(i,ik)=(i-1)*step
        if(rad(i,ik).gt.Rps(ik)) exit
      enddo
      NRps(ik)=i
      if(NRps(ik).ge.Nrmax) stop 'NRps>Nrmax at input_pseudopotential_KY'
      NRloc(ik)=NRps(ik)
      Rloc(ik)=Rps(ik)
      radnl(:,ik)=rad(:,ik)

      do l=0,Mlps(ik)
        do i=1,NRps(ik)
          vloctbl(i,ik)=vpp(i-1,Lref(ik))
!          dvloctbl(i,ik)=dvpp(i-1,Lref(ik))
          udVtbl(i,l,ik)=(vpp(i-1,l)-vpp(i-1,Lref(ik)))*upp(i-1,l)
!          dudVtbl(i,l,ik)=(dvpp(i-1,l)-dvpp(i-1,Lref(ik)))*upp(i-1,l) + (vpp(i-1,l)-vpp(i-1,Lref(ik)))*dupp(i-1,l)
        enddo
        if (inorm(l,ik) == 0) cycle
        udVtbl(1:NRps(ik),l,ik)=udVtbl(1:NRps(ik),l,ik)/anorm(l,ik)
!        dudVtbl(1:NRps(ik),l,ik)=dudVtbl(1:NRps(ik),l,ik)/anorm(l,ik)
      enddo
    enddo

  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  CALL MPI_BCAST(Zps,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Mlps,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Rps,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NRps,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(NRloc,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Rloc,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(anorm,(Lmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(inorm,(Lmax+1)*NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rad,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(radnl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(vloctbl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(dvloctbl,Nrmax*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(udVtbl,Nrmax*(Lmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(dudVtbl,Nrmax*(Lmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!  CALL MPI_BCAST(Mass,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  if(myrank == 0)write(*,*)'input_pseudopotential_KY is complete'

  return
End Subroutine PSE_input_pseudopotential_KY
