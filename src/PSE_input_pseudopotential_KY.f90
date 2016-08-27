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
