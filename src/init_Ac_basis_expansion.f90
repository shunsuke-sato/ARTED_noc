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
subroutine init_Ac_basis_expansion
  use global_variables
  implicit none
  integer :: iter
  real(8) :: tt,xx
  real(8) :: f0_1,f0_2,omega_1,omega_2,tpulse_1,tpulse_2,T1_T2
  integer :: nt_exp, nlines, io, it
  real(8),allocatable :: tt_exp(:), Eexp(:), Aexp(:)
  real(8),allocatable :: Eexp_org(:)
  real(8) :: Eexp_max, texp_ave, t_exp_sigma, cut_sigma, ss
  real(8) :: Ac_tmp
  real(8),allocatable :: tt_exp_pump(:), Aexp_pump(:), Aexp_pump_org(:)
  real(8),allocatable :: Eexp_pump(:)
  real(8) :: Aexp_pump_max, Eexp_pump_max, texp_ave_pump, t_exp_sigma_pump
  real(8) :: cut_sigma_pump
  integer :: nt_exp_pump
  real(8),parameter :: t_offset = 0d0/0.024189d0, tpump_center = 25d0/0.024189d0
  real(8) :: amat(3,3), bmat(3,3), yvec(3), xvec(3)
  real(8) :: x0, x1, x2, y0, y1, y2

  if(myrank == 0)write(*,"(A)")"== Start: Initialization of vector potential."

  allocate(Actot_BE(0:Nt+2),javt_BE(0:Nt+1,3))

  f0_1=5.338d-9*sqrt(IWcm2_1)      ! electric field in a.u.
  omega_1=omegaev_1/(2d0*Ry)  ! frequency in a.u.
  tpulse_1=tpulsefs_1/0.024189d0 ! pulse duration in a.u.
  f0_2=5.338d-9*sqrt(IWcm2_2)      ! electric field in a.u.
  omega_2=omegaev_2/(2d0*Ry)  ! frequency in a.u.
  tpulse_2=tpulsefs_2/0.024189d0 ! pulse duration in a.u.
  T1_T2=T1_T2fs/0.024189d0 ! pulse duration in a.u.
  javt_BE=0.d0
  Actot_BE = 0d0 

! 'cos2cos'                          ! laser_type 'cos2cos', 'cos4cos' or 'impulse'
  select case(laser_type)
  case('impulse')
    Actot_BE = dAc
  case('cos2cos')
! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
    do iter=0,Nt+2
      tt=iter*dt
      xx = tt - 0.5d0*tpulse_1
      if (abs(xx)<0.5d0*tpulse_1) then
        Actot_BE(iter)=-f0_1/omega_1*(cos(pi*xx/tpulse_1))**2*sin(omega_1*(1d0+chirp_1*xx)*xx+phi_CEP_1*2d0*pi)
      end if
    enddo
! probe laser
    do iter=0,Nt+2
      tt=iter*dt
      xx = tt - 0.5d0*tpulse_1 - T1_T2
      if (abs(xx)<0.5d0*tpulse_2) then
        Actot_BE(iter)=Actot_BE(iter) &
          &-f0_2/omega_2*(cos(pi*xx/tpulse_2))**2*sin(omega_2*(1d0+chirp_2*xx)*xx+phi_CEP_2*2d0*pi)
      endif
    enddo
  case('cos4cos')
! pulse shape : A(t)=f0/omega*sin(Pi t/T)**4 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
    do iter=0,Nt+2
      tt=iter*dt
      if (tt<tpulse_1) then
        Actot_BE(iter)=-f0_1/omega_1*(cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1))**4 &
        *sin(omega_1*(1d0+chirp_1*(tt-0.5d0*tpulse_1))*(tt-0.5d0*tpulse_1)+phi_CEP_1*2d0*pi)
      end if
    enddo
! probe laser
    do iter=0,Nt+2
      tt=iter*dt
      if ( (tt-0.5d0*tpulse_1 - T1_T2 >-0.5*tpulse_2) .and. (tt-0.5d0*tpulse_1 - T1_T2 < 0.5*tpulse_2) ) then
        Actot_BE(iter)=Actot_BE(iter) &
          &-f0_2/omega_2*(cos(pi*(tt-(0.5d0*tpulse_1+T1_T2))/tpulse_2))**4&
          *sin(omega_2*(1d0+chirp_2*(tt-(0.5d0*tpulse_1+T1_T2)))*(tt-(0.5d0*tpulse_1+T1_T2))&
          +phi_CEP_2*2d0*pi)
      endif
    enddo
  case('cos_2_4')
! pulse shape : A(t)=f0/omega*sin(Pi t/T)**4 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
    do iter=0,Nt+2
      tt=iter*dt
      if (tt<tpulse_1) then
        Actot_BE(iter)=-f0_1/omega_1*(cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1))**2&
            *sin(omega_1*(1d0+chirp_1*(tt-0.5d0*tpulse_1))*(tt-0.5d0*tpulse_1)+phi_CEP_1*2d0*pi)
      end if
    enddo
! probe laser
    do iter=0,Nt+2
      tt=iter*dt
      if ( (tt-0.5d0*tpulse_1 - T1_T2 >-0.5*tpulse_2) .and. (tt-0.5d0*tpulse_1 - T1_T2 < 0.5*tpulse_2) ) then
        Actot_BE(iter)=Actot_BE(iter) &
          &-f0_2/omega_2*(cos(pi*(tt-(0.5d0*tpulse_1+T1_T2))/tpulse_2))**4&
          *sin(omega_2*(1d0+chirp_2*(tt-(0.5d0*tpulse_1+T1_T2)))*(tt-(0.5d0*tpulse_1+T1_T2)) &
          +phi_CEP_2*2d0*pi)
      endif
    enddo
  case('ge_exp')

! probe pulse
    if(myrank == 0)then
      nlines = 0
      open(131,file="exp_probe_field.dat")
      read(131,*)
      nlines = nlines + 1
      do
        read(131,*,iostat=io)
        IF (io/=0) EXIT
        nlines = nlines + 1
      end do
      close(131)
      nt_exp = nlines -1
    end if
    call MPI_BCAST(nt_exp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)    
    allocate(tt_exp(nt_exp))
    allocate(Eexp(nt_exp))
    allocate(Eexp_org(nt_exp))
    allocate(Aexp(nt_exp))
    if(myrank == 0)then
      open(131,file="exp_probe_field.dat")
      read(131,*)
      do it = 1, nt_exp
        read(131,*)tt_exp(it),Eexp(it)
      end do
      close(131)
      tt_exp = tt_exp/0.024189d0
    end if
    call MPI_BCAST(tt_exp,nt_exp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Eexp,nt_exp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    Eexp_max = maxval(abs(Eexp))
    Eexp = (Eexp/Eexp_max)

    texp_ave = sum(Eexp**2*tt_exp)/sum(Eexp**2)
    t_exp_sigma = sum(Eexp**2*(tt_exp-texp_ave)**2)/sum(Eexp**2)
    if(myrank == 0)write(*,*)'t_exp_sigma',sqrt(t_exp_sigma)*0.024189d0

    Eexp = Eexp*f0_2
    Eexp_org = Eexp

    cut_sigma = (sqrt(t_exp_sigma)*4d0)**10
    Eexp = Eexp * exp(-0.5d0*(tt_exp-texp_ave)**10/cut_sigma)

    Aexp = 0d0
    ss = 0d0
    ss = 0.5d0*Eexp(1)*(tt_exp(2)-tt_exp(1))
    do it = 2, nt_exp
      ss = ss + 0.5d0*Eexp(it)*(tt_exp(it)-tt_exp(it-1))
      Aexp(it) = -ss
      ss = ss + 0.5d0*Eexp(it)*(tt_exp(it)-tt_exp(it-1))
    end do

    Aexp = Aexp * exp(-0.5d0*(tt_exp-texp_ave)**10/cut_sigma)
    
    if(myrank == 0)then
      open(132,file="test_field.dat")
      do it = 1, nt_exp
        write(132,"(999e26.16e3)")tt_exp(it), Eexp_org(it), Eexp(it), aexp(it)
      end do
      close(132)

      open(132,file="test_ac_field.dat")
      do it = 1, nt_exp-1
        write(132,"(999e26.16e3)")0.5d0*(tt_exp(it)+tt_exp(it+1)) &
            ,-(Aexp(it+1)-Aexp(it))/(tt_exp(it+1)-tt_exp(it))
      end do
      close(132)
    end if


    if(myrank == 0)then
      nlines = 0
      open(131,file="exp_pump_A_field.dat")
      read(131,*)
      nlines = nlines + 1
      do
        read(131,*,iostat=io)
        IF (io/=0) EXIT
        nlines = nlines + 1
      end do
      close(131)
      nt_exp_pump = nlines -1
    end if
    call MPI_BCAST(nt_exp_pump,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)    
    allocate(tt_exp_pump(nt_exp_pump), Aexp_pump(nt_exp_pump), Aexp_pump_org(nt_exp_pump))
    if(myrank == 0)then
      open(131,file="exp_pump_A_field.dat")
      read(131,*)
      do it = 1, nt_exp_pump
        read(131,*)tt_exp_pump(it),Aexp_pump(it)
!        write(*,*)tt_exp_pump(it),Aexp_pump(it)
      end do
      close(131)
      tt_exp_pump = tt_exp_pump/0.024189d0
    end if
    call MPI_BCAST(tt_exp_pump,nt_exp_pump,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Aexp_pump,nt_exp_pump,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    Aexp_pump_max = maxval(abs(Aexp_pump))
    Aexp_pump = Aexp_pump/Aexp_pump_max
    texp_ave_pump = sum(Aexp_pump**2*tt_exp_pump)/sum(Aexp_pump**2)
    t_exp_sigma_pump = sum(Aexp_pump**2*(tt_exp_pump-texp_ave_pump)**2)/sum(Aexp_pump**2)
    if(myrank == 0)write(*,*)'t_exp_sigma_pump',sqrt(t_exp_sigma_pump)*0.024189d0
    cut_sigma_pump = (sqrt(t_exp_sigma_pump)*6d0)

    allocate(Eexp_pump(nt_exp_pump-1))
    do it = 1, nt_exp_pump-1
      Eexp_pump(it) = -(Aexp_pump(it+1)-Aexp_pump(it))/(tt_exp_pump(it+1)-tt_exp_pump(it))
    end do
    Eexp_pump_max = maxval(abs(Eexp_pump))
!    write(*,*)"Eexp_pump_max",Eexp_pump_max
    Aexp_pump = Aexp_pump*f0_1/Eexp_pump_max
    Aexp_pump_org = Aexp_pump

! window (24-25 fs)
    do it = 1, nt_exp_pump
      xx = abs(tt_exp_pump(it) - texp_ave_pump)*0.024189d0
      if(xx > tpump_center*0.024189d0)then
        Aexp_pump(it) = 0d0
      else if(tpump_center*0.024189d0 -1d0 <xx)then
        ss = xx-(tpump_center*0.024189d0 -1d0)
        Aexp_pump(it) = Aexp_pump(it) * cos(0.5d0*pi*ss)**2
      end if
    end do

!    Aexp_pump = Aexp_pump * exp(-0.5d0*((tt_exp_pump-texp_ave_pump)/cut_sigma_pump)**12)

    if(myrank == 0)then    
      open(132,file="test_ac_pump_field.dat")
      do it = 1, nt_exp_pump
        write(132,"(999e26.16e3)")tt_exp_pump(it), Aexp_pump(it), Aexp_pump_org(it)
!        write(*,"(999e26.16e3)")tt_exp_pump(it), Aexp_pump(it), Aexp_pump_org(it)
      end do
      close(132)
    end if
!    stop

! pulse shape : A(t)=f0/omega*sin(Pi t/T)**4 *cos (omega t+phi_CEP*2d0*pi) 
! pump laser
    do iter=0,Nt+2
      tt=iter*dt - t_offset
      xx = tt  - tpump_center
      if (abs(xx) < tpump_center) then
        ss = xx +  texp_ave_pump
        do it = 1, nt_exp_pump
          if(ss < tt_exp_pump(it))then
            x0 = tt_exp_pump(it-1)-tt_exp_pump(it-1)
            x1 = tt_exp_pump(it)-tt_exp_pump(it-1)
            x2 = tt_exp_pump(it+1)-tt_exp_pump(it-1)
            y0 = Aexp_pump(it-1)
            y1 = Aexp_pump(it)
            y2 = Aexp_pump(it+1)
            if(it-1 <= 0 .or. nt_exp_pump < it +1)then
              write(*,*)"Eror",it
            end if

            amat(1,1) = x0**2; amat(1,2) = x0; amat(1,3) = 1d0
            amat(2,1) = x1**2; amat(2,2) = x1; amat(2,3) = 1d0
            amat(3,1) = x2**2; amat(3,2) = x2; amat(3,3) = 1d0
            call sub_mat_inv(amat, bmat)
            yvec(1) = y0; yvec(2) = y1; yvec(3) = y2

            xvec = matmul(bmat, yvec)

            Ac_tmp = xvec(1)*(ss-tt_exp_pump(it-1))**2 &
                   + xvec(2)*(ss-tt_exp_pump(it-1)) &
                   + xvec(3)

!            Ac_tmp = Aexp_pump(it) !debug
            if(myrank == 0)write(*,*)it, Ac_tmp

            Actot_BE(iter) = Ac_tmp
            exit
          end if
        end do


!        Actot_BE(iter)=-f0_1/omega_1*(cos(pi*(tt-0.5d0*tpulse_1)/tpulse_1))**2&
!            *sin(omega_1*(1d0+chirp_1*(tt-0.5d0*tpulse_1))*(tt-0.5d0*tpulse_1)+phi_CEP_1*2d0*pi)
      end if
    enddo
! probe laser
    do iter=0,Nt+2
      tt=iter*dt - t_offset
 !     ss = tt -0.5d0*tpulse_1 + texp_ave - T1_T2]
     ss = tt - tpump_center - T1_T2
      if(ss > tt_exp(1) .and. ss < tt_exp(nt_exp))then
        do it = 1, nt_exp
          if(ss < tt_exp(it))then
            
            Ac_tmp = Aexp(it-1)*(tt_exp(it)-ss)/(tt_exp(it)-tt_exp(it-1)) &
                + Aexp(it)*(ss-tt_exp(it-1))/(tt_exp(it)-tt_exp(it-1))
            exit
          end if
        end do
        Actot_BE(iter)=Actot_BE(iter) + Ac_tmp
      end if
    end do
  case default
    err_message='error in init_Ac'
    call err_finalize
  end select

  write(*,*)"maxval(abs(Actot_BE))",maxval(abs(Actot_BE)) ! debug
  if(myrank == 0)write(*,"(A)")"== End: Initialization of vector potential."

  return

contains
  subroutine sub_mat_inv(a, b)
    implicit none
    real(8),intent(in) :: a(3,3)
    real(8),intent(out) :: b(3,3)
    real(8) :: detA
    
    detA=a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3) &
        -a(1,3)*a(2,2)*a(3,1)-a(2,3)*a(3,2)*a(1,1)-a(3,3)*a(1,2)*a(2,1)
    
    b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
    b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
    b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
    
    b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
    b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
    b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
    
    b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
    b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
    b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    
    b = b/detA
    
    return
  end subroutine sub_mat_inv
end subroutine init_Ac_basis_expansion
