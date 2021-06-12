!
!  Copyright 2018 S.A. Sato
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
subroutine dt_evolve_macro_field
  use global_variables
  use ms_maxwell_ks_variables
  implicit none


  call dt_evolve_macro_field_explicit_multi_step
!  call dt_evolve_macro_field_explicit
!  call dt_evolve_macro_field_implicit

contains
!--------------------------------------------------------------
  subroutine dt_evolve_macro_field_explicit_multi_step
    implicit none
    real(8) :: clap0, clap1
    real(8) :: Lap_Ac(nx_s:nx_e)
    real(8) :: jt_tmp(Mx), af(Mx), bf(Mx), cf(Mx)
    integer :: ix
    integer :: iter
    real(8) :: xx
    
    clap0 = -2d0/dx_m**2
    clap1 = 1d0/dx_m**2

    af = (jt_m-2d0*jt_old_m+jt_old2_m)/dt**2
    bf = 0.5d0*(3d0*jt_m-4d0*jt_old_m+jt_old2_m)/dt
    cf = jt_m

    do iter = 0, nt_internal_m-1
      xx = dt_m*iter
      jt_tmp = 0.5d0*af*xx**2+bf*xx+cf


      Lap_Ac(nx_s) = clap0*Ac_m(nx_s) + clap1*Ac_m(nx_s+1)
      do ix = nx_s+1,nx_e-1
        Lap_Ac(ix) = clap0*Ac_m(ix) + clap1*(Ac_m(ix+1)+Ac_m(ix-1))
      end do
      Lap_Ac(nx_e) = clap0*Ac_m(nx_e) + clap1*Ac_m(nx_e-1)

      Ac_m_n = 2d0*Ac_m - Ac_m_o + (dt_m*clight)**2*Lap_Ac
      Ac_m_n(1:Mx) = Ac_m_n(1:Mx) -4d0*pi*dt_m**2*jt_tmp(1:Mx)

      Ac_m_o = Ac_m
      Ac_m = Ac_m_n
      
    end do


  end subroutine dt_evolve_macro_field_explicit_multi_step
!--------------------------------------------------------------
  subroutine dt_evolve_macro_field_explicit
    implicit none
    real(8) :: clap0, clap1
    real(8) :: Lap_Ac(nx_s:nx_e)
    integer :: ix
    
    clap0 = -2d0/dx_m**2
    clap1 = 1d0/dx_m**2

    Lap_Ac(nx_s) = clap0*Ac_m(nx_s) + clap1*Ac_m(nx_s+1)
    do ix = nx_s+1,nx_e-1
      Lap_Ac(ix) = clap0*Ac_m(ix) + clap1*(Ac_m(ix+1)+Ac_m(ix-1))
    end do
    Lap_Ac(nx_e) = clap0*Ac_m(nx_e) + clap1*Ac_m(nx_e-1)

    Ac_m_n = 2d0*Ac_m - Ac_m_o + (dt*clight)**2*Lap_Ac
    Ac_m_n(1:Mx) = Ac_m_n(1:Mx) -4d0*pi*dt**2*jt_m(1:Mx)
  end subroutine dt_evolve_macro_field_explicit
!--------------------------------------------------------------
  subroutine dt_evolve_macro_field_implicit
    implicit none
    real(8),parameter :: eps_error = 1d-28
    real(8) :: xvec(nx_s:nx_e), bvec(nx_s:nx_e)
    real(8) :: rvec(nx_s:nx_e), pvec(nx_s:nx_e)
    real(8) :: Axvec(nx_s:nx_e), Apvec(nx_s:nx_e)
    real(8) :: rr_p, rr_p_old, alpha, beta
    real(8) :: pAp_p
    real(8) :: clap0, clap1
    real(8) :: coeff0, coeff1
    real(8) :: bnorm, xnorm
    integer :: iter

    clap0 = -2d0/dx_m**2
    clap1 = 1d0/dx_m**2

    coeff0 = 1d0 -0.5d0*(dt*clight)**2*clap0
    coeff1 = -0.5d0*(dt*clight)**2*clap1

    call operate(Ac_m_o,bvec, coeff0, coeff1)
    bvec = -bvec + 2d0*Ac_m 
    bvec(1:Mx) = bvec(1:Mx) - 4d0*pi*dt**2*jt_m(1:Mx)
    bnorm = sum(bvec**2)

! initial guess
    call dt_evolve_macro_field_explicit
    xvec = Ac_m_n
    xnorm = sum(xvec**2)

    call operate(xvec,Axvec, coeff0, coeff1)
    rvec = bvec-Axvec
    pvec = rvec

    rr_p = sum(rvec**2)
    rr_p_old = rr_p
    iter = 0
    if(sqrt(rr_p/xnorm) < eps_error)then
      return
    end if

    do iter = 1, 10
      call operate(pvec,Apvec, coeff0, coeff1)
      pAp_p = sum(pvec*Apvec)

      alpha = rr_p/pAp_p
      xvec = xvec + alpha*pvec
      xnorm = sum(xvec**2)
      rvec = rvec - alpha*Apvec
      rr_p_old = rr_p
      rr_p = sum(rvec**2)
      if(sqrt(rr_p/xnorm) < eps_error)exit

      beta = rr_p/rr_p_old
      pvec = rvec + beta*pvec

    end do
    Ac_m_n = xvec

  end subroutine dt_evolve_macro_field_implicit
!--------------------------------------------------------------
  subroutine operate(xvec, Axvec, coeff0, coeff1)
    implicit none
    real(8),intent(in)  :: xvec(nx_s:nx_e), coeff0, coeff1
    real(8),intent(out) :: Axvec(nx_s:nx_e)
    integer :: ix

    Axvec(nx_s) = coeff0*xvec(nx_s) + coeff1*xvec(nx_s+1)
    do ix = nx_s+1,nx_e-1
      Axvec(ix) = coeff0*xvec(ix) + coeff1*(xvec(ix+1)+xvec(ix-1))
    end do
    Axvec(nx_e) = coeff0*xvec(nx_e) + coeff1*xvec(nx_e-1)


  end subroutine operate


end subroutine dt_evolve_macro_field

