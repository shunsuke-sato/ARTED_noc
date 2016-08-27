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
subroutine PSE_band_calculation
  use global_variables
  use PSE_variables
  use PSE_band_calc_variables
  implicit none
  integer :: isymp,id,ik,ib,iter_scf
  real(8) :: kk,dk,vec(3,2)

  if(myrank == 0) write(*,*)'Start band calculation'

  if(myrank == 0)then
    open(55,file = 'input_band_map.dat')
    read(55,*)Nsym_point, NB_band_calc
  end if
  
  call MPI_BCAST(Nsym_point,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(NB_band_calc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(b0_band_vec(Nsym_point),band_Cvec_d(3,Nsym_point))
  allocate(Name_sym_point(Nsym_point),Ndist_band(2:Nsym_point))

  if(myrank == 0)then
    read(55,*) b0_band_vec(1),band_Cvec_d(1,1),band_Cvec_d(2,1),band_Cvec_d(3,1),Name_sym_point(1)
    do isymp = 2,Nsym_point
      read(55,*) b0_band_vec(isymp),band_Cvec_d(1,isymp),band_Cvec_d(2,isymp),band_Cvec_d(3,isymp) &
        &,Ndist_band(isymp),Name_sym_point(isymp)
    end do
    close(55)
  end if

  call MPI_BCAST(b0_band_vec,Nsym_point,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(band_Cvec_d,3*Nsym_point,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Name_sym_point,10*Nsym_point,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Ndist_band,Nsym_point-1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  NK = Nsym_point + sum(Ndist_band(:))
  NB = NB_band_calc
  call NK_split
  

  deallocate(kAc_Rvec,kAc0_Rvec,kAc_Cvec,kAc0_Cvec)
  allocate(kAc_Rvec(3,NK),kAc0_Rvec(3,NK))
  allocate(kAc_Cvec(3,NK),kAc0_Cvec(3,NK))

  deallocate(ekr)
  allocate(ekr(Nps,NI,NK_s:NK_e)) ! sato

  deallocate(zu,zu_GS,zu_GS0)
  allocate(zu(NL,NB_TD,NK_s:NK_e),zu_GS(NL,NB,NK_s:NK_e),zu_GS0(NL,NB,NK_s:NK_e))


  deallocate(esp,esp_l)
  allocate(esp(NB,NK),esp_l(NB,NK))

  deallocate(esp,esp_l)
  allocate(esp(NB,NK),esp_l(NB,NK))

  deallocate(zutmp_diag,za_diag)
  allocate(zutmp_diag(NL,NB),za_diag(NB,NB))
  deallocate(esp_var)
  allocate(esp_var(NB,NK))


  kAc0_Cvec(1,1)=b0_band_vec(1)*band_Cvec_d(1,1)
  kAc0_Cvec(2,1)=b0_band_vec(1)*band_Cvec_d(2,1)
  kAc0_Cvec(3,1)=b0_band_vec(1)*band_Cvec_d(3,1)
  ik=1
  do isymp = 2,Nsym_point
    
    do id=1,Ndist_band(isymp)
      ik = ik + 1
      kAc0_Cvec(1,ik) = b0_band_vec(isymp-1)*band_Cvec_d(1,isymp-1)*dble(Ndist_band(isymp)+1-id)/dble(Ndist_band(isymp)+1) &
        & +b0_band_vec(isymp)*band_Cvec_d(1,isymp)*dble(id)/dble(Ndist_band(isymp)+1) 
        
      kAc0_Cvec(2,ik) = b0_band_vec(isymp-1)*band_Cvec_d(2,isymp-1)*dble(Ndist_band(isymp)+1-id)/dble(Ndist_band(isymp)+1) &
        & +b0_band_vec(isymp)*band_Cvec_d(2,isymp)*dble(id)/dble(Ndist_band(isymp)+1) 
        
      kAc0_Cvec(3,ik) = b0_band_vec(isymp-1)*band_Cvec_d(3,isymp-1)*dble(Ndist_band(isymp)+1-id)/dble(Ndist_band(isymp)+1) &
        & +b0_band_vec(isymp)*band_Cvec_d(3,isymp)*dble(id)/dble(Ndist_band(isymp)+1) 
    end do

    ik = ik + 1
    kAc0_Cvec(1,ik)=b0_band_vec(isymp)*band_Cvec_d(1,isymp)
    kAc0_Cvec(2,ik)=b0_band_vec(isymp)*band_Cvec_d(2,isymp)
    kAc0_Cvec(3,ik)=b0_band_vec(isymp)*band_Cvec_d(3,isymp)
    
  end do

  kAc_Cvec(:,:)=kAc0_Cvec(:,:)


  call init_wf

  if(Myrank == 0)write(*,'(A)')'Eigan states calculation start'
  do iter_scf=1,100
    if(Myrank == 0)write(*,'(A,2X,I5)')'iter_scf = ',iter_scf
    call Gram_Schmidt
    call PSE_Conjugate_Gradient(20)
    call Gram_Schmidt
    call PSE_subspace_diag


    if(Myrank == 0)then
      write(*,'(4(I3,e16.6e3))')(ib,esp(ib,NK_s),ib=1,NB)
    end if
  end do
  if(Myrank == 0)write(*,'(A)')'End Eigan states calculation start'

  if(Myrank == 0)write(*,'(A)')'Write band map'
  if(myrank == 0)then
    open(55,file = trim(SYSname)//'_band_map.dat')
    kk = 0d0
    ik = 1
    write(55,'(A,2x,A)')'#',Name_sym_point(1)
    write(55,'(999e26.16e3)')kk,kAc_Cvec(:,1),esp(:,1)

    do isymp = 2,Nsym_point
      vec(:,1) = b0_band_vec(isymp-1)*band_Cvec_d(:,isymp-1)
      vec(:,2) = b0_band_vec(isymp-1)*band_Cvec_d(:,isymp)
      dk = sqrt(sum((vec(:,1)-vec(:,2))**2))/dble(Ndist_band(isymp)+1)
      if(Ndist_band(isymp) == 0)dk = 0d0
      do id=1,Ndist_band(isymp)
        ik = ik + 1
        kk = kk + dk
        write(55,'(999e26.16e3)')kk,kAc_Cvec(:,ik),esp(:,ik)
      end do
      ik = ik + 1
      kk = kk + dk
      write(55,'(A,2x,A)')'#',Name_sym_point(isymp)
      write(55,'(999e26.16e3)')kk,kAc_Cvec(:,ik),esp(:,ik)

    end do

    close(55)
  end if
    

  if(myrank == 0) write(*,*)'Complete band calculation'

  return
end subroutine PSE_band_calculation
