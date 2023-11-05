subroutine read_occupation_from_file
  use global_variables
  use PSE_variables
  implicit none
  integer :: ib, ik
  real(8) :: occ_tmp(NB)

  occ_tmp = 0d0
  if(myrank == 0)then
    open(101,file='occs')
    do ib = 1, NB
      read(101,*)occ_tmp(ib)
    end do
    close(101)
  end if

  call MPI_BCAST(occ_tmp,NB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  do ik = NK_s, NK_e
    occ(1:NB,ik) = occ_tmp(1:NB)/dble(NK)
  end do

end subroutine read_occupation_from_file
