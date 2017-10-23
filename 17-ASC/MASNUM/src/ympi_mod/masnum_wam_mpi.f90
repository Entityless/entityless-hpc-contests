!###############################################################################
!-------------------------------------------------------------------------------

  program masnum_wam_mpi

!-------------------------------------------------------------------------------

  use wammpi_mod
  use ympi_mod !, only: myid
 
  implicit none

  integer wav_count0, wav_count1, wav_count_rate, wav_count_max
  real :: elapse_time

!-------------------------------------------------------------------------------

  call wammpi_init

! Modified by Zhenya.Song, for timer, 2016/04/03
  if (myid == 0) then
    call system_clock(wav_count0, wav_count_rate, wav_count_max)
    if (wav_count_rate ==0) then
     write(6, '(a33)') '--- No system clock available ---'
     stop
    endif
  endif
! End Modify


  !Modified by Huang Chenghuan, 20170422.
  call init_main(myid) ! for the implemention of C code.

  call readwi_mpi

! Modified by Zhenya.Song, for timer, 2016/04/03
  if (myid == 0) then
    call system_clock(wav_count1, wav_count_rate, wav_count_max)
  endif
  if (myid == 0) then
    if (wav_count1 > wav_count0) then
     elapse_time = (wav_count1 - wav_count0 ) * 1.0 / wav_count_rate
    else
     elapse_time = (wav_count_max + wav_count1 - wav_count0) * 1.0 / wav_count_rate
    endif
    open(unit=11, file="elapse_time.txt", access = "append", status = "old")
    write(11, *) "Elapsed Time is ", elapse_time, " seconds."
    close(11)
    !print *, "Elapsed Time is ", elapse_time, " seconds."  !Modified by Huang Chenghuan, 20170422.
  endif
! End Modify
  
  call ympi_final

!-------------------------------------------------------------------------------

  end program masnum_wam_mpi

!-------------------------------------------------------------------------------
!###############################################################################
