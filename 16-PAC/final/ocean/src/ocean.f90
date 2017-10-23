program ocean
    use mod_data
    use mpdata_adiff
#ifdef USE_MPI
    use mpi
#endif
    implicit none
    integer start_time, stop_time, clock_rate, clock_max
    integer ierr, myrank
    real*8 t1

#ifndef USE_MPI
! General procedure.
! *************************************************************
    write (*,*) "[INIT] Initializing data... NO MPI"
    call allocate_init_data()

    write (*,*) "[INIT] Done!"
    write (*,*) 

    write (*,*) "[RUN]  Start calculation..."
    call system_clock(start_time,clock_rate,clock_max)

    call mpdata_adiff_tile(LBi, UBi, LBj, UBj, MYDATA%oHz, MYDATA%Huon, MYDATA%Hvom, MYDATA%W, MYDATA%Ta, &
                         & MYDATA%Uind, MYDATA%Dn, MYDATA%Dm, MYDATA%Ua)

    call system_clock(stop_time,clock_rate,clock_max)
    t1 = real(stop_time - start_time)/real(clock_rate)

    write(*,'(a,f11.3,a)') '        Calculation time:',t1,'(s)'

    write (*,*) "[RUN]  Done!"
    write (*,*) 

    write (*,*) "[VERIFY] Start verification."
    call verify_data()
    write (*,*) "[VERIFY] Done!"
    call deallocate_data()

#else
! MPI version framwork
!**************************************************************
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if (myrank.eq.0) write (*,*) "[INIT] Initializing data... USE MPI"
    if (myrank.eq.0) call allocate_init_data()

! Should be totally implemented for multi-process version.
! Would not counted in to the calculation time.
    call distribute_init_data() ! Not implemented.

! A dummy call only used to make the cache dirty.
    call flush_cache(t1)

    if (myrank.eq.0) write (*,*) "[INIT] Done!"
    if (myrank.eq.0) write (*,*) 

    if (myrank.eq.0) write (*,*) "[RUN]  Start calculation..."
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call system_clock(start_time,clock_rate,clock_max)

! If programmer needs to use distributed version, just use this style function call.
    call mpdata_adiff_tile_mpi() ! Not implemented.

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call system_clock(stop_time,clock_rate,clock_max)
    t1 = real(stop_time - start_time)/real(clock_rate)

    if (myrank.eq.0) write(*,'(a,f8.3,a)') '        Calculation time:',t1,'(s)'

    if (myrank.eq.0) write (*,*) "[RUN]  Done!"
    if (myrank.eq.0) write (*,*) 

    if (myrank.eq.0) write (*,*) "[VERIFY] Start verification."
! Needs to be implemented to gather the distributed results to MYDATA%Ua on proc0.
! Would not counted in to the calculation time.
    call gather_data()

    if (myrank.eq.0) call verify_data()
    if (myrank.eq.0) write (*,*) "[VERIFY] Done!"

    if (myrank.eq.0) call deallocate_data()

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)

#endif
end program
