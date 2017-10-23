subroutine update_communication(e_comm, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(4), intent(inout) :: e_comm(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    integer :: ia, ic
    !integer :: req(4)
    !integer :: stat(mpi_status_size,4)
!    if(allocated(dnee1_3d)) then
!        deallocate(dnee1_3d, upee1_3d, lfee1_3d, rtee1_3d)
!    endif
!    allocate(dnee1_3d(kl, jl, ixl, halo), upee1_3d(kl, jl, ixl, halo), &
!        lfee1_3d(kl, jl, halo, iyl), rtee1_3d(kl, jl, halo, iyl))
!    if(allocated(dnee_3d)) then
!       deallocate(dnee_3d, upee_3d, lfee_3d, rtee_3d)
!    endif
!    allocate(dnee_3d(kl, jl, ixl, halo), upee_3d(kl, jl, ixl, halo),&
!        lfee_3d(kl, jl, halo, iyl), rtee_3d(kl, jl, halo, iyl))

    if(allocated(dnee1_3d)) then
        !deallocate(dnee1_3d, upee1_3d, lfee1_3d, rtee1_3d)
    else
        allocate(dnee1_3d(kl, jl, my_hor_len), upee1_3d(kl, jl, my_hor_len))
        allocate(lfee1_3d(kl, jl, my_ver_len), rtee1_3d(kl, jl, my_ver_len))

        dnee1_3d = small
        upee1_3d = small
        lfee1_3d = small
        rtee1_3d = small
        
        !allocate(lfee1_3d(kl, jl, iyl), rtee1_3d(kl, jl, iyl))
        allocate(lsreq(mypebox%nl_send))
        allocate(lrreq(mypebox%nl_recv))
        allocate(rsreq(mypebox%nr_send))
        allocate(rrreq(mypebox%nr_recv))
    endif

    if(allocated(dnee_3d)) then
        !deallocate(dnee1_3d, upee1_3d, lfee1_3d, rtee1_3d)
    else
        allocate(dnee_3d(kl, jl, my_hor_len), upee_3d(kl, jl, my_hor_len), &
            lfee_3d(kl, jl, iyl), rtee_3d(kl, jl, iyl))
        !not my_ver_len here!!!
    endif

    !actually, only halosize = 1 will enter this function

    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif

!    do i = 1, mypebox%nl_send 
!      ncount = kl*jl*(mypebox%left_send(i, 3)-mypebox%left_send(i, 2)+1)
!      if(ncount < 0 .or. mypebox%left_send(i, 1) < 0 .or. mypebox%left_send(i, 1) > 63)then
!        print * , "isend ****, myid = ", myid, i, "left_send: ", mypebox%left_send(i, :)
!      endif
!    end do

    
    ! --- From down to up.
    !call mpi_sendrecv(upee1_3d, ncount, mpi_real, mypebox%up, 100, dnee_3d , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! ! --- From up to down.
    !call mpi_sendrecv(dnee1_3d, ncount, mpi_real, mypebox%dn, 200, upee_3d , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    ! --- Send to up.

    ic = 2
    do ia = my_hor_start, my_hor_end 
      if(nsp(ia, ic) == 0.0) cycle
      dnee1_3d(:, :, ia - my_hor_start + 1) = e_comm(:, :, ia, ic)
    end do 

    ic = iyl - 1
    do ia = my_hor_start, my_hor_end 
      if(nsp(ia, ic) == 0.0) cycle
      upee1_3d(:, :, ia - my_hor_start + 1) = e_comm(:, :, ia, ic)
    end do 

    !dnee1_3d(:, :, :) = e_comm(:, :, my_hor_start : my_hor_end, 2)
    !upee1_3d(:, :, :) = e_comm(:, :, my_hor_start : my_hor_end, iyl - 1)

    lfee1_3d(:, :, :) = e_comm(:, :,  2, my_ver_start : my_ver_end)
    rtee1_3d(:, :, :) = e_comm(:, :,  ixl-1, my_ver_start : my_ver_end)

    ia = 2
    do ic = my_ver_start, my_ver_end 
      if(nsp(ia, ic) == 0.0) cycle
      lfee1_3d(:, :, ic - my_ver_start + 1) = e_comm(:, :,  ia, ic)
    end do 
    
    ia = ixl-1
    do ic = my_ver_start, my_ver_end 
      if(nsp(ia, ic) == 0.0) cycle
      rtee1_3d(:, :, ic - my_ver_start + 1) = e_comm(:, :,  ia, ic)
    end do 


!    ncount = kl*jl*my_hor_len
!    if(mypebox%up /= MPI_PROC_NULL) call mpi_isend(upee1_3d, ncount, mpi_real, mypebox%up, 100, mpi_comm_ympi, updn_req(1), ierr)
!   !  !--- Send to down.
!    if(mypebox%dn /= MPI_PROC_NULL) call mpi_isend(dnee1_3d, ncount, mpi_real, mypebox%dn, 200, mpi_comm_ympi, updn_req(2), ierr)    
!   
!    ncount = kl*jl*my_hor_len
!       !  !--- Receive from dnee.
!    if(mypebox%dn /= MPI_PROC_NULL) call mpi_irecv(dnee_3d,  ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, updn_req(3), ierr)
!   !  !--- Receive from up.
!    if(mypebox%up /= MPI_PROC_NULL) call mpi_irecv(upee_3d,  ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, updn_req(4), ierr)

    ncount = kl*jl*my_hor_len
    if(mod(mypebox%ver_id, 2) == 0)then
      if(mypebox%up /= MPI_PROC_NULL) call mpi_isend(upee1_3d, ncount, mpi_real, mypebox%up, 100, mpi_comm_ympi, updn_req(1), ierr)
      if(mypebox%dn /= MPI_PROC_NULL) call mpi_irecv(dnee_3d,  ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, updn_req(2), ierr)
      if(mypebox%dn /= MPI_PROC_NULL) call mpi_isend(dnee1_3d, ncount, mpi_real, mypebox%dn, 200, mpi_comm_ympi, updn_req(3), ierr)
      if(mypebox%up /= MPI_PROC_NULL) call mpi_irecv(upee_3d,  ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, updn_req(4), ierr)
    else
      if(mypebox%dn /= MPI_PROC_NULL) call mpi_irecv(dnee_3d,  ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, updn_req(2), ierr)
      if(mypebox%up /= MPI_PROC_NULL) call mpi_isend(upee1_3d, ncount, mpi_real, mypebox%up, 100, mpi_comm_ympi, updn_req(1), ierr)
      if(mypebox%up /= MPI_PROC_NULL) call mpi_irecv(upee_3d,  ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, updn_req(4), ierr)
      if(mypebox%dn /= MPI_PROC_NULL) call mpi_isend(dnee1_3d, ncount, mpi_real, mypebox%dn, 200, mpi_comm_ympi, updn_req(3), ierr)
    endif

    if(mod(mypebox%hor_id, 2) == 0)then
      do i = 1, mypebox%nr_send 
        ncount = kl*jl*(mypebox%right_send(i, 3)-mypebox%right_send(i, 2)+1)
        call mpi_isend(rtee1_3d(1, 1, mypebox%right_send(i, 2)),ncount, mpi_real, mypebox%right_send(i, 1), 101,mpi_comm_ympi, rsreq(i), ierr                         )
      end do
      do i = 1, mypebox%nl_recv 
        ncount = kl*jl*(mypebox%left_recv(i, 3)-mypebox%left_recv(i, 2)+1)
        call mpi_irecv(lfee_3d(1, 1, mypebox%left_recv(i, 2)),ncount, mpi_real, mypebox%left_recv(i, 1), 101,mpi_comm_ympi, lrreq(i), ierr                         )
      end do
      do i = 1, mypebox%nl_send 
        ncount = kl*jl*(mypebox%left_send(i, 3)-mypebox%left_send(i, 2)+1)
        call mpi_isend(lfee1_3d(1, 1, mypebox%left_send(i, 2)),ncount, mpi_real, mypebox%left_send(i, 1), 101,mpi_comm_ympi, lsreq(i), ierr                         )
      end do
      do i = 1, mypebox%nr_recv 
        ncount = kl*jl*(mypebox%right_recv(i, 3)-mypebox%right_recv(i, 2)+1)
        call mpi_irecv(rtee_3d(1, 1, mypebox%right_recv(i, 2)),ncount, mpi_real, mypebox%right_recv(i, 1), 101,mpi_comm_ympi, rrreq(i), ierr                         )
      end do
    else
      do i = 1, mypebox%nl_recv 
        ncount = kl*jl*(mypebox%left_recv(i, 3)-mypebox%left_recv(i, 2)+1)
        call mpi_irecv(lfee_3d(1, 1, mypebox%left_recv(i, 2)),ncount, mpi_real, mypebox%left_recv(i, 1), 101,mpi_comm_ympi, lrreq(i), ierr                         )
      end do
      do i = 1, mypebox%nr_send 
        ncount = kl*jl*(mypebox%right_send(i, 3)-mypebox%right_send(i, 2)+1)
        call mpi_isend(rtee1_3d(1, 1, mypebox%right_send(i, 2)),ncount, mpi_real, mypebox%right_send(i, 1), 101,mpi_comm_ympi, rsreq(i), ierr                         )
      end do
      do i = 1, mypebox%nr_recv 
        ncount = kl*jl*(mypebox%right_recv(i, 3)-mypebox%right_recv(i, 2)+1)
        call mpi_irecv(rtee_3d(1, 1, mypebox%right_recv(i, 2)),ncount, mpi_real, mypebox%right_recv(i, 1), 101,mpi_comm_ympi, rrreq(i), ierr                         )
      end do
      do i = 1, mypebox%nl_send 
        ncount = kl*jl*(mypebox%left_send(i, 3)-mypebox%left_send(i, 2)+1)
        call mpi_isend(lfee1_3d(1, 1, mypebox%left_send(i, 2)),ncount, mpi_real, mypebox%left_send(i, 1), 101,mpi_comm_ympi, lsreq(i), ierr                         )
      end do
    endif

!    do i = 1, mypebox%nl_send 
!      ncount = kl*jl*(mypebox%left_send(i, 3)-mypebox%left_send(i, 2)+1)
!      call mpi_isend(lfee1_3d(1, 1, mypebox%left_send(i, 2)),ncount, mpi_real, mypebox%left_send(i, 1), 101,mpi_comm_ympi, lsreq(i), ierr                         )
!    end do
!
!
!    do i = 1, mypebox%nl_recv 
!      ncount = kl*jl*(mypebox%left_recv(i, 3)-mypebox%left_recv(i, 2)+1)
!      call mpi_irecv(lfee_3d(1, 1, mypebox%left_recv(i, 2)),ncount, mpi_real, mypebox%left_recv(i, 1), 101,mpi_comm_ympi, lrreq(i), ierr                         )
!    end do
!
!    do i = 1, mypebox%nr_send 
!      ncount = kl*jl*(mypebox%right_send(i, 3)-mypebox%right_send(i, 2)+1)
!      call mpi_isend(rtee1_3d(1, 1, mypebox%right_send(i, 2)),ncount, mpi_real, mypebox%right_send(i, 1), 101,mpi_comm_ympi, rsreq(i), ierr                         )
!    end do
!
!    do i = 1, mypebox%nr_recv 
!      ncount = kl*jl*(mypebox%right_recv(i, 3)-mypebox%right_recv(i, 2)+1)
!      call mpi_irecv(rtee_3d(1, 1, mypebox%right_recv(i, 2)),ncount, mpi_real, mypebox%right_recv(i, 1), 101,mpi_comm_ympi, rrreq(i), ierr                         )
!    end do

end subroutine update_communication

subroutine update_wait(e_comm, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(4), intent(inout) :: e_comm(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
            ! --- Send to up.
    if(mypebox%up /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(1),stat,ierr)
      !--- Send to down.
    if(mypebox%dn /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(2),stat,ierr)
      !--- Receive from dnee.
    if(mypebox%dn /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(3),stat,ierr)
      !--- Receive from up.
    if(mypebox%up /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(4),stat,ierr)
    
    !if(mypebox%dn /= MPI_PROC_NULL)e_comm(:, :, my_hor_start : my_hor_end, 1) = dnee_3d
    !if(mypebox%up /= MPI_PROC_NULL)e_comm(:, :, my_hor_start : my_hor_end, iyl) = upee_3d


    if(mypebox%nl_send > 0)then
      call mpi_waitall(mypebox%nl_send, lsreq, stat_l_send(:, 1:mypebox%nl_send), ierr)
    endif

    if(mypebox%nl_recv > 0)then
      call mpi_waitall(mypebox%nl_recv, lrreq, stat_l_recv(:, 1:mypebox%nl_recv), ierr)
      !e_comm(:, :, 1, :) = lfee_3d
    endif

    if(mypebox%nr_send > 0)then
      call mpi_waitall(mypebox%nr_send, rsreq, stat_r_send(:, 1:mypebox%nr_send), ierr)
    endif

    if(mypebox%nr_recv > 0)then
      call mpi_waitall(mypebox%nr_recv, rrreq, stat_r_recv(:, 1:mypebox%nr_recv), ierr)
      !e_comm(:, :, ixl, :) = rtee_3d
    endif


!    if(mypebox%nr > 0)then
!      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
!      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
!      e_comm(:, :, ixl, :) = rtee_3d
!    endif
!    if(mypebox%nl > 0)then
!      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
!      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
!      e_comm(:, :, 1, :) = lfee_3d
!    endif

end subroutine update_wait

subroutine update_copy_back(e_comm, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(4), intent(inout) :: e_comm(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif

!    if(mypebox%up /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(1),stat,ierr)
!      !--- Send to down.
!    if(mypebox%dn /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(2),stat,ierr)
!      !--- Receive from dnee.
!    if(mypebox%dn /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(3),stat,ierr)
!      !--- Receive from up.
!    if(mypebox%up /= MPI_PROC_NULL) call mpi_waitall(1, updn_req(4),stat,ierr)
    
    if(mypebox%dn /= MPI_PROC_NULL)e_comm(:, :, my_hor_start : my_hor_end, 1) = dnee_3d
    if(mypebox%up /= MPI_PROC_NULL)e_comm(:, :, my_hor_start : my_hor_end, iyl) = upee_3d


!    if(mypebox%nl_send > 0)then
!      call mpi_waitall(mypebox%nl_send, lsreq, stat_l_send(:, 1:mypebox%nl_send), ierr)
!    endif

    if(mypebox%nl_recv > 0)then
      !call mpi_waitall(mypebox%nl_recv, lrreq, stat_l_recv(:, 1:mypebox%nl_recv), ierr)
      e_comm(:, :, 1, :) = lfee_3d
    endif

!    if(mypebox%nr_send > 0)then
!      call mpi_waitall(mypebox%nr_send, rsreq, stat_r_send(:, 1:mypebox%nr_send), ierr)
!    endif

    if(mypebox%nr_recv > 0)then
      !call mpi_waitall(mypebox%nr_recv, rrreq, stat_r_recv(:, 1:mypebox%nr_recv), ierr)
      e_comm(:, :, ixl, :) = rtee_3d
    endif


  end subroutine update_copy_back

!-------------------------------------------------------------------------------
!###############################################################################
