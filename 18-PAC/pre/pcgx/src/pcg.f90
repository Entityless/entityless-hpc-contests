 program pcg_sample

   use kinds_mod
   use netcdf_mod
   use blocks
   use distribution
   use domain
   use domain_size
   use constants
   use boundary
   use global_reductions
   use gather_scatter
   use broadcast
   use grid
   use io
   use time_management
   use exit_mod
   use communicate, only: my_task, master_task
   use initial, only: initialize_pop
   use solver_pcg_mod
   use timers, only: get_timer, timer_start, timer_stop,timer_print_all

   implicit none
   character(len=9) :: outname
   integer (int_kind) :: ncid,n
   integer (int_kind), save :: pcg_timer

   call initialize_pop

   lprecond = .false.
   solv_max_iters = 1000
   solv_ncheck    = 10
   solv_convrg    = 3121084.94389626
!-----------------------------------------------------------------------
!
!  compute initial residual and initialize S
!
!-----------------------------------------------------------------------
   write(outname,'(i9)')my_task
   if (my_task<10) outname = '00'//trim(adjustl(outname//'.nc'))
   if (my_task>=10 .and. my_task<100) outname = '0'//trim(adjustl(outname//'.nc'))
   if (my_task>=100 .and. my_task<1000) outname = trim(adjustl(outname//'.nc'))

   call open_nc(NcID, './data_in/'//outname, "read")
   call readnc(NcID, 'A0', A0)
   call readnc(NcID, 'AN', AN)
   call readnc(NcID, 'AE', AE)
   call readnc(NcID, 'ANE', ANE)
   call readnc(NcID, 'RCALCT_B', RCALCT_B)
   call readnc(NcID, 'X', X0)
   call readnc(NcID, 'B', B)
   call close_nc(NcID)

   call get_timer(pcg_timer,'pcg',1,distrb_clinic%nprocs)
   call timer_start(pcg_timer)

   do n =1,50
      X = X0
      call pcg(X,B)
   enddo

   call timer_stop(pcg_timer)
!-----------------------------------------------------------------------
!  Result validation with original pcg
!-----------------------------------------------------------------------
   solv_ncheck    = 1
   if (my_task .eq. 1) then
      write(6,*) 'Being verified...'
   endif
   call pcg(X,B)

   call timer_print_all(stats=.true.)

   call exit_POP(sigExit,'Successful completion of pcg run')

 end program pcg_sample
