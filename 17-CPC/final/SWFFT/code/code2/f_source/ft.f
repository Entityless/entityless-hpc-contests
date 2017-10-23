
      program main_fft

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none

      include 'mpif.h'
      integer i, ierr
      
c---------------------------------------------------------------------
c u0, u1, u2 are the main arrays in the problem. 
c Depending on the decomposition, these arrays will have different 
c dimensions. To accomodate all possibilities, we allocate them as 
c one-dimensional arrays and pass them to subroutines for different 
c views
c  - u0 contains the initial (transformed) initial condition
c  - u1 and u2 are working arrays
c---------------------------------------------------------------------

      double complex, allocatable ::  u0(:),u1(:),u2(:)
      double precision,allocatable :: twiddle(:)

      double complex pad1(3), pad2(3), pad3(3)
c      integer :: 
      double complex :: dbg_dc_val
      double precision :: dbg_d_val

      integer iter
      double precision total_time, mflops
      logical verified
      character class
      call MPI_Init(ierr)

c---------------------------------------------------------------------
c Run the entire problem once to make sure all data is touched. 
c This reduces variable startup costs, which is important for such a 
c short benchmark. 
c---------------------------------------------------------------------
      do i = 1, t_max
         call timer_clear(i)
      end do

      call timer_start(T_init)
      call setup()
      ntdivnp=((nx*ny)/np_min)*nz
      ntotal_f=1.d0*nx*ny*nz
      allocate(u(nx),u0(ntdivnp),u1(ntdivnp),u2(ntdivnp))
      allocate(sums(0:niter_default))  
      allocate(twiddle(ntdivnp))
      call compute_indexmap(twiddle, dims(1,3), dims(2,3), dims(3,3))
      call compute_initial_conditions(u1, dims(1,1), dims(2,1), 
     >                                dims(3,1))
      call fft_init (dims(1,1))
      call fft(1, u1, u0)
      call timer_stop(T_init)
      if (me .eq. 0) then
         print *,'Initialization time =',timer_read(T_init)
      endif

c---------------------------------------------------------------------
c Start over from the beginning. Note that all operations must
c be timed, in contrast to other benchmarks. 
c---------------------------------------------------------------------
      
      do i = 1, t_max
         call timer_clear(i)
      end do
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call timer_start(t_total)

      call compute_indexmap(twiddle, dims(1,3), dims(2,3), dims(3,3))
      call compute_initial_conditions(u1, dims(1,1), dims(2,1), 
     >                                dims(3,1))
      call fft_init (dims(1,1))

      call fft(1, u1, u0)

      if(me .eq. 0) then
         print *, "<check> fortran"
         print *, niter, fftblockpad_default, maxdim
         print *, fftblock, layout_type, me
         print *, np1, np2, np
         print *, ntdivnp
         print *, transblockpad, transblock, fftblockpad
         print *, nx, ny, nz
         print *, "</check> fortran"
      endif
      call c_param_init(u0, u1, u2, twiddle, dims, niter,
     >         fftblockpad_default, fftblock, maxdim, layout_type,
     >         u, me, np1, np2, ntdivnp,
     >         np, transblockpad,
     >         transblock, fftblockpad,
     >         nx, ny, nz,
     >         xstart, ystart, zstart,
     >         xend, yend, zend, sums)

      if(me .eq. 0) then
         dbg_dc_val = (2.5, 22.1)
         dbg_d_val = -7.9
         print *, "<fortran dbg: "
         print *, dbg_d_val, real(dbg_dc_val), imag(dbg_dc_val)
         print *, dconjg(dbg_dc_val)
         print *, dbg_dc_val * dbg_d_val
         print *, dbg_dc_val * (-2.0, -2.0)
         print *, dbg_dc_val + (-2.0, -2.0)
         print *, dbg_dc_val - (-2.0, -2.0)
         print *, "</fortran dbg: "
         call fortran_c_dbg(dbg_dc_val, dbg_d_val)
      endif

      call c_fft_iter()

c      do iter = 1, niter
c         call evolve(u0, u1, twiddle, dims(1,1), dims(2,1), dims(3,1))
c!         if(mod(iter,100) .eq. 0) then
c         call fft(-1, u1, u2)
c!         if (timers_enabled) call synchup()
c            call checksum(iter, u2, dims(1,1), dims(2,1), dims(3,1))
c!         endif
c      end do

      call c_finalize()

      call timer_stop(t_total)
      total_time = timer_read(t_total)

      if( total_time .ne. 0. ) then
         mflops = 1.0d-6*ntotal_f *
     >             (14.8157+7.19641*log(ntotal_f)
     >          +  (5.23518+7.21113*log(ntotal_f))*niter)
     >                 /total_time
      else
         mflops = 0.0
      endif
      if (me .eq. 0) then
         call print_results('FT', class, nx, ny, nz, niter, np_min, np,
     >     total_time, mflops, '          floating point', verified, 
     >     npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)
      endif
C       if (me .eq. 0) print *, "total time = ", total_time
      if (timers_enabled) call print_timers()
      deallocate(u,u0,u1,u2)
      deallocate(sums)
      deallocate(twiddle)   
      call MPI_Finalize(ierr)
      end

      subroutine mmpmmpmmp(i, allchk)
      use global
      implicit none
      integer :: i
      double complex allchk

      sums(i) = allchk

      print *, sums(i)

      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine evolve(u0, u1, twiddle, d1, d2, d3)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c evolve u0 -> u1 (t time steps) in fourier space
c---------------------------------------------------------------------
      use global
      implicit none
c      include 'global.h'
      integer d1, d2, d3
      double precision exi
      double complex u0(d1,d2,d3)
      double complex u1(d1,d2,d3)
      double precision twiddle(d1,d2,d3)
      integer i, j, k

      do k = 1, d3
         do j = 1, d2
            do i = 1, d1
               u0(i,j,k) = u0(i,j,k)*(twiddle(i,j,k))
               u1(i,j,k) = u0(i,j,k)
            end do
         end do
      end do

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine compute_initial_conditions(u0, d1, d2, d3)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c Fill in array u0 with initial conditions from 
c random number generator 
c---------------------------------------------------------------------
      use global
      implicit none
c      include 'global.h'
      integer d1, d2, d3
      double complex u0(d1, d2, d3)
      integer k
      double precision x0, start, an, dummy
      
c---------------------------------------------------------------------
c 0-D and 1-D layouts are easy because each processor gets a contiguous
c chunk of the array, in the Fortran ordering sense. 
c For a 2-D layout, it's a bit more complicated. We always
c have entire x-lines (contiguous) in processor. 
c We can do ny/np1 of them at a time since we have
c ny/np1 contiguous in y-direction. But then we jump
c by z-planes (nz/np2 of them, total). 
c For the 0-D and 1-D layouts we could do larger chunks, but
c this turns out to have no measurable impact on performance. 
c---------------------------------------------------------------------


      start = seed                                    
c---------------------------------------------------------------------
c Jump to the starting element for our first plane.
c---------------------------------------------------------------------
      call ipow46(a, 2*nx, (zstart(1)-1)*ny + (ystart(1)-1), an)
      dummy = randlc(start, an)
      call ipow46(a, 2*nx, ny, an)
      
c---------------------------------------------------------------------
c Go through by z planes filling in one square at a time.
c---------------------------------------------------------------------
      do k = 1, dims(3, 1) ! nz/np2
         x0 = start
         call vranlc(2*nx*dims(2, 1), x0, a, u0(1, 1, k))
         if (k .ne. dims(3, 1)) dummy = randlc(start, an)
      end do

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine ipow46(a, exp_1, exp_2, result)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c compute a^exponent mod 2^46
c---------------------------------------------------------------------

      implicit none
      double precision a, result, dummy, q, r
      integer exp_1, exp_2, n, n2, ierr
      external randlc
      double precision randlc
      logical  two_pow
c---------------------------------------------------------------------
c Use
c   a^n = a^(n/2)*a^(n/2) if n even else
c   a^n = a*a^(n-1)       if n odd
c---------------------------------------------------------------------
      result = 1
      if (exp_2 .eq. 0 .or. exp_1 .eq. 0) return
      q = a
      r = 1
      n = exp_1
      two_pow = .true.

      do while (two_pow)
         n2 = n/2
         if (n2 * 2 .eq. n) then
            dummy = randlc(q, q)
            n = n2
         else
            n = n * exp_2
            two_pow = .false.
         endif
      end do

      do while (n .gt. 1)
         n2 = n/2
         if (n2 * 2 .eq. n) then
            dummy = randlc(q, q) 
            n = n2
         else
            dummy = randlc(r, q)
            n = n-1
         endif
      end do
      dummy = randlc(r, q)
      result = r
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine setup

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
c      include 'global.h'

      integer ierr, i, j, fstatus
      debug = .FALSE.
      
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)

      if (.not. convertdouble) then
         dc_type = MPI_DOUBLE_COMPLEX
      else
         dc_type = MPI_COMPLEX
      endif


      if (me .eq. 0) then
         write(*, 1000)

         open (unit=2,file='timer.flag',status='old',iostat=fstatus)
         timers_enabled = .false.
         if (fstatus .eq. 0) then
            timers_enabled = .true.
            close(2)
         endif

         open (unit=2,file='inputft.data',status='old', iostat=fstatus)

         if (fstatus .eq. 0) then
            write(*,233) 
 233        format(' Reading from input file inputft.data')
            read (2,*) nx,ny,nz,maxdim
            read (2,*) niter_default
            read (2,*) np_min
!            read (2,*) niter
            read (2,*) layout_type
            read (2,*) np1, np2
            close(2)
            niter=niter_default
c---------------------------------------------------------------------
c check to make sure input data is consistent
c---------------------------------------------------------------------

    
c---------------------------------------------------------------------
c 1. product of processor grid dims must equal number of processors
c---------------------------------------------------------------------

            if (np1 * np2 .ne. np) then
               write(*, 238)
 238           format(' np1 and np2 given in input file are not valid.')
               write(*, 239) np1*np2, np
 239           format(' Product is ', i5, ' and should be ', i5)
               call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

c---------------------------------------------------------------------
c 2. layout type must be valid
c---------------------------------------------------------------------

            if (layout_type .ne. layout_0D .and.
     >          layout_type .ne. layout_1D .and.
     >          layout_type .ne. layout_2D) then
               write(*, 240)
 240           format(' Layout type specified in inputft.data is 
     >                  invalid ')
               call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

c---------------------------------------------------------------------
c 3. 0D layout must be 1x1 grid
c---------------------------------------------------------------------

            if (layout_type .eq. layout_0D .and.
     >            (np1 .ne.1 .or. np2 .ne. 1)) then
               write(*, 241)
 241           format(' For 0D layout, both np1 and np2 must be 1 ')
               call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
c---------------------------------------------------------------------
c 4. 1D layout must be 1xN grid
c---------------------------------------------------------------------

            if (layout_type .eq. layout_1D .and. np1 .ne. 1) then
               write(*, 242)
 242           format(' For 1D layout, np1 must be 1 ')
               call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

         else
            write(*,234) 
            if (np .eq. 1) then
               np1 = 1
               np2 = 1
               layout_type = layout_0D
            else if (np .le. nz) then
               np1 = 1
               np2 = np
               layout_type = layout_1D
            else
               np1 = nz
               np2 = np/nz
               layout_type = layout_2D
            endif
         endif

         if (np .lt. np_min) then
            write(*, 10) np_min
 10         format(' Error: Compiled for ', I5, ' processors. ')
            write(*, 11) np
 11         format(' Only ',  i5, ' processors found ')
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
         endif

 234     format(' No input file inputft.data. Using compiled defaults')
         write(*, 1001) nx, ny, nz
         write(*, 1002) niter
         write(*, 1004) np
         write(*, 1005) np1, np2
!         if (np .ne. np_min) write(*, 1006) np_min

         if (layout_type .eq. layout_0D) then
            write(*, 1010) '0D'
         else if (layout_type .eq. layout_1D) then
            write(*, 1010) '1D'
         else
            write(*, 1010) '2D'
         endif

 1000 format(//,' CPC -- FFT Benchmark',/)
 1001    format(' Size                : ', i4, 'x', i4, 'x', i4)
 1002    format(' Iterations          : ', 7x, i7)
 1004    format(' Number of processes : ', 7x, i7)
 1005    format(' Processor array     : ', 5x, i4, 'x', i4)
 1006    format(' WARNING: compiled for ', i5, ' processes. ',
     >          ' Will not verify. ')
 1010    format(' Layout type         : ', 9x, A5)
      endif


c---------------------------------------------------------------------
c Broadcast parameters 
c---------------------------------------------------------------------
      call MPI_BCAST(np1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(np2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(layout_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     &               ierr)
      call MPI_BCAST(niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(maxdim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(niter_default, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &              ierr)
      call MPI_BCAST(np_min, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    
c      call MPI_BCAST(timers_enabled, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
c     &               ierr)

      if (np1 .eq. 1 .and. np2 .eq. 1) then
        layout_type = layout_0D
      else if (np1 .eq. 1) then
         layout_type = layout_1D
      else
         layout_type = layout_2D
      endif

      if (layout_type .eq. layout_0D) then
         do i = 1, 3
            dims(1, i) = nx
            dims(2, i) = ny
            dims(3, i) = nz
         end do
      else if (layout_type .eq. layout_1D) then
         dims(1, 1) = nx
         dims(2, 1) = ny
         dims(3, 1) = nz

         dims(1, 2) = nx
         dims(2, 2) = ny
         dims(3, 2) = nz

         dims(1, 3) = nz
         dims(2, 3) = nx
         dims(3, 3) = ny
      else if (layout_type .eq. layout_2D) then
         dims(1, 1) = nx
         dims(2, 1) = ny
         dims(3, 1) = nz

         dims(1, 2) = ny
         dims(2, 2) = nx
         dims(3, 2) = nz

         dims(1, 3) = nz
         dims(2, 3) = nx
         dims(3, 3) = ny

      endif
      do i = 1, 3
         dims(2, i) = dims(2, i) / np1
         dims(3, i) = dims(3, i) / np2
      end do


c---------------------------------------------------------------------
c Determine processor coordinates of this processor
c Processor grid is np1xnp2. 
c Arrays are always (n1, n2/np1, n3/np2)
c Processor coords are zero-based. 
c---------------------------------------------------------------------
      me2 = mod(me, np2)  ! goes from 0...np2-1
      me1 = me/np2        ! goes from 0...np1-1
c---------------------------------------------------------------------
c Communicators for rows/columns of processor grid. 
c commslice1 is communicator of all procs with same me1, ranked as me2
c commslice2 is communicator of all procs with same me2, ranked as me1
c mpi_comm_split(comm, color, key, ...)
c---------------------------------------------------------------------
      call MPI_Comm_split(MPI_COMM_WORLD, me1, me2, commslice1, ierr)
      call MPI_Comm_split(MPI_COMM_WORLD, me2, me1, commslice2, ierr)
!      if (timers_enabled) call synchup()

      if (debug) print *, 'proc coords: ', me, me1, me2

c---------------------------------------------------------------------
c Determine which section of the grid is owned by this
c processor. 
c---------------------------------------------------------------------
      if (layout_type .eq. layout_0d) then

         do i = 1, 3
            xstart(i) = 1
            xend(i)   = nx
            ystart(i) = 1
            yend(i)   = ny
            zstart(i) = 1
            zend(i)   = nz
         end do

      else if (layout_type .eq. layout_1d) then

         xstart(1) = 1
         xend(1)   = nx
         ystart(1) = 1
         yend(1)   = ny
         zstart(1) = 1 + me2 * nz/np2
         zend(1)   = (me2+1) * nz/np2

         xstart(2) = 1
         xend(2)   = nx
         ystart(2) = 1
         yend(2)   = ny
         zstart(2) = 1 + me2 * nz/np2
         zend(2)   = (me2+1) * nz/np2

         xstart(3) = 1
         xend(3)   = nx
         ystart(3) = 1 + me2 * ny/np2
         yend(3)   = (me2+1) * ny/np2
         zstart(3) = 1
         zend(3)   = nz

      else if (layout_type .eq. layout_2d) then

         xstart(1) = 1
         xend(1)   = nx
         ystart(1) = 1 + me1 * ny/np1
         yend(1)   = (me1+1) * ny/np1
         zstart(1) = 1 + me2 * nz/np2
         zend(1)   = (me2+1) * nz/np2

         xstart(2) = 1 + me1 * nx/np1
         xend(2)   = (me1+1)*nx/np1
         ystart(2) = 1
         yend(2)   = ny
         zstart(2) = zstart(1)
         zend(2)   = zend(1)

         xstart(3) = xstart(2)
         xend(3)   = xend(2)
         ystart(3) = 1 + me2 *ny/np2
         yend(3)   = (me2+1)*ny/np2
         zstart(3) = 1
         zend(3)   = nz
      endif

c---------------------------------------------------------------------
c Set up info for blocking of ffts and transposes.  This improves
c performance on cache-based systems. Blocking involves
c working on a chunk of the problem at a time, taking chunks
c along the first, second, or third dimension. 
c
c - In cffts1 blocking is on 2nd dimension (with fft on 1st dim)
c - In cffts2/3 blocking is on 1st dimension (with fft on 2nd and 3rd dims)

c Since 1st dim is always in processor, we'll assume it's long enough 
c (default blocking factor is 16 so min size for 1st dim is 16)
c The only case we have to worry about is cffts1 in a 2d decomposition. 
c so the blocking factor should not be larger than the 2nd dimension. 
c---------------------------------------------------------------------

      fftblock = fftblock_default
      fftblockpad = fftblockpad_default

      if (layout_type .eq. layout_2d) then
         if (dims(2, 1) .lt. fftblock) fftblock = dims(2, 1)
         if (dims(2, 2) .lt. fftblock) fftblock = dims(2, 2)
         if (dims(2, 3) .lt. fftblock) fftblock = dims(2, 3)
      endif
      
      if (fftblock .ne. fftblock_default) fftblockpad = fftblock+3

      return
      end

      
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine compute_indexmap(twiddle, d1, d2, d3)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2 
c for time evolution exponent. 
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
c      include 'global.h'
      integer d1, d2, d3
      integer i, j, k, ii, ii2, jj, ij2, kk
      double precision ap, twiddle(d1, d2, d3)

c---------------------------------------------------------------------
c this function is very different depending on whether 
c we are in the 0d, 1d or 2d layout. Compute separately. 
c basically we want to convert the fortran indices 
c   1 2 3 4 5 6 7 8 
c to 
c   0 1 2 3 -4 -3 -2 -1
c The following magic formula does the trick:
c mod(i-1+n/2, n) - n/2
c---------------------------------------------------------------------

      ap = - 4.d0 * alpha * pi *pi

      if (layout_type .eq. layout_0d) then ! xyz layout
         do i = 1, dims(1,3)
            ii =  mod(i+xstart(3)-2+nx/2, nx) - nx/2
            ii2 = ii*ii
            do j = 1, dims(2,3)
               jj = mod(j+ystart(3)-2+ny/2, ny) - ny/2
               ij2 = jj*jj+ii2
               do k = 1, dims(3,3)
                  kk = mod(k+zstart(3)-2+nz/2, nz) - nz/2
                  twiddle(i,j,k) = dexp(ap*dfloat(kk*kk+ij2))
               end do
            end do
         end do
      else if (layout_type .eq. layout_1d) then ! zxy layout 
         do i = 1,dims(2,3)
            ii =  mod(i+xstart(3)-2+nx/2, nx) - nx/2
            ii2 = ii*ii
            do j = 1,dims(3,3)
               jj = mod(j+ystart(3)-2+ny/2, ny) - ny/2
               ij2 = jj*jj+ii2
               do k = 1,dims(1,3)
                  kk = mod(k+zstart(3)-2+nz/2, nz) - nz/2
                  twiddle(k,i,j) = dexp(ap*dfloat(kk*kk+ij2))
               end do
            end do
         end do
      else if (layout_type .eq. layout_2d) then ! zxy layout
         do i = 1,dims(2,3)
            ii =  mod(i+xstart(3)-2+nx/2, nx) - nx/2
            ii2 = ii*ii
            do j = 1, dims(3,3)
               jj = mod(j+ystart(3)-2+ny/2, ny) - ny/2
               ij2 = jj*jj+ii2
               do k =1,dims(1,3)
                  kk = mod(k+zstart(3)-2+nz/2, nz) - nz/2
                  twiddle(k,i,j) = dexp(ap*dfloat(kk*kk+ij2))
               end do
            end do
         end do
      else
         print *, ' Unknown layout type ', layout_type
         stop
      endif

      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine print_timers()

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer i, ierr
c      include 'global.h'
      include 'mpinpb.h'
      character*25 tstrings(T_max+2)
      double precision t1(T_max+2), tsum(T_max+2),
     >                 tming(T_max+2), tmaxg(T_max+2)
      data tstrings / '          total ', 
     >                '          vnpsetup ', 
     >                '            fft ', 
     >                '         evolve ', 
     >                '       checksum ', 
     >                '         fftlow ', 
     >                '        fftcopy ', 
     >                '      transpose ', 
     >                ' transpose1_loc ', 
     >                ' transpose1_glo ', 
     >                ' transpose1_fin ', 
     >                ' transpose2_loc ', 
     >                ' transpose2_glo ', 
     >                ' transpose2_fin ', 
     >                '           sync ',
     >                '           init ',
     >                '        totcomp ',
     >                '        totcomm ' /

      do i = 1, t_max
         t1(i) = timer_read(i)
      end do
      t1(t_max+2) = t1(t_transxzglo) + t1(t_transxyglo) + t1(t_synch)
      t1(t_max+1) = t1(t_total) - t1(t_max+2)

      call MPI_Reduce(t1, tsum,  t_max+2, MPI_DOUBLE_PRECISION, 
     >                MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(t1, tming, t_max+2, MPI_DOUBLE_PRECISION, 
     >                MPI_MIN, 0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(t1, tmaxg, t_max+2, MPI_DOUBLE_PRECISION, 
     >                MPI_MAX, 0, MPI_COMM_WORLD, ierr)

      if (me .ne. 0) return
      write(*, 800) np
      do i = 1, t_max+2
         if (tsum(i) .ne. 0.0d0) then
            write(*, 810) i, tstrings(i), tming(i), tmaxg(i), tsum(i)/np
         endif
      end do
 800  format(' nprocs =', i6, 19x, 'minimum', 5x, 'maximum', 
     >       5x, 'average')
 810  format(' timer ', i2, '(', A16, ') :', 3(2X,F10.4))
      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine fft(dir, x1, x2)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
c      include 'global.h'
      integer dir
      double complex x1(ntdivnp), x2(ntdivnp)

      double complex scratch(fftblockpad_default*maxdim*2)

c---------------------------------------------------------------------
c note: args x1, x2 must be different arrays
c note: args for cfftsx are (direction, layout, xin, xout, scratch)
c       xin/xout may be the same and it can be somewhat faster
c       if they are
c note: args for transpose are (layout1, layout2, xin, xout)
c       xin/xout must be different
c---------------------------------------------------------------------

      if (dir .eq. 1) then
         if (layout_type .eq. layout_0d) then
c            print *, "dir .eq. 1, 0d"
            call cffts1(1, dims(1,1), dims(2,1), dims(3,1), 
     >                  x1, x1, scratch)
            call cffts2(1, dims(1,2), dims(2,2), dims(3,2), 
     >                  x1, x1, scratch)
            call cffts3(1, dims(1,3), dims(2,3), dims(3,3), 
     >                  x1, x2, scratch)
         else if (layout_type .eq. layout_1d) then
c            print *, "dir .eq. 1, 1d"
            call cffts1(1, dims(1,1), dims(2,1), dims(3,1), 
     >                  x1, x1, scratch)
            call cffts2(1, dims(1,2), dims(2,2), dims(3,2), 
     >                  x1, x1, scratch)
            if (timers_enabled) call timer_start(T_transpose)
            call transpose_xy_z(2, 3, x1, x2)
            if (timers_enabled) call timer_stop(T_transpose)
            call cffts1(1, dims(1,3), dims(2,3), dims(3,3), 
     >                  x2, x2, scratch)
         else if (layout_type .eq. layout_2d) then
c            print *, "dir .eq. 1, 2d"
            call cffts1(1, dims(1,1), dims(2,1), dims(3,1), 
     >                  x1, x1, scratch)
            if (timers_enabled) call timer_start(T_transpose)
            call transpose_x_y(1, 2, x1, x2)
            if (timers_enabled) call timer_stop(T_transpose)
            call cffts1(1, dims(1,2), dims(2,2), dims(3,2), 
     >                  x2, x2, scratch)
            if (timers_enabled) call timer_start(T_transpose)
            call transpose_x_z(2, 3, x2, x1)
            if (timers_enabled) call timer_stop(T_transpose)
            call cffts1(1, dims(1,3), dims(2,3), dims(3,3), 
     >                  x1, x2, scratch)
         endif
      else
         if (layout_type .eq. layout_0d) then
c            print *, "dir .eq. -1, 0d"
            call cffts3(-1, dims(1,3), dims(2,3), dims(3,3), 
     >                  x1, x1, scratch)
            call cffts2(-1, dims(1,2), dims(2,2), dims(3,2), 
     >                  x1, x1, scratch)
            call cffts1(-1, dims(1,1), dims(2,1), dims(3,1), 
     >                  x1, x2, scratch)
         else if (layout_type .eq. layout_1d) then
c            print *, "dir .eq. -1, 1d"
            call cffts1(-1, dims(1,3), dims(2,3), dims(3,3), 
     >                  x1, x1, scratch)
            if (timers_enabled) call timer_start(T_transpose)
            call transpose_x_yz(3, 2, x1, x2)
            if (timers_enabled) call timer_stop(T_transpose)
            call cffts2(-1, dims(1,2), dims(2,2), dims(3,2), 
     >                  x2, x2, scratch)
            call cffts1(-1, dims(1,1), dims(2,1), dims(3,1), 
     >                  x2, x2, scratch)
         else if (layout_type .eq. layout_2d) then
c            print *, "dir .eq. -1, 2d"
            call cffts1(-1, dims(1,3), dims(2,3), dims(3,3), 
     >                  x1, x1, scratch)
            if (timers_enabled) call timer_start(T_transpose)
            call transpose_x_z(3, 2, x1, x2)
            if (timers_enabled) call timer_stop(T_transpose)
            call cffts1(-1, dims(1,2), dims(2,2), dims(3,2), 
     >                  x2, x2, scratch)
            if (timers_enabled) call timer_start(T_transpose)
            call transpose_x_y(2, 1, x2, x1)
            if (timers_enabled) call timer_stop(T_transpose)
            call cffts1(-1, dims(1,1), dims(2,1), dims(3,1), 
     >                  x1, x2, scratch)
         endif
      endif
      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine cffts1(is, d1, d2, d3, x, xout, y)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none

c      include 'global.h'
      integer is, d1, d2, d3, logd1
      double complex x(d1,d2,d3)
      double complex xout(d1,d2,d3)
      double complex y(fftblockpad, d1, 2) 
      integer i, j, k, jj

      logd1 = ilog2(d1)

      do k = 1, d3
         do jj = 0, d2 - fftblock, fftblock
            if (timers_enabled) call timer_start(T_fftcopy)
            do j = 1, fftblock
               do i = 1, d1
                  y(j,i,1) = x(i,j+jj,k)
               enddo
            enddo
            if (timers_enabled) call timer_stop(T_fftcopy)
            
            if (timers_enabled) call timer_start(T_fftlow)
            call cfftz (is, logd1, d1, y, y(1,1,2))
            if (timers_enabled) call timer_stop(T_fftlow)

            if (timers_enabled) call timer_start(T_fftcopy)
            do j = 1, fftblock
               do i = 1, d1
                  xout(i,j+jj,k) = y(j,i,1)
               enddo
            enddo
            if (timers_enabled) call timer_stop(T_fftcopy)
         enddo
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine cffts2(is, d1, d2, d3, x, xout, y)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none

c      include 'global.h'
      integer is, d1, d2, d3, logd2
      double complex x(d1,d2,d3)
      double complex xout(d1,d2,d3)
      double complex y(fftblockpad, d2, 2) 
      integer i, j, k, ii

      logd2 = ilog2(d2)

      do k = 1, d3
        do ii = 0, d1 - fftblock, fftblock
           if (timers_enabled) call timer_start(T_fftcopy)
           do j = 1, d2
              do i = 1, fftblock
                 y(i,j,1) = x(i+ii,j,k)
              enddo
           enddo
           if (timers_enabled) call timer_stop(T_fftcopy)

           if (timers_enabled) call timer_start(T_fftlow)
           call cfftz (is, logd2, d2, y, y(1, 1, 2))
           if (timers_enabled) call timer_stop(T_fftlow)

           if (timers_enabled) call timer_start(T_fftcopy)
           do j = 1, d2
              do i = 1, fftblock
                 xout(i+ii,j,k) = y(i,j,1)
              enddo
           enddo
           if (timers_enabled) call timer_stop(T_fftcopy)
        enddo
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine cffts3(is, d1, d2, d3, x, xout, y)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none

c      include 'global.h'
      integer is, d1, d2, d3, logd3
      double complex x(d1,d2,d3)
      double complex xout(d1,d2,d3)
      double complex y(fftblockpad, d3, 2) 
      integer i, j, k, ii

      logd3 = ilog2(d3)

      do j = 1, d2
        do ii = 0, d1 - fftblock, fftblock
           if (timers_enabled) call timer_start(T_fftcopy)
           do k = 1, d3
              do i = 1, fftblock
                 y(i,k,1) = x(i+ii,j,k)
              enddo
           enddo
           if (timers_enabled) call timer_stop(T_fftcopy)

           if (timers_enabled) call timer_start(T_fftlow)
           call cfftz (is, logd3, d3, y, y(1, 1, 2))
           if (timers_enabled) call timer_stop(T_fftlow)

           if (timers_enabled) call timer_start(T_fftcopy)
           do k = 1, d3
              do i = 1, fftblock
                 xout(i+ii,j,k) = y(i,k,1)
              enddo
           enddo
           if (timers_enabled) call timer_stop(T_fftcopy)
        enddo
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine fft_init (n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c compute the roots-of-unity array that will be used for subsequent FFTs. 
c---------------------------------------------------------------------
      use global
      implicit none
c      include 'global.h'

      integer m,n,nu,ku,i,j,ln
      double precision t, ti
        

c---------------------------------------------------------------------
c   Initialize the U array with sines and cosines in a manner that permits
c   stride one access at each FFT iteration.
c---------------------------------------------------------------------
      nu = n
      m = ilog2(n)
      u(1) = m
      ku = 2
      ln = 1

      do j = 1, m
         t = pi / ln
         
         do i = 0, ln - 1
            ti = i * t
            u(i+ku) = dcmplx (cos (ti), sin(ti))
         enddo
         
         ku = ku + ln
         ln = 2 * ln
      enddo
      
      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine cfftz (is, m, n, x, y)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   Computes NY N-point complex-to-complex FFTs of X using an algorithm due
c   to Swarztrauber.  X is both the input and the output array, while Y is a 
c   scratch array.  It is assumed that N = 2^M.  Before calling CFFTZ to 
c   perform FFTs, the array U must be initialized by calling CFFTZ with IS 
c   set to 0 and M set to MX, where MX is the maximum value of M for any 
c   subsequent call.
c---------------------------------------------------------------------
      use global
      implicit none
c      include 'global.h'

      integer is,m,n,i,j,l,mx
      double complex x, y
        
      dimension x(fftblockpad,n), y(fftblockpad,n)
c---------------------------------------------------------------------
c   Check if input parameters are invalid.
c---------------------------------------------------------------------
      mx = u(1)
      if ((is .ne. 1 .and. is .ne. -1) .or. m .lt. 1 .or. m .gt. mx)    
     >  then
        write (*, 1)  is, m, mx
 1      format ('CFFTZ: Either U has not been initialized, or else'/    
     >    'one of the input parameters is invalid', 3I5)
        stop
      endif

c---------------------------------------------------------------------
c   Perform one variant of the Stockham FFT.
c---------------------------------------------------------------------
      do l = 1, m, 2
        call fftz2 (is, l, m, n, fftblock, fftblockpad, u, x, y)
        if (l .eq. m) goto 160
        call fftz2 (is, l + 1, m, n, fftblock, fftblockpad, u, y, x)
      enddo

      goto 180

c---------------------------------------------------------------------
c   Copy Y to X.
c---------------------------------------------------------------------
 160  do j = 1, n
        do i = 1, fftblock
          x(i,j) = y(i,j)
        enddo
      enddo

 180  continue

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine fftz2 (is, l, m, n, ny, ny1, u, x, y)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   Performs the L-th iteration of the second variant of the Stockham FFT.
c---------------------------------------------------------------------
        use global, only:me
      implicit none

      integer is,k,l,m,n,ny,ny1,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22
      double complex u,x,y,u1,x11,x21
      dimension u(n), x(ny1,n), y(ny1,n)

c      print *, real(u(1)), imag(u(1)), u(1)

c      print *, is, l, m, n, ny, ny1

c      y(1, 1) = (1.5, 2.1)
c      y(2, 2) = (21.5, 22.1)

c      print *, "f y(1, 1) = ", y(1, 1)
c      print *, "f y(2, 2) = ", y(2, 2)
c      call c_fftz2 (is, l, m, n, ny, ny1, u(1), x(1, 1), y(1, 1))
c      return

c---------------------------------------------------------------------
c   Set initial parameters.
c---------------------------------------------------------------------

      n1 = n / 2
      lk = 2 ** (l - 1)
      li = 2 ** (m - l)
      lj = 2 * lk
      ku = li + 1

      do i = 0, li - 1
        i11 = i * lk + 1
        i12 = i11 + n1
        i21 = i * lj + 1
        i22 = i21 + lk
        if (is .ge. 1) then
          u1 = u(ku+i)
        else
          u1 = dconjg (u(ku+i))
        endif

c---------------------------------------------------------------------
c   This loop is vectorizable.
c   ni zhe ge da pian zi!
c---------------------------------------------------------------------
        do k = 0, lk - 1
          do j = 1, ny
c            print *, "i, k, j", i, k, j
            x11 = x(j,i11+k)
            x21 = x(j,i12+k)
            y(j,i21+k) = x11 + x21
            y(j,i22+k) = u1 * (x11 - x21)
          enddo
        enddo
      enddo




c      print *, "</y(1, 1) = ", y(1, 1)

      return
      end

      subroutine fortran_starstar(left_val, right_val, ret)
      implicit none
      integer left_val, right_val, ret
      ret = left_val ** right_val
      return
      end

      subroutine fortran_starstar_m1(left_val, right_val, ret)
      implicit none
      integer left_val, right_val, ret
      ret = left_val ** (right_val - 1)
      return
      end

      subroutine fortran_two_starstar_m1(m, li)
      implicit none
      integer m, li
      li = 2 ** (m - 1)
      !print *, "(- w -) (m - 1) = ", (m - 1)
      !print *, "(- w -) 2 ** (m - l) = ", 2 ** (m - l)
      return
      end

      subroutine fftz2_dbg (is, l, m, n, ny, ny1, u, x, y)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   Performs the L-th iteration of the second variant of the Stockham FFT.
c---------------------------------------------------------------------
        use global, only:me
      implicit none

      integer is,k,l,m,n,ny,ny1,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22
      double complex u,x,y,u1,x11,x21
      dimension u(n), x(ny1,n), y(ny1,n)

      integer ret

c---------------------------------------------------------------------
c   Set initial parameters.
c---------------------------------------------------------------------

      n1 = n / 2
      lk = 2 ** (l - 1)
      li = 2 ** (m - l)
      lj = 2 * lk
      ku = li + 1

      print *, is, l, m, n, ny, ny1
      print *, n1, lk, li, lj, ku

      do i = 0, li - 1
        i11 = i * lk + 1
        i12 = i11 + n1
        i21 = i * lj + 1
        i22 = i21 + lk
        if (is .ge. 1) then
          u1 = u(ku+i)
        else
          u1 = dconjg (u(ku+i))
        endif

c---------------------------------------------------------------------
c   This loop is vectorizable.
c   ni zhe ge da pian zi!
c---------------------------------------------------------------------
        do k = 0, lk - 1
          do j = 1, ny
            x11 = x(j,i11+k)
            x21 = x(j,i12+k)
            print *, i, k, j
            print *, "x11 = ", x11
            print *, "x11 = ", x21
            print *, "u1 = ", u1
            y(j,i21+k) = x11 + x21
            y(j,i22+k) = u1 * (x11 - x21)
            print *, "y21 = ", y(j,i21+k)
            print *, "y22 = ", y(j,i22+k)
          enddo
        enddo
      enddo




c      print *, "</y(1, 1) = ", y(1, 1)

      return
      end

c---------------------------------------------------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer function ilog2(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
      integer n, nn, lg
      if (n .eq. 1) then
         ilog2=0
         return
      endif
      lg = 1
      nn = 2
      do while (nn .lt. n)
         nn = nn*2
         lg = lg+1
      end do
      ilog2 = lg
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_yz(l1, l2, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer l1, l2
      double complex xin(ntdivnp), xout(ntdivnp)

      call transpose2_local(dims(1,l1),dims(2, l1)*dims(3, l1),
     >                          xin, xout)

      call transpose2_global(xout, xin)

      call transpose2_finish(dims(1,l1),dims(2, l1)*dims(3, l1),
     >                          xin, xout)

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_xy_z(l1, l2, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer l1, l2
      double complex xin(ntdivnp), xout(ntdivnp)

      call transpose2_local(dims(1,l1)*dims(2, l1),dims(3, l1),
     >                          xin, xout)
      call transpose2_global(xout, xin)
      call transpose2_finish(dims(1,l1)*dims(2, l1),dims(3, l1),
     >                          xin, xout)

      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose2_local(n1, n2, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
      integer n1, n2
      double complex xin(n1, n2), xout(n2, n1)
      
      double complex z(transblockpad, transblock)

      integer i, j, ii, jj

      if (timers_enabled) call timer_start(T_transxzloc)

c---------------------------------------------------------------------
c If possible, block the transpose for cache memory systems. 
c How much does this help? Example: R8000 Power Challenge (90 MHz)
c Blocked version decreases time spend in this routine 
c from 14 seconds to 5.2 seconds on 8 nodes class A.
c---------------------------------------------------------------------

      if (n1 .lt. transblock .or. n2 .lt. transblock) then
         if (n1 .ge. n2) then 
            do j = 1, n2
               do i = 1, n1
                  xout(j, i) = xin(i, j)
               end do
            end do
         else
            do i = 1, n1
               do j = 1, n2
                  xout(j, i) = xin(i, j)
               end do
            end do
         endif
      else
         do j = 0, n2-1, transblock
            do i = 0, n1-1, transblock
               
c---------------------------------------------------------------------
c Note: compiler should be able to take j+jj out of inner loop
c---------------------------------------------------------------------
               do jj = 1, transblock
                  do ii = 1, transblock
                     z(jj,ii) = xin(i+ii, j+jj)
                  end do
               end do
               
               do ii = 1, transblock
                  do jj = 1, transblock
                     xout(j+jj, i+ii) = z(jj,ii)
                  end do
               end do
               
            end do
         end do
      endif
      if (timers_enabled) call timer_stop(T_transxzloc)

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose2_global(xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
      double complex xin(ntdivnp)
      double complex xout(ntdivnp) 
      integer ierr

!      if (timers_enabled) call synchup()

      if (timers_enabled) call timer_start(T_transxzglo)
      call mpi_alltoall(xin, ntdivnp/np, dc_type,
     >                  xout, ntdivnp/np, dc_type,
     >                  commslice1, ierr)
      if (timers_enabled) call timer_stop(T_transxzglo)

      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose2_finish(n1, n2, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer n1, n2, ioff
      double complex xin(n2, n1/np2, 0:np2-1), xout(n2*np2, n1/np2)
      
      integer i, j, p

      if (timers_enabled) call timer_start(T_transxzfin)
      do p = 0, np2-1
         ioff = p*n2
         do j = 1, n1/np2
            do i = 1, n2
               xout(i+ioff, j) = xin(i, j, p)
            end do
         end do
      end do
      if (timers_enabled) call timer_stop(T_transxzfin)

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_z(l1, l2, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer l1, l2
      double complex xin(ntdivnp), xout(ntdivnp)

      call transpose_x_z_local(dims(1,l1),dims(2,l1),dims(3,l1),
     >                         xin, xout)
      call transpose_x_z_global(dims(1,l1),dims(2,l1),dims(3,l1), 
     >                          xout, xin)
      call transpose_x_z_finish(dims(1,l2),dims(2,l2),dims(3,l2), 
     >                          xin, xout)
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_z_local(d1, d2, d3, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer d1, d2, d3
      double complex xin(d1,d2,d3)
      double complex xout(d3,d2,d1)
      integer block1, block3
      integer i, j, k, kk, ii, i1, k1

      double complex buf(transblockpad, maxdim)
      if (timers_enabled) call timer_start(T_transxzloc)
      if (d1 .lt. 32) goto 100
      block3 = d3
      if (block3 .eq. 1)  goto 100
      if (block3 .gt. transblock) block3 = transblock
      block1 = d1
      if (block1*block3 .gt. transblock*transblock) 
     >          block1 = transblock*transblock/block3
c---------------------------------------------------------------------
c blocked transpose
c---------------------------------------------------------------------
      do j = 1, d2
         do kk = 0, d3-block3, block3
            do ii = 0, d1-block1, block1
               
               do k = 1, block3
                  k1 = k + kk
                  do i = 1, block1
                     buf(k, i) = xin(i+ii, j, k1)
                  end do
               end do

               do i = 1, block1
                  i1 = i + ii
                  do k = 1, block3
                     xout(k+kk, j, i1) = buf(k, i)
                  end do
               end do

            end do
         end do
      end do
      goto 200
      

c---------------------------------------------------------------------
c basic transpose
c---------------------------------------------------------------------
 100  continue
      
      do j = 1, d2
         do k = 1, d3
            do i = 1, d1
               xout(k, j, i) = xin(i, j, k)
            end do
         end do
      end do

c---------------------------------------------------------------------
c all done
c---------------------------------------------------------------------
 200  continue

      if (timers_enabled) call timer_stop(T_transxzloc)
      return 
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_z_global(d1, d2, d3, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
      integer d1, d2, d3
      double complex xin(d3,d2,d1)
      double complex xout(d3,d2,d1) ! not real layout, but right size
      integer ierr

!      if (timers_enabled) call synchup()

c---------------------------------------------------------------------
c do transpose among all  processes with same 1-coord (me1)
c---------------------------------------------------------------------
      if (timers_enabled)call timer_start(T_transxzglo)
      call mpi_alltoall(xin, d1*d2*d3/np2, dc_type,
     >                  xout, d1*d2*d3/np2, dc_type,
     >                  commslice1, ierr)
      if (timers_enabled) call timer_stop(T_transxzglo)
      return
      end
      

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_z_finish(d1, d2, d3, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer d1, d2, d3
      double complex xin(d1/np2, d2, d3, 0:np2-1)
      double complex xout(d1,d2,d3)
      integer i, j, k, p, ioff
      if (timers_enabled) call timer_start(T_transxzfin)
c---------------------------------------------------------------------
c this is the most straightforward way of doing it. the
c calculation in the inner loop doesn't help. 
c      do i = 1, d1/np2
c         do j = 1, d2
c            do k = 1, d3
c               do p = 0, np2-1
c                  ii = i + p*d1/np2
c                  xout(ii, j, k) = xin(i, j, k, p)
c               end do
c            end do
c         end do
c      end do
c---------------------------------------------------------------------

      do p = 0, np2-1
         ioff = p*d1/np2
         do k = 1, d3
            do j = 1, d2
               do i = 1, d1/np2
                  xout(i+ioff, j, k) = xin(i, j, k, p)
               end do
            end do
         end do
      end do
      if (timers_enabled) call timer_stop(T_transxzfin)
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_y(l1, l2, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer l1, l2
      double complex xin(ntdivnp), xout(ntdivnp)

c---------------------------------------------------------------------
c xy transpose is a little tricky, since we don't want
c to touch 3rd axis. But alltoall must involve 3rd axis (most 
c slowly varying) to be efficient. So we do
c (nx, ny/np1, nz/np2) -> (ny/np1, nz/np2, nx) (local)
c (ny/np1, nz/np2, nx) -> ((ny/np1*nz/np2)*np1, nx/np1) (global)
c then local finish. 
c---------------------------------------------------------------------


      call transpose_x_y_local(dims(1,l1),dims(2,l1),dims(3,l1),
     >                         xin, xout)
      call transpose_x_y_global(dims(1,l1),dims(2,l1),dims(3,l1), 
     >                          xout, xin)
      call transpose_x_y_finish(dims(1,l2),dims(2,l2),dims(3,l2), 
     >                          xin, xout)

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_y_local(d1, d2, d3, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer d1, d2, d3
      double complex xin(d1, d2, d3)
      double complex xout(d2, d3, d1)
      integer i, j, k
      if (timers_enabled) call timer_start(T_transxyloc)

      do k = 1, d3
         do i = 1, d1
            do j = 1, d2
               xout(j,k,i)=xin(i,j,k)
            end do
         end do
      end do
      if (timers_enabled) call timer_stop(T_transxyloc)
      return 
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_y_global(d1, d2, d3, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
      integer d1, d2, d3
c---------------------------------------------------------------------
c array is in form (ny/np1, nz/np2, nx)
c---------------------------------------------------------------------
      double complex xin(d2,d3,d1)
      double complex xout(d2,d3,d1) ! not real layout but right size
      integer ierr

!      if (timers_enabled) call synchup()

c---------------------------------------------------------------------
c do transpose among all processes with same 1-coord (me1)
c---------------------------------------------------------------------
      if (timers_enabled) call timer_start(T_transxyglo)
      call mpi_alltoall(xin, d1*d2*d3/np1, dc_type,
     >                  xout, d1*d2*d3/np1, dc_type,
     >                  commslice2, ierr)
      if (timers_enabled) call timer_stop(T_transxyglo)

      return
      end
      

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine transpose_x_y_finish(d1, d2, d3, xin, xout)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      integer d1, d2, d3
      double complex xin(d1/np1, d3, d2, 0:np1-1)
      double complex xout(d1,d2,d3)
      integer i, j, k, p, ioff
      if (timers_enabled) call timer_start(T_transxyfin)
c---------------------------------------------------------------------
c this is the most straightforward way of doing it. the
c calculation in the inner loop doesn't help. 
c      do i = 1, d1/np1
c         do j = 1, d2
c            do k = 1, d3
c               do p = 0, np1-1
c                  ii = i + p*d1/np1
c note order is screwy bcz we have (ny/np1, nz/np2, nx) -> (ny, nx/np1, nz/np2)
c                  xout(ii, j, k) = xin(i, k, j, p)
c               end do
c            end do
c         end do
c      end do
c---------------------------------------------------------------------

      do p = 0, np1-1
         ioff = p*d1/np1
         do k = 1, d3
            do j = 1, d2
               do i = 1, d1/np1
                  xout(i+ioff, j, k) = xin(i, k, j, p)
               end do
            end do
         end do
      end do
      if (timers_enabled) call timer_stop(T_transxyfin)
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine checksum(i, u1, d1, d2, d3)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
      integer i, d1, d2, d3
      double complex u1(d1, d2, d3)
      integer j, q,r,s, ierr
      double complex chk,allchk
      chk = (0.0,0.0)
      do j=1,1024
         q = mod(j, nx)+1
         if (q .ge. xstart(1) .and. q .le. xend(1)) then
            r = mod(3*j,ny)+1
            if (r .ge. ystart(1) .and. r .le. yend(1)) then
               s = mod(5*j,nz)+1
               if (s .ge. zstart(1) .and. s .le. zend(1)) then
                  chk=chk+u1(q-xstart(1)+1,r-ystart(1)+1,s-zstart(1)+1)
               end if
            end if
         end if
      end do
      chk = chk/ntotal_f

      if (timers_enabled) call timer_start(T_synch)
      call MPI_Reduce(chk, allchk, 1, dc_type, MPI_SUM, 
     >                0, MPI_COMM_WORLD, ierr)      
      if (timers_enabled) call timer_stop(T_synch)
      if (me .eq. 0) then
            write (*, 30) i, allchk
 30         format (' T =',I5,5X,'Checksum =',1P2D22.12)
      endif

      sums(i) = allchk
      if(me .eq. 0)  write(*,*) sums(i)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine synchup

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use global
      implicit none
      include 'mpinpb.h'
      integer ierr
      call timer_start(T_synch)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call timer_stop(T_synch)
      return
      end


