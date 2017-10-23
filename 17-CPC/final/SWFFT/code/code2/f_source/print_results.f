
      subroutine print_results(name, class, n1, n2, n3, niter, 
     >               nprocs_compiled, nprocs_total,
     >               t, mops, optype, verified, npbversion, 
     >               compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)
      use global  , only : nx,ny,nz,np,sums    
      implicit none
      character*2 name
      character*1 class
      integer n1, n2, n3, niter, nprocs_compiled, nprocs_total, j
      double precision t, mops
      character optype*24, size*15
      logical verified
      character*(*) npbversion, compiletime, 
     >              cs1, cs2, cs3, cs4, cs5, cs6, cs7

         write (*, 2) 
 2       format(//, ' ', A2, ' Benchmark Completed.')


         if ((n2 .eq. 0) .and. (n3 .eq. 0)) then
            if (name(1:2) .eq. 'EP') then
               write(size, '(f15.0)' ) 2.d0**n1
               j = 15
               if (size(j:j) .eq. '.') j = j - 1
               write (*,42) size(1:j)
 42            format(' Size            = ',9x, a15)
            else
               write (*,44) n1
 44            format(' Size            = ',12x, i12)
            endif
         else
            write (*, 4) n1,n2,n3
 4          format(' Size            =  ',9x, i4,'x',i4,'x',i4)
         endif

         write (*, 5) niter
 5       format(' Iterations      = ', 12x, i12)
         
         write (*, 6) t
 6       format(' Time in seconds = ',12x, f12.2)
         
         write (*,7) nprocs_total
 7       format(' Total processes = ', 12x, i12)
         
         write (*,8) nprocs_compiled
 8       format(' Compiled procs  = ', 12x, i12)

         write (*,9) mops
 9       format(' Mop/s total     = ',12x, f12.2)

         write (*,10) mops/float( nprocs_total )
 10      format(' Mop/s/process   = ', 12x, f12.2)        
         
         write(*, 11) optype
 11      format(' Operation type  = ', a24)

        open(155,file='fftout',form='unformatted')
!        open(155,file='fftout')
        write(155) nx,ny,nz,niter,np
        write(155) sums
        close(155)

         return
         end

