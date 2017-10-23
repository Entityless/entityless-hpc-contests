module global
      include 'npbparams.h'

      integer nx,ny,nz,maxdim
      integer niter_default
      integer np_min
      integer ntdivnp
      double precision ntotal_f

      integer np1, np2, np
     
      integer layout_type
      integer layout_0D, layout_1D, layout_2D
      parameter (layout_0D = 0, layout_1D = 1, layout_2D = 2)



      integer fftblock_default, fftblockpad_default
      parameter (fftblock_default=16, fftblockpad_default=18)
      integer transblock, transblockpad
      parameter(transblock=32, transblockpad=34)
      
      integer fftblock, fftblockpad

      integer me, me1, me2
      integer commslice1, commslice2



      integer dims(3, 3)
      integer xstart(3), ystart(3), zstart(3)
      integer xend(3), yend(3), zend(3)

      integer T_total, T_setup, T_fft, T_evolve, T_checksum,   &
             T_fftlow, T_fftcopy, T_transpose, &
             T_transxzloc, T_transxzglo, T_transxzfin, &
             T_transxyloc, T_transxyglo, T_transxyfin, &
             T_synch, T_init, T_max
      parameter (T_total = 1, T_setup = 2, T_fft = 3, &
                T_evolve = 4, T_checksum = 5, &
                T_fftlow = 6, T_fftcopy = 7, T_transpose = 8,&
                T_transxzloc = 9, T_transxzglo = 10, T_transxzfin = 11, &
                T_transxyloc = 12, T_transxyglo = 13, &
                T_transxyfin = 14,  T_synch = 15, T_init = 16,&
                T_max = 16)



      logical timers_enabled


      external timer_read
      double precision timer_read
      external ilog2
      integer ilog2

      external randlc
      double precision randlc


      logical debug, debugsynch

      double precision seed, a, pi, alpha
      parameter (seed = 314159265.d0, a = 1220703125.d0, &
       pi = 3.141592653589793238d0, alpha=1.0d-6)

      double complex,allocatable ::u(:)


!      double complex sums(niter_default/100)
      double complex,allocatable ::sums(:)

      integer niter
end module
