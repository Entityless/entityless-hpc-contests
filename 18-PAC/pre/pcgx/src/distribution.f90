!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module distribution

!BOP
! !MODULE: distribution
!
! !DESCRIPTION:
! This module provides data types and routines for distributing
! blocks across processors.
!
! !REVISION HISTORY:
! CVS:$Id: distribution.F90,v 1.11 2003/12/23 22:11:40 pwjones Exp $
! CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use communicate
   use broadcast
   use blocks
   use exit_mod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: & !
      nblock_mic_pp = 1 ,&! block# per MIC process/rank
      nproc_cpu_pn = 28 ,&! MPI process# on CPU per node
      nproc_mic_pn = 0 ,&! MPI process# on MIC per node
      nproc_pn = 28 ,&! MPI process# per node
      nnode = 2 ! total node number

   real(r8), parameter, public :: & !
      k_cpu_mic = 3.0_r8 ! quick hack on block distribution

! !PUBLIC TYPES:

   type, public :: distrb ! distribution data type
      integer (int_kind) :: &
         nprocs ,&! number of processors in this dist
         communicator ! communicator to use in this dist

      integer (int_kind), dimension(:), pointer :: &
         proc ,&! processor location for this block
         local_block ! block position in local array on proc

      integer (int_kind), dimension(:), pointer :: &
         local_block_ids ! global id for blocks in local process

      integer (int_kind) :: &
         local_block_num ! number of blocks for local process
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution, &
             create_distribution_tiny, &
             create_local_block_ids, &
             create_local_block_ids_tiny, &
             init_mic_proc_flag

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 function create_distribution(dist_type, nprocs, work_per_block)

! !DESCRIPTION:
! This routine determines the distribution of blocks across processors
! by call the appropriate subroutine based on distribution type
! requested. Currently only two distributions are supported:
! 2-d Cartesian distribution (cartesian) and a load-balanced
! distribution (balanced) based on an input amount of work per
! block.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      dist_type ! method for distributing blocks
                            ! either cartesian or balanced

   integer (int_kind), intent(in) :: &
      nprocs ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distribution ! resulting structure describing
                            ! distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
! select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('cartesian')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('balanced')

      create_distribution = create_distrb_balanced(nprocs, &
                                                   work_per_block)

   ! ljm tuning
   case('cart_cpuonly')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('cart_micmixed')

      create_distribution = create_distrb_mic(nprocs, work_per_block)

   case default

      call exit_POP(sigAbort,'distribution: unknown distribution type')

   end select

!-----------------------------------------------------------------------
!EOC

 end function create_distribution

!***********************************************************************
!BOP
! !IROUTINE: init_mic_proc_flag
! !INTERFACE:

 subroutine init_mic_proc_flag()

! !DESCRIPTION:
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC

   lmic_proc = (my_task >= nproc_cpu_pn * nnode)
   !lmic_proc = (my_task == nproc_cpu_pn * nnode)
   !lmic_proc = .false.
   lmic_trace = .false.

!EOC

 end subroutine init_mic_proc_flag


!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids
! !INTERFACE:

 subroutine create_local_block_ids(block_ids, distribution)

! !DESCRIPTION:
! This routine determines which blocks in an input distribution are
! located on the local processor and creates an array of block ids
! for all local blocks.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution ! input distribution for which local
                             ! blocks required

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), pointer :: &
      block_ids ! array of block ids for every block
                             ! that resides on the local processor

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bid, bcount ! dummy counters

!-----------------------------------------------------------------------
!
! first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%proc)
      if (distribution%proc(n) == my_task+1) bcount = bcount + 1
   end do

   !ljm tuning
   write(6,*) 'bcount ', my_task, bcount
   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
! now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      do n=1,size(distribution%proc)
         if (distribution%proc(n) == my_task+1) then
            block_ids(distribution%local_block(n)) = n
         endif
      end do
   endif

!EOC

 end subroutine create_local_block_ids

!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids_tiny
! !INTERFACE:

 subroutine create_local_block_ids_tiny(block_ids, distribution)

! !DESCRIPTION:
! This routine determines which blocks in an input distribution are
! located on the local processor and creates an array of block ids
! for all local blocks.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution ! input distribution for which local
                             ! blocks required

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), pointer :: &
      block_ids ! array of block ids for every block
                             ! that resides on the local processor

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bid, bcount ! dummy counters

!-----------------------------------------------------------------------
!
! first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = nblocks_tot_tiny
   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
! now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      do n=1,bcount
         block_ids(n) = n + nblocks_tot ! global id
      end do
   endif

!EOC

 end subroutine create_local_block_ids_tiny


!***********************************************************************
!BOP
! !IROUTINE: create_distrb_cart
! !INTERFACE:

 function create_distrb_cart(nprocs, work_per_block)

! !DESCRIPTION:
! This function creates a distribution of blocks across processors
! using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_cart ! resulting structure describing Cartesian
                          ! distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n ,&! dummy loop indices
      iblock, jblock, nblck ,&!
      is, ie, js, je ,&! start, end block indices for each proc
      local_block ,&! block location on this processor
      nprocs_x ,&! num of procs in x for global domain
      nprocs_y ,&! num of procs in y for global domain
      nblocks_x_loc ,&! num of blocks per processor in x
      nblocks_y_loc ! num of blocks per processor in y

   type (distrb) :: dist ! temp hold distribution

!----------------------------------------------------------------------
!
! create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
! try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

   ! ljm tuning
   if (my_task == nproc_cpu_pn*nnode) then ! first MIC process
   !if (my_task == master_task) then
    write(6,*) 'dist_cart: nprocs,nprocs_x,nprocs_y',nprocs,nprocs_x,nprocs_y
    write(6,*) '         : nblocks_tot,nblocks_x,nblocks_y',nblocks_tot,nblocks_x,nblocks_y
   endif
!----------------------------------------------------------------------
!
! allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc (nblocks_tot), &
             dist%local_block(nblocks_tot))

!----------------------------------------------------------------------
!
! distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   nblocks_x_loc = (nblocks_x-1)/nprocs_x + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_y + 1

   do j=1,nprocs_y
   do i=1,nprocs_x
      n = (j-1)*nprocs_x + i

      is = (i-1)*nblocks_x_loc + 1
      ie = i *nblocks_x_loc
      if (ie > nblocks_x) ie = nblocks_x
      js = (j-1)*nblocks_y_loc + 1
      je = j *nblocks_y_loc
      if (je > nblocks_y) je = nblocks_y

      local_block = 0
      do jblock = js,je
      do iblock = is,ie
         nblck = (jblock - 1)*nblocks_x + iblock
         if (work_per_block(nblck) /= 0) then
            local_block = local_block + 1
            dist%proc(nblck) = n
            dist%local_block(nblck) = local_block
         else
            dist%proc(nblck) = 0
            dist%local_block(nblck) = 0
         endif
      end do
      end do
   end do
   end do

!----------------------------------------------------------------------

   create_distrb_cart = dist ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_cart

!**********************************************************************
!BOP
! !IROUTINE: create_distribution_tiny
! !INTERFACE:

 subroutine create_distribution_tiny(dist)

! !DESCRIPTION:
! This subroutine append the tiny blocks that derived from the normal blocks
! distributed to this process.
!
! !REVISION HISTORY:
! same as module

! !INPUT/OUTPUT PARAMETERS:

   type (distrb), intent(inout) :: &
      dist ! resulting structure describing Cartesian

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      n ! loop indices

   integer (int_kind),dimension(:),allocatable :: &
      local_proc,local_block ! tmp buffer for reallocate

   allocate(local_proc(nblocks_tot), &
            local_block(nblocks_tot))

   local_proc(1:nblocks_tot) = dist%proc(1:nblocks_tot)
   local_block(1:nblocks_tot) = dist%local_block(1:nblocks_tot)
   deallocate(dist%proc, dist%local_block)
   allocate(dist%proc(nblocks_tot+nblocks_tot_tiny), &
            dist%local_block(nblocks_tot+nblocks_tot_tiny))
   dist%proc(1:nblocks_tot) = local_proc(1:nblocks_tot)
   dist%local_block(1:nblocks_tot) = local_block(1:nblocks_tot)
   deallocate(local_proc, local_block)

!----------------------------------------------------------------------
! Reaplace the normal(macro) blocks with tiny blocks ?
! Or keep the original normal blocks info somewhere?
   do n=1,nblocks_tot_tiny
      dist%proc(n+nblocks_tot) = my_task + 1
      dist%local_block(n+nblocks_tot) = n
   end do

!----------------------------------------------------------------------
!EOC

 end subroutine create_distribution_tiny


!***********************************************************************
!BOP
! !IROUTINE: create_distrb_mic
! !INTERFACE:

 function create_distrb_mic(nprocs, work_per_block)

! !DESCRIPTION:
! This function creates a distribution of blocks across processors
! using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_mic ! resulting structure describing Cartesian
                          ! distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------
   character (char_len) :: outstring

   integer (int_kind) :: &
      i, j ,&! dummy loop indices
      n, n_mic ,&! dummy loop indices
      nprocs_cpu ,&! total num of cpu processes
      nblocks_mic ,&! dummy counter
      nbmap_pos ,&! stack top pointer to mic_block_map
      nbmap_delta ,&! adjustment for stack top pointer
      max_micblock_pp ,&! dummy counter
      max_micproc_pp ,&! max block# to [1,max_micproc_pp] MIC ranks,
      max_cpuproc_pp ,&! the following got max-1 block#
      cur_nmicblocks_pp ,&! block# for MIC rank(s) in current tile
      cur_ncpublocks_pp ,&! block# for CPU rank in current tile
      cur_nmicblocks_counter,&! dummy counter
      iblock, jblock, nblck ,&!
      nblocks_eff ,&! num of non-null blocks
      nblocks_avgeff ,&! num of non-null blocks per tile
      quota_tile,quota_mic ,&! accumualted diff. from max load per tile
      ncols ,&! num of olumns that contains at least one non-null block
      i_bcount ,&! num of non-null blocks in current column
      nstretch_y ,&! num of blocks per tile in y
      is, ie, istep ,&! start, end block indices for each proc
      local_block ,&! block location on this processor
      local_block_cpu ,&! block location on this processor
      local_block_mic ,&! block location on this processor
      nprocs_x ,&! num of procs in x for global domain
      nprocs_y ,&! num of procs in y for global domain
      nblocks_loc ,&! num of blocks per processor
      nblocks_x_loc ,&! num of blocks per processor in x
      nblocks_y_loc ! num of blocks per processor in y

   logical (log_kind):: &
      leff_col ,&! scan order in x direction(interleaving)
      reverse_order ,&! scan order in x direction(interleaving)
      rectshape ! scan order in x direction(interleaving)

   integer (int_kind),dimension(:),allocatable :: &
      mic_block_map ! block location on this processor

   integer (int_kind),dimension(:,:),allocatable :: &
      dist_proc ! block location on this processor

   type (distrb) :: dist ! temp hold distribution

!----------------------------------------------------------------------
!
! create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
! try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs
   nprocs_cpu = nproc_cpu_pn * nnode

   call proc_decomposition(nprocs_cpu, nprocs_x, nprocs_y)

   ! avoid 0 block distrb
   nblocks_x_loc = (nblocks_x-1)/nprocs_x + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_y + 1
   nblck = nblocks_x_loc*nprocs_x * &! virtual total size,small is better
           nblocks_y_loc*nprocs_y
   ! test tranpose
   nblocks_x_loc = (nblocks_x-1)/nprocs_y + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_x + 1
   if (nblck > nblocks_x_loc * nblocks_y_loc) then
      nblck = nprocs_x
      nprocs_x = nprocs_y
      nprocs_y = nblck
   endif
! iblock = (nblocks_x-1)/nblocks_x_loc + 1
! jblock = (nblocks_y-1)/nblocks_y_loc + 1
! ! swap nprocs_x and nprocs_y
! if (iblock < nprocs_x .or. jblock < nprocs_y) then
! nblck = nprocs_x
! nprocs_x = nprocs_y
! nprocs_y = nblck
! endif

   !if (my_task == nproc_cpu_pn*nnode) then ! first MIC process
   if (my_task == master_task) then
    write(6,*) 'dist_cart: nprocs,nprocs_x,nprocs_y',nprocs_cpu,nprocs_x,nprocs_y
    write(6,*) '         : nblocks_tot,nblocks_x,nblocks_y',nblocks_tot,nblocks_x,nblocks_y
   endif
!----------------------------------------------------------------------
!
! allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc (nblocks_tot), &
             dist%local_block(nblocks_tot))
   dist%proc = 0
   dist%local_block = 0;

   nblocks_x_loc = (nblocks_x-1)/nprocs_x + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_y + 1
   nblocks_loc = nblocks_x_loc * nblocks_y_loc
!----------------------------------------------------------------------
!
! collect the histogram of block distribution in advance
!----------------------------------------------------------------------
   if (nproc_mic_pn /= 0) then
   nblocks_mic = nblock_mic_pp * nproc_mic_pn * nnode
   max_micblock_pp = (nblock_mic_pp * nproc_mic_pn - 1)/nproc_cpu_pn + 1
   max_micproc_pp = nblock_mic_pp*nproc_mic_pn - (max_micblock_pp-1)*nproc_cpu_pn
   allocate (mic_block_map(nblocks_mic))
   mic_block_map = 0
   nbmap_pos = 1
   else
   nblocks_mic = 0
   max_micblock_pp = 0
   max_micproc_pp = 0
   endif
! nblocks_eff = 0
! do j=1,nblocks_tot
! if (work_per_block(j) /= 0) then
! nblocks_eff = nblocks_eff + 1
! endif
! end do
   nblocks_eff = count(work_per_block /= 0)
   nblocks_avgeff = (nblocks_eff - 1)/nprocs_cpu + 1
   max_cpuproc_pp = nblocks_eff - (nblocks_avgeff-1)*nprocs_cpu
   quota_tile = nprocs_cpu - max_cpuproc_pp
   quota_mic = nproc_mic_pn - max_micproc_pp
   nstretch_y = floor(sqrt(real(nblocks_avgeff))) ! previous use floor
   if (nstretch_y < nblocks_avgeff/nstretch_y) then
      nstretch_y = nstretch_y + 1
   endif
   if ((nblocks_avgeff/nstretch_y)*nstretch_y == nblocks_avgeff .and. &
       quota_tile > 0 .and. &
       nblocks_avgeff > 2) then
      rectshape = .true.
   else
      rectshape = .false.
   endif

   !if (my_task == nproc_cpu_pn*nnode) then ! first MIC process
   if (my_task == master_task) then
    write(6,*) 'dist parameters:',nblocks_eff,nblocks_avgeff,max_cpuproc_pp,nstretch_y
    write(6,*) '               :',max_micblock_pp,max_micproc_pp,quota_tile,quota_mic,rectshape
   endif
   n = 1 ! started processor
   n_mic = 1 ! started mic processor
   local_block = 0 ! init counter
   local_block_cpu = 0 ! init counter
   local_block_mic = 0 ! init counter
   cur_nmicblocks_counter = 0
   cur_nmicblocks_pp = max_micblock_pp
   cur_ncpublocks_pp = nblocks_avgeff
   reverse_order = .false.

   ! stats on effective columns to avoid over-distrb
   ncols = 0
   do j=1,nblocks_y,nstretch_y
   do i=1,nblocks_x
      leff_col = .false.
      effcol_loop: do jblock = j,min(j+nstretch_y-1,nblocks_y)
         nblck = (jblock - 1)*nblocks_x + i
         if (work_per_block(nblck) /= 0) then
            leff_col = .true.
            exit effcol_loop
         endif
      end do effcol_loop
      if (leff_col) &
         ncols = ncols + 1
   end do
   end do

   do j=1,nblocks_y,nstretch_y
   if (reverse_order) then
      is = nblocks_x
      ie = 1
      istep = -1
   else
      is = 1
      ie = nblocks_x
      istep = 1
   endif
   reverse_order = .not. reverse_order
   do i=is,ie,istep
      i_bcount = 0
      do jblock = j,min(j+nstretch_y-1,nblocks_y)
         nblck = (jblock - 1)*nblocks_x + i
         if (work_per_block(nblck) /= 0) then
            i_bcount = i_bcount + 1
         endif
      end do

   !if (my_task == nproc_cpu_pn*nnode) &! first MIC process
   !if (my_task == master_task) &
   ! write(6,*) 'j,i,i_bc,cnp,lblk,quot:',j,i,i_bcount,cur_ncpublocks_pp,local_block,quota_tile

      ! indented with diff. proc: discard cur column or go ahead
      if (rectshape .and. &
          i_bcount >= cur_ncpublocks_pp - local_block &
         ) then
         ! can discard
         if (i_bcount > cur_ncpublocks_pp - local_block .and. &
             quota_tile >= cur_ncpublocks_pp - local_block &
            ) then
            quota_tile = quota_tile - (cur_ncpublocks_pp - local_block)
            ! decrease mic blocks in current tile
            nbmap_delta = min(max_micblock_pp, cur_ncpublocks_pp-local_block)
            nbmap_delta = min(nbmap_delta, quota_mic)
            quota_mic = quota_mic - nbmap_delta
            ! replace mic with cpu proc
            do jblock=nbmap_pos - nbmap_delta, nbmap_pos-1
               nblck = mic_block_map(jblock)
               local_block_cpu = local_block_cpu + 1
               dist%proc(nblck) = n
               dist%local_block(nblck) = local_block_cpu
            end do
            nbmap_pos = nbmap_pos - nbmap_delta
   !if (my_task == nproc_cpu_pn*nnode) &! first MIC process
   !if (my_task == master_task) &
   ! write(6,*) 'def nbmap_pos:',nbmap_pos

            n = n + 1
            local_block = 0
            local_block_cpu = 0
            cur_nmicblocks_counter = 0
            ! block# for current tile: MIC and CPU
            cur_nmicblocks_pp = min(nblocks_mic-nbmap_pos+1,max_micblock_pp)
   !if (my_task == nproc_cpu_pn*nnode) &! first MIC process
   !if (my_task == master_task) &
   ! write(6,*) 'reset cur_nmicblocks_pp:',cur_nmicblocks_pp

            cur_ncpublocks_pp = nblocks_avgeff
         endif
      endif ! rectshape .and. ,etc
      do jblock = j,min(j+nstretch_y-1,nblocks_y)
         nblck = (jblock - 1)*nblocks_x + i
         if (work_per_block(nblck) /= 0) then
            ! distrb to mic
            if (cur_nmicblocks_counter < cur_nmicblocks_pp .and. &
                (jblock /= j .and. jblock /= min(j+nstretch_y-1,nblocks_y) .and. &
                 local_block_cpu > nstretch_y &
                 .or. &
                 cur_ncpublocks_pp - local_block <= nstretch_y &
                 .or. &
                 ncols<=nprocs_cpu-n+2) &
               ) then
              local_block_mic = local_block_mic + 1
              dist%proc(nblck) = n_mic + nprocs_cpu
              dist%local_block(nblck) = local_block_mic
              cur_nmicblocks_counter = cur_nmicblocks_counter + 1
              if (rectshape .and. &
                  nbmap_pos <= nblocks_mic) then
                 mic_block_map(nbmap_pos) = nblck
                 nbmap_pos = nbmap_pos + 1
   !if (my_task == nproc_cpu_pn*nnode) & ! first MIC process
   !if (my_task == master_task) &
   ! write(6,*) 'inc nbmap_pos:',nbmap_pos

              endif
              if (local_block_mic == nblock_mic_pp) then
                n_mic = n_mic + 1
                local_block_mic = 0
              endif
            else
            ! distrb to cpu
              local_block_cpu = local_block_cpu + 1
              dist%proc(nblck) = n
              dist%local_block(nblck) = local_block_cpu
            endif
            local_block = local_block + 1
            if (local_block == cur_ncpublocks_pp) then
               n = n + 1
               local_block = 0
               local_block_cpu = 0
               cur_nmicblocks_counter = 0
               ! block# for current tile: MIC and CPU
               if (.not. rectshape) then
                  if (mod(n-1,nproc_cpu_pn)+1 <= max_micproc_pp) then
                    cur_nmicblocks_pp = max_micblock_pp
                  else
                    cur_nmicblocks_pp = max_micblock_pp - 1
                  endif
                  if (n <= max_cpuproc_pp) then
                    cur_ncpublocks_pp = nblocks_avgeff
                  else
                    cur_ncpublocks_pp = nblocks_avgeff - 1
                  endif
               else
                  ! MIC
                    cur_nmicblocks_pp = min(nblocks_mic-nbmap_pos+1,max_micblock_pp)
                  ! CPU
                  cur_ncpublocks_pp = nblocks_avgeff
               endif
            endif ! switch to next proc
         else ! empty/null block
            dist%proc(nblck) = 0
            dist%local_block(nblck) = 0
         endif
      end do ! jblock
      ! check if over-distrb
      if (i_bcount > 0) then
         ncols = ncols - 1
         if (rectshape .and. ncols <= nprocs_cpu - n+1 .and. &
             n<nprocs_cpu) then
               n = n + 1
               local_block = 0
               local_block_cpu = 0
               cur_nmicblocks_counter = 0
               ! block# for current tile: MIC and CPU
               ! MIC
               cur_nmicblocks_pp = min(nblocks_mic-nbmap_pos+1,max_micblock_pp)
               ! CPU
               cur_ncpublocks_pp = nblocks_avgeff
         endif ! switch to next proc
      endif ! check over-distrb
   end do
   end do

   !distrb blocks for MIC
   if (nproc_mic_pn /= 0) then
     if (rectshape) then
        do j=1,nproc_mic_pn*nnode
        do i=1,nblock_mic_pp
           nblck = mic_block_map((j-1)*nblock_mic_pp+i)
           if (nblck /= 0) then ! null task at the tail
           dist%proc(nblck) = nprocs_cpu + j
           dist%local_block(nblck) = i
           endif
        end do
        end do
     endif
     deallocate(mic_block_map)
   endif

   !if (my_task == nproc_cpu_pn*nnode) then ! first MIC process
   if (my_task == master_task) then
      allocate(dist_proc(nblocks_x,nblocks_y))
      dist_proc = reshape(dist%proc, (/nblocks_x,nblocks_y/))
      do j = 1,nblocks_y
         print *,'y ',j, (dist_proc(i,j), i = 1,nblocks_x)
      end do
      deallocate(dist_proc)
   endif
   ! tuning debug
   !call exit_POP(sigAbort,'dist_map completed.')

!----------------------------------------------------------------------

   create_distrb_mic = dist ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_mic

!**********************************************************************
!BOP
! !IROUTINE: create_distrb_balanced
! !INTERFACE:

 function create_distrb_balanced(nprocs, work_per_block)

! !DESCRIPTION:
! This function distributes blocks across processors in a
! load-balanced manner based on the amount of work per block.
! A rake algorithm is used in which the blocks are first distributed
! in a Cartesian distribution and then a rake is applied in each
! Cartesian direction.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_balanced ! resulting structure describing
                              ! load-balanced distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,k,n ,&! dummy loop indices
      pid ,&! dummy for processor id
      local_block ,&! local block position on processor
      max_work ,&! max amount of work in any block
      nprocs_x ,&! num of procs in x for global domain
      nprocs_y ! num of procs in y for global domain

   integer (int_kind), dimension(:), allocatable :: &
      priority ,&! priority for moving blocks
      work_tmp ,&! work per row or column for rake algrthm
      proc_tmp ,&! temp processor id for rake algrthm
      block_count ! counter to determine local block indx

   type (distrb) :: dist ! temp hold distribution

!----------------------------------------------------------------------
!
! first set up as Cartesian distribution
! retain the Cartesian distribution if nblocks_tot = nprocs
! to avoid processors with no work
!
!----------------------------------------------------------------------

   dist = create_distrb_cart(nprocs, work_per_block)
   if (nblocks_tot == nprocs) then
      create_distrb_balanced = dist ! return the result
      return
   endif

!----------------------------------------------------------------------
!
! now re-distribute blocks using a rake in each direction
!
!----------------------------------------------------------------------

   max_work = maxval(work_per_block)

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

!----------------------------------------------------------------------
!
! load-balance using a rake algorithm in the x-direction first
!
!----------------------------------------------------------------------

   allocate(priority(nblocks_tot))

   !*** set highest priority such that eastern-most blocks
   !*** and blocks with the least amount of work are
   !*** moved first

   do j=1,nblocks_y
   do i=1,nblocks_x
      n=(j-1)*nblocks_x + i
      if (work_per_block(n) > 0) then
         priority(n) = (max_work + 1)*(nblocks_x + i) - &
                       work_per_block(n)
      else
         priority(n) = 0
      endif
   end do
   end do

   allocate(work_tmp(nprocs_x), &
            proc_tmp(nprocs_x))

   do j=1,nprocs_y

      work_tmp(:) = 0
      do i=1,nprocs_x
         pid = (j-1)*nprocs_x + i
         proc_tmp(i) = pid
         do n=1,nblocks_tot
            if (dist%proc(n) == pid) then
               work_tmp(i) = work_tmp(i) + work_per_block(n)
            endif
         end do
      end do

      call rake (work_tmp, proc_tmp, work_per_block, priority, dist)

   end do

   deallocate(work_tmp, proc_tmp)

!----------------------------------------------------------------------
!
! use a rake algorithm in the y-direction now
!
!----------------------------------------------------------------------

   !*** set highest priority for northern-most blocks

   do j=1,nblocks_y
   do i=1,nblocks_x
      n=(j-1)*nblocks_x + i
      if (work_per_block(n) > 0) then
         priority(n) = (max_work + 1)*(nblocks_y + j) - &
                       work_per_block(n)
      else
         priority(n) = 0
      endif
   end do
   end do

   allocate(work_tmp(nprocs_y), &
            proc_tmp(nprocs_y))

   do i=1,nprocs_x

      work_tmp(:) = 0
      do j=1,nprocs_y
         pid = (j-1)*nprocs_x + i
         proc_tmp(j) = pid
         do n=1,nblocks_tot
            if (dist%proc(n) == pid) then
               work_tmp(j) = work_tmp(j) + work_per_block(n)
            endif
         end do
      end do

      call rake (work_tmp, proc_tmp, work_per_block, priority, dist)

   end do

   deallocate(work_tmp, proc_tmp)
   deallocate(priority)

!----------------------------------------------------------------------
!
! reset local_block info based on new distribution
!
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   do pid=1,nprocs
      local_block = 0
      do n=1,nblocks_tot
         if (dist%proc(n) == pid) then
            local_block = local_block + 1
            dist%local_block(n) = local_block
            proc_tmp(pid) = proc_tmp(pid) + 1
         endif
      end do
   end do

   if (minval(proc_tmp) < 1) then
      call exit_POP(sigAbort,'Load-balanced distribution failed')
   endif

   deallocate(proc_tmp)

!----------------------------------------------------------------------

   create_distrb_balanced = dist ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_balanced

!**********************************************************************
!BOP
! !IROUTINE: proc_decomposition
! !INTERFACE:

 subroutine proc_decomposition(nprocs, nprocs_x, nprocs_y)

! !DESCRIPTION:
! This subroutine attempts to find an optimal (nearly square)
! 2d processor decomposition for a given number of processors.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nprocs ! total number or processors

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      nprocs_x, nprocs_y ! number of procs in each dimension

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iguess, jguess ! guesses for nproc_x,y

   real (r4) :: &
      square ! square root of nprocs

!----------------------------------------------------------------------
!
! start with an initial guess that is closest to square decomp
!
!----------------------------------------------------------------------

   square = sqrt(real(nprocs))
   nprocs_x = 0
   nprocs_y = 0

   iguess = nint(square)

!----------------------------------------------------------------------
!
! try various decompositions to find the best
!
!----------------------------------------------------------------------

   proc_loop: do
      jguess = nprocs/iguess

      if (iguess*jguess == nprocs) then ! valid decomp

         !***
         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         !***

         if (mod(nblocks_x,iguess) == 0 .and. &
             mod(nblocks_y,jguess) == 0) then
            nprocs_x = iguess
            nprocs_y = jguess
            exit proc_loop

         !***
         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         !***

         else if (mod(nblocks_x,jguess) == 0 .and. &
                mod(nblocks_y,iguess) == 0) then
            nprocs_x = jguess
            nprocs_y = iguess
            exit proc_loop

         !***
         !*** A valid decomposition, but keep searching for
         !*** a better one
         !***

         else
            if (nprocs_x == 0) then
               nprocs_x = iguess
               nprocs_y = jguess
            endif
            iguess = iguess - 1
            if (iguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         endif

      else ! invalid decomp - keep trying

         iguess = iguess - 1
         if (iguess == 0) then
            exit proc_loop
         else
            cycle proc_loop
         endif
      endif
   end do proc_loop

   if (nprocs_x == 0) then
      call exit_POP(sigAbort,'Unable to find 2d processor config')
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine proc_decomposition

!**********************************************************************
!BOP
! !IROUTINE: rake
! !INTERFACE:

 subroutine rake (proc_work, proc_id, block_work, priority, dist)

! !DESCRIPTION:
! This subroutine performs a rake algorithm to distribute the work
! along a vector of processors. In the rake algorithm, a work
! threshold is first set. Then, moving from left to right, work
! above that threshold is raked to the next processor in line.
! The process continues until the end of the vector is reached
! and then the threshold is reduced by one for a second rake pass.
! In this implementation, a priority for moving blocks is defined
! such that the rake algorithm chooses the highest priority
! block to be moved to the next processor. This can be used
! for example to always choose the eastern-most block or to
! ensure a block does not stray too far from its neighbors.
!
! !REVISION HISTORY:
! same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent(inout), dimension(:) :: &
      proc_work ,&! amount of work per processor
      priority ! priority for moving a given block

   integer (int_kind), intent(in), dimension(:) :: &
      block_work ,&! amount of work per block
      proc_id ! global processor number

   type (distrb), intent(inout) :: &
      dist ! distribution to change

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n, m, np1, &
      iproc, inext, &
      nprocs, nblocks, &
      last_priority, last_loc, &
      residual, &
      work_mean, work_max, work_diff, &
      iter, niters, itransfer, ntransfers, &
      min_priority

!----------------------------------------------------------------------
!
! initialization
!
!----------------------------------------------------------------------

   nprocs = size(proc_work)
   nblocks = size(block_work)

   !*** mean work per processor

   work_mean = sum(proc_work)/nprocs + 1
   work_max = maxval(proc_work)
   residual = mod(work_mean,nprocs)

   min_priority = 1000000
   do n=1,nprocs
      iproc = proc_id(n)
      do i=1,nblocks
         if (dist%proc(i) == iproc) then
            min_priority = min(min_priority,priority(i))
         endif
      end do
   end do

!----------------------------------------------------------------------
!
! do two sets of transfers
!
!----------------------------------------------------------------------

   transfer_loop: do

!----------------------------------------------------------------------
!
! do rake across the processors
!
!----------------------------------------------------------------------

      ntransfers = 0
       do n=1,nprocs
          if (n < nprocs) then
             np1 = n+1
          else
             np1 = 1
          endif
         iproc = proc_id(n)
         inext = proc_id(np1)

         if (proc_work(n) > work_mean) then !*** pass work to next
            work_diff = proc_work(n) - work_mean

            rake1: do while (work_diff > 1)

               !*** attempt to find a block with the required
               !*** amount of work and with the highest priority
               !*** for moving (eg boundary blocks first)

               last_priority = 0
               last_loc = 0
               do i=1,nblocks
                  if (dist%proc(i) == iproc) then
                     if (priority(i) > last_priority ) then
                        last_priority = priority(i)
                        last_loc = i
                     endif
                  endif
               end do
               if (last_loc == 0) exit rake1 ! could not shift work

               ntransfers = ntransfers + 1
               dist%proc(last_loc) = inext
               if (np1 == 1) priority(last_loc) = min_priority
               work_diff = work_diff - block_work(last_loc)

               proc_work(n ) = proc_work(n )-block_work(last_loc)
               proc_work(np1) = proc_work(np1)+block_work(last_loc)
            end do rake1
         endif

      end do

!----------------------------------------------------------------------
!
! increment work_mean by one and repeat
!
!----------------------------------------------------------------------

      work_mean = work_mean + 1
      if (ntransfers == 0 .or. work_mean > work_max) exit transfer_loop

   end do transfer_loop

!----------------------------------------------------------------------
!EOC

end subroutine rake

!***********************************************************************

end module distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
