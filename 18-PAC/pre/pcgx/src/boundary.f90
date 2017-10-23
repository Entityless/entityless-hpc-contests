!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: boundary

 module boundary

! !DESCRIPTION:
!  This module contains data types and routines for updating ghost cell
!  boundaries using MPI calls
!
! !REVISION HISTORY:
!  CVS:$Id: boundary.F90,v 1.13 2004/01/07 19:56:32 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use communicate
   use constants
   use blocks
   use distribution
   use exit_mod
   use domain_size
   !use timers

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: bndy
     integer (int_kind) :: &
       communicator       ,&! communicator to use for update messages
       nmsg_ew_snd        ,&! number of messages to send for e-w update
       nmsg_ns_snd        ,&! number of messages to send for n-s update
       nmsg_ew_rcv        ,&! number of messages to recv for e-w update
       nmsg_ns_rcv        ,&! number of messages to recv for n-s update
       maxblocks_ew_snd   ,&! max num blocks involved in east-west sends
       maxblocks_ew_rcv   ,&! max num blocks involved in east-west recvs
       maxblocks_ns_snd   ,&! max num blocks involved in north-south sends
       maxblocks_ns_rcv   ,&! max num blocks involved in north-south recvs
       nlocal_ew,nlocal_ew_tiny ,&! num local copies for east-west bndy update
       nlocal_ns,nlocal_ns_tiny   ! num local copies for east-west bndy update

     integer (int_kind), dimension(:), pointer :: &
       nblocks_ew_snd     ,&! num blocks in each east-west send msg
       nblocks_ns_snd     ,&! num blocks in each north-south send msg
       nblocks_ew_rcv     ,&! num blocks in each east-west recv msg
       nblocks_ns_rcv     ,&! num blocks in each north-south recv msg
       ew_snd_proc        ,&! dest   proc for east-west send message
       ew_rcv_proc        ,&! source proc for east-west recv message
       ns_snd_proc        ,&! dest   proc for north-south send message
       ns_rcv_proc        ,&! source proc for north-south recv message
       local_ew_src_block ,&! source block for each local east-west copy
       local_ew_dst_block ,&! dest   block for each local east-west copy
       local_ns_src_block ,&! source block for each local north-south copy
       local_ns_dst_block   ! dest   block for each local north-south copy

     integer (int_kind), dimension(:,:), pointer :: &
       local_ew_src_add   ,&! starting source address for local e-w copies
       local_ew_dst_add   ,&! starting dest   address for local e-w copies
       local_ns_src_add   ,&! starting source address for local n-s copies
       local_ns_dst_add     ! starting dest   address for local n-s copies

     integer (int_kind), dimension(:,:), pointer :: &
       ew_src_block       ,&! source block for sending   e-w bndy msg
       ew_dst_block       ,&! dest   block for receiving e-w bndy msg
       ns_src_block       ,&! source block for sending   n-s bndy msg
       ns_dst_block         ! dest   block for sending   n-s bndy msg

     integer (int_kind), dimension(:,:,:), pointer :: &
       ew_src_add         ,&! starting source address for e-w msgs
       ew_dst_add         ,&! starting dest   address for e-w msgs
       ns_src_add         ,&! starting source address for n-s msgs
       ns_dst_add           ! starting dest   address for n-s msgs

     integer (int_kind), dimension(:), pointer :: &
       local_ew_src_block_tiny       ,&! source block for sending   e-w bndy msg
       local_ew_dst_block_tiny       ,&! dest   block for receiving e-w bndy msg
       local_ns_src_block_tiny       ,&! source block for sending   n-s bndy msg
       local_ns_dst_block_tiny         ! dest   block for sending   n-s bndy msg

     integer (int_kind), dimension(:,:), pointer :: &
       local_ew_src_add_tiny         ,&! starting source address for e-w msgs
       local_ew_dst_add_tiny         ,&! starting dest   address for e-w msgs
       local_ns_src_add_tiny         ,&! starting source address for n-s msgs
       local_ns_dst_add_tiny           ! starting dest   address for n-s msgs
   end type bndy

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_boundary,  &
             destroy_boundary, &
             update_ghost_cells, &
             update_ghost_cells0, &
             update_ghost_cells1

   interface update_ghost_cells  ! generic interface
     module procedure boundary_2d_dbl,  &
                      boundary_2d_real, &
                      boundary_2d_int,  &
                      boundary_3d_dbl,  &
                      boundary_3d_real, &
                      boundary_3d_int,  &
                      boundary_4d_dbl,  &
                      boundary_4d_real, &
                      boundary_4d_int
   end interface

  interface update_ghost_cells0  ! generic interface
     module procedure boundary_2d_dbl0
   end interface

  interface update_ghost_cells1  ! generic interface
     module procedure boundary_2d_dbl1
   end interface

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  global boundary buffers for tripole boundary
!
!-----------------------------------------------------------------------

   type (distrb) :: &
      loc_dist       ! distribution of blocks across procs

   integer (int_kind), dimension(:,:), allocatable :: &
      tripole_ibuf,  &
      tripole_ighost

   real (r4), dimension(:,:), allocatable :: &
      tripole_rbuf,  &
      tripole_rghost

   real (r8), dimension(:,:), allocatable :: &
      tripole_dbuf,  &
      tripole_dghost

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: create_boundary
! !INTERFACE:

 subroutine create_boundary(newbndy, dist, &
                            ns_bndy_type, ew_bndy_type, &
                            nx_global, ny_global)

! !DESCRIPTION:
!  This routine creates a boundary type with info necessary for
!  performing a boundary (ghost cell) update based on the input block
!  distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      dist       ! distribution of blocks across procs

   character (*), intent(in) :: &
      ns_bndy_type,             &! type of boundary to use in ns dir
      ew_bndy_type               ! type of boundary to use in ew dir

   integer (int_kind), intent(in) :: &
      nx_global, ny_global       ! global extents of domain

! !OUTPUT PARAMETERS:

   type (bndy), intent(out) :: &
      newbndy    ! a new boundary type with info for updates

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,n,ii,iii, ierr,               &! dummy counters
      iblock_src  , jblock_src  ,  &! i,j index of source block
      iblock_dst  , jblock_dst  ,  &! i,j index of dest   block
      iblock_north, jblock_north,  &! i,j index of north neighbor block
      iblock_south, jblock_south,  &! i,j index of south neighbor block
      iblock_east , jblock_east ,  &! i,j index of east  neighbor block
      iblock_west , jblock_west ,  &! i,j index of west  neighbor block
      src_block_loc,               &! local block location of source
      dst_block_loc,               &! local block location of dest
      imsg_ew_snd, imsg_ew_rcv,    &! counters for ew send/recv
      imsg_ns_snd, imsg_ns_rcv,    &! counters for ns send/recv
      nprocs,                      &! num of processors involved
      nblocks,                     &! total number of blocks
      bid, pid,                    &! block and processor locators
      iblk, imsg,                  &!
      iloc_ew, iloc_ns,            &!
      iloc_ew_tiny, iloc_ns_tiny,  &
      src_proc, dst_proc            ! src,dst processor for message

   logical (log_kind) :: &
      lalloc_tripole      ! flag for allocating tripole buffers

   integer (int_kind), dimension(:), allocatable :: &
      ew_snd_count,    &! array for counting blocks in each message
      ew_rcv_count,    &! array for counting blocks in each message
      ns_snd_count,    &! array for counting blocks in each message
      ns_rcv_count,    &! array for counting blocks in each message
      msg_ew_snd  ,    &! msg counter for each active processor
      msg_ew_rcv  ,    &! msg counter for each active processor
      msg_ns_snd  ,    &! msg counter for each active processor
      msg_ns_rcv        ! msg counter for each active processor

   type (block) ::     &
      src_block,       &! block info for source      block
      dst_block         ! block info for destination block
   integer (int_kind), dimension(:), allocatable :: &
      east_boulder_tiny_blocks, &
      west_boulder_tiny_blocks, &
      north_boulder_tiny_blocks, &
      south_boulder_tiny_blocks 
!-----------------------------------------------------------------------
!
!  Initialize some useful variables and return if this task not
!  in the current distribution.
!
!-----------------------------------------------------------------------
   !ljm tuning
   loc_dist = dist

   nprocs = dist%nprocs

   if (my_task >= nprocs) return

   if (.not. lmic_proc) then
   nblocks = size(dist%proc(:))
   else  ! lmic_proc
   nblocks = nblocks_tot
   endif  ! lmic_proc
   lalloc_tripole = .false.
   newbndy%communicator = dist%communicator

!-----------------------------------------------------------------------
!
!  Count the number of messages to send/recv from each processor
!  and number of blocks in each message.  These quantities are
!  necessary for allocating future arrays.
!
!-----------------------------------------------------------------------

   allocate (ew_snd_count(nprocs), ew_rcv_count(nprocs), &
             ns_snd_count(nprocs), ns_rcv_count(nprocs))

   ew_snd_count = 0
   ew_rcv_count = 0
   ns_snd_count = 0
   ns_rcv_count = 0

   block_loop1: do n=1,nblocks
      src_proc  = dist%proc(n)
      src_block = get_block(n,n)

      iblock_src = src_block%iblock  ! i,j index of this block in
      jblock_src = src_block%jblock  !   block cartesian decomposition

      !*** compute cartesian i,j block indices for each neighbor
      !*** use zero if off the end of closed boundary
      !*** use jnorth=nblocks_y and inorth < 0 for tripole boundary
      !***   to make sure top boundary communicated to all top
      !***   boundary blocks

      select case(ew_bndy_type)
      case ('cyclic')
         iblock_east = mod(iblock_src,nblocks_x) + 1
         iblock_west = iblock_src - 1
         if (iblock_west == 0) iblock_west = nblocks_x
         jblock_east = jblock_src
         jblock_west = jblock_src
      case ('closed')
         iblock_east = iblock_src + 1
         iblock_west = iblock_src - 1
         if (iblock_east > nblocks_x) iblock_east = 0
         if (iblock_west < 1        ) iblock_west = 0
         jblock_east = jblock_src
         jblock_west = jblock_src
      case default
         call exit_POP(sigAbort, 'Unknown east-west boundary type')
      end select

      select case(ns_bndy_type)
      case ('cyclic')
         jblock_north = mod(jblock_src,nblocks_y) + 1
         jblock_south = jblock_src - 1
         if (jblock_south == 0) jblock_south = nblocks_y
         iblock_north = iblock_src
         iblock_south = iblock_src
      case ('closed')
         jblock_north = jblock_src + 1
         jblock_south = jblock_src - 1
         if (jblock_north > nblocks_y) jblock_north = 0
         if (jblock_south < 1        ) jblock_south = 0
         iblock_north = iblock_src
         iblock_south = iblock_src
      case ('tripole')
         lalloc_tripole = .true.
         jblock_north = jblock_src + 1
         jblock_south = jblock_src - 1
         iblock_north = iblock_src
         iblock_south = iblock_src
         if (jblock_south < 1        ) jblock_south = 0
         if (jblock_north > nblocks_y) then
            jblock_north = nblocks_y
            iblock_north = -iblock_src
         endif
      case default
         call exit_POP(sigAbort, 'Unknown north-south boundary type')
      end select

      !***
      !*** if any boundary is closed boundary, create a local
      !*** copy pseudo-message to fill ghost cells
      !***

      if (iblock_east == 0) &
         call increment_message_counter(ew_snd_count, ew_rcv_count, &
                                        0, src_proc)
      if (iblock_west == 0) &
         call increment_message_counter(ew_snd_count, ew_rcv_count, &
                                        0, src_proc)
      if (jblock_north == 0) &
         call increment_message_counter(ns_snd_count, ns_rcv_count, &
                                        0, src_proc)
      if (jblock_south == 0) &
         call increment_message_counter(ns_snd_count, ns_rcv_count, &
                                        0, src_proc)

      !***
      !*** now look through all the blocks for the neighbors
      !*** of the source block and check whether a message is
      !*** required for communicating with the neighbor
      !***

      do k=1,nblocks
         dst_block = get_block(k,k)

         iblock_dst = dst_block%iblock  !*** i,j block index of
         jblock_dst = dst_block%jblock  !*** potential neighbor block

         dst_proc = dist%proc(k)  ! processor that holds dst block

         !***
         !*** if this block is an eastern neighbor
         !*** increment message counter
         !***

         if (iblock_dst == iblock_east .and. &
             jblock_dst == jblock_east) then

            call increment_message_counter(ew_snd_count, ew_rcv_count, &
                                           src_proc, dst_proc)
         endif

         !***
         !*** if this block is an western neighbor
         !*** increment message counter
         !***

         if (iblock_dst == iblock_west .and. &
             jblock_dst == jblock_west) then

            call increment_message_counter(ew_snd_count, ew_rcv_count, &
                                           src_proc, dst_proc)
         endif

         !***
         !*** if this block is an northern neighbor
         !*** find out whether a message is required
         !*** for tripole, must communicate with all
         !*** north row blocks (triggered by iblock_dst <0)
         !***

         if ((iblock_dst == iblock_north .or. iblock_north < 0) .and. &
              jblock_dst == jblock_north) then

            call increment_message_counter(ns_snd_count, ns_rcv_count, &
                                           src_proc, dst_proc)
         endif

         !***
         !*** if this block is an southern neighbor
         !*** find out whether a message is required
         !***

         if (iblock_dst == iblock_south .and. &
             jblock_dst == jblock_south) then

            call increment_message_counter(ns_snd_count, ns_rcv_count, &
                                           src_proc, dst_proc)
         endif

      end do  ! search for dest blocks
   end do block_loop1

   !*** if messages are received from the same processor
   !*** the message is actually a local copy - count them
   !*** and reset to zero

   newbndy%nlocal_ew = ew_rcv_count(my_task+1)
   newbndy%nlocal_ns = ns_rcv_count(my_task+1)
   ew_snd_count(my_task+1) = 0
   ew_rcv_count(my_task+1) = 0
   ns_snd_count(my_task+1) = 0
   ns_rcv_count(my_task+1) = 0

   !*** now count the number of actual messages to be
   !*** sent and received

   newbndy%nmsg_ew_snd = count(ew_snd_count /= 0)
   newbndy%nmsg_ns_snd = count(ns_snd_count /= 0)
   newbndy%nmsg_ew_rcv = count(ew_rcv_count /= 0)
   newbndy%nmsg_ns_rcv = count(ns_rcv_count /= 0)

   !*** find the maximum number of blocks sent in any one
   !*** message to use as an array size parameter

   newbndy%maxblocks_ew_snd = maxval(ew_snd_count)
   newbndy%maxblocks_ew_rcv = maxval(ew_rcv_count)
   newbndy%maxblocks_ns_snd = maxval(ns_snd_count)
   newbndy%maxblocks_ns_rcv = maxval(ns_rcv_count)

   !***
   !*** create buffers for tracking which message
   !*** is sent/received from which processor
   !***

   allocate(msg_ew_snd(nprocs), msg_ew_rcv(nprocs), &
            msg_ns_snd(nprocs), msg_ns_rcv(nprocs))

   msg_ew_snd = 0
   msg_ew_rcv = 0
   msg_ns_snd = 0
   msg_ns_rcv = 0

   !***
   !*** assign a location in buffer for each message to a
   !*** different processor. scramble the processor order
   !*** using current task id as offset to prevent all
   !*** processors sending to the same processor at the
   !*** same time
   !***

   imsg_ew_snd = 0
   imsg_ew_rcv = 0
   imsg_ns_snd = 0
   imsg_ns_rcv = 0

   do n=1,nprocs
      dst_proc = modulo(my_task+n,nprocs) + 1
      if (ew_snd_count(dst_proc) /= 0) then
         imsg_ew_snd = imsg_ew_snd + 1
         msg_ew_snd(dst_proc) = imsg_ew_snd
      endif
      if (ew_rcv_count(dst_proc) /= 0) then
         imsg_ew_rcv = imsg_ew_rcv + 1
         msg_ew_rcv(dst_proc) = imsg_ew_rcv
      endif
      if (ns_snd_count(dst_proc) /= 0) then
         imsg_ns_snd = imsg_ns_snd + 1
         msg_ns_snd(dst_proc) = imsg_ns_snd
      endif
      if (ns_rcv_count(dst_proc) /= 0) then
         imsg_ns_rcv = imsg_ns_rcv + 1
         msg_ns_rcv(dst_proc) = imsg_ns_rcv
      endif
   end do

   deallocate(ew_snd_count, ew_rcv_count, ns_snd_count, ns_rcv_count)

!-----------------------------------------------------------------------
!
!  allocate buffers and arrays necessary for boundary comms
!
!-----------------------------------------------------------------------

   allocate (newbndy%nblocks_ew_snd(newbndy%nmsg_ew_snd), &
             newbndy%nblocks_ns_snd(newbndy%nmsg_ns_snd), &
             newbndy%nblocks_ew_rcv(newbndy%nmsg_ew_rcv), &
             newbndy%nblocks_ns_rcv(newbndy%nmsg_ns_rcv))

   allocate (newbndy%local_ew_src_block(newbndy%nlocal_ew), &
             newbndy%local_ew_dst_block(newbndy%nlocal_ew), &
             newbndy%local_ns_src_block(newbndy%nlocal_ns), &
             newbndy%local_ns_dst_block(newbndy%nlocal_ns), &
             newbndy%local_ew_src_add(2,newbndy%nlocal_ew), &
             newbndy%local_ew_dst_add(2,newbndy%nlocal_ew), &
             newbndy%local_ns_src_add(2,newbndy%nlocal_ns), &
             newbndy%local_ns_dst_add(2,newbndy%nlocal_ns))

   allocate ( &
     newbndy%ew_snd_proc (newbndy%nmsg_ew_snd), &
     newbndy%ew_rcv_proc (newbndy%nmsg_ew_rcv), &
     newbndy%ns_snd_proc (newbndy%nmsg_ns_snd), &
     newbndy%ns_rcv_proc (newbndy%nmsg_ns_rcv), &
     newbndy%ew_src_block(newbndy%maxblocks_ew_snd,newbndy%nmsg_ew_snd), &
     newbndy%ew_dst_block(newbndy%maxblocks_ew_rcv,newbndy%nmsg_ew_rcv), &
     newbndy%ns_src_block(newbndy%maxblocks_ns_snd,newbndy%nmsg_ns_snd), &
     newbndy%ns_dst_block(newbndy%maxblocks_ns_rcv,newbndy%nmsg_ns_rcv), &
     newbndy%ew_src_add(2,newbndy%maxblocks_ew_snd,newbndy%nmsg_ew_snd), &
     newbndy%ew_dst_add(2,newbndy%maxblocks_ew_rcv,newbndy%nmsg_ew_rcv), &
     newbndy%ns_src_add(2,newbndy%maxblocks_ns_snd,newbndy%nmsg_ns_snd), &
     newbndy%ns_dst_add(3,newbndy%maxblocks_ns_rcv,newbndy%nmsg_ns_rcv))

   newbndy%nblocks_ew_snd = 0
   newbndy%nblocks_ns_snd = 0
   newbndy%nblocks_ew_rcv = 0
   newbndy%nblocks_ns_rcv = 0

   newbndy%ew_snd_proc = 0
   newbndy%ew_rcv_proc = 0
   newbndy%ns_snd_proc = 0
   newbndy%ns_rcv_proc = 0

   newbndy%local_ew_src_block = 0
   newbndy%local_ew_dst_block = 0
   newbndy%local_ns_src_block = 0
   newbndy%local_ns_dst_block = 0
   newbndy%local_ew_src_add = 0
   newbndy%local_ew_dst_add = 0
   newbndy%local_ns_src_add = 0
   newbndy%local_ns_dst_add = 0

   newbndy%ew_src_block = 0
   newbndy%ew_dst_block = 0
   newbndy%ns_src_block = 0
   newbndy%ns_dst_block = 0
   newbndy%ew_src_add = 0
   newbndy%ew_dst_add = 0
   newbndy%ns_src_add = 0
   newbndy%ns_dst_add = 0

!-----------------------------------------------------------------------
!
!  now set up indices into buffers and address arrays
!
!-----------------------------------------------------------------------

   allocate (ew_snd_count(newbndy%nmsg_ew_snd), &
             ew_rcv_count(newbndy%nmsg_ew_rcv), &
             ns_snd_count(newbndy%nmsg_ns_snd), &
             ns_rcv_count(newbndy%nmsg_ns_rcv))

   ew_snd_count = 0
   ew_rcv_count = 0
   ns_snd_count = 0
   ns_rcv_count = 0

   iloc_ew = 0
   iloc_ns = 0
   
   !call MPI_BARRIER(MPI_COMM_OCN, ierr)
   !call exit_POP(sigAbort,'debug pcg-zero-R...')
   if (lmic_proc) then
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !local copy for tiny blocks which occurs inside the normal block
      ! on mic(not on normal block boundarys), xiaobin
      !North
      !13 >14 >15 >16
      ! 9 >10 >11 >12
      ! 5 > 6 > 7 > 8
      ! 1 > 2 > 3 > 4
      newbndy%nlocal_ew_tiny = 2*nblocks_tot_tiny*(nx_bkfactor-1)/nx_bkfactor &
                               + newbndy%nlocal_ew*ny_bkfactor
      newbndy%nlocal_ns_tiny = 2*nblocks_tot_tiny*(ny_bkfactor-1)/ny_bkfactor &
                               + newbndy%nlocal_ns*nx_bkfactor

      !write(6,*)'mytask+1, nlocal_ew_tiny, nlocal_ns_tiny: ',my_task+1, newbndy%nlocal_ew_tiny, newbndy%nlocal_ns_tiny

      allocate (newbndy%local_ew_src_block_tiny(newbndy%nlocal_ew_tiny), &
                newbndy%local_ew_dst_block_tiny(newbndy%nlocal_ew_tiny), &
                newbndy%local_ns_src_block_tiny(newbndy%nlocal_ns_tiny), &
                newbndy%local_ns_dst_block_tiny(newbndy%nlocal_ns_tiny), &
                newbndy%local_ew_src_add_tiny(2,newbndy%nlocal_ew_tiny), &
                newbndy%local_ew_dst_add_tiny(2,newbndy%nlocal_ew_tiny), &
                newbndy%local_ns_src_add_tiny(2,newbndy%nlocal_ns_tiny), &
                newbndy%local_ns_dst_add_tiny(2,newbndy%nlocal_ns_tiny))
      newbndy%local_ew_src_block_tiny = 0
      newbndy%local_ew_dst_block_tiny = 0
      newbndy%local_ns_src_block_tiny = 0
      newbndy%local_ns_dst_block_tiny = 0
      newbndy%local_ew_src_add_tiny = 0
      newbndy%local_ew_dst_add_tiny = 0
      newbndy%local_ns_src_add_tiny = 0
      newbndy%local_ns_dst_add_tiny = 0
      iloc_ew_tiny = 0
      iloc_ns_tiny = 0

      !??check nblocks_tot_tiny?? right or wrong here?
      !write(6,*)'my_task+1,nblocks_tot_tiny',my_task+1,nblocks_tot_tiny
      do k=1, nblocks_tot_tiny !number of internal tiny block boundarys
         if (mod(k, nx_bkfactor) /= 0) then
            iloc_ew_tiny = iloc_ew_tiny + 2
            newbndy%local_ew_src_block_tiny(iloc_ew_tiny-1) = k
            newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-1) = block_size_x + 1
            newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-1) = 1
            newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-1) = k+1
            newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-1) = 1
            newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-1) = 1

            ! src is dst, dst is scr
            newbndy%local_ew_src_block_tiny(iloc_ew_tiny) = k+1
            newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny) = nghost + 1
            newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny) = 1
            newbndy%local_ew_dst_block_tiny(iloc_ew_tiny) = k
            newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny) = nghost + block_size_x + 1
            newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny) = 1
         end if
         if (mod(k, nx_bkfactor*ny_bkfactor) <= nx_bkfactor*(ny_bkfactor-1) .and.&
             mod(k, nx_bkfactor*ny_bkfactor) /= 0) then
            iloc_ns_tiny = iloc_ns_tiny + 2
            newbndy%local_ns_src_block_tiny(iloc_ns_tiny-1) = k
            newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-1) = 1
            newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-1) = block_size_y + 1
            newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-1) = k + nx_bkfactor
            newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-1) = 1
            newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-1) = 1

            newbndy%local_ns_src_block_tiny(iloc_ns_tiny) = k + nx_bkfactor
            newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny) = 1
            newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny) = nghost + 1
            newbndy%local_ns_dst_block_tiny(iloc_ns_tiny) = k
            newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny) = 1
            newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny) = nghost + block_size_y + 1
         endif
         !write(6,*)'my_task+1,iloc_ew_tiny,iloc_ns_tiny',my_task+1,iloc_ew_tiny,iloc_ns_tiny
      enddo
      !call MPI_BARRIER(MPI_COMM_OCN, ierr)
      !call exit_POP(sigAbort,'debug create_boundary...')

      allocate(east_boulder_tiny_blocks(ny_bkfactor))
      allocate(west_boulder_tiny_blocks(ny_bkfactor))
      allocate(north_boulder_tiny_blocks(nx_bkfactor))
      allocate(south_boulder_tiny_blocks(nx_bkfactor))
      east_boulder_tiny_blocks = 0
      west_boulder_tiny_blocks = 0
      north_boulder_tiny_blocks= 0
      south_boulder_tiny_blocks= 0

      do ii = 1, ny_bkfactor
         east_boulder_tiny_blocks(ii) = nx_bkfactor*ii
         west_boulder_tiny_blocks(ii) = nx_bkfactor*ii-nx_bkfactor+1
      end do

      do ii = 1, nx_bkfactor
         north_boulder_tiny_blocks(ii) = ii + (ny_bkfactor-1)*nx_bkfactor
         south_boulder_tiny_blocks(ii) = ii 
      end do
      !write(6,*)'east_boulder_tiny_blocks: ',east_boulder_tiny_blocks(:)
      !write(6,*)'west_boulder_tiny_blocks: ',west_boulder_tiny_blocks(:)
      !write(6,*)'north_boulder_tiny_blocks: ',north_boulder_tiny_blocks(:)
      !write(6,*)'south_boulder_tiny_blocks: ',south_boulder_tiny_blocks(:)
   end if
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !call MPI_BARRIER(MPI_COMM_OCN, ierr)
   !call exit_POP(sigAbort,'debug create_boundary...')


!write(6,*)'create_bound1,mytask+1,loc_ew,iloc_ns',my_task+1,iloc_ew_tiny,iloc_ns_tiny

!-----------------------------------------------------------------------
!
!  repeat loop through blocks but this time, determine all the
!  required message information for each message or local copy
!
!-----------------------------------------------------------------------

   block_loop2: do n=1,nblocks

      src_proc  = dist%proc(n)    ! processor location for this block
      src_block = get_block(n,n)  ! block info for this block

      iblock_src = src_block%iblock  ! i,j index of this block in
      jblock_src = src_block%jblock  !   block cartesian decomposition

      if (src_proc /= 0) then
         src_block_loc = dist%local_block(n)  ! local block location
      else
         src_block_loc = 0  ! block is a land block
      endif

      !*** compute cartesian i,j block indices for each neighbor
      !*** use zero if off the end of closed boundary
      !*** use jnorth=nblocks_y and inorth < 0 for tripole boundary
      !***   to make sure top boundary communicated to all top
      !***   boundary blocks

      select case(ew_bndy_type)
      case ('cyclic')
         iblock_east = mod(iblock_src,nblocks_x) + 1
         iblock_west = iblock_src - 1
         if (iblock_west == 0) iblock_west = nblocks_x
         jblock_east = jblock_src
         jblock_west = jblock_src
      case ('closed')
         iblock_east = iblock_src + 1
         iblock_west = iblock_src - 1
         if (iblock_east > nblocks_x) iblock_east = 0
         if (iblock_west < 1        ) iblock_west = 0
         jblock_east = jblock_src
         jblock_west = jblock_src
      case default
         call exit_POP(sigAbort, 'Unknown east-west boundary type')
      end select

      select case(ns_bndy_type)
      case ('cyclic')
         jblock_north = mod(jblock_src,nblocks_y) + 1
         jblock_south = jblock_src - 1
         if (jblock_south == 0) jblock_south = nblocks_y
         iblock_north = iblock_src
         iblock_south = iblock_src
      case ('closed')
         jblock_north = jblock_src + 1
         jblock_south = jblock_src - 1
         if (jblock_north > nblocks_y) jblock_north = 0
         if (jblock_south < 1        ) jblock_south = 0
         iblock_north = iblock_src
         iblock_south = iblock_src
      case ('tripole')
         jblock_north = jblock_src + 1
         jblock_south = jblock_src - 1
         iblock_north = iblock_src
         iblock_south = iblock_src
         if (jblock_south < 1        ) jblock_south = 0
         if (jblock_north > nblocks_y) then
            jblock_north = nblocks_y
            iblock_north = -iblock_src
         endif
      case default
         call exit_POP(sigAbort, 'Unknown north-south boundary type')
      end select

      !***
      !*** if blocks are at closed boundary, create local copy
      !*** pseudo-message to fill ghost cells
      !***

      if (src_block_loc /= 0 .and. src_proc == my_task+1) then
         if (iblock_east == 0) then
            iloc_ew = iloc_ew + 1
            newbndy%local_ew_src_block(iloc_ew) = 0
            newbndy%local_ew_src_add(1,iloc_ew) = 0
            newbndy%local_ew_src_add(2,iloc_ew) = 0
            newbndy%local_ew_dst_block(iloc_ew) = src_block_loc
            newbndy%local_ew_dst_add(1,iloc_ew) = src_block%ie + 1
            newbndy%local_ew_dst_add(2,iloc_ew) = 1
            if (lmic_proc) then
               iloc_ew_tiny = iloc_ew_tiny + ny_bkfactor
               newbndy%local_ew_src_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
               newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
               newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
               newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                               (src_block_loc-1)*nx_bkfactor*ny_bkfactor+east_boulder_tiny_blocks
               newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = nghost + block_size_x + 1
               newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
            end if
         else if (iblock_west == 0) then
            iloc_ew = iloc_ew + 1
            newbndy%local_ew_src_block(iloc_ew) = 0
            newbndy%local_ew_src_add(1,iloc_ew) = 0
            newbndy%local_ew_src_add(2,iloc_ew) = 0
            newbndy%local_ew_dst_block(iloc_ew) = src_block_loc
            newbndy%local_ew_dst_add(1,iloc_ew) = 1
            newbndy%local_ew_dst_add(2,iloc_ew) = 1
            if (lmic_proc) then
               iloc_ew_tiny = iloc_ew_tiny + ny_bkfactor
               newbndy%local_ew_src_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
               newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
               newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
               newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+west_boulder_tiny_blocks
               newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
               newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
            end if
         else if (jblock_north == 0) then
            iloc_ns = iloc_ns + 1
            newbndy%local_ns_src_block(iloc_ns) = 0
            newbndy%local_ns_src_add(1,iloc_ns) = 0
            newbndy%local_ns_src_add(2,iloc_ns) = 0
            newbndy%local_ns_dst_block(iloc_ns) = src_block_loc
            newbndy%local_ns_dst_add(1,iloc_ns) = 1
            newbndy%local_ns_dst_add(2,iloc_ns) = src_block%je + 1
            if (lmic_proc) then
               iloc_ns_tiny = iloc_ns_tiny + nx_bkfactor
               newbndy%local_ns_src_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
               newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
               newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
               newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+north_boulder_tiny_blocks
               newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
               newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = nghost + block_size_y + 1
            end if
         else if (jblock_south == 0) then
            iloc_ns = iloc_ns + 1
            newbndy%local_ns_src_block(iloc_ns) = 0
            newbndy%local_ns_src_add(1,iloc_ns) = 0
            newbndy%local_ns_src_add(2,iloc_ns) = 0
            newbndy%local_ns_dst_block(iloc_ns) = src_block_loc
            newbndy%local_ns_dst_add(1,iloc_ns) = 1
            newbndy%local_ns_dst_add(2,iloc_ns) = 1
            if (lmic_proc) then
               iloc_ns_tiny = iloc_ns_tiny + nx_bkfactor
               newbndy%local_ns_src_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
               newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
               newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
               newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks
               newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
               newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
               !write(6,"(A,3I5)")'create_bound1.1,mytask+1,iloc_ns,dst_block',my_task,iloc_ns_tiny,(src_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks
            end if
         endif
      endif

      !***
      !*** now search through blocks looking for neighbors to
      !*** the source block
      !***

      do k=1,nblocks

         dst_proc      = dist%proc(k)  ! processor holding dst block

         !***
         !*** compute the rest only if this block is not a land block
         !***

         if (dst_proc /= 0) then

            dst_block = get_block(k,k)  ! block info for this block

            iblock_dst = dst_block%iblock  ! i,j block index in 
            jblock_dst = dst_block%jblock  ! Cartesian block decomposition

            dst_block_loc = dist%local_block(k)  ! local block location

            !***
            !*** if this block is an eastern neighbor
            !*** determine send/receive addresses
            !***

            if (iblock_dst == iblock_east .and. &
                jblock_dst == jblock_east) then

               if (src_proc == my_task+1 .and. &
                   src_proc == dst_proc) then
                  !*** local copy from one block to another
                  iloc_ew = iloc_ew + 1
                  newbndy%local_ew_src_block(iloc_ew) = src_block_loc
                  newbndy%local_ew_src_add(1,iloc_ew) = src_block%ie - &
                                                     nghost + 1
                  newbndy%local_ew_src_add(2,iloc_ew) = 1
                  newbndy%local_ew_dst_block(iloc_ew) = dst_block_loc
                  newbndy%local_ew_dst_add(1,iloc_ew) = 1
                  newbndy%local_ew_dst_add(2,iloc_ew) = 1

                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ew_tiny = iloc_ew_tiny + ny_bkfactor
                     newbndy%local_ew_src_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+east_boulder_tiny_blocks
                     newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                                block_size_x + 1
                     newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1

                     newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+west_boulder_tiny_blocks
                     newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
                     newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
                  end if

               else if (src_proc == 0 .and. dst_proc == my_task+1) then
                  !*** source block is all land so treat as local copy
                  !*** with source block zero to fill ghost cells with 
                  !*** zeroes
                  iloc_ew = iloc_ew + 1
                  newbndy%local_ew_src_block(iloc_ew) = 0
                  newbndy%local_ew_src_add(1,iloc_ew) = 0
                  newbndy%local_ew_src_add(2,iloc_ew) = 0
                  newbndy%local_ew_dst_block(iloc_ew) = dst_block_loc
                  newbndy%local_ew_dst_add(1,iloc_ew) = 1
                  newbndy%local_ew_dst_add(2,iloc_ew) = 1
                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ew_tiny = iloc_ew_tiny + ny_bkfactor
                     !local copy for normal block boundary
                     newbndy%local_ew_src_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
                     newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
                     newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0

                     newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+west_boulder_tiny_blocks
                     newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
                     newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
                  end if
               else if (src_proc == my_task+1 .and. &
                        dst_proc /= my_task+1) then
                  !*** an actual message must be sent
                  imsg = msg_ew_snd(dst_proc)
                  ew_snd_count(imsg) = ew_snd_count(imsg) + 1
                  iblk = ew_snd_count(imsg)
                  newbndy%ew_snd_proc (     imsg) = dst_proc
                  newbndy%ew_src_block(iblk,imsg) = src_block_loc
                  newbndy%ew_src_add(1,iblk,imsg) = src_block%ie - &
                                                    nghost + 1
                  newbndy%ew_src_add(2,iblk,imsg) = 1
                  newbndy%nblocks_ew_snd(imsg) = &
                  newbndy%nblocks_ew_snd(imsg) + 1
               else if (dst_proc == my_task+1 .and. &
                        src_proc /= my_task+1) then
                  !*** must receive a message
                  imsg = msg_ew_rcv(src_proc)
                  ew_rcv_count(imsg) = ew_rcv_count(imsg) + 1
                  iblk = ew_rcv_count(imsg)
                  newbndy%ew_rcv_proc (     imsg) = src_proc
                  newbndy%ew_dst_block(iblk,imsg) = dst_block_loc
                  newbndy%ew_dst_add(:,iblk,imsg) = 1
                  newbndy%nblocks_ew_rcv(imsg) = &
                  newbndy%nblocks_ew_rcv(imsg) + 1
               endif

            endif ! east neighbor

            !***
            !*** if this block is a western neighbor
            !*** determine send/receive addresses
            !***

            if (iblock_dst == iblock_west .and. &
                jblock_dst == jblock_west) then

               if (src_proc == my_task+1 .and. &
                   src_proc == dst_proc) then
                  !*** perform a local copy
                  iloc_ew = iloc_ew + 1
                  newbndy%local_ew_src_block(iloc_ew) = src_block_loc
                  newbndy%local_ew_src_add(1,iloc_ew) = nghost + 1
                  newbndy%local_ew_src_add(2,iloc_ew) = 1
                  newbndy%local_ew_dst_block(iloc_ew) = dst_block_loc
                  newbndy%local_ew_dst_add(1,iloc_ew) = dst_block%ie + 1
                  newbndy%local_ew_dst_add(2,iloc_ew) = 1

                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ew_tiny = iloc_ew_tiny + ny_bkfactor
                     newbndy%local_ew_src_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+west_boulder_tiny_blocks
                     newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                                nghost + 1
                     newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1

                     newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+east_boulder_tiny_blocks
                     newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                                nghost + block_size_x + 1
                     newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
                  end if

               else if (src_proc == 0 .and. dst_proc == my_task+1) then
                  !*** neighbor is a land block so zero ghost cells
                  iloc_ew = iloc_ew + 1
                  newbndy%local_ew_src_block(iloc_ew) = 0
                  newbndy%local_ew_src_add(1,iloc_ew) = 0
                  newbndy%local_ew_src_add(2,iloc_ew) = 0
                  newbndy%local_ew_dst_block(iloc_ew) = dst_block_loc
                  newbndy%local_ew_dst_add(1,iloc_ew) = dst_block%ie + 1
                  newbndy%local_ew_dst_add(2,iloc_ew) = 1
                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ew_tiny = iloc_ew_tiny + ny_bkfactor
                     !local copy for normal block boundary
                     newbndy%local_ew_src_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
                     newbndy%local_ew_src_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0
                     newbndy%local_ew_src_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 0

                     newbndy%local_ew_dst_block_tiny(iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+east_boulder_tiny_blocks
                     newbndy%local_ew_dst_add_tiny(1,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = &
                                                                nghost + block_size_x + 1
                     newbndy%local_ew_dst_add_tiny(2,iloc_ew_tiny-ny_bkfactor+1 : iloc_ew_tiny) = 1
                  end if

               else if (src_proc == my_task+1 .and. &
                        dst_proc /= my_task+1) then
                  !*** message must be sent
                  imsg = msg_ew_snd(dst_proc)
                  ew_snd_count(imsg) = ew_snd_count(imsg) + 1
                  iblk = ew_snd_count(imsg)
                  newbndy%ew_snd_proc (     imsg) = dst_proc
                  newbndy%ew_src_block(iblk,imsg) = src_block_loc
                  newbndy%ew_src_add(1,iblk,imsg) = nghost + 1
                  newbndy%ew_src_add(2,iblk,imsg) = 1
                  newbndy%nblocks_ew_snd(imsg) = &
                  newbndy%nblocks_ew_snd(imsg) + 1
               else if (dst_proc == my_task+1 .and. &
                        src_proc /= my_task+1) then
                  !*** message must be received
                  imsg = msg_ew_rcv(src_proc)
                  ew_rcv_count(imsg) = ew_rcv_count(imsg) + 1
                  iblk = ew_rcv_count(imsg)
                  newbndy%ew_rcv_proc (     imsg) = src_proc
                  newbndy%ew_dst_block(iblk,imsg) = dst_block_loc
                  newbndy%ew_dst_add(1,iblk,imsg) = dst_block%ie + 1
                  newbndy%ew_dst_add(2,iblk,imsg) = 1
                  newbndy%nblocks_ew_rcv(imsg) = &
                  newbndy%nblocks_ew_rcv(imsg) + 1
               endif

            endif ! west neighbor

            !***
            !*** if this block is a northern neighbor
            !***  compute send/recv addresses
            !*** for tripole, must communicate with all
            !*** north row blocks (triggered by iblock_dst <0)
            !***

            if ((iblock_dst == iblock_north .or. iblock_north < 0) .and. &
                 jblock_dst == jblock_north) then

               if (src_proc == my_task+1 .and. &
                   src_proc == dst_proc) then
                  !*** local copy
                  iloc_ns = iloc_ns + 1
                  newbndy%local_ns_src_block(iloc_ns) = src_block_loc
                  newbndy%local_ns_src_add(1,iloc_ns) = 1
                  newbndy%local_ns_src_add(2,iloc_ns) = src_block%je - &
                                                        nghost + 1
                  newbndy%local_ns_dst_block(iloc_ns) = dst_block_loc
                  newbndy%local_ns_dst_add(1,iloc_ns) = 1
                  newbndy%local_ns_dst_add(2,iloc_ns) = 1

                  if (iblock_north < 0) then !*** tripole boundary

                     newbndy%local_ns_dst_block(iloc_ns) = -dst_block_loc
                     !*** copy nghost+1 northern rows of physical
                     !*** domain into global north tripole buffer
                     newbndy%local_ns_src_add(1,iloc_ns) = &
                                        src_block%i_glob(nghost+1)
                     newbndy%local_ns_src_add(2,iloc_ns) = &
                                        dst_block%je - nghost

                     !*** copy out of tripole ghost cell buffer
                     !*** over-write the last row of the destination
                     !*** block to enforce for symmetry for fields
                     !*** located on domain boundary
                     newbndy%local_ns_dst_add(1,iloc_ns) = &
                                          dst_block%i_glob(nghost+1)
                     newbndy%local_ns_dst_add(2,iloc_ns) = & 
                                          dst_block%je
                  endif

                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ns_tiny = iloc_ns_tiny + nx_bkfactor
                     newbndy%local_ns_src_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+north_boulder_tiny_blocks
                     newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
                     newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = block_size_y + 1

                     newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks
                     newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1

                     newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1

                     if (iblock_north < 0) then !*** tripole boundary

                        newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              -1*((dst_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks)
                        !*** copy nghost+1 northern rows of physical
                        !*** domain into global north tripole buffer
                        iii=0
                        do ii = iloc_ns_tiny-nx_bkfactor+1 , iloc_ns_tiny
                           newbndy%local_ns_src_add_tiny(1, ii) = &
                                              src_block%i_glob(nghost+1)+iii*block_size_x
                           newbndy%local_ns_src_add_tiny(2, ii) = &
                                              block_size_y
                          ! newbndy%local_ns_src_add(2, ii) = &
                          !                    dst_block%je - nghost
                           !*** copy out of tripole ghost cell buffer
                           !*** over-write the last row of the destination
                           !*** block to enforce for symmetry for fields
                           !*** located on domain boundary
                           newbndy%local_ns_dst_add_tiny(1, ii) = &
                                                dst_block%i_glob(nghost+1)+iii*block_size_x
                           newbndy%local_ns_dst_add_tiny(2, ii) = &
                                                block_size_y+2*nghost
                           iii = iii+1
                        enddo
                     endif
                  endif

               else if (src_proc == 0 .and. dst_proc == my_task+1) then
                  !*** source is land block so zero ghost cells
                  iloc_ns = iloc_ns + 1
                  newbndy%local_ns_src_block(iloc_ns) = 0
                  newbndy%local_ns_src_add(1,iloc_ns) = 0
                  newbndy%local_ns_src_add(2,iloc_ns) = 0
                  newbndy%local_ns_dst_block(iloc_ns) = dst_block_loc
                  newbndy%local_ns_dst_add(1,iloc_ns) = 1
                  newbndy%local_ns_dst_add(2,iloc_ns) = 1
                  if (iblock_north < 0) then !*** tripole boundary
                     newbndy%local_ns_dst_block(iloc_ns) = -dst_block_loc
                     !*** replace i addresses with global i location
                     !*** for copies into and out of global buffer
                     newbndy%local_ns_dst_add(1,iloc_ns) = &
                                             dst_block%i_glob(nghost+1)
                     newbndy%local_ns_dst_add(2,iloc_ns) = dst_block%je
                  endif

                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ns_tiny = iloc_ns_tiny + nx_bkfactor
                     newbndy%local_ns_src_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
                     newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
                     newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0

                     newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks
                     newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1

                     newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
                  !write(6,"(A,12I5)")'xiaobin createbound1.21,task,iloc_ns,dst_block,dst_block_loc-1',my_task,iloc_ns_tiny,newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny),dst_block_loc-1,nx_bkfactor*ny_bkfactor,south_boulder_tiny_blocks(:)

                     if (iblock_north < 0) then !*** tripole boundary
                        newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              -1*((dst_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks)
                        !*** replace i addresses with global i location
                        !*** for copies into and out of global buffer
                        iii = 0
                        do ii = iloc_ns_tiny-nx_bkfactor+1 , iloc_ns_tiny
                           newbndy%local_ns_dst_add_tiny(1,ii) = &
                                                   dst_block%i_glob(nghost+1)+iii*block_size_x
                           newbndy%local_ns_dst_add_tiny(2,ii) = block_size_y+nghost
                           !newbndy%local_ns_dst_add_tiny(2,ii) = dst_block%je
                           iii = iii+1
                        end do
                     endif
                  end if
               else if (src_proc == my_task+1 .and. &
                        dst_proc /= my_task+1) then
                  !*** message must be sent
                  imsg = msg_ns_snd(dst_proc)
                  ns_snd_count(imsg) = ns_snd_count(imsg) + 1
                  iblk = ns_snd_count(imsg)
                  newbndy%ns_snd_proc (     imsg) = dst_proc
                  newbndy%ns_src_block(iblk,imsg) = src_block_loc
                  newbndy%ns_src_add(1,iblk,imsg) = 1
                  newbndy%ns_src_add(2,iblk,imsg) = src_block%je - &
                                                    nghost + 1
                  newbndy%nblocks_ns_snd(imsg) = &
                  newbndy%nblocks_ns_snd(imsg) + 1
                  if (iblock_north < 0) then !*** tripole boundary
                     !*** need extra ghost cell for U points
                     newbndy%ns_src_add(2,iblk,imsg) = src_block%je - nghost
                  endif
               else if (dst_proc == my_task+1 .and. &
                        src_proc /= my_task+1) then
                  !*** message must be received
                  imsg = msg_ns_rcv(src_proc)
                  ns_rcv_count(imsg) = ns_rcv_count(imsg) + 1
                  iblk = ns_rcv_count(imsg)
                  newbndy%ns_rcv_proc (     imsg) = src_proc
                  newbndy%ns_dst_block(iblk,imsg) = dst_block_loc
                  newbndy%ns_dst_add(1,iblk,imsg) = 1
                  newbndy%ns_dst_add(2,iblk,imsg) = 1
                  newbndy%nblocks_ns_rcv(imsg) = &
                  newbndy%nblocks_ns_rcv(imsg) + 1
                  if (iblock_north < 0) then !*** tripole
                     newbndy%ns_dst_block(iblk,imsg) = -dst_block_loc
                     !*** upon receiving message, store in global 
                     !*** tripole buffer for src, then copy out of
                     !*** ghost cell buffer once global buffer filled
                     !*** i address for storing in global buffer
                     newbndy%ns_dst_add(1,iblk,imsg) = &
                                            src_block%i_glob(nghost+1)
                     !*** addresses for copying out of ghost buffer
                     newbndy%ns_dst_add(2,iblk,imsg) = dst_block%je
                     newbndy%ns_dst_add(3,iblk,imsg) = &
                                       dst_block%i_glob(nghost+1)
                  endif
               endif

            endif ! north neighbor

            !***
            !*** if this block is a southern neighbor
            !*** determine send/receive addresses
            !***

            if (iblock_dst == iblock_south .and. &
                jblock_dst == jblock_south) then

               if (src_proc == my_task+1 .and. &
                   src_proc == dst_proc) then
                  !*** local copy
                  iloc_ns = iloc_ns + 1
                  newbndy%local_ns_src_block(iloc_ns) = src_block_loc
                  newbndy%local_ns_src_add(1,iloc_ns) = 1
                  newbndy%local_ns_src_add(2,iloc_ns) = nghost + 1
                  newbndy%local_ns_dst_block(iloc_ns) = dst_block_loc
                  newbndy%local_ns_dst_add(1,iloc_ns) = 1
                  newbndy%local_ns_dst_add(2,iloc_ns) = dst_block%je + 1
                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ns_tiny = iloc_ns_tiny + nx_bkfactor
                     newbndy%local_ns_src_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (src_block_loc-1)*nx_bkfactor*ny_bkfactor+south_boulder_tiny_blocks
                     newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
                     newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = nghost + 1

                     newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+north_boulder_tiny_blocks
                     newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
                     newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                              nghost + block_size_y + 1

                  end if
               else if (src_proc == 0 .and. dst_proc == my_task+1) then
                  !*** neighbor is a land block so zero ghost cells
                  iloc_ns = iloc_ns + 1
                  newbndy%local_ns_src_block(iloc_ns) = 0
                  newbndy%local_ns_src_add(1,iloc_ns) = 0
                  newbndy%local_ns_src_add(2,iloc_ns) = 0
                  newbndy%local_ns_dst_block(iloc_ns) = dst_block_loc
                  newbndy%local_ns_dst_add(1,iloc_ns) = 1
                  newbndy%local_ns_dst_add(2,iloc_ns) = dst_block%je + 1
                  ! local copy for tiny block xiaobin
                  if (lmic_proc) then
                     iloc_ns_tiny = iloc_ns_tiny + nx_bkfactor
                     newbndy%local_ns_src_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
                     newbndy%local_ns_src_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0
                     newbndy%local_ns_src_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 0

                     newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                                (dst_block_loc-1)*nx_bkfactor*ny_bkfactor+north_boulder_tiny_blocks
                     newbndy%local_ns_dst_add_tiny(1,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = 1
                     newbndy%local_ns_dst_add_tiny(2,iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny) = &
                                                                nghost + block_size_y + 1
                  !write(6,"(A,6I5)")'xiaobin createbound1.3,task,iloc_ns,dst_block',my_task,iloc_ns_tiny,newbndy%local_ns_dst_block_tiny(iloc_ns_tiny-nx_bkfactor+1 : iloc_ns_tiny)
                  end if
               else if (src_proc == my_task+1 .and. &
                        dst_proc /= my_task+1) then
                  !*** message must be sent
                  imsg = msg_ns_snd(dst_proc)
                  ns_snd_count(imsg) = ns_snd_count(imsg) + 1
                  iblk = ns_snd_count(imsg)
                  newbndy%ns_snd_proc (     imsg) = dst_proc
                  newbndy%ns_src_block(iblk,imsg) = src_block_loc
                  newbndy%ns_src_add(1,iblk,imsg) = 1
                  newbndy%ns_src_add(2,iblk,imsg) = nghost+1
                  newbndy%nblocks_ns_snd(imsg) = &
                  newbndy%nblocks_ns_snd(imsg) + 1
               else if (dst_proc == my_task+1 .and. &
                        src_proc /= my_task+1) then
                  !*** message must be received
                  imsg = msg_ns_rcv(src_proc)
                  ns_rcv_count(imsg) = ns_rcv_count(imsg) + 1
                  iblk = ns_rcv_count(imsg)
                  newbndy%ns_rcv_proc (     imsg) = src_proc
                  newbndy%ns_dst_block(iblk,imsg) = dst_block_loc
                  newbndy%ns_dst_add(1,iblk,imsg) = 1
                  newbndy%ns_dst_add(2,iblk,imsg) = dst_block%je + 1
                  newbndy%nblocks_ns_rcv(imsg) = &
                  newbndy%nblocks_ns_rcv(imsg) + 1
               endif
            endif ! south neighbor

         endif  ! not a land block

      end do
   end do block_loop2

   !write(6,*)'create_bound2,mytask+1,iloc_ew,iloc_ns',my_task+1, iloc_ew_tiny, iloc_ns_tiny
   !if (lmic_proc) then
   !   do i =1, newbndy%nlocal_ew_tiny
   !      write(6,"(A,7I5)")'task,ew_src_block,add1,2,dst_block,add1,2',my_task,newbndy%local_ew_src_block_tiny(i),newbndy%local_ew_src_add_tiny(1,i),newbndy%local_ew_src_add_tiny(2,i),newbndy%local_ew_dst_block_tiny(i),newbndy%local_ew_dst_add_tiny(1,i),newbndy%local_ew_dst_add_tiny(2,i)
   !   enddo
   !   do i =1, newbndy%nlocal_ns_tiny
   !      write(6,"(A,7I5)")'task,ns_src_block,add1,2,dst_block,add1,2',my_task,newbndy%local_ns_src_block_tiny(i),newbndy%local_ns_src_add_tiny(1,i),newbndy%local_ns_src_add_tiny(2,i),newbndy%local_ns_dst_block_tiny(i),newbndy%local_ns_dst_add_tiny(1,i),newbndy%local_ns_dst_add_tiny(2,i)
   !   enddo
   !endif
!-----------------------------------------------------------------------

   deallocate(ew_snd_count, ew_rcv_count, ns_snd_count, ns_rcv_count)
   deallocate(msg_ew_snd, msg_ew_rcv, msg_ns_snd, msg_ns_rcv)

!-----------------------------------------------------------------------
!
!  if necessary, create tripole boundary buffers for each
!  common data type.  the ghost cell buffer includes an
!  extra row for the last physical row in order to enforce
!  symmetry conditions on variables at U points.  the other buffer
!  contains an extra row for handling y-offset for north face or
!  northeast corner points.
!
!-----------------------------------------------------------------------

   if (lalloc_tripole .and. .not. allocated(tripole_ibuf)) then
      allocate(tripole_ibuf  (nx_global,nghost+1), &
               tripole_rbuf  (nx_global,nghost+1), &
               tripole_dbuf  (nx_global,nghost+1), &
               tripole_ighost(nx_global,nghost+1), &
               tripole_rghost(nx_global,nghost+1), &
               tripole_dghost(nx_global,nghost+1))
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine create_boundary

!***********************************************************************
!BOP
! !IROUTINE: destroy_boundary
! !INTERFACE:

 subroutine destroy_boundary(in_bndy)

! !DESCRIPTION:
!  This routine destroys a boundary by deallocating all memory
!  associated with the boundary and nullifying pointers.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
     in_bndy          ! boundary structure to be destroyed

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  reset all scalars
!
!-----------------------------------------------------------------------

   in_bndy%communicator      = 0
   in_bndy%nmsg_ew_snd       = 0
   in_bndy%nmsg_ns_snd       = 0
   in_bndy%nmsg_ew_rcv       = 0
   in_bndy%nmsg_ns_rcv       = 0
   in_bndy%maxblocks_ew_snd  = 0
   in_bndy%maxblocks_ew_rcv  = 0
   in_bndy%maxblocks_ns_snd  = 0
   in_bndy%maxblocks_ns_rcv  = 0
   in_bndy%nlocal_ew         = 0
   in_bndy%nlocal_ns         = 0

!-----------------------------------------------------------------------
!
!  deallocate all pointers
!
!-----------------------------------------------------------------------

   deallocate(in_bndy%nblocks_ew_snd,     &
              in_bndy%nblocks_ns_snd,     &
              in_bndy%nblocks_ew_rcv,     &
              in_bndy%nblocks_ns_rcv,     &
              in_bndy%ew_snd_proc,        &
              in_bndy%ew_rcv_proc,        &
              in_bndy%ns_snd_proc,        &
              in_bndy%ns_rcv_proc,        &
              in_bndy%local_ew_src_block, &
              in_bndy%local_ew_dst_block, &
              in_bndy%local_ns_src_block, &
              in_bndy%local_ns_dst_block, &
              in_bndy%local_ew_src_add,   &
              in_bndy%local_ew_dst_add,   &
              in_bndy%local_ns_src_add,   &
              in_bndy%local_ns_dst_add,   &
              in_bndy%ew_src_block,       &
              in_bndy%ew_dst_block,       &
              in_bndy%ns_src_block,       &
              in_bndy%ns_dst_block,       &
              in_bndy%ew_src_add,         &
              in_bndy%ew_dst_add,         &
              in_bndy%ns_src_add,         &
              in_bndy%ns_dst_add )

!-----------------------------------------------------------------------
!EOC

 end subroutine destroy_boundary

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

 subroutine boundary_2d_dbl(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) ::     &
      this_block        ! block information for current block

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r8), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   real (r8) :: &
      xavg               ! scalar for enforcing symmetry at U pts

   !logical (log_kind), save :: first_call = .true.
   !integer (int_kind), save :: bndy_2d_local, bndy_2d_recv, &
   !                            bndy_2d_send, bndy_2d_wait, bndy_2d_final

   if (.not. lmic_proc) then
!   !ljm tuning
!   if (lmic_trace) then
!   do j=1,loc_dist%local_block_num
!      i = loc_dist%local_block_ids(j)
!      if (i>nblocks_tot) then ! tiny-block
!        this_block = get_block(i,j) 
!        write(6,*) 'pcg-zero:RHS:',this_block%block_id,i-nblocks_tot,sum(ARRAY(:,:,j))
!      else
!        write(6,*) 'pcg-zero:RHS:',i,0,sum(ARRAY(:,:,j))
!      endif
!   enddo
!   call MPI_BARRIER(MPI_COMM_OCN, ierr)
!   call exit_POP(sigAbort,'debug pcg-zero-R...')
!   endif

!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   !if (first_call) then
   !  first_call = .false.
   !  call get_timer(bndy_2d_local, 'BNDY_2D_LOCAL')
   !  call get_timer(bndy_2d_recv,  'BNDY_2D_RECV')
   !  call get_timer(bndy_2d_send,  'BNDY_2D_SEND')
   !  call get_timer(bndy_2d_wait,  'BNDY_2D_WAIT')
   !  call get_timer(bndy_2d_final, 'BNDY_2D_FINAL')
   !endif

   allocate(buf_ew_snd(nghost, ny_block, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_dbuf)) nx_global = size(tripole_dbuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ew_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n)
         buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:buff_ew_snd:',my_task,n,sum(buf_ew_snd(:,:,:,n))
      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_local)
   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = c0
      endif
   end do
!   if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:update_local',my_task,sum(ARRAY(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------


   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr) !8.7
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:buff_ns_snd:',my_task,n,sum(buf_ns_snd(:,:,:,n))
      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) tripole_dbuf = c0

   !call timer_start(bndy_2d_local)
   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = c0
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_dghost(nx_global,1) = tripole_dbuf(nx_global,nghost+1)
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_dghost(i,1+j) = isign* &
                                 tripole_dbuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(i,:)
               end do
            endif

            !*** now copy the remaining ghost cell values

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(i,:)
               end do
            endif

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

   else   ! lmic_proc
      call boundary_2d_dbl_tiny(ARRAY, in_bndy, grid_loc, field_type)
!xiaobin
!      do i=1,nblocks_tot_tiny
!         write(6,"(A,2I5,4F10.3)")'task,nlocal_tiny,4boulders',my_task,i,sum(ARRAY(1:block_size_x+2*nghost,1,i)),&
!         sum(ARRAY(1:block_size_x+2*nghost,block_size_y+2*nghost,i)),sum(ARRAY(1,1:block_size_y+2*nghost,i)),&
!         sum(ARRAY(block_size_x+2*nghost,1:block_size_y+2*nghost,i))
!      enddo
   endif  ! lmic_proc
!-----------------------------------------------------------------------
!xiaobintest
 !call MPI_BARRIER(MPI_COMM_OCN, ierr)
 !call exit_POP(sigAbort, 'xiaobin test boundary2ddbl')
 end subroutine boundary_2d_dbl

  subroutine boundary_2d_dbl0(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) ::     &
      this_block        ! block information for current block

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r8), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   real (r8) :: &
      xavg               ! scalar for enforcing symmetry at U pts

   !logical (log_kind), save :: first_call = .true.
   !integer (int_kind), save :: bndy_2d_local, bndy_2d_recv, &
   !                            bndy_2d_send, bndy_2d_wait, bndy_2d_final

   if (.not. lmic_proc) then
!   !ljm tuning
!   if (lmic_trace) then
!   do j=1,loc_dist%local_block_num
!      i = loc_dist%local_block_ids(j)
!      if (i>nblocks_tot) then ! tiny-block
!        this_block = get_block(i,j) 
!        write(6,*) 'pcg-zero:RHS:',this_block%block_id,i-nblocks_tot,sum(ARRAY(:,:,j))
!      else
!        write(6,*) 'pcg-zero:RHS:',i,0,sum(ARRAY(:,:,j))
!      endif
!   enddo
!   call MPI_BARRIER(MPI_COMM_OCN, ierr)
!   call exit_POP(sigAbort,'debug pcg-zero-R...')
!   endif

!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   !if (first_call) then
   !  first_call = .false.
   !  call get_timer(bndy_2d_local, 'BNDY_2D_LOCAL')
   !  call get_timer(bndy_2d_recv,  'BNDY_2D_RECV')
   !  call get_timer(bndy_2d_send,  'BNDY_2D_SEND')
   !  call get_timer(bndy_2d_wait,  'BNDY_2D_WAIT')
   !  call get_timer(bndy_2d_final, 'BNDY_2D_FINAL')
   !endif

   allocate(buf_ew_snd(nghost, ny_block, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_dbuf)) nx_global = size(tripole_dbuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ew_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n)
         buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:buff_ew_snd:',my_task,n,sum(buf_ew_snd(:,:,:,n))
      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_local)
   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = c0
      endif
   end do
!   if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:update_local',my_task,sum(ARRAY(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------


   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:buff_ns_snd:',my_task,n,sum(buf_ns_snd(:,:,:,n))
      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) tripole_dbuf = c0

   !call timer_start(bndy_2d_local)
   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = c0
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_dghost(nx_global,1) = tripole_dbuf(nx_global,nghost+1)
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_dghost(i,1+j) = isign* &
                                 tripole_dbuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(i,:)
               end do
            endif

            !*** now copy the remaining ghost cell values

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(i,:)
               end do
            endif

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

   else   ! lmic_proc
      call boundary_2d_dbl_tiny(ARRAY, in_bndy, grid_loc, field_type)
!xiaobin
!      do i=1,nblocks_tot_tiny
!         write(6,"(A,2I5,4F10.3)")'task,nlocal_tiny,4boulders',my_task,i,sum(ARRAY(1:block_size_x+2*nghost,1,i)),&
!         sum(ARRAY(1:block_size_x+2*nghost,block_size_y+2*nghost,i)),sum(ARRAY(1,1:block_size_y+2*nghost,i)),&
!         sum(ARRAY(block_size_x+2*nghost,1:block_size_y+2*nghost,i))
!      enddo
   endif  ! lmic_proc
!-----------------------------------------------------------------------
!xiaobintest
 !call MPI_BARRIER(MPI_COMM_OCN, ierr)
 !call exit_POP(sigAbort, 'xiaobin test boundary2ddbl')
 end subroutine boundary_2d_dbl0

  subroutine boundary_2d_dbl1(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) ::     &
      this_block        ! block information for current block

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r8), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   real (r8) :: &
      xavg               ! scalar for enforcing symmetry at U pts

   !logical (log_kind), save :: first_call = .true.
   !integer (int_kind), save :: bndy_2d_local, bndy_2d_recv, &
   !                            bndy_2d_send, bndy_2d_wait, bndy_2d_final

   if (.not. lmic_proc) then
!   !ljm tuning
!   if (lmic_trace) then
!   do j=1,loc_dist%local_block_num
!      i = loc_dist%local_block_ids(j)
!      if (i>nblocks_tot) then ! tiny-block
!        this_block = get_block(i,j) 
!        write(6,*) 'pcg-zero:RHS:',this_block%block_id,i-nblocks_tot,sum(ARRAY(:,:,j))
!      else
!        write(6,*) 'pcg-zero:RHS:',i,0,sum(ARRAY(:,:,j))
!      endif
!   enddo
!   call MPI_BARRIER(MPI_COMM_OCN, ierr)
!   call exit_POP(sigAbort,'debug pcg-zero-R...')
!   endif

!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   !if (first_call) then
   !  first_call = .false.
   !  call get_timer(bndy_2d_local, 'BNDY_2D_LOCAL')
   !  call get_timer(bndy_2d_recv,  'BNDY_2D_RECV')
   !  call get_timer(bndy_2d_send,  'BNDY_2D_SEND')
   !  call get_timer(bndy_2d_wait,  'BNDY_2D_WAIT')
   !  call get_timer(bndy_2d_final, 'BNDY_2D_FINAL')
   !endif

   allocate(buf_ew_snd(nghost, ny_block, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_dbuf)) nx_global = size(tripole_dbuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ew_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n)
         buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:buff_ew_snd:',my_task,n,sum(buf_ew_snd(:,:,:,n))
      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_local)
   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = c0
      endif
   end do
!   if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:update_local',my_task,sum(ARRAY(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------


   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl:buff_ns_snd:',my_task,n,sum(buf_ns_snd(:,:,:,n))
      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) tripole_dbuf = c0

   !call timer_start(bndy_2d_local)
   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = c0
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_dghost(nx_global,1) = tripole_dbuf(nx_global,nghost+1)
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_dghost(i,1+j) = isign* &
                                 tripole_dbuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(i,:)
               end do
            endif

            !*** now copy the remaining ghost cell values

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_dghost(i,:)
               end do
            endif

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

   else   ! lmic_proc
      call boundary_2d_dbl_tiny(ARRAY, in_bndy, grid_loc, field_type)
!xiaobin
!      do i=1,nblocks_tot_tiny
!         write(6,"(A,2I5,4F10.3)")'task,nlocal_tiny,4boulders',my_task,i,sum(ARRAY(1:block_size_x+2*nghost,1,i)),&
!         sum(ARRAY(1:block_size_x+2*nghost,block_size_y+2*nghost,i)),sum(ARRAY(1,1:block_size_y+2*nghost,i)),&
!         sum(ARRAY(block_size_x+2*nghost,1:block_size_y+2*nghost,i))
!      enddo
   endif  ! lmic_proc
!-----------------------------------------------------------------------
!xiaobintest
 !call MPI_BARRIER(MPI_COMM_OCN, ierr)
 !call exit_POP(sigAbort, 'xiaobin test boundary2ddbl')
 end subroutine boundary_2d_dbl1

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

 subroutine boundary_2d_dbl_tiny(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) ::     &
      this_block        ! block information for current block

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      ib_src_tiny,ie_src_tiny,jb_src_tiny,je_src_tiny, &!
      ib_dst_tiny,ie_dst_tiny,jb_dst_tiny,je_dst_tiny, &!
      iblock,jblock,iblock2,jblock2,&!
      src_tiny_block,dst_tiny_block,&! local block number for tiny blocks
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r8), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   real (r8) :: &
      xavg               ! scalar for enforcing symmetry at U pts
   integer :: omp_get_thread_num
   !logical (log_kind), save :: first_call = .true.
   !integer (int_kind), save :: bndy_2d_local, bndy_2d_recv, &
   !                            bndy_2d_send, bndy_2d_wait, bndy_2d_final

!-----------------------------------------------------------------------
!
!  update tiny-blocks within the same macro-block
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   rewrite to expose more parallelism to OpenMP threading
!-----------------------------------------------------------------------
!   ! east-west
!   !DIR$ IVDEP
!   !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(k,jblock,iblock,n,j,i)
!   do k=1,nblocks_tot_tiny
!      do jblock=1,ny_bkfactor
!      do iblock=1,nx_bkfactor-1
!         n = iblock + (jblock-1)*nx_bkfactor + k-1
!         do j=nghost+1,ny_block-nghost
!         !DIR$ UNROLL
!         do i=1,nghost
!         ARRAY(i,j,n+1) = ARRAY(i+block_size_x,j,n) ! west=>east
!         ARRAY(nghost+i+block_size_x,j,n) = ARRAY(nghost+i,j,n+1) ! east=>west
!         end do ! i
!         end do ! j
!      end do ! iblock
!      end do ! jblock
!   end do !k
!   !$OMP END PARALLEL DO
!   ! north-south
!   !$OMP PARALLEL DO PRIVATE(k,n,j,i)
!   do k=1,nblocks_tot_tiny
!      do n=k,k-1+(ny_bkfactor-1)*nx_bkfactor
!         !n = iblock + (jblock-1)*nx_bkfactor + k-1
!         do j=1,nghost
!         do i=1,nx_block
!         ARRAY(i,j,n+nx_bkfactor) = ARRAY(i,j+block_size_y,n) ! north=>south
!         ARRAY(i,nghost+j+block_size_y,n) = ARRAY(i,nghost+j,n+nx_bkfactor) ! south=>north
!         end do ! i
!         end do ! j
!      end do ! n
!      !$OMP END PARALLEL DO
!   end do !k
!   !$OMP END PARALLEL DO

!   !DIR$ IVDEP
!   do k=1,nblocks_tot_tiny,nx_bkfactor*ny_bkfactor
!       ! north-south
!      !$OMP PARALLEL DO PRIVATE(n,j,i)
!      do n=k,k-1+(ny_bkfactor-1)*nx_bkfactor
!         !n = iblock + (jblock-1)*nx_bkfactor + k-1
!         do j=1,nghost
!         do i=1,nx_block
!         ARRAY(i,j,n+nx_bkfactor) = ARRAY(i,j+block_size_y,n) ! north=>south
!         ARRAY(i,nghost+j+block_size_y,n) = ARRAY(i,nghost+j,n+nx_bkfactor) ! south=>north
!         end do ! i
!         end do ! j
!         write(6,"(A,7I5)")'orignaltest,task,ns_src_block,add1,2,dst_block,add1,2',my_task,n,1,1+block_size_y,n+nx_bkfactor,1,1
!         write(6,"(A,7I5)")'orignaltest,task,ns_src_block,add1,2,dst_block,add1,2',my_task,n+nx_bkfactor,1,nghost+1,n,1,nghost+1+block_size_y
!      end do ! n
!      !$OMP END PARALLEL DO
!  end do !k

!   !ljm tuning
!   if (lmic_trace) then
!   do j=1,loc_dist%local_block_num
!      i = loc_dist%local_block_ids(j)
!      if (i>nblocks_tot) then ! tiny-block
!        this_block = get_block(i,j) 
!        write(6,*) 'pcg-zero:RHS:',this_block%block_id,i-nblocks_tot,sum(ARRAY(:,:,j))
!      else
!        write(6,*) 'pcg-zero:RHS:',i,0,sum(ARRAY(:,:,j))
!      endif
!   enddo
!   call MPI_BARRIER(MPI_COMM_OCN, ierr)
!   call exit_POP(sigAbort,'debug pcg-zero-R...')
!   endif

!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   !if (first_call) then
   !  first_call = .false.
   !  call get_timer(bndy_2d_local, 'BNDY_2D_LOCAL')
   !  call get_timer(bndy_2d_recv,  'BNDY_2D_RECV')
   !  call get_timer(bndy_2d_send,  'BNDY_2D_SEND')
   !  call get_timer(bndy_2d_wait,  'BNDY_2D_WAIT')
   !  call get_timer(bndy_2d_final, 'BNDY_2D_FINAL')
   !endif

   allocate(buf_ew_snd(nghost, ny_block_mac, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block_mac, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_dbuf)) nx_global = size(tripole_dbuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block_mac*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ew_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !write(6,*) 'bound_dbl:IRECV:rcv',my_task,in_bndy%nmsg_ew_rcv
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block_mac*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         !ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n) ! macro-block
         !if (my_task == nproc_cpu_pn*nnode) &
         !   write(6,*) 'boundary_2d_dbl_tiny:snd_info:',my_task,src_block,ib_src
         !buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
         if (ib_src <= 2*nghost) then
            iblock = 1
         else
            iblock = nx_bkfactor
         endif
         !iblock=(ib_src-nghost-2)/block_size_x+1
         ib_src_tiny = ib_src-(iblock-1)*block_size_x
         ie_src_tiny = ib_src_tiny + nghost - 1
         !if (my_task == nproc_cpu_pn*nnode) &
         !   write(6,*) 'boundary_2d_dbl_tiny:snd_tiny:',my_task,iblock,ib_src_tiny
         src_tiny_block = iblock + (1-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ew_snd(:,1:nghost,i,n) = &
            ARRAY(ib_src_tiny:ie_src_tiny,1:nghost,src_tiny_block)
         !$OMP PARALLEL DO PRIVATE(jblock,src_tiny_block)
         do jblock=1,ny_bkfactor
            src_tiny_block = iblock + (jblock-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
            buf_ew_snd(:,nghost+(jblock-1)*block_size_y+1:nghost+jblock*block_size_y,i,n) = &
               ARRAY(ib_src_tiny:ie_src_tiny,nghost+1:nghost+block_size_y,src_tiny_block)
         end do ! jblock
         !$OMP END PARALLEL DO
         src_tiny_block = iblock + (ny_bkfactor-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ew_snd(:,nghost+(ny_bkfactor)*block_size_y+1:ny_block_mac,i,n) = &
            ARRAY(ib_src_tiny:ie_src_tiny,nghost+block_size_y+1:ny_block,src_tiny_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl_tiny:buff_ew_snd:',my_task,n,sum(buf_ew_snd(:,:,:,n))
      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !write(6,*) 'bound_dbl:ISEND:ew',my_task,in_bndy%nmsg_ew_snd
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_local)
!write(6,*)'sum local_ew_src_block_tiny',sum(in_bndy%local_ew_src_block_tiny(:))
!write(6,*)'sum local_ew_src_add_tiny1',sum(in_bndy%local_ew_src_add_tiny(1,:))
!write(6,*)'sum local_ew_dst_add_tiny1',sum(in_bndy%local_ew_dst_add_tiny(1,:))
   !$OMP PARALLEL DO PRIVATE(n,src_block,dst_block,ib_src,ib_dst,ie_src,ie_dst) schedule(static,2)
   do n=1,in_bndy%nlocal_ew_tiny
      src_block = in_bndy%local_ew_src_block_tiny(n)
      dst_block = in_bndy%local_ew_dst_block_tiny(n)
      !write(6,"(A,4I5)")'a my_task,thread,src,dst_block:',my_task,omp_get_thread_num(),src_block,dst_block
      ib_src = in_bndy%local_ew_src_add_tiny(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add_tiny(1,n)
      ie_dst = ib_dst + nghost - 1
      !write(6,*)'my_task+1',my_task+1,'src_block,ib_src,ie_src',src_block,ib_src,ie_src,'dst_block,ib_dst,ie_dst',dst_block,ib_dst,ie_dst

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = c0
      endif
   end do
   !$OMP END PARALLEL DO

!   if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl_tiny:update_local',my_task,sum(ARRAY(nghost+1:nx_block-nghost,nghost+1:ny_block-nghost,1:nx_bkfactor*ny_bkfactor))
   
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)
   !write(6,*) 'bound_dbl:WAITALL:ew_rcv',my_task,in_bndy%nmsg_ew_rcv
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      !ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
      iblock=(ib_dst-nghost-2)/block_size_x+1
      ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
      ie_dst_tiny = ib_dst_tiny + nghost - 1

      do jblock=1,ny_bkfactor
         dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
         ARRAY(ib_dst_tiny:ie_dst_tiny,1:ny_block,dst_tiny_block) = &
            buf_ew_rcv(1:nghost,(jblock-1)*block_size_y+1:2*nghost+jblock*block_size_y,k,n)
      end do

   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)
   !write(6,*) 'bound_dbl:WAITALL:ew_snd',my_task,in_bndy%nmsg_ew_snd
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block_mac, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block_mac, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_recv)
   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block_mac*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !write(6,*) 'bound_dbl:IRECV:ns_rcv',my_task,in_bndy%nmsg_ns_rcv
   !call timer_stop(bndy_2d_recv)

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_send)
   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block_mac*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         !buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
         if (jb_src <= 2*nghost) then
            jblock = 1
         else
            jblock = ny_bkfactor
         endif
         !jblock=(jb_src-nghost-2)/block_size_y+1
         jb_src_tiny = jb_src-(jblock-1)*block_size_y
         je_src_tiny = jb_src_tiny + nghost
         src_tiny_block = 1 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ns_snd(1:nghost,:,i,n) = &
            ARRAY(1:nghost,jb_src_tiny:je_src_tiny,src_tiny_block)
         do iblock=1,nx_bkfactor
         src_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ns_snd(nghost+(iblock-1)*block_size_x+1:nghost+iblock*block_size_x,:,i,n) = &
            ARRAY(nghost+1:nghost+block_size_x,jb_src_tiny:je_src_tiny,src_tiny_block)
         end do
         src_tiny_block = nx_bkfactor + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ns_snd(nghost+(nx_bkfactor)*block_size_x+1:nx_block_mac,:,i,n) = &
            ARRAY(nghost+block_size_x+1:nx_block,jb_src_tiny:je_src_tiny,src_tiny_block)
      end do

!      if (my_task == nproc_cpu_pn*nnode) &
!      write(6,*) 'bound_2d_dbl_tiny:buff_ns_snd:',my_task,n,sum(buf_ns_snd(:,:,:,n))
      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_dbl, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !write(6,*) 'bound_dbl:ISEND:ns_snd',my_task,in_bndy%nmsg_ns_snd
   !call timer_stop(bndy_2d_send)

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) tripole_dbuf = c0

   !$OMP PARALLEL DO PRIVATE(n,src_block,dst_block,ib_src,ib_dst,ie_src,ie_dst)
   do n=1,in_bndy%nlocal_ns_tiny
      src_block = in_bndy%local_ns_src_block_tiny(n)
      dst_block = in_bndy%local_ns_dst_block_tiny(n)
      !write(6,"(A,4I5)")'b my_task,thread,src,dst_block:',my_task,omp_get_thread_num(),src_block,dst_block
      if (dst_block > 0) then ! straight local copy
         jb_src = in_bndy%local_ns_src_add_tiny(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add_tiny(2,n)
         je_dst = jb_dst + nghost - 1
         !write(6,*)'my_task+1',my_task+1,'src_block,ib_src,ie_src',src_block,ib_src,ie_src,'dst_block,ib_dst,ie_dst',dst_block,ib_dst,ie_dst

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = c0
         endif
      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add_tiny(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add_tiny(1,n)
         !ie_dst = ib_dst + (nx_block-2*nghost) - 1
         ie_dst = ib_dst + block_size_x - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do
   !$OMP END PARALLEL DO
   
   !call timer_start(bndy_2d_local)
!   do n=1,in_bndy%nlocal_ns
!      src_block = in_bndy%local_ns_src_block(n)
!      dst_block = in_bndy%local_ns_dst_block(n)
!
!      if (dst_block > 0) then ! straight local copy
!
!         jb_src = in_bndy%local_ns_src_add(2,n)
!         je_src = jb_src + nghost - 1
!         jb_dst = in_bndy%local_ns_dst_add(2,n)
!         je_dst = jb_dst + nghost - 1
!write(6,"(A)")'xiaobintest0,task'
!         if (src_block /= 0) then
!write(6,"(A)")'xiaobintest0.1,task'
!            !ARRAY(:,jb_dst:je_dst,dst_block) = &
!            !ARRAY(:,jb_src:je_src,src_block)
!            jblock=(jb_src-nghost-2)/block_size_y+1
!            jb_src_tiny = jb_src-(jblock-1)*block_size_y
!            je_src_tiny = jb_src_tiny + nghost - 1
!            jblock2=(jb_dst-nghost-2)/block_size_y+1
!            jb_dst_tiny = jb_dst-(jblock2-1)*block_size_y
!            je_dst_tiny = jb_dst_tiny + nghost - 1
!            do iblock=1,nx_bkfactor
!               src_tiny_block = iblock + (jblock-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
!               dst_tiny_block = iblock + (jblock2-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
!               do j=1,nghost
!               do i=1,nx_block
!               ARRAY(i,jb_dst_tiny+j-1,dst_tiny_block) = &
!               ARRAY(i,jb_src_tiny+j-1,src_tiny_block)
!               end do ! i
!write(6,"(A,5I5)")'xiaobintest1,task,dst,src',my_task,dst_tiny_block,jb_dst_tiny+j-1,src_tiny_block,jb_src_tiny+j-1
!               end do ! j
!            end do
!         else
!            !ARRAY(:,jb_dst:je_dst,dst_block) = c0
!            jblock=(jb_dst-nghost-2)/block_size_y+1
!            jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
!            je_dst_tiny = jb_dst_tiny + nghost - 1
!            do iblock=1,nx_bkfactor
!               dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
!               ARRAY(:,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = c0
!               write(6,"(A,7I5)")'orignaltest,task,ns_src_block,add1,2,dst_block,add1,2',my_task,0,0,0,dst_tiny_block,jb_dst_tiny,1
!            end do
!         endif
!
!      else  !north boundary tripole grid - copy into global north buffer
!write(6,"(A,I5)")'xiaobintest1.2,task',my_task
!
!         jb_src = in_bndy%local_ns_src_add(2,n)
!         je_src = jb_src + nghost ! need nghost+1 rows for tripole
!
!         !*** determine start, end addresses of physical domain
!         !*** for both global buffer and local block
!
!         ib_dst = in_bndy%local_ns_src_add(1,n)
!         ie_dst = ib_dst + (nx_block_mac-2*nghost) - 1
!         if (ie_dst > nx_global) ie_dst = nx_global
!         ib_src = nghost + 1
!         ie_src = ib_src + ie_dst - ib_dst
!         if (src_block /= 0) then
!            !tripole_dbuf(ib_dst:ie_dst,:) = &
!            !      ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
!            jblock=(jb_src-nghost-2)/block_size_y+1
!            jb_src_tiny = jb_src-(jblock-1)*block_size_y
!            je_src_tiny = jb_src_tiny + nghost
!            iblock2=(ie_src-nghost-2)/block_size_x+1
!            ie_src_tiny = ie_src-(iblock2-1)*block_size_x
!            do iblock=1,iblock2-1
!            src_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
!            tripole_dbuf(ib_dst:ib_dst+block_size_x-1,:) = &
!               ARRAY(nghost+1:nghost+block_size_x,jb_src_tiny:je_src_tiny,src_tiny_block)
!            ib_dst = ib_dst + block_size_x
!            end do
!            src_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
!            tripole_dbuf(ib_dst:ie_dst,:) = &
!               ARRAY(nghost+1:ie_src_tiny,jb_src_tiny:je_src_tiny,src_tiny_block)
!         endif
!      endif
!   end do
   !call timer_stop(bndy_2d_local)

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)
   !write(6,*) 'bound_dbl:WAITALL:ns_rcv',my_task,in_bndy%nmsg_ns_rcv
   !call timer_stop(bndy_2d_wait)

   !call timer_start(bndy_2d_final)
   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      !??? src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         !ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
         jblock=(jb_dst-nghost-2)/block_size_y+1
         jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
         je_dst_tiny = jb_dst_tiny + nghost - 1
   
         do iblock=1,nx_bkfactor
            dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(1:nx_block,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               buf_ns_rcv((iblock-1)*block_size_x+1:2*nghost+iblock*block_size_x,1:nghost,k,n)
         end do
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_dbuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do
   !call timer_stop(bndy_2d_final)

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_dbuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_dghost(nx_global,1) = tripole_dbuf(nx_global,nghost+1)
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_dghost(:,1) = tripole_dbuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_dbuf(i     ,nghost+1)) + &
                       abs(tripole_dbuf(ib_dst,nghost+1)))
            tripole_dghost(i     ,1) = sign(xavg, &
                                            tripole_dbuf(i,nghost+1))
            tripole_dghost(ib_dst,1) = sign(xavg, &
                                            tripole_dbuf(ib_dst,nghost+1))
         end do
         tripole_dbuf(:,nghost+1) = tripole_dghost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_dghost(i,1+j) = isign* &
                                 tripole_dbuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary
               jblock=(jb_dst-nghost-2)/block_size_y+1
               jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
               je_dst_tiny = jb_dst_tiny + nghost

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               dst_tiny_block = 1 + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(1:nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_dghost(nx_global-nghost+1:nx_global,:)
               ie_src = ib_src + nx_block_mac - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block_mac - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               iblock=(ie_dst-nghost-2)/block_size_x+1
               ie_dst_tiny = ie_dst-(iblock-1)*block_size_x
               dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ie_dst_tiny+1:ie_dst_tiny+nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_dghost(1:nghost,:)
            endif

            !*** now copy the remaining ghost cell values

            !ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
            !   tripole_dghost(ib_src:ie_src,:)
            iblock=(ib_dst-nghost-2)/block_size_x+1
            ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
            iblock2=(ie_dst-nghost-2)/block_size_x+1
            ie_dst_tiny = ie_dst-(iblock2-1)*block_size_x
            do i=iblock,iblock2-1
               dst_tiny_block = i + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ib_dst_tiny:nx_block,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_dghost(ib_src:ib_src+(nx_block-ib_dst_tiny),:)
               ib_src = ib_src + (block_size_x-ib_dst_tiny+1)
               ib_dst_tiny = 1
            end do
            dst_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(ib_dst_tiny:ie_dst_tiny,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary
               jblock=(jb_dst-nghost-2)/block_size_y+1
               jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
               je_dst_tiny = jb_dst_tiny + nghost

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               dst_tiny_block = 1 + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(1:nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_dghost(nx_global-nghost+1:nx_global,:)
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               iblock=(ie_dst-nghost-2)/block_size_x+1
               ie_dst_tiny = ie_dst-(iblock-1)*block_size_x
               dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ie_dst_tiny+1:ie_dst_tiny+nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_dghost(1:nghost,:)
            endif

            !ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
            !   tripole_dghost(ib_src:ie_src,:)
            iblock=(ib_dst-nghost-2)/block_size_x+1
            ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
            iblock2=(ie_dst-nghost-2)/block_size_x+1
            ie_dst_tiny = ie_dst-(iblock2-1)*block_size_x
            do i=iblock,iblock2-1
               dst_tiny_block = i + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ib_dst_tiny:nx_block,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_dghost(ib_src:ib_src+(nx_block-ib_dst_tiny),:)
               ib_src = ib_src + (block_size_x-ib_dst_tiny+1)
               ib_dst_tiny = 1
            end do
            dst_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(ib_dst_tiny:ie_dst_tiny,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               tripole_dghost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)
   !write(6,*) 'bound_dbl:WAITALL:ns_snd',my_task,in_bndy%nmsg_ns_snd
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------

 end subroutine boundary_2d_dbl_tiny


!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

 subroutine boundary_2d_real(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   real (r4), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r4), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   real (r4) :: &
      xavg               ! scalar for enforcing symmetry at U pts

   if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   allocate(buf_ew_snd(nghost, ny_block, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_rbuf)) nx_global = size(tripole_rbuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_real,  &
                     in_bndy%ew_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n)
         buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
      end do

      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_real, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = c0
      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)

   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
   end do
   end do

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_real,  &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
      end do

      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_real, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_rbuf)) tripole_rbuf = c0

   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = c0
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_rbuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)

   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_rbuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_rbuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_rghost(:,1) = tripole_rbuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_rbuf(i     ,nghost+1)) + &
                       abs(tripole_rbuf(ib_dst,nghost+1)))
            tripole_rghost(i     ,1) = sign(xavg, &
                                            tripole_rbuf(i,nghost+1))
            tripole_rghost(ib_dst,1) = sign(xavg, &
                                            tripole_rbuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_rghost(nx_global,1) = tripole_rbuf(nx_global,nghost+1)
         tripole_rbuf(:,nghost+1) = tripole_rghost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_rghost(:,1) = tripole_rbuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_rbuf(i     ,nghost+1)) + &
                       abs(tripole_rbuf(ib_dst,nghost+1)))
            tripole_rghost(i     ,1) = sign(xavg, &
                                            tripole_rbuf(i,nghost+1))
            tripole_rghost(ib_dst,1) = sign(xavg, &
                                            tripole_rbuf(ib_dst,nghost+1))
         end do
         tripole_rbuf(:,nghost+1) = tripole_rghost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_rghost(i,1+j) = isign* &
                                 tripole_rbuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(i,:)
               end do
            endif

            !*** now copy the remaining ghost cell values

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_rghost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(i,:)
               end do
            endif

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_rghost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

   else   ! lmic_proc
      call boundary_2d_real_tiny(ARRAY, in_bndy, grid_loc, field_type)
   endif  ! lmic_proc
!-----------------------------------------------------------------------

end subroutine boundary_2d_real

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

 subroutine boundary_2d_real_tiny(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   real (r4), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) ::     &
      this_block        ! block information for current block

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      ib_src_tiny,ie_src_tiny,jb_src_tiny,je_src_tiny, &!
      ib_dst_tiny,ie_dst_tiny,jb_dst_tiny,je_dst_tiny, &!
      iblock,jblock,iblock2,jblock2,&!
      src_tiny_block,dst_tiny_block,&! local block number for tiny blocks
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   real (r4), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   real (r4) :: &
      xavg               ! scalar for enforcing symmetry at U pts

         !if (my_task == nproc_cpu_pn*nnode) &
         call exit_POP(sigAbort,'boundary_2d_real_tiny:NOT IMPL.')
!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   allocate(buf_ew_snd(nghost, ny_block, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_rbuf)) nx_global = size(tripole_rbuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_real,  &
                     in_bndy%ew_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n)
         buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
      end do

      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_real, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = c0
      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)

   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
   end do
   end do

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_real,  &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
      end do

      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_real, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_rbuf)) tripole_rbuf = c0

   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = c0
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_rbuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)

   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_rbuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_rbuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_rghost(:,1) = tripole_rbuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_rbuf(i     ,nghost+1)) + &
                       abs(tripole_rbuf(ib_dst,nghost+1)))
            tripole_rghost(i     ,1) = sign(xavg, &
                                            tripole_rbuf(i,nghost+1))
            tripole_rghost(ib_dst,1) = sign(xavg, &
                                            tripole_rbuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_rghost(nx_global,1) = tripole_rbuf(nx_global,nghost+1)
         tripole_rbuf(:,nghost+1) = tripole_rghost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_rghost(:,1) = tripole_rbuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_rbuf(i     ,nghost+1)) + &
                       abs(tripole_rbuf(ib_dst,nghost+1)))
            tripole_rghost(i     ,1) = sign(xavg, &
                                            tripole_rbuf(i,nghost+1))
            tripole_rghost(ib_dst,1) = sign(xavg, &
                                            tripole_rbuf(ib_dst,nghost+1))
         end do
         tripole_rbuf(:,nghost+1) = tripole_rghost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_rghost(i,1+j) = isign* &
                                 tripole_rbuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(i,:)
               end do
            endif

            !*** now copy the remaining ghost cell values

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_rghost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_rghost(i,:)
               end do
            endif

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_rghost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------

end subroutine boundary_2d_real_tiny

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

 subroutine boundary_2d_int(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   integer (int_kind), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   integer (int_kind) :: &
      xavg               ! scalar for enforcing symmetry at U pts

   if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   allocate(buf_ew_snd(nghost, ny_block, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_ibuf)) nx_global = size(tripole_ibuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_integer, &
                     in_bndy%ew_rcv_proc(n)-1,                  &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n),   &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n)
         buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
      end do

      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_integer, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         ARRAY(ib_dst:ie_dst,:,dst_block) = &
         ARRAY(ib_src:ie_src,:,src_block)
      else
         ARRAY(ib_dst:ie_dst,:,dst_block) = 0
      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)

   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
   end do
   end do

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_integer,   &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
      end do

      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_integer, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_ibuf)) tripole_ibuf = c0

   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            ARRAY(:,jb_dst:je_dst,dst_block) = &
            ARRAY(:,jb_src:je_src,src_block)
         else
            ARRAY(:,jb_dst:je_dst,dst_block) = 0
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_ibuf(ib_dst:ie_dst,:) = &
                  ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
         endif

      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)

   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_ibuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_ibuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_ighost(:,1) = tripole_ibuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_ibuf(i     ,nghost+1)) + &
                       abs(tripole_ibuf(ib_dst,nghost+1)))
            tripole_ighost(i     ,1) = sign(xavg, &
                                          tripole_ibuf(i,nghost+1))
            tripole_ighost(ib_dst,1) = sign(xavg, &
                                          tripole_ibuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_ighost(nx_global,1) = tripole_ibuf(nx_global,nghost+1)
         tripole_ibuf(:,nghost+1) = tripole_ighost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_ighost(:,1) = tripole_ibuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_ibuf(i     ,nghost+1)) + &
                       abs(tripole_ibuf(ib_dst,nghost+1)))
            tripole_ighost(i     ,1) = sign(xavg, &
                                          tripole_ibuf(i,nghost+1))
            tripole_ighost(ib_dst,1) = sign(xavg, &
                                          tripole_ibuf(ib_dst,nghost+1))
         end do
         tripole_ibuf(:,nghost+1) = tripole_ighost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_ighost(i,1+j) = isign* &
                                 tripole_ibuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_ighost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_ighost(i,:)
               end do
            endif

            !*** now copy the remaining ghost cell values

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_ighost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               do i=1,nghost
                  ARRAY(i,jb_dst:je_dst,dst_block) = &
                     tripole_ighost(nx_global-nghost+i,:)
               end do
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               do i=1,nghost
                  ARRAY(ie_dst+i,jb_dst:je_dst,dst_block) = &
                     tripole_ighost(i,:)
               end do
            endif

            ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
               tripole_ighost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

   else   ! lmic_proc
      call boundary_2d_int_tiny(ARRAY, in_bndy, grid_loc, field_type)
   endif  ! lmic_proc
!-----------------------------------------------------------------------

end subroutine boundary_2d_int

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

 subroutine boundary_2d_int_tiny(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 2d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

   include 'mpif.h'   ! MPI Fortran include file

! !INPUT PARAMETERS:

   type (bndy), intent(in) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)
! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      i,j,k,m,n,                   &! dummy loop indices
      ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
      ib_dst,ie_dst,jb_dst,je_dst, &!
      ib_src_tiny,ie_src_tiny,jb_src_tiny,je_src_tiny, &!
      ib_dst_tiny,ie_dst_tiny,jb_dst_tiny,je_dst_tiny, &!
      iblock,jblock,iblock2,jblock2,&!
      src_tiny_block,dst_tiny_block,&! local block number for tiny blocks
      nx_global,                   &! global domain size in x
      src_block,                   &! local block number for source
      dst_block,                   &! local block number for destination
      bufsize,                     &! buffer size for send/recv buffers
      xoffset, yoffset,            &! address shifts for tripole
      isign,                       &! sign factor for tripole grids
      ierr                          ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      snd_request,              &! MPI request ids
      rcv_request                ! MPI request ids

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status,               &! MPI status flags
      rcv_status                 ! MPI status flags

   integer (int_kind), dimension(:,:,:,:), allocatable :: &
      buf_ew_snd,       &! message buffer for east-west sends
      buf_ew_rcv,       &! message buffer for east-west recvs
      buf_ns_snd,       &! message buffer for north-south sends
      buf_ns_rcv         ! message buffer for north-south recvs

   integer (int_kind) :: &
      xavg               ! scalar for enforcing symmetry at U pts

!-----------------------------------------------------------------------
!
!  update tiny-blocks within the same macro-block
!
!-----------------------------------------------------------------------
   do k=1,nblocks_tot_tiny,nx_bkfactor*ny_bkfactor
      ! east-west
      do jblock=1,ny_bkfactor
      do iblock=1,nx_bkfactor-1
         m = iblock + (jblock-1)*nx_bkfactor + k-1
         n = m + 1
         ARRAY(1:nghost,:,n) = ARRAY(nx_block-2*nghost+1:nx_block-nghost,:,m) ! west=>east
         ARRAY(nx_block-nghost+1:nx_block,:,m) = ARRAY(nghost+1:2*nghost,:,n) ! east=>west
      end do ! iblock
      end do ! jblock
       ! north-south
      do jblock=1,ny_bkfactor-1
      do iblock=1,nx_bkfactor
         m = iblock + (jblock-1)*nx_bkfactor + k-1
         n = m + nx_bkfactor
         ARRAY(:,1:nghost,n) = ARRAY(:,ny_block-2*nghost+1:ny_block-nghost,m) ! north=>south
         ARRAY(:,ny_block-nghost+1:ny_block,m) = ARRAY(:,nghost+1:2*nghost,n) ! south=>north
      end do ! iblock
      end do ! jblock
  end do !k

!-----------------------------------------------------------------------
!
!  allocate buffers for east-west sends and receives
!
!-----------------------------------------------------------------------

   allocate(buf_ew_snd(nghost, ny_block_mac, &
                       in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block_mac, &
                       in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

   allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

   if (allocated(tripole_ibuf)) nx_global = size(tripole_ibuf,dim=1)

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_rcv

      bufsize = ny_block_mac*nghost*in_bndy%nblocks_ew_rcv(n)

      call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_integer, &
                     in_bndy%ew_rcv_proc(n)-1,                  &
                     mpitag_bndy_2d + in_bndy%ew_rcv_proc(n),   &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !write(6,*) 'bound_int:IRECV:rcv',my_task,in_bndy%nmsg_ew_rcv

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ew_snd

      bufsize = ny_block_mac*nghost*in_bndy%nblocks_ew_snd(n)

      do i=1,in_bndy%nblocks_ew_snd(n)
         ib_src    = in_bndy%ew_src_add(1,i,n)
         ie_src    = ib_src + nghost - 1
         src_block = in_bndy%ew_src_block(i,n) ! macro-block
         !if (my_task == nproc_cpu_pn*nnode) &
         !   write(6,*) 'boundary_2d_dbl_tiny:snd_info:',my_task,src_block,ib_src
         !buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
         iblock=(ib_src-nghost-2)/block_size_x+1
         ib_src_tiny = ib_src-(iblock-1)*block_size_x
         ie_src_tiny = ib_src_tiny + nghost - 1
         !if (my_task == nproc_cpu_pn*nnode) &
         !   write(6,*) 'boundary_2d_dbl_tiny:snd_tiny:',my_task,iblock,ib_src_tiny
         src_tiny_block = iblock + (1-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ew_snd(:,1:nghost+block_size_y,i,n) = &
            ARRAY(ib_src_tiny:ie_src_tiny,1:nghost+block_size_y,src_tiny_block)
         do jblock=2,ny_bkfactor-1
            src_tiny_block = iblock + (jblock-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
            buf_ew_snd(:,nghost+(jblock-1)*block_size_y+1:nghost+jblock*block_size_y,i,n) = &
               ARRAY(ib_src_tiny:ie_src_tiny,nghost+1:nghost+block_size_y,src_tiny_block)
         end do
         src_tiny_block = iblock + (ny_bkfactor-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ew_snd(:,nghost+(ny_bkfactor-1)*block_size_y+1:ny_block_mac,i,n) = &
            ARRAY(ib_src_tiny:ie_src_tiny,nghost+1:ny_block,src_tiny_block)
      end do

      call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_integer, &
                     in_bndy%ew_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !write(6,*) 'bound_int:ISEND:ew',my_task,in_bndy%nmsg_ew_snd

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!  also initialize ghost cells to zero
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nlocal_ew
      src_block = in_bndy%local_ew_src_block(n)
      dst_block = in_bndy%local_ew_dst_block(n)

      ib_src = in_bndy%local_ew_src_add(1,n)
      ie_src = ib_src + nghost - 1
      ib_dst = in_bndy%local_ew_dst_add(1,n)
      ie_dst = ib_dst + nghost - 1

      if (src_block /= 0) then
         !ARRAY(ib_dst:ie_dst,:,dst_block) = &
         !ARRAY(ib_src:ie_src,:,src_block)
         iblock=(ib_src-nghost-2)/block_size_x+1
         ib_src_tiny = ib_src-(iblock-1)*block_size_x
         ie_src_tiny = ib_src_tiny + nghost - 1
         iblock2=(ib_dst-nghost-2)/block_size_x+1
         ib_dst_tiny = ib_dst-(iblock2-1)*block_size_x
         ie_dst_tiny = ib_dst_tiny + nghost - 1
         do jblock=1,ny_bkfactor
            src_tiny_block = iblock + (jblock-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
            dst_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(ib_dst_tiny:ie_dst_tiny,:,dst_tiny_block) = &
            ARRAY(ib_src_tiny:ie_src_tiny,:,src_tiny_block)
         end do
      else
         !ARRAY(ib_dst:ie_dst,:,dst_block) = 0
         iblock=(ib_dst-nghost-2)/block_size_x+1
         ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
         ie_dst_tiny = ib_dst_tiny + nghost - 1
         do jblock=1,ny_bkfactor
            dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(ib_dst_tiny:ie_dst_tiny,:,dst_tiny_block) = 0
         end do
      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)
   !write(6,*) 'bound_int:WAITALL:ew_rcv',my_task,in_bndy%nmsg_ew_rcv

   do n=1,in_bndy%nmsg_ew_rcv
   do k=1,in_bndy%nblocks_ew_rcv(n)
      dst_block = in_bndy%ew_dst_block(k,n)

      ib_dst = in_bndy%ew_dst_add(1,k,n)
      ie_dst = ib_dst + nghost - 1

      !ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
      iblock=(ib_dst-nghost-2)/block_size_x+1 ! Fix bug:further minus 1 or 2
      ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
      ie_dst_tiny = ib_dst_tiny + nghost - 1

      do jblock=1,ny_bkfactor
         dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
         ARRAY(ib_dst_tiny:ie_dst_tiny,1:ny_block,dst_tiny_block) = &
            buf_ew_rcv(1:nghost,(jblock-1)*block_size_y+1:2*nghost+jblock*block_size_y,k,n)
      end do

   end do
   end do

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)
   !write(6,*) 'bound_int:WAITALL:ew_snd',my_task,in_bndy%nmsg_ew_snd

   deallocate(buf_ew_snd, buf_ew_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------
!
!  now exchange north-south boundary info
!
!-----------------------------------------------------------------------

   allocate(buf_ns_snd(nx_block_mac, nghost+1, &
                       in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block_mac, nghost+1, &
                       in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

   allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

!-----------------------------------------------------------------------
!
!  post receives
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_rcv

      bufsize = nx_block_mac*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

      call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_integer,   &
                     in_bndy%ns_rcv_proc(n)-1,                &
                     mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                     in_bndy%communicator, rcv_request(n), ierr)
   end do
   !write(6,*) 'bound_int:IRECV:ns_rcv',my_task,in_bndy%nmsg_ns_rcv

!-----------------------------------------------------------------------
!
!  fill send buffer and post sends
!
!-----------------------------------------------------------------------

   do n=1,in_bndy%nmsg_ns_snd

      bufsize = nx_block_mac*(nghost+1)*in_bndy%nblocks_ns_snd(n)

      do i=1,in_bndy%nblocks_ns_snd(n)
         jb_src    = in_bndy%ns_src_add(2,i,n)
         je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
         src_block = in_bndy%ns_src_block(i,n)
         !buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
         jblock=(jb_src-nghost-2)/block_size_y+1
         jb_src_tiny = jb_src-(jblock-1)*block_size_y
         je_src_tiny = jb_src_tiny + nghost - 1
         src_tiny_block = 1 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ns_snd(1:nghost+block_size_x,:,i,n) = &
            ARRAY(1:nghost+block_size_x,jb_src_tiny:je_src_tiny,src_tiny_block)
         do iblock=2,nx_bkfactor-1
         src_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ns_snd(nghost+(iblock-1)*block_size_x+1:nghost+iblock*block_size_x,:,i,n) = &
            ARRAY(nghost+1:nghost+block_size_x,jb_src_tiny:je_src_tiny,src_tiny_block)
         end do
         src_tiny_block = nx_bkfactor + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
         buf_ns_snd(nghost+(nx_bkfactor-1)*block_size_x+1:nx_block_mac,:,i,n) = &
            ARRAY(nghost+1:nx_block,jb_src_tiny:je_src_tiny,src_tiny_block)
      end do

      call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_integer, &
                     in_bndy%ns_snd_proc(n)-1, &
                     mpitag_bndy_2d + my_task + 1, &
                     in_bndy%communicator, snd_request(n), ierr)
   end do
   !write(6,*) 'bound_int:ISEND:ns_snd',my_task,in_bndy%nmsg_ns_snd

!-----------------------------------------------------------------------
!
!  do local copies while waiting for messages to complete
!
!-----------------------------------------------------------------------

   if (allocated(tripole_ibuf)) tripole_ibuf = 0

   do n=1,in_bndy%nlocal_ns
      src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%local_ns_dst_block(n)

      if (dst_block > 0) then ! straight local copy

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost - 1
         jb_dst = in_bndy%local_ns_dst_add(2,n)
         je_dst = jb_dst + nghost - 1

         if (src_block /= 0) then
            !ARRAY(:,jb_dst:je_dst,dst_block) = &
            !ARRAY(:,jb_src:je_src,src_block)
            jblock=(jb_src-nghost-2)/block_size_y+1
            jb_src_tiny = jb_src-(jblock-1)*block_size_y
            je_src_tiny = jb_src_tiny + nghost - 1
            jblock2=(jb_dst-nghost-2)/block_size_y+1
            jb_dst_tiny = jb_dst-(jblock2-1)*block_size_y
            je_dst_tiny = jb_dst_tiny + nghost - 1
            do iblock=1,nx_bkfactor
               src_tiny_block = iblock + (jblock-1)*nx_bkfactor+(src_block-1)*nx_bkfactor*ny_bkfactor
               dst_tiny_block = iblock + (jblock2-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(:,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               ARRAY(:,jb_src_tiny:je_src_tiny,src_tiny_block)
            end do
         else
            !ARRAY(:,jb_dst:je_dst,dst_block) = 0
            jblock=(jb_dst-nghost-2)/block_size_y+1
            jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
            je_dst_tiny = jb_dst_tiny + nghost - 1
            do iblock=1,nx_bkfactor
               dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(:,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = 0
            end do
         endif

      else  !north boundary tripole grid - copy into global north buffer

         jb_src = in_bndy%local_ns_src_add(2,n)
         je_src = jb_src + nghost ! need nghost+1 rows for tripole

         !*** determine start, end addresses of physical domain
         !*** for both global buffer and local block

         ib_dst = in_bndy%local_ns_src_add(1,n)
         ie_dst = ib_dst + (nx_block_mac-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            !tripole_ibuf(ib_dst:ie_dst,:) = &
            !      ARRAY(ib_src:ie_src,jb_src:je_src,src_block)
            jblock=(jb_src-nghost-2)/block_size_y+1
            jb_src_tiny = jb_src-(jblock-1)*block_size_y
            je_src_tiny = jb_src_tiny + nghost
            iblock2=(ie_src-nghost-2)/block_size_x+1
            ie_src_tiny = ie_src-(iblock2-1)*block_size_x
            do iblock=1,iblock2-1
            src_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
            tripole_ibuf(ib_dst:ib_dst+block_size_x-1,:) = &
               ARRAY(nghost+1:nghost+block_size_x,jb_src_tiny:je_src_tiny,src_tiny_block)
            ib_dst = ib_dst + block_size_x
            end do
            src_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
            tripole_ibuf(ib_dst:ie_dst,:) = &
               ARRAY(nghost+1:ie_src_tiny,jb_src_tiny:je_src_tiny,src_tiny_block)
         endif

      endif
   end do

!-----------------------------------------------------------------------
!
!  wait for receives to finish and then unpack the recv buffer into
!  ghost cells
!
!-----------------------------------------------------------------------

   call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)
   !write(6,*) 'bound_int:WAITALL:ns_rcv',my_task,in_bndy%nmsg_ns_rcv

   do n=1,in_bndy%nmsg_ns_rcv
   do k=1,in_bndy%nblocks_ns_rcv(n)
      !??? src_block = in_bndy%local_ns_src_block(n)
      dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

      if (dst_block > 0) then  ! normal receive
         jb_dst = in_bndy%ns_dst_add(2,k,n)
         je_dst = jb_dst + nghost - 1

         !ARRAY(:,jb_dst:je_dst,dst_block) = buf_ns_rcv(:,1:nghost,k,n)
         jblock=(jb_dst-nghost-2)/block_size_y+1
         jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
         je_dst_tiny = jb_dst_tiny + nghost - 1
   
         do iblock=1,nx_bkfactor
            dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+(dst_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(1:nx_block,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               buf_ns_rcv((iblock-1)*block_size_x+1:2*nghost+iblock*block_size_x,1:nghost,k,n)
         end do
      else ! northern tripole bndy: copy into global tripole buffer

         !*** determine start,end of physical domain for both
         !*** global buffer and local buffer

         ib_dst = in_bndy%ns_dst_add(1,k,n)
         ie_dst = ib_dst + (nx_block-2*nghost) - 1
         if (ie_dst > nx_global) ie_dst = nx_global
         ib_src = nghost + 1
         ie_src = ib_src + ie_dst - ib_dst
         if (src_block /= 0) then
            tripole_ibuf(ib_dst:ie_dst,:) = &
            buf_ns_rcv(ib_src:ie_src,:,k,n)
         endif
      endif
   end do
   end do

!-----------------------------------------------------------------------
!
!  take care of northern boundary in tripole case
!
!-----------------------------------------------------------------------

   if (allocated(tripole_ibuf)) then

      select case (grid_loc)
      case (field_loc_center)   ! cell center location
         xoffset = 1
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_ighost(:,1) = tripole_ibuf(:,nghost+1)
      case (field_loc_NEcorner)   ! cell corner (velocity) location
         xoffset = 0
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global - i
            if (ib_dst == 0) ib_dst = nx_global
            xavg = p5*(abs(tripole_ibuf(i     ,nghost+1)) + &
                       abs(tripole_ibuf(ib_dst,nghost+1)))
            tripole_ighost(i     ,1) = sign(xavg, &
                                          tripole_ibuf(i,nghost+1))
            tripole_ighost(ib_dst,1) = sign(xavg, &
                                          tripole_ibuf(ib_dst,nghost+1))
         end do
         !*** catch nx_global point
         tripole_ighost(nx_global,1) = tripole_ibuf(nx_global,nghost+1)
         tripole_ibuf(:,nghost+1) = tripole_ighost(:,1)
      case (field_loc_Eface)   ! cell center location
         xoffset = 0
         yoffset = 1
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain (mostly for symmetry enforcement)
         tripole_ighost(:,1) = tripole_ibuf(:,nghost+1)
      case (field_loc_Nface)   ! cell corner (velocity) location
         xoffset = 1
         yoffset = 0
         !*** enforce symmetry
         !*** first row of ghost cell buffer is actually the last
         !*** row of physical domain
         do i = 1, nx_global/2
            ib_dst = nx_global + 1 - i
            xavg = p5*(abs(tripole_ibuf(i     ,nghost+1)) + &
                       abs(tripole_ibuf(ib_dst,nghost+1)))
            tripole_ighost(i     ,1) = sign(xavg, &
                                          tripole_ibuf(i,nghost+1))
            tripole_ighost(ib_dst,1) = sign(xavg, &
                                          tripole_ibuf(ib_dst,nghost+1))
         end do
         tripole_ibuf(:,nghost+1) = tripole_ighost(:,1)
      case default
         call exit_POP(sigAbort, 'Unknown location in boundary_2d')
      end select

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case default
         call exit_POP(sigAbort, 'Unknown field type in boundary')
      end select

      !*** copy source (physical) cells into ghost cells
      !*** global source addresses are:
      !*** nx_global + xoffset - i
      !*** ny_global + yoffset - j
      !*** in the actual tripole buffer, the indices are:
      !*** nx_global + xoffset - i = ib_src - i
      !*** ny_global + yoffset - j - (ny_global - nghost) + 1 =
      !***    nghost + yoffset +1 - j = jb_src - j

      ib_src = nx_global + xoffset
      jb_src = nghost + yoffset + 1

      do j=1,nghost
      do i=1,nx_global
         tripole_ighost(i,1+j) = isign* &
                                 tripole_ibuf(ib_src-i, jb_src-j)
      end do
      end do

      !*** copy out of global ghost cell buffer into local
      !*** ghost cells

      do n=1,in_bndy%nlocal_ns
         dst_block = in_bndy%local_ns_dst_block(n)

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%local_ns_dst_add(2,n)
            je_dst = jb_dst + nghost
            ib_src = in_bndy%local_ns_dst_add(1,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary
               jblock=(jb_dst-nghost-2)/block_size_y+1
               jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
               je_dst_tiny = jb_dst_tiny + nghost

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               dst_tiny_block = 1 + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(1:nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_ighost(nx_global-nghost+1:nx_global,:)
               ie_src = ib_src + nx_block_mac - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block_mac - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               iblock=(ie_dst-nghost-2)/block_size_x+1
               ie_dst_tiny = ie_dst-(iblock-1)*block_size_x
               dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ie_dst_tiny+1:ie_dst_tiny+nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_ighost(1:nghost,:)
            endif

            !*** now copy the remaining ghost cell values

            !ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
            !   tripole_ighost(ib_src:ie_src,:)
            iblock=(ib_dst-nghost-2)/block_size_x+1
            ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
            iblock2=(ie_dst-nghost-2)/block_size_x+1
            ie_dst_tiny = ie_dst-(iblock2-1)*block_size_x
            do i=iblock,iblock2-1
               dst_tiny_block = i + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ib_dst_tiny:nx_block,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_ighost(ib_src:ib_src+(nx_block-ib_dst_tiny),:)
               ib_src = ib_src + (block_size_x-ib_dst_tiny+1)
               ib_dst_tiny = 1
            end do
            dst_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(ib_dst_tiny:ie_dst_tiny,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               tripole_ighost(ib_src:ie_src,:)
         endif

      end do

      do n=1,in_bndy%nmsg_ns_rcv
      do k=1,in_bndy%nblocks_ns_rcv(n)
         dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

         if (dst_block < 0) then
            dst_block = -dst_block

            jb_dst = in_bndy%ns_dst_add(2,k,n)
            je_dst = jb_dst + nghost ! last phys row incl for symmetry
            ib_src = in_bndy%ns_dst_add(3,k,n)
            !*** ib_src is glob address of 1st point in physical
            !*** domain.  must now adjust to properly copy
            !*** east-west ghost cell info in the north boundary
               jblock=(jb_dst-nghost-2)/block_size_y+1
               jb_dst_tiny = jb_dst-(jblock-1)*block_size_y
               je_dst_tiny = jb_dst_tiny + nghost

            if (ib_src == 1) then  ! western boundary
               !*** impose cyclic conditions at western boundary
               !*** then set up remaining indices to copy rest
               !*** of domain from tripole ghost cell buffer
               dst_tiny_block = 1 + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(1:nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_ighost(nx_global-nghost+1:nx_global,:)
               ie_src = ib_src + nx_block - nghost - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = nghost + 1
               ie_dst = ib_dst + (ie_src - ib_src)
            else
               ib_src = ib_src - nghost
               ie_src = ib_src + nx_block - 1
               if (ie_src > nx_global) ie_src = nx_global
               ib_dst = 1
               ie_dst = ib_dst + (ie_src - ib_src)
            endif
            if (ie_src == nx_global) then ! eastern boundary
               !*** impose cyclic conditions in ghost cells
               iblock=(ie_dst-nghost-2)/block_size_x+1
               ie_dst_tiny = ie_dst-(iblock-1)*block_size_x
               dst_tiny_block = iblock + (jblock-1)*nx_bkfactor+ (dst_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ie_dst_tiny+1:ie_dst_tiny+nghost,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_ighost(1:nghost,:)
            endif

            !ARRAY(ib_dst:ie_dst,jb_dst:je_dst,dst_block) = &
            !   tripole_ighost(ib_src:ie_src,:)
            iblock=(ib_dst-nghost-2)/block_size_x+1
            ib_dst_tiny = ib_dst-(iblock-1)*block_size_x
            iblock2=(ie_dst-nghost-2)/block_size_x+1
            ie_dst_tiny = ie_dst-(iblock2-1)*block_size_x
            do i=iblock,iblock2-1
               dst_tiny_block = i + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
               ARRAY(ib_dst_tiny:nx_block,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
                  tripole_ighost(ib_src:ib_src+(nx_block-ib_dst_tiny),:)
               ib_src = ib_src + (block_size_x-ib_dst_tiny+1)
               ib_dst_tiny = 1
            end do
            dst_tiny_block = iblock2 + (jblock-1)*nx_bkfactor+ (src_block-1)*nx_bkfactor*ny_bkfactor
            ARRAY(ib_dst_tiny:ie_dst_tiny,jb_dst_tiny:je_dst_tiny,dst_tiny_block) = &
               tripole_ighost(ib_src:ie_src,:)
         endif


      end do
      end do

   endif

!-----------------------------------------------------------------------
!
!  wait for sends to complete and deallocate arrays
!
!-----------------------------------------------------------------------

   !call timer_start(bndy_2d_wait)
   call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)
   !write(6,*) 'bound_int:WAITALL:ns_snd',my_task,in_bndy%nmsg_ns_snd
   !call timer_stop(bndy_2d_wait)

   deallocate(buf_ns_snd, buf_ns_rcv)
   deallocate(snd_request, rcv_request, snd_status, rcv_status)

!-----------------------------------------------------------------------

end subroutine boundary_2d_int_tiny

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

subroutine boundary_3d_dbl(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 3d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   real (r8), dimension(:,:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      k,m                           ! dummy loop indices

!-----------------------------------------------------------------------

   m = size(ARRAY,3)
   do k = 1, m
      call boundary_2d_dbl(ARRAY(:,:,k,:),in_bndy,grid_loc,field_type)
   end do

!-----------------------------------------------------------------------

end subroutine boundary_3d_dbl

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

subroutine boundary_3d_real(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 3d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   real (r4), dimension(:,:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      k,m                           ! dummy loop indices

!-----------------------------------------------------------------------

   m = size(ARRAY,3)
   do k = 1, m
      call boundary_2d_real(ARRAY(:,:,k,:),in_bndy,grid_loc,field_type)
   end do

!-----------------------------------------------------------------------

end subroutine boundary_3d_real

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

subroutine boundary_3d_int(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 3d horizontal arrays of integer.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   integer (int_kind), dimension(:,:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      k,m                           ! dummy loop indices

!-----------------------------------------------------------------------

   m = size(ARRAY,3)
   do k = 1, m
      call boundary_2d_int(ARRAY(:,:,k,:),in_bndy,grid_loc,field_type)
   end do

!-----------------------------------------------------------------------

end subroutine boundary_3d_int

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

subroutine boundary_4d_dbl(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 3d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   real (r8), dimension(:,:,:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      k,l,m,n                      ! dummy loop indices

!-----------------------------------------------------------------------

   l = size(ARRAY,dim=3)
   n = size(ARRAY,dim=4)
   do k=1,l
   do m=1,n
      call boundary_2d_dbl(ARRAY(:,:,k,m,:),in_bndy,grid_loc,field_type)
   end do
   end do

!-----------------------------------------------------------------------

end subroutine boundary_4d_dbl

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

subroutine boundary_4d_real(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 3d horizontal arrays of single precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   real (r4), dimension(:,:,:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      k,l,m,n                      ! dummy loop indices

!-----------------------------------------------------------------------

   l = size(ARRAY,dim=3)
   n = size(ARRAY,dim=4)
   do k=1,l
   do m=1,n
     call boundary_2d_real(ARRAY(:,:,k,m,:),in_bndy,grid_loc,field_type)
   end do
   end do

!-----------------------------------------------------------------------

end subroutine boundary_4d_real

!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells
! !INTERFACE:

subroutine boundary_4d_int(ARRAY, in_bndy, grid_loc, field_type)

! !DESCRIPTION:
!  This routine updates ghost cells for an input array and is a
!  member of a group of routines under the generic interface
!  update\_ghost\_cells.  This routine is the specific interface
!  for 3d horizontal arrays of double precision.
!
! !REVISION HISTORY:
!  same as module

! !USER:

! !INPUT/OUTPUT PARAMETERS:

   type (bndy), intent(inout) :: &
      in_bndy                 ! boundary update structure for the array

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      grid_loc                   ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   integer (int_kind), dimension(:,:,:,:,:), intent(inout) :: &
      ARRAY              ! array containing horizontal slab to update

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::           &
      k,l,m,n                      ! dummy loop indices

!-----------------------------------------------------------------------

   l = size(ARRAY,dim=3)
   n = size(ARRAY,dim=4)

   do k=1,l
   do m=1,n
      call boundary_2d_int(ARRAY(:,:,k,m,:),in_bndy,grid_loc,field_type)
   end do
   end do

!-----------------------------------------------------------------------

end subroutine boundary_4d_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: increment_message_counter
! !INTERFACE:

   subroutine increment_message_counter(snd_counter, rcv_counter, &
                                        src_proc, dst_proc)

! !DESCRIPTION:
!  This is a utility routine to increment the arrays for counting
!  whether messages are required.  It is used only for creating
!  boundary structures for updating ghost cells.

! !REVISION HISTORY:
!  Same as module.

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      src_proc,              & ! source processor for communication
      dst_proc                 ! destination processor for communication

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), dimension(:), intent(inout) :: &
      snd_counter,       &! array for counting messages to be sent
      rcv_counter         ! array for counting messages to be received

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if destination all land (proc = 0), then no send is necessary, 
!  so do the rest only for proc /= 0
!
!-----------------------------------------------------------------------

   if (dst_proc /= 0) then

      !*** if the current processor is the source, must send
      !*** data (local copy if dst_proc = src_proc)

      if (src_proc == my_task + 1) snd_counter(dst_proc) = &
                                   snd_counter(dst_proc) + 1

      !*** if the current processor is the destination, must
      !*** receive data (local copy if dst_proc = src_proc)

      if (dst_proc == my_task + 1) then

         if (src_proc /= 0) then  
            !*** the source block has ocean points so
            !*** increment the number of messages from
            !*** the source process

            rcv_counter(src_proc) = rcv_counter(src_proc) + 1

         else
            !*** the source block has no ocean points so
            !*** count this as a local copy in order to
            !*** fill ghost cells with zeroes

            rcv_counter(dst_proc) = rcv_counter(dst_proc) + 1

         endif
      endif
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine increment_message_counter

!***********************************************************************

end module boundary

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
