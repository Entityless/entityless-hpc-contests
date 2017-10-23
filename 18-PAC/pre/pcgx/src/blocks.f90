!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module blocks

!BOP
! !MODULE: blocks
!
! !DESCRIPTION:
! This module contains data types and tools for decomposing a global
! horizontal domain into a set of blocks. It contains a data type
! for describing each block and contains routines for creating and
! querying the block decomposition for a global domain.
!
! !REVISION HISTORY:
! CVS:$Id: blocks.F90,v 1.10 2003/12/23 22:11:40 pwjones Exp $
! CVS:$Name: POP_2_0_1 $
!
! !USES:

   use kinds_mod
   use exit_mod
   use domain_size
   use communicate, only: my_task, master_task

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: block ! block data type
      integer (int_kind) :: &
         block_id ,&! global block number
         local_id ,&! local address of block in current distrib
         ib, ie, jb, je ,&! begin,end indices for physical domain
         iblock, jblock ! cartesian i,j position for bloc

      integer (int_kind), dimension(:), pointer :: &
         i_glob, j_glob ! global domain location for each point
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_blocks ,&
             create_blocks_mac ,&
             create_blocks_tiny ,&
             destroy_blocks ,&
             get_block ,&
             get_block_parameter

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: &
      nghost = 2 ! number of ghost cells around each block

   integer (int_kind), parameter, public :: &! size of block domain in
      nx_block = block_size_x + 2*nghost, &! x,y dir including ghost
      ny_block = block_size_y + 2*nghost ! cells

   integer (int_kind), parameter, public :: &
      nx_bkfactor = 1 ,&! size of macro-block(nx_bkfactor*nx_block,ny_bkfactor*ny_block)
      ny_bkfactor = 1 ! for cache blocking opt

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      nblocks_tot ,&! total number of blocks in decomposition
      nblocks_x ,&! tot num blocks in i direction
      nblocks_y ! tot num blocks in j direction

   logical (log_kind), public :: &!
      lmic_trace, &! flag to debug trace
      lmic_proc ! flag to MIC process using tiny blocks

   integer (int_kind), parameter, public :: &
      block_size_x_mac = block_size_x*nx_bkfactor, &! size of macro-block in first horizontal dimension
      block_size_y_mac = block_size_y*ny_bkfactor ! size of macro-block in second horizontal dimension

   integer (int_kind), parameter, public :: &! size of macro-block domain in
      nx_block_mac = block_size_x_mac + 2*nghost, &! x,y dir including ghost
      ny_block_mac = block_size_y_mac + 2*nghost ! cells

   integer (int_kind), public :: &
      nblocks_tot_tiny ! total number of tiny blocks, the macro-blocks are in the first
                          !

!EOP
!BOC
!-----------------------------------------------------------------------
!
! module private data
!
!-----------------------------------------------------------------------

   type (block), dimension(:), allocatable :: &
      all_blocks ! block information for all blocks in domain

   integer (int_kind), dimension(:,:), allocatable, target :: &
      i_global, &! global i index for each point in each block
      j_global, &! global j index for each point in each block
      i_global_tiny, &! suppl. global i index for each point in tiny block
      j_global_tiny ! suppl. global j index for each point in tiny block

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: create_blocks
! !INTERFACE:

 subroutine create_blocks(nx_global, ny_global, ew_boundary_type, &
                                                ns_boundary_type)

! !DESCRIPTION:
! This subroutine decomposes the global domain into blocks and
! fills the data structures with all the necessary block information.
!
! !REVISION HISTORY:
! same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nx_global, ny_global ! global domain size in x,y

   character (*), intent(in) :: &
      ew_boundary_type, &! type of boundary in logical east-west dir
      ns_boundary_type ! type of boundary in logical north-south dir

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, ip1, j, jp1, n ,&! loop indices
      iblock, jblock ,&! block loop indices
      is, ie, js, je ! temp start, end indices

!----------------------------------------------------------------------
!
! compute number of blocks and cartesian decomposition
! if the requested block size does not divide the global domain
! size evenly, add additional block space to accomodate padding
!
!----------------------------------------------------------------------

   nblocks_x = (nx_global-1)/block_size_x + 1
   nblocks_y = (ny_global-1)/block_size_y + 1
   nblocks_tot = nblocks_x*nblocks_y

!----------------------------------------------------------------------
!
! allocate block arrays
!
!----------------------------------------------------------------------

   allocate(all_blocks(nblocks_tot))
   allocate(i_global(nx_block,nblocks_tot), &
            j_global(ny_block,nblocks_tot))

!----------------------------------------------------------------------
!
! fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   n = 0
   do jblock=1,nblocks_y
      js = (jblock-1)*block_size_y + 1
      je = js + block_size_y - 1
      if (js > ny_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: ny_block too large?')
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nblocks_x
         n = n + 1 ! global block id

         is = (iblock-1)*block_size_x + 1
         ie = is + block_size_x - 1
         if (is > nx_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: nx_block too large?')
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = n
         all_blocks(n)%iblock = iblock
         all_blocks(n)%jblock = jblock
         all_blocks(n)%ib = nghost + 1
         all_blocks(n)%jb = nghost + 1
         all_blocks(n)%ie = nx_block - nghost ! default value
         all_blocks(n)%je = ny_block - nghost ! default value

         do j=1,ny_block
            j_global(j,n) = js - nghost + j - 1


            !*** southern ghost cells

            if (j_global(j,n) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) + ny_global
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select
            endif

            !*** padding required

            if (j_global(j,n) > ny_global + nghost) then
               j_global(j,n) = 0 ! padding

            !*** northern ghost cells

            else if (j_global(j,n) > ny_global) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) - ny_global
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = -j_global(j,n)
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global(j,n) == ny_global .and. &
                     j > all_blocks(n)%jb) then
               all_blocks(n)%je = j ! last physical point in padded domain
            endif
         end do

         all_blocks(n)%j_glob => j_global(:,n)

         do i=1,nx_block
            i_global(i,n) = is - nghost + i - 1

            !*** western ghost cells

            if (i_global(i,n) < 1) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) + nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (i_global(i,n) > nx_global + nghost) then
               i_global(i,n) = 0

            !*** eastern ghost cells

            else if (i_global(i,n) > nx_global) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) - nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global(i,n) == nx_global .and. &
                     i > all_blocks(n)%ib) then
               all_blocks(n)%ie = i
            endif
         end do

         all_blocks(n)%i_glob => i_global(:,n)

      end do
   end do

!EOC
!----------------------------------------------------------------------

end subroutine create_blocks

!***********************************************************************
!BOP
! !IROUTINE: create_blocks_mac
! !INTERFACE:

 subroutine create_blocks_mac(nx_global, ny_global, ew_boundary_type, &
                                                ns_boundary_type)

! !DESCRIPTION:
! This subroutine decomposes the global domain into blocks and
! fills the data structures with all the necessary block information.
!
! !REVISION HISTORY:
! same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nx_global, ny_global ! global domain size in x,y

   character (*), intent(in) :: &
      ew_boundary_type, &! type of boundary in logical east-west dir
      ns_boundary_type ! type of boundary in logical north-south dir

!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, ip1, j, jp1, n ,&! loop indices
      iblock, jblock ,&! block loop indices
      is, ie, js, je ! temp start, end indices

!----------------------------------------------------------------------
!
! compute number of blocks and cartesian decomposition
! if the requested block size does not divide the global domain
! size evenly, add additional block space to accomodate padding
!
!----------------------------------------------------------------------

   nblocks_x = (nx_global-1)/block_size_x_mac + 1
   nblocks_y = (ny_global-1)/block_size_y_mac + 1
   nblocks_tot = nblocks_x*nblocks_y

!----------------------------------------------------------------------
!
! allocate block arrays
!
!----------------------------------------------------------------------

   allocate(all_blocks(nblocks_tot))
   allocate(i_global(nx_block_mac,nblocks_tot), &
            j_global(ny_block_mac,nblocks_tot))

!----------------------------------------------------------------------
!
! fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   n = 0
   do jblock=1,nblocks_y
      js = (jblock-1)*block_size_y_mac + 1
      je = js + block_size_y_mac - 1
      if (js > ny_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: ny_block too large?')
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nblocks_x
         n = n + 1 ! global block id

         is = (iblock-1)*block_size_x_mac + 1
         ie = is + block_size_x_mac - 1
         if (is > nx_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: nx_block too large?')
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = n
         all_blocks(n)%iblock = iblock
         all_blocks(n)%jblock = jblock
         all_blocks(n)%ib = nghost + 1
         all_blocks(n)%jb = nghost + 1
         all_blocks(n)%ie = nx_block_mac - nghost ! default value
         all_blocks(n)%je = ny_block_mac - nghost ! default value

         do j=1,ny_block_mac
            j_global(j,n) = js - nghost + j - 1


            !*** southern ghost cells

            if (j_global(j,n) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) + ny_global
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select
            endif

            !*** padding required

            if (j_global(j,n) > ny_global + nghost) then
               j_global(j,n) = 0 ! padding

            !*** northern ghost cells

            else if (j_global(j,n) > ny_global) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) - ny_global
               case ('closed')
                  j_global(j,n) = 0
               case ('tripole')
                  j_global(j,n) = -j_global(j,n)
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global(j,n) == ny_global .and. &
                     j > all_blocks(n)%jb) then
               all_blocks(n)%je = j ! last physical point in padded domain
            endif
         end do

         all_blocks(n)%j_glob => j_global(:,n)

         do i=1,nx_block_mac
            i_global(i,n) = is - nghost + i - 1

            !*** western ghost cells

            if (i_global(i,n) < 1) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) + nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (i_global(i,n) > nx_global + nghost) then
               i_global(i,n) = 0

            !*** eastern ghost cells

            else if (i_global(i,n) > nx_global) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) - nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global(i,n) == nx_global .and. &
                     i > all_blocks(n)%ib) then
               all_blocks(n)%ie = i
            endif
         end do

         all_blocks(n)%i_glob => i_global(:,n)

      end do
   end do

!EOC
!----------------------------------------------------------------------

end subroutine create_blocks_mac

!***********************************************************************
!BOP
! !IROUTINE: create_blocks_mac
! !INTERFACE:

 subroutine create_blocks_tiny(nx_global, ny_global, &
                               ew_boundary_type, &
                               ns_boundary_type, &
                               nblocks_local, block_ids )
! !DESCRIPTION:
! This subroutine decomposes the global domain into blocks and
! fills the data structures with all the necessary block information.
!
! !REVISION HISTORY:
! same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      nx_global, ny_global ! global domain size in x,y

   character (*), intent(in) :: &
      ew_boundary_type, &! type of boundary in logical east-west dir
      ns_boundary_type ! type of boundary in logical north-south dir

   integer (int_kind), intent(in) :: &
      nblocks_local ! num of distrb blocks to local MIC process

   integer (int_kind), dimension(:), intent(in) :: &
      block_ids ! array of macro block ids for every block
                             !that resides on the local process
!EOP
!BOC
!----------------------------------------------------------------------
!
! local variables
!
!----------------------------------------------------------------------

   type (block), dimension(:), allocatable :: &
      backup_blocks ! block information for all blocks in domain

   integer (int_kind) :: &
      i, ip1, j, jn, n ,&! loop indices
      iblock, jblock ,&! block loop indices
      k ,&! normal block loop index
      block_id,base_i,base_j,&!
      is, ie, js, je ! temp start, end indices

   nblocks_tot_tiny = nblocks_local * nx_bkfactor * ny_bkfactor

!----------------------------------------------------------------------
!
! allocate block arrays
!
!----------------------------------------------------------------------

   allocate(backup_blocks(nblocks_tot))
   backup_blocks = all_blocks
   deallocate(all_blocks)
   allocate(all_blocks(nblocks_tot + nblocks_tot_tiny))
   all_blocks(1:nblocks_tot) = backup_blocks(1:nblocks_tot)
   deallocate(backup_blocks)
   ! TODO: let tiny blocks in the same row/column share the same {j,i}_glob
   ! array
   allocate(i_global_tiny(nx_block,nblocks_tot_tiny), &
            j_global_tiny(ny_block,nblocks_tot_tiny))
!----------------------------------------------------------------------
!
! fill block data structures for all tiny blocks in domain
! append to normal blocks
!----------------------------------------------------------------------

   n = nblocks_tot
   do k=1,nblocks_local
   block_id = block_ids(k)
   base_j = all_blocks(block_id)%j_glob(all_blocks(block_id)%jb)
   base_i = all_blocks(block_id)%i_glob(all_blocks(block_id)%ib)
! if(all_blocks(block_id)%ie < nx_block_mac-nghost .or. &
! all_blocks(block_id)%je < ny_block_mac-nghost ) &
! write(6,*) 'create_block_tiny:rootblock:',my_task, block_id,&
! all_blocks(block_id)%ib,all_blocks(block_id)%jb,&
! all_blocks(block_id)%ie,all_blocks(block_id)%je

   do jblock=1,ny_bkfactor
      js = (jblock-1)*block_size_y + base_j
      je = js + block_size_y - 1
      if (js > ny_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: ny_block too large?')
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nx_bkfactor
         n = n + 1 ! global block id

         is = (iblock-1)*block_size_x + base_i
         ie = is + block_size_x - 1
         if (is > nx_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: nx_block too large?')
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = block_id ! macro-block id
         all_blocks(n)%local_id = n - nblocks_tot ! tiny-block id
         all_blocks(n)%iblock = iblock
         all_blocks(n)%jblock = jblock
         all_blocks(n)%ib = nghost + 1
         all_blocks(n)%jb = nghost + 1
         all_blocks(n)%ie = nx_block - nghost ! default value
         all_blocks(n)%je = ny_block - nghost ! default value

         jn = n - nblocks_tot
         do j=1,ny_block
            j_global_tiny(j,jn) = js - nghost + j - 1


            !*** southern ghost cells

            if (j_global_tiny(j,jn) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global_tiny(j,jn) = j_global_tiny(j,jn) + ny_global
               case ('closed')
                  j_global_tiny(j,jn) = 0
               case ('tripole')
                  j_global_tiny(j,jn) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select
            endif

            !*** padding required

            if (j_global_tiny(j,jn) > ny_global + nghost) then
               j_global_tiny(j,jn) = 0 ! padding

            !*** northern ghost cells

            else if (j_global_tiny(j,jn) > ny_global) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global_tiny(j,jn) = j_global_tiny(j,jn) - ny_global
               case ('closed')
                  j_global_tiny(j,jn) = 0
               case ('tripole')
                  j_global_tiny(j,jn) = -j_global_tiny(j,jn)
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global_tiny(j,jn) == ny_global .and. &
                     j > all_blocks(n)%jb) then
               all_blocks(n)%je = j ! last physical point in padded domain
            endif
         end do

         all_blocks(n)%j_glob => j_global_tiny(:,jn)

         do i=1,nx_block
            i_global_tiny(i,jn) = is - nghost + i - 1

            !*** western ghost cells

            if (i_global_tiny(i,jn) < 1) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global_tiny(i,jn) = i_global_tiny(i,jn) + nx_global
               case ('closed')
                  i_global_tiny(i,jn) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (i_global_tiny(i,jn) > nx_global + nghost) then
               i_global_tiny(i,jn) = 0

            !*** eastern ghost cells

            else if (i_global_tiny(i,jn) > nx_global) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global_tiny(i,jn) = i_global_tiny(i,jn) - nx_global
               case ('closed')
                  i_global_tiny(i,jn) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global_tiny(i,jn) == nx_global .and. &
                     i > all_blocks(n)%ib) then
               all_blocks(n)%ie = i
            endif
         end do

         all_blocks(n)%i_glob => i_global_tiny(:,jn)

      end do ! loop iblock
   end do ! loop jblock
   end do ! loop k

!EOC
!----------------------------------------------------------------------

end subroutine create_blocks_tiny


!***********************************************************************
!BOP
! !IROUTINE: get_block
! !INTERFACE:

 function get_block(block_id,local_id)

! !DESCRIPTION:
! This function returns the block data structure for the block
! associated with the input block id.
!
! !REVISION HISTORY:
! same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block_id, &! global block id for requested block info
      local_id ! local block id to assign to this block

! !OUTPUT PARAMETERS:

   type (block) :: &
      get_block ! block information returned for requested block

!EOP
!BOC
!----------------------------------------------------------------------
!
! check for valid id. if valid, return block info for requested block
!
!----------------------------------------------------------------------
   if (.not. lmic_proc) then

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call exit_POP(sigAbort,'get_block: invalid block_id')
   endif

   else ! lmic_proc

   if (block_id < 1 .or. block_id > nblocks_tot+nblocks_tot_tiny) then
      call exit_POP(sigAbort,'get_block: invalid block_id')
   endif

   endif ! lmic_proc

   get_block = all_blocks(block_id)
   get_block%local_id = local_id
!----------------------------------------------------------------------
!EOC

 end function get_block

!**********************************************************************
!BOP
! !IROUTINE: get_block_parameter
! !INTERFACE:

 subroutine get_block_parameter(block_id, local_id, ib, ie, jb, je, &
                                iblock, jblock, i_glob, j_glob)

! !DESCRIPTION:
! This routine returns requested parts of the block data type
! for the block associated with the input block id
!
! !REVISION HISTORY:
! same as module
!
! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block_id ! global block id for which parameters are requested

! !OUTPUT PARAMETERS:

   !(optional) parts of block data type to extract if requested

   integer (int_kind), intent(out), optional :: &
      local_id ,&! local id assigned to block in current distrb
      ib, ie, jb, je ,&! begin,end indices for physical domain
      iblock, jblock ! cartesian i,j position for bloc

   integer (int_kind), dimension(:), pointer, optional :: &
      i_glob, j_glob ! global domain location for each point

!EOP
!BOC
!----------------------------------------------------------------------
!
! extract each component of data type if requested
!
!----------------------------------------------------------------------

   if (.not. lmic_proc) then

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call exit_POP(sigAbort,'get_block_parameter: invalid block_id')
   endif

   else ! lmic_proc

   if (block_id < 1 .or. block_id > nblocks_tot+nblocks_tot_tiny) then
      call exit_POP(sigAbort,'get_block_parameter: invalid block_id')
   endif

   endif ! lmic_proc

   if (present(local_id)) local_id = all_blocks(block_id)%local_id
   if (present(ib )) ib = all_blocks(block_id)%ib
   if (present(ie )) ie = all_blocks(block_id)%ie
   if (present(jb )) jb = all_blocks(block_id)%jb
   if (present(je )) je = all_blocks(block_id)%je
   if (present(iblock )) iblock = all_blocks(block_id)%iblock
   if (present(jblock )) jblock = all_blocks(block_id)%jblock
   if (present(i_glob )) i_glob = all_blocks(block_id)%i_glob
   if (present(j_glob )) j_glob = all_blocks(block_id)%j_glob

!----------------------------------------------------------------------
!EOC

 end subroutine get_block_parameter

!**********************************************************************
!BOP
! !IROUTINE: destroy_blocks
! !INTERFACE:

 subroutine destroy_blocks

! !DESCRIPTION:
! This subroutine deallocates the array with block information.
!
! !REVISION HISTORY:
! same as module
!EOP
!----------------------------------------------------------------------
!BOC

   deallocate(all_blocks)

!EOC
!----------------------------------------------------------------------

 end subroutine destroy_blocks

!***********************************************************************

 end module blocks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
