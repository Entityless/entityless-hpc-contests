!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: gather_scatter

 module gather_scatter

! !DESCRIPTION:
!  This module contains routines for gathering data to a single
!  processor from a distributed array, scattering data from a
!  single processor to a distributed array and changing distribution
!  of blocks of data (eg from baroclinic to barotropic and back).
!
! !REVISION HISTORY:
!  CVS: $Id: gather_scatter.F90,v 1.6 2003/12/23 22:35:00 pwjones Exp $
!  CVS: $Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use communicate
   use constants
   use blocks
   use distribution
   use domain
   use domain_size
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: gather_global,      &
             scatter_global,     & 
             redistribute_blocks

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface gather_global
     module procedure gather_global_dbl,  &
                      gather_global_real, &
                      gather_global_int
   end interface 

   interface scatter_global
     module procedure scatter_global_dbl,  &
                      scatter_global_real, &
                      scatter_global_int
   end interface 

   interface redistribute_blocks
     module procedure redistribute_blocks_dbl,  &
                      redistribute_blocks_real, &
                      redistribute_blocks_int
   end interface 

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: gather_global
! !INTERFACE:

 subroutine gather_global_dbl(ARRAY_G, ARRAY, dst_task, src_dist)

! !DESCRIPTION:
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific inteface for double precision arrays 
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).


! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (r8), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r8), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

   if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_dbl, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(src_dist%local_block_num), &
              snd_status (MPI_STATUS_SIZE, src_dist%local_block_num))

     nsends = 0
     do n=1,src_dist%local_block_num
         src_block = n
         nsends = nsends + 1
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_dbl, dst_task, 3*mpitag_gs+src_dist%local_block_ids(n), &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

  else  ! lmic_proc
     call gather_global_dbl_tiny(ARRAY_G, ARRAY, dst_task, src_dist)
     !write(6,*) 'gather_global_dbl_tiny:done',my_task
  endif ! lmic_proc
!-----------------------------------------------------------------------

 end subroutine gather_global_dbl

!***********************************************************************
!BOP
! !IROUTINE: gather_global
! !INTERFACE:

 subroutine gather_global_dbl_tiny(ARRAY_G, ARRAY, dst_task, src_dist)

! !DESCRIPTION:
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific inteface for double precision arrays 
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).


! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (r8), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind) :: &
     k,local_id         ,&! dummy loop indices
     jblock,iblock      ,&! dummy loop indices
     js,je,is,ie          ! 

   real (r8), dimension(:,:,:), allocatable :: &
     MACRO_ARRAY      ! buffer for receiving blocks

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r8), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

   ! tuning debug
   call exit_POP(sigAbort,'gather_global_dbl_tiny:MIC is dst_task.')

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_dbl, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     !allocate(snd_request(src_dist%local_block_num), &
     allocate(snd_request(nblocks_clinic_mac), &
              snd_status (MPI_STATUS_SIZE, nblocks_clinic_mac))
     allocate(MACRO_ARRAY(nx_block_mac,ny_block_mac,nblocks_clinic_mac))

     local_id = 0
     do k=1,nblocks_clinic_mac
       do jblock=1,ny_bkfactor
       js = (jblock-1)*block_size_y + 1
       je = js + ny_block - 1
       do iblock=1,nx_bkfactor
          local_id = local_id + 1
          is = (iblock-1)*block_size_x + 1
          ie = is + nx_block - 1
          MACRO_ARRAY(is:ie,js:je,k) = ARRAY(1:nx_block,1:ny_block,local_id)
       end do
       end do
     end do

     nsends = 0
     local_id = 0
     do k=1,src_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(src_dist%local_block_ids(k),k) ! tiny-block
         n = this_block%block_id ! id of root macro-block

         nsends = nsends + 1
         call MPI_ISEND(MACRO_ARRAY(1,1,local_id), nx_block_mac*ny_block_mac, &
                     mpi_dbl, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)

     deallocate(MACRO_ARRAY)

     deallocate(snd_request, snd_status)

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_dbl_tiny


!***********************************************************************

 subroutine gather_global_real(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   real (r4), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r4), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

   if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_real, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,src_dist%local_block_num

         nsends = nsends + 1
         src_block = n
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_real, dst_task, 3*mpitag_gs+src_dist%local_block_ids(n), &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

  else  ! lmic_proc
     call gather_global_real_tiny(ARRAY_G, ARRAY, dst_task, src_dist)
     !write(6,*) 'gather_global_real_tiny:done',my_task
  endif ! lmic_proc
!-----------------------------------------------------------------------

 end subroutine gather_global_real

!***********************************************************************

 subroutine gather_global_real_tiny(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   real (r4), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind) :: &
     k,local_id         ,&! dummy loop indices
     jblock,iblock      ,&! dummy loop indices
     js,je,is,ie          ! 

   real (r4), dimension(:,:,:), allocatable :: &
     MACRO_ARRAY      ! buffer for receiving blocks

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r4), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

   ! tuning debug
   call exit_POP(sigAbort,'gather_global_real_tiny:MIC is dst_task.')

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_real, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     !allocate(snd_request(src_dist%local_block_num), &
     allocate(snd_request(nblocks_clinic_mac), &
              snd_status (MPI_STATUS_SIZE, nblocks_clinic_mac))
     allocate(MACRO_ARRAY(nx_block_mac,ny_block_mac,nblocks_clinic_mac))

     local_id = 0
     do k=1,nblocks_clinic_mac
       do jblock=1,ny_bkfactor
       js = (jblock-1)*block_size_y + 1
       je = js + ny_block - 1
       do iblock=1,nx_bkfactor
          local_id = local_id + 1
          is = (iblock-1)*block_size_x + 1
          ie = is + nx_block - 1
          MACRO_ARRAY(is:ie,js:je,k) = ARRAY(1:nx_block,1:ny_block,local_id)
       end do
       end do
     end do

     nsends = 0
     local_id = 0
     do k=1,src_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(src_dist%local_block_ids(k),k) ! tiny-block
         n = this_block%block_id ! id of root macro-block

         nsends = nsends + 1
         call MPI_ISEND(MACRO_ARRAY(1,1,local_id), nx_block_mac*ny_block_mac, &
                     mpi_real, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)

     deallocate(MACRO_ARRAY)

     deallocate(snd_request, snd_status)

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_real_tiny


!***********************************************************************

 subroutine gather_global_int(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_integer, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,src_dist%local_block_num

         nsends = nsends + 1
         src_block = n
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_integer, dst_task, 3*mpitag_gs+src_dist%local_block_ids(n), &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_int

!***********************************************************************

 subroutine gather_global_int_tiny(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind) :: &
     k,local_id         ,&! dummy loop indices
     jblock,iblock      ,&! dummy loop indices
     js,je,is,ie          ! 

   integer (int_kind), dimension(:,:,:), allocatable :: &
     MACRO_ARRAY      ! buffer for receiving blocks

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

   ! tuning debug
   call exit_POP(sigAbort,'gather_global_int_tiny:MIC is dst_task.')

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_integer, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     !allocate(snd_request(src_dist%local_block_num), &
     allocate(snd_request(nblocks_clinic_mac), &
              snd_status (MPI_STATUS_SIZE, nblocks_clinic_mac))
     allocate(MACRO_ARRAY(nx_block_mac,ny_block_mac,nblocks_clinic_mac))

     local_id = 0
     do k=1,nblocks_clinic_mac
       do jblock=1,ny_bkfactor
       js = (jblock-1)*block_size_y + 1
       je = js + ny_block - 1
       do iblock=1,nx_bkfactor
          local_id = local_id + 1
          is = (iblock-1)*block_size_x + 1
          ie = is + nx_block - 1
          MACRO_ARRAY(is:ie,js:je,k) = ARRAY(1:nx_block,1:ny_block,local_id)
       end do
       end do
     end do

     nsends = 0
     local_id = 0
     do k=1,src_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(src_dist%local_block_ids(k),k) ! tiny-block
         n = this_block%block_id ! id of root macro-block

         nsends = nsends + 1
         call MPI_ISEND(MACRO_ARRAY(1,1,local_id), nx_block_mac*ny_block_mac, &
                     mpi_integer, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)

     deallocate(MACRO_ARRAY)

     deallocate(snd_request, snd_status)

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_int_tiny

!EOC
!***********************************************************************
!BOP
! !IROUTINE: scatter_global
! !INTERFACE:

 subroutine scatter_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

! !DESCRIPTION:
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.  Note that this routine only
!  is guaranteed to scatter correct values in the physical domain
!  on each block.  Ghost cells are not filled correctly for tripole
!  boundary conditions due to the complexity of dealing with field
!  locations.  It is wise to call the boundary update routine after
!  each scatter call.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r8), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (r8), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

  if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         !ljm tuning
         !write(6,*) 'scatter_dbl:0sending:',dst_dist%proc(n)-1,n,nx_block,ny_block
         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_dbl, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,dst_dist%local_block_num
         dst_block = n
         this_block = get_block(dst_dist%local_block_ids(n),n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     !allocate (rcv_request(nblocks_tot), &
     !          rcv_status(MPI_STATUS_SIZE, nblocks_tot))
     allocate (rcv_request(dst_dist%local_block_num), &
               rcv_status(MPI_STATUS_SIZE, dst_dist%local_block_num))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,dst_dist%local_block_num
         nrecvs = nrecvs + 1
         dst_block = n
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_dbl, src_task, 3*mpitag_gs+dst_dist%local_block_ids(n), &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
     end do
     !write(6,*) 'scatter_dbl:MPI_IRECV:',my_task,nrecvs !dst_dist%local_block_ids(n)

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)
     !write(6,*) 'scatter_dbl:MPI_WAITALL:',my_task

     deallocate(rcv_request, rcv_status)
   endif

  else  ! lmic_proc
     call scatter_global_dbl_tiny(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)
  endif ! lmic_proc
!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl

!***********************************************************************
!BOP
! !IROUTINE: scatter_global
! !INTERFACE:

 subroutine scatter_global_dbl_tiny(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

! !DESCRIPTION:
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.  Note that this routine only
!  is guaranteed to scatter correct values in the physical domain
!  on each block.  Ghost cells are not filled correctly for tripole
!  boundary conditions due to the complexity of dealing with field
!  locations.  It is wise to call the boundary update routine after
!  each scatter call.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r8), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     root_block,&! block info for root macro block
     this_block  ! block info for current block

   integer (int_kind) :: &
     k,local_id         ,&! dummy loop indices
     jblock,iblock      ,&! dummy loop indices
     js,je,is,ie          ! 

   real (r8), dimension(:,:,:), allocatable :: &
     MACRO_ARRAY      ! buffer for receiving blocks

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (r8), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

   !ljm tuning
   call exit_POP(sigAbort, 'scatter_global_dbl_tiny:MIC is NOT assumed to be sender')
     !*** send non-local blocks away

     allocate (msg_buffer(nx_block_mac,ny_block_mac))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n) ! look at the macro block

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block_mac
            do i=1,nx_block_mac
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block_mac
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block_mac
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block_mac*ny_block_mac, &
                       mpi_dbl, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks
     ! tiny blocks
     local_id = 0
     do k=1,dst_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(dst_dist%local_block_ids(k),k) ! tiny-block
         root_block = get_block(this_block%block_id,local_id) ! root macro-block

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (root_block%iblock > 1         .and. &
             root_block%iblock < nblocks_x .and. &
             root_block%jblock > 1         .and. &
             root_block%jblock < nblocks_y) then

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do
           end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (root_block%jblock /= nblocks_y) then

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do
           end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do
           end do

         endif
     end do ! do k

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     !allocate (rcv_request(dst_dist%local_block_num), &
     allocate (rcv_request(nblocks_clinic_mac), &
               rcv_status(MPI_STATUS_SIZE, nblocks_clinic_mac))
     allocate(MACRO_ARRAY(nx_block_mac,ny_block_mac,nblocks_clinic_mac))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     local_id = 0
     do k=1,dst_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(dst_dist%local_block_ids(k),k) ! tiny-block
         n = this_block%block_id ! id of root macro-block

         nrecvs = nrecvs + 1
         call MPI_IRECV(MACRO_ARRAY(1,1,local_id), nx_block_mac*ny_block_mac, &
                       mpi_dbl, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
     end do
         !write(6,*) 'scatter_dbl_tiny:MPI_IRECV:',my_task,n,local_id

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)
         !write(6,*) 'scatter_dbl_tiny:MPI_WAITALL:',my_task

     local_id = 0
     do k=1,nblocks_clinic_mac
       do jblock=1,ny_bkfactor
       js = (jblock-1)*block_size_y + 1
       je = js + ny_block - 1
       do iblock=1,nx_bkfactor
          local_id = local_id + 1
          is = (iblock-1)*block_size_x + 1
          ie = is + nx_block - 1
          ARRAY(1:nx_block,1:ny_block,local_id) = MACRO_ARRAY(is:ie,js:je,k)
       end do
       end do
     end do
     deallocate(MACRO_ARRAY)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl_tiny


!***********************************************************************

 subroutine scatter_global_real(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r4), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (r4), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

  if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_real, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,dst_dist%local_block_num
         dst_block = n
         this_block = get_block(dst_dist%local_block_ids(n),n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,dst_dist%local_block_num
         nrecvs = nrecvs + 1
         dst_block = n
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_real, src_task, 3*mpitag_gs+dst_dist%local_block_ids(n), &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

  else  ! lmic_proc
     call scatter_global_real_tiny(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)
  endif ! lmic_proc
!-----------------------------------------------------------------------

 end subroutine scatter_global_real

!***********************************************************************

 subroutine scatter_global_real_tiny(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r4), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     root_block,&! block info for root macro block
     this_block  ! block info for current block

   integer (int_kind) :: &
     k,local_id         ,&! dummy loop indices
     jblock,iblock      ,&! dummy loop indices
     js,je,is,ie          ! 

   real (r4), dimension(:,:,:), allocatable :: &
     MACRO_ARRAY      ! buffer for receiving blocks

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (r4), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block_mac,ny_block_mac))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block_mac
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block_mac
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block_mac*ny_block_mac, &
                       mpi_real, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     ! tiny blocks
     local_id = 0
     do k=1,dst_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(dst_dist%local_block_ids(k),k) ! tiny-block
         root_block = get_block(this_block%block_id,local_id) ! root macro-block

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (root_block%iblock > 1         .and. &
             root_block%iblock < nblocks_x .and. &
             root_block%jblock > 1         .and. &
             root_block%jblock < nblocks_y) then

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do
           end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do
           end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do
           end do

         endif
     end do ! do k

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     !allocate (rcv_request(dst_dist%local_block_num), &
     allocate (rcv_request(nblocks_clinic_mac), &
               rcv_status(MPI_STATUS_SIZE, nblocks_clinic_mac))
     allocate(MACRO_ARRAY(nx_block_mac,ny_block_mac,nblocks_clinic_mac))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     local_id = 0
     do k=1,dst_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(dst_dist%local_block_ids(k),k) ! tiny-block
         n = this_block%block_id ! id of root macro-block

         nrecvs = nrecvs + 1
         call MPI_IRECV(MACRO_ARRAY(1,1,local_id), nx_block_mac*ny_block_mac, &
                       mpi_real, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     local_id = 0
     do k=1,nblocks_clinic_mac
       do jblock=1,ny_bkfactor
       js = (jblock-1)*block_size_y + 1
       je = js + ny_block - 1
       do iblock=1,nx_bkfactor
          local_id = local_id + 1
          is = (iblock-1)*block_size_x + 1
          ie = is + nx_block - 1
          ARRAY(1:nx_block,1:ny_block,local_id) = MACRO_ARRAY(is:ie,js:je,k)
       end do
       end do
     end do
     deallocate(MACRO_ARRAY)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_real_tiny


!***********************************************************************

 subroutine scatter_global_int(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

   if (.not. lmic_proc) then
!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = 0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_integer, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,dst_dist%local_block_num
         dst_block = n
         this_block = get_block(dst_dist%local_block_ids(n),n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,dst_dist%local_block_num
         nrecvs = nrecvs + 1
         dst_block = n
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_integer, src_task, 3*mpitag_gs+dst_dist%local_block_ids(n), &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

!     if (my_task == nproc_cpu_pn*nnode) &
!         write(6,*) 'KMT(macro):sum:',my_task,sum(ARRAY(:,:,1))

     deallocate(rcv_request, rcv_status)
   endif

  else  ! lmic_proc
     call scatter_global_int_tiny(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)
  endif ! lmic_proc
!-----------------------------------------------------------------------

 end subroutine scatter_global_int

!***********************************************************************

 subroutine scatter_global_int_tiny(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     root_block,&! block info for root macro block
     this_block  ! block info for current block

   integer (int_kind) :: &
     k,local_id         ,&! dummy loop indices
     jblock,iblock      ,&! dummy loop indices
     js,je,is,ie          ! 

   integer (int_kind), dimension(:,:,:), allocatable :: &
     MACRO_ARRAY      ! buffer for receiving blocks

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block_mac,ny_block_mac))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n) ! look at the macro block

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block_mac
            do i=1,nx_block_mac
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block_mac
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block_mac
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block_mac
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block_mac*ny_block_mac, &
                       mpi_integer, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks
     ! tiny blocks
     local_id = 0
     do k=1,dst_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(dst_dist%local_block_ids(k),k) ! tiny-block
         root_block = get_block(this_block%block_id,local_id) ! root macro-block

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (root_block%iblock > 1         .and. &
             root_block%iblock < nblocks_x .and. &
             root_block%jblock > 1         .and. &
             root_block%jblock < nblocks_y) then

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do
           end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (root_block%jblock /= nblocks_y) then

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do
           end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

           do n=k,k+nx_bkfactor*ny_bkfactor-1
            dst_block = n
            this_block = get_block(dst_dist%local_block_ids(n),n) ! tiny-block
            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do
           end do

         endif
     end do ! do k

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     !allocate (rcv_request(dst_dist%local_block_num), &
     allocate (rcv_request(nblocks_clinic_mac), &
               rcv_status(MPI_STATUS_SIZE, nblocks_clinic_mac))
     allocate(MACRO_ARRAY(nx_block_mac,ny_block_mac,nblocks_clinic_mac))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     local_id = 0
     do k=1,dst_dist%local_block_num,nx_bkfactor*ny_bkfactor
         local_id = local_id+1
         this_block = get_block(dst_dist%local_block_ids(k),k) ! tiny-block
         n = this_block%block_id ! id of root macro-block

         nrecvs = nrecvs + 1
         call MPI_IRECV(MACRO_ARRAY(1,1,local_id), nx_block_mac*ny_block_mac, &
                       mpi_integer, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
     end do
         !write(6,*) 'scatter_int_tiny:MPI_IRECV:',my_task,n,local_id

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)
         !write(6,*) 'scatter_int_tiny:MPI_WAITALL:',my_task

!     if (my_task == nproc_cpu_pn*nnode) &
!         write(6,*) 'KMT(tiny):sum:',my_task,sum(MACRO_ARRAY(:,:,1))

     local_id = 0
     do k=1,nblocks_clinic_mac
       do jblock=1,ny_bkfactor
       js = (jblock-1)*block_size_y + 1
       je = js + ny_block - 1
       do iblock=1,nx_bkfactor
          local_id = local_id + 1
          is = (iblock-1)*block_size_x + 1
          ie = is + nx_block - 1
          ARRAY(1:nx_block,1:ny_block,local_id) = MACRO_ARRAY(is:ie,js:je,k)
       end do
       end do
     end do
     deallocate(MACRO_ARRAY)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_int_tiny

!EOC
!***********************************************************************
!BOP
! !IROUTINE: redistribute_blocks
! !INTERFACE:

 subroutine redistribute_blocks_dbl(DST_ARRAY, dst_dist, &
                                    SRC_ARRAY, src_dist)

! !DESCRIPTION:
!  This subroutine converts an array distributed in a one decomposition
!  an array distributed in a different decomposition
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
     src_dist    ,&! info on distribution of blocks for source array
     dst_dist      ! info on distribution of blocks for dest   array

   real (r8), dimension(:,:,:), intent(in) :: &
     SRC_ARRAY     ! array containing field in source distribution

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
     DST_ARRAY     ! array containing field in dest distribution

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,n, &
     src_task,       &! task where the block currently resides
     dst_task,       &! task where the block needs to end up
     src_blck,       &! blck where the block currently resides
     dst_blck,       &! blck where the block needs to end up
     num_sends,      &! number of messages sent from this task
     num_recvs,      &! number of messages received by this task
     ierr             ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request,    &! request array for receives
     snd_request      ! request array for sends

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status,     &! status array for receives
     snd_status       ! status array for sends

!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcv_request(nblocks_tot), &
             snd_request(nblocks_tot), &
             rcv_status(MPI_STATUS_SIZE, nblocks_tot), &
             snd_status(MPI_STATUS_SIZE, nblocks_tot))

   rcv_request = 0
   snd_request = 0
   rcv_status = 0
   snd_status = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   num_recvs = 0
   do n=1,dst_dist%local_block_num
     src_task = src_dist%proc(dst_dist%local_block_ids(n)) - 1
     if (src_task /= my_task) then

       num_recvs = num_recvs + 1
       dst_blck = n

       call MPI_IRECV(DST_ARRAY(1,1,dst_blck), nx_block*ny_block, &
                      mpi_dbl, src_task, 3*mpitag_gs+dst_dist%local_block_ids(n), &
                      MPI_COMM_OCN, rcv_request(num_recvs), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  now determine which sends are required and post the sends
!
!-----------------------------------------------------------------------

   num_sends = 0
   do n=1,src_dist%local_block_num

     dst_task = dst_dist%proc(src_dist%local_block_ids(n)) - 1

     if (dst_task /= my_task) then

       num_sends = num_sends + 1
       src_blck = n

       call MPI_ISEND(SRC_ARRAY(1,1,src_blck), nx_block*ny_block, &
                      mpi_dbl, dst_task, 3*mpitag_gs+src_dist%local_block_ids(n), &
                      MPI_COMM_OCN, snd_request(num_sends), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  if blocks are local, simply copy the proper buffers
!
!-----------------------------------------------------------------------

   do n=1,src_dist%local_block_num
     i = src_dist%local_block_ids(n)
     if (dst_dist%proc(i) == my_task+1) then

       DST_ARRAY(:,:,dst_dist%local_block(i)) = &
       SRC_ARRAY(:,:,n)

     endif
   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (num_sends /= 0) &
     call MPI_WAITALL(num_sends, snd_request, snd_status, ierr)
   if (num_recvs /= 0) &
     call MPI_WAITALL(num_recvs, rcv_request, rcv_status, ierr)

   deallocate (rcv_request, snd_request, rcv_status, snd_status)

!-----------------------------------------------------------------------

 end subroutine redistribute_blocks_dbl

!***********************************************************************

 subroutine redistribute_blocks_real(DST_ARRAY, dst_dist, &
                                     SRC_ARRAY, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine converts an array distributed in a baroclinic
!  data decomposition to an array in a barotropic decomposition
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(in) :: &
     src_dist      ! info on distribution of blocks for source array

   real (r4), dimension(:,:,:), intent(in) :: &
     SRC_ARRAY     ! array containing field in source distribution

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(inout) :: &
     dst_dist      ! info on dist of blocks for destination array

   real (r4), dimension(:,:,:), intent(inout) :: &
     DST_ARRAY     ! array containing field in dest distribution

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,n, &
     src_task,       &! task where the block currently resides
     dst_task,       &! task where the block needs to end up
     src_blck,       &! blck where the block currently resides
     dst_blck,       &! blck where the block needs to end up
     num_sends,      &! number of messages sent from this task
     num_recvs,      &! number of messages received by this task
     ierr             ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request,    &! request array for receives
     snd_request      ! request array for sends

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status,     &! status array for receives
     snd_status       ! status array for sends

!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcv_request(nblocks_tot), &
             snd_request(nblocks_tot), &
             rcv_status(MPI_STATUS_SIZE, nblocks_tot), &
             snd_status(MPI_STATUS_SIZE, nblocks_tot))

   rcv_request = 0
   snd_request = 0
   rcv_status = 0
   snd_status = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   num_recvs = 0
   do n=1,dst_dist%local_block_num
     src_task = src_dist%proc(dst_dist%local_block_ids(n)) - 1
     if (src_task /= my_task) then

       num_recvs = num_recvs + 1
       dst_blck = n

       call MPI_IRECV(DST_ARRAY(1,1,dst_blck), nx_block*ny_block, &
                      mpi_real, src_task, 3*mpitag_gs+dst_dist%local_block_ids(n), &
                      MPI_COMM_OCN, rcv_request(num_recvs), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  now determine which sends are required and post the sends
!
!-----------------------------------------------------------------------

   num_sends = 0
   do n=1,src_dist%local_block_num

     dst_task = dst_dist%proc(src_dist%local_block_ids(n)) - 1

     if (dst_task /= my_task) then

       num_sends = num_sends + 1
       src_blck = n

       call MPI_ISEND(SRC_ARRAY(1,1,src_blck), nx_block*ny_block, &
                      mpi_real, dst_task, 3*mpitag_gs+src_dist%local_block_ids(n), &
                      MPI_COMM_OCN, snd_request(num_sends), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  if blocks are local, simply copy the proper buffers
!
!-----------------------------------------------------------------------

   do n=1,src_dist%local_block_num
     i = src_dist%local_block_ids(n)
     if (dst_dist%proc(i) == my_task+1) then

       DST_ARRAY(:,:,dst_dist%local_block(i)) = &
       SRC_ARRAY(:,:,n)

     endif
   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (num_sends /= 0) &
     call MPI_WAITALL(num_sends, snd_request, snd_status, ierr)
   if (num_recvs /= 0) &
     call MPI_WAITALL(num_recvs, rcv_request, rcv_status, ierr)

   deallocate (rcv_request, snd_request, rcv_status, snd_status)

!-----------------------------------------------------------------------

 end subroutine redistribute_blocks_real

!***********************************************************************

 subroutine redistribute_blocks_int(DST_ARRAY, dst_dist, &
                                    SRC_ARRAY, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine converts an array distributed in a baroclinic
!  data decomposition to an array in a barotropic decomposition
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(in) :: &
     src_dist      ! info on distribution of blocks for source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     SRC_ARRAY     ! array containing field in source distribution

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(inout) :: &
     dst_dist      ! info on dist of blocks for destination array

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     DST_ARRAY     ! array containing field in dest distribution

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,n, &
     src_task,       &! task where the block currently resides
     dst_task,       &! task where the block needs to end up
     src_blck,       &! blck where the block currently resides
     dst_blck,       &! blck where the block needs to end up
     num_sends,      &! number of messages sent from this task
     num_recvs,      &! number of messages received by this task
     ierr             ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request,    &! request array for receives
     snd_request      ! request array for sends

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status,     &! status array for receives
     snd_status       ! status array for sends

!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcv_request(nblocks_tot), &
             snd_request(nblocks_tot), &
             rcv_status(MPI_STATUS_SIZE, nblocks_tot), &
             snd_status(MPI_STATUS_SIZE, nblocks_tot))

   rcv_request = 0
   snd_request = 0
   rcv_status = 0
   snd_status = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   num_recvs = 0
   do n=1,dst_dist%local_block_num
     src_task = src_dist%proc(dst_dist%local_block_ids(n)) - 1
     if (src_task /= my_task) then

       num_recvs = num_recvs + 1
       dst_blck = n

       call MPI_IRECV(DST_ARRAY(1,1,dst_blck), nx_block*ny_block, &
                      mpi_integer, src_task, 3*mpitag_gs+dst_dist%local_block_ids(n), &
                      MPI_COMM_OCN, rcv_request(num_recvs), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  now determine which sends are required and post the sends
!
!-----------------------------------------------------------------------

   num_sends = 0
   do n=1,src_dist%local_block_num

     dst_task = dst_dist%proc(src_dist%local_block_ids(n)) - 1

     if (dst_task /= my_task) then

       num_sends = num_sends + 1
       src_blck = n

       call MPI_ISEND(SRC_ARRAY(1,1,src_blck), nx_block*ny_block, &
                      mpi_integer, dst_task, 3*mpitag_gs+src_dist%local_block_ids(n), &
                      MPI_COMM_OCN, snd_request(num_sends), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  if blocks are local, simply copy the proper buffers
!
!-----------------------------------------------------------------------

   do n=1,src_dist%local_block_num
     i = src_dist%local_block_ids(n)
     if (dst_dist%proc(i) == my_task+1) then

       DST_ARRAY(:,:,dst_dist%local_block(i)) = &
       SRC_ARRAY(:,:,n)

     endif
   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (num_sends /= 0) &
     call MPI_WAITALL(num_sends, snd_request, snd_status, ierr)
   if (num_recvs /= 0) &
     call MPI_WAITALL(num_recvs, rcv_request, rcv_status, ierr)

   deallocate (rcv_request, snd_request, rcv_status, snd_status)

!-----------------------------------------------------------------------

 end subroutine redistribute_blocks_int

!EOC
!***********************************************************************

 end module gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
