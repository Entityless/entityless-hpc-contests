!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module domain_size

!BOP
! !MODULE: domain_size
!
! !DESCRIPTION:
! This module contains parameters for the global model domain size
! decomposition block size. It is used by the domain and block
! modules for decomposing the model domain across processors.
!
! !REVISION HISTORY:
! CVS:$Id: domain.F90,v 1.16 2002/11/18 17:35:33 jfd Exp $
! CVS:$Name: $

! !USES:

   use kinds_mod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: & ! model size parameters
      nx_global = 3600 ,&! extent of horizontal axis in i direction
      ny_global = 2400 ,&! extent of horizontal axis in j direction
      km = 40 ,&! number of vertical levels
      nt = 2 ! total number of tracers

   integer (int_kind), parameter, public :: &
      block_size_x = 180, &! size of block in first horizontal dimension
      block_size_y = 100 ! size of block in second horizontal dimension

   !*** The model will inform the user of the correct
   !*** values for theparameters below. A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary. A value that
   !*** is too low will cause the code to exit.
   !*** A good initial guess is found using
   !*** max=(nx_global/block_size_x)*(ny_global/block_size_y)/
   !*** num_procs

   integer (int_kind), parameter, public :: &
      max_blocks_clinic = 9, &! max number of blocks per processor
      max_blocks_tropic = 9 ! in each distribution

!EOP
!BOC
!EOC
!***********************************************************************

 end module domain_size

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
