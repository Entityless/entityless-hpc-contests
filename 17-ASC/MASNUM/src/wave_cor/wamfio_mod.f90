!###############################################################################
!-------------------------------------------------------------------------------

  module wamfio_mod

!-------------------------------------------------------------------------------

  use time_mod
  use netcdf_mod
  use wamvar_mod
  use wamcpl_mod
  
  implicit none

!-------------------------------------------------------------------------------

  public wamfio_mod_init
  public my_getwind, get_wind, inprst, outrst, output
  private outwav_t, outmix, settopog
  private

!-------------------------------------------------------------------------------

  real, allocatable, dimension(:, :) :: windx1, windx2, windy1, windy2

  integer :: wflag = 0
  integer :: mod_init = 0
  integer :: recd = 0
  integer :: wflag_qbln = 0, wflag_ncep = 0
  double precision :: dwtime1 = -1.d0, dwtime2 = -1.d0

  real, allocatable :: v2(:, :), v3(:, :, :)
  integer, allocatable :: iv2(:, :)

  real, allocatable, dimension(:, :, :) :: wcoe
  integer, allocatable, dimension(:, :) :: widxx, widxy
  real, allocatable, dimension(:, :) :: wu, wv
  real, allocatable, dimension(:) :: lon, lat
  integer*2, allocatable :: ivar(:, :)

!-------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------

  include 'wamfio_mod_init.inc'
  include 'settopog.inc'
  include 'setwave.inc'
  include 'nlweight.inc'

  include 'get_wind.inc'

  include 'io_rest.inc'
  include 'output.inc'

  include 'outwav_t.inc'
  include 'outeng_t.inc'

  include 'outmix.inc'

  include 'outmix_tau.inc'
  include 'outmix_bv.inc'
  include 'outmix_bvl.inc'
 
  include 'outmix_wit.inc'  !shenhj 2012-09-23

!-------------------------------------------------------------------------------

  end module wamfio_mod

!-------------------------------------------------------------------------------
!###############################################################################
