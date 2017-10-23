!###############################################################################
!-------------------------------------------------------------------------------

  module wamcpl_mod

!-------------------------------------------------------------------------------

  use wamvar_mod

  implicit none

!-------------------------------------------------------------------------------

  public intact, mixture, mixture_limit, set_ice, set_uv, set_uv0, init_wamcpl
  public zyyz, taubb11, taubb12, taubb22, taubb33
  public ux, uy, uxx, uxy, uyx, uyy
  public bv, bvl, bh1, bh2
  public bv_wtv, bv_wtd                   !shenhj 2012-09-22
  public mixture_wit                      !shenhj 2012-09-23

  private

!-------------------------------------------------------------------------------

! --- Background 2D currents: ux, uy
!     uxy=dux/dy, uxx=dux/dx, uyy=duy/dy, uyx=duy/dx

  real, allocatable :: ux(:, :)
  real, allocatable :: uy(:, :)
  real, allocatable :: uxx(:, :)
  real, allocatable :: uxy(:, :)
  real, allocatable :: uyx(:, :)
  real, allocatable :: uyy(:, :)

!-------------------------------------------------------------------------------

! --- For wave induced mixing & Reynolds stresses.

  real, allocatable :: zyyz(:)

! --- Wave induced mixing (m^2/s)
  real, pointer :: bv(:, :, :)

! --- Wave induced Reynold stresses (m^2/s^2)
  real, pointer :: taubb11(:, :, :)
  real, pointer :: taubb12(:, :, :)
  real, pointer :: taubb22(:, :, :)
  real, pointer :: taubb33(:, :, :)

  real, pointer :: bvl(:, :, :)
  real, pointer :: bh1(:, :, :)
  real, pointer :: bh2(:, :, :)
  
  real, allocatable :: bv_wtv(:, :, :)           !shenhj 2012-09-22
  real, allocatable :: bv_wtd(:, :, :)           !shenhj 2012-09-22

!-------------------------------------------------------------------------------

  integer :: init_ice = 0
  double precision :: last_time_ice = -99999.d0

  real, allocatable :: noicensp(:, :)
  real, allocatable :: icensp(:, :)

!-------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------

  include 'set_uv.inc'
  include 'intact.inc'

  include 'mixture.inc'
  include 'mixture_limit.inc'
  include 'mixture_wit.inc'              !shenhj 2012-09-22
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: init_wamcpl

  subroutine init_wamcpl(mixflag)

!-------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------

  integer, intent(in) :: mixflag

!-------------------------------------------------------------------------

! --- Background 2D currents: ux, uy
!     uxy=dux/dy, uxx=dux/dx, uyy=duy/dy, uyx=duy/dx

  allocate(ux(ixs:ixl, iys:iyl) )
  allocate(uy(ixs:ixl, iys:iyl) )
  allocate(uxx(ixs:ixl, iys:iyl))
  allocate(uxy(ixs:ixl, iys:iyl))
  allocate(uyx(ixs:ixl, iys:iyl))
  allocate(uyy(ixs:ixl, iys:iyl))

!-------------------------------------------------------------------------------
! --- For wave induced mixing & Reynolds stresses.

  allocate(zyyz(kb)             )

  if(mixflag == 1)then
    allocate(bv(ixs:ixl, iys:iyl, kb)     )
    allocate(taubb11(ixs:ixl, iys:iyl, kb))
    allocate(taubb12(ixs:ixl, iys:iyl, kb))
    allocate(taubb22(ixs:ixl, iys:iyl, kb))
    allocate(taubb33(ixs:ixl, iys:iyl, kb))
    allocate(ea(kl,jl,ixs:ixl, iys:iyl));   ea = 0.0
  elseif(mixflag == 2)then
    allocate(bvl(ixs:ixl, iys:iyl, kb)    )
    allocate(bh1(ixs:ixl, iys:iyl, kb)    )
    allocate(bh2(ixs:ixl, iys:iyl, kb)    )
    allocate(taubb11(ixs:ixl, iys:iyl, kb))
    allocate(taubb12(ixs:ixl, iys:iyl, kb))
    allocate(taubb22(ixs:ixl, iys:iyl, kb))
    allocate(taubb33(ixs:ixl, iys:iyl, kb))
    allocate(ea(kl,jl,ixs:ixl, iys:iyl));   ea = 0.0
  elseif(mixflag >= 3)then
    allocate(bv(ixs:ixl, iys:iyl, kb)     )
    allocate(bvl(ixs:ixl, iys:iyl, kb)    )
    allocate(bh1(ixs:ixl, iys:iyl, kb)    )
    allocate(bh2(ixs:ixl, iys:iyl, kb)    )
    allocate(taubb11(ixs:ixl, iys:iyl, kb))
    allocate(taubb12(ixs:ixl, iys:iyl, kb))
    allocate(taubb22(ixs:ixl, iys:iyl, kb))
    allocate(taubb33(ixs:ixl, iys:iyl, kb))
    allocate(ea(kl,jl,ixs:ixl, iys:iyl));   ea = 0.0
    allocate(bv_wtv(ixs:ixl, iys:iyl, kb) )    !shenhj 2012-09-22
    allocate(bv_wtd(ixs:ixl, iys:iyl, kb) )    !shenhj 2012-09-22
  endif

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------

  end subroutine init_wamcpl

!-------------------------------------------------------------------------------

  end module wamcpl_mod

!-------------------------------------------------------------------------------
!###############################################################################
