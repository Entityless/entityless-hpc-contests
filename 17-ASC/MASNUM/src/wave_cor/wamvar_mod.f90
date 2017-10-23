!###############################################################################
!-------------------------------------------------------------------------------

  module wamvar_mod

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  public wamvar_mod_init

  public

!-------------------------------------------------------------------------------

  integer, allocatable :: pebox(:, :);
  integer, allocatable :: peboxex(:, :);
  integer :: ixs
  integer :: ixl
  integer :: ix1
  integer :: ix2

  integer :: iys
  integer :: iyl
  integer :: iy1
  integer :: iy2
  
  integer :: my_hor_start, my_hor_end, my_hor_len
  integer :: my_ver_start, my_ver_end, my_ver_len
  integer :: ver_inner_start, ver_inner_end, hor_inner_start, hor_inner_end

  integer :: kb

! Saved index and intermediate values
! huangh223, 2017-02-25
  integer :: has_c_inited = 0
  integer :: has_propagat_inited = 0
  integer :: has_c_acce_kernel_inited = 0
  integer, allocatable :: idxs(:, :, :, :, :)
  real, allocatable :: tmp_values(:, :, :, :, :)
  
!-------------------------------------------------------------------------------

  integer, parameter :: kl    = 25
  integer, parameter :: kld   = 30
  integer, parameter :: klp1  = kl + 1
  integer, parameter :: kldp1 = kld + 1
  integer, parameter :: jl    = 12 
  integer, parameter :: jlp1  = jl + 1

  integer :: used1 = 0, dd
  integer:: inputxflag = 2, inputyflag = 2

!-------------------------------------------------------------------------------

  real, parameter :: acu       = 0.            ! --- Theoretical coefficient for wave-current interaction
!  real, parameter :: beta0     = 1.12         ! --- Coefficient for wind input.
  real, parameter :: beta0     = 1.0           ! --- Coefficient for wind input.
  real, parameter :: beta10    = beta0*0.25*1.25*0.001
  real, parameter :: rs        = 6367451.637   ! --- The global (Earth) radius.
  real, parameter :: pi        = 3.1415926     ! --- Pi.
  real, parameter :: pi2       = pi/2.0        ! --- Pi/2.
  real, parameter :: zpi       = 2.0*pi        ! --- 2*Pi.
  real, parameter :: pid180    = pi/180.       ! --- Pi/180.
  real, parameter :: g         = 9.81          ! --- Acceleration of gravity.
  real, parameter :: gc2       = 0.877**2 * g  ! --- Acceleration of gravity.
  real, parameter :: gg        = g * g         ! --- Acceleration of gravity.
  real, parameter :: tztz      = 1.099314      ! --- Coefficient for zero-crossing wave period.
  real, parameter :: d1        = 0.000132      ! --- Coefficient in wave-breaking dissipation formula
  real, parameter :: d2        = 2.61          ! --- Coefficient in wave-breaking dissipation formula
  real, parameter :: pwk       = 1.21          ! --- The constant for discretion of wave-number
  !real, parameter :: alog10pwk = alog10(pwk)   ! --- alog10(pwk)
  real, parameter :: alog10pwk = 0.08278537031645008150039994248605
  real, parameter :: wkmax     = 0.6894        ! --- Maximum of wave-umber amplitude
  real, parameter :: wkmin     = 0.0071        ! --- Minimum of wave-number amplitude
  real, parameter :: wfmax     = 0.413         ! --- Maximum of wave frequency
  real, parameter :: wfmin     = 0.042         ! --- Minimum of wave frequency
  real, parameter :: ads       = 1.0           ! --- Theoretical coefficient for wave-breaking dissipation
  real, parameter :: abo       = 1.0           ! --- Theoretical coefficient for bottom dissipation
  real, parameter :: p         = 0.025         ! --- Coefficient for growth spectrum limiter
  real, parameter :: cksp      = 14.0          ! --- Coefficient for estimation of wave-number range
  real, parameter :: cksa      = 4.5           ! --- Coefficient for estimation of wave-number range
  real, parameter :: small     = 0.000001      ! --- Small value.

  character(len=12), parameter ::  rstfile = 'wave_rest.nc'

!-------------------------------------------------------------------------------

  real, allocatable :: x(:)                    ! --- For longitude: rang from 0~360.
  real, allocatable :: rx(:)                   ! --- Real longitude.
  real, allocatable :: y(:)                    ! --- For latitude
  real, allocatable :: rslat(:)                ! --- For rs * cos(lat)
  real, allocatable :: d(:, :)                 ! --- For topography: the bottom depth (m) for topography
  real, allocatable :: nsp(:, :)               ! --- Computing mask for scalar variables: 
                                               !     0 for land, 1 for water and 2 for open boundary.
  real, allocatable :: deltx(:), delty(:)      ! --- For grid interval.
  real, allocatable :: dddx(:, :), dddy(:, :)  ! --- For d(d)/d(x) and d(d) / d(y)

  real, allocatable :: e(:, :, :, :)           ! --- Wave spectrum, before the influence of source function & propagation.
  real, allocatable :: ee(:, :, :, :)          ! --- Wave spectrum, after the influence of source function & propagation.
  real, allocatable :: ea(:, :, :, :)          ! --- Time averaged ee.
  !real, allocatable :: em(:, :, :, :)          ! --- Huang Chenghuan 20170319

  real, allocatable :: sein(:, :)              ! --- (kl, jl)
  real, allocatable :: seds(:, :)              ! --- (kl, jl)
  real, allocatable :: sebo(:, :)              ! --- (kl, jl)

  real, pointer :: pein(:, :)                  ! --- (ixl, iyl)
  real, pointer :: peds(:, :)                  ! --- (ixl, iyl)
  real, pointer :: pebo(:, :)                  ! --- (ixl, iyl)
  !real, allocatable :: grids(:, :)             ! --- (ixl, iyl) 2013/4/11 18:07:13

!-------------------------------------------------------------------------------

  real, pointer :: wx(:, :)      ! --- The wind along longitude
  real, pointer :: wy(:, :)      ! --- The wind along latitude
  real, allocatable :: w(:, :)       ! --- The wind speed.

!-------------------------------------------------------------------------------

  real, allocatable :: grolim(:)               ! --- Growth spectrum limiter
  real, allocatable :: thet(:)                 ! --- Directional discrete for 12 bands (300 resolution)
  real, allocatable :: wk(:)                   ! --- Wave-number discrete for 25 bands on a logarithmic scale
  real, allocatable :: wkh(:)                  ! --- coefficient for computing high-frenquency spectrum
  real, allocatable :: dwk(:)                  ! --- Discrete bands of wave-number
                                               
  real, allocatable :: wp(:, :, :)             ! --- weighting factor of K+ for decomposing wave-wave transfer
  real, allocatable :: wm(:, :, :)             ! --- weighting factor of K- for decomposing wave-wave transfer
  real, allocatable :: wks17(:)                ! --- 17/2 exponential of wave-number, that is K**(17/2)
                                               
  real, allocatable :: se(:, :)                ! --- Increment of spectrum (m4/s)
  real, allocatable :: dse(:, :)               ! --- Deferential of se with respect to spectrum (s-1)
  real, allocatable :: enh(:, :)               ! --- shallow water factor for nonlinear wave-wave transfer

!-------------------------------------------------------------------------------

  real, allocatable :: ae(:, :)                ! --- The zero-order moment of the spectrum
  real, allocatable :: awf(:, :)               ! --- The frequency first-order moment of the spectrum
  real, allocatable :: asi(:, :)               ! --- The frequency negative first-order moment of the spectrum
  real, allocatable :: awk(:, :)               ! --- The frequency negative second-order moment of the spectrum
  real, allocatable :: ark(:, :)               ! --- spectral mean of wave-number
  real, allocatable :: hb(:,:)                 ! --- The Hs limiter through wave breaking
  real, allocatable :: hbb(:, :)               ! --- The Hs limiter through wave breaking
                                               
  real, allocatable :: fconst0(:, :, :)        ! --- computing mask for threshold of wave-number 
  real, allocatable :: wf(:, :, :)             ! --- wave frequency
  real, allocatable :: ccg(:, :, :)            ! --- group velocity
  real, allocatable :: dwf(:, :, :)            ! --- discrete bands of wave frequency multiplied by half of directional discrete
  real, allocatable :: ssbos(:, :, :)          ! --- Huang Chenghuan, 20170423
                                              
  integer, allocatable :: ks0(:, :)            ! --- threshold of wave-number determined by wind speed or spectral mean of wave-number
  integer, allocatable :: kpmt0(:, :)          ! --- threshold of wave-number determined by wind speed 
  integer, allocatable :: kakt0(:, :)          ! --- threshold of wave-number determined by spectral mean of wave-number
                                               
  integer, allocatable :: jp1(:, :)            ! --- Discrete-interaction configurations of  nonlinear wave-wave transfer
  integer, allocatable :: jp2(:, :)            ! --- Discrete-interaction configurations of  nonlinear wave-wave transfer
  integer, allocatable :: jm1(:, :)            ! --- Discrete-interaction configurations of  nonlinear wave-wave transfer
  integer, allocatable :: jm2(:, :)            ! --- Discrete-interaction configurations of  nonlinear wave-wave transfer
                                               ! --- angular discrete-interaction configurations of nonlinear wave-wave transfer                                               

  integer, allocatable :: ikp(:)               ! --- wave-number discrete-interaction configurations of nonlinear wave-wave transfer 
  integer, allocatable :: ikp1(:)              ! --- wave-number discrete-interaction configurations of nonlinear wave-wave transfer 
  integer, allocatable :: ikm(:)               ! --- wave-number discrete-interaction configurations of nonlinear wave-wave transfer 
  integer, allocatable :: ikm1(:)              ! --- wave-number discrete-interaction configurations of nonlinear wave-wave transfer 

!-------------------------------------------------------------------------------

  real, pointer :: h1_3(:,:)     ! --- hs: significant wave height (m)
  real, pointer :: aet(:,:)      ! --- th: mean wave direction (Deg)
  real, pointer :: tpf(:,:)      ! --- tp: spectrum peak wave period (s)
  real, pointer :: ape(:,:)      ! --- tz: zero-crossing wave period (s)

!-------------------------------------------------------------------------------

! --- Input parameters. ( in the file of ctlparams)

  character(len=100) :: wind_path  ! --- Path for wind.
  character(len=100) :: data_path  ! --- Path for model setting.
  character(len=14 ) :: title      ! --- Symbal for model output.

!-------------------------------------------------------------------------------

  real    :: delttm       ! --- Length of integral time step, in minutes.
  real    :: lonref       ! --- The real longitude for lon=0 used in model.
  integer :: istime(6)    ! --- Integral start time
  integer :: ietime(6)    ! --- Integral end time
  integer :: cools_days   ! --- The time (days) for cool start.
  integer :: wndfreq      ! --- The frequency of wind data, in hours.
  integer :: wndtype      ! --- Wind type: 0 for same grid with model
                          !                1 for GFS wind data, no interp.
  integer :: outflag      ! --- The method of output. 
                          !     0, only wave output; 
                          !     1 wave & mix together; 
                          !     2 for wave & mix seperately.
                          !     3 Only output wave variables into file 
                          !       multi-records, everyday 1 file.
                          !     4 Only output wave variables into file 
                          !       multi-records, one run 1 file
  integer :: wiofreq      ! --- The output frequency for wave results (hour).
  integer :: ciofreq      ! --- The output frequency for current coef.s (hour).
  integer :: rstfreq      ! --- The output frequency for model restart (hour).

!-------------------------------------------------------------------------------

  integer :: glbflag      ! --- global flag, 0 for global model, 1 for regional.
  integer :: mixflag      ! --- Flag for mixing method, 0 for none, 
                          !                             1 for unlimited, 
                          !                             2 for limited depth.
                          !                             3 for both methods.

  integer :: flageng, flagmpi = -1
                            
!-------------------------------------------------------------------------------

! --- For variables about time.

  integer :: itime(6)
  double precision :: dtime, dtime0, dtimeend
  double precision :: nxttim                      !jiangxj 2012.6.1
  character(len=14) :: ctime, cistime, cietime

! --- ia for x-direction, ic for y-direction, it for time steps, itend for
!     number of all steps, iprint for print interval.

  integer :: iwiofreq, iciofreq, irstfreq
  integer :: it, itend, ia, ic, j, k, key, number

  real :: cong, al11, al21, al31, al12, al22, al13, al23
  real :: deltt, deltt5, deltth

!-------------------------------------------------------------------------------

  namelist/ctlparams/data_path, wind_path, title, cistime, cietime, cools_days, &
                     delttm, wndfreq, wndtype, outflag, wiofreq, ciofreq, rstfreq

  namelist/ctlparams1/title, istime, ietime, cools_days, delttm, lonref,       &
                      glbflag, wndfreq, wndtype, outflag, wiofreq, iwiofreq,   &
                      ciofreq, iciofreq, rstfreq, irstfreq, ixs, ixl, iys, iyl
                       
!-------------------------------------------------------------------------------

  contains
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: wamvar_mod_init

  subroutine wamvar_mod_init

!-------------------------------------------------------------------------------
  
  implicit none

!-------------------------------------------------------------------------------

  ixs=1
  ix1=ixs+1
  ix2=ixl-1

  iys=1
  iy1=iys+1
  iy2=iyl-1

! For saved indexs and intermediate values
! huangh223, 2017-02-25
  allocate(idxs(5, kl, jl, ixs:ixl, iys:iyl))
  allocate(tmp_values(6, kl, jl, ixs:ixl, iys:iyl))
  
!-------------------------------------------------------------------------------

  allocate(x(ixs:ixl))
  allocate(rx(ixs:ixl))
  allocate(y(iys:iyl))
  allocate(rslat(iys:iyl))
  allocate(d(ixs:ixl, iys:iyl))
  allocate(nsp(ixs:ixl, iys:iyl))

  allocate(deltx(ixs:ixl), delty(iys:iyl))
  allocate(dddx(ixs:ixl, iys:iyl), dddy(ixs:ixl, iys:iyl))

  allocate(e(kl, jl, ixs:ixl, iys:iyl))
  allocate(ee(kl, jl, ixs:ixl, iys:iyl))
  !allocate(em(kl,jl, ixs:ixl, iys:iyl))

  allocate(sein(kl, jl), seds(kl, jl), sebo(kl, jl))
  allocate(pein(ixl, iyl), peds(ixl, iyl), pebo(ixl, iyl)) 
  !2013/4/11 18:05:44, grids(ixl, iyl))

!-------------------------------------------------------------------------------

! --- For wind data.

  allocate(wx(ixs:ixl, iys:iyl))
  allocate(wy(ixs:ixl, iys:iyl))
  allocate(w(ixs:ixl, iys:iyl))

!-------------------------------------------------------------------------------

  allocate(grolim(kl)  )
  allocate(thet(jlp1)  )
  allocate(wk(kldp1)   )
!  allocate(wkh(kld)    ) ! yinxq
!  allocate(dwk(kld)    ) ! yinxq
  allocate(wkh(kldp1)    ) ! yinxq
  allocate(dwk(kldp1)    ) ! yinxq
  allocate(wp(kl, 2, 2)  )
  allocate(wm(kl, 2, 2)  )
  allocate(wks17(kl)   )
  allocate(se(klp1, jl) )
  allocate(dse(klp1, jl))
  allocate(enh(ixs:ixl, iys:iyl))

!-------------------------------------------------------------------------------

  allocate(ae(ixs:ixl, iys:iyl)        )
  allocate(awf(ixs:ixl, iys:iyl)       )
  allocate(asi(ixs:ixl, iys:iyl)       )
  allocate(awk(ixs:ixl, iys:iyl)       )
  allocate(ark(ixs:ixl, iys:iyl)       )
  allocate(hb(ixs:ixl, iys:iyl)        )
  allocate(hbb(ixs:ixl, iys:iyl)       )
  allocate(fconst0(kl,ixs:ixl, iys:iyl))
  allocate(wf(kldp1,ixs:ixl, iys:iyl) )
  allocate(ccg(kldp1,ixs:ixl, iys:iyl))
  allocate(dwf(kldp1,ixs:ixl, iys:iyl))
  allocate(ssbos(kldp1,ixs:ixl, iys:iyl))

  allocate(ks0(ixs:ixl, iys:iyl))
  allocate(kpmt0(ixs:ixl, iys:iyl))
  allocate(kakt0(ixs:ixl, iys:iyl))
  allocate(jp1(2,jl)     )
  allocate(jp2(2,jl)     )
  allocate(jm1(2,jl)     )
  allocate(jm2(2,jl)     )
  allocate(ikp(kl)       )
  allocate(ikp1(kl)      )
  allocate(ikm(kl)       )
  allocate(ikm1(kl)      )

!-------------------------------------------------------------------------------

! --- hs: significant wave height (m)
  allocate(h1_3(ixs:ixl, iys:iyl))
  allocate(aet(ixs:ixl, iys:iyl) )
  allocate(tpf(ixs:ixl, iys:iyl) )
  allocate(ape(ixs:ixl, iys:iyl) )

!-------------------------------------------------------------------------------

  e = small
  ee = small

  h1_3 = 0.0
  aet = 0.0
  tpf = 0.0
  ape = 0.0

  pein = 0.0
  pebo = 0.0
  peds = 0.0

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------
  
  end subroutine wamvar_mod_init

!-------------------------------------------------------------------------------

  end module wamvar_mod

!-------------------------------------------------------------------------------
!###############################################################################
