!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: wamfio_mod_init

  subroutine wamfio_mod_init(lon, lat, mask, dep, top)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  real, optional, intent(in) :: lon(:), lat(:), mask(:, :), dep(:), top(:, :)

  integer :: ncid

!-------------------------------------------------------------------------------

! --- Input the parameters used in this model.

  title      = 'glb500'
  istime     = [2005, 1, 1, 0, 0, 0]
  ietime     = [2005, 1, 30, 0, 0, 0]
  cools_days = 5
  glbflag    = 1
  mixflag    = 0
  
  open(11, file='ctlparams')
  read(11, nml=ctlparams)
  close(11)

!-------------------------------------------------------------------------------

! --- Common parameters for wind data & module init.

  wflag    = 0
  dwtime1  = -1.d0
  dwtime2  = -1.d0
  mod_init = 1
  recd     = 0

  istime = 0; ietime = 0; mixflag = 0
  read(cistime, '(i4,i2,i2)')istime(1:3)
  read(cietime, '(i4,i2,i2)')ietime(1:3)

! --- mixflag = 1: tau, bv in one file.
! --- mixflag = 2: tau, bvl in one file.
! --- mixflag = 3: tau, bv, bvl in one file.
! --- mixflag = 4: bv in one file.
! --- mixflag = 5: tau in one file; bv in one file; bvl in one file.

  if(ciofreq /= 0)mixflag = 1

!-------------------------------------------------------------------------------

  if(present(lon))then
    ixl = size(lon)
    iyl = size(lat)
    kb = size(dep)
  else
    call open_nc(ncid, trim(data_path), 'r')
    ixl = get_dimension_len(ncid, 'lon')
    iyl = get_dimension_len(ncid, 'lat')
    kb = get_dimension_len(ncid, 'zyyz')
    call close_nc(ncid)
  endif

!-------------------------------------------------------------------------------

  allocate(windx1(ixl, iyl), windx2(ixl, iyl))
  allocate(windy1(ixl, iyl), windy2(ixl, iyl))
  allocate(v2(ixl, iyl), v3(ixl, iyl, kb), iv2(ixl, iyl))

!-------------------------------------------------------------------------------

! --- Initialize for variables.

  call wamvar_mod_init

  call init_wamcpl(mixflag)

  if(present(lon))then
    
    rx = lon; y = lat; zyyz = dep; d = top; nsp = mask

  else

    call open_nc(ncid, trim(data_path), 'r')
    call readnc(ncid, 'lon', rx)
    call readnc(ncid, 'lat', y)
    call readnc(ncid, 'zyyz', zyyz)
    call readnc(ncid, 'depyyz', d)
    call readnc(ncid, 'nspyyz', nsp)
    call close_nc(ncid)

  endif

  call settopog
  call setwave
  call nlweight
  call set_uv0

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine wamfio_mod_init

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
