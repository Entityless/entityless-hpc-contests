!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_ice

  subroutine set_ice

!-------------------------------------------------------------------------------

  use netcdf_mod
  use time_mod
  
  implicit none

!-------------------------------------------------------------------------------

  integer :: im, jm, ncid, i, j
  real, allocatable :: xice(:), yice(:), xweight(:), yweight(:)
  integer, allocatable :: ice(:, :), ilist(:, :), jlist(:, :)
  logical :: ext
  real :: a, b

  integer :: iceitime(6)
  
!-------------------------------------------------------------------------------

!  write(6, *)ctime, 'ice'
!  if(sum(itime(4:5)) /= 0 .or. itime(3) /= 1)return
!  if(last_time_ice >= dtime)return  !yinxq 2011/10/31 8:51:17
  if(last_time_ice >= dtime)return

!-------------------------------------------------------------------------------

  inquire(file='ice_clim_mask.nc', exist=ext)
  if(.not.ext)return

!-------------------------------------------------------------------------------

! --- Initial for seting ice.

  if(init_ice == 0)then
    allocate(icensp(ixs:ixl, iys:iyl), noicensp(ixs:ixl, iys:iyl))
    noicensp = nsp
    init_ice = 1
  endif

!-------------------------------------------------------------------------------

! --- Input ice mask, 0 for land/ice, 1 for water.

  call open_nc(ncid, 'ice_clim_mask.nc', 'r')
  im = get_dimension_len(ncid, 'x')
  jm = get_dimension_len(ncid, 'y')
  allocate(xice(im), yice(jm), ice(im, jm),     &
          ilist(ixs:ixl, 2), jlist(iys:iyl, 2), &
          xweight(ixs:ixl), yweight(iys:iyl)    )
  call readnc(ncid, 'x', xice)
  call readnc(ncid, 'y', yice)
  call readnc(ncid, 'ice_clim_mask', ice, itime(2))
  call close_nc(ncid)

  do ia = ixs, ixl
    call find_i1i2a(xice, rx(ia), im, ilist(ia, 1), ilist(ia, 2), xweight(ia))
  enddo
  do ic = iys, iyl
    call find_i1i2a(yice, y(ic), jm, jlist(ic, 1), jlist(ic, 2), yweight(ic))
  enddo

  do ia = ixs, ixl
  do ic = iys, iyl
    a = xweight(ia); b = yweight(ic)
    icensp(ia, ic) = a     * b     * ice(ilist(ia, 1), jlist(ic, 1)) &
                   + a     * (1-b) * ice(ilist(ia, 1), jlist(ic, 2)) &
                   + (1-a) * (1-b) * ice(ilist(ia, 2), jlist(ic, 2)) &
                   + (1-a) * b     * ice(ilist(ia, 2), jlist(ic, 1))
    if(icensp(ia, ic) < 1)icensp(ia, ic) = 0
  enddo
  enddo

  !write(6, *)'set_ice: ', minval(icensp), maxval(icensp)

  deallocate(xice, yice, ice, ilist, jlist, xweight, yweight)

!-------------------------------------------------------------------------------

! --- Update mask in nsp.

  nsp = icensp * noicensp
!  last_time_ice = dtime

  iceitime = datevec(dtime)
  iceitime(2:3) = [iceitime(2) + 1, 1]
  if(iceitime(2) > 12)iceitime(1:2) = [iceitime(1)+1, 1]
  last_time_ice = datenum(iceitime)

! --- Reset time steps.

  call reset_time_step

  !write(6, *)'set_ice is okay.', ctime

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine set_ice

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_uv

  subroutine set_uv

  implicit none

!  call reset_time_step

  return

  end subroutine set_uv

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_uv0

  subroutine set_uv0

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  ux = 0.0
  uy = 0.0
  uxx = 0.0
  uyx = 0.0
  uxy = 0.0
  uyy = 0.0

  call reset_time_step
  !write(6, *)'set_uv0 is okay.'

!-------------------------------------------------------------------------------

  end subroutine set_uv0

!-------------------------------------------------------------------------------

  subroutine reset_time_step

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  real, allocatable :: tim(:, :)
  real :: maxtimestep, prop_time_step, iii

!-------------------------------------------------------------------------------

  allocate(tim(ixs:ixl, iys:iyl))

  tim = -10.
  do ia = ixs, ixl
  do ic = iys, iyl
    if(nsp(ia,ic) /= 1)cycle
    tim(ia, ic) = min(deltx(ia) * rslat(ic) / (ccg(1, ia, ic) + ux(ia, ic)), &
                      delty(ic) * rs        / (ccg(1, ia, ic) + uy(ia, ic))  )
  enddo
  enddo

  tim = tim * pi / 1.08e4 ! 10800 = 180(deg) * 60(m)

  maxtimestep = minval(tim, mask=tim>0.)
  
  deallocate(tim)

!-------------------------------------------------------------------------------
! ----------delttm is defined by cltparams !jiangxj
!  prop_time_step = 60.
!  write(6, *)'maxtimestep, prop_time_step = ', maxtimestep, prop_time_step
!   do while (maxtimestep < prop_time_step)
!    iii = int(60. / prop_time_step)
!     iii = nint(60. / prop_time_step)       !jiangxj
!     iii = iii + 1
!     prop_time_step = 60. / iii
!     write(6, *)'maxtimestep, prop_time_step = ', maxtimestep, prop_time_step
!   enddo
  
!  prop_time_step = 1.2554928e-2   !shenhj
  
!shenhj  prop_time_step = 0.002
! jiangxj temporary changed for wang's project
!  prop_time_step = 0.015

!  delttm = prop_time_step

  iii = max(1, nint(60. / maxtimestep))
  if( 60. / iii > maxtimestep)then
    iii = iii + 1
  endif
  
! Modified by Zhenya.Song, 2015/12/27
!  delttm = 60. / iii
  delttm = min(delttm,60. / iii)
! End modify
  !write(6, *)'delttm, maxtimestep = ', delttm, maxtimestep
    
!-------------------------------------------------------------------------------

  deltt    = delttm  * 60.  ! yinxq
  deltt5   = delttm  * 30.  ! yinxq
  iwiofreq = wiofreq * 60. / delttm
  iciofreq = ciofreq * 60. / delttm
  irstfreq = rstfreq * 60. / delttm
  number = (1-key) * cools_days * 1440. / delttm

!-------------------------------------------------------------------------------

!  write(6, *)'***************************************************************'
!  write(6, *)'data_path = ', trim(data_path)
!  write(6, *)'wind_path = ', trim(wind_path)
!  write(6, *)'---------------------------------------------------------------'
!  write(6, nml=ctlparams1)
!  write(6, *)'---------------------------------------------------------------'

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine reset_time_step

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: find_i1i2a

  subroutine find_i1i2a(x, xo, im, i1, i2, a)

  implicit none

  integer, intent(in) :: im
  real, intent(in) :: x(im), xo
  integer, intent(out) :: i1, i2
  real, intent(out) :: a

  integer :: i, flag

!-------------------------------------------------------------------------------

  if(x(1) >= x(2))then ! --- X is sorted ascending

    if(xo >= x(1))then
      i1 = 1; i2 = 1
    elseif(xo <= x(im))then
      i1 = im; i2 = im
    else
      do i = 1, im-1
        if(x(i) > xo .and. x(i+1) <= xo)exit
      enddo
      i1 = i + 1; i2 = i
    endif

!-------------------------------------------------------------------------------

  elseif(x(1) < x(2))then ! --- X is sorted descending

    if(xo <= x(1))then
      i1 = 1; i2 = 1
    elseif(xo >= x(im))then
      i1 = im; i2 = im
    else
      do i = 1, im-1
        if(x(i) <= xo .and. x(i+1) > xo)exit
      enddo
      i1 = i; i2 = i + 1
    endif

  endif

  if(i1 /= i2)then
    a = (x(i2) - xo) / (x(i2) - x(i1))
  else
    a = 0.5
  endif

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine find_i1i2a

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
