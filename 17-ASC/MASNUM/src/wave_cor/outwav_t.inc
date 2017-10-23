!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outwav_t

  subroutine outwav_t(filename, rec)

  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(in) :: rec

!-------------------------------------------------------------------------------

  integer :: ncid, jd
  integer*2, parameter :: ivland = nf_fill_int2

  logical :: ext
  integer :: outrecord, timerec, ittt
  double precision, allocatable :: timealready(:)
  
!-------------------------------------------------------------------------------

  if(mod_init /= 1)call wamfio_mod_init()
  
  jd = datenum(1950, 1, 1, 0, 0, 0)

!-------------------------------------------------------------------------------

  inquire(file=filename, exist=ext)
  
!  if(rec == 1)then
  if(.not.ext)then

    outrecord = 1
    
    call open_nc(ncid, filename, 'c')

    call dimension_define(ncid, 'lon', ixl, 'lon', nf_real)
    call dimension_define(ncid, 'lat', iyl, 'lat', nf_real)
    call dimension_define(ncid, 'time', 0, 'time', nf_double)    !jiangxj
    call set_attribute(ncid, 'units', 'degrees_north', 'lat')
    call set_attribute(ncid, 'units', 'degrees_east', 'lon')
    call set_attribute(ncid, 'modulo', '', 'lon')
    call set_attribute(ncid, 'units', 'Days since 1950-01-01 00:00:0.0.', 'time')
    call set_attribute(ncid, 'Start_time', ctime)

    call variable_define(ncid, 'windx', nf_int2, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'windy', nf_int2, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'hs',    nf_int2, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'tp',    nf_int2, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'tz',    nf_int2, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'th',    nf_int2, ['lon', 'lat', 'time'])

    call set_attribute(ncid, 'missing_value', ivland, 'windx')
    call set_attribute(ncid, 'missing_value', ivland, 'windy')
    call set_attribute(ncid, 'missing_value', ivland, 'hs')
    call set_attribute(ncid, 'missing_value', ivland, 'tp')
    call set_attribute(ncid, 'missing_value', ivland, 'tz')
    call set_attribute(ncid, 'missing_value', ivland, 'th')

    call set_attribute(ncid, 'scale_factor', 0.01, 'windx')
    call set_attribute(ncid, 'scale_factor', 0.01, 'windy')
    call set_attribute(ncid, 'scale_factor', 0.01, 'hs')
    call set_attribute(ncid, 'scale_factor', 0.01, 'tp')
    call set_attribute(ncid, 'scale_factor', 0.01, 'tz')
    call set_attribute(ncid, 'scale_factor', 0.1, 'th')

    call set_attribute(ncid, 'units', 'm/s', 'windx')
    call set_attribute(ncid, 'units', 'm/s', 'windy')
    call set_attribute(ncid, 'units', 'm'  , 'hs')
    call set_attribute(ncid, 'units', 's'  , 'tp')
    call set_attribute(ncid, 'units', 's'  , 'tz')
    call set_attribute(ncid, 'units', 'deg', 'th')

    call set_attribute(ncid, 'longname', 'Zonal Wind Velocity '     , 'windx')
    call set_attribute(ncid, 'longname', 'Meridional Wind Velocity ', 'windy')
    call set_attribute(ncid, 'longname', 'Significant wave height'  , 'hs')
    call set_attribute(ncid, 'longname', 'Mean wave direction'      , 'th')
    call set_attribute(ncid, 'longname', 'Spectrum peak wave period', 'tp')
    call set_attribute(ncid, 'longname', 'Zero-crossing wave period', 'tz')

    call end_define(ncid)
  
    call writenc(ncid, 'lon', rx)
    call writenc(ncid, 'lat', y)

  else

    call open_nc(ncid, filename, 'w')
    
    timerec = get_dimension_len(ncid, 'time')
    allocate(timealready(timerec))
    call readnc(ncid, 'time', timealready)
    
    outrecord = 1
    do ittt = 1, timerec !-1
    	if(dtime-jd > timealready(ittt))then ! .and. dtime-jd <= timealready(ittt+1))then
    	  outrecord = ittt + 1
    	endif
    enddo
    
    deallocate(timealready)
    
    write(6, *)'outrecord = ', outrecord

  endif

!-------------------------------------------------------------------------------

  call writenc(ncid, 'time', dtime-jd, outrecord)
  call setland_iv2(wx, 0.01)
  call writenc(ncid, 'windx', iv2, outrecord)
  call setland_iv2(wy, 0.01)
  call writenc(ncid, 'windy', iv2, outrecord)
  call setland_iv2(h1_3, 0.01)
  call writenc(ncid, 'hs', iv2, outrecord)
  call setland_iv2(tpf, 0.01)
  call writenc(ncid, 'tp', iv2, outrecord)
  call setland_iv2(ape, 0.01)
  call writenc(ncid, 'tz', iv2, outrecord)
  call setland_iv2(aet, 0.1)
  call writenc(ncid, 'th', iv2, outrecord)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outwav_t

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setland_iv2

  subroutine setland_iv2(var, scal)

  implicit none

  real, intent(in) :: var(:, :), scal
  integer :: i, j

  do j = iys, iyl
  do i = ixs, ixl
!    if(nsp(i, j) == 0)then
    if(nsp(i, j) /= 1)then
      iv2(i, j) = nf_fill_int2
    else
      iv2(i, j) = var(i, j) / scal
    endif
  enddo
  enddo

  return

  end subroutine setland_iv2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setland_v2

  subroutine setland_v2(var)

  implicit none

  real, intent(in) :: var(:, :)
  integer :: i, j

  do j = iys, iyl
  do i = ixs, ixl
    if(nsp(i, j) == 0)then
      v2(i, j) = nf_fill_real
    else
      v2(i, j) = var(i, j)
    endif
  enddo
  enddo

  return

  end subroutine setland_v2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setland_v3

  subroutine setland_v3(var)

  implicit none

  real, intent(in) :: var(:, :, :)

  integer :: i, j


  do j = iys, iyl
  do i = ixs, ixl
    if(nsp(i, j) == 0)then
      v3(i, j, :) = nf_fill_real
    else
      v3(i, j, :) = var(i, j, :)
    endif
  enddo
  enddo

  return

  end subroutine setland_v3

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
