!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outeng_t

  subroutine outeng_t(filename, rec)

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

    call variable_define(ncid, 'pein',    nf_real, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'pebo',    nf_real, ['lon', 'lat', 'time'])
    call variable_define(ncid, 'peds',    nf_real, ['lon', 'lat', 'time'])

    call set_attribute(ncid, 'missing_value', nf_fill_real, 'pein')
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'pebo')
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'peds')

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

  endif

!-------------------------------------------------------------------------------

  call writenc(ncid, 'time', dtime-jd, outrecord)

  call setland_v2(pein)
  call writenc(ncid, 'pein', v2, outrecord)

  call setland_v2(pebo)
  call writenc(ncid, 'pebo', v2, outrecord)

  call setland_v2(peds)
  call writenc(ncid, 'peds', v2, outrecord)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outeng_t

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
