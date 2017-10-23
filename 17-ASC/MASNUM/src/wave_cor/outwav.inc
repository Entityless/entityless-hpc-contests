!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outwav

  subroutine outwav(filename)

  implicit none

  character(len=*), intent(in) :: filename

!-------------------------------------------------------------------------------

  integer :: ncid

!-------------------------------------------------------------------------------

  if(mod_init /= 1)call wamfio_mod_init()

!-------------------------------------------------------------------------------

  call open_nc(ncid, filename, 'c')

  call dimension_define(ncid, 'lon', ixl, 'lon', nf_real)
  call dimension_define(ncid, 'lat', iyl, 'lat', nf_real)
  call set_attribute(ncid, 'units', 'degrees_north', 'lat')
  call set_attribute(ncid, 'units', 'degrees_east', 'lon')
  call set_attribute(ncid, 'modulo', '', 'lon')
  call set_attribute(ncid, 'ctime', ctime)

  call variable_define(ncid, 'windx', nf_real, ['lon', 'lat'])
  call variable_define(ncid, 'windy', nf_real, ['lon', 'lat'])

  call variable_define(ncid, 'hs', nf_real, ['lon', 'lat'])
  call variable_define(ncid, 'tp', nf_real, ['lon', 'lat'])
  call variable_define(ncid, 'tz', nf_real, ['lon', 'lat'])
  call variable_define(ncid, 'th', nf_real, ['lon', 'lat'])

  call set_attribute(ncid, 'missing_value', vland, 'windx')
  call set_attribute(ncid, 'missing_value', vland, 'windy')
  call set_attribute(ncid, 'missing_value', vland, 'hs')
  call set_attribute(ncid, 'missing_value', vland, 'tp')
  call set_attribute(ncid, 'missing_value', vland, 'tz')
  call set_attribute(ncid, 'missing_value', vland, 'th')

  call end_define(ncid)

!-------------------------------------------------------------------------------

  call writenc(ncid, 'lon', rx)
  call writenc(ncid, 'lat', y)

  call setland_v2(wx)
  call writenc(ncid, 'windx', v2)
  call setland_v2(wy)
  call writenc(ncid, 'windy', v2)

  call setland_v2(h1_3)
  call writenc(ncid, 'hs', v2)
  call setland_v2(tpf)
  call writenc(ncid, 'tp', v2)
  call setland_v2(ape)
  call writenc(ncid, 'tz', v2)
  call setland_v2(aet)
  call writenc(ncid, 'th', v2)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outwav

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------

