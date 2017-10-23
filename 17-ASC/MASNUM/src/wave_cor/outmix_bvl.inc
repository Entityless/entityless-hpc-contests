!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix_bvl

  subroutine outmix_bvl(filename)

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
  call dimension_define(ncid, 'dep', kb, 'dep', nf_real)
  call set_attribute(ncid, 'units', 'degrees_north', 'lat')
  call set_attribute(ncid, 'units', 'degrees_east', 'lon')
  call set_attribute(ncid, 'modulo', '', 'lon')
  call set_attribute(ncid, 'ctime', ctime)


  call variable_define(ncid, 'bvl'   , nf_real, ['lon', 'lat', 'dep'])
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'bvl')

  call variable_define(ncid, 'bh1'   , nf_real, ['lon', 'lat', 'dep'])
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'bh1')

  call variable_define(ncid, 'bh2'   , nf_real, ['lon', 'lat', 'dep'])
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'bh2')

  call end_define(ncid)

!-------------------------------------------------------------------------------

  call writenc(ncid, 'lon', rx)
  call writenc(ncid, 'lat', y)
  call writenc(ncid, 'dep', zyyz)

  call setland_v3(bvl)
  call writenc(ncid, 'bvl'   , v3)
  call setland_v3(bh1)
  call writenc(ncid, 'bh1'   , v3)
  call setland_v3(bh2)
  call writenc(ncid, 'bh2'   , v3)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outmix_bvl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
