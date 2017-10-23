!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix_tau

  subroutine outmix_tau(filename)

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

  call variable_define(ncid, 'tau11', nf_real, ['lon', 'lat', 'dep'])
  call variable_define(ncid, 'tau12', nf_real, ['lon', 'lat', 'dep'])
  call variable_define(ncid, 'tau22', nf_real, ['lon', 'lat', 'dep'])
  call variable_define(ncid, 'tau33', nf_real, ['lon', 'lat', 'dep'])

  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau11')
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau12')
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau22')
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau33')

  call end_define(ncid)

!-------------------------------------------------------------------------------

  call writenc(ncid, 'lon', rx)
  call writenc(ncid, 'lat', y)
  call writenc(ncid, 'dep', zyyz)

  call setland_v3(taubb11)
  call writenc(ncid, 'tau11', v3)
  call setland_v3(taubb12)
  call writenc(ncid, 'tau12', v3)
  call setland_v3(taubb22)
  call writenc(ncid, 'tau22', v3)
  call setland_v3(taubb33)
  call writenc(ncid, 'tau33', v3)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outmix_tau

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
