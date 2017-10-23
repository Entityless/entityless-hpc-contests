!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckShenhj: outmix_wit

  subroutine outmix_wit(filename)

!-------------------------------------------------------------------------------
  
  implicit none

!-------------------------------------------------------------------------------
  
  character(len=*), intent(in) :: filename
  
!-------------------------------------------------------------------------------

  integer :: ncid
  
!-------------------------------------------------------------------------------

	if(mod_init /= 1)call wamfio_mod_init
	
!-------------------------------------------------------------------------------

  call open_nc(ncid, filename, 'c')

  call dimension_define(ncid, 'lon', ixl, 'lon', nf_real)
  call dimension_define(ncid, 'lat', iyl, 'lat', nf_real)
  call dimension_define(ncid, 'dep', kb, 'dep', nf_real)
  call set_attribute(ncid, 'units', 'degrees_north', 'lat')
  call set_attribute(ncid, 'units', 'degrees_east', 'lon')
  call set_attribute(ncid, 'modulo', '', 'lon')
  call set_attribute(ncid, 'ctime', ctime)  
 
  call variable_define(ncid, 'bv_wtv'   , nf_real, ['lon', 'lat', 'dep'])
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'bv_wtv')

  call variable_define(ncid, 'bv_wtd'   , nf_real, ['lon', 'lat', 'dep'])
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'bv_wtd')

  call end_define(ncid)

!-------------------------------------------------------------------------------

  call writenc(ncid, 'lon', rx)
  call writenc(ncid, 'lat', y)
  call writenc(ncid, 'dep', zyyz)
 
  call setland_v3(bv_wtv)
  call writenc(ncid, 'bv_wtv'   , v3)
  
  call setland_v3(bv_wtd)
  call writenc(ncid, 'bv_wtd'   , v3)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outmix_wit

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
