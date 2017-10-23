!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix

  subroutine outmix(filename)

  implicit none

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

  call variable_define(ncid, 'tau11', nf_real, ['lon', 'lat', 'dep'])
  call variable_define(ncid, 'tau12', nf_real, ['lon', 'lat', 'dep'])
  call variable_define(ncid, 'tau22', nf_real, ['lon', 'lat', 'dep'])
  call variable_define(ncid, 'tau33', nf_real, ['lon', 'lat', 'dep'])

  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau11')
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau12')
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau22')
  call set_attribute(ncid, 'missing_value', nf_fill_real, 'tau33')

  if(mixflag == 1 .or. mixflag == 3)then

    call variable_define(ncid, 'bv'   , nf_real, ['lon', 'lat', 'dep'])
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'bv')

  endif

  if(mixflag == 2 .or. mixflag == 3)then

    call variable_define(ncid, 'bvl'   , nf_real, ['lon', 'lat', 'dep'])
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'bvl')

    call variable_define(ncid, 'bh1'   , nf_real, ['lon', 'lat', 'dep'])
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'bh1')

    call variable_define(ncid, 'bh2'   , nf_real, ['lon', 'lat', 'dep'])
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'bh2')

    call variable_define(ncid, 'bv_wtv'   , nf_real, ['lon', 'lat', 'dep'])
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'bv_wtv')

    call variable_define(ncid, 'bv_wtd'   , nf_real, ['lon', 'lat', 'dep'])
    call set_attribute(ncid, 'missing_value', nf_fill_real, 'bv_wtd')

  endif

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

  if(mixflag == 1 .or. mixflag == 3)then

    call setland_v3(bv)
    call writenc(ncid, 'bv'   , v3)

  endif

  if(mixflag == 2 .or. mixflag == 3)then

    call setland_v3(bvl)
    call writenc(ncid, 'bvl'   , v3)
    call setland_v3(bh1)
    call writenc(ncid, 'bh1'   , v3)
    call setland_v3(bh2)
    call writenc(ncid, 'bh2'   , v3)

    call setland_v3(bv_wtv)
    call writenc(ncid, 'bv_wtv'   , v3)
    call setland_v3(bv_wtd)
    call writenc(ncid, 'bv_wtd'   , v3)

  endif

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outmix

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
