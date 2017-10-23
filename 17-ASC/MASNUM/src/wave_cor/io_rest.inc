!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outrst

  subroutine outrst(filename)

  implicit none

  character(len=*), intent(in) :: filename

  integer :: ncid

!-------------------------------------------------------------------------------

  if(mod_init /= 1)call wamfio_mod_init()
! Modified by Zhenya.Song, 2016/04/03  
!  if(mod(it-number, irstfreq) /= 0)return
  if (irstfreq == 0) then
    return
  else
   if(mod(it-number, irstfreq) /= 0) return
  endif
! End Modify

!-------------------------------------------------------------------------------

  call open_nc(ncid, filename, 'c')

  call dimension_define(ncid, 'kk', kl, 'kk', nf_real)
  call dimension_define(ncid, 'jj', jl, 'jj', nf_real)
  call dimension_define(ncid, 'lon', ixl, 'lon', nf_real)
  call dimension_define(ncid, 'lat', iyl, 'lat', nf_real)

  call variable_define(ncid, 'ee', nf_real, ['kk', 'jj', 'lon', 'lat'])

  call set_attribute(ncid, 'units', 'degrees_north', 'lat')
  call set_attribute(ncid, 'units', 'degrees_east', 'lon')
  call set_attribute(ncid, 'modulo', '', 'lon')
  call set_attribute(ncid, 'ctime', ctime)

  call end_define(ncid)

!  call writenc(ncid, 'kk', x) ! wk(1:kl)
!  call writenc(ncid, 'jj', x) ! 1:13
  call writenc(ncid, 'lon', rx)
  call writenc(ncid, 'lat', y)
  call writenc(ncid, 'ee', ee)

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outrst

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: inprst

  subroutine inprst(filename, key)

  implicit none

  integer, intent(out) :: key
  character(len=*), intent(in) :: filename

  integer :: ncid
  logical :: ext

!-------------------------------------------------------------------------------

  if(mod_init /= 1)call wamfio_mod_init()

  inquire(file=filename, exist=ext)
  if(ext)then
    key = 1
  else
    key = 0
    return
  endif

!-------------------------------------------------------------------------------

  call open_nc(ncid, filename, 'r')
  call readnc(ncid, 'ee', ee)
  call get_attribute(ncid, 'ctime', ctime)
  call close_nc(ncid)

  write(*, *)'Restart time is :', ctime
  read(ctime, '(i4.4,5i2.2)')istime

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine inprst

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
