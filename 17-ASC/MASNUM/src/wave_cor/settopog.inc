!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: settopog

  subroutine settopog

!-------------------------------------------------------------------------------

  use netcdf_mod

  implicit none

!-------------------------------------------------------------------------------

  integer :: i, j, ncid

!-------------------------------------------------------------------------------

!  if(mod_init /= 1)call wamfio_mod_init

! --- Read in x, y, rslat, d, nsp

!  call open_nc(ncid, trim(data_path)//'wamyyz.nc', 'r')
!  call readnc(ncid, 'lon', rx)
!  call readnc(ncid, 'lat', y)
!  call readnc(ncid, 'zyyz', zyyz)
!  call readnc(ncid, 'depyyz', d)
!  call readnc(ncid, 'nspyyz', nsp)
!  call close_nc(ncid)

!-------------------------------------------------------------------------------

!  x = rx; lonref = 0.0
!!  if(minval(rx) < 0)then ! ympi
!  if(minval(rx) < 0 .or. maxval(rx) > 360)then
!    lonref = minval(x)
!    x = x - lonref
!  endif

  x = rx; lonref = 0.0
  do ia = ixs, ixl
  	if(ia < ixl/2)then
  	  if(x(ia) > x(ia+3))x(ia) = x(ia) - 360.
  	else
  		if(x(ia) < x(ia-3))x(ia) = x(ia) + 360.
  	endif
  enddo

  if(minval(x) < 0 .or. maxval(x) > 360)then
    lonref = minval(x)
    x = x - lonref
  endif

!-------------------------------------------------------------------------------

  do ic = iys, iyl
    rslat(ic) = rs * cosd(y(ic))
  enddo

!-------------------------------------------------------------------------------

  if((x(ixs) + 360. -x(ixl)) < (x(ixs+1) - x(ixs)))then
    glbflag = 0 ! --- For global model.
  else
    glbflag = 1 ! --- For regional model.
  endif
  if(flagmpi /= -1)glbflag = 1  ! ympi

!-------------------------------------------------------------------------------

!  dx = x(2) - x(1)
!  dy = y(2) - y(1)

  do ia = ixs, ixl
    if(ia == ixs)then
      deltx(ia) = x(ia+1) - x(ia)
    elseif(ia == ixl)then
      deltx(ia) = x(ia) - x(ia-1)
    else
      deltx(ia) = (x(ia+1) - x(ia-1)) * 0.5
    endif
  enddo

  do ic = iys, iyl
    if(ic == iys)then
      delty(ic) = y(ic+1) - y(ic)
    elseif(ic == iyl)then
      delty(ic) = y(ic) - y(ic-1)
    else
      delty(ic) = (y(ic+1) - y(ic-1)) * 0.5
    endif
  enddo

  do ia = ixs, ixl
  do ic = iys, iyl

    if(glbflag == 0)then ! --- for global mdel. 2010.04.09
      if(ia == ixs)then
        dddx(ia, ic)=(d(ia+1,ic)-d(ixl-1,ic))/(2.*deltx(ia))
      elseif(ia == ixl)then
        dddx(ia, ic)=(d(ixs+1,ic)-d(ia-1,ic))/(2.*deltx(ia))
      else
        dddx(ia, ic)=(d(ia+1,ic)-d(ia-1,ic))/(2.*deltx(ia))
      endif
    else
      if(ia == ixs)then
        dddx(ia, ic)=(d(ia+1,ic)-d(ia,ic))/(deltx(ia))
      elseif(ia == ixl)then
        dddx(ia, ic)=(d(ia,ic)-d(ia-1,ic))/(deltx(ia))
      else
        dddx(ia, ic)=(d(ia+1,ic)-d(ia-1,ic))/(2.*deltx(ia))
      endif
    endif

    if(ic == iys)then
      dddy(ia, ic)=(d(ia,ic+1)-d(ia,ic))/(delty(ic))
    elseif(ic == iyl)then
      dddy(ia, ic)=(d(ia,ic)-d(ia,ic-1))/(delty(ic))
    else
      dddy(ia, ic)=(d(ia,ic+1)-d(ia,ic-1))/(2. * delty(ic))
    endif

  enddo
  enddo

!-------------------------------------------------------------------------------

! --- Checking the setting of longitude.

  if(flagmpi == -1)then ! ympi
  
	  if(maxval(x) > 360)then
	    write(6, *)'***STOP*** Please check the input longitude:'
	    write(6, *)'It should between -360 & 360, and max_lon - min_lon <= 360'
	    stop
	  endif

	  if(glbflag == 0 .and. rx(ixl) - rx(ixs) /= 360)then
	    write(6, *)'***STOP*** Please check the input longitude:'
	    write(6, *)'For global model, max_lon - min_lon = 360'
	    stop
	  endif
	  
  endif

  do ia = ixs, ixl-1
!    if(rx(ia) >= rx(ia+1))then
    if(x(ia) >= x(ia+1))then
      write(6, *)'***STOP*** Please check the input longitude:'
      write(6, *)'It should be increasing along index.'
      write(6, *)'longitude = ', x
      write(6, *)'Wrong lon: ', x(ia), x(ia+1), ia
!      write(6, *)'longitude = ', rx
!      write(6, *)'Wrong lon: ', rx(ia), rx(ia+1), ia
      stop
    endif
  enddo

!-------------------------------------------------------------------------------

! --- Checking the setting of latitude.

  do ic = iys, iyl-1
    if(y(ic) >= y(ic+1))then
      write(6, *)'***STOP*** Please check the input latitude:'
      write(6, *)'It should be increasing along index.'
      write(6, *)'latitude = ', y
      write(6, *)'Wrong lat: ', y(ic), y(ic+1), ic
      stop
    endif
  enddo

!-------------------------------------------------------------------------------

  return

  do ia = ixs, ixl
    do ic = iys, iyl
      if(abs(y(ic)) > 65)then
        nsp(ia, ic) = 0
      endif
    enddo
  enddo

  write(6, *)'cut 65'
  
!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine settopog

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
