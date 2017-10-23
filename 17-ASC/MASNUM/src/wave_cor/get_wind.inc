!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_wind

  subroutine my_getwind(my_dtime,my_itime)

  implicit none
  double precision, intent(in) :: my_dtime
  integer, intent(inout)  :: my_itime(6)
  integer :: ii, i, j
  real :: a, b

!-------------------------------------------------------------------------------

  if(mod_init /= 1)call wamfio_mod_init()

!-------------------------------------------------------------------------------

  if(dwtime1 == -1.d0 .and. dwtime2 == -1.d0)then

    !write(6, *)'WND: ', datestr(my_itime)
    ii = wndfreq * int(my_itime(4)/wndfreq)
    dwtime2 = datenum([my_itime(1), my_itime(2), my_itime(3), ii, 0, 0])

    ! --- read wind2
    if(wndtype == 0)then
      call read_wind(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 1)then
      call read_wind_gfs(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 2)then
      call read_wind_qbln_interp(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 3)then
      call read_wind_ncep_interp(dwtime2, windx2, windy2, wflag)
    endif

    if(wflag /= 0)then
      write(6, *)'WND: ', 'There is no wind data input.'
      stop
    endif

  endif

!-------------------------------------------------------------------------------

  if(my_dtime >= dwtime2 .and. wflag == 0)then

    !write(6, *)'WND: ', datestr(my_itime)
    dwtime1 = dwtime2
    windx1 = windx2
    windy1 = windy2

    ii = wndfreq * int(my_itime(4)/wndfreq) + wndfreq
    dwtime2 = datenum([my_itime(1), my_itime(2), my_itime(3), ii, 0, 0])

    ! --- read wind2
    if(wndtype == 0)then
      call read_wind(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 1)then
      call read_wind_gfs(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 2)then
      call read_wind_qbln_interp(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 3)then
      call read_wind_ncep_interp(dwtime2, windx2, windy2, wflag)
    endif

    if(wflag /= 0)then
      dwtime2 = dwtime1
      windx2 = windx1
      windy2 = windy1
    endif

  endif

!-------------------------------------------------------------------------------

  if(wflag /= 0)then
    wx = windx1
    wy = windy1
    w = sqrt(wx*wx + wy*wy)
    return
  endif

  b = (my_dtime - dwtime1)/(dwtime2 - dwtime1)
  a = 1 - b

  wx = windx1 * a + windx2 * b
  wy = windy1 * a + windy2 * b
  w = sqrt(abs(wx*wx + wy*wy))

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine my_getwind

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_wind

  subroutine get_wind

  implicit none

  integer :: ii, i, j
  real :: a, b

!-------------------------------------------------------------------------------

  if(mod_init /= 1)call wamfio_mod_init()

!-------------------------------------------------------------------------------

  if(dwtime1 == -1.d0 .and. dwtime2 == -1.d0)then

    !write(6, *)'WND: ', datestr(itime)
    ii = wndfreq * int(itime(4)/wndfreq)
    dwtime2 = datenum([itime(1), itime(2), itime(3), ii, 0, 0])

    ! --- read wind2
    if(wndtype == 0)then
      call read_wind(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 1)then
      call read_wind_gfs(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 2)then
      call read_wind_qbln_interp(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 3)then
      call read_wind_ncep_interp(dwtime2, windx2, windy2, wflag)
    endif

    if(wflag /= 0)then
      write(6, *)'WND: ', 'There is no wind data input.'
      stop
    endif

  endif

!-------------------------------------------------------------------------------

  if(dtime >= dwtime2 .and. wflag == 0)then

    !write(6, *)'WND: ', datestr(itime)
    dwtime1 = dwtime2
    windx1 = windx2
    windy1 = windy2

    ii = wndfreq * int(itime(4)/wndfreq) + wndfreq
    dwtime2 = datenum([itime(1), itime(2), itime(3), ii, 0, 0])

    ! --- read wind2
    if(wndtype == 0)then
      call read_wind(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 1)then
      call read_wind_gfs(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 2)then
      call read_wind_qbln_interp(dwtime2, windx2, windy2, wflag)
    elseif(wndtype == 3)then
      call read_wind_ncep_interp(dwtime2, windx2, windy2, wflag)
    endif

    if(wflag /= 0)then
      dwtime2 = dwtime1
      windx2 = windx1
      windy2 = windy1
    endif

  endif

!-------------------------------------------------------------------------------

  if(wflag /= 0)then
    wx = windx1
    wy = windy1
    w = sqrt(wx*wx + wy*wy)
    return
  endif

  b = (dtime - dwtime1)/(dwtime2 - dwtime1)
  a = 1 - b

  wx = windx1 * a + windx2 * b
  wy = windy1 * a + windy2 * b
  w = sqrt(abs(wx*wx + wy*wy))

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine get_wind

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: read_wind

  subroutine read_wind(dwtime, windx, windy, flag)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  double precision, intent(in) :: dwtime
  real, intent(out) :: windx(:, :), windy(:, :)
  integer, intent(out) :: flag

!-------------------------------------------------------------------------------

  integer :: ncid, rec, iwtime(6), maxrec
  character(len=14) :: cwtime
  character(len=100) :: filename
  logical :: ext

!-------------------------------------------------------------------------------

  iwtime = datevec(dwtime)
  cwtime = datestr(iwtime)
  filename = trim(wind_path)//'wind_'//cwtime(1:6)//'.nc'

  !write(6, *)'WND: ', trim(filename)

  inquire(file=filename, exist=ext)
  flag = 0
  if(.not. ext)then
    flag = 1
    return
  endif

!-------------------------------------------------------------------------------

!  rec = (dwtime - datenum([iwtime(1), iwtime(2), iwtime(3), 0, 0, 0])) * 24.0 &
  rec = (dwtime - datenum([iwtime(1), iwtime(2), 1, 0, 0, 0])) * 24.0 &
      / wndfreq + 1
      
  !write(6, *)'WND: ', datestr(dwtime), ' ~ ', &
  !           datestr([iwtime(1), iwtime(2), 1, 0, 0, 0]), wndfreq, rec

  call open_nc(ncid, filename, 'r')
  maxrec = get_dimension_len(ncid, 'time')
  if(rec <= maxrec)then
	  call readnc(ncid, 'windx', windx, rec)
	  call readnc(ncid, 'windy', windy, rec)
	else
	  flag = 1
	  write(6, *)'WND: ', 'maxrec, rec', maxrec, rec
	endif
  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine read_wind

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: read_wind_gfs

  subroutine read_wind_gfs(dwtime, windx, windy, flag)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  double precision, intent(in) :: dwtime
  real, intent(out) :: windx(:, :), windy(:, :)
  integer, intent(out) :: flag

!-------------------------------------------------------------------------------

  integer :: ncid, rec, iwtime(6), maxrec
  character(len=14) :: cwtime
  character(len=100) :: filename
  logical :: ext
  
  integer :: witime(6)
  character(len=10) :: wctime

!-------------------------------------------------------------------------------

  iwtime = datevec(dwtime)
  cwtime = datestr(dtime0)
!  filename = trim(wind_path)//'gfs.wind.'//cwtime(1:8)//'00.nc'
  filename = trim(wind_path) 
  !write(6, *)'WND: ', trim(filename)
  
  inquire(file=filename, exist=ext)
  flag = 0
  if(.not. ext)then
    flag = 1
    return
  endif

!-------------------------------------------------------------------------------

  call open_nc(ncid, filename, 'r')
  maxrec = get_dimension_len(ncid, 'time')

  witime = 0
  call readnc(ncid, 'time', wctime, 1)
  read(wctime, '(i4.4,3i2.2)')witime(1:4)
  
  rec = (dwtime - datenum(witime)) * 24. / wndfreq + 1
  !write(6, *)'WND: ', datestr(dwtime), ' ~ ', datestr(witime), wndfreq, rec

  if(rec <= maxrec)then
	  call readnc(ncid, 'uwnd', v2(1:ix2, :), rec)
	  v2(ixl, :) = v2(1, :); windx = v2 * 0.01
	  call readnc(ncid, 'vwnd', v2(1:ix2, :), rec)
	  v2(ixl, :) = v2(1, :); windy = v2 * 0.01
	else
	  flag = 1
	  write(6, *)'WND: ', 'maxrec, rec', maxrec, rec
	endif

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine read_wind_gfs

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: read_wind_qbln_interp

  subroutine read_wind_qbln_interp(dwtime, windx, windy, flag)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  double precision, intent(in) :: dwtime
  real, intent(out) :: windx(:, :), windy(:, :)
  integer, intent(out) :: flag

!-------------------------------------------------------------------------------

  integer :: ncid, rec, iwtime(6), maxrec
  character(len=14) :: cwtime
  character(len=100) :: filename
  logical :: ext

  integer, parameter :: im = 720, jm = 353
  real, parameter :: wdx = 0.5, wdy = 0.5, lon0 = 0.5, lat0 = -88.

  integer :: i, j, i1, i2
  real :: x1, x2, x3, x4, y1, y2, y3, y4, w1, w2, w3, w4, xx, yy

!-------------------------------------------------------------------------------

  iwtime = datevec(dwtime)
  cwtime = datestr(iwtime)
  filename = trim(wind_path)//'wind'//cwtime(1:6)//'.nc'

  !write(6, *)'WND: ', trim(filename)

  inquire(file=filename, exist=ext)
  flag = 0
  if(.not. ext)then
    flag = 1
    return
  endif

!-------------------------------------------------------------------------------

  if(wflag_qbln == 0)then

    allocate(lon(im), lat(jm))
    allocate(wu(im, jm), wv(im, jm))
    allocate(widxx(ixl, 2), widxy(iyl, 2), wcoe(ixl, iyl, 4))

    call open_nc(ncid, filename, 'r')
    call readnc(ncid, 'lon', lon)
    call readnc(ncid, 'lat', lat)
    call close_nc(ncid)
!    lon = lon - lonref

    !write(6, *)'input lon, lat ok.'

    do i = ixs, ixl
!      xx = x(i)
      xx = rx(i)
      if(xx >= 360.0)xx = xx - 360.
      if(xx < 0.0)xx = xx + 360. ! yinxq
      i1 = floor((xx - lon0) / wdx) + 1
      i2 = floor((xx - lon0) / wdx) + 2
      if(i1 < 1 )i1 = im
      if(i2 > im)i2 = 1
      widxx(i, :) = [i1, i2]
    enddo

    do i = iys, iyl
      i1 = floor((y(i) - lat0) / wdy) + 1
      i2 = floor((y(i) - lat0) / wdy) + 2
      if(i1 < 1 )i1 = 1
      if(i2 < 1 )i2 = 1
      if(i1 > jm)i1 = jm
      if(i2 > jm)i2 = jm
      widxy(i, :) = [i1, i2]
    enddo
    
    !write(6, *)'Set index ok.'
    
    do j = iys, iyl
    do i = ixs, ixl

      !xx = x(i)
      xx = rx(i)
      yy = y(j)
      x1 = lon(widxx(i, 1))
      x2 = lon(widxx(i, 2))
      y1 = lat(widxy(j, 1))
      y2 = lat(widxy(j, 2))

      if(x1 > x2)x1 = x1 - 360
      if(xx < x1)xx = xx + 360
      if(xx > x2)xx = xx - 360

      if(y1 == y2)then
        w1 = 0.5 * (x2 - xx) / (x2 - x1)
        w2 = w1; w3 = 0.5 - w1; w4 = w3
      else
        call bilinear_4points(xx, yy, x1, x1, x2, x2,                &
                                      y1, y2, y2, y1, w1, w2, w3, w4 )
      endif
      
      wcoe(i, j, :) = [w1, w2, w3, w4]

!      write(6, *)x1, x2, xx, i
!      write(6, *)y1, y2, yy, j
!      write(6, *)w1, w2, w3, w4
      
    enddo
    enddo

    wflag_qbln = 1
    !write(6, *)'WND: ', datestr(dwtime), ' initial for interp. '

  endif

!-------------------------------------------------------------------------------

!  rec = (dwtime - datenum([iwtime(1), iwtime(2), iwtime(3), 0, 0, 0])) * 24.0 &
  rec = (dwtime - datenum([iwtime(1), iwtime(2), 1, 0, 0, 0])) * 24.0 &
      / wndfreq + 1

  if(iwtime(2) == 2 .and. iwtime(3) == 29)then
    rec = (datenum([iwtime(1), iwtime(2), 28, iwtime(4:6)])      &
          - datenum([iwtime(1), iwtime(2), 1, 0, 0, 0])) * 24.0  &
        / wndfreq + 1
  endif

  !write(6, *)'WND: ', datestr(dwtime), ' ~ ', &
  !           datestr([iwtime(1), iwtime(2), 1, 0, 0, 0]), wndfreq, rec

  call open_nc(ncid, filename, 'r')

  maxrec = get_dimension_len(ncid, 'time')
  if(rec <= maxrec)then
	  call readnc(ncid, 'windu', wu, rec)
	  call readnc(ncid, 'windv', wv, rec)
	else
	  !flag = 1
	  write(6, *)'WND: ', 'maxrec, rec', maxrec, rec
    rec = maxrec
	  call readnc(ncid, 'windu', wu, rec)
	  call readnc(ncid, 'windv', wv, rec)
	endif

  call close_nc(ncid)

!-------------------------------------------------------------------------------

  do j = iys, iyl
    do i = ixs, ixl

      windx(i, j) = wu(widxx(i, 1), widxy(j, 1)) * wcoe(i, j, 1) &
                  + wu(widxx(i, 1), widxy(j, 2)) * wcoe(i, j, 2) &
                  + wu(widxx(i, 2), widxy(j, 2)) * wcoe(i, j, 3) &
                  + wu(widxx(i, 2), widxy(j, 1)) * wcoe(i, j, 4)

      windy(i, j) = wv(widxx(i, 1), widxy(j, 1)) * wcoe(i, j, 1) &
                  + wv(widxx(i, 1), widxy(j, 2)) * wcoe(i, j, 2) &
                  + wv(widxx(i, 2), widxy(j, 2)) * wcoe(i, j, 3) &
                  + wv(widxx(i, 2), widxy(j, 1)) * wcoe(i, j, 4)
  
    enddo
  enddo

  !write(6, *)'WND: ', datestr(dwtime), ' interp is over. '

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine read_wind_qbln_interp

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: read_wind_ncep_interp

  subroutine read_wind_ncep_interp(dwtime, windx, windy, flag)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  double precision, intent(in) :: dwtime
  real, intent(out) :: windx(:, :), windy(:, :)
  integer, intent(out) :: flag

!-------------------------------------------------------------------------------

  integer :: ncid, rec, iwtime(6), maxrec
  character(len=14) :: cwtime
  !modify by xiefan 2016/4/12 路径长度有时候会超过100个字符，因此改为200
  !character(len=100) :: filename, filename1
  character(len=200) :: filename, filename1
  !end modify
  logical :: ext

  integer :: im, jm
  real :: wdx, wdy, lon0, lat0

  integer :: i, j, i1, i2
  real :: x1, x2, x3, x4, y1, y2, y3, y4, w1, w2, w3, w4, xx, yy

!-------------------------------------------------------------------------------

  iwtime = datevec(dwtime)
  cwtime = datestr(iwtime)
  filename  = trim(wind_path)//'uwnd.sig995.'//cwtime(1:4)//'.nc'
  filename1 = trim(wind_path)//'vwnd.sig995.'//cwtime(1:4)//'.nc'
 
  
  !write(6, *)'WND: ', trim(filename)
  inquire(file=filename, exist=ext)
  flag = 0
  if(.not. ext)then
    flag = 1
    return
  endif

!-------------------------------------------------------------------------------

  if(wflag_ncep == 0)then

    call open_nc(ncid, filename, 'r')

    im = get_dimension_len(ncid, 'lon')
    jm = get_dimension_len(ncid, 'lat')

    allocate(lon(im), lat(jm))
    allocate(wu(im, jm), wv(im, jm), ivar(im, jm))
    allocate(widxx(ixl, 2), widxy(iyl, 2), wcoe(ixl, iyl, 4))

    call readnc(ncid, 'lon', lon)
    call readnc(ncid, 'lat', lat)
    call close_nc(ncid)
    !lon = lon - lonref

    wdx = lon(2) - lon(1)
    wdy = lat(2) - lat(1)
    lon0 = lon(1)
    lat0 = lat(1)

    !write(6, *)'input lon, lat ok.'

    do i = ixs, ixl
!      xx = x(i)
      xx = rx(i)
      if(xx >= 360.0)xx = xx - 360.
      i1 = floor((xx - lon0) / wdx) + 1
      i2 = floor((xx - lon0) / wdx) + 2
      if(i1 < 1 )i1 = im
      if(i2 > im)i2 = 1
      widxx(i, :) = [i1, i2]
    enddo

    do i = iys, iyl
      !i1 = floor((y(i) - lat0) / wdy) + 1
      !i2 = floor((y(i) - lat0) / wdy) + 2
      i2 = floor((y(i) - lat0) / wdy) + 1; i1 = i2 + 1
      if(i1 < 1 )i1 = 1
      if(i2 < 1 )i2 = 1
      if(i1 > jm)i1 = jm
      if(i2 > jm)i2 = jm
      widxy(i, :) = [i1, i2]
    enddo
    
    !write(6, *)'Set index ok.'
    
    do j = iys, iyl
    do i = ixs, ixl

      xx = rx(i)
      yy = y(j)
      x1 = lon(widxx(i, 1))
      x2 = lon(widxx(i, 2))
      y1 = lat(widxy(j, 1))
      y2 = lat(widxy(j, 2))

      if(x1 > x2)x1 = x1 - 360
      if(xx < x1)xx = xx + 360
      if(xx > x2)xx = xx - 360

      if(y1 == y2)then
        w1 = 0.5 * (x2 - xx) / (x2 - x1)
        w2 = w1; w3 = 0.5 - w1; w4 = w3
      else
        call bilinear_4points(xx, yy, x1, x1, x2, x2,                &
                                      y1, y2, y2, y1, w1, w2, w3, w4 )
      endif
      
      wcoe(i, j, :) = [w1, w2, w3, w4]
           
    enddo
    enddo

    wflag_ncep = 1
    !write(6, *)'WND: ', datestr(dwtime), ' initial for interp. '

  endif

!-------------------------------------------------------------------------------

  rec = (dwtime - datenum([iwtime(1), 1, 1, 0, 0, 0])) * 24.0 &
      / wndfreq + 1

  !write(6, *)'WND: ', datestr(dwtime), ' ~ ', &
  !           datestr([iwtime(1), 1, 1, 0, 0, 0]), wndfreq, rec

!-------------------------------------------------------------------------------

  call open_nc(ncid, filename, 'r')
  maxrec = get_dimension_len(ncid, 'time')
  if(rec <= maxrec)then
	  call readnc(ncid, 'uwnd', ivar, rec)
    wu = ivar * 0.01 + 225.45
	else
	  flag = 1
	  write(6, *)'WND: ', 'maxrec, rec', maxrec, rec
	endif
  call close_nc(ncid)

  call open_nc(ncid, filename1, 'r')
  maxrec = get_dimension_len(ncid, 'time')
  if(rec <= maxrec)then
	  call readnc(ncid, 'vwnd', ivar, rec)
    wv = ivar * 0.01 + 225.45
	else
	  flag = 1
	  write(6, *)'WND: ', 'maxrec, rec', maxrec, rec
	endif
  call close_nc(ncid)

!-------------------------------------------------------------------------------

  do j = iys, iyl
    do i = ixs, ixl

      windx(i, j) = wu(widxx(i, 1), widxy(j, 1)) * wcoe(i, j, 1) &
                  + wu(widxx(i, 1), widxy(j, 2)) * wcoe(i, j, 2) &
                  + wu(widxx(i, 2), widxy(j, 2)) * wcoe(i, j, 3) &
                  + wu(widxx(i, 2), widxy(j, 1)) * wcoe(i, j, 4)

      windy(i, j) = wv(widxx(i, 1), widxy(j, 1)) * wcoe(i, j, 1) &
                  + wv(widxx(i, 1), widxy(j, 2)) * wcoe(i, j, 2) &
                  + wv(widxx(i, 2), widxy(j, 2)) * wcoe(i, j, 3) &
                  + wv(widxx(i, 2), widxy(j, 1)) * wcoe(i, j, 4)
  
    enddo
  enddo

  !write(6, *)'WND: ', datestr(dwtime), ' interp is over. '

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine read_wind_ncep_interp

!--------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!--------------------------------------------------------------------------------
!*DeckYinxq : bilinear_4points

  subroutine bilinear_4points(x, y, x1, x2, x3, x4, &
                                    y1, y2, y3, y4, &
                                    w1, w2, w3, w4, &
                                    f1, f2, f3, f4, f)

!**********************************************
!----------------------------------------------
! --- Compute the weights for bilinear interpolation.
!     Copied and modified from gomc27_test87f.f, subroutine BLINT.
!     (Xi,Yi,fi) = data grid & values surounding model point (x,y).
!     w# = weights at the model grid points.
!     Recopy=ied from Xiaohua Lin's program.
!                    --- Xunqiang Yin, 2005/08/22
! --- Find 4-points.
!**********************************************
!----------------------------------------------
!
!     x1 < x < x2, x4 < x < x3
!     y1 < y < y4, y2 < y < y3
!
!          (2)---(3)
!           |     |
!           |     |
!          (1)---(4)
!
!----------------------------------------------
!**********************************************

  implicit none

  real, intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y
  real, intent(out) :: w1,w2,w3,w4

  real, optional, intent(in) :: f1, f2, f3, f4
  real, optional, intent(out) :: f

  real :: a1, a2, a3, a4, b1, b2, b3, b4, A, B, C, t, s

  a1 = x1 - x2 + x3 - x4
  a2 = -x1 + x4
  a3 = -x1 + x2
  a4 = x1 - x

  b1 = y1 - y2 + y3 - y4
  b2 = -y1 + y4
  b3 = -y1 + y2
  b4 = y1 - y

  A = a3 * b1 - a1 * b3
  B = b2 * a3 + b1 * a4 - a1 * b4 - a2 * b3
  C = -a2 * b4 + a4 * b2
  if(ABS(A*C).gt.0.002*B**2) then
    t = (-B - sqrt(B*B-4.*A*C)) / (2.*A)
  else
    t = C / ABS(B)
  endif

  A = a2 * b1 - a1 * b2
  B = b3 * a2 + b1 * a4 - a1 * b4 - a3 * b2
  C = -a3 * b4 + a4 * b3
  if(ABS(A*C).gt.0.002*B**2) then
    s = (-B + sqrt(B*B-4.*A*C)) / (2.*A)
  else
    s = - C / ABS(B)
  endif

  w1 = (1. - t) * (1. - s)
  w2 = t * (1. - s)
  w3 = s * t
  w4 = (1. - t) * s

  if(present(f))then
    f = f1 * w1 + f2 * w2 + f3 * w3 + f4 * w4
  endif

  return

  end subroutine bilinear_4points

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
