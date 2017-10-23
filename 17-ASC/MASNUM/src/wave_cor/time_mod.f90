!###############################################################################
!-------------------------------------------------------------------------------
!*DeckYinxq: time_mod 2.1

  module time_mod

!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------
!                                                Copyright (C) 2007 Xunqiang Yin
!                                                MODULE NAME : time_mod
!                                                OLD VERSION : 2005-11-04
!                                                PRESENT VERSION : 2007-10-16
!
! --- USAGE : To convert time among different types: datestr, datevec, datenum.
! --- DEPEND: None
!
! --- NOTE for describing of subroutine / function :
!  A. The parameters bracketed with [], means optional parameter.
!  B. The describe for the parameters of subroutine / function, started with:
!   * It means input prameter;
!   # It means output prameter;
!   @ It means input and output prameter(it will be changed inside).
!
!-------------------------------------------------------------------------------
! ******************************************************************************
! ***                       INTERFACE DESCRIBE                               ***
! ******************************************************************************
!-------------------------------------------------------------------------------
!
!  1. datenum : Get the date in a double precision number.
!
!     datenum(year, month, day, hour, minute, second)
!     datenum(datev[, flag])
!
!    # datenum      = a double precision number --- the counted time from
!                     reference_date. (in days, hours, minutes or seconds)
!    * datev   = a vector contains : yyyy, mm, dd, HH, MM, SS.
!    * flag = the unit type for datenum ( 1 = days, 2 = hours,
!                              3 = minutes & 4 = seconds.)
!
!-------------------------------------------------------------------------------
!
!  2. datevec : Get the date in a vector.
!
!     datevec(daten[, flag])
!
!    # datevec      = a vector contains : yyyy, mm, dd, HH, MM, SS.
!    * daten   = a double precision number --- the counted time from
!                     reference_date. (in days, hours, minutes or seconds)
!    * flag = the unit type for datenum ( 1 = days, 2 = hours,
!                              3 = minutes & 4 = seconds.)
!
!-------------------------------------------------------------------------------
!
!  3. datestr : Get the date in a string.
!
!     datestr(datev)
!     datestr(daten)
!
!    # datestr      = a string contains a date with the formate of datestr_form.
!    * datev   = a vector contains : yyyy, mm, dd, HH, MM, SS.
!    * daten   = a double precision number --- the counted time from
!                     reference_date. (in days, hours, minutes or seconds)
!
!
!-------------------------------------------------------------------------------
!
!  4. init_time_mod : Initialize for time_mod (Set refference date for counting).
!
!     init_time_mod(reference_date_in)
!
!   * reference_date_in = New refference date for count. Default is
!                         [0, 0, 0, 0, 0, 0] (0000-00-00 00:00:00).
!
!                                     --- Written by Xunqiang Yin, 2007-10-16
!
!
!-------------------------------------------------------------------------------
!
! --- For leap year 0, it should check 03.01, some time it will add 1 more day.
!                                                --- Xunqiang Yin, 2007-12-06
!
!-------------------------------------------------------------------------------

  public init_time_mod, datenum, datevec, datestr, dtime_now
  private datenum_days
  private

!-------------------------------------------------------------------------------

  interface datenum
    module procedure datenum1, datenum2, datenum3
  end interface datenum

!-------------------------------------------------------------------------------

  interface datestr
    module procedure datestr0, datestr1
  end interface datestr

!-------------------------------------------------------------------------------

  interface datevec
    module procedure datevec0, datevec1
  end interface datevec

!-------------------------------------------------------------------------------

  integer, parameter :: mon_days(12) = [31, 28, 31, 30, 31, 30, &
                                        31, 31, 30, 31, 30, 31  ]

  integer :: datev_default(6) = [0, 0, 0, 0, 0, 0]

!-------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: init_date_func_mod

  subroutine init_time_mod(datev)

  implicit none

  integer, intent(in) :: datev(6)

  datev_default = datev

  return

  end subroutine init_time_mod

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datenum

  double precision function datenum3(year, month, day, hour, minute, second)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  integer, intent(in) :: year, month, day, hour, minute, second

!-------------------------------------------------------------------------------

  datenum3 = datenum1([year, month, day, hour, minute, second], 1)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datenum3

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datenum

  double precision function datenum2(datev)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  integer, intent(in) :: datev(6)

!-------------------------------------------------------------------------------

  datenum2 = datenum1(datev, 1)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datenum2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datenum

  double precision function datenum1(datev, flag)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  integer, intent(in) :: datev(6), flag

!-------------------------------------------------------------------------------

  datenum1 = datenum_days(datev) - datenum_days(datev_default)

!-------------------------------------------------------------------------------

  if(flag == 1)datenum1 = datenum1
  if(flag == 2)datenum1 = datenum1 * 24.d0
  if(flag == 3)datenum1 = datenum1 * 1440.d0
  if(flag == 4)datenum1 = datenum1 * 86400.d0

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datenum1

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datenum_days

  double precision function datenum_days(datev)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  integer, intent(in) :: datev(6)

!-------------------------------------------------------------------------------

  integer :: leap, year, mon, flag

!-------------------------------------------------------------------------------
! --- for years

  datenum_days = 365 * datev(1)

  leap = 0
  do year = 0, datev(1)-1, 1
    if((mod(year, 4)==0 .and. mod(year, 100) /= 0) .or. (mod(year, 400)==0))then
      leap = leap + 1
    endif
  enddo

  datenum_days = datenum_days + leap

!-------------------------------------------------------------------------------
! --- for month

  do mon = 1, datev(2)-1, 1
    datenum_days = datenum_days + mon_days(mon)
  enddo

  leap = 0
  if(datev(2) > 2)then
    if((mod(year, 4)==0 .and. mod(year, 100) /= 0) .or. (mod(year, 400)==0))then
      leap = 1
    endif
  endif

  datenum_days = datenum_days + leap

!-------------------------------------------------------------------------------

! --- for day
  datenum_days = datenum_days + datev(3)

!-------------------------------------------------------------------------------

! --- for time : hour, minute, second

  datenum_days = datenum_days + datev(4) / 24.d0 + datev(5) / 1440.d0 &
                              + datev(6) / 86400.d0

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datenum_days

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datevec1

  function datevec1(daten, flag)

  implicit none

  double precision, intent(in) :: daten
  integer, intent(in) :: flag
  integer :: datevec1(6)

!-------------------------------------------------------------------------------

  if(flag == 1)datevec1 = datevec0(daten)
  if(flag == 2)datevec1 = datevec0(daten/24.d0)
  if(flag == 3)datevec1 = datevec0(daten/1440.d0)
  if(flag == 4)datevec1 = datevec0(daten/86400.d0)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datevec1

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datevec0

  function datevec0(daten)

  implicit none

  double precision, intent(in) :: daten
  integer :: datevec0(6)

  integer :: rest_days, daysFeb
  double precision :: time_in_1day, daten1

  integer*4 :: dseconds

!-------------------------------------------------------------------------------
! --- Initial value.

  daten1 = daten + datenum_days(datev_default)
  datevec0 = 0

!-------------------------------------------------------------------------------
! --- For year

  datevec0(1) = int(daten1/365.d0)
  rest_days = daten1 - datenum_days(datevec0)
  do while (rest_days < 0)
    datevec0(1) = datevec0(1) - 1
    rest_days = daten1 - datenum_days(datevec0)
  enddo

!-------------------------------------------------------------------------------
! --- For mon & day

  daysFeb = 28
  if((mod(datevec0(1), 4)==0 .and. mod(datevec0(1), 100) /= 0) &
     .or. (mod(datevec0(1), 400)==0))daysFeb = 29

  datevec0(2) = 1
  do while (rest_days > max(mon_days(datevec0(2)), daysFeb))
    datevec0(2) = datevec0(2) + 1
    rest_days = rest_days - max(mon_days(datevec0(2)-1), daysFeb)
  enddo

  datevec0(3) = daten1 - datenum_days(datevec0)

  if(datevec0(3) <= 0 .and. datevec0(2) == 1)then
    datevec0(1) = datevec0(1) - 1
    datevec0(2) = 12
    datevec0(3) = 0
    datevec0(3) = daten1 - datenum_days(datevec0)
  endif

!-------------------------------------------------------------------------------
! --- For time

!!  time_in_1day = daten1 - int(daten1)
!!
!!  datevec0(4) = mod(int(time_in_1day*24.d0   ), 24)
!!  datevec0(5) = mod(int(time_in_1day*1440.d0 ), 60)
!!  datevec0(6) = mod(ceiling(time_in_1day*86400.d0), 60)

  time_in_1day = daten1 - int(daten1)
!  dseconds = ceiling(time_in_1day*86400.d0)   !jiangxj
  dseconds = nint(time_in_1day*86400.d0)       !jiangxj
  if(dseconds >= 86400)then
    datevec0(4:5) = 0
    datevec0(3) = datevec0(3) + 1
    if(datevec0(3) > max(mon_days(datevec0(2)-1), daysFeb))then
      datevec0(3) = 1
      datevec0(2) = datevec0(2) + 1
      if(datevec0(2) > 12)then
        datevec0(2) = 1
        datevec0(1) = datevec0(1) + 1
      endif
    endif
  else
    datevec0(4) = dseconds / 3600
    datevec0(5) = (dseconds - datevec0(4)*3600)/ 60
    datevec0(6) = dseconds - datevec0(4)*3600 - datevec0(5)*60
  endif

  if(datevec0(6) == 1)then
    if(daten - datenum_days(datevec0) < 5.787037037037037d-006)&
    datevec0(6) = 0
  endif

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datevec0

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datestr

  character(len=14) function datestr0(datev)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  integer, intent(in) :: datev(6)

!-------------------------------------------------------------------------------

  write(datestr0, '(i4.4, 5i2.2)')datev

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datestr0

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: datestr1

  character(len=14) function datestr1(daten)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  double precision, intent(in) :: daten

  integer :: datev(6)

!-------------------------------------------------------------------------------

  datev = datevec(daten)
  datestr1 = datestr0(datev)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end function datestr1

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: dtime_now

  double precision function dtime_now()

  integer :: val(8)

  call date_and_time(values=val)

  dtime_now = datenum([val(1:3), val(5:7)])

  return

  end function dtime_now

!-------------------------------------------------------------------------------

  end module time_mod

!-------------------------------------------------------------------------------
!###############################################################################
