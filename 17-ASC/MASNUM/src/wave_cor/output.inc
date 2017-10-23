!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: output

  subroutine output

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  integer :: ncid
  character(len=80) :: filename, file1, file2, file3, file4, fileeng
  character(len=80) :: timstr, gamastr
  logical :: chkflag
  real :: gama_wt

!-------------------------------------------------------------------------------

!  if(mod_init /= 1)call wamfio_mod_init()

!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------

! --- output wave variables into file multi-records, 
!     if outflag =  1  , one file every year,
!                   2  , one file every month,
!                   3  , one file every day,
!                  else, one file every run.

      !chkflag = sum(itime(4:5)) == 0
      !filename = trim(title)//'_wav_'//ctime(1:8)//'.nc'
      !fileeng = trim(title)//'_eng_'//ctime(1:8)//'.nc'
      !if (dtime >= dtimeend - 1.0) call outwav_t(trim(filename), recd)
      !return

  if(iwiofreq > 0 .and. mod(it-number, iwiofreq) == 0)then
   
    chkflag = .false.
    if(outflag == 1)then
      ! --- Output one file by yearly.
      chkflag = sum(itime(4:5)) == 0 .and. itime(3) == 1 .and. itime(2) == 1
      filename = trim(title)//'_wav_'//ctime(1:4)//'.nc'
      fileeng = trim(title)//'_eng_'//ctime(1:4)//'.nc'
    elseif(outflag == 2)then
      ! --- Output one file by monthly.
      chkflag = sum(itime(4:5)) == 0 .and. itime(3) == 1
      filename = trim(title)//'_wav_'//ctime(1:6)//'.nc'
      fileeng = trim(title)//'_eng_'//ctime(1:6)//'.nc'
    elseif(outflag == 3)then
    ! --- Output one file by daily.
      chkflag = sum(itime(4:5)) == 0
      filename = trim(title)//'_wav_'//ctime(1:8)//'.nc'
      fileeng = trim(title)//'_eng_'//ctime(1:8)//'.nc'
    else
      filename = trim(title)//'_wav_'//ctime//'.nc'
      fileeng = trim(title)//'_eng_'//ctime//'.nc'
    endif

    if(chkflag)recd = 0
    recd = recd + 1

    if (dtime >= dtimeend - 1.0) call outwav_t(trim(filename), recd)
    !if(flageng == 1)call outeng_t(trim(fileeng), recd)

  endif



!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------

! --- Output mixing variables. If iciofreq <= 0, no mixing variables output.
  
  if(iciofreq <= 0)return

  ! return any way
  ! don't spend time on IO while testing others
  return;
  
  if(mod(it-number, iciofreq) == 0)then

    if(it == number)then
      ea = ee
    else
      ea = ea / float(iciofreq)
    endif

    if(mixflag <= 3)then

      if(mixflag == 1 .or. mixflag == 3)call mixture
      if(mixflag == 2 .or. mixflag == 3)call mixture_limit

      if(ciofreq < 24)then
        filename = trim(title)//'_mix_'//ctime(1:14)//'.nc'
      else
        filename = trim(title)//'_mix_'//ctime(1:8)//'.nc'
      endif

      call intact

      call outmix(trim(filename))

    else

      !if(mixflag == 4 .or. mixflag == 5)call mixture
      !if(mixflag == 5)call mixture_limit

      if(ciofreq < 24)then
        timstr = ctime(1:14)//'.nc'
      else
        timstr = ctime(1:8)//'.nc'
      endif

      file1 = trim(title)//'_tau_'//trim(timstr)
      file2 = trim(title)//'_bvv_'//trim(timstr)
      file3 = trim(title)//'_bvl_'//trim(timstr)

      if(mixflag == 4)then

        call mixture
        call outmix_bv(trim(file2))

      elseif(mixflag == 5)then

        call intact
        call mixture
        call mixture_limit
        
        call outmix_tau(trim(file1))
        call outmix_bv(trim(file2))
        call outmix_bvl(trim(file3))

        do gama_wt = 1.0, 3.0, 1.0
          write(gamastr, "('gama_', f3.1, '_')")gama_wt
          file4 = trim(title)//'_wit_'//trim(gamastr)//trim(timstr)
          call mixture_wit(gama_wt)
          call outmix_wit(trim(file4))
        enddo
                
      endif

    endif

    ea = 0.0
    ea = ea + ee

  else

    ea = ea + ee

  endif

!-------------------------------------------------------------------------------
!*******************************************************************************
!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine output

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
