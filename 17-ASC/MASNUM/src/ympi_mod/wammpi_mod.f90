!###############################################################################
!-------------------------------------------------------------------------------

  module wammpi_mod

!-------------------------------------------------------------------------------

  use time_mod
  use netcdf_mod
  
  use ympi_mod
  use wamvar_mod
  
  use wamfio_mod, only : wamfio_mod_init, get_wind, my_getwind
  use wamcor_mod, only : setspec, c_init_once, propagat_init_once, mean1
  use wamcpl_mod, only : set_uv, set_ice
  use wamcpl_mod, only : zyyz, bv, bvl, bh1, bh2
  use wamcpl_mod, only : taubb11, taubb12, taubb22, taubb33
  use wamcpl_mod, only : mixture, mixture_limit, intact
  use wamcpl_mod, only : mixture, mixture_limit, intact

  use wamcpl_mod, only : bv_wtv, bv_wtd, mixture_wit  !shenhj 2012-09-23
    
  implicit none

!-------------------------------------------------------------------------------

  public wammpi_init, readwi_mpi, ympi_final

  private
    
!-------------------------------------------------------------------------------

  real, allocatable :: v2(:, :), v3(:, :, :)
  integer, allocatable :: iv2(:, :)
  integer :: recd
! Modified by Zhenya.Song, for point check, 2016/04/06
  !Modified by Huang Chenghuan, 20170422.
  !will be changed in C function.
  integer :: ipoint
  integer :: jpoint
! End modify
  logical :: ext
!-------------------------------------------------------------------------------

  contains
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: readwi_mpi

  subroutine readwi_mpi

!-------------------------------------------------------------------------------

  implicit none
  include 'mpif.h'

  real*8 :: stt, edt = 0.0, ttstt, ttedt
  real*8 :: get_wind_total
  real*8 :: set_uv_total
  real*8 :: set_ice_total
  real*8 :: propagat_halo_total
  real*8 :: propagat_inner_total
  real*8 :: timestep_total
  real*8 :: setspec2_total
  real*8 :: comm1_total
  real*8 :: comm2_total
  real*8 :: wait1_total
  real*8 :: wait2_total
  real*8 :: mean1_time_total
  real*8 :: output_time_total
  real*8 :: copy_back_total
  real*8 :: copy_back2_total
  real*8 :: ppg_inner_wait_total
  real*8 :: others_time_total
  real*8 :: ppg_inner_wait2_total

!  integer :: point_check_file_open_flag = 0
  
  double precision :: my_dtime
  integer :: my_itime(6)

  integer :: ncount, i, command!, print_interval

  real*8 :: collected_times(50)

  get_wind_total = 0.0
  set_uv_total = 0.0
  set_ice_total = 0.0
  propagat_halo_total = 0.0
  propagat_inner_total = 0.0
  timestep_total = 0.0
  setspec2_total = 0.0
  comm1_total = 0.0
  comm2_total = 0.0
  wait1_total = 0.0
  wait2_total = 0.0
  mean1_time_total = 0.0
  output_time_total = 0.0
  copy_back_total = 0.0
  copy_back2_total = 0.0
  ppg_inner_wait_total = 0.0
  ppg_inner_wait2_total = 0.0
  others_time_total = 0.0

!-------------------------------------------------------------------------------
  if (myid == 0) then
    open(371,file="point_check.txt",status="replace")
    close(371)
    print *, "delttm", delttm
  endif


  !point_check_file_open_flag = 0
!  point_check_flag = 0
!  if (mypebox%i1 <= ipoint .and. mypebox%i2 >= ipoint .and. mypebox%j1 <= jpoint .and. mypebox%j2 >= jpoint) then
!      point_check_flag = 1
!      !open(371,file="point_check.txt",access="append")
!      !write(371, *) "Time = ", ctime(1:12),", lon = ", x(ipoint-mypebox%ei1+1), ", lat = ", y(jpoint-mypebox%ej1+1), ", hs = ", h1_3(ipoint-mypebox%ei1+1, jpoint-mypebox%ej1+1)
!      !close(371)
!  endif

! --- Input restart file or set for cool restart.

  !call inprst(rstfile, key)

  !write(6, *) myid, 'inp';

!-------------------------------------------------------------------------------

  dtime0 = datenum(istime)
  dtimeend = datenum(ietime)
  number = (1-key) * cools_days * 1440. / delttm
  itend = number + (datenum(ietime) - dtime0 + 1) * 1440. / delttm

!-------------------------------------------------------------------------------

  key = 0
  dtime = dtime0
  dd=x(ixs+1)-x(ixs)
   do ia=ixs,ixl-1
     if(x(ia+1)-x(ia)/=dd) then
       inputxflag=min(inputxflag,1)
      endif
     if(x(ia+1)<x(ia)) then
       inputxflag=min(inputxflag,0)
      endif
   enddo
   dd=y(iys+1)-y(iys)
    do ic=iys,iyl-1
     if(y(ic+1)-y(ic)/=dd) then
       inputyflag=min(inputyflag,1)
      endif
     if(y(ic+1)<y(ic)) then
       inputyflag=min(inputyflag,0)
      endif
    enddo

!   if (has_c_inited .eq. 0) then
!         call c_init_once()
!         has_c_inited = has_c_inited + 1 
!     endif
  
!      if (has_propagat_inited .eq. 0) then
!          call propagat_init_once()
!       has_propagat_inited = has_propagat_inited + 1
!      endif

!   if (has_c_acce_kernel_inited .eq. 0) then 
! !   call init_c_acce_kernel( &
! !       hor_inner_start, hor_inner_end, ver_inner_start, ver_inner_end, &
! !       hor_inner_start, hor_inner_end, ver_inner_start, ver_inner_end, &
! !                my_hor_start, my_hor_end, my_ver_start, my_ver_end &
! !     )
!         call init_c_acce_kernel( &
!                 my_hor_start + 1, my_hor_end - 1, my_ver_start + 1, my_ver_end - 1, &
!                 halosize * 2 + 1, ixl - halosize * 2, halosize * 2 + 1, iyl - halosize * 2, &
!                 my_hor_start, my_hor_end, my_ver_start, my_ver_end, &
!                 ipoint, jpoint &
!             )
!         e(:, :, :, :) = ee(:, :, :, :)
!     has_c_acce_kernel_inited = has_c_acce_kernel_inited + 1
!   endif

  ttstt = MPI_Wtime()

!print *, "glbflag = ", glbflag, myid

  !Modified by Huang Chenghuan, 20170422.
  !will be changed in init_c_acce_kernel.
  ipoint = 225
  jpoint = 50


  do it = 0, itend 

!-------------------------------------------------------------------------------

! --- Set time & key for cool start.
    dtime = dtime + key * delttm / 1440.
    if(it >= number)key = 1
!    dtime = dtime0 + key * (it - number) * delttm / 1440.
    itime = datevec(dtime)
    ctime = datestr(itime)

    if(dtime > dtimeend)then
      !print *, "direct quit!!!"
      exit
    endif

    command = itend / 50
    command = command + 2
    !ncount = it / 1.8
    if(myid==0 .and. mod(it, command) == 0)write(6, *)'READWI: ', ctime(1:14), delttm, it - number
!  print*, __FILE__, __LINE__
!-------------------------------------------------------------------------------

! --- Read in wind data: wx, wy, w is the wind used in model.

    !stt = MPI_Wtime()
    !if(it == 0) call get_wind
    if (it==0) call get_wind
    !edt = MPI_Wtime() - stt
    get_wind_total = get_wind_total + edt

    !write(6, *) myid, 'get success';

    if(key == 0 .and. it == 0)then
      if(myid==0)write(6, *)'READWI: ', 'Cool start'
      call setspec(1)
    endif

!-------------------------------------------------------------------------------

!modify

  stt = MPI_Wtime()
  call set_uv
  edt = MPI_Wtime() - stt
  set_uv_total = set_uv_total + edt

  inquire(file='ice_clim_mask.nc', exist=ext)
  if(ext) then 
    stt = MPI_Wtime()
    call set_ice
    edt = MPI_Wtime() - stt
    set_ice_total = set_ice_total + edt

    stt = MPI_Wtime()
    call mpi_set_timesteps  ! mpi, yinxq, 2011/10/13 11:45:51
    edt = MPI_Wtime() - stt
    timestep_total = timestep_total + edt
  endif
!end modify

!-------------------------------------------------------------------------------

! --- For wave propagation.

    !write(6, *) myid, 'pro_pre'
    
    stt = MPI_Wtime()
    
    edt = MPI_Wtime() - stt
    set_uv_total = set_uv_total + edt

    stt = MPI_Wtime()
!    if(mypebox%up >= 0) then
!        call c_propagat(my_hor_start, my_hor_end, iyl - 1, iyl - 1)
!    endif
!    if(mypebox%dn >= 0) then
!        call c_propagat(my_hor_start, my_hor_end, 2, 2)
!    endif
!    if(mypebox%nr > 0) call c_propagat(ixl - 1, ixl - 1, ver_inner_start, ver_inner_end)
!    if(mypebox%nl > 0) call c_propagat(2, 2, ver_inner_end, ver_inner_end)

!    call c_propagat(ixs,ixl,halosize+1,halosize*2)
!    call c_propagat(ixs,ixl,iyl-halosize*2+1,iyl-halosize)
!    call c_propagat(halosize+1,halosize*2,iys,iyl)
!    call c_propagat(ixl-halosize*2+1,ixl-halosize,iys,iyl)
!    if(mypebox%dn == MPI_PROC_NULL) call c_propagat(ixs,ixl,1,halosize)
!    if(mypebox%up == MPI_PROC_NULL) call c_propagat(ixs,ixl,iyl-halosize+1,iyl)
!    if(mypebox%nr <= 0) call c_propagat(ixl-halosize+1,ixl,iys,iyl)
!    if(mypebox%nl <= 0) call c_propagat(1,halosize,iys,iyl)
!    call c_propagat(ixs,ixl,halosize*2+1,halosize*3)
!    call c_propagat(ixs,ixl,iyl-halosize*3+1,iyl-halosize*2)
!    call c_propagat(halosize*2+1,halosize*3,iys,iyl)
!    call c_propagat(ixl-halosize*3+1,ixl-halosize*2,iys,iyl)

    if(it == 0)then
      call c_propagat(my_hor_start, my_hor_end, my_ver_start, my_ver_start)
      do i = my_ver_start + 1, my_ver_end - 1
        call c_propagat(my_hor_start, my_hor_start, i, i)
        call c_propagat(my_hor_end, my_hor_end, i, i)
      end do 
      call c_propagat(my_hor_start, my_hor_end, my_ver_end, my_ver_end)
    endif


    call collect_halo_ppg_marine_finalize()
    !!call check_halo_ppg_nsp(1111)
    !!call check_inner_ppg_nsp(1111)
    edt = MPI_Wtime() - stt
    propagat_halo_total = propagat_halo_total + edt         

    stt = MPI_Wtime()
    if(glbflag == 1)then
      !print *, "setspec", myid
      !call setspec(2)
      call c_setspec2(my_hor_start, my_hor_end, my_ver_start, my_ver_end)
    elseif(glbflag == 2)then
      ! --- nest from ...
      !call wamnst_io
    endif
    edt = MPI_Wtime() - stt
    setspec2_total = setspec2_total + edt


    stt = MPI_Wtime()
!     call c_implsch(ixs, ixl, iys + halosize, iys + halosize * 2 - 1);
!     call c_implsch(ixs, ixl, iyl - 2 * halosize + 1, iyl - halosize);
!     call c_implsch(ixs + halosize, ixs + 2 * halosize - 1, iys, iyl);
!     call c_implsch(ixl - 2 * halosize + 1, ixl - halosize, iys, iyl);
!     if(mypebox%dn == MPI_PROC_NULL) call c_implsch(ixs, ixl, iys, iys + halosize - 1);
!     if(mypebox%up == MPI_PROC_NULL) call c_implsch(ixs, ixl, iyl - halosize + 1, iyl);
!     if(mypebox%nr <= 0) call c_implsch(ixl - halosize + 1, ixl, iys, iyl);
!     if(mypebox%nl <= 0) call c_implsch(ixs, ixs + halosize - 1, iys, iyl);

!     call collect_halo_ips_marine_finalize()
    !!call check_halo_ppg_nsp(2222)
    !!call check_inner_ppg_nsp(2222)

    call c_implsch_inner()
    edt = MPI_Wtime() - stt
    ppg_inner_wait_total = ppg_inner_wait_total + edt


    ! if(halosize == 1)call updatev(ee, ixl, iyl, kl, jl, halosize)
    stt = MPI_Wtime()
    if(halosize == 1) call update_communication(e, ixl, iyl, kl, jl, halosize)
    edt = MPI_Wtime() - stt
    comm1_total = comm1_total + edt
    !!call check_halo_ppg_nsp(3333)
    !!call check_inner_ppg_nsp(3333)

    !!!!!call check_halo_ppg_nsp(4444)
    !!call check_inner_ppg_nsp(4444)
    
    stt = MPI_Wtime()
    if(it == 0)then
        command = 1
    else
        command = 3
    endif
    call inner_propagat_begin(command);
    edt = MPI_Wtime() - stt
    propagat_inner_total = propagat_inner_total + edt
  
    stt = MPI_Wtime()
    if(halosize == 1) call update_wait(e, ixl, iyl, kl, jl, halosize)
    edt = MPI_Wtime() - stt
    wait1_total = wait1_total + edt
    !!!!!call check_halo_ppg_nsp(5555)
    !!call check_inner_ppg_nsp(5555)

    !stt = MPI_Wtime()
    !call get_wind
    if (it<itend) then
        my_dtime = dtime + key * delttm / 1440.
        my_itime = datevec(my_dtime)
        if (my_dtime <= dtimeend)call my_getwind(my_dtime,my_itime)
    endif
    !edt = MPI_Wtime() - stt
    !get_wind_total = get_wind_total + edt

    stt = MPI_Wtime()

    if(it == 0)then
        command = 1
    else
        command = 3
    endif
    call inner_propagat_stop(command);
    edt = MPI_Wtime() - stt
    ppg_inner_wait_total = ppg_inner_wait_total + edt
    !!call check_halo_ppg_nsp(6666)
    !!call check_inner_ppg_nsp(6666)

    stt = MPI_Wtime()
    if(halosize == 1) call update_copy_back(e, ixl, iyl, kl, jl, halosize)
    edt = MPI_Wtime() - stt
    copy_back_total = copy_back_total + edt
    !!call check_halo_ppg_nsp(7777)
    !!call check_inner_ppg_nsp(7777)

!-------------------------------------------------------------------------------

! --- For reginal model: setspec or nesting.


!-------------------------------------------------------------------------------

    ! mpi, yinxq, 2011/10/13 11:45:51
   ! if(halosize == 1)call updatev(ee, ixl, iyl, kl, jl, halosize)

    stt = MPI_Wtime()

    !call c_getem()
    !call wait_mean1_inner()
    !call mean1_halo
    !!call check_halo_ppg_nsp(8888)
    !!call check_inner_ppg_nsp(8888)

    edt = MPI_Wtime() - stt
    mean1_time_total = mean1_time_total + edt


    stt = MPI_Wtime()
    call mean1_inner
    edt = MPI_Wtime() - stt
    mean1_time_total = mean1_time_total + edt
    stt = MPI_Wtime()
    call wait_mean1_inner
    edt = MPI_Wtime() - stt
    mean1_time_total = mean1_time_total + edt

    ! mpi, yinxq, 2011/10/13 11:45:51
!    stt = MPI_Wtime()
!    call updatev(ee, ixl, iyl, kl, jl, halosize)
!    edt = MPI_Wtime() - stt
!    updatev_wait_total = updatev_wait_total + edt

    stt = MPI_Wtime()
    if(halosize == 1) call update_communication(ee, ixl, iyl, kl, jl, halosize)
    edt = MPI_Wtime() - stt
    comm2_total = comm2_total + edt
    !!call check_halo_ppg_nsp(9999)
    !!call check_inner_ppg_nsp(9999)

    stt = MPI_Wtime()
    command = 2
    call inner_propagat_begin(command);
    edt = MPI_Wtime() - stt
    propagat_inner_total = propagat_inner_total + edt

  
    stt = MPI_Wtime()
    if(halosize == 1) call update_wait(ee, ixl, iyl, kl, jl, halosize)
    edt = MPI_Wtime() - stt
    wait2_total = wait2_total + edt
    !!call check_halo_ppg_nsp(2233)
    !!call check_inner_ppg_nsp(2233)

    ! Modified by Zhenya.Song, for point check, 2016/04/06
    stt = MPI_Wtime()
    if (mypebox%i1 <= ipoint .and. mypebox%i2 >= ipoint .and. mypebox%j1 <= jpoint .and. mypebox%j2 >= jpoint) then
      if(it == 0) open(371,file="point_check.txt",access="append")
      !open(371,file="point_check.txt",access="append")
      write(371, *) "Time = ", ctime(1:12),", lon = ", x(ipoint-mypebox%ei1+1), ", lat = ", y(jpoint-mypebox%ej1+1), ", hs = ", h1_3(ipoint-mypebox%ei1+1, jpoint-mypebox%ej1+1)
      !close(371)
    endif
    edt = MPI_Wtime() - stt
    set_ice_total = set_ice_total + edt
    ! End modify  

    stt = MPI_Wtime()
    command = 2
    call inner_propagat_stop(command);
    edt = MPI_Wtime() - stt
    ppg_inner_wait2_total = ppg_inner_wait2_total + edt

    stt = MPI_Wtime()
    if(halosize == 1) call update_copy_back(ee, ixl, iyl, kl, jl, halosize)
    edt = MPI_Wtime() - stt
    copy_back2_total = copy_back2_total + edt
    !!call check_halo_ppg_nsp(12450)
    !!call check_inner_ppg_nsp(12450)

!-------------------------------------------------------------------------------

! --- Compute for wave parameters.

!-------------------------------------------------------------------------------  

!-------------------------------------------------------------------------------

! --- Output & save restart.




    if(key == 0)cycle

!     if(myid==0)write(6, *)'output'
!     stt = MPI_Wtime()
!     call output
!     edt = MPI_Wtime() - stt
!     output_time_total = output_time_total + edt
!     if(myid==0)write(6, *)'output over'

!    if(itime(2) == 1 .and. itime(3) == 1)then
!    if(itime(3) == 1)then
!      call outrst('restart_'//ctime(1:12)//'.nc')
!    else
!      call outrst('restart_'//ctime(1:8)//'.nc')
!    endif
    

!-------------------------------------------------------------------------------

  enddo

!     stt = MPI_Wtime()
!     call output
!     edt = MPI_Wtime() - stt
!     output_time_total = output_time_total + edt

  it = it - 1
  ttedt = MPI_Wtime() - ttstt
  others_time_total = others_time_total + ttedt

  if (mypebox%i1 <= ipoint .and. mypebox%i2 >= ipoint .and. mypebox%j1 <= jpoint .and. mypebox%j2 >= jpoint) then
      close(371)
  endif

  collected_times(1) = get_wind_total ! / it
  collected_times(2) = set_uv_total ! / it
  collected_times(3) = set_ice_total ! / it
  collected_times(4) = propagat_halo_total ! / it
  collected_times(5) = propagat_inner_total ! / it
  collected_times(6) = timestep_total ! / it
  collected_times(7) = setspec2_total ! / it
  collected_times(8) = mean1_time_total ! / it
  collected_times(9) = comm1_total ! / it
  collected_times(10) = wait1_total ! / it
  collected_times(11) = comm2_total ! / it
  collected_times(12) = wait2_total ! / it
  collected_times(13) = output_time_total ! / it
  collected_times(14) = copy_back_total ! / it
  collected_times(15) = ppg_inner_wait_total ! / it
  collected_times(16) = copy_back2_total ! / it
  collected_times(17) = ppg_inner_wait2_total ! / it
  collected_times(18) = others_time_total ! / it

  call c_athread_finalize()
  call c_time_print(collected_times)


  !Modified by Huang Chenghuan, 20170422.
  !since Propagat and Implsch is combined, we get the Propagat and Implsch time from C function
  !Confirmed by the ASC TechSupport, the time below is only for reference
  if(myid==0)then
      open(unit=11, file="elapse_time.txt", status = "replace")
      write(11,*)"Propagat Time is ",collected_times(19)," seconds."
      write(11,*)"Implsch Time is ",collected_times(20)," seconds."
      write(11,*)"DMA Time is ",collected_times(21)," seconds."
      close(11)
      print *, "file written"
  endif

  !call outrst('restart_'//ctime(1:12)//'.nc')

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine readwi_mpi

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: wammpi_init

  subroutine wammpi_init

  implicit none
        
  integer :: ncid, i, i1, i2
  character(len=80) :: filename

  real, allocatable :: lon1(:), lat1(:), mask1(:, :), dep(:), top1(:, :)
  real, allocatable :: lon2(:), lat2(:)
  
  real :: addlon
  
!-------------------------------------------------------------------------------

  call ympi_init

!-------------------------------------------------------------------------------

! --- Set parameters & pebox.

  halosize = 1
  
  open(11, file='ctlparams')
  read(11, nml=ctlparams)
  close(11)
  
  filename = trim(data_path)//'wamyyz.nc'

  call init_pebox(filename, 'lon', 'lat', 'nspyyz', npe, myid)
                    
!-------------------------------------------------------------------------------

! --- Input topo, x, y, depyyz, ...

  ixl = mypebox%isize
  iyl = mypebox%jsize

  !write(6, *)'read topo ...'

  call opennc(ncid, trim(data_path)//'wamyyz.nc', 'r')
  kb = get_dimension_len(ncid, 'zyyz')
  allocate(lon2(im), lat2(jm))
  allocate(lon1(ixl), lat1(iyl), mask1(ixl, iyl), dep(kb), top1(ixl, iyl))
  call readnc(ncid, 'lon', lon2)
  call readnc(ncid, 'lat', lat2)
  call readnc(ncid, 'zyyz', dep)
  call readnc(ncid, 'nspyyz', mask1(is:ie, js:je), locs=loc2)
  call readnc(ncid, 'depyyz', top1(is:ie, js:je), locs=loc2)
  call closenc(ncid)

  !write(6, *)'read topo is okay.'

  lat1(:) = lat2(mypebox%ej1:mypebox%ej2)
  
  if(flagxcycle == 1)then
      
    lon1(:) = lon2(mypebox%ei1:mypebox%ei2)
    
  else
      
      do i = mypebox%ei1, mypebox%ei2
          
        i1 = i - mypebox%ei1 + 1
        
        if(i < 1)then
            i2 = i + (im - ixoverlay)
            addlon = -360.
            !write(6, *)lon2(i2), i2, i, im, ixoverlay
        elseif(i > im)then
            i2 = i - (im - ixoverlay)
            addlon = 360.
        else
            i2 = i        
            addlon = 0.
        endif
        
          lon1(i1) = lon2(i2) + addlon

      enddo
      
  endif

  call updatev(mask1, ixl, iyl, halosize)
  call updatev(top1, ixl, iyl, halosize)
  
!-------------------------------------------------------------------------------

! --- wamfio_mod_init

  call wamfio_mod_init(lon1, lat1, mask1, dep, top1)

! --- Set timestep

  call mpi_set_timesteps

  deallocate(lon1, lat1, mask1, dep, top1, lon2, lat2)

  call init_main(myid)

  call c_init_once()
  
  call propagat_init_once()

  call init_c_acce_kernel( &
          my_hor_start + 1, my_hor_end - 1, my_ver_start + 1, my_ver_end - 1, &
          halosize * 2 + 1, ixl - halosize * 2, halosize * 2 + 1, iyl - halosize * 2, &
          my_hor_start, my_hor_end, my_ver_start, my_ver_end, &
          ipoint, jpoint &
      )
  e(:, :, :, :) = ee(:, :, :, :)

!  print*, __FILE__, __LINE__
  
!-------------------------------------------------------------------------------

  end subroutine wammpi_init
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: mpi_set_timesteps

  subroutine mpi_set_timesteps
      
  implicit none
      
! --- Set timestep

  call get_mpimin(delttm)

  deltt    = delttm  * 60.  ! yinxq
  deltt5   = delttm  * 30.  ! yinxq
  iwiofreq = wiofreq * 60. / delttm
  iciofreq = ciofreq * 60. / delttm
  irstfreq = rstfreq * 60. / delttm
  number   = (1-key) * cools_days * 1440. / delttm

  !if(myid==0)write(6, nml=ctlparams1)
  !write(6, *)'mpi_set_timesteps:', delttm
  
  end subroutine mpi_set_timesteps

  include '../wave_cor/output.inc'
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix

  subroutine outmix(filename)

  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: filename

!-------------------------------------------------------------------------------

  integer :: ncid
  
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real(4),allocatable :: recvBuf(:)
  real(4),allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real(4),allocatable :: v3_global(:,:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
!-------------------------------------------------------------------------------

  if(myid == 0)then
    
      call open_nc(ncid, filename, 'c')
    
      call dimension_define(ncid, 'lon', im, 'lon', nf_real)
      call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
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
        
      endif
    
      call end_define(ncid)

      call writenc(ncid, 'lon', lon)
      call writenc(ncid, 'lat', lat)
      call writenc(ncid, 'dep', zyyz)
    
      call close_nc(ncid)

  endif

!-------------------------------------------------------------------------------

 !-------------------------------------------------------------------------------
  
  !分配sendBuf
  sendLengthNoOverlay = (ie - is + 1) * (je - js + 1) * kb
  sendLength = sendLengthNoOverlay
  if(ixoverlay /= 0 .and. loc2(1) == 1) then
    sendLength = sendLength + (ixoverlay + 1)*(je - js + 1) * kb
  end if
  allocate(sendBuf(sendLength))

  totalCounts = 0
  !write(6,*) 'i am in tau my id is:',myid,'im is',im,'jm is',jm

!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
if(myid == 0) then
    do i=1,npe
        length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
        length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
        counts(i) = length_i * length_j*kb
        if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
            counts(i) = counts(i) + (ixoverlay + 1)*length_j*kb
        end if
        totalCounts = totalCounts + counts(i)
        if(i == 1) then
            startPos(i) = 0
        else
            startPos(i) = startPos(i-1) + counts(i-1)
        end if
    end do
    !debug
    !write(*,*) 'total coutns:',totalCounts
    !分配recvBuf
    allocate(recvBuf(totalCounts))
    allocate(v3_global(im,jm,kb))
end if

!A部分结果
    call setland_v3(taubb11)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau11', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!B部分结果
    call setland_v3(taubb12)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau12', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!C部分结果
    call setland_v3(taubb22)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau22', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!D部分结果
    call setland_v3(taubb33)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau33', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!E部分结果    
    if(mixflag == 1 .or. mixflag == 3)then
        call setland_v3(bv)
        call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
        if(ixoverlay /= 0 .and. loc2(1) == 1)then 
            call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
        end if
        !调用gatherV
        call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

        if(myid == 0) then !0号进程将recvBuf中的数据转换广且写兿
            do i=1,npe
                i_start = pebox(2,i) 
                i_end = pebox(3,i)
                j_start = pebox(4,i)
                j_end = pebox(5,i)
                if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                    length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                    tmp = (ixoverlay + 1)*length_j*kb
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                    call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
                else
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
                end if
            end do
            call open_nc(ncid, filename, 'w')
            call writenc(ncid, 'bv', v3_global(:,:,:), locs=[1,1,1])
            call close_nc(ncid)
        end if

    endif
!F部分结果
    if(mixflag == 2 .or. mixflag == 3)then

        call setland_v3(bvl)
        call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
        if(ixoverlay /= 0 .and. loc2(1) == 1)then 
            call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
        end if
        !调用gatherV
        call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

        if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
            do i=1,npe
                i_start = pebox(2,i) 
                i_end = pebox(3,i)
                j_start = pebox(4,i)
                j_end = pebox(5,i)
                if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                    length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                    tmp = (ixoverlay + 1)*length_j*kb
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                    call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
                else
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
                end if
            end do
            call open_nc(ncid, filename, 'w')
            call writenc(ncid, 'bvl', v3_global(:,:,:), locs=[1,1,1])
            call close_nc(ncid)
        end if
        
        call setland_v3(bh1)
        call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
        if(ixoverlay /= 0 .and. loc2(1) == 1)then 
            call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
        end if
        !调用gatherV
        call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

        if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
            do i=1,npe
                i_start = pebox(2,i) 
                i_end = pebox(3,i)
                j_start = pebox(4,i)
                j_end = pebox(5,i)
                if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                    length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                    tmp = (ixoverlay + 1)*length_j*kb
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                    call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
                else
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
                end if
            end do
            call open_nc(ncid, filename, 'w')
            call writenc(ncid, 'bh1', v3_global(:,:,:), locs=[1,1,1])
            call close_nc(ncid)
        end if
    
        call setland_v3(bh2)
        call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
        if(ixoverlay /= 0 .and. loc2(1) == 1)then 
            call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
        end if
        !调用gatherV
        call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

        if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
            do i=1,npe
                i_start = pebox(2,i) 
                i_end = pebox(3,i)
                j_start = pebox(4,i)
                j_end = pebox(5,i)
                if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                    length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                    tmp = (ixoverlay + 1)*length_j*kb
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                    call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
                else
                    call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
                end if
            end do
            call open_nc(ncid, filename, 'w')
            call writenc(ncid, 'bh2', v3_global(:,:,:), locs=[1,1,1])
            call close_nc(ncid)
        end if
      endif
      
      deallocate(sendBuf)
      if(myid == 0) then
        deallocate(recvBuf)
        deallocate(v3_global)
      end if


!-------------------------------------------------------------------------------

  return


!-------------------------------------------------------------------------------

  end subroutine outmix

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix_wit

  subroutine outmix_wit(filename)

  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: filename

!-------------------------------------------------------------------------------

  integer :: ncid
  
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real(4),allocatable :: recvBuf(:)
  real(4),allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real(4),allocatable :: v3_global(:,:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
!-------------------------------------------------------------------------------

  if(myid == 0)then
    
      call open_nc(ncid, filename, 'c')
    
      call dimension_define(ncid, 'lon', im, 'lon', nf_real)
      call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
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

      call writenc(ncid, 'lon', lon)
      call writenc(ncid, 'lat', lat)
      call writenc(ncid, 'dep', zyyz)
    
      call close_nc(ncid)

  endif

!-------------------------------------------------------------------------------

   !分配sendBuf
  sendLengthNoOverlay = (ie - is + 1) * (je - js + 1) * kb
  sendLength = sendLengthNoOverlay
  if(ixoverlay /= 0 .and. loc2(1) == 1) then
    sendLength = sendLength + (ixoverlay + 1)*(je - js + 1) * kb
  end if
  allocate(sendBuf(sendLength))

  totalCounts = 0


!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
    if(myid == 0) then
        do i=1,npe
            length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
            length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
            counts(i) = length_i * length_j*kb
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                counts(i) = counts(i) + (ixoverlay + 1)*length_j*kb
            end if
            totalCounts = totalCounts + counts(i)
            if(i == 1) then
                startPos(i) = 0
            else
                startPos(i) = startPos(i-1) + counts(i-1)
            end if
        end do
        !debug
        !write(*,*) 'total coutns:',totalCounts
        !分配recvBuf
        allocate(recvBuf(totalCounts))
        allocate(v3_global(im,jm,kb))
    end if
    
    !分配完毕，开始gather
    !A部分结果
    call setland_v3(bv_wtv)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'bv_wtv', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
    
    !B部分结果
    call setland_v3(bv_wtd)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'bv_wtd', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
    
    deallocate(sendBuf)
    if(myid == 0) then
        deallocate(recvBuf)
        deallocate(v3_global)
    end if
    
!-------------------------------------------------------------------------------

!  call opennc(ncid, filename, 'w')
  
!  call setland_v3(bv_wtv)                                              
!  call writenc(ncid, 'bv_wtv'   ,v3(is:ie, js:je, :), locs=loc3)       
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then                            
!      call writenc(ncid, 'bv_wtv', v3(is:is+ixoverlay, js:je, :), &      
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])              
!  endif                                                                
                                                                       
!  call setland_v3(bv_wtd)                                              
!  call writenc(ncid, 'bv_wtd'   , v3(is:ie, js:je, :), locs=loc3)      
! if(ixoverlay /= 0 .and. loc2(1) == 1)then                            
!      call writenc(ncid, 'bv_wtd', v3(is:is+ixoverlay, js:je, :), &      
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])              
!  endif                                                                
      
 ! call closenc(ncid)

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outmix_wit

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix_tau

  subroutine outmix_tau(filename)

  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: filename

!-------------------------------------------------------------------------------

  integer :: ncid
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real,allocatable :: recvBuf(:)
  real,allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real,allocatable :: v3_global(:,:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
!-------------------------------------------------------------------------------

  if(myid == 0)then
      
      call open_nc(ncid, filename, 'c')
    
      call dimension_define(ncid, 'lon', im, 'lon', nf_real)
      call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
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
    
      call writenc(ncid, 'lon', lon)
      call writenc(ncid, 'lat', lat)
      call writenc(ncid, 'dep', zyyz)
    
      call close_nc(ncid)

  endif
    
!-------------------------------------------------------------------------------

 !分配sendBuf
  sendLengthNoOverlay = (ie - is + 1) * (je - js + 1) * kb
  sendLength = sendLengthNoOverlay
  if(ixoverlay /= 0 .and. loc2(1) == 1) then
    sendLength = sendLength + (ixoverlay + 1)*(je - js + 1) * kb
  end if
  allocate(sendBuf(sendLength))

  totalCounts = 0


!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
if(myid == 0) then
    do i=1,npe
        length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
        length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
        counts(i) = length_i * length_j*kb
        if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
            counts(i) = counts(i) + (ixoverlay + 1)*length_j*kb
        end if
        totalCounts = totalCounts + counts(i)
        if(i == 1) then
            startPos(i) = 0
        else
            startPos(i) = startPos(i-1) + counts(i-1)
        end if
    end do
    !debug
    !write(*,*) 'total coutns:',totalCounts
    !分配recvBuf
    allocate(recvBuf(totalCounts))
    allocate(v3_global(im,jm,kb))
end if

!A部分结果
    call setland_v3(taubb11)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau11', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!B部分结果
    call setland_v3(taubb12)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau12', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!C部分结果
    call setland_v3(taubb22)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau22', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!D部分结果
    call setland_v3(taubb33)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tau33', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
    
    deallocate(sendBuf)
    if(myid == 0) then
        deallocate(recvBuf)
        deallocate(v3_global)
    end if

!  call opennc(ncid, filename, 'w')

!  call setland_v3(taubb11)
!  call writenc(ncid, 'tau11'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'tau11', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call setland_v3(taubb12)
!  call writenc(ncid, 'tau12'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'tau12', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call setland_v3(taubb22)
!  call writenc(ncid, 'tau22'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'tau22', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call setland_v3(taubb33)
!  call writenc(ncid, 'tau33'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'tau33', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call closenc(ncid)

!-------------------------------------------------------------------------------

  return


!-------------------------------------------------------------------------------

  end subroutine outmix_tau

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix_bvl

  subroutine outmix_bvl(filename)

  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: filename

!-------------------------------------------------------------------------------

  integer :: ncid
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real,allocatable :: recvBuf(:)
  real,allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real,allocatable :: v3_global(:,:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp

!-------------------------------------------------------------------------------

  if(myid == 0)then
      
      call open_nc(ncid, filename, 'c')
    
      call dimension_define(ncid, 'lon', im, 'lon', nf_real)
      call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
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
    
      call writenc(ncid, 'lon', lon)
      call writenc(ncid, 'lat', lat)
      call writenc(ncid, 'dep', zyyz)
    
      call close_nc(ncid)

  endif
    
!-------------------------------------------------------------------------------

      !分配sendBuf
  sendLengthNoOverlay = (ie - is + 1) * (je - js + 1) * kb
  sendLength = sendLengthNoOverlay
  if(ixoverlay /= 0 .and. loc2(1) == 1) then
    sendLength = sendLength + (ixoverlay + 1)*(je - js + 1) * kb
  end if
  allocate(sendBuf(sendLength))

  totalCounts = 0


!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
if(myid == 0) then
    do i=1,npe
        length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
        length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
        counts(i) = length_i * length_j*kb
        if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
            counts(i) = counts(i) + (ixoverlay + 1)*length_j*kb
        end if
        totalCounts = totalCounts + counts(i)
        if(i == 1) then
            startPos(i) = 0
        else
            startPos(i) = startPos(i-1) + counts(i-1)
        end if
    end do
    !debug
    !write(*,*) 'total coutns:',totalCounts
    !分配recvBuf
    allocate(recvBuf(totalCounts))
    allocate(v3_global(im,jm,kb))
end if

!A部分结果
    call setland_v3(bvl)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'bvl', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!B部分结濿
    call setland_v3(bh1)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'bh1', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
!C部分结果
    call setland_v3(bh2)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'bh2', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
    
    deallocate(sendBuf)
    if(myid == 0) then
        deallocate(recvBuf)
        deallocate(v3_global)
    end if
!  call opennc(ncid, filename, 'w')

!  call setland_v3(bvl)
!  call writenc(ncid, 'bvl'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'bvl', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call setland_v3(bh1)
!  call writenc(ncid, 'bh1'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'bh1', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call setland_v3(bh2)
!  call writenc(ncid, 'bh2'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'bh2', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call closenc(ncid)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outmix_bvl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outmix_bv

  subroutine outmix_bv(filename)

  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: filename

!-------------------------------------------------------------------------------

  integer :: ncid
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real,allocatable :: recvBuf(:)
  real,allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real,allocatable :: v3_global(:,:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
!-------------------------------------------------------------------------------

  if(myid == 0)then
      
      call open_nc(ncid, filename, 'c')
    
      call dimension_define(ncid, 'lon', im, 'lon', nf_real)
      call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
      call dimension_define(ncid, 'dep', kb, 'dep', nf_real)
      call set_attribute(ncid, 'units', 'degrees_north', 'lat')
      call set_attribute(ncid, 'units', 'degrees_east', 'lon')
      call set_attribute(ncid, 'modulo', '', 'lon')
      call set_attribute(ncid, 'ctime', ctime)
    
      call variable_define(ncid, 'bv'   , nf_real, ['lon', 'lat', 'dep'])
      call set_attribute(ncid, 'missing_value', nf_fill_real, 'bv')
    
      call end_define(ncid)
    
      call writenc(ncid, 'lon', lon)
      call writenc(ncid, 'lat', lat)
      call writenc(ncid, 'dep', zyyz)
    
      call close_nc(ncid)

  endif
    
!-------------------------------------------------------------------------------

       !分配sendBuf
  sendLengthNoOverlay = (ie - is + 1) * (je - js + 1) * kb
  sendLength = sendLengthNoOverlay
  if(ixoverlay /= 0 .and. loc2(1) == 1) then
    sendLength = sendLength + (ixoverlay + 1)*(je - js + 1) * kb
  end if
  allocate(sendBuf(sendLength))

  totalCounts = 0


!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
if(myid == 0) then
    do i=1,npe
        length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
        length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
        counts(i) = length_i * length_j*kb
        if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
            counts(i) = counts(i) + (ixoverlay + 1)*length_j*kb
        end if
        totalCounts = totalCounts + counts(i)
        if(i == 1) then
            startPos(i) = 0
        else
            startPos(i) = startPos(i-1) + counts(i-1)
        end if
    end do
    !debug
    !write(*,*) 'total coutns:',totalCounts
    !分配recvBuf
    allocate(recvBuf(totalCounts))
    allocate(v3_global(im,jm,kb))
end if

!A部分结果
    call setland_v3(bv)
    call array3DTo1DReal(v3(is:ie, js:je, :),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array3DTo1DReal(v3(is:is+ixoverlay, js:je, :),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kb
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v3_global(i_start:i_end,j_start:j_end,:))
                call array1DTo3DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v3_global(im-ixoverlay:im,j_start:j_end,:))
            else
                call array1DTo3DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v3_global(i_start:i_end,j_start:j_end,:))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'bv', v3_global(:,:,:), locs=[1,1,1])
        call close_nc(ncid)
    end if
    
    deallocate(sendBuf)
      if(myid == 0) then
        deallocate(recvBuf)
        deallocate(v3_global)
    end if
!  call opennc(ncid, filename, 'w')

!  call setland_v3(bv)
!  call writenc(ncid, 'bv'   , v3(is:ie, js:je, :), locs=loc3)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'bv', v3(is:is+ixoverlay, js:je, :), &
!                   locs = [im-ixoverlay, loc3(2), loc3(3)])
!  endif

!  call closenc(ncid)

!-------------------------------------------------------------------------------

  return
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------

  end subroutine outmix_bv

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outeng_t

  subroutine outeng_t(filename, rec)

  implicit none
  include 'mpif.h'
  
  character(len=*), intent(in) :: filename
  integer, intent(in) :: rec

!-------------------------------------------------------------------------------

  integer :: ncid, jd
  integer*2, parameter :: ivland = nf_fill_int2

  logical :: ext
  integer :: outrecord, timerec, ittt
  double precision, allocatable :: timealready(:)
  
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real(4),allocatable :: recvBuf(:)
  real(4),allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real(4),allocatable :: v2_global(:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
  
  real(8) :: recvTimeBuf(npe)
  real(8) :: sendTime
!-------------------------------------------------------------------------------

  jd = datenum(1950, 1, 1, 0, 0, 0)

!-------------------------------------------------------------------------------

  if(myid == 0)then

      inquire(file=filename, exist=ext)
    
      if(.not.ext)then
        
        outrecord = 1
    
        call open_nc(ncid, filename, 'c')
    
        call dimension_define(ncid, 'lon', im, 'lon', nf_real)
        call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
        call dimension_define(ncid, 'time', 0, 'time', nf_real)
        call set_attribute(ncid, 'units', 'degrees_north', 'lat')
        call set_attribute(ncid, 'units', 'degrees_east', 'lon')
        call set_attribute(ncid, 'modulo', '', 'lon')
        call set_attribute(ncid, 'units', 'Days since 1950-01-01 00:00:0.0.', 'time')
        call set_attribute(ncid, 'Start_time', ctime)
    
        call variable_define(ncid, 'pein',    nf_real, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'pebo',    nf_real, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'peds',    nf_real, ['lon', 'lat', 'time'])
    
        call set_attribute(ncid, 'missing_value', nf_fill_real, 'pein')
        call set_attribute(ncid, 'missing_value', nf_fill_real, 'pebo')
        call set_attribute(ncid, 'missing_value', nf_fill_real, 'peds')
    
        call end_define(ncid)
      
        call writenc(ncid, 'lon', lon)
        call writenc(ncid, 'lat', lat)
        
        call close_nc(ncid)
    
      else
    
        call open_nc(ncid, filename, 'w')
        
        timerec = get_dimension_len(ncid, 'time')
        allocate(timealready(timerec))
        call readnc(ncid, 'time', timealready)
        
        outrecord = 1
        do ittt = 1, timerec !-1
            if(dtime-jd > timealready(ittt))then
              outrecord = ittt + 1
            endif
        enddo
        
        deallocate(timealready)
        
        call close_nc(ncid)
    
      endif
      
    endif ! if(myid == 0)then

  call bcast_to_all(outrecord)

!-------------------------------------------------------------------------------

  write(6, *)'eng '
  !call opennc(ncid, filename, 'w')

  !call writenc(ncid, 'time', dtime-jd, outrecord)
  !call closenc(ncid)
  sendTime = dtime - jd
  call mpi_gather(sendTime,1,MPI_DOUBLE_PRECISION,recvTimeBuf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(myid == 0) then
      call open_nc(ncid, filename, 'w')
    do i=1,npe
        call writenc(ncid, 'time', recvTimeBuf(i), outrecord)
    end do
    call close_nc(ncid)
  end if
  
  !分配sendBuf
    sendLengthNoOverlay = (ie - is + 1) * (je - js + 1)
    sendLength = sendLengthNoOverlay
    if(ixoverlay /= 0 .and. loc2(1) == 1)then
        sendLength = sendLength + (ixoverlay + 1)*(je - js + 1)
    end if
    allocate(sendBuf(sendLength))

    totalCounts = 0

!  if(ixoverlay /= 0 .and. loc2(1) == 1)then 
!      call writenc(ncid, 'windx', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
    if(myid == 0) then
        do i=1,npe
            length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
            length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
            counts(i) = length_i * length_j
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                counts(i) = counts(i) + (ixoverlay + 1)*length_j
            end if
            totalCounts = totalCounts + counts(i)
            if(i == 1) then
                startPos(i) = 0
            else
                startPos(i) = startPos(i-1) + counts(i-1)
            end if
        end do
        !debug
        !write(*,*) 'total coutns:',totalCounts
        !分配recvBuf
        allocate(recvBuf(totalCounts))
        allocate(v2_global(im,jm))
    end if

    call opennc(ncid, filename, 'w')
    call writenc(ncid, 'time', dtime-jd, outrecord)
    call closenc(ncid)
!各个进程首先重置iv2的值，然后把iv2的数据放到sendBuf丿

!  A部分结果
    call setland_v2(pein)
    call array2DTo1DReal(v2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DReal(v2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'pein', v2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
!  B部分结果
    call setland_v2(pebo)
    call array2DTo1DReal(v2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DReal(v2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'pebo', v2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
    
!  C部分结果
    call setland_v2(peds)
    call array2DTo1DReal(v2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DReal(v2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),v2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),v2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),v2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'peds', v2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
    
    deallocate(sendBuf)
    if(myid == 0) then
        deallocate(recvBuf)
        deallocate(v2_global)
    end if
  
!  call setland_v2(pein)
!  call writenc(ncid, 'pein', v2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'pein', v2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call setland_v2(pebo)
!  call writenc(ncid, 'pebo', v2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'pebo', v2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call setland_v2(peds)
!  call writenc(ncid, 'peds', v2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'peds', v2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call closenc(ncid)

!  write(6, *)'eng over'

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine outeng_t

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outwav_t

  subroutine outwav_t(filename, rec)

  implicit none
  include 'mpif.h'
  
  character(len=*), intent(in) :: filename
  integer, intent(in) :: rec

!-------------------------------------------------------------------------------

  integer :: ncid, jd
  integer*2, parameter :: ivland = nf_fill_int2

  logical :: ext
  integer :: outrecord, timerec, ittt
  double precision, allocatable :: timealready(:)
  
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  integer,allocatable :: recvBuf(:)
  integer,allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  integer,allocatable :: iv2_global(:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
  
  real(8) :: recvTimeBuf(npe)
  real(8) :: sendTime
!-------------------------------------------------------------------------------

  jd = datenum(1950, 1, 1, 0, 0, 0)

!-------------------------------------------------------------------------------

  if(myid == 0)then

      inquire(file=filename, exist=ext)
  
    !  if(rec == 1)then
      if(.not.ext)then

        outrecord = 1
        call open_nc(ncid, filename, 'c')
        call dimension_define(ncid, 'lon', im, 'lon', nf_real)
        call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
        call dimension_define(ncid, 'time', 0, 'time', nf_double)
        call set_attribute(ncid, 'units', 'degrees_north', 'lat')
        call set_attribute(ncid, 'units', 'degrees_east', 'lon')
        call set_attribute(ncid, 'modulo', '', 'lon')
        call set_attribute(ncid, 'units', 'Days since 1950-01-01 00:00:0.0.', 'time')
        call set_attribute(ncid, 'Start_time', ctime)
    
        call variable_define(ncid, 'windx', nf_int2, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'windy', nf_int2, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'hs',    nf_int2, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'tp',    nf_int2, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'tz',    nf_int2, ['lon', 'lat', 'time'])
        call variable_define(ncid, 'th',    nf_int2, ['lon', 'lat', 'time'])
    
        call set_attribute(ncid, 'missing_value', ivland, 'windx')
        call set_attribute(ncid, 'missing_value', ivland, 'windy')
        call set_attribute(ncid, 'missing_value', ivland, 'hs')
        call set_attribute(ncid, 'missing_value', ivland, 'tp')
        call set_attribute(ncid, 'missing_value', ivland, 'tz')
        call set_attribute(ncid, 'missing_value', ivland, 'th')
    
        call set_attribute(ncid, 'scale_factor', 0.01, 'windx')
        call set_attribute(ncid, 'scale_factor', 0.01, 'windy')
        call set_attribute(ncid, 'scale_factor', 0.01, 'hs')
        call set_attribute(ncid, 'scale_factor', 0.01, 'tp')
        call set_attribute(ncid, 'scale_factor', 0.01, 'tz')
        call set_attribute(ncid, 'scale_factor', 0.1, 'th')
    
        call set_attribute(ncid, 'units', 'm/s', 'windx')
        call set_attribute(ncid, 'units', 'm/s', 'windy')
        call set_attribute(ncid, 'units', 'm'  , 'hs')
        call set_attribute(ncid, 'units', 's'  , 'tp')
        call set_attribute(ncid, 'units', 's'  , 'tz')
        call set_attribute(ncid, 'units', 'deg', 'th')
    
        call set_attribute(ncid, 'longname', 'Zonal Wind Velocity '     , 'windx')
        call set_attribute(ncid, 'longname', 'Meridional Wind Velocity ', 'windy')
        call set_attribute(ncid, 'longname', 'Significant wave height'  , 'hs')
        call set_attribute(ncid, 'longname', 'Mean wave direction'      , 'th')
        call set_attribute(ncid, 'longname', 'Spectrum peak wave period', 'tp')
        call set_attribute(ncid, 'longname', 'Zero-crossing wave period', 'tz')
    
        call end_define(ncid)
      
        call writenc(ncid, 'lon', lon)
        call writenc(ncid, 'lat', lat)
        
        call close_nc(ncid)
    
      else

        call open_nc(ncid, filename, 'w')
        
        timerec = get_dimension_len(ncid, 'time')
        allocate(timealready(timerec))
        call readnc(ncid, 'time', timealready)
        
        outrecord = 1
        do ittt = 1, timerec !-1
            if(dtime-jd > timealready(ittt))then 
              outrecord = ittt + 1
            endif
        enddo
        
        deallocate(timealready)
        
        call close_nc(ncid)
    
      endif
      
    endif ! if(myid == 0)then

  call bcast_to_all(outrecord)

  if(myid==0)write(6, *)'outrecord = ', outrecord

!-------------------------------------------------------------------------------

!分配sendBuf
sendLengthNoOverlay = (ie - is + 1) * (je - js + 1)
sendLength = sendLengthNoOverlay
if(ixoverlay /= 0 .and. loc2(1) == 1)then
    sendLength = sendLength + (ixoverlay + 1)*(je - js + 1)
end if
allocate(sendBuf(sendLength))

totalCounts = 0

!  if(ixoverlay /= 0 .and. loc2(1) == 1)then 
!      call writenc(ncid, 'windx', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
if(myid == 0) then
    do i=1,npe
        length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
        length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
        counts(i) = length_i * length_j
        if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
            counts(i) = counts(i) + (ixoverlay + 1)*length_j
        end if
        totalCounts = totalCounts + counts(i)
        if(i == 1) then
            startPos(i) = 0
        else
            startPos(i) = startPos(i-1) + counts(i-1)
        end if
    end do
    !debug
    !write(*,*) 'total coutns:',totalCounts
    !分配recvBuf
    allocate(recvBuf(totalCounts))
    allocate(iv2_global(im,jm))
end if

    !call opennc(ncid, filename, 'w')
    !call writenc(ncid, 'time', dtime-jd, outrecord)
    !call closenc(ncid)
    sendTime = dtime-jd
    call mpi_gather(sendTime,1,MPI_DOUBLE_PRECISION,recvTimeBuf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(myid == 0) then
        call open_nc(ncid, filename, 'w')
        do i=1,npe
            call writenc(ncid, 'time', recvTimeBuf(i), outrecord)
        end do
        call close_nc(ncid) 
    end if
!各个进程首先重置iv2的值，然后把iv2的数据放到sendBuf丿

!  A部分结果
    call setland_iv2(wx, 0.01)
    call array2DTo1DInt(iv2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DInt(iv2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),iv2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DInt(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),iv2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),iv2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'windx', iv2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
!  B部分结果
    call setland_iv2(wy, 0.01)
    call array2DTo1DInt(iv2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DInt(iv2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),iv2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DInt(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),iv2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),iv2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'windy', iv2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if

!  C部分结果
    call setland_iv2(h1_3, 0.01)
    call array2DTo1DInt(iv2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DInt(iv2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),iv2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DInt(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),iv2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),iv2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'hs', iv2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
! D部分结果    
    call setland_iv2(tpf, 0.01)
    call array2DTo1DInt(iv2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DInt(iv2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),iv2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DInt(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),iv2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),iv2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tp', iv2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if

! E部分结果
    call setland_iv2(ape, 0.01)
    call array2DTo1DInt(iv2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DInt(iv2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),iv2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DInt(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),iv2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),iv2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'tz', iv2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
    
! F部分结果
    call setland_iv2(aet, 0.1)
    call array2DTo1DInt(iv2(is:ie,js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc2(1) == 1)then 
        call array2DTo1DInt(iv2(is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_INTEGER,recvBuf,counts,startPos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),iv2_global(i_start:i_end,j_start:j_end))
                call array1DTo2DInt(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),iv2_global(im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo2DInt(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),iv2_global(i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'th', iv2_global(1:im,1:jm), locs=[1,1], recnum=outrecord)
        call close_nc(ncid)
    end if
    deallocate(sendBuf)
    if(myid == 0) then
        deallocate(recvBuf)
        deallocate(iv2_global)
    end if
    
!  call opennc(ncid, filename, 'w')
!  call writenc(ncid, 'time', dtime-jd, outrecord)

!  call setland_iv2(wx, 0.01)
!  call writenc(ncid, 'windx', iv2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'windx', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

! call setland_iv2(wy, 0.01)
!  call writenc(ncid, 'windy', iv2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'windy', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call setland_iv2(h1_3, 0.01)
!  call writenc(ncid, 'hs', iv2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'hs', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call setland_iv2(tpf, 0.01)
!  call writenc(ncid, 'tp', iv2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'tp', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call setland_iv2(ape, 0.01)
!  call writenc(ncid, 'tz', iv2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'tz', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call setland_iv2(aet, 0.1)
!  call writenc(ncid, 'th', iv2(is:ie, js:je), locs=loc2, recnum=outrecord)
!  if(ixoverlay /= 0 .and. loc2(1) == 1)then
!      call writenc(ncid, 'th', iv2(is:is+ixoverlay, js:je), &
!                   locs = [im-ixoverlay, loc2(2)], recnum=outrecord)
!  endif

!  call closenc(ncid)

!-------------------------------------------------------------------------------

  return


!-------------------------------------------------------------------------------

  end subroutine outwav_t

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: outrst

  subroutine outrst(filename)

  implicit none
   include 'mpif.h'
   
  character(len=*), intent(in) :: filename

  integer :: ncid
  !自定义的变量
  integer :: startPos(npe)
  integer :: counts(npe)
  real(4),allocatable :: recvBuf(:)
  real(4),allocatable :: sendBuf(:)
  integer :: sendLength
  integer :: totalCounts !一共要发送的整数个数
  integer :: i
  integer :: length_i,length_j
  integer :: i_start,i_end,j_start,j_end
  integer :: ierr
  real(4),allocatable :: ee_global(:,:,:,:)
  integer :: sendLengthNoOverlay
  integer :: tmp
!-------------------------------------------------------------------------------
! Modified by Zhenya.Song, 2016/04/03
!  if(mod(it-number, irstfreq) /= 0)return
  if(irstfreq == 0) then
    return
  else
    if(mod(it-number, irstfreq) /= 0)return
  endif

!-------------------------------------------------------------------------------

  if(myid == 0)then
    
      call open_nc(ncid, filename, 'c')
    
      call dimension_define(ncid, 'kk', kl, 'kk', nf_real)
      call dimension_define(ncid, 'jj', jl, 'jj', nf_real)
      call dimension_define(ncid, 'lon', im, 'lon', nf_real)
      call dimension_define(ncid, 'lat', jm, 'lat', nf_real)
    
      call variable_define(ncid, 'ee', nf_real, ['kk', 'jj', 'lon', 'lat'])
    
      call set_attribute(ncid, 'units', 'degrees_north', 'lat')
      call set_attribute(ncid, 'units', 'degrees_east', 'lon')
      call set_attribute(ncid, 'modulo', '', 'lon')
      call set_attribute(ncid, 'ctime', ctime)
    
      call end_define(ncid)

        !  call writenc(ncid, 'kk', x) ! wk(1:kl)
        !  call writenc(ncid, 'jj', x) ! 1:13
    call writenc(ncid, 'lon', lon)
    call writenc(ncid, 'lat', lat)
      call close_nc(ncid)

  endif

  !分配sendBuf
  sendLengthNoOverlay = kl * jl * (ie - is + 1) * (je - js + 1)
  sendLength = sendLengthNoOverlay
  if(ixoverlay /= 0 .and. loc4(3) == 1)then
    sendLength = sendLength + kl * jl * (ixoverlay + 1)*(je - js + 1)
  end if
  allocate(sendBuf(sendLength))

  totalCounts = 0


!0号进程计算startPos以及counts并且分配recvBuf,长度计算信息从pebox里面取得
if(myid == 0) then
    do i=1,npe
        length_i = pebox(3,i) - pebox(2,i) + 1 !ie - is + 1
        length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
        counts(i) = length_i * length_j*kl*jl
        if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
            counts(i) = counts(i) + (ixoverlay + 1)*length_j*kl*jl
        end if
        totalCounts = totalCounts + counts(i)
        if(i == 1) then
            startPos(i) = 0
        else
            startPos(i) = startPos(i-1) + counts(i-1)
        end if
    end do
    !debug
    !write(*,*) 'total coutns:',totalCounts
    !分配recvBuf
    allocate(recvBuf(totalCounts))
end if


    call array4DTo1DReal(ee(:, :, is:ie, js:je),sendBuf(1:sendLengthNoOverlay))
    if(ixoverlay /= 0 .and. loc4(3) == 1)then 
        call array4DTo1DReal(ee(:, :, is:is+ixoverlay, js:je),sendBuf(sendLengthNoOverlay+1:sendLength))
    end if
    !调用gatherV
    call mpi_gatherv(sendBuf,sendLength,MPI_REAL,recvBuf,counts,startPos,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    if(myid == 0) then !0号进程将recvBuf中的数据转换并且写入
        allocate(ee_global(kl,jl,im,jm))
        do i=1,npe
            i_start = pebox(2,i) 
            i_end = pebox(3,i)
            j_start = pebox(4,i)
            j_end = pebox(5,i)
            if(ixoverlay /= 0 .and. pebox(2,i) == 1) then
                length_j = pebox(5,i) - pebox(4,i) + 1 !je - js + 1
                tmp = (ixoverlay + 1)*length_j*kl*jl
                call array1DTo4DReal(recvBuf(startPos(i)+1:startPos(i)+counts(i)-tmp),ee_global(:,:,i_start:i_end,j_start:j_end))
                call array1DTo4DReal(recvBuf(startPos(i)+counts(i)-tmp+1:startPos(i)+counts(i)),ee_global(:,:,im-ixoverlay:im,j_start:j_end))
            else
                call array1DTo4DReal(recvBuf(startPos(i)+1:startPos(i) + counts(i) + 1),ee_global(:,:,i_start:i_end,j_start:j_end))
            end if
        end do
        call open_nc(ncid, filename, 'w')
        call writenc(ncid, 'ee', ee_global(:,:,:,:), locs=[1,1,1,1])
        call close_nc(ncid)
        deallocate(ee_global)
    end if
  !call opennc(ncid, filename, 'w')

  !call writenc(ncid, 'ee', ee(:, :, is:ie, js:je), locs=loc4)
  
! --- How to deal with x=360, 0?

  !if(ixoverlay /= 0 .and. loc4(3) == 1)then
  !    call writenc(ncid, 'ee', ee(:, :, is:is+ixoverlay, js:je), &
  !               locs = [loc4(1), loc4(2), im-ixoverlay, loc4(4)])
  !endif
  
  !call closenc(ncid)
  
  deallocate(sendBuf)
  if(myid == 0) then
    deallocate(recvBuf)
  end if
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

  inquire(file=filename, exist=ext)
  if(ext)then
    key = 1
  else
    key = 0
    return 
  endif

!-------------------------------------------------------------------------------

  call opennc(ncid, filename, 'r')
  call readnc(ncid, 'ee', ee(:, :, is:ie, js:je), locs=loc4)
  call get_attribute(ncid, 'ctime', ctime)
  call closenc(ncid)

  write(6, *)'Restart time is :', ctime
  read(ctime, '(i4.4,5i2.2)')istime
  istime(5:6) = 0

  call updatev(ee, ixl, iyl, kl, jl, halosize)

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine inprst

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setland_iv2

  subroutine setland_iv2(var, scal)

  implicit none

  real, intent(in) :: var(:, :), scal
  integer :: i, j

  if(.not. allocated(iv2))allocate(iv2(ixl, iyl))

  do j = iys, iyl
  do i = ixs, ixl
!    if(nsp(i, j) == 0)then
    if(nsp(i, j) /= 1)then
      iv2(i, j) = nf_fill_int2
    else
      iv2(i, j) = var(i, j) / scal
    endif
  enddo
  enddo

  return

  end subroutine setland_iv2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setland_v2

  subroutine setland_v2(var)

  implicit none

  real, intent(in) :: var(:, :)
  integer :: i, j

  if(.not. allocated(v2))allocate(v2(ixl, iyl))

  do j = iys, iyl
  do i = ixs, ixl
    if(nsp(i, j) == 0)then
      v2(i, j) = nf_fill_real
    else
      v2(i, j) = var(i, j)
    endif
  enddo
  enddo

  return

  end subroutine setland_v2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setland_v3

  subroutine setland_v3(var)

  implicit none

  real, intent(in) :: var(:, :, :)

  integer :: i, j

  if(.not. allocated(v3))allocate(v3(ixl, iyl, kb))
  
  do j = iys, iyl
  do i = ixs, ixl
    if(nsp(i, j) == 0)then
      v3(i, j, :) = nf_fill_real
    else
      v3(i, j, :) = var(i, j, :)
    endif
  enddo
  enddo

  return

  end subroutine setland_v3

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
  !array util中的函数，后面可以独立一个模块出板
  subroutine array2DTo1DInt(array2D,array1D)
        integer :: array2D(:,:)
        integer :: array1D(:)
        integer :: array_size(2),rowNum,colNum
        integer :: curOffset
        integer :: i
        
        array_size = shape(array2D)
        rowNum = array_size(1)
        colNum = array_size(2)
        !debug代码
        !write(*,*) 'row num:',rowNum,'col num:',colNum
        curOffset = 1
        do i=1,colNum
            array1D(curOffset:curOffset+rowNum-1) = array2D(1:rowNum,i)
            curOffset = curOffset + rowNum
        end do
        !debug代码
        !do i=1,rowNum*colNum
        !    write(*,*) 'array1D:',i,':',array1D(i)
        !end do
    end subroutine array2DTo1DInt
    
    subroutine array2DTo1DReal(array2D,array1D)
        real :: array2D(:,:)
        real :: array1D(:)
        integer :: array_size(2),rowNum,colNum
        integer :: curOffset
        integer :: i
        
        array_size = shape(array2D)
        rowNum = array_size(1)
        colNum = array_size(2)
        !debug代码
        !write(*,*) 'row num:',rowNum,'col num:',colNum
        curOffset = 1
        do i=1,colNum
            array1D(curOffset:curOffset+rowNum-1) = array2D(1:rowNum,i)
            curOffset = curOffset + rowNum
        end do
        !debug代码
        !do i=1,rowNum*colNum
        !    write(*,*) 'array1D:',i,':',array1D(i)
        !end do
    end subroutine array2DTo1DReal
    
    !注意array1D和array2D数组要预先定义好大小,且长度要刚刚奿
    subroutine array1DTo2DInt(array1D,array2D)
        integer :: array2D(:,:),array1D(:)
        integer :: array_size(2),rowNum,colNum
        integer :: curOffset
        integer :: i
        
        array_size = shape(array2D)
        rowNum = array_size(1)
        colNum = array_size(2)
        curOffset = 1
        do i=1,colNum
            array2D(1:rowNum,i) = array1D(curOffset:curOffset+rowNum-1)
            curOffset = curOffset + rowNum
        end do
    end subroutine array1DTo2DInt
    
    subroutine array1DTo2DReal(array1D,array2D)
        real :: array2D(:,:),array1D(:)
        integer :: array_size(2),rowNum,colNum
        integer :: curOffset
        integer :: i
        
        array_size = shape(array2D)
        rowNum = array_size(1)
        colNum = array_size(2)
        curOffset = 1
        do i=1,colNum
            array2D(1:rowNum,i) = array1D(curOffset:curOffset+rowNum-1)
            curOffset = curOffset + rowNum
        end do
    end subroutine array1DTo2DReal
    
    
    subroutine printArray1DInt(array1D)
        integer :: array1D(:)
        integer :: length(1)
        integer :: i
        length = shape(array1D)
        
        do i=1,length(1)
            write(*,*) 'array1D  i:',i,array1D(i)
        end do
    end subroutine printArray1DInt
    
    subroutine printArray1DReal(array1D)
        real :: array1D(:)
        integer :: length(1)
        integer :: i
        length = shape(array1D)
        
        do i=1,length(1)
            write(*,*) 'array1D  i:',i,array1D(i)
        end do
    end subroutine printArray1DReal
    
    subroutine printArray2DInt(array2D)
        integer :: array2D(:,:)
        integer :: length(2)
        integer :: i,j
        length = shape(array2D)
        
        do i=1,length(2)
            do j=1,length(1)
                write(*,*) 'array2D  col:',i,'row',j,array2D(j,i)
            end do
        end do
    end subroutine printArray2DInt
    
    subroutine printArray2DReal(array2D)
        real :: array2D(:,:)
        integer :: length(2)
        integer :: i,j
        length = shape(array2D)
        
        do i=1,length(2)
            do j=1,length(1)
                write(*,*) 'array2D  col:',i,'row',j,array2D(j,i)
            end do
        end do
    end subroutine printArray2DReal
    
    subroutine array1DTo4DReal(array1D,array4D)
        real(4) :: array1D(:),array4D(:,:,:,:)
        integer :: array_size(4)
        integer :: curOffset
        integer :: i,j,k
        array_size = shape(array4D)
        curOffset = 1
        do i=1,array_size(4)
            do j=1,array_size(3)
                do k=1,array_size(2)
                    array4D(1:array_size(1),k,j,i) = array1D(curOffset:curOffset+array_size(1)-1)
                    curOffset = curOffset + array_size(1)
                end do
            end do
        end do
    end subroutine array1DTo4DReal
    
    subroutine array4DTo1DReal(array4D,array1D)
        real(4) :: array1D(:),array4D(:,:,:,:)
        integer :: array_size(4)
        integer :: curOffset
        integer :: i,j,k
        array_size = shape(array4D)
        curOffset = 1
        do i=1,array_size(4)
            do j=1,array_size(3)
                do k=1,array_size(2)
                    array1D(curOffset:curOffset+array_size(1)-1) = array4D(1:array_size(1),k,j,i)
                    curOffset = curOffset + array_size(1)
                end do
            end do
        end do
    end subroutine array4DTo1DReal
    
    subroutine array1DTo3DReal(array1D,array3D)
        real(4) :: array1D(:),array3D(:,:,:)
        integer :: array_size(3)
        integer :: curOffset
        integer :: i,j
        array_size = shape(array3D)
        curOffset = 1
        do i=1,array_size(3)
            do j=1,array_size(2)
                array3D(1:array_size(1),j,i) = array1D(curOffset:curOffset+array_size(1)-1)
                curOffset = curOffset + array_size(1)
            end do
        end do
    end subroutine array1DTo3DReal
    
    subroutine array3DTo1DReal(array3D,array1D)
        real(4) :: array1D(:),array3D(:,:,:)
        integer :: array_size(3)
        integer :: curOffset
        integer :: i,j,k
        array_size = shape(array3D)
        curOffset = 1
        do i=1,array_size(3)
            do j=1,array_size(2)
                array1D(curOffset:curOffset+array_size(1)-1) = array3D(1:array_size(1),j,i)
                curOffset = curOffset + array_size(1)
            end do
        end do
    end subroutine array3DTo1DReal
!-----------------------------------------
  
  
  end module wammpi_mod

!-------------------------------------------------------------------------------
!###############################################################################
