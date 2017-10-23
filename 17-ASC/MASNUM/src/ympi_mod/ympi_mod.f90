!###############################################################################
!-------------------------------------------------------------------------------

  module ympi_mod

!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------
!                                                Copyright (C) 2005 Xunqiang Yin
!                                                MODULE NAME : ympi_mod
!                                                Current VERSION : 2011/10/19
!
! --- USAGE : Tools for parallel computation using MPI.
! --- DEPEND: netcdf_mod developed by Xunqiang Yin.
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
!  1. subroutine ympi_init
!     --- Initialize this module.
!
!-------------------------------------------------------------------------------
!
!  2. subroutine ympi_final
!     --- Finalize this module.
!
!-------------------------------------------------------------------------------
!
!  3. subroutine init_pebox(filename, xname, yname, mname, n, myid)
!     * character(len=*) :: filename = the name of grid file.
!     * character(len=*) :: xname = the name of x-direction.
!     * character(len=*) :: yname = the name of y-direction.
!     * character(len=*) :: mname = the name of mask arrary.
!     * integer          :: n = number of all PEs.
!     * integer          :: myid = index of the current PE.
!
!-------------------------------------------------------------------------------
!
!  4. subroutine opennc(ncid, filename, action)
!     # integer :: ncid = unit number of opened netcdf file.
!     * character(len=*) :: filename = the name of netcdf file.
!     * character :: action = 'create', 'define', 'write' or 'read'.
!
!-------------------------------------------------------------------------------
!
!  5. subroutine closenc(ncid)
!     * integer :: ncid = unit number of opened netcdf file.
!
!-------------------------------------------------------------------------------
!
!  6. subroutine get_mpimin(x)
!     @ integer/real/double precision :: x = the value need to get min.
!
!-------------------------------------------------------------------------------
!
!  7. subroutine bcast_to_all(x)
!     @ integer/real/double precision :: x = the value need to be broadcasted.
!
!-------------------------------------------------------------------------------
!
!  8. subroutine updatev(ee, ixl, iyl[, kl[, jl]], halo[, flag])
!     @ real(4/8) :: ee = the array need to be updated, could be 2d,3d & 4d.
!     * integer :: ixl, iyl, kl, jl = is for the shape of ee.
!     * integer :: halo = the size of halo in spatial partition.
!     * flag = the flag for 4d direction of 4d array.
!       --- If it is 4d & flag is given, the direction of dimensions is 
!           (ixl, iyl, kl, jl), otherwise, it is (kl, jl, ixl, iyl)
!       --- For 3d array, it is defined as ee(ixl, iyl, kl)
!
!-------------------------------------------------------------------------------
!
! --- Useful variables:
!
!     lon, lat, im, jm, halosize, flagxcycle, ixoverlay
!     myid, npe, mypebox, loc2, loc3, loc4, is, ie, js, je
!
!-------------------------------------------------------------------------------
!
!                                                --- Xunqiang Yin, 2011/10/19
!                                                 E-Mail: XunqiangYin@gmail.com
!
!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------
  
  use wamvar_mod!, only : pebox, peboxex, my_hor_start, my_hor_end, my_hor_len, my_ver_start, my_ver_end, my_ver_len, ver_inner_start, ver_inner_end, hor_inner_start, hor_inner_end

  use netcdf_mod
  
  implicit none

!-------------------------------------------------------------------------------

  public ympi_init, ympi_final, init_pebox
  public get_mpimin, bcast_to_all, updatev, opennc, closenc, runbyturn
  public lon, lat, im, jm, halosize, flagxcycle, ixoverlay
  public myid, npe, mypebox, loc2, loc3, loc4, is, ie, js, je
  !private

!-------------------------------------------------------------------------------

  include 'mpif.h'
  integer :: myid, npe, mpi_comm_ympi
  integer :: loc2(2), loc3(3), loc4(4), is, ie, js, je
  integer, allocatable :: qrs(:), qrr(:), qls(:), qlr(:), tgl(:), tgr(:)
  integer, allocatable :: statr(:, :), statl(:, :)
  integer, allocatable :: stat_r_send(:, :), stat_r_recv(:, :)
  integer, allocatable :: stat_l_send(:, :), stat_l_recv(:, :)
  integer :: im, jm, flagxcycle, halosize, ixoverlay !,ia,ic,j,k
  real(4), allocatable :: lon(:), lat(:), mask(:, :)
  real(4), allocatable :: dnee1_3d(:,:,:), dnee_3d(:,:,:)
  real(4), allocatable :: upee1_3d(:,:,:), upee_3d(:,:,:)
  real(4), allocatable :: lfee1_3d(:,:,:), lfee_3d(:,:,:)
  real(4), allocatable :: rtee1_3d(:,:,:), rtee_3d(:,:,:)
  integer :: updn_req(4)
  integer, allocatable  :: lsreq(:), lrreq(:)
  integer, allocatable  :: rsreq(:), rrreq(:)
  !real(4), allocatable :: bsend_buffer(:)



!-------------------------------------------------------------------------------

  type pebox_type
    integer :: myid, i1, i2, j1, j2   ! 1,2,3,4,5
    integer :: up, dn, nr, nl, np     ! 7, 6, 8, 9, 10
    integer :: ei1, ei2, ej1, ej2     ! 11, 12, 13, 14
    integer :: halo, isize, jsize     ! 15, 16
    integer, allocatable :: left(:, :)    ! (nl, 5): id, rj1, rj2, sj1, sj2
    integer, allocatable :: right(:, :)   ! (nr, 5): id, rj1, rj2, sj1, sj2
    integer :: nr_send, nr_recv, nl_send, nl_recv
    integer, allocatable :: left_send(:, :)  !(nl, 3): id, sj1, sj2
    integer, allocatable :: right_send(:, :) !(nr, 3): id, sj1, sj2
    integer, allocatable :: left_recv(:, :)  !(nl, 3): id, rj1, rj2
    integer, allocatable :: right_recv(:, :) !(nr, 3): id, rj1, rj2
    integer :: hor_id, ver_id
  end type pebox_type
  type(pebox_type) :: mypebox
  
!-------------------------------------------------------------------------------

  interface updatev
    module procedure update2d_sgl, update3d_sgl, update4d_sgl, update4d_sgle, &
                     update2d_dbl, update3d_dbl, update4d_dbl, update4d_dble
  end interface updatev

  interface get_mpimin
    module procedure get_mpimini, get_mpiminr, get_mpimind
  end interface get_mpimin

  interface bcast_to_all
    module procedure bcast_to_alli, bcast_to_allr, bcast_to_alld
  end interface bcast_to_all
  
!-------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: runbyturn

  subroutine runbyturn(flag)
    integer, intent(in) :: flag
    integer :: ierr, i = 0
    integer :: stat(mpi_status_size)
    if(flag == 0 .and. myid /= 0)then
      ! recv message from myid-1, tag=229
      call mpi_recv(i, 1, mpi_integer, myid-1, 229, mpi_comm_ympi, stat, ierr)
    endif
    if(flag == 1 .and. myid < npe-1)then
      ! send message to myid+1, tag=229
      call mpi_send(i, 1, mpi_integer, myid+1, 229, mpi_comm_ympi, ierr)
    endif
  end subroutine runbyturn

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: opennc

  subroutine opennc(ncid, filename, action)
    character(len=*), intent(in) :: filename, action
    integer, intent(in) :: ncid
    integer :: flag, ierr
    integer :: stat(mpi_status_size)
    if(myid /= 0)then
      ! recv message from myid-1, tag=229
      call mpi_recv(flag, 1, mpi_integer, myid-1, 229, mpi_comm_ympi, stat, ierr)
    endif
    call open_nc(ncid, filename, action)
  end subroutine opennc

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: closenc

  subroutine closenc(ncid)
    integer, intent(in) :: ncid
    integer :: flag=0, ierr
    integer :: stat(mpi_status_size)
    call close_nc(ncid)
    if(myid < npe-1)then
      ! send message to myid+1, tag=229
      call mpi_send(flag, 1, mpi_integer, myid+1, 229, mpi_comm_ympi, ierr           )
    endif
  end subroutine closenc

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: ympi_init

  subroutine ympi_init
    integer :: ierr
    character(len=10) :: str
    ! --- Initial mpi, get myid & npe
    call mpi_init(ierr)
    mpi_comm_ympi = mpi_comm_world
    call mpi_comm_rank(mpi_comm_ympi, myid, ierr)
    call mpi_comm_size(mpi_comm_ympi, npe, ierr)
    !print *, myid;
  end subroutine ympi_init

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: ympi_final

  subroutine ympi_final
    integer :: rc
    close(6)
    call mpi_finalize(rc); stop
  end subroutine ympi_final

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_mpimini

  subroutine get_mpimini(x)
    integer, intent(inout) :: x
    integer :: y
    integer :: ierr
    y = x
    !call mpi_reduce(y, x, 1, mpi_integer, mpi_min, 0, mpi_comm_ympi, ierr)
    !call mpi_bcast(x, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
  call MPI_ALLREDUCE( y, x, 1, mpi_integer,mpi_min,mpi_comm_ympi, ierr)
  end subroutine get_mpimini
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_mpiminr

  subroutine get_mpiminr(x)
    real(4), intent(inout) :: x
    real(4) :: y
    integer :: ierr
    y = x
    !call mpi_reduce(y, x, 1, mpi_real, mpi_min, 0, mpi_comm_ympi, ierr)
    !call mpi_bcast(x, 1, mpi_real, 0, mpi_comm_ympi, ierr)
  call MPI_ALLREDUCE( y, x, 1, mpi_real,mpi_min,mpi_comm_ympi, ierr)
  end subroutine get_mpiminr

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_mpimind

  subroutine get_mpimind(x)
    real(8), intent(inout) :: x
    real(8) :: y
    integer :: ierr
    y = x
    !call mpi_reduce(y, x, 1, mpi_double_precision, mpi_min, 0, mpi_comm_ympi, ierr)
    !call mpi_bcast(x, 1, mpi_double_precision, 0, mpi_comm_ympi, ierr)
  call MPI_ALLREDUCE( y, x, 1, mpi_double_precision,mpi_min,mpi_comm_ympi, ierr)
  end subroutine get_mpimind
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: bcast_to_alli

  subroutine bcast_to_alli(ii)
    integer, intent(inout) :: ii
    integer :: ierr
    call mpi_bcast(ii, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
  end subroutine bcast_to_alli
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: bcast_to_allr

  subroutine bcast_to_allr(ii)
    real(4), intent(inout) :: ii
    integer :: ierr
    call mpi_bcast(ii, 1, mpi_real, 0, mpi_comm_ympi, ierr)
  end subroutine bcast_to_allr
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: bcast_to_alld

  subroutine bcast_to_alld(ii)
    real(8), intent(inout) :: ii
    integer :: ierr
    call mpi_bcast(ii, 1, mpi_double_precision, 0, mpi_comm_ympi, ierr)
  end subroutine bcast_to_alld
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update2d_sgl

  include 'mpi_update_mod.inc'
  subroutine update2d_sgl(ee, ixl, iyl, halo)
    integer, intent(in) :: ixl, iyl, halo
    real(4), intent(inout) :: ee(ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(4), dimension(ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    call mpi_barrier(mpi_comm_ympi, ierr)
    dnee1(:, :) = ee(:, halosize+1:halosize*2)
    upee1(:, :) = ee(:, iyl-halosize*2+1:iyl-halosize)
    ncount = ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr)
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr)
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, iyl-halosize+1:iyl) = upee
    lfee1(:, :) = ee(halosize+1:halosize*2, :)
    rtee1(:, :) = ee(ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, mypebox%left(i, 4)), ncount, mpi_real, mypebox%left(i, 1), 101, mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, mypebox%left(i, 4)), ncount, mpi_real, mypebox%left(i, 1), 201, mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr, 1 !myrightno
      ! --- Send to right.
      ncount = halosize*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, mypebox%right(i, 4)), ncount, mpi_real, mypebox%right(i, 1), 201, mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv( rtee(1, mypebox%right(i, 4)), ncount, mpi_real, mypebox%right(i, 1), 101, mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(1:halosize, :) = lfee
    endif
  end subroutine update2d_sgl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update3d_sgl

  subroutine update3d_sgl(ee, ixl, iyl, kl, halo)
    integer, intent(in) :: ixl, iyl, kl, halo
    real(4), intent(inout) :: ee(ixl, iyl, kl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(4), dimension(kl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(kl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do i = 1, ixl
      dnee1(k, i, :) = ee(i, halosize+1:halosize*2, k)
      upee1(k, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k)
    enddo
    enddo
    ncount = kl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
      enddo
      enddo
    endif
    do k = 1, kl
    do i = 1, iyl
      lfee1(k, :, i) = ee(halosize+1:halosize*2        , i, k)
      rtee1(k, :, i) = ee(ixl-halosize*2+1:ixl-halosize, i, k)
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k) = rtee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(1:halosize, i, k) = lfee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
    endif
  end subroutine update3d_sgl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_sgle

  subroutine update4d_sgle(ee, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(4), intent(inout) :: ee(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(4), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    dnee1(:, :, :, :) = ee(:, :, :, halosize+1:halosize*2)
    upee1(:, :, :, :) = ee(:, :, :, iyl-halosize*2+1:iyl-halosize)
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, :, :, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, :, :, iyl-halosize+1:iyl) = upee
    lfee1(:, :, :, :) = ee(:, :, halosize+1:halosize*2, :)
    rtee1(:, :, :, :) = ee(:, :, ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(:, :, ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(:, :, 1:halosize, :) = lfee
    endif
  end subroutine update4d_sgle

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_sgl

  subroutine update4d_sgl(ee, ixl, iyl, kl, jl, halo, flag)
    integer, intent(in) :: ixl, iyl, kl, jl, halo, flag
    real(4), intent(inout) :: ee(ixl, iyl, kl, jl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(4), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, ixl
      dnee1(k, j, i, :) = ee(i, halosize+1:halosize*2, k, j)
      upee1(k, j, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k, j)
    enddo
    enddo
    enddo
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
      enddo
      enddo
      enddo
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, iyl
      lfee1(k, j, :, i) = ee(halosize+1:halosize*2        , i, k, j)
      rtee1(k, j, :, i) = ee(iyl-halosize*2+1:iyl-halosize, i, k, j)
    enddo
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k, j) = rtee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(1:halosize, i, k, j) = lfee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
  end subroutine update4d_sgl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update2d_dbl

  subroutine update2d_dbl(ee, ixl, iyl, halo)
    integer, intent(in) :: ixl, iyl, halo
    real(8), intent(inout) :: ee(ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(8), dimension(ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    call mpi_barrier(mpi_comm_ympi, ierr)
    dnee1(:, :) = ee(:, halosize+1:halosize*2)
    upee1(:, :) = ee(:, iyl-halosize*2+1:iyl-halosize)
    ncount = ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr)
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr)
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, iyl-halosize+1:iyl) = upee
    lfee1(:, :) = ee(halosize+1:halosize*2, :)
    rtee1(:, :) = ee(ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, mypebox%left(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101, mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, mypebox%left(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201, mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr, 1 !myrightno
      ! --- Send to right.
      ncount = halosize*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, mypebox%right(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201, mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv( rtee(1, mypebox%right(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101, mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(1:halosize, :) = lfee
    endif
  end subroutine update2d_dbl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update3d_dbl

  subroutine update3d_dbl(ee, ixl, iyl, kl, halo)
    integer, intent(in) :: ixl, iyl, kl, halo
    real(8), intent(inout) :: ee(ixl, iyl, kl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(8), dimension(kl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(kl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do i = 1, ixl
      dnee1(k, i, :) = ee(i, halosize+1:halosize*2, k)
      upee1(k, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k)
    enddo
    enddo
    ncount = kl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
      enddo
      enddo
    endif
    do k = 1, kl
    do i = 1, iyl
      lfee1(k, :, i) = ee(halosize+1:halosize*2        , i, k)
      rtee1(k, :, i) = ee(ixl-halosize*2+1:ixl-halosize, i, k)
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k) = rtee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(1:halosize, i, k) = lfee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
    endif
  end subroutine update3d_dbl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_dble

  subroutine update4d_dble(ee, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(8), intent(inout) :: ee(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(8), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    dnee1(:, :, :, :) = ee(:, :, :, halosize+1:halosize*2)
    upee1(:, :, :, :) = ee(:, :, :, iyl-halosize*2+1:iyl-halosize)
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, :, :, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, :, :, iyl-halosize+1:iyl) = upee
    lfee1(:, :, :, :) = ee(:, :, halosize+1:halosize*2, :)
    rtee1(:, :, :, :) = ee(:, :, ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(:, :, ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
     call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(:, :, 1:halosize, :) = lfee
    endif
  end subroutine update4d_dble

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_dbl

  subroutine update4d_dbl(ee, ixl, iyl, kl, jl, halo, flag)
    integer, intent(in) :: ixl, iyl, kl, jl, halo, flag
    real(8), intent(inout) :: ee(ixl, iyl, kl, jl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(8), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, ixl
      dnee1(k, j, i, :) = ee(i, halosize+1:halosize*2, k, j)
      upee1(k, j, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k, j)
    enddo
    enddo
    enddo
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
      enddo
      enddo
      enddo
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, iyl
      lfee1(k, j, :, i) = ee(halosize+1:halosize*2        , i, k, j)
      rtee1(k, j, :, i) = ee(iyl-halosize*2+1:iyl-halosize, i, k, j)
    enddo
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k, j) = rtee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(1:halosize, i, k, j) = lfee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
  end subroutine update4d_dbl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_pebox !input parameters is needed before this sub.

  subroutine init_pebox(filename, xname, yname, mname, n, myid)
    integer, intent(in) :: n, myid
    character(len=*), intent(in) :: filename, xname, yname, mname
    integer :: ncid, i, ierr, im1, x, j
    !integer, allocatable :: pebox(:, :), peleft(:, :, :), peright(:, :, :)
    integer, allocatable :: peleft(:, :, :), peright(:, :, :)
    integer, allocatable :: peleftsend(:, :, :), perightsend(:, :, :)
    integer, allocatable :: peleftrecv(:, :, :), perightrecv(:, :, :)
    integer :: ncount

    allocate(pebox(17, n), peleft(n, n, 5), peright(n, n, 5))
    allocate(peleftsend(n, n, 3), perightsend(n, n, 3))
    allocate(peleftrecv(n, n, 3), perightrecv(n, n, 3))
    allocate(peboxex(7, n))

    if(myid == 0)then
      ! --- Input lon, lat & mask.
      call open_nc(ncid, filename, 'r')
      im = get_dimension_len(ncid, xname)
      jm = get_dimension_len(ncid, yname)
      allocate(lon(im), lat(jm), mask(im, jm))
      call readnc(ncid, xname, lon)
      call readnc(ncid, yname, lat)
      call readnc(ncid, mname, mask)
      call close_nc(ncid)
      ! --- Check xcycle & ixoverlay for global model.
      flagxcycle = 1 ! --- For regional model.
      !if((lon(1) + 360. - lon(im)) < (lon(2) - lon(1)))flagxcycle = 0
      if((2*lon(1) + 360. - lon(im) - lon(2)) < 0)flagxcycle = 0 ! --- global
      ixoverlay = 0
      if(flagxcycle == 0)then
        do i = 1, im
          if(lon(i) > lon(im) - 360)then
            exit
          else
            ixoverlay = i
          endif
        enddo
      endif
      im1 = im - ixoverlay
      ! --- Set pebox for all processors.
      call set_pebox(im1, jm, n, flagxcycle, halosize, lon(1:im-1), lat, mask(1:im1, :), pebox, peleft, peright, peleftsend, perightsend, peleftrecv, perightrecv, peboxex)
    endif
  ! --- Bcast pebox, peleft, peright & some parameters.
    call mpi_bcast(pebox  , 17*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peleft , 5*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peright, 5*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(im, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(jm, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(flagxcycle, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(ixoverlay, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peleftsend, 3*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peleftrecv, 3*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(perightsend, 3*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(perightrecv, 3*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peboxex  , 7*n, mpi_integer, 0, mpi_comm_ympi, ierr)

    mypebox%myid  = pebox(1, myid+1)
    mypebox%i1    = pebox(2, myid+1)
    mypebox%i2    = pebox(3, myid+1)
    mypebox%j1    = pebox(4, myid+1)
    mypebox%j2    = pebox(5, myid+1)
    mypebox%dn    = pebox(6, myid+1)
    mypebox%up    = pebox(7, myid+1)
    mypebox%nr    = pebox(8, myid+1)
    mypebox%nl    = pebox(9, myid+1)
    mypebox%np    = pebox(10, myid+1)
    mypebox%ei1   = pebox(11, myid+1)
    mypebox%ei2   = pebox(12, myid+1)
    mypebox%ej1   = pebox(13, myid+1)
    mypebox%ej2   = pebox(14, myid+1)
    mypebox%halo  = pebox(15, myid+1)
    mypebox%isize = pebox(16, myid+1)
    mypebox%jsize = pebox(17, myid+1)

    mypebox%nr_send = peboxex(4, myid+1)
    mypebox%nr_recv = peboxex(5, myid+1)
    mypebox%nl_send = peboxex(2, myid+1)
    mypebox%nl_recv = peboxex(3, myid+1)
    !mypebox%hor_id  = peboxex(6, myid+1)
    !mypebox%ver_id  = peboxex(7, myid+1)

    allocate(mypebox%left(mypebox%nl, 5), mypebox%right(mypebox%nr, 5))
    mypebox%left(:, :) = peleft(myid+1, 1:mypebox%nl, 1:5)
    mypebox%right(:, :) = peright(myid+1, 1:mypebox%nr, 1:5)

    allocate(mypebox%left_send(mypebox%nl_send, 3), mypebox%left_recv(mypebox%nl_recv, 3))
    allocate(mypebox%right_send(mypebox%nr_send, 3), mypebox%right_recv(mypebox%nr_recv, 3))
    mypebox%left_send(:, :) = peleftsend(myid+1, 1:mypebox%nl_send, 1:3)
    mypebox%left_recv(:, :) = peleftrecv(myid+1, 1:mypebox%nl_recv, 1:3)
    mypebox%right_send(:, :) = perightsend(myid+1, 1:mypebox%nr_send, 1:3)
    mypebox%right_recv(:, :) = perightrecv(myid+1, 1:mypebox%nr_recv, 1:3)

    my_hor_start = 1
    hor_inner_start = 1
    if(mypebox%nl > 0)then
      my_hor_start = my_hor_start + 1
      hor_inner_start = my_hor_start + 1
    endif

    my_hor_end = mypebox%isize
    hor_inner_end = mypebox%isize
    if(mypebox%nr > 0)then
      my_hor_end = my_hor_end - 1
      hor_inner_end = my_hor_end - 1
    endif

    my_ver_start = 1
    ver_inner_start = 1
    my_ver_end = mypebox%jsize
    ver_inner_end = mypebox%jsize

    if(mypebox%up >= 0)then
      my_ver_end = my_ver_end - 1
      ver_inner_end = my_ver_end - 1
    endif

    if(mypebox%dn >= 0)then
      my_ver_start = my_ver_start + 1
      ver_inner_start = my_ver_start + 1
    endif

    my_ver_len = my_ver_end - my_ver_start + 1
    my_hor_len = my_hor_end - my_hor_start + 1

    do i = 1, n
      if(i == myid + 1)then
        print *, "myid = ", myid, mypebox%i1, mypebox%i2, mypebox%j1, mypebox%j2!, mypebox%hor_id, mypebox%ver_id
        print *, mypebox%nl_send, mypebox%nl_recv, mypebox%nr_send, mypebox%nr_recv, mypebox%dn, mypebox%up
        print *, my_ver_start, my_ver_end, my_hor_start, my_hor_end 
        print *, ver_inner_start, ver_inner_end, hor_inner_start, hor_inner_end 
        !do j = 1, mypebox%nl_send 
        !  print *, "left_send", j, mypebox%left_send(j, :)
        !end do 
        !do j = 1, mypebox%nl_recv 
        !  print *, "left_recv", j, mypebox%left_recv(j, :)
        !end do 
        !do j = 1, mypebox%nr_send 
        !  print *, "right_send", j, mypebox%right_send(j, :)
        !end do 
        !do j = 1, mypebox%nr_recv 
        !  print *, "right_recv", j, mypebox%right_recv(j, :)
        !end do 
        !do j = 1, mypebox%nl
        !    print *, "nl        ", j, mypebox%left(j, :)
        !end do 
        !do j = 1, mypebox%nr
        !    print *, "nr        ", j, mypebox%right(j, :)
        !end do 
      endif
      call mpi_barrier(mpi_comm_ympi, ierr)
    end do

    loc2 = [mypebox%i1, mypebox%j1]
    loc3 = [mypebox%i1, mypebox%j1, 1]
    loc4 = [1, 1, mypebox%i1, mypebox%j1]
    is = mypebox%i1 - mypebox%ei1 + 1
    ie = mypebox%i2 - mypebox%ei1 + 1
    js = mypebox%j1 - mypebox%ej1 + 1
    je = mypebox%j2 - mypebox%ej1 + 1
    if(mypebox%nr>0)then
      allocate(qrs(mypebox%nr), qrr(mypebox%nr), tgr(mypebox%nr), statr(mpi_status_size, mypebox%nr))
      qrs = 0; qrr = 0;
    endif
    if(mypebox%nl>0)then
      allocate(qls(mypebox%nl), qlr(mypebox%nl), tgl(mypebox%nl), statl(mpi_status_size, mypebox%nl))
      qlr = 0; qls = 0
    endif
    !do x = 0, n

    if(mypebox%nr_send > 0)then
      allocate(stat_r_send(mpi_status_size, mypebox%nr_send))
    endif
    if(mypebox%nr_recv > 0)then
      allocate(stat_r_recv(mpi_status_size, mypebox%nr_recv))
    endif
    if(mypebox%nl_send > 0)then
      allocate(stat_l_send(mpi_status_size, mypebox%nl_send))
    endif
    if(mypebox%nl_recv > 0)then
      allocate(stat_l_recv(mpi_status_size, mypebox%nl_recv))
    endif

    if(myid == 0)then
      write(6, *)'myid  = ', mypebox%myid      
      write(6, *)'i1    = ', mypebox%i1        
      write(6, *)'i2    = ', mypebox%i2        
      write(6, *)'j1    = ', mypebox%j1        
      write(6, *)'j2    = ', mypebox%j2        
      write(6, *)'dn    = ', mypebox%dn        
      write(6, *)'up    = ', mypebox%up        
      write(6, *)'nr    = ', mypebox%nr        
      write(6, *)'nl    = ', mypebox%nl        
      write(6, *)'np    = ', mypebox%np        
      write(6, *)'ei1   = ', mypebox%ei1       
      write(6, *)'ei2   = ', mypebox%ei2       
      write(6, *)'ej1   = ', mypebox%ej1       
      write(6, *)'ej2   = ', mypebox%ej2       
      write(6, *)'halo  = ', mypebox%halo 
      write(6, *)'isize = ', mypebox%isize     
      write(6, *)'jsize = ', mypebox%jsize     
      write(6, "('left  (id, rj1, rj2, sj1, sj2) : ', 5i5)")mypebox%left(:, :)
      write(6, "('right (id, rj1, rj2, sj1, sj2) : ', 5i5)")mypebox%right(:, :)
      print *, my_hor_start, my_hor_end, my_hor_len, hor_inner_start, hor_inner_end 
      print *, my_ver_start, my_ver_end, my_ver_len, ver_inner_start, ver_inner_end 
    endif
    write(6, *) mypebox%myid,mypebox%i1,mypebox%i2,mypebox%j1,mypebox%j2,mypebox%ei1,mypebox%ei2,mypebox%ej1,mypebox%ej2
    !enddo
    !deallocate(pebox, peleft, peright)

!  do i = 1, mypebox%nl_send 
!      ncount = 25*12*(mypebox%left_send(i, 3)-mypebox%left_send(i, 2)+1)
!      if(ncount < 0 .or. mypebox%left_send(i, 1) < 0 .or. mypebox%left_send(i, 1) > 63)then
!        print * , "isend, myid = ", myid, i, "left_send: ", mypebox%left_send(i, :)
!      endif
!  end do

  !check avail

  end subroutine init_pebox

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_pebox

  subroutine set_pebox(im, jm, n, flagxcycle, halosize, lon, lat, mask, pebox, peleft, peright, peleftsend, perightsend, peleftrecv, perightrecv, peboxex)
    !implicit none
    integer, intent(in) :: im, jm, n, flagxcycle, halosize
    real(4), intent(in) :: lon(im), lat(jm), mask(im, jm)
    integer, intent(out) :: pebox(17, n), peleft(n, n, 5), peright(n, n, 5)
    integer, intent(out) :: peboxex(7, n), peleftsend(n, n, 3), perightsend(n, n, 3), peleftrecv(n, n, 3), perightrecv(n, n, 3)
    real(4) :: avepoint, ttpoint
    integer, allocatable :: mask_water_only(:,:)
    integer :: i1, i2, i, n1, n2, ii, j1, j2, j, jj, nn,iii,jjj
    integer :: nr, nl, nx, ny, nxmore
    integer :: i_send_start, i_send_end, i_recv_start, i_recv_end
    integer :: nn_send_start, nn_send_end, nn_recv_start, nn_recv_end
    integer :: have_send_flag, have_recv_flag
    integer :: nr_send, nr_recv, nl_send, nl_recv 

    !modified by Huang Chenghuan, 20170418
    integer, allocatable :: mask_walked(:)
    integer :: nnn, sum_col, sum_dst, sum_last, sum_cur, cur_row, cur_col, diff_cur, diff_last, cur_id_in_col, row_start
    integer :: ver_id, hor_id
    integer :: total_max

    !initialize the mask_water_only
    allocate(mask_water_only(im,jm))
    do jjj=1,jm
        do iii=1,im
                if(mask(iii,jjj).ne.1) then
                        mask_water_only(iii,jjj) = 0
                else
                        mask_water_only(iii,jjj) = 1
                endif
        enddo
    enddo

    !modified by Huang Chenghuan, 20170418
    allocate(mask_walked(n))
    mask_walked = 0
  
   ! write(6, *)'x is not ok'
!    ! --- Along x-direction.
    ttpoint = sum(mask_water_only)
    avepoint = sum(mask_water_only) / float(n)

    !Huang Chenghuan 2018-03-31 21:04:24
    !妈蛋国桥你写的我看不懂啊，这缩进也没法看
    nx = sqrt(n * im / float(jm))
    if(nx < 1) then
      nx = 1
    endif
    ny = n / nx
    nxmore = n - nx * ny
    print *, "nx, ny, nxmore ", nx, ny, nxmore

    !当前遍历到的进程列
    cur_col = 1
    !当前进程列的起始格点列
    jj = 1
    !未结算的所有列的进程数量
    iii = n
    !当前进程列的起点进程号
    jjj = 1
    !未结算的所有列的海洋格点数量之和
    sum_col = ttpoint

    do j = 1, nx
      !这一进程列的进程数量
      ii = ny
      if(j <= nxmore)then
        !假如nxmore = 1，那么只有j=1有padding
        ii = ii + 1
      endif

      sum_cur = 0
      sum_last = 0
      sum_dst = (sum_col * ii) / iii

      print *, "nx, expect ", nx, ii, sum_dst, sum_col
      
      !如果只剩下最后一进程列，那么就不用算了
      if(j == nx)then
        print *, "nx, range ", nx, jjj, jjj + ii - 1, jj, im, sum(mask_water_only(jj : im, :))
        do nn = jjj, jjj + ii - 1
          pebox(1 : 3, nn) = [nn - 1, jj, im]
        end do
        exit
      endif

      !开始遍历
      do cur_row = jj, im
        sum_cur = sum_cur + sum(mask_water_only(cur_row, :))
        if(sum_cur >= sum_dst)then
          exit
        endif
        sum_last = sum_cur
      end do

      !结束遍历，开始确定要选那边
      diff_cur = sum_cur - sum_dst
      diff_last = sum_dst - sum_last

      !判断是否需要回退一步
      if(diff_cur > diff_last)then
        cur_row = cur_row - 1
      endif

      !设定pebox
      print *, "nx, range", nx, jjj, jjj + ii - 1, jj, cur_row, sum(mask_water_only(jj : cur_row, :))
      do nn = jjj, jjj + ii - 1
        pebox(1 : 3, nn) = [nn - 1, jj, cur_row]
      end do

      sum_col = sum_col - sum(mask_water_only(jj : cur_row, :))
      
      jj = cur_row + 1
      jjj = jjj + ii
      iii = iii - ii
    end do

      !开始寻找这一列的终点x

    !Do some bug fix on the code of Guoqiao Ye.
    !Huang Chenghuan, 20170418
    !on the vertical direction, in one column

    print *, "begin divide"
    total_max = 0
    hor_id = 1
    do nn = 1, n
        if(mask_walked(nn) /= 0) cycle
        !now this is the bottom of this column
        !if this is the last proc, then do nothing!
!         if(nn == n) exit

        !list all procs in this column
        do nnn = nn + 1, n
            if(pebox(2, nn) /= pebox(2, nnn)) then
                !nnn = nnn - 1
                exit
            endif
            !if(nnn == n) exit !to avoid modification of nnn to n + 1
        end do 
        nnn = nnn - 1
        !now, nnn is the top proc in this column

        !tag all the procs as walked.
        !并不是一个个进程看，而是一列列看
        mask_walked(nn : nnn) = 1

!         在删除了叶国桥的代码后这段不算数
!         if(nnn == nn) cycle !no need to do anything!
        if(nnn == nn)then
          pebox(4, nn) = 1
          pebox(5, nn) = jm
!           print *, "only one", pebox(2, nn), pebox(3, nn)
          cycle
        endif

        !print *, "enter stage 2, nn =", nn, " nnn =", nnn

        !如果当前列只有一个进程，那么就不用分了

        sum_col = sum(mask_water_only(pebox(2, nn) : pebox(3, nn), :))
        sum_dst = sum_col / (nnn - nn + 1)
        sum_cur = 0
        sum_last = 0
        row_start = 1

!         print *, "in this column, the total marine point is ", sum_col
!         print *, "now, dst is ", sum_dst, ", col:", pebox(2, nn), pebox(3, nn)

        cur_id_in_col = 0
        !do cur_row = 1, jm
        cur_row = 1
        do while(.true.)
            if(cur_row == jm + 1)exit
            sum_cur = sum_cur + sum(mask_water_only(pebox(2, nn) : pebox(3, nn), cur_row))
            if(sum_cur >= sum_dst)then

                !set a proc
                diff_cur = sum_cur - sum_dst
                diff_last = sum_dst - sum_last
                print *, "sum_cur >= sum_dst, cur =", sum_cur, "last =", sum_last, "diff = ", diff_cur, diff_last, cur_row
                

                if(diff_cur < diff_last)then
                    pebox(4, nn + cur_id_in_col) = row_start
                    pebox(5, nn + cur_id_in_col) = cur_row
                    sum_col = sum_col - sum_cur
!                     total_max = max(total_max, sum_cur)
!                     print *, sum_cur
                else
                    pebox(4, nn + cur_id_in_col) = row_start
                    pebox(5, nn + cur_id_in_col) = cur_row - 1
                    sum_col = sum_col - sum_last

                    !黄承欢于2018-03-31 22:33:55
                    !非常严重的bug，如果没有加入这一行
                     cur_row = cur_row - 1
!                     total_max = max(total_max, sum_last)
!                     print *, sum_last
                endif

                print *, pebox(2, nn + cur_id_in_col), pebox(3, nn + cur_id_in_col), pebox(4, nn + cur_id_in_col), pebox(5, nn + cur_id_in_col)
                print *, sum(mask_water_only(pebox(2, nn + cur_id_in_col) : pebox(3, nn + cur_id_in_col), pebox(4, nn + cur_id_in_col) : pebox(5, nn + cur_id_in_col)))

                !print *, "proc ", nn + cur_id_in_col, "set to ", pebox(4, nn + cur_id_in_col), "~", pebox(5, nn + cur_id_in_col)

                row_start = pebox(5, nn + cur_id_in_col) + 1

                !黄承欢于2018-04-01 00:40:55
                cur_row = pebox(5, nn + cur_id_in_col)

                !reset the dst, in case of lll or hhh
                sum_dst = sum_col / (nnn - nn + 1 - (1 + cur_id_in_col))

                !print *, "row_start set to ", row_start, "sum_dst: ", sum_dst

                peboxex(7, nn + cur_id_in_col) = cur_id_in_col
                peboxex(6, nn + cur_id_in_col) = hor_id

                cur_id_in_col = cur_id_in_col + 1

                sum_cur = 0

                if(cur_id_in_col == nnn - nn)then
                    !prepare to exit, because the last row of the next proc must be jm
                    pebox(4, nn + cur_id_in_col) = row_start
                    pebox(5, nn + cur_id_in_col) = jm
!                     print *, sum_cur
!                     total_max = max(total_max, sum_cur)
                    !print *, "for the last, ", pebox(4, nn + cur_id_in_col), pebox(5, nn + cur_id_in_col)
                    exit
                end if
            endif
            sum_last = sum_cur
            cur_row = cur_row + 1
        end do 

        hor_id = hor_id + 1
    end do 

    do nn = 1, n
      sum_cur = sum(mask_water_only(pebox(2, nn) : pebox(3, nn), pebox(4, nn) : pebox(5, nn)))
      total_max = max(total_max, sum_cur)
      print *, sum_cur, pebox(2, nn) , pebox(3, nn), pebox(4, nn) , pebox(5, nn)
    end do

    print *, "total max = ", total_max
    print *, "avepoint = ", avepoint
    print *, "imba eff = ", avepoint / total_max
    print *, "end divide"

    !

  ! --- Ids of up & down.
    pebox(6, :) = MPI_PROC_NULL
    pebox(7, :) = MPI_PROC_NULL
    do nn = 1, n
      i1 = pebox(2, nn)
      i2 = pebox(3, nn)
      j1 = pebox(4, nn)
      j2 = pebox(5, nn)
      do i = 1, n
        if(i1 == pebox(2, i) .and. i2 == pebox(3, i) &
                             .and. j1 == pebox(5, i) + 1)then
          pebox(6, nn) = pebox(1, i)  ! my down
        endif
        if(i1 == pebox(2, i) .and. i2 == pebox(3, i) &
                             .and. j2 == pebox(4, i) - 1)then
          pebox(7, nn) = pebox(1, i)  ! my up
        endif
      enddo
    enddo
  ! --- halosize
    do nn = 1, n
      i1 = pebox(2, nn)
      i2 = pebox(3, nn)
      j1 = pebox(4, nn)
      j2 = pebox(5, nn)
      if(flagxcycle == 1)then
        pebox(11, nn) = max(i1 - halosize, 1)
        pebox(12, nn) = min(i2 + halosize, im)
      else
        pebox(11, nn) = i1 - halosize
        pebox(12, nn) = i2 + halosize
      endif
      pebox(13, nn) = max(1, j1 - halosize)
      pebox(14, nn) = min(jm, j2 + halosize)
      pebox(15, nn) = halosize
      pebox(16, nn) = pebox(12, nn) - pebox(11, nn) + 1 ! isize
      pebox(17, nn) = pebox(14, nn) - pebox(13, nn) + 1 ! jsize
    enddo
  ! --- Ids, js and je of right & left.
    peright = MPI_PROC_NULL; peleft = MPI_PROC_NULL; nl = 0; nr = 0
    do nn = 1, n
      i1 = pebox(2, nn)
      i2 = pebox(3, nn)
      j1 = pebox(4, nn)
      j2 = pebox(5, nn)
      nn_send_start = j1
      nn_send_end = j2
      nn_recv_start = pebox(13, nn)
      nn_recv_end = pebox(14, nn)

      nr = 0 ! ---- My right.
      nr_send = 0
      nr_recv = 0
      do i = 1, n
        if(flagxcycle == 1)then
          if(i2 /= pebox(2, i) - 1)cycle
        else
          if(i2 /= pebox(2, i) - 1 .and. i2 /= im)cycle
          if(i2 == im .and. pebox(2, i) /= 1)cycle
        endif

        have_send_flag = 0
        have_recv_flag = 0

        i_send_start = pebox(4, i)
        i_send_end = pebox(5, i)
        i_recv_start = pebox(13, i)
        i_recv_end = pebox(14, i)
        
        !send to i
        if(max(nn_send_start, i_recv_start) - min(nn_send_end, i_recv_end) <= 0)then
            have_send_flag = 1
            nr_send = nr_send + 1
            perightsend(nn, nr_send, 1) = pebox(1, i)
            perightsend(nn, nr_send, 2) = max(nn_send_start, i_recv_start)
            perightsend(nn, nr_send, 3) = min(nn_send_end, i_recv_end)
            perightsend(nn, nr_send, 2 : 3) = perightsend(nn, nr_send, 2 : 3) - nn_send_start + 1
        endif
        !recv from i
        if(max(nn_recv_start, i_send_start) - min(nn_recv_end, i_send_end) <= 0)then
            have_recv_flag = 1
            nr_recv = nr_recv + 1
            perightrecv(nn, nr_recv, 1) = pebox(1, i)
            perightrecv(nn, nr_recv, 2) = max(nn_recv_start, i_send_start)
            perightrecv(nn, nr_recv, 3) = min(nn_recv_end, i_send_end)
            perightrecv(nn, nr_recv, 2 : 3) = perightrecv(nn, nr_recv, 2 : 3) - nn_recv_start + 1
        endif

        if(max(j2, pebox(5, i)) - min(j1, pebox(4, i)) &
            <= j2 + pebox(5, i) - j1 - pebox(4, i))then
          nr = nr + 1
          peright(nn, nr, 1) = pebox(1, i)
          peright(nn, nr, 2) = max(pebox(4, i), j1) !- j1 + 1
          peright(nn, nr, 3) = min(pebox(5, i), j2) !- j1 + 1
          peright(nn, nr, 4) = max(1 , peright(nn, nr, 2) - halosize) &
                             - pebox(13, nn) + 1
          peright(nn, nr, 5) = min(jm, peright(nn, nr, 3) + halosize) &
                             - pebox(13, nn) + 1
        endif
      enddo

      nl = 0 ! ---- My left.
      nl_send = 0
      nl_recv = 0
      do i = 1, n
        if(flagxcycle == 1)then
          if(i1 /= pebox(3, i) + 1)cycle
        else
          if(i1 /= pebox(3, i) + 1 .and. i1 /= 1)cycle
          if(i1 == 1 .and. pebox(3, i) /= im)cycle
        endif

        have_send_flag = 0
        have_recv_flag = 0

        i_send_start = pebox(4, i)
        i_send_end = pebox(5, i)
        i_recv_start = pebox(13, i)
        i_recv_end = pebox(14, i)

        !send to i
        if(max(nn_send_start, i_recv_start) - min(nn_send_end, i_recv_end) <= 0)then
            have_send_flag = 1
            nl_send = nl_send + 1
            peleftsend(nn, nl_send, 1) = pebox(1, i)
            peleftsend(nn, nl_send, 2) = max(nn_send_start, i_recv_start)
            peleftsend(nn, nl_send, 3) = min(nn_send_end, i_recv_end)
            peleftsend(nn, nl_send, 2 : 3) = peleftsend(nn, nl_send, 2 : 3) - nn_send_start + 1
        endif
        !recv from i
        if(max(nn_recv_start, i_send_start) - min(nn_recv_end, i_send_end) <= 0)then
            have_recv_flag = 1
            nl_recv = nl_recv + 1
            peleftrecv(nn, nl_recv, 1) = pebox(1, i)
            peleftrecv(nn, nl_recv, 2) = max(nn_recv_start, i_send_start)
            peleftrecv(nn, nl_recv, 3) = min(nn_recv_end, i_send_end)
            peleftrecv(nn, nl_recv, 2 : 3) = peleftrecv(nn, nl_recv, 2 : 3) - nn_recv_start + 1
        endif

        if(max(j2, pebox(5, i)) - min(j1, pebox(4, i)) &
            <= j2 + pebox(5, i) - j1 - pebox(4, i))then
          nl = nl + 1
          peleft(nn, nl, 1) = pebox(1, i)
          peleft(nn, nl, 2) = max(pebox(4, i), j1) !- j1 + 1
          peleft(nn, nl, 3) = min(pebox(5, i), j2) !- j1 + 1
          peleft(nn, nl, 4) = max(1 , peleft(nn, nl, 2) - halosize) &
                            - pebox(13, nn) + 1
          peleft(nn, nl, 5) = min(jm, peleft(nn, nl, 3) + halosize) &
                            - pebox(13, nn) + 1
        endif
      enddo

      peboxex(1, nn) = nn
      peboxex(2, nn) = nl_send 
      peboxex(3, nn) = nl_recv
      peboxex(4, nn) = nr_send 
      peboxex(5, nn) = nr_recv

!      print *, "nn =", nn
!      !if(nn == 1)then
!        print *, "nl_send =", nl_send 
!        do i = 1, nl_send 
!            print *, "   ", peleftsend(nn, i, 1 : 3)
!        end do 
!
!        print *, "nl_recv =", nl_recv 
!        do i = 1, nl_recv 
!            print *, "   ", peleftrecv(nn, i, 1 : 3)
!        end do 
!
!        print *, "nr_send =", nr_send 
!        do i = 1, nr_send 
!            print *, "   ", perightsend(nn, i, 1 : 3)
!        end do 
!
!        print *, "nr_recv =", nr_recv 
!        do i = 1, nr_recv 
!            print *, "   ", perightrecv(nn, i, 1 : 3)
!        end do 
    !endif

      pebox(8, nn) = nr
      pebox(9, nn) = nl
      pebox(10, nn) = sum(mask_water_only(i1:i2, j1:j2))
    enddo
    !do i = 1, n
      !write(11, '(17(i5,x))')pebox(:, i)
      !write(12, '(<n>(i5,x))')peleft(i, :, 1)
      !write(13, '(<n>(i5,x))')peright(i, :, 1)
    !enddo
  end subroutine set_pebox
!-------------------------------------------------------------------------------

  end module ympi_mod

!-------------------------------------------------------------------------------
!###############################################################################
