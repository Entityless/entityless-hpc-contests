MODULE mpdata_adiff
!
    implicit none
    integer, parameter :: vec_len = 16
CONTAINS

#ifndef USE_MPI

! General implementation funcions.
!
!***********************************************************************
SUBROUTINE mpdata_adiff_tile (LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
!***********************************************************************
!
    USE mod_data, ONLY: ND, N
    implicit none
!
!  Imported variable declarations.
!
    integer, intent(in) :: LBi, UBi, LBj, UBj
!
    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)

!    integer, allocatable :: Uindind(:)
    integer, allocatable :: Ucounter(:)
    integer, allocatable :: Ustp(:)
    integer, allocatable :: Ustp2(:)
    real*8, allocatable :: new_dn(:), new_dm(:)

    integer :: max_ucounter, ilen, vec_count, vec_end
    integer, allocatable :: true_false_arr(:), false_true_arr(:)
    integer, allocatable :: tmp_ctr(:), ft_stp(:), ft_stp2(:), place_counter(:)
    integer, allocatable :: u_offset(:)
    real*8, allocatable :: vec_ucoeff(:), vec_udm(:), vec_udn(:)


!
!  Local variable declarations.
!
    integer :: Istr, Iend, Jstr, Jend, Upos, i, indindlen, l, indval, j, k, cur_ind

    INTEGER:: count1, count_rate1, count_max1
    INTEGER:: count2, count_rate2, count_max2

    CALL SYSTEM_CLOCK(count1, count_rate1, count_max1)

    Istr = LBi+2
    Iend = UBi-2
    Jstr = LBj+2
    Jend = UBj-2

    ilen = Iend - Istr + 1
    vec_count = ilen / vec_len
    vec_end = vec_count * vec_len + Istr - 1

    print *, Istr, Iend, Jstr, Jend, N

    indindlen = Iend - Istr + 1

    allocate(Ucounter(indindlen))
    allocate(Ustp(indindlen))
    allocate(Ustp2(indindlen))
!    allocate(Uindind(ND))
    allocate(new_dn(ND))
    allocate(new_dm(ND))

    Ucounter = 0
    do l = 1, ND
        indval = Uind(l) - Istr + 1
        if(indval > 0 .and. indval < indindlen + 1) Ucounter(indval) = Ucounter(indval) + 1
    enddo

    Ustp(1) = 1
    Ustp2(1) = 1

    max_ucounter = Ucounter(1)

    do i = 2, indindlen
        Ustp(i) = Ustp(i - 1) + Ucounter(i - 1)
        Ustp2(i) = Ustp(i)

        if(max_ucounter < Ucounter(i) .and. i <= vec_count * vec_len) max_ucounter = Ucounter(i)
    end do

    do l = 1, ND
        indval = Uind(l) - Istr + 1
        !Uindind(Ustp2(indval)) = l
        new_dm(Ustp2(indval)) = sin(dm(l))
        new_dn(Ustp2(indval)) = dn(l)
        Ustp2(indval) = Ustp2(indval) + 1
    end do


    !******************************************
    allocate(ft_stp(max_ucounter + 1))
    allocate(ft_stp2(max_ucounter + 1))
    allocate(tmp_ctr(max_ucounter + 1))
    allocate(false_true_arr(vec_count * vec_len + 1))
    allocate(true_false_arr(vec_count * vec_len + 1))
    allocate(u_offset(vec_count * vec_len + 1))
    allocate(place_counter(vec_count * vec_len + 1))

    tmp_ctr = 0

    do i = 1, vec_count * vec_len
        tmp_ctr(Ucounter(i) + 1) = tmp_ctr(Ucounter(i) + 1) + 1
    end do

    do i = 1, max_ucounter + 1
        print *, i - 1, tmp_ctr(i)
    end do
    print *, "max_ucounter", max_ucounter

    ft_stp(1) = 1
    ft_stp2(1) = 1
    do i = 2, max_ucounter + 1
        ft_stp(i) = ft_stp(i - 1) + tmp_ctr(i - 1)
        ft_stp2(i) = ft_stp(i)
    end do

!    do i = 1, max_ucounter + 1
!        do j = ft_stp(i), ft_stp(i) + tmp_ctr(i) - 1
!            place_counter(j) = i - 1
!        end do
!    end do

    do i = 1, vec_count * vec_len
        indval = Ucounter(i) + 1
        false_true_arr(ft_stp2(indval)) = i
        place_counter(ft_stp2(indval)) = Ucounter(i)
        true_false_arr(i) = ft_stp2(indval)

        ft_stp2(indval) = ft_stp2(indval) + 1
    end do

    !do i = 1, max_ucounter + 1
    !    print *, "ft_stp(i): ", i, ft_stp(i), ft_stp2(i)
    !end do

    u_offset(1) = 1
    do i = vec_len + 1, vec_count * vec_len + 1, vec_len
        u_offset(i) = u_offset(i - vec_len) + place_counter(i - 1) * vec_len
    end do

    !print *, u_offset

    allocate(vec_ucoeff(u_offset(vec_count * vec_len + 1)))
    allocate(vec_udn(u_offset(vec_count * vec_len + 1)))
    allocate(vec_udm(u_offset(vec_count * vec_len + 1)))

    do i = 1, vec_count * vec_len, vec_len
        do j = 1, place_counter(i + vec_len - 1)
            do k = 1, vec_len
                cur_ind = u_offset(i) + (j - 1) * vec_len + k - 1 !in vec...
                indval = false_true_arr(i + k - 1) !in false_true_arr, the ori ind
                !some with coeff = 0, or "out of range"
                if(place_counter(i + k - 1) < j) then
                    vec_ucoeff(cur_ind) = 0.0
                    vec_udn(cur_ind) = 1.0
                    vec_udm(cur_ind) = 0.0
                else
                    l = Ustp(indval) + j - 1
                    vec_ucoeff(cur_ind) = 1.0
                    vec_udn(cur_ind) = new_dn(l)
                    vec_udm(cur_ind) = new_dm(l)
                endif
            end do
        end do
    end do

    !print *, Uindind

    CALL SYSTEM_CLOCK(count2, count_rate2, count_max2)
    print *, "it: ", (count2-count1)*1.0/count_rate2
    CALL SYSTEM_CLOCK(count1, count_rate1, count_max1)
    !$omp parallel private(i)

    do i = 1, 1
        call stencil (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, new_dn, new_dm, Ua, Ustp, Ucounter, true_false_arr, false_true_arr, vec_count * vec_len, place_counter, u_offset, u_offset(vec_count * vec_len + 1), vec_ucoeff, vec_udm, vec_udn)
        !$omp barrier
        !$omp master
        CALL SYSTEM_CLOCK(count2, count_rate2, count_max2)
        print *, "tt: ", (count2-count1)*1.0/count_rate2
        count1 = count2
        !$omp end master
    end do

    !$omp end parallel
    return
END SUBROUTINE mpdata_adiff_tile

!
!***********************************************************************
SUBROUTINE stencil (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, new_dn, new_dm, Ua, Ustp, Ucounter, true_false_arr, false_true_arr, vec_total_len, place_counter, u_offset, vec_usz, vec_ucoeff, vec_udm, vec_udn)
!***********************************************************************
!  Imported variable declarations.
!
    USE mod_data, ONLY: N, ND, dt
    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: LBi, UBi, LBj, UBj
    integer, intent(in) :: vec_total_len, vec_usz

    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    real*8, intent(in) :: new_dn(ND)
    real*8, intent(in) :: new_dm(ND)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)

    integer, intent(in) :: Ucounter(UBi - LBi + 1)
    integer, intent(in) :: Ustp(UBi - LBi + 1)
    integer, intent(in) :: true_false_arr(vec_total_len)
    integer, intent(in) :: false_true_arr(vec_total_len)
    integer, intent(in) :: place_counter(vec_total_len)
    integer, intent(in) :: u_offset(vec_total_len)
    real*8, intent(in) :: vec_ucoeff(vec_total_len)
    real*8, intent(in) :: vec_udm(vec_total_len)
    real*8, intent(in) :: vec_udn(vec_total_len)

!
!  Local variable declarations.
!
    integer :: i, j, k, l, m, i_vec, i_inner
    integer :: total_jk_count, my_jk_count, my_start_jk, my_end_jk
    integer :: my_end_j, my_end_k, my_start_j, my_start_k
    integer :: tmp_jstr, tmp_jend, tid, nthreads

    real*8, parameter :: eps = 1.0E-14

    real*8 :: A, B, Um, Vm, X, Y
    real*8 :: AA, BB, AB
    real*8 :: XX, YY, XY
    real*8 :: sig_alfa, sig_beta, sig_gama
    real*8 :: sig_a, sig_b, sig_c

    !real*8, dimension(LBi:UBi,N) :: C
    !real*8, dimension(LBi:UBi,N) :: Wm
    real*8 :: Cv, Wmv

    INTEGER,EXTERNAL :: OMP_GET_NUM_THREADS
    INTEGER,EXTERNAL :: OMP_GET_THREAD_NUM

    integer :: ilen, vec_count, vec_end, left_end, cur_ind

    !dir$ attributes align:64 :: ta_i0j1k0
    !dir$ attributes align:64 :: ta_i1j1k0
    !dir$ attributes align:64 :: ta_i0j0k1
    !dir$ attributes align:64 :: ta_i1j0k1
    !dir$ attributes align:64 :: ta_i0j1k1
    !dir$ attributes align:64 :: ta_i1j1k1
    !dir$ attributes align:64 :: ta_i0j2k1
    !dir$ attributes align:64 :: ta_i1j2k1
    !dir$ attributes align:64 :: ta_i0j1k2
    !dir$ attributes align:64 :: ta_i1j1k2
    !dir$ attributes align:64 :: huon_i1j1k1
    !dir$ attributes align:64 :: hvom_i0j1k1
    !dir$ attributes align:64 :: hvom_i1j1k1
    !dir$ attributes align:64 :: hvom_i0j2k1
    !dir$ attributes align:64 :: hvom_i1j2k1
    !dir$ attributes align:64 :: w_i0j1k0
    !dir$ attributes align:64 :: w_i1j1k0
    !dir$ attributes align:64 :: w_i0j1k1
    !dir$ attributes align:64 :: w_i1j1k1
    !dir$ attributes align:64 :: ohz_i0j0k1
    !dir$ attributes align:64 :: ohz_i1j0k1
    !dir$ attributes align:64 :: ohz_i0j1k1
    !dir$ attributes align:64 :: ohz_i1j1k1
    !dir$ attributes align:64 :: ohz_i0j2k1
    !dir$ attributes align:64 :: ohz_i1j2k1
    !dir$ attributes align:64 :: zero_ones
    !dir$ attributes align:64 :: res_val

    
    real*8 :: ta_i0j1k0(vec_len), ta_i1j1k0(vec_len), ta_i0j0k1(vec_len), ta_i1j0k1(vec_len), ta_i0j1k1(vec_len)
    real*8 :: ta_i1j1k1(vec_len), ta_i0j2k1(vec_len), ta_i1j2k1(vec_len), ta_i0j1k2(vec_len), ta_i1j1k2(vec_len)

    real*8 :: huon_i1j1k1(vec_len)

    real*8 :: hvom_i0j1k1(vec_len), hvom_i1j1k1(vec_len), hvom_i0j2k1(vec_len), hvom_i1j2k1(vec_len)

    real*8 :: w_i0j1k0(vec_len), w_i1j1k0(vec_len), w_i0j1k1(vec_len), w_i1j1k1(vec_len)

    real*8 :: ohz_i0j0k1(vec_len), ohz_i1j0k1(vec_len), ohz_i0j1k1(vec_len), ohz_i1j1k1(vec_len), ohz_i0j2k1(vec_len), ohz_i1j2k1(vec_len)

    real*8 :: zero_ones(vec_len), res_val(vec_len)

    !dir$ attributes align:64 :: Av
    !dir$ attributes align:64 :: Bv
    !dir$ attributes align:64 :: Umv
    !dir$ attributes align:64 :: Vmv
    !dir$ attributes align:64 :: Xv
    !dir$ attributes align:64 :: Yv
    !dir$ attributes align:64 :: AAv
    !dir$ attributes align:64 :: BBv
    !dir$ attributes align:64 :: ABv
    !dir$ attributes align:64 :: XXv
    !dir$ attributes align:64 :: YYv
    !dir$ attributes align:64 :: XYv
    !dir$ attributes align:64 :: sig_alfav
    !dir$ attributes align:64 :: sig_betav
    !dir$ attributes align:64 :: sig_gamav
    !dir$ attributes align:64 :: sig_av
    !dir$ attributes align:64 :: sig_bv
    !dir$ attributes align:64 :: sig_cv
    !dir$ attributes align:64 :: Cvv
    !dir$ attributes align:64 :: Wmvv
    !dir$ attributes align:64 :: abs_um

    real*8 :: Av(vec_len), Bv(vec_len), Umv(vec_len), Vmv(vec_len), Xv(vec_len), Yv(vec_len), AAv(vec_len), BBv(vec_len), ABv(vec_len), XXv(vec_len), YYv(vec_len), XYv(vec_len)
    real*8 :: sig_alfav(vec_len), sig_betav(vec_len), sig_gamav(vec_len), sig_av(vec_len), sig_bv(vec_len), sig_cv(vec_len)
    real*8 :: Cvv(vec_len), Wmvv(vec_len), abs_um(vec_len)

    !dir$ attributes align:64 :: val1
    !dir$ attributes align:64 :: val2
    !dir$ attributes align:64 :: val3
    !dir$ attributes align:64 :: val4
    !dir$ attributes align:64 :: val5
    !dir$ attributes align:64 :: val6
    !dir$ attributes align:64 :: val7
    !dir$ attributes align:64 :: val8
    !dir$ attributes align:64 :: val9
    !dir$ attributes align:64 :: val10
    !dir$ attributes align:64 :: val11
    real*8 :: val1(vec_len), val2(vec_len), val3(vec_len), val4(vec_len), val5(vec_len), val6(vec_len), val7(vec_len), val8(vec_len), val9(vec_len), val10(vec_len), val11(vec_len)

    !dir$ attributes align:64 :: yicsi
    !dir$ attributes align:64 :: yicba
    !dir$ attributes align:64 :: yicshiliu
    !dir$ attributes align:64 :: yicsaner
    !dir$ attributes align:64 :: yicsidt
    !dir$ attributes align:64 :: yicbadt
    !dir$ attributes align:64 :: yicsanerdt
    !dir$ attributes align:64 :: epsv
    !dir$ attributes align:64 :: cv_adjust
    !dir$ attributes align:64 :: zeros
    !dir$ attributes align:64 :: ones
    !dir$ attributes align:64 :: twos
    !dir$ attributes align:64 :: w_adjust
    !dir$ attributes align:64 :: dtdtr
    real*8 :: yicsi(vec_len), yicba(vec_len), yicshiliu(vec_len), yicsaner(vec_len)
    real*8 :: yicsidt(vec_len), yicbadt(vec_len), yicsanerdt(vec_len)
    real*8 :: epsv(vec_len), cv_adjust(vec_len), zeros(vec_len), ones(vec_len), twos(vec_len)
    real*8 :: w_adjust(vec_len), dtdtr(vec_len)

    !dir$ attributes align:64 :: tmp_vec_ua
    !dir$ attributes align:64 :: tmp_vec_val11
    !dir$ attributes align:64 :: tmp_vec_val12
    !dir$ attributes align:64 :: tmp_vec_ua_false
    !dir$ attributes align:64 :: tmp_vec_val11_false
    !dir$ attributes align:64 :: tmp_vec_val12_false
    real*8, allocatable :: tmp_vec_ua(:)
    real*8, allocatable :: tmp_vec_val11(:)
    real*8, allocatable :: tmp_vec_val12(:)
    real*8, allocatable :: tmp_vec_ua_false(:)
    real*8, allocatable :: tmp_vec_val11_false(:)
    real*8, allocatable :: tmp_vec_val12_false(:)

    !dir$ attributes align:64 :: tmp_ohz_j0k1
    !dir$ attributes align:64 :: tmp_ohz_j1k1
    !dir$ attributes align:64 :: tmp_ohz_j2k1
    !dir$ attributes align:64 :: tmp_huon
    !dir$ attributes align:64 :: tmp_hvom_j1k1
    !dir$ attributes align:64 :: tmp_hvom_j2k1
    !dir$ attributes align:64 :: tmp_w_j1k0
    !dir$ attributes align:64 :: tmp_w_j1k1
    !dir$ attributes align:64 :: tmp_ta_j1k0
    !dir$ attributes align:64 :: tmp_ta_j0k1
    !dir$ attributes align:64 :: tmp_ta_j1k1
    !dir$ attributes align:64 :: tmp_ta_j2k1
    !dir$ attributes align:64 :: tmp_ta_j1k2
    real*8, allocatable :: tmp_ohz_j0k1(:), tmp_ohz_j1k1(:), tmp_ohz_j2k1(:)
    real*8, allocatable :: tmp_huon(:)
    real*8, allocatable :: tmp_hvom_j1k1(:), tmp_hvom_j2k1(:)
    real*8, allocatable :: tmp_w_j1k0(:), tmp_w_j1k1(:)
    real*8, allocatable :: tmp_ta_j1k0(:), tmp_ta_j0k1(:), tmp_ta_j1k1(:), tmp_ta_j2k1(:), tmp_ta_j1k2(:)

!    !dec$ attributes fastmem :: tmp_vec_ua
!    !dec$ attributes fastmem :: tmp_vec_val11
!    !dec$ attributes fastmem :: tmp_vec_val12
!    !dec$ attributes fastmem :: tmp_vec_ua_false
!    !dec$ attributes fastmem :: tmp_vec_val11_false
!    !dec$ attributes fastmem :: tmp_vec_val12_false

    yicsi = 0.25
    yicba = 0.125
    yicshiliu = 0.0625
    yicsaner = 0.03125
    yicsidt = 0.25 * dt
    yicbadt = 0.125 * dt
    yicsanerdt = 0.03125 * dt
    epsv = eps
    zeros = 0.0
    ones = 1.0
    twos = 2.0
    dtdtr = 1.0_8 / (dt * dt)

    allocate(tmp_vec_ua(vec_total_len))
    allocate(tmp_vec_val11(vec_total_len))
    allocate(tmp_vec_val12(vec_total_len))
    allocate(tmp_vec_ua_false(vec_total_len))
    allocate(tmp_vec_val11_false(vec_total_len))
    allocate(tmp_vec_val12_false(vec_total_len))

    ilen = Iend - Istr + 1
    vec_count = ilen / vec_len
!    vec_count = 0
    vec_end = vec_count * vec_len + Istr - 1

    allocate(tmp_ohz_j0k1(ilen + 1))
    allocate(tmp_ohz_j1k1(ilen + 1))
    allocate(tmp_ohz_j2k1(ilen + 1))
    allocate(tmp_huon(ilen))
    allocate(tmp_hvom_j1k1(ilen + 1))
    allocate(tmp_hvom_j2k1(ilen + 1))
    allocate(tmp_w_j1k0(ilen + 1))
    allocate(tmp_w_j1k1(ilen + 1))
    allocate(tmp_ta_j1k0(ilen + 1))
    allocate(tmp_ta_j0k1(ilen + 1))
    allocate(tmp_ta_j1k1(ilen + 1))
    allocate(tmp_ta_j2k1(ilen + 1))
    allocate(tmp_ta_j1k2(ilen + 1))

    tid = OMP_GET_THREAD_NUM()
    nthreads = OMP_GET_NUM_THREADS()

    total_jk_count = N * (Jend - Jstr + 1)
    my_jk_count = total_jk_count / nthreads
    my_start_jk = my_jk_count * tid

    if(mod(total_jk_count, nthreads) > tid)then
        my_jk_count = my_jk_count + 1
        my_start_jk = my_start_jk + tid
    ELSE
        my_start_jk = my_start_jk + mod(total_jk_count, nthreads)
    endif

    my_start_jk = my_start_jk + 1 !critical in fortran!!!!!
    my_end_jk = my_start_jk + my_jk_count - 1

    my_start_k = (my_start_jk - 1) / (Jend - Jstr + 1)
    my_start_k = my_start_k + 1
    my_start_j = mod(my_start_jk - 1, (Jend - Jstr + 1)) + Jstr

    my_end_k = (my_end_jk - 1) / (Jend - Jstr + 1)
    my_end_k = my_end_k + 1
    my_end_j = mod(my_end_jk - 1, (Jend - Jstr + 1)) + Jstr

    DO k = my_start_k, my_end_k
        tmp_jstr = Jstr
        if(k == my_start_k) tmp_jstr = my_start_j
        tmp_jend = Jend
        if(k == my_end_k) tmp_jend = my_end_j
        if(k /= 1 .and. k /= N)then
            cv_adjust(1 : vec_len) = yicshiliu(1 : vec_len)
        else
            cv_adjust(1 : vec_len) = yicsi(1 : vec_len)
        endif
        w_adjust(1 : vec_len) = ones(1 : vec_len)
        if(k == N) w_adjust(1 : vec_len) = zeros(1 : vec_len)

        tmp_hvom_j2k1(1 : ilen + 1) = hvom(Istr - 1 : Iend, tmp_jstr, k)
        tmp_ohz_j1k1(1 : ilen + 1) = ohz(Istr - 1 : Iend, tmp_jstr - 1, k)
        tmp_ohz_j2k1(1 : ilen + 1) = ohz(Istr - 1 : Iend, tmp_jstr, k)
        tmp_ta_j1k1(1 : ilen + 1) = ta(Istr - 1 : Iend, tmp_jstr - 1, k)
        tmp_ta_j2k1(1 : ilen + 1) = ta(Istr - 1 : Iend, tmp_jstr, k)

        DO j = tmp_jstr, tmp_jend
            tmp_huon(1 : ilen) = huon(Istr : Iend, j, k)
            tmp_w_j1k0(1 : ilen + 1) = w(Istr - 1 : Iend, j, k - 1)
            tmp_w_j1k1(1 : ilen + 1) = w(Istr - 1 : Iend, j, k)
            tmp_hvom_j1k1(1 : ilen + 1) = tmp_hvom_j2k1(1 : ilen + 1)
            tmp_hvom_j2k1(1 : ilen + 1) = hvom(Istr - 1 : Iend, j + 1, k)
            tmp_ohz_j0k1(1 : ilen + 1) = tmp_ohz_j1k1(1 : ilen + 1)
            tmp_ohz_j1k1(1 : ilen + 1) = tmp_ohz_j2k1(1 : ilen + 1)
            tmp_ohz_j2k1(1 : ilen + 1) = ohz(Istr - 1 : Iend, j + 1, k)
            if(k > 1) tmp_ta_j1k0(1 : ilen + 1) = ta(Istr - 1 : Iend, j, k - 1)
            tmp_ta_j0k1(1 : ilen + 1) = tmp_ta_j1k1(1 : ilen + 1)
            tmp_ta_j1k1(1 : ilen + 1) = tmp_ta_j2k1(1 : ilen + 1)
            tmp_ta_j2k1(1 : ilen + 1) = ta(Istr - 1 : Iend, j + 1, k)
            if(k < N) tmp_ta_j1k2(1 : ilen + 1) = ta(Istr - 1 : Iend, j, k + 1)

            do i_vec = Istr, vec_end, vec_len

                !memory copy, for aligned computation
!                ta_i0j0k1(1 : vec_len) = ta(i_vec - 1 : i_vec + vec_len - 2, j - 1, k)
!                ta_i1j0k1(1 : vec_len) = ta(i_vec : i_vec + vec_len - 1, j - 1, k)
!                ta_i0j1k1(1 : vec_len) = ta(i_vec - 1 : i_vec + vec_len - 2, j, k)
!                ta_i1j1k1(1 : vec_len) = ta(i_vec : i_vec + vec_len - 1, j, k)
!                ta_i0j2k1(1 : vec_len) = ta(i_vec - 1 : i_vec + vec_len - 2, j + 1, k)
!                ta_i1j2k1(1 : vec_len) = ta(i_vec : i_vec + vec_len - 1, j + 1, k)
!
!                hvom_i0j1k1(1 : vec_len) = hvom(i_vec - 1 : i_vec + vec_len - 2, j, k)
!                hvom_i1j1k1(1 : vec_len) = hvom(i_vec : i_vec + vec_len - 1, j, k)
!                hvom_i0j2k1(1 : vec_len) = hvom(i_vec - 1 : i_vec + vec_len - 2, j + 1, k)
!                hvom_i1j2k1(1 : vec_len) = hvom(i_vec : i_vec + vec_len - 1, j + 1, k)
!
!                ohz_i0j0k1(1 : vec_len) = ohz(i_vec - 1 : i_vec + vec_len - 2, j - 1, k)
!                ohz_i1j0k1(1 : vec_len) = ohz(i_vec : i_vec + vec_len - 1, j - 1, k)
!                ohz_i0j1k1(1 : vec_len) = ohz(i_vec - 1 : i_vec + vec_len - 2, j, k)
!                ohz_i1j1k1(1 : vec_len) = ohz(i_vec : i_vec + vec_len - 1, j, k)
!                ohz_i0j2k1(1 : vec_len) = ohz(i_vec - 1 : i_vec + vec_len - 2, j + 1, k)
!                ohz_i1j2k1(1 : vec_len) = ohz(i_vec : i_vec + vec_len - 1, j + 1, k)
!
!                w_i0j1k1(1 : vec_len) = w(i_vec - 1 : i_vec + vec_len - 2, j, k) * w_adjust(1 : vec_len)
!                w_i1j1k1(1 : vec_len) = w(i_vec : i_vec + vec_len - 1, j, k) * w_adjust(1 : vec_len)
!
!                huon_i1j1k1(1 : vec_len) = huon(i_vec : i_vec + vec_len - 1, j, k)

                ta_i0j0k1(1 : vec_len) = tmp_ta_j0k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                ta_i1j0k1(1 : vec_len) = tmp_ta_j0k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                ta_i0j1k1(1 : vec_len) = tmp_ta_j1k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                ta_i1j1k1(1 : vec_len) = tmp_ta_j1k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                ta_i0j2k1(1 : vec_len) = tmp_ta_j2k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                ta_i1j2k1(1 : vec_len) = tmp_ta_j2k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)

                hvom_i0j1k1(1 : vec_len) = tmp_hvom_j1k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                hvom_i1j1k1(1 : vec_len) = tmp_hvom_j1k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                hvom_i0j2k1(1 : vec_len) = tmp_hvom_j2k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                hvom_i1j2k1(1 : vec_len) = tmp_hvom_j2k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)

                ohz_i0j0k1(1 : vec_len) = tmp_ohz_j0k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                ohz_i1j0k1(1 : vec_len) = tmp_ohz_j0k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                ohz_i0j1k1(1 : vec_len) = tmp_ohz_j1k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                ohz_i1j1k1(1 : vec_len) = tmp_ohz_j1k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                ohz_i0j2k1(1 : vec_len) = tmp_ohz_j2k1(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                ohz_i1j2k1(1 : vec_len) = tmp_ohz_j2k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)

                w_i0j1k1(1 : vec_len) = tmp_w_j1k1(i_vec - Istr + 1 : i_vec - Istr + vec_len) * w_adjust(1 : vec_len)
                w_i1j1k1(1 : vec_len) = tmp_w_j1k1(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1) * w_adjust(1 : vec_len)

                huon_i1j1k1(1 : vec_len) = tmp_huon(i_vec - Istr + 1 : i_vec - Istr + vec_len)

                if(k > 1) then
                    ta_i0j1k0(1 : vec_len) = tmp_ta_j1k0(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                    ta_i1j1k0(1 : vec_len) = tmp_ta_j1k0(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                    w_i0j1k0(1 : vec_len) = tmp_w_j1k0(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                    w_i1j1k0(1 : vec_len) = tmp_w_j1k0(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                else
                    ta_i0j1k0(1 : vec_len) = ta_i0j1k1(1 : vec_len)
                    ta_i1j1k0(1 : vec_len) = ta_i1j1k1(1 : vec_len)
                    w_i0j1k0(1 : vec_len) = zeros(1 : vec_len)
                    w_i1j1k0(1 : vec_len) = zeros(1 : vec_len)
                endif

                if(k < N) then
                    ta_i0j1k2(1 : vec_len) = tmp_ta_j1k2(i_vec - Istr + 1 : i_vec - Istr + vec_len)
                    ta_i1j1k2(1 : vec_len) = tmp_ta_j1k2(i_vec - Istr + 2 : i_vec - Istr + vec_len + 1)
                else
                    ta_i0j1k2(1 : vec_len) = ta_i0j1k1(1 : vec_len)
                    ta_i1j1k2(1 : vec_len) = ta_i1j1k1(1 : vec_len)
                endif

                Wmvv(1 : vec_len) = (w_i0j1k0(1 : vec_len) + w_i1j1k0(1 : vec_len) + w_i0j1k1(1 : vec_len) + w_i1j1k1(1 : vec_len)) * yicsidt(1 : vec_len)

                val1(1 : vec_len) = ones(1 : vec_len) / (ta_i0j1k1(1 : vec_len) + ta_i1j1k1(1 : vec_len) + epsv(1 : vec_len))

                Cvv(1 : vec_len) = (ta_i1j1k2(1 : vec_len) - ta_i1j1k0(1 : vec_len) + ta_i0j1k2(1 : vec_len) - ta_i0j1k0(1 : vec_len)) * val1(1 : vec_len) * cv_adjust(1 : vec_len)

                Av(1 : vec_len) = (ta_i1j1k1(1 : vec_len) - ta_i0j1k1(1 : vec_len)) * val1(1 : vec_len)
                Bv(1 : vec_len) = yicsaner(1 : vec_len) * (ta_i1j2k1(1 : vec_len) - ta_i1j0k1(1 : vec_len) + ta_i0j2k1(1 : vec_len) - ta_i0j0k1(1 : vec_len)) * val1(1 : vec_len)
                ABv(1 : vec_len) = Av(1 : vec_len) * Bv(1 : vec_len)

                val6(1 : vec_len) = abs(Av(1 : vec_len))
                val5(1 : vec_len) = ones(1 : vec_len) - val6(1 : vec_len)
                AAv(1 : vec_len) = Av(1 : vec_len) * Av(1 : vec_len)
                val4(1 : vec_len) = ones(1 : vec_len) - AAv(1 : vec_len)

                sig_alfav(1 : vec_len) = ones(1 : vec_len) / (val5(1 : vec_len) + epsv(1 : vec_len))
                sig_betav(1 : vec_len) = -Av(1 : vec_len) / (val4(1 : vec_len) * val5(1 : vec_len) + epsv(1 : vec_len))
                val7(1 : vec_len) = abs(AAv(1 : vec_len) * Av(1 : vec_len))
                sig_gamav(1 : vec_len) = twos(1 : vec_len) * val7(1 : vec_len) / (val5(1 : vec_len) * val4(1 : vec_len) * (ones(1 : vec_len) - val7(1 : vec_len)) + epsv(1 : vec_len))
                val8(1 : vec_len) = ones(1 : vec_len) - abs(ABv(1 : vec_len))
                sig_av(1 : vec_len) = -Bv(1 : vec_len) / (val5(1 : vec_len) * val8(1 : vec_len) + epsv(1 : vec_len))
                val9(1 : vec_len) = abs(Bv(1 : vec_len))
                sig_bv(1 : vec_len) = ABv(1 : vec_len) / (val5(1 : vec_len) * (ones(1 : vec_len) - AAv(1 : vec_len) * val9(1 : vec_len)) + epsv(1 : vec_len)) &
                                    * (val9(1 : vec_len) / (val8(1 : vec_len) + epsv(1 : vec_len)) + twos(1 : vec_len) * Av(1 : vec_len) / (val4(1 : vec_len) - epsv(1 : vec_len)))
                BBv(1 : vec_len) = Bv(1 : vec_len) * Bv(1 : vec_len)
                val10(1 : vec_len) = BBv(1 : vec_len) * val6(1 : vec_len)
                sig_cv(1 : vec_len) = val10(1 : vec_len) / (val5(1 : vec_len) * (ones(1 : vec_len) - val10(1 : vec_len)) * val8(1 : vec_len) + epsv(1 : vec_len))

                Umv(1 : vec_len) = yicbadt(1 : vec_len) * huon_i1j1k1(1 : vec_len) * (ohz_i0j1k1(1 : vec_len) + ohz_i1j1k1(1 : vec_len))

                Vmv(1 : vec_len) = yicsanerdt(1 : vec_len) * (hvom_i0j1k1(1 : vec_len) * (ohz_i0j1k1(1 : vec_len) + ohz_i0j0k1(1 : vec_len)) + hvom_i0j2k1(1 : vec_len) * (ohz_i0j2k1(1 : vec_len) + &
                                    ohz_i0j1k1(1 : vec_len)) + hvom_i1j1k1(1 : vec_len) * (ohz_i1j1k1(1 : vec_len) + ohz_i1j0k1(1 : vec_len)) + hvom_i1j2k1(1 : vec_len) * (ohz_i1j2k1(1 : vec_len) + ohz_i1j1k1(1 : vec_len)))

                abs_um(1 : vec_len) = abs(Umv(1 : vec_len))
                Xv(1 : vec_len) = Umv(1 : vec_len) * Umv(1 : vec_len) - abs_um(1 : vec_len)
                Xv(1 : vec_len) = Xv(1 : vec_len) * Av(1 : vec_len)
                val2(1 : vec_len) = Umv(1 : vec_len) * Vmv(1 : vec_len)
                val3(1 : vec_len) = Cvv(1 : vec_len) * Wmvv(1 : vec_len)
                Xv(1 : vec_len) = val2(1 : vec_len) * Bv(1 : vec_len) + Xv(1 : vec_len)
                Xv(1 : vec_len) = val3(1 : vec_len) * Umv(1 : vec_len) + Xv(1 : vec_len)
                Xv(1 : vec_len) = zeros(1 : vec_len) - Xv(1 : vec_len)

                Yv(1 : vec_len) = Vmv(1 : vec_len) * Vmv(1 : vec_len) - abs(Vmv(1 : vec_len))
                Yv(1 : vec_len) = Yv(1 : vec_len) * Bv(1 : vec_len)
                Yv(1 : vec_len) = val2(1 : vec_len) * Av(1 : vec_len) + Yv(1 : vec_len)
                val11(1 : vec_len) = val3(1 : vec_len) * Vmv(1 : vec_len)
                Yv(1 : vec_len) = val11(1 : vec_len) + Yv(1 : vec_len)
                Yv(1 : vec_len) = zeros(1 : vec_len) - Yv(1 : vec_len)

                XXv(1 : vec_len) = Xv(1 : vec_len) * Xv(1 : vec_len)
                !YYv(1 : vec_len) = Yv(1 : vec_len) * Yv(1 : vec_len)
                !XYv(1 : vec_len) = Xv(1 : vec_len) * Yv(1 : vec_len)

                res_val(1 : vec_len) = sig_alfav(1 : vec_len) * Xv(1 : vec_len) + sig_betav(1 : vec_len) * XXv(1 : vec_len) + sig_gamav(1 : vec_len) * XXv(1 : vec_len) * Xv(1 : vec_len) + &
                                        sig_av(1 : vec_len) * Xv(1 : vec_len) * Yv(1 : vec_len) + sig_bv(1 : vec_len) * XXv(1 : vec_len) * Yv(1 : vec_len) + sig_cv(1 : vec_len) * Xv(1 : vec_len) * Yv(1 : vec_len) * Yv(1 : vec_len)

                res_val(1 : vec_len) = min(abs(res_val(1 : vec_len)), abs(Umv(1 : vec_len)) * sign(ones(1 : vec_len), res_val(1 : vec_len)))

                !res_val(1 : vec_len) = abs(res_val(1 : vec_len))

                tmp_vec_ua(i_vec - Istr + 1 : i_vec - Istr + vec_len) = res_val(1 : vec_len)
                tmp_vec_val11(i_vec - Istr + 1 : i_vec - Istr + vec_len) = val11(1 : vec_len)
                tmp_vec_val12(i_vec - Istr + 1 : i_vec - Istr + vec_len) = Umv(1 : vec_len) * Wmvv(1 : vec_len)
            end do

            do i = 1, vec_total_len
                tmp_vec_ua_false(true_false_arr(i)) = tmp_vec_ua(i)
                tmp_vec_val11_false(true_false_arr(i)) = tmp_vec_val11(i)
                tmp_vec_val12_false(true_false_arr(i)) = tmp_vec_val12(i)
            end do
!            tmp_vec_ua_false(1 : vec_total_len) = tmp_vec_ua(1 : vec_total_len)!wrong!
!            tmp_vec_val11_false(1 : vec_total_len) = tmp_vec_val11(1 : vec_total_len)!wrong!
!            tmp_vec_val12_false(1 : vec_total_len) = tmp_vec_val12(1 : vec_total_len)!wrong!

            do i = 1, vec_total_len, vec_len
                do m = 1, place_counter(i + vec_len - 1)
                    cur_ind = u_offset(i) + (m - 1) * vec_len!in vec...
                    tmp_vec_ua_false(i : i + vec_len - 1) = tmp_vec_ua_false(i : i + vec_len - 1) + (abs(tmp_vec_ua_false(i : i + vec_len - 1)) ** vec_udn(cur_ind : cur_ind + vec_len - 1) * vec_ucoeff(cur_ind : cur_ind + vec_len - 1) * tmp_vec_val12_false(i : i + vec_len - 1) + &
                                                            abs(vec_udm(cur_ind : cur_ind + vec_len - 1) * tmp_vec_val11_false(i : i + vec_len - 1))) * dtdtr(1 : vec_len)
                end do
            end do

            do i = 1, vec_total_len
                ua(Istr - 1 + false_true_arr(i), j, k) = tmp_vec_ua_false(i)
            end do
!            ua(Istr : vec_end, j, k) = tmp_vec_ua_false(1 : vec_total_len)!wrong!

            where((ta(Istr - 1 : vec_end - 1, j, k) .le. 0.0) .or. (ta(Istr : vec_end, j, k) .le. 0.0))
                ua(Istr : vec_end, j, k) = 0.0
            end where

            DO i = vec_end + 1, Iend
                IF ((Ta(i-1,j,k).le.0.0).or.(Ta(i,j,k).le.0.0)) THEN
                    Ua(i,j,k)=0.0
                ELSE
                    if(k == 1)THEN
                        Cv = 0.25*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i-1,j,k+1)-Ta(i-1,j,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
                        Wmv =0.25*dt*(W(i-1,j,k)+W(i,j,k))
                    else if(k == N)then
                        Cv =0.25*  ((Ta(i,j,k)-Ta(i,j,k-1))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
                        Wmv=0.25*dt*( W(i-1,j,k-1)+ W(i,j,k-1))
                    else
                        Cv =0.0625*(Ta(i,j,k+1)-Ta(i,j,k-1)+ &
                                &  Ta(i-1,j,k+1)-Ta(i-1,j,k-1))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
                        Wmv=0.25*dt*((W(i-1,j,k-1)+ W(i-1,j,k))+(W(i,j,k)+ W(i,j,k-1)))
                    endif
                    A=(Ta(i,j,k)-Ta(i-1,j,k))/(Ta(i,j,k)+Ta(i-1,j,k)+eps)
                    B=0.03125*(Ta(i,j+1,k)-Ta(i,j-1,k)+Ta(i-1,j+1,k)-Ta(i-1,j-1,k))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
    !
                    Um=0.125*dt*Huon(i,j,k)*(oHz(i-1,j,k)+oHz(i,j,k))
                    Vm=0.03125*dt*(Hvom(i-1,j,k)*(oHz(i-1,j,k)+oHz(i-1,j-1,k))+Hvom(i-1,j+1,k)*(oHz(i-1,j+1,k)+ & 
                             &     oHz(i-1,j,k))+Hvom(i,j,k)*(oHz(i,j,k)+oHz(i,j-1,k))+Hvom(i,j+1,k)*(oHz(i,j+1,k)+oHz(i,j,k)))
    !
                    X=(ABS(Um)-Um*Um)*A-B*Um*Vm-Cv*Um*Wmv
                    Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-Cv*Vm*Wmv
    !
                    AA=A*A
                    BB=B*B
                    AB=A*B

                    XX=X*X
                    YY=Y*Y
                    XY=X*Y

                    sig_alfa=1.0/(1.0-ABS(A)+eps)
                    sig_beta=-A/((1.0-ABS(A))*(1.0-AA)+eps)
                    sig_gama=2.0*ABS(AA*A)/((1.0-ABS(A))*(1.0-AA)*(1.0-ABS(AA*A))+eps)
                    sig_a=-B/((1.0-ABS(A))*(1.0-ABS(AB))+eps)
                    sig_b=AB/((1.0-ABS(A))*(1.0-AA*ABS(B))+eps)*(ABS(B)/(1.0-ABS(AB)+eps)+2.0*A/(1.0-AA+eps))
                    sig_c=ABS(A)*BB/((1.0-ABS(A))*(1.0-BB*ABS(A))*(1.0-ABS(AB))+eps)

                    Ua(i,j,k)=sig_alfa*X+sig_beta*XX+sig_gama*XX*X+sig_a*XY+sig_b*XX*Y+sig_c*X*YY

                    Ua(i,j,k)=MIN(ABS(Ua(i,j,k)), ABS(Um)*SIGN(1.0,Ua(i,j,k)))
                    !Ua(i,j,k) = abs(Ua(i,j,k))

                    do m = 1, Ucounter(i - Istr + 1)
                        l = Ustp(i - Istr + 1) + m - 1
                        Ua(i,j,k)=Ua(i,j,k)+(abs(Ua(i,j,k))**new_dn(l)*Um*Wmv+ABS(new_dm(l)*Vm*Cv*Wmv))/(dt * dt)
                    enddo
                END IF
            END DO

        END DO
    END DO
!    !$omp end do nowait

    RETURN
END SUBROUTINE stencil 

#else 
! MPI version functions. Not implemented at all.
!
!***********************************************************************
SUBROUTINE distribute_init_data ()
!***********************************************************************
! Distribute the data to all processes from proc0.
! Would not be counted into calculation time.
! Don't do other work inside this subroutine.
END SUBROUTINE
!***********************************************************************
!

!
!***********************************************************************
SUBROUTINE mpdata_adiff_tile_mpi ()
!***********************************************************************
! Distributed implementation for calculation.
END SUBROUTINE
!***********************************************************************
!

!
!***********************************************************************
SUBROUTINE gather_data ()
!***********************************************************************
! Gather the distributed output data Ua to MYDATA%Ua proc0.
! Would not be counted into calculation time.
! Don't do other work inside this subroutine.
END SUBROUTINE
!***********************************************************************
!

#endif

END MODULE mpdata_adiff
