module mod_data
    implicit none
    type T_DATA
        real*8, pointer :: oHz(:,:,:)
        real*8, pointer :: Huon(:,:,:)
        real*8, pointer :: Hvom(:,:,:)
        real*8, pointer :: W(:,:,:)
        real*8, pointer :: Ta(:,:,:)
        integer, pointer :: Uind(:)
        real*8, pointer :: Dn(:)
        real*8, pointer :: Dm(:)
        real*8, pointer :: Ua(:,:,:)
    end type

    type(T_DATA) :: MYDATA

    integer :: LBi, UBi, LBj, UBj, N, ND
    real*8, parameter :: dt=86.0

contains
subroutine allocate_init_data()
    integer :: i, j
    integer :: seed
    real*8, allocatable :: rndtmp(:)

    seed = 35*16
    N = 300
    ND = 2500
    LBi = 10
    UBi = 960
    LBj = 28
    UBj = 2016

    allocate (MYDATA % oHz(LBi:UBi,LBj:UBj,N))
    allocate (MYDATA % Huon(LBi:UBi,LBj:UBj,N))
    allocate (MYDATA % Hvom(LBi:UBi,LBj:UBj,N))
    allocate (MYDATA % W(LBi:UBi,LBj:UBj,N))
    allocate (MYDATA % Ta(LBi:UBi,LBj:UBj,N))
    allocate (MYDATA % Uind(1:ND))
    allocate (MYDATA % Dn(1:ND))
    allocate (MYDATA % Dm(1:ND))
    allocate (MYDATA % Ua(LBi:UBi,LBj:UBj,N))
    allocate (rndtmp(1:ND))

    call random_seed(seed)

    call random_number(MYDATA%oHz)
    call random_number(MYDATA%Huon)
    call random_number(MYDATA%Hvom)
    call random_number(MYDATA%W)
    call random_number(MYDATA%Ta)
    call random_number(MYDATA%Dn)
    call random_number(MYDATA%Dm)
    call random_number(MYDATA%Ua)
    call random_number(rndtmp)

    MYDATA%Uind(1) = rndtmp(1)*0.1*(UBi-LBi-4)+LBi+2
    DO i=2, ND
    IF(rndtmp(i).gt.0.75) THEN
        MYDATA%Uind(i) = MYDATA%Uind(i-1-MOD(int(rndtmp(i-1)*0.5*i),4))
    ELSE
        MYDATA%Uind(i) = MYDATA%Uind(i-1)+rndtmp(i)*10+0.5
        IF(MYDATA%Uind(i).gt.UBi-2) THEN
            MYDATA%Uind(i) = MOD(MYDATA%Uind(i), UBi-2)+LBi+2
        END IF
    END IF
    END DO
    deallocate (rndtmp)
end subroutine

subroutine flush_cache(t1)
    real*8, allocatable, dimension(:) :: dummyi
    real*8 :: t1

    allocate (dummyi(104857600))
    call random_number(dummyi)
    t1 = dummyi(104857600)
    deallocate (dummyi)

end subroutine

subroutine verify_data()
    real*8, allocatable, dimension(:,:,:) :: Ub
    real*8 :: vmax
    integer :: i,j,k
    allocate (Ub(LBi:UBi,LBj:UBj,N))

    call mpdata_adiff_tile_ver(LBi, UBi, LBj, UBj, MYDATA%oHz, MYDATA%Huon, MYDATA%Hvom, MYDATA%W, MYDATA%Ta, &
                         & MYDATA%Uind, MYDATA%Dn, MYDATA%Dm, Ub)
    vmax = 0.0
    do k=1, N
    do j=LBj+2, UBj-2
    do i=LBi+2, UBi-2
        if(ABS(Ub(i,j,k).gt.1E-6)) THEN
        if(ABS((MYDATA%Ua(i,j,k)-Ub(i,j,k))/Ub(i,j,k)).gt.vmax) THEN
            vmax = (MYDATA%Ua(i,j,k)-Ub(i,j,k))/Ub(i,j,k)
        endif
        elseif (ABS(MYDATA%Ua(i,j,k)-Ub(i,j,k)).gt.vmax) THEN
            vmax = Ub(i,j,k)-MYDATA%Ua(i,j,k)
        endif
    enddo
    enddo
    enddo
    write (*,'(a, f, \)') "          Max error: ", vmax
    if (vmax.le.1E-6) THEN
        write (*,*) "       PASS!"
    else
        write (*,*) "       FAILED..."
    endif
    deallocate (Ub)
end subroutine

subroutine mpdata_adiff_tile_ver (LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
    implicit none
    integer, intent(in) :: LBi, UBi, LBj, UBj
    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)

    integer :: Istr, Iend, Jstr, Jend

    Istr = LBi+2 
    Iend = UBi-2 
    Jstr = LBj+2 
    Jend = UBj-2 

    call stencil_ver (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
end subroutine

subroutine deallocate_data()
    deallocate (MYDATA % oHz)
    deallocate (MYDATA % Huon)
    deallocate (MYDATA % Hvom)
    deallocate (MYDATA % W)
    deallocate (MYDATA % Ta)
    deallocate (MYDATA % Uind)
    deallocate (MYDATA % Dn)
    deallocate (MYDATA % Dm)
    deallocate (MYDATA % Ua)
end subroutine

end module
