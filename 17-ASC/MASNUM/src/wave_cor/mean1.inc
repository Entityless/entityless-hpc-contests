!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: mean1
! --- To compute the characters of wave for output.

  subroutine mean1

!-------------------------------------------------------------------------------

  implicit none

  integer :: i, j, k, k1, i1
  real :: aets, aetc, thmax, akmax, eemax, eformax, dwkk, dwfk
  real :: wfk, wfk1, wsk, wsk1, wkk, wkk1, eef0, eekj, eekj1, eekjth
  real :: sinth, costh, aett, chbh

  real :: hmax, iahm, ichm

!-------------------------------------------------------------------------------

  hmax=0.
  iahm=1
  ichm=1

!-------------------------------------------------------------------------------

  do 100 ic=iys,iyl
  do 100 ia=ixs,ixl

    ae(ia,ic)=0.
    asi(ia,ic)=0.
    awf(ia,ic)=0.
    awk(ia,ic)=0.
    ark(ia,ic)=0.
    ape(ia,ic)=0.
    aet(ia,ic)=0.
    hb(ia,ic)=0.
    hbb(ia,ic)=0.
    h1_3(ia,ic)=0.
    aets=0.0
    aetc=0.0

    if(nsp(ia,ic).ne.1) cycle

    thmax=0.0
    akmax=0.0
    eemax=-999
    eformax=-999
    do 200 k=1,kl
      k1=k+1
      i=k-kl+1
      i1=i+1
      dwkk=dwk(k)
      dwfk=dwf(k,ia,ic)
      wfk=wf(k,ia,ic)
      wfk1=wf(k1,ia,ic)
      !      wfk=fr(k)
      !      wfk1=fr(k1)
      wsk=zpi*wfk
      wsk1=zpi*wfk1
      wkk=wk(k)
      wkk1=wk(k1)
      eef0=0.0
      do 300 j=1,jl
        if (k.lt.kl) then
          eekj =ee(k,j,ia,ic)
          eekj1=ee(k1,j,ia,ic)
        else
          eekj =ee(kl,j,ia,ic)*wkh(i)
          eekj1=ee(kl,j,ia,ic)*wkh(i1)
        endif
        eekjth=eekj*deltth
        eef0=eef0+eekjth
        sinth=sin(thet(j))
        costh=cos(thet(j))
        !      atfs=atfs+sinth*eekjth
        !      atfc=atfs+costh*eekjth
        ae(ia,ic)=ae(ia,ic)+(eekj+eekj1)*dwkk
        asi(ia,ic)=asi(ia,ic)+(eekj/wfk**2+eekj1/wfk1**2)*dwkk
        ape(ia,ic)=ape(ia,ic)+(eekj*wsk**2+eekj1*wsk1**2)*dwkk
        awf(ia,ic)=awf(ia,ic)+(eekj*wfk+eekj1*wfk1)*dwkk
        awk(ia,ic)=awk(ia,ic)+(eekj*wkk+eekj1*wkk1)*dwkk
        !      ark(ia,ic)=ark(ia,ic)+(eekj/sqrt(wkk)+eekj1/sqrt(wkk1))*dwkk
        aets=aets+(eekj+eekj1)*wk(k)*sinth*dwkk
        aetc=aetc+(eekj+eekj1)*wk(k)*costh*dwkk
      300      continue
    !      ef(k,ia,ic)=zpi*wkk/ccg(k,ia,ic)*eef0
    200      continue
    !^^^^^^^^^^^^^^^^^^^^^^^
    ape(ia,ic)=tztz*zpi/sqrt(ape(ia,ic)/ae(ia,ic))
    tpf(ia,ic)=(asi(ia,ic)/ae(ia,ic))*(awf(ia,ic)/ae(ia,ic))
    !        tpf(ia,ic)=asi(ia,ic)*awf(ia,ic)/ae(ia,ic)**2

    if(abs(aetc).lt.0.000001)aetc=0.00001

    aett=atan2d(aets,aetc)
    if (aett.lt.0.) aett=360.+aett
    aet(ia,ic)=aett
    awk(ia,ic)=awk(ia,ic)/ae(ia,ic)
    h1_3(ia,ic)=4.*sqrt(ae(ia,ic))
    !^^^^^^^^^^^^^^^^^^^^^^^
    !        hb(ia,ic)=zpi/awk(ia,ic)*0.12*
    !     &          tanh(d(ia,ic)*awk(ia,ic))/1.6726
    !^^^^^^^^^^^^^^^^^^^^^^^
    hb(ia,ic)=zpi/awk(ia,ic)*0.142*                                   &
    &            tanh(d(ia,ic)*awk(ia,ic))
    hbb(ia,ic)=0.78125*d(ia,ic)/1.6726
    !      hbb(ia,ic)=0.78125*d(ia,ic)/1.5864792
    if(hb(ia,ic).gt.hbb(ia,ic)) hb(ia,ic)=hbb(ia,ic)

    if(h1_3(ia,ic).gt.hb(ia,ic)) then
      chbh=(hb(ia,ic)/h1_3(ia,ic))**2
      !      do 6783 j=1,nth
      !      do 6784 k=1,nwk
      do 6783 j=1,jl
      do 6784 k=1,kl
        ee(k,j,ia,ic)=max(chbh*ee(k,j,ia,ic), small)
      6784      continue
      6783      continue
      h1_3(ia,ic)=hb(ia,ic)
    endif

	! --- This constrain is from Hua F. et al, 2004. Advances in marine science
	tpf(ia, ic) = min(1.41*ape(ia, ic), tpf(ia, ic)) 

    if(h1_3(ia,ic).gt.hmax) then
      hmax=h1_3(ia,ic)
      iahm=ia
      ichm=ic
    endif

!    if(h1_3(ia,ic).gt.h1_3maxt(ia,ic)) then
!      h1_3maxt(ia,ic)=h1_3(ia,ic)
!      tahmt(ia,ic)=ape(ia,ic)
!      yhmt(ia,ic)=itime(1)
!      monhmt(ia,ic)=itime(2)
!      dhmt(ia,ic)=itime(3)
!      hhmt(ia,ic)=itime(4)
!      mhmt(ia,ic)=itime(5)
!    endif

  100      continue

!-------------------------------------------------------------------------------

  do ia=ixs,ixl
    do ic=iys,iyl
      if(nsp(ia,ic).ne.0)then
        if(h1_3(ia,ic).le.0.05)h1_3(ia,ic)=0.05
      endif
    enddo
  enddo

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine mean1

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
