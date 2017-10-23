!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: intact
!¼ÆËã²¨ÖÂÀ×ÅµÓ¦Á¦
!--   The following subroutines are created for wave-current interaction  ---

  subroutine intact

  implicit none

  integer :: kh, k, k1, i, i1, j
  real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, theta0, sinth, costh, ekj, ekj1

  do kh=1,kb
  do ic=iys,iyl
  do ia=ixs,ixl
    taubb11(ia,ic,kh)=0.0
    taubb12(ia,ic,kh)=0.0
    taubb22(ia,ic,kh)=0.0
    taubb33(ia,ic,kh)=0.0
  enddo
  enddo
  enddo

  do 100 kh=1,kb
  do 100 ic=iys,iyl
  do 100 ia=ixs,ixl

    if(nsp(ia,ic).ne.1) cycle

    do 200 k=1,kld

      k1=k+1
      i=k-kl+1
      i1=i+1
      dwkk=dwk(k)
      wfk=wf(k,ia,ic)
      wfk1=wf(k1,ia,ic)
      wsk=zpi*wfk
      wsk1=zpi*wfk1
      wkk=wk(k)
      wkk1=wk(k1)

      do 200 j=1,jl

        theta0=thet(j)
        sinth=sin(theta0)
        costh=cos(theta0)

        if (k.lt.kl) then
          ekj=ea(k,j,ia,ic)
          ekj1=ea(k1,j,ia,ic)
        else
          ekj=ea(kl,j,ia,ic)*wkh(i)
          ekj1=ea(kl,j,ia,ic)*wkh(i1)
        endif

        taubb11(ia,ic,kh)=taubb11(ia,ic,kh)                               &
        &                      -wsk**2*costh**2*(ekj+ekj1)                 &
        &                      *exp(2*wkk*zyyz(kh))*dwkk
        taubb12(ia,ic,kh)=taubb12(ia,ic,kh)                               &
        &                      -wsk**2*sinth*costh*(ekj+ekj1)              &
        &                      *exp(2*wkk*zyyz(kh))*dwkk
        taubb22(ia,ic,kh)=taubb22(ia,ic,kh)                               &
        &                      -wsk**2*sinth**2*(ekj+ekj1)                 &
        &                      *exp(2*wkk*zyyz(kh))*dwkk
        taubb33(ia,ic,kh)=taubb33(ia,ic,kh)                               &
        &                      -wsk**2*(ekj+ekj1)                          &
        &                      *exp(2*wkk*zyyz(kh))*dwkk

        !========
        if(abs(zyyz(kh)).gt.d(ia,ic))then
          taubb11(ia,ic,kh)=0.0
          taubb12(ia,ic,kh)=0.0
          taubb22(ia,ic,kh)=0.0
          taubb33(ia,ic,kh)=0.0
        endif
        !=========
    200      continue

  100      continue

  return

  end subroutine intact

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
