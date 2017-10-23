!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: mixture

  subroutine mixture

!  use wamvar_mod, only : ixl,iyl,kb, iys, ixs, kld, kl, jl, zpi
!  use wamvar_mod, only : nsp, dwk, bv, wf, wk, ea, wkh, zyyz, d

  implicit none

  integer :: kh, ic, ia, i, i1, j, k, k1
  real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1, bv1, bv2, bv3

  do 100 kh=1,kb
  do 100 ic=iys,iyl
  do 100 ia=ixs,ixl

    bv1 = 0.0; bv2 = 0.0; bv3 = 0.0;
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
        if (k.lt.kl) then
          ekj=ea(k,j,ia,ic)
          ekj1=ea(k1,j,ia,ic)
        else
          ekj=ea(kl,j,ia,ic)*wkh(i)
          ekj1=ea(kl,j,ia,ic)*wkh(i1)
        endif

        bv1=bv1+(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
        bv2 = bv2+wsk**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
        bv3=bv3 +wkk*wsk**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
    200      continue

    bv(ia,ic,kh)=bv1/sqrt(bv2)*bv3

    !=======
    if(abs(zyyz(kh)).gt.d(ia,ic))then
      bv(ia,ic,kh)=0.0
    endif
    !========
  100      continue

  return

  end subroutine mixture

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
