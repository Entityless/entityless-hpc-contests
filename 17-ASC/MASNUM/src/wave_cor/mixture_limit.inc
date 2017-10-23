!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: mixture_limit

  subroutine mixture_limit

!-------------------------------------------------------------------------------

!  use wamvar_mod, only : ixl,iyl,kb, iys, ixs, kld, kl, jl, zpi
!  use wamvar_mod, only : nsp, dwk, bv, wf, wk, ea, wkh, zyyz, d

  implicit none

!-------------------------------------------------------------------------------

  integer :: kh, ic, ia, i, i1, j, k, k1
  real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1
  real ::  bv1, bv2, bh11, bh12, bh21, bh22

  real :: theta0, sinth, costh, temp

!-------------------------------------------------------------------------------

	do 100 kh = 1, kb
	do 100 ic = iys, iyl
	do 100 ia = ixs, ixl

    bv1 = 0.0; bv2 = 0.0; 
    bh11 = 0.0; bh12 = 0.0; 
    bh21 = 0.0; bh22 = 0.0; 

    if(nsp(ia, ic).ne.1) cycle

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

!-------------------------------------------------------------------------------

        if(wkk*d(ia,ic).lt.0.14.and.abs(zyyz(kh)).le.d(ia,ic))then

          !bv1=bv1+       (ekj+ekj1)*(d(ia,ic)+zyyz(kh))**2/(d(ia,ic))**2*dwkk
          !bv2=bv2+wsk**2*(ekj+ekj1)*(d(ia,ic)+zyyz(kh))**2/(d(ia,ic))**2*dwkk	
          temp = (ekj+ekj1)*(d(ia,ic)+zyyz(kh))**2/(d(ia,ic))**2*dwkk
          bv1 = bv1 + temp
          bv2 = bv2 + wsk**2 * temp

          !bh11=bh11+       costh**2*(ekj+ekj1)*(1.0+wkk**2*(d(ia,ic)
          !            +zyyz(kh))**2)/(wkk*d(ia,ic))**2*dwkk
          !bh12=bh12+wsk**2*costh**2*(ekj+ekj1)*(1.0+wkk**2*(d(ia,ic)
          !           +zyyz(kh))**2)/(wkk*d(ia,ic))**2*dwkk
          temp = costh**2*(ekj+ekj1)*(1.0+wkk**2*(d(ia,ic)+zyyz(kh))**2) &
               /(wkk*d(ia,ic))**2*dwkk
          bh11 = bh11 + temp
          bh12 = bh12 + wsk**2 * temp

          !bh21=bh21+       sinth**2*(ekj+ekj1)*(1.0+wkk**2*(d(ia,ic)
          !                +zyyz(kh))**2)/(wkk*d(ia,ic))**2*dwkk
          !bh22=bh22+wsk**2*sinth**2*(ekj+ekj1)*(1.0+wkk**2*(d(ia,ic)
          !                +zyyz(kh))**2)/(wkk*d(ia,ic))**2*dwkk
          temp = sinth**2*(ekj+ekj1)*(1.0+wkk**2*(d(ia,ic)+zyyz(kh))**2) &
               /(wkk*d(ia,ic))**2*dwkk
          bh21 = bh21 + temp
          bh22 = bh22 + wsk**2 * temp

        endif

!-------------------------------------------------------------------------------

        if(wkk*d(ia,ic).ge.0.14 .and. wkk*d(ia,ic).lt.4.0 &
           .and. abs(zyyz(kh)).le.d(ia,ic)                )then

          !bv1=bv1+       (ekj+ekj1)*(sinh(wkk*(d(ia,ic)+zyyz(kh))))**2
          !/(sinh(wkk*d(ia,ic)))**2*dwkk
          !bv2=bv2+wsk**2*(ekj+ekj1)*(sinh(wkk*(d(ia,ic)+zyyz(kh))))**2
          !/(sinh(wkk*d(ia,ic)))**2*dwkk
          temp = (ekj+ekj1)*(sinh(wkk*(d(ia,ic)+zyyz(kh))))**2 &
               / (sinh(wkk*d(ia,ic)))**2*dwkk
          bv1 = bv1 + temp
          bv2 = bv2 + wsk**2 * temp


          !bh11=bh11+       costh**2*(ekj+ekj1)*(cosh(wkk*(d(ia,ic)
          !+zyyz(kh))))**2/(sinh(wkk*d(ia,ic)))**2*dwkk
          !bh12=bh12+wsk**2*costh**2*(ekj+ekj1)*(cosh(wkk*(d(ia,ic)
          !+zyyz(kh))))**2/(sinh(wkk*d(ia,ic)))**2*dwkk
          temp = costh**2*(ekj+ekj1)*(cosh(wkk*(d(ia,ic)+zyyz(kh))))**2 &
               / (sinh(wkk*d(ia,ic)))**2*dwkk
          bh11 = bh11 + temp
          bh12 = bh12 + wsk**2 * temp

          !bh21=bh21+       sinth**2*(ekj+ekj1)*(cosh(wkk*(d(ia,ic)
          !+zyyz(kh))))**2/(sinh(wkk*d(ia,ic)))**2*dwkk
          !bh22=bh22+wsk**2*sinth**2*(ekj+ekj1)*(cosh(wkk*(d(ia,ic)
          !+zyyz(kh))))**2/(sinh(wkk*d(ia,ic)))**2*dwkk
          temp = sinth**2*(ekj+ekj1)*(cosh(wkk*(d(ia,ic)+zyyz(kh))))**2 &
               / (sinh(wkk*d(ia,ic)))**2*dwkk
          bh21 = bh21 + temp
          bh22 = bh22 + wsk**2 * temp

        endif

!-------------------------------------------------------------------------------

        if(wkk*d(ia,ic).ge.4.0.and.abs(zyyz(kh)).le.d(ia,ic))then

          !bv1=bv1+       (ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          !bv2=bv2+wsk**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          temp = (ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          bv1 = bv1 + temp
          bv2 = bv2 + wsk**2 * temp


          !bh11=bh11+       costh**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          !bh12=bh12+wsk**2*costh**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          temp = costh**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          bh11 = bh11 + temp     
          bh12 = bh12 + wsk**2 * temp

          !bh21=bh21+       sinth**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          !bh22=bh22+wsk**2*sinth**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          temp = sinth**2*(ekj+ekj1)*exp(2*wkk*zyyz(kh))*dwkk
          bh21 = bh21 + temp      
          bh22 = bh22 + wsk**2 * temp

        endif

!-------------------------------------------------------------------------------

      200	continue

	    bvl(ia, ic, kh) = (bv1 * bv2  )**0.5
	    bh1(ia, ic, kh) = (bh11 * bh12)**0.5
	    bh2(ia, ic, kh) = (bh21 * bh22)**0.5

      if(abs(zyyz(kh)) .gt. d(ia, ic))then
        bvl(ia, ic, kh) = 0.0
        bh1(ia, ic, kh) = 0.0
        bh2(ia, ic, kh) = 0.0
      endif

  100	continue

!-------------------------------------------------------------------------------

	return

!-------------------------------------------------------------------------------

	end subroutine mixture_limit

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
