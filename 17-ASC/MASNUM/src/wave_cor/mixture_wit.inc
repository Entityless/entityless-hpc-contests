!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: mixture_wit

  subroutine mixture_wit(gama_wt)

!-------------------------------------------------------------------------------

  implicit none

!-------------------------------------------------------------------------------

  real, intent(in) :: gama_wt
  integer :: kh, ic, ia, i, i1, j, k, k1
  real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1, dwv
  real, allocatable :: dwv3d(:, :, :)
  real :: theta0, sinth, costh, temp
  real :: amp2, beta_wt, sigma_wt

!-------------------------------------------------------------------------------

  allocate(dwv3d(ixs:ixl, iys:iyl, kb))

!-------------------------------------------------------------------------------
	
	beta_wt = 1.0
!	gama_wt = 2.0
  sigma_wt = 0.5
  
!-------------------------------------------------------------------------------

	do 100 kh = 1, kb
	do 100 ic = iys, iyl
	do 100 ia = ixs, ixl

    dwv = 0.0

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
          temp = (ekj+ekj1)*wkk**2*wsk**2*(0.5+(wkk*(d(ia,ic)+zyyz(kh)))**2)**2 &
		       / (wkk*d(ia,ic))**2*dwkk
          dwv = dwv + temp
        endif

!-------------------------------------------------------------------------------

        if(wkk*d(ia,ic).ge.0.14 .and. wkk*d(ia,ic).lt.4.0 &
           .and. abs(zyyz(kh)).le.d(ia,ic))then
          temp = (ekj+ekj1)*wkk**2*wsk**2*cosh(2*wkk*(d(ia,ic)+zyyz(kh))) &
		       / (sinh(wkk*d(ia,ic)))**2*dwkk
          dwv = dwv + temp
        endif

!-------------------------------------------------------------------------------

        if(wkk*d(ia,ic).ge.4.0.and.abs(zyyz(kh)).le.d(ia,ic))then
          temp = (ekj+ekj1)*wkk**2*wsk**2*2*exp(2*wkk*zyyz(kh))*dwkk
          dwv = dwv + temp
        endif

!-------------------------------------------------------------------------------

      200	continue

      dwv =  dwv**0.5
      dwv3d(ia, ic, kh) = dwv
      
      amp2 = (h1_3(ia,ic)/4.0)**2
      bv_wtv(ia, ic, kh) = (7.0/8)**0.5*beta_wt*amp2*dwv3d(ia, ic, 1)  &
                           *(dwv3d(ia, ic, kh)/dwv3d(ia, ic, 1))**gama_wt
      bv_wtd(ia, ic, kh) = bv_wtv(ia, ic, kh)/sigma_wt

      if(abs(zyyz(kh)) .gt. d(ia, ic))then
        bv_wtv(ia, ic, kh) = 0.0
        bv_wtd(ia, ic, kh) = 0.0
      endif

  100	continue

  deallocate(dwv3d)

!-------------------------------------------------------------------------------

	return

!-------------------------------------------------------------------------------

	end subroutine mixture_wit

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
