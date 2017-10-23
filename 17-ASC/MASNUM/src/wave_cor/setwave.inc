!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setwave

  subroutine setwave

  implicit none

  integer :: j, k
  real :: wh, di, wkk, dk, tanhdk, wfk, wsk, deltts
  real :: cgro, yy

  deltth=zpi/float(jl)
  do 777 j=1,jlp1
    thet(j)=(j-1)*deltth
  777  continue

  wh=sqrt((1./pwk)**7)
  wkh(1)=1.
  do 444 k=1,kldp1
    wk(k)=wkmin*(pwk**(k-1)) !discretion of wave number
    if (k.le.kld) dwk(k)=(pwk-1.)*(wk(k)**2)*deltth/2
    if (k.ge.2) wkh(k)=wkh(k-1)*wh
  444  continue

  do 900 ic=iys,iyl
  do 900 ia=ixs,ixl
    di=d(ia,ic)
    do 800 k=1,kldp1
      wkk=wk(k)
      dk=di*wkk
      tanhdk=1.
      if (dk.lt.4.) tanhdk=tanh(dk)
      wfk=sqrt(g*wkk*tanhdk)/zpi
      wsk=wfk*zpi
      wf(k,ia,ic)=wfk
      if (dk.gt.4.) then
        ccg(k,ia,ic)=0.5*wsk/wkk
      else
        if (dk.lt.0.14) then
          ccg(k,ia,ic)=sqrt(g*di)
        else
          ccg(k,ia,ic)=0.5*wsk*(1.+2.*dk/sinh(2.*dk))/wkk
        endif
      endif
    800   continue
  900   continue

  do 910 ic=iys,iyl
  do 910 ia=ixs,ixl
  do 910 k=1,kld
    dwf(k,ia,ic)=(wf(k+1,ia,ic)-wf(k,ia,ic))*deltth/2
  910 continue

  deltts=delttm*60.
  cgro=0.0000091*p*deltts
  do 666 k=1,kl
    grolim(k)=cgro/wk(k)**4
  666   continue

! yinxq: 2011.05.06
!  do ic = iys, iyl
!  do ia = ixs, ixl
! 	  yy=((ic-1)*0.5-90.)*pi/180.0
!	  grids(ia,ic)=rs*0.5*pi/180.*rs*cos(yy)*0.5*pi/180.
!  enddo
!  enddo ! 2013/4/11 18:05:20
! yinxq: 2011.05.06

  return

  end subroutine setwave

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
