!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: inter

  subroutine inter(u1,u2,v1,v2,aa,bb,cc,dd,xt,yt,value)

  implicit none

  real, intent(in) :: u1,u2,v1,v2,aa,bb,cc,dd,xt,yt
  real, intent(out) :: value

  real :: dox, doy, yvalue, r, q, phia, phib, phic, phid

  dox=u2-u1
  doy=v2-v1

  if(dox.eq.0.0.and.doy.eq.0.0)then
    yvalue=(aa+bb+cc+dd)/4.0
  endif

  if(dox.eq.0.0.and.doy.ne.0.0)then
    r=(yt-v1)/doy
    yvalue=r*(cc+dd)/2.0+(1-r)*(aa+bb)/2.0
  endif

  if(dox.ne.0.0.and.doy.eq.0.0)then
    q=(xt-u1)/dox
    yvalue=q*(bb+dd)/2.0+(1-q)*(aa+cc)/2.0
  endif

  if(dox.ne.0.0.and.doy.ne.0.0)then
    q=(xt-u1)/dox
    r=(yt-v1)/doy
    phia=(1-q)*(1-r)
    phib=q*(1-r)
    phic=(1-q)*r
    phid=q*r

    yvalue=aa*phia
    yvalue=yvalue+bb*phib
    yvalue=yvalue+cc*phic
    yvalue=yvalue+dd*phid
  endif

  value=yvalue

  return

  end subroutine inter

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
