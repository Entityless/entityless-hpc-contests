!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: nlweight

  subroutine nlweight

  implicit none

  integer :: mr, j, js, j1, j2, is1, icl1, icl2, k, ikn
  real :: alamd, con, deltha, delphi, delphi2, cl1, cl2, ch, acl1, acl2
!  real :: cl11, cl21, alog10pwk, wkk
  real :: cl11, cl21, wkk
  real :: wklp, wklm, wkp, wkp1, wklap, wklap1, wklam, wklam1, wkm, wkm1

  !*      1. setting initial value

  1000      continue
  alamd=0.25
  con=7.862532087
  cong=con*sqrt(g)
  deltha=deltth*180/pi
  delphi=-11.48
  delphi2=33.56

  !*      2.      computation of weights f angular grid

  2000      continue
  cl1=delphi/deltha
  cl2=delphi2/deltha
  ic=1
  do 2001 mr=1,2
  do 2002 j=1,jl
  js=j
  ch=ic*cl1
  j1=jafu(ch,js,jlp1)
  is1=1
  if (ch.lt.0) is1=-1
  j2=j1+is1
  if (j2.gt.jl) j2=1
  if (j2.lt.1) j2=jl
  jp1(mr,j)=j1
  jp2(mr,j)=j2
  js=j
  ch=ic*cl2
  j1=jafu(ch,js,jlp1)
  is1=1
  if (ch.lt.0) is1=-1
  j2=j1+is1
  if (j2.gt.jl) j2=1
  if (j2.lt.1) j2=jl
  jm1(mr,j)=j1
  jm2(mr,j)=j2
  2002      continue
  ic=-1
  2001      continue
  icl1=cl1
  cl1=cl1-icl1
  icl2=cl2
  cl2=cl2-icl2
  acl1=abs(cl1)
  acl2=abs(cl2)
  cl11=1.-acl1
  cl21=1.-acl2
  al11=1.+alamd
  al21=1.-alamd
  al31=al11*al21
  al12=al11**2
  al22=al21**2
  al13=al11**3
  al23=al21**3

  !*      indices and weights f frquency grid

  2100      continue
!yinxq:  alog10pwk=alog10(pwk)
  do 2101 k=1,kl
  wkk=wk(k)
  wks17(k)=sqrt(wkk**17)
  wklp=wkk*al12
  wklm=wkk*al22
  ikn=alog10(al12)/alog10pwk
  ikn=k+ikn
  ikp(k)=ikn
  wkp=wk(ikp(k))
  ikp1(k)=ikp(k)+1
  wkp1=wk(ikp1(k))
  wklap=(wklp-wkp)/(wkp1-wkp)
  wklap1=1.-wklap
  if(ikp(k).le.kl)go to 2104
  ikp(k)=klp1
  ikp1(k)=klp1
  2104      continue
  wp(k,1,1)=wklap1*cl11
  wp(k,1,2)=wklap1*acl1
  wp(k,2,1)=wklap*cl11
  wp(k,2,2)=wklap*acl1
  if(wk(1).lt.wklm)go to 2103
  ikm(k)=1
  ikm1(k)=1
  wklam=0.
  wklam1=0.
  go to 2102
  2103      continue
  ikn=alog10(al22)/alog10pwk
  ikn=k+ikn-1
  if(ikn.lt.1) ikn=1
  ikm(k)=ikn
  wkm=wk(ikm(k))
  ikm1(k)=ikm(k)+1
  wkm1=wk(ikm1(k))
  wklam=(wklm-wkm)/(wkm1-wkm)
  wklam1=1.-wklam
  2102      continue
!  wm(k,1,1)=wklam1*cl11
!  wm(k,1,2)=wklam1*acl1
!  wm(k,2,1)=wklam*cl11
!  wm(k,2,2)=wklam*acl1
   wm(k,1,1)=wklam1*cl21      !jiangxj 2012.4.27
   wm(k,1,2)=wklam1*acl2      !jiangxj 2012.4.27
   wm(k,2,1)=wklam*cl21       !jiangxj 2012.4.27
   wm(k,2,2)=wklam*acl2       !jiangxj 2012.4.27
  2101      continue
  return
  end subroutine nlweight

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: jafu

  integer function jafu(cl,j,ian)

  implicit none

  real, intent(in) :: cl
  integer, intent(in) :: j,ian

  integer :: idph, ja

  idph=cl
  ja=j+idph
  if(ja.le.0)ja=ian+ja-1
  if(ja.ge.ian)ja=ja-ian+1
  jafu=ja

  return

  end function jafu

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
