!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------

 ! 计算谱矩：用于计算源函数
 subroutine mean2_async(ia_start,ia_end,ic_start,ic_end)

  implicit none

!-------------------------------------------------------------------------------
  integer :: ia_start,ia_end,ic_start,ic_end
  integer :: k, k1, i, i1, j
  real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1

!-------------------------------------------------------------------------------

  do 100 ic=ic_start,ic_end
  do 100 ia=ia_start,ia_end

    ae(ia,ic)=0.
    asi(ia,ic)=0.
    awf(ia,ic)=0.
    awk(ia,ic)=0.
    ark(ia,ic)=0.

    if(nsp(ia,ic).ne.1) cycle

    do 200 k=1,kld
      k1=k+1
      i=k-kl+1
      i1=i+1
      dwkk=dwk(k)
      wfk=wf(k,ia,ic)
      wfk1=wf(k1,ia,ic)
      !      wfk=fr(k)
      !      wfk1=fr(k1)
      wsk=zpi*wfk
      wsk1=zpi*wfk1
      wkk=wk(k)
      wkk1=wk(k1)
      do 200 j=1,jl
        if (k.lt.kl) then
          ekj=e(k,j,ia,ic)
          ekj1=e(k1,j,ia,ic)
        else
          ekj=e(kl,j,ia,ic)*wkh(i)
          ekj1=e(kl,j,ia,ic)*wkh(i1)
        endif
        ae(ia,ic)=ae(ia,ic)+(ekj+ekj1)*dwkk
        awf(ia,ic)=awf(ia,ic)+(ekj*wfk+ekj1*wfk1)*dwkk
        asi(ia,ic)=asi(ia,ic)+(ekj/wsk+ekj1/wsk1)*dwkk
        awk(ia,ic)=awk(ia,ic)+(ekj*wkk+ekj1*wkk1)*dwkk
        ark(ia,ic)=ark(ia,ic)+(ekj/sqrt(wkk)+ekj1/sqrt(wkk1))*dwkk
    200      continue

    asi(ia,ic)=ae(ia,ic)/asi(ia,ic)
    awf(ia,ic)=awf(ia,ic)/ae(ia,ic)
    awk(ia,ic)=awk(ia,ic)/ae(ia,ic)
    ark(ia,ic)=(ae(ia,ic)/ark(ia,ic))**2

  100      continue

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine mean2_async
 
 
 
 
 
  subroutine mean2

  implicit none

!-------------------------------------------------------------------------------

  integer :: k, k1, i, i1, j
  real :: dwkk, wfk, wfk1, wsk, wsk1, wkk, wkk1, ekj, ekj1

!-------------------------------------------------------------------------------

  do 100 ic=iys,iyl
  do 100 ia=ixs,ixl

    ae(ia,ic)=0.
    asi(ia,ic)=0.
    awf(ia,ic)=0.
    awk(ia,ic)=0.
    ark(ia,ic)=0.

    if(nsp(ia,ic).ne.1) cycle

    do 200 k=1,kld
      k1=k+1
      i=k-kl+1
      i1=i+1
      dwkk=dwk(k)
      wfk=wf(k,ia,ic)
      wfk1=wf(k1,ia,ic)
      !      wfk=fr(k)
      !      wfk1=fr(k1)
      wsk=zpi*wfk
      wsk1=zpi*wfk1
      wkk=wk(k)
      wkk1=wk(k1)
      do 200 j=1,jl
        if (k.lt.kl) then
          ekj=e(k,j,ia,ic)
          ekj1=e(k1,j,ia,ic)
        else
          ekj=e(kl,j,ia,ic)*wkh(i)
          ekj1=e(kl,j,ia,ic)*wkh(i1)
        endif
        ae(ia,ic)=ae(ia,ic)+(ekj+ekj1)*dwkk
        awf(ia,ic)=awf(ia,ic)+(ekj*wfk+ekj1*wfk1)*dwkk
        asi(ia,ic)=asi(ia,ic)+(ekj/wsk+ekj1/wsk1)*dwkk
        awk(ia,ic)=awk(ia,ic)+(ekj*wkk+ekj1*wkk1)*dwkk
        ark(ia,ic)=ark(ia,ic)+(ekj/sqrt(wkk)+ekj1/sqrt(wkk1))*dwkk
    200      continue

    asi(ia,ic)=ae(ia,ic)/asi(ia,ic)
    awf(ia,ic)=awf(ia,ic)/ae(ia,ic)
    awk(ia,ic)=awk(ia,ic)/ae(ia,ic)
    ark(ia,ic)=(ae(ia,ic)/ark(ia,ic))**2

  100      continue

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine mean2

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
