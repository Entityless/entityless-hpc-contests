!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: setspec

! --- Usage: if n = 1, this subroutine will prepare for initial conditions
!            if n = 2, this subroutine will prepare for OBC.

  subroutine setspec(n)

  implicit none

!-------------------------------------------------------------------------------

  integer, intent(in) :: n

  real :: vx, vy, ww, xj, xj0, arlfa, wsj, wkj, &
          theta0, sinth, costh, wk0, wf0, ws0, wl, sigma, alpha
  integer :: j, k

  real, parameter :: gama=3.3, sq3=1.7320508075688772935274463415059

!-------------------------------------------------------------------------------

  do 100 ia=my_hor_start,my_hor_end 
  do 100 ic=my_ver_start,my_ver_end 

    if(nsp(ia,ic).lt.n) cycle
    !if(nsp(ia,ic).lt.n) goto 100
    !      if(nsp(ia,ic).ne.1) goto 100

    vx=wx(ia,ic)
    vy=wy(ia,ic)

    ww=vx**2+vy**2
    w(ia,ic)=sqrt(ww)

    if (w(ia,ic).le.0.) w(ia,ic)=0.9
    !      xj0=0.5*ro(ia,ic)
    !      if (xj0.lt.25000.) xj0=25000.
    xj0=200.*1000.
    xj=g*xj0/(w(ia,ic)**2)
    arlfa=(0.076*(xj**(-0.4)))/pi
    wsj=22.*(xj**(-0.33))*g/w(ia,ic)
    wkj=wsj**2/g

    do 200 j=1,jl
      theta0=thet(j)
      costh=cos(theta0)
      sinth=sin(theta0)
      do 500 k=1,kl
        wk0=wk(k)
        wf0=wf(k,ia,ic)
        ws0=zpi*wf0
        wl=vx*costh+vy*sinth

        if (wl.gt.0) then
          if (ws0.le.wsj) then
            sigma=0.07
          else
            sigma=0.09
          endif
          alpha=arlfa/wk0**4*exp(-1.25*(wkj/wk0)**2)                     &
          &             *gama**(exp(-0.5*((1.-ws0/wsj)/sigma)**2))       &
          &             *(wl/w(ia,ic))**2
        else
          alpha=0.0
        endif

!        e(k,j,ia,ic)=alpha
!        ee(k,j,ia,ic)=max(alpha, small)
        !if(myid == 5 .and. k > 23 .and. j > 10)then
        !  print *, "f:"
        !  print *, "  k, j, ia, ic: ", k, j, ia, ic
        !  
        !endif
        e(k,j,ia,ic)=max(alpha, small)

      500      continue
    200      continue
  100      continue

!-------------------------------------------------------------------------------

  return

!-------------------------------------------------------------------------------

  end subroutine setspec

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
