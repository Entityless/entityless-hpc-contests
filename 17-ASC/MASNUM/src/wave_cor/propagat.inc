!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: propagat

! huangh223, 2017-02-25
! Many indexs and intermediate values remains unchanged in iterations,
! so they should be stored and reused instead of recalculated.

subroutine propagat_init_once
  implicit none
  integer :: j, k, i, iyyz, ixx, jyy, iwk, iwk1, i1, jth
  real :: x0, y0, d0, dddx0, dddy0, duxdx0, duxdy0, duydy0, duydx0, th0
  real :: sinth, costh, wk0, ws0, dk0, cg, cgx, cgy 
  real :: xx, yy, x1, x2, y1, y2, dsidd
  real :: ssr1, ssr2, ssrwk, ssrth, wks, fien, fien1, wk1, wk2, dtth0, ths
  real :: dx, dy, th1, th2
  integer :: iii
  integer :: speread_coor
  
  real :: q_4, r_4, q, r;

  real :: pi_inv_degree, rs_inv
  real :: rslat_inv(iys:iyl)
  real :: wk_inv(kldp1), wk0_inv

  real :: spread_distance_x
  real :: spread_distance_y

  pi_inv_degree = 180.0 / pi
  rs_inv = 1.0 / rs
  rslat_inv = 1.0 / rslat
  wk_inv = 1.0 / wk

  spread_distance_x = x(ixl) - x(ixs)
  spread_distance_x = spread_distance_x / (ixl - ixs)
  spread_distance_x = spread_distance_x * 1.001
  spread_distance_y = y(iyl) - y(iys)
  spread_distance_y = spread_distance_y / (iyl - iys)
  spread_distance_y = spread_distance_y * 1.001
  
  do 100 ia=ixs,ixl
  do 100 ic=iys,iyl

    if(nsp(ia,ic).ne.1) cycle

    x0=x(ia)
    y0=y(ic)
    d0=d(ia,ic)

    dx = deltx(ia)
    dy = delty(ic)

    dddx0 = dddx(ia, ic) ! yinxq
    dddy0 = dddy(ia, ic) ! yinxq

    duxdx0=uxx(ia,ic)
    duxdy0=uxy(ia,ic)
    duydx0=uyx(ia,ic)
    duydy0=uyy(ia,ic)

    do 200 j=1,jl
      th0=thet(j)
      sinth=sin(th0)
      costh=cos(th0)
      do 500 k=1,kl
        wk0=wk(k)
        wk0_inv=wk_inv(k)
        ws0=zpi*wf(k,ia,ic)
        dk0=d0*wk0
        cg=ccg(k,ia,ic)
        cgx=cg*costh
        cgy=cg*sinth

        !******  1.  "the calculation of wave engery-current spreading"
        xx=x0-deltt*(cgx+ux(ia,ic))*rslat_inv(ic)*pi_inv_degree
        yy=y0-deltt*(cgy+uy(ia,ic))*rs_inv*pi_inv_degree

        if(glbflag == 0)then
	        if(xx < x(ixs))xx = xx + 360.
	        if(xx > x(ixl))xx = xx - 360.
        else
	        if(xx.lt.x(ixs))xx=x(ixs)
	        if(xx.ge.x(ixl))xx=x(ixl)
        endif

	ixx = 0;
	!call binary_search(xx, x, ixs, ixl, ixx);

        !if(xx < x(ixs) .or. xx > x(ixl)) cycle

        !always == or <

        speread_coor = (xx - x(ixs)) / spread_distance_x
        speread_coor = speread_coor + ixs
        speread_coor = min(speread_coor, ixl - 1)

        if(xx >= x(speread_coor) .and. xx <= x(speread_coor + 1))then
          ixx = speread_coor
        else
          do iii = speread_coor + 1, ixl - 1
            if(xx >= x(iii) .and. xx <= x(iii+1))then
              ixx = iii; exit
            endif
          enddo
          if(ixx == 0)then
            do iii = speread_coor - 1, ixl - 1, -1
              if(xx >= x(iii) .and. xx <= x(iii+1))then
                ixx = iii; exit
              endif
            enddo
          endif
        endif

        !do iii = ixs, ixl-1
        !  if(xx >= x(iii) .and. xx <= x(iii+1))then
        !  	ixx = iii; exit
        !  endif
        !enddo

        if(ixx < 1 .or. ixx > ixl)cycle
        
        if(yy < y(iys))yy = y(iys)
        if(yy > y(iyl))yy = y(iyl)
        jyy = 0

	!call binary_search(yy, y, iys, iyl, jyy);


        speread_coor = (yy - y(iys)) / spread_distance_y
        speread_coor = speread_coor + iys
        speread_coor = min(speread_coor, iyl - 1)

        if(yy >= y(speread_coor) .and. yy <= y(speread_coor + 1))then
          jyy = speread_coor
        else
          do iii = speread_coor + 1, iyl - 1
            if(yy >= y(iii) .and. yy <= y(iii+1))then
              jyy = iii; exit
            endif
          enddo
          if(jyy == 0)then
            do iii = speread_coor - 1, iyl - 1, -1
              if(yy >= y(iii) .and. yy <= y(iii+1))then
                jyy = iii; exit
              endif
            enddo
          endif
        endif


        !do iii = iys, iyl-1
        !  if(yy >= y(iii) .and. yy <= y(iii+1))then
        !  	jyy = iii; exit
        !  endif
        !enddo

        if(jyy < 0 .or. jyy > iyl)cycle

!-------------------------------------------------------------------------------

        x1=x(ixx)
        x2=x(ixx + 1)
        y1=y(jyy)
        y2=y(jyy + 1)

        !******  2.  "the effect of refraction caused by topography and current"

        if (dk0.lt.40.) then
          dsidd=0.5*g/cosh(dk0)*wk0**2/ws0/cosh(dk0)
        else
          dsidd=0.
        endif
		
        ssr1=(dsidd*dddx0+wk0*costh*duxdx0+wk0*sinth*duydx0)*pi_inv_degree
        ssr2=(dsidd*dddy0+wk0*costh*duxdy0+wk0*sinth*duydy0)*pi_inv_degree
        ssrwk=-(ssr1*costh*rslat_inv(ic)+ssr2*sinth*rs_inv)
        ssrth=(ssr1*sinth*rslat_inv(ic)-ssr2*costh*rs_inv)*wk0_inv
        ssrth=ssrth-cg*costh*tand(y(ic))*rs_inv !?yinxq
        wks=wk0-deltt*ssrwk
        if (wks.lt.0.) wks=0.
        if (wks.le.wkmin) then
          iwk=1
          iwk1=1
          fien=0.
          fien1=1.
          wk1=0.
          wk2=wk(iwk1)
        else
          if (wks.lt.wk(kld)) then
            !call binary_search(wks, wk, 1, kld, iwk);
	    
			do iyyz=1,kld
              if(wks.ge.wk(iyyz).and.wks.lt.wk(iyyz+1))iwk=iyyz
            enddo

            iwk1=iwk+1

            wk1=wk(iwk)
            wk2=wk(iwk+1)
            if (iwk.lt.kl) then
              iwk1=iwk+1
              fien=1.
              fien1=1.
            else
              i=iwk-kl+1
              i1=i+1
              fien=wkh(i)
              fien1=wkh(i1)
              iwk=kl
              iwk1=kl
            endif
          else
            wks=wk(kld)
            iwk=kl
            iwk1=kl
            i=kld-kl+1
            i1=i+1
            fien=wkh(i)
            fien1=wkh(i1)
          endif
        endif
		
        dtth0=deltt*ssrth
        ths=th0-dtth0
        if (ths.ge.(-1.0)*zpi.and.ths.lt.(2.0)*zpi) then
          goto 1200
        else
          ths=th0
        endif
        1200      continue
        if (ths.ge.zpi) ths=ths-zpi
        if (ths.lt. 0.) ths=ths+zpi
		
        jth=int(ths/deltth)+1
        if(jth.eq.jlp1)jth=jl

        th1=thet(jth)
        th2=thet(jth+1)

        !******  3.  "determing the wave energy at the physical space point (xx,yy)
        !******                  and the wave space point (wks,ths)"

        !*  Store indexs and intermediate values
		idxs(1, k, j, ia, ic) = ixx
		idxs(2, k, j, ia, ic) = jyy
		idxs(3, k, j, ia, ic) = jth
		idxs(4, k, j, ia, ic) = iwk
		idxs(5, k, j, ia, ic) = iwk1
		
		call get_bilinear_interpolation_qr(wk1, wk2, th1, th2, wks, ths, q_4, r_4)
		call get_bilinear_interpolation_qr( x1,  x2,  y1,  y2,  xx,  yy,   q,   r)
		
		tmp_values(1, k, j, ia, ic) = q_4
		tmp_values(2, k, j, ia, ic) = r_4
		tmp_values(3, k, j, ia, ic) = q
		tmp_values(4, k, j, ia, ic) = r
		tmp_values(5, k, j, ia, ic) = fien
		tmp_values(6, k, j, ia, ic) = fien1

      500      continue
    200      continue
  100      continue
  
  call c_propagat_init(idxs, tmp_values)
  
end subroutine propagat_init_once