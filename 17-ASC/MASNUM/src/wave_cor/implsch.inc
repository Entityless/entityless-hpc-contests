!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------

subroutine c_init_once()
	implicit none

    integer :: ia, ic, k
    real :: ssbo, wk0, dk, d0, sbo


    sbo=0.038*2./g

    do ic=iys,iyl
        do ia=ixs,ixl
            if(nsp(ia,ic).ne.1) cycle
            d0=d(ia,ic)
            do k=1,kl
                wk0=wk(k)
                dk=d0*wk0
                if (dk.ge.30.)then
                    ssbo=0.
                else
                    ssbo=-abo*sbo*wk0/sinh(2.*dk)
                endif
                ssbos(k, ia, ic) = ssbo
            end do 
        end do 
    end do 

	call c_implement_init_once( &
		ixs, ixl, iys, iyl, wx, wy, wk, nsp, kpmt0, kakt0, ks0, &
		fconst0, awk, ae, asi, awf, ark, wf, dwk, wkh, e, w, ee, dse, &
		ssbos, deltt, deltt5, grolim, pein, peds, pebo, sein, sebo, seds, &
		uxx, uxy, uyx, uyy, thet, ccg, d, enh, wp, wm, ikp, ikp1, ikm, &
		ikm1, wks17, jp1, jp2, jm1, jm2, cong, al31, al21, al13, al23, al11, &
		x, y , &
    ape, aet, hb, hbb, h1_3, dwf, tpf, &
    deltth ,&
    glbflag ,&
    ix2 &
	)

	return
end subroutine c_init_once

  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
