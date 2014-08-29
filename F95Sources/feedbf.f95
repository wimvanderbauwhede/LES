module module_feedbf

contains

      subroutine feedbf(km,jm,im,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f, &
      g,h)
      use common_sn ! create_new_include_statements() line 102
        real(kind=4), intent(In) :: alpha
        real(kind=4), intent(In) :: beta
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1) , intent(In) :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: cmask1
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(In) :: dmask1
        real(kind=4), intent(In) :: dt
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fz
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: h
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: usum
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: wsum
! 
#ifdef WV_DEBUG
    print *, 'F95 UVWSUMSUM after bondv1:',sum(usum)+sum(vsum)+sum(wsum)
    print *, 'F95 USUMSUM after bondv1:',sum(usum)
    print *, 'F95 VSUMSUM after bondv1:',sum(vsum)
    print *, 'F95 WSUMSUM after bondv1:',sum(wsum)

#endif
! 
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        usum(i,j,k) = (usum(i,j,k)+u(i,j,k))*bmask1(i,j,k)
        vsum(i,j,k) = (vsum(i,j,k)+v(i,j,k))*cmask1(i,j,k)
        wsum(i,j,k) = (wsum(i,j,k)+w(i,j,k))*dmask1(i,j,k)

        f1x = alpha*usum(i,j,k)*dt
        f1y = alpha*vsum(i,j,k)*dt
        f1z = alpha*wsum(i,j,k)*dt
        f2x = beta*u(i,j,k)*bmask1(i,j,k)
        f2y = beta*v(i,j,k)*cmask1(i,j,k)
        f2z = beta*w(i,j,k)*dmask1(i,j,k)

        fx(i,j,k) = f1x+f2x
        fy(i,j,k) = f1y+f2y
        fz(i,j,k) = f1z+f2z
       end do
       end do
       end do

       do k = 1,km
       do j = 1,jm
       do i = 1,im
         f(i,j,k) = f(i,j,k)+fx(i,j,k)
         g(i,j,k) = g(i,j,k)+fy(i,j,k)
         h(i,j,k) = h(i,j,k)+fz(i,j,k)
       end do
       end do
       end do
#ifdef WV_DEBUG
    print *, 'F95 FGHSUM after feedbf:',sum(f)+sum(g)+sum(h)
    print *, 'F95 FSUM after feedbf:',sum(f)
    print *, 'F95 GSUM after feedbf:',sum(g)
    print *, 'F95 HSUM after feedbf:',sum(h)
    print *, 'F95 UVWSUMSUM after feedbf:',sum(usum)+sum(vsum)+sum(wsum)
        print *, 'F95 USUMSUM after feedbf:',sum(usum)
    print *, 'F95 VSUMSUM after feedbf:',sum(vsum)
    print *, 'F95 WSUMSUM after feedbf:',sum(wsum)

print *, 'F95 UVWSUM after feedbf:', sum(u)+sum(v)+sum(w)
    print *, 'F95 USUM after feedbf:', sum(u)
    print *, 'F95 VSUM after feedbf:', sum(v)
    print *, 'F95 WSUM after feedbf:', sum(w)

#endif
! 
      return
      end subroutine feedbf




end module module_feedbf

