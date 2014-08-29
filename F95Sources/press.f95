module module_press
#ifdef USE_NETCDF_OUTPUT
      use module_LES_write_netcdf
#endif      
      use module_bondFG ! add_module_decls() line 156
      use module_boundp ! add_module_decls() line 156
contains

      subroutine press(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,f,g,h,dt,cn1,cn2l,p,cn2s,cn3l,cn3s,cn4l,cn4s, &
      n,nmax,data20,usum,vsum,wsum)
      use common_sn ! create_new_include_statements() line 102
        real(kind=4), dimension(ip,jp,kp) , intent(In) :: cn1
        real(kind=4), dimension(ip) , intent(In) :: cn2l
        real(kind=4), dimension(ip) , intent(In) :: cn2s
        real(kind=4), dimension(jp) , intent(In) :: cn3l
        real(kind=4), dimension(jp) , intent(In) :: cn3s
        real(kind=4), dimension(kp) , intent(In) :: cn4l
        real(kind=4), dimension(kp) , intent(In) :: cn4s
        character(len=70), intent(In) :: data20
        real(kind=4), intent(In) :: dt
        real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
        real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzn
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: h
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
        integer, intent(In) :: n
        integer, intent(In) :: nmax
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(InOut) :: p
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(Out) :: rhs
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: usum
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: wsum
        integer :: nn
! 
! --only used mac method
      real, parameter  :: pjuge = 0.0001
      integer, parameter  :: nmaxp = 50 ! WV was 50
      real, parameter  :: omega = 1.
!      

      call bondfg(km,jm,f,im,g,h)
! 
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        rhs(i,j,k) = (-u(i-1,j,k)+u(i,j,k))/dx1(i) +(-v(i,j-1,k)+v(i,j,k))/dy1(j) +(-w(i,j,k-1)+w(i, &
      j,k))/dzn(k)
! --stretch
        rhs(i,j,k) = (f(i,j,k)-f(i-1,j,k))/dx1(i) +(g(i,j,k)-g(i,j-1,k))/dy1(j) +(h(i,j,k)-h(i,j, &
      k-1))/dzn(k)+rhs(i,j,k)/dt
! 
      end do
      end do
      end do
! 
!    rhs=0
      rhsav = 0.0
      area = 0.0
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        rhsav = rhsav+dx1(i)*dy1(j)*dzn(k)*rhs(i,j,k)
        area = area +dx1(i)*dy1(j)*dzn(k)
      end do
      end do
      end do
      rhsav = rhsav/area
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        rhs(i,j,k) = rhs(i,j,k)-rhsav
      end do
      end do
      end do
! --SOR
      do l = 1,nmaxp
        sor = 0.0
        do nrd = 0,1
        do k = 1,km
        do j = 1,jm
        do i = 1+mod(k+j+nrd,2),im,2
          reltmp = omega*(cn1(i,j,k) *(cn2l(i)*p(i+1,j,k) +cn2s(i)*p(i-1,j,k) +cn3l(j)*p(i,j+1, &
      k) +cn3s(j)*p(i,j-1,k) +cn4l(k)*p(i,j,k+1) +cn4s(k)*p(i,j,k-1) -rhs(i,j,k))-p(i,j,k))
! 
          p(i,j,k) = p(i,j,k) +reltmp
          sor = sor+reltmp*reltmp
      end do
      end do
      end do
        call boundp1(km,jm,p,im)

      end do
        call boundp2(jm,im,p,km)
#ifndef NO_IO
#ifdef VERBOSE
! --check
!      if ((mod(n-1,10) == 0).and.(mod(l,20) == 0)) then
!        print *, 'time step, iteration step, conv =',n,l,sor
!      end if
#endif      
#endif
!
        if (sor < pjuge) goto 510                         !Break
      end do
  510 call noop

      pav = 0.0
      pco = 0.0
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        pav = pav+p(i,j,k)*dx1(i)*dy1(j)*dzn(k)
        pco = pco+dx1(i)*dy1(j)*dzn(k)
      end do
      end do
      end do
! 
      pav = pav/pco
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        p(i,j,k) = p(i,j,k)-pav
      end do
      end do
      end do
! 
!      print *, "F95: P_SUM_ADJ=",sum(p)
      call boundp1(km,jm,p,im)
      ! print *, "F95: P_SUM_1=",sum(p)
      call boundp2(jm,im,p,km)
!      print *, "F95: P_SUM_BOUND=",sum(p)

! 
#ifndef NO_IO
#ifdef VERBOSE
      if (mod(n-1,20) == 0) then
         write(6,*) '=mac= time step, iteration step, conv =',n,l,sor
      end if
! --check
      if (mod(n-1,20) == 0) then
      do k = 1,km,10
      write(6,*) 'Inflow=',k,'u=',u(0,jm/2,k),v(0,jm/2,k) ,w(0,jm/2,k)
      end do
      do k = 1,km,10
      write(6,*) 'Urbanflow=',k,'u=',u(im/2,jm/2,k),v(im/2,jm/2,k) ,w(im/2,jm/2,k)
      end do
! 
      cflu = 0.
      cflv = 0.
      cflw = 0.
      do k = 1,km
      do j = 1,jm
      do i = 1,im
       cflu = amax1(cflu,abs(u(i,j,k)*dt/dx1(i)))
       cflv = amax1(cflv,abs(v(i,j,k)*dt/dy1(j)))
       cflw = amax1(cflw,abs(w(i,j,k)*dt/dzn(k)))
      end do 
      end do
      end do
      end if

      if (mod(n-1,20) == 0) then
      write(6,*) 'Check_CFL,u*dt/dx,v*dt/dy,w*dt/dz=',cflu,cflv,cflw
      end if
#endif
#endif
      if(mod(n,1000) == 0.or.n == nmax) then
        nn = n/1000
        print *, 'timestep: ',nn,' pressure at centre: ',p(ip/2,jp/2,kp/2), &
                'vel at centre: ', &
                u(ip/2,jp/2,kp/2),v(ip/2,jp/2,kp/2),w(ip/2,jp/2,kp/2)
#ifdef USE_NETCDF_OUTPUT
        call write_to_netcdf_file(p,u,v,w,usum,vsum,wsum,nn)
#endif
#ifndef NO_IO
        open(unit=20,file=data20,form='unformatted',status='unknown')
        write(20) (((u(i,j,k),i=1,im),j=1,jm),k=1,km)
        write(20) (((v(i,j,k),i=1,im),j=1,jm),k=1,km)
        write(20) (((w(i,j,k),i=1,im),j=1,jm),k=1,km)
        write(20) (((p(i,j,k),i=1,im),j=1,jm),k=1,km)
        write(20) (((usum(i,j,k),i=1,im),j=1,jm),k=1,km)
        write(20) (((vsum(i,j,k),i=1,im),j=1,jm),k=1,km)
        write(20) (((wsum(i,j,k),i=1,im),j=1,jm),k=1,km)
      close(unit=20)
#endif
      end if

! 
      return
      end subroutine press




end module module_press

