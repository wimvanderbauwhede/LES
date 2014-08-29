module module_ifdata

      use module_bondv1 
      use module_velFG
      use module_les 
      use module_boundp 
#if IFBF == 1        
      use module_feedbfm 
      use module_feedbf 
#endif
#if IADAM == 1
      use module_adam
#endif
contains
      subroutine zero_arrays(cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9, &
              dfu1, dfv1, dfw1, &
              diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,  &
              f,g,h, &
              nou1,nou2,nou3,nou4,nou5,nou6,nou7,nou8,nou9 &
              )  
      use common_sn       
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov9
        real(kind=4), dimension(0:ip,jp,kp) , intent(Out) :: dfu1
        real(kind=4), dimension(ip,0:jp,kp) , intent(Out) :: dfv1
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: dfw1
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu9
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: h
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou9

          do k = 0,kp+1
          do j = 0,jp+2
          do i = -1,ip+2
              cov1(i,j,k) = 0.0
              cov5(i,j,k) = 0.0
              diu1(i,j,k) = 0.0
              diu5(i,j,k) = 0.0
              nou1(i,j,k) = 0.0
              nou5(i,j,k) = 0.0
          end do
          end do
          end do 
          do k = 0,kp+1
          do j = 0,jp+2
          do i = 0,ip+2
              cov2(i,j,k) = 0.0
              cov3(i,j,k) = 0.0
              cov4(i,j,k) = 0.0
              cov6(i,j,k) = 0.0
              cov7(i,j,k) = 0.0
              cov8(i,j,k) = 0.0
              cov9(i,j,k) = 0.0
              diu2(i,j,k) = 0.0
              diu3(i,j,k) = 0.0
              diu4(i,j,k) = 0.0
              diu6(i,j,k) = 0.0
              diu7(i,j,k) = 0.0
              diu8(i,j,k) = 0.0
              diu9(i,j,k) = 0.0
              nou2(i,j,k) = 0.0
              nou3(i,j,k) = 0.0
              nou4(i,j,k) = 0.0
              nou6(i,j,k) = 0.0
              nou7(i,j,k) = 0.0
              nou8(i,j,k) = 0.0
              nou9(i,j,k) = 0.0
          end do
          end do
          end do 
          do k = 1,kp
          do j = 1,jp
          do i = 0,ip
              dfu1(i,j,k) = 0.0
          end do
          end do
          end do 
          do k = 1,kp
          do j = 0,jp
          do i = 1,ip
              dfv1(i,j,k) = 0.0
          end do
          end do
          end do 
          do k = 1,kp
          do j = 1,jp
          do i = 1,ip
              dfw1(i,j,k) = 0.0
          end do
          end do
          end do 

          do k = 0,kp
          do j = 0,jp
          do i = 0,ip
              f(i,j,k) = 0.0
              g(i,j,k) = 0.0
              h(i,j,k) = 0.0
          end do
          end do
          end do 

      end subroutine

      subroutine ifdata( &
#if ICAL == 1
      data30,data31, fold,gold,hold,fghold, time &
#endif
      n,u,im,jm,km,v,w,p,usum,vsum,wsum, &
      delx1,dx1,dy1,dzn,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h,z2,dt, &
      dxs,cov1,cov2,cov3,dfu1,vn,cov4,cov5,cov6,dfv1,cov7,cov8,cov9,dfw1,dzs,nou1,nou5,nou9,nou2, &
      nou3,nou4,nou6,nou7,nou8,bmask1,cmask1,dmask1,alpha,beta,fx,fy,fz,amask1,zbm)
      use common_sn ! create_new_include_statements() line 102
#if ICAL == 1
        character(len=70), intent(In) :: data30
        character(len=70), intent(In) :: data31
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: fghold
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: fold
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: gold
        real(kind=4), dimension(ip,jp,kp) , intent(InOut) :: hold
        real(kind=4), intent(InOut) :: time
#endif
        real(kind=4), intent(In) :: alpha
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(Out) :: amask1
        real(kind=4), intent(In) :: beta
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1) , intent(InOut) :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: cmask1
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: cov9
        real(kind=4), dimension(kp) , intent(Out) :: delx1
        real(kind=4), dimension(0:ip,jp,kp) , intent(Out) :: dfu1
        real(kind=4), dimension(ip,0:jp,kp) , intent(Out) :: dfv1
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: dfw1
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: diu9
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(InOut) :: dmask1
        real(kind=4), intent(In) :: dt
        real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
        real(kind=4), dimension(0:ip) , intent(In) :: dxs
        real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzn
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzs
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: fz
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: h
        integer, intent(InOut) :: im
        integer, intent(InOut) :: jm
        integer, intent(InOut) :: km
        integer, intent(InOut) :: n
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2) , intent(Out) :: nou9
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(InOut) :: p
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1) , intent(Out) :: sm
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: usum
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: v
        real(kind=4), intent(In) :: vn
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(InOut) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: wsum
        real(kind=4), dimension(kp+2) , intent(In) :: z2
        real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1) , intent(InOut) :: zbm
#if IADAM == 1
        character(len=70) :: data21dummy
        integer :: n1,n2
#endif

    if ((im/=ip) .or. (jm/=jp) .or. (km/=kp)) then
            print *, "im,km,km is different from ip,jp,kp, aborting!"
            call exit(-1)
    end if
!        
! 
! 
!      ical =  0;initial start(grid write)
!           =  1;continuous data read,start
!           = 10;continuous data write
!           
#if ICAL == 1
!        if(ical == 1) then
        open(unit=30,file=data30,form='unformatted',status='unknown')
        read(30) n,time
        read(30) (((u(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(30) (((v(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(30) (((w(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(30) (((p(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(30) (((usum(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(30) (((vsum(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(30) (((wsum(i,j,k),i=1,im),j=1,jm),k=1,km) 
        close(30)

        open(unit=31,file=data31,form='unformatted',status='unknown')
        read(31) (((fold(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(31) (((gold(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(31) (((hold(i,j,k),i=1,im),j=1,jm),k=1,km)
        read(31) (((fghold(i,j,k),i=1,im),j=1,jm),k=1,km)
        close(31)
!        end if
#endif
! WV: I added this routine to explicitly set all arrays to zero

        call zero_arrays( &
              cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9, &
              dfu1, dfv1, dfw1, &
              diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,  &
              f,g,h, &
              nou1,nou2,nou3,nou4,nou5,nou6,nou7,nou8,nou9 &
         )

        call bondv1(jm,u,z2,dzn,v,w,km,n,im,dt,dxs)

        call boundp1(km,jm,p,im)

        call boundp2(jm,im,p,km)
        call velfg(km,jm,im,dx1,cov1,cov2,cov3,dfu1,diu1,diu2,dy1,diu3,dzn,vn,f,cov4,cov5,cov6,dfv1, &
      diu4,diu5,diu6,g,cov7,cov8,cov9,dfw1,diu7,diu8,diu9,dzs,h,nou1,u,nou5,v,nou9,w,nou2,nou3, &
      nou4,nou6,nou7,nou8)
#if IFBF == 1        
!        if(ifbf == 1) then
        call feedbf(km,jm,im,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f,g, &
      h)
        call feedbfm(km,jm,im,amask1,bmask1,cmask1,dmask1,zbm,z2,dzn)
!        endif
#endif
        call les(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h)

! --adam
! WV iadam is not defined!
#if IADAM == 1
! WV        if(iadam.eq.1) then
            n1=1
            n2=2
            data21dummy=""
          call adam(n1,n2,data21dummy,fold,im,jm,km,gold,hold,fghold,f,g,h)
! WV        end if
#endif
! 
      end subroutine ifdata

end module module_ifdata

