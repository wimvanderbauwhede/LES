module module_init

      use module_feedbfm ! add_module_decls() line 156
contains

      subroutine init(km,jm,im,u,v,w,p,cn2s,dxs,cn2l,cn3s,dys,cn3l,dzs,cn4s,cn4l,cn1,amask1, &
      bmask1,cmask1,dmask1,zbm,z2,dzn)
      use common_sn ! create_new_include_statements() line 102
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(Out) :: amask1
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1) , intent(Out) :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(Out) :: cmask1
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: cn1
        real(kind=4), dimension(ip) , intent(Out) :: cn2l
        real(kind=4), dimension(ip) , intent(Out) :: cn2s
        real(kind=4), dimension(jp) , intent(Out) :: cn3l
        real(kind=4), dimension(jp) , intent(Out) :: cn3s
        real(kind=4), dimension(kp) , intent(Out) :: cn4l
        real(kind=4), dimension(kp) , intent(Out) :: cn4s
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(Out) :: dmask1
        real(kind=4), dimension(0:ip) , intent(In) :: dxs
        real(kind=4), dimension(0:jp) , intent(In) :: dys
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzn
        real(kind=4), dimension(-1:kp+2) , intent(In) :: dzs
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(Out) :: p
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(Out) :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(Out) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(Out) :: w
        real(kind=4), dimension(kp+2) , intent(In) :: z2
        real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1) , intent(InOut) :: zbm
! 
! WV: The original boundary was 0,km;0,jm;0,im. This does not init the boundary values,so I changed it to the dimensions of u,v,w,p
!      do k = 0,km
!      do j = 0,jm
!      do i = 0,im
!        u(i,j,k) = 0.0
!        v(i,j,k) = 0.0
!        w(i,j,k) = 0.0
!        p(i,j,k) = 0.0
!      end do
!      end do
!      end do
!print *, 'Zero arrays'
      do k = 0,kp+1
      do j = 0,jp+2
      do i = 0,ip+2
          p(i,j,k) = 0.0
      end do
      end do
      end do

      do k = 0,kp+1
      do j = -1,jp+1
      do i = 0,ip+1
        u(i,j,k) = 0.0
        v(i,j,k) = 0.0
      end do
      end do
      end do 

      do k = -1,kp+1
      do j = -1,jp+1
      do i = 0,ip+1
        w(i,j,k) =  0.0
      end do
      end do
      end do
!print *, 'call feedbfm'
! --check
#if IFBF == 1
      call feedbfm(km,jm,im,amask1,bmask1,cmask1,dmask1,zbm,z2,dzn)
!      if(ifbf == 1) call feedbfm(km,jm,im,amask1,bmask1,cmask1,dmask1,zbm,z2,dzn)
#endif
!print *, 'Parameter settings for solving Poisson equation'
! =====================================================
! 
!      Parameter settings for solving Poisson equation
! 
! =====================================================
      do i = 1,im
      cn2s(i) = 2./(dxs(i-1)*(dxs(i-1)+dxs(i)))
      cn2l(i) = 2./(dxs(i)*(dxs(i-1)+dxs(i)))
      end do
      do j = 1,jm
      cn3s(j) = 2./(dys(j-1)*(dys(j-1)+dys(j)))
      cn3l(j) = 2./(dys(j)*(dys(j-1)+dys(j)))
      end do
! 
      do k = 1,km
          dz1 = dzs(k-1)
          dz2 = dzs(k)
          cn4s(k) = 2./(dz1*(dz1+dz2))
          cn4l(k) = 2./(dz2*(dz1+dz2))
      do j = 1,jm
      do i = 1,im
          cn1(i,j,k) = &
      ! ((dxs(i)+dxs(i-1)) /(dxs(i-1)+dxs(i)))* 
      2./(dxs(i-1)*dxs(i))  + &
      ! ((dys(j)+dys(j-1))/(dys(j-1)+dys(j)))*
      2./(dys(j-1)*dys(j)) + &
      2./(dz1*dz2) ! *((dz1+dz2)/(dz1+dz2))
          cn1(i,j,k) = 1./cn1(i,j,k)
      end do
      end do
      end do
      return
      end subroutine init

end module module_init

