module module_anime

contains

      subroutine anime(n,n0,nmax,km,jm,im,dxl,dx1,dyl,dy1,z2,data22,data23,u,w,v,amask1)
      use common_sn ! create_new_include_statements() line 102
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1) , intent(In) :: amask1
        character(len=70), intent(In) :: data22
        character(len=70), intent(In) :: data23
        real(kind=4), dimension(-1:ip+1) , intent(In) :: dx1
        real(kind=4), dimension(0:ip) , intent(In) :: dxl
        real(kind=4), dimension(0:jp+1) , intent(In) :: dy1
        real(kind=4), dimension(0:jp) , intent(In) :: dyl
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
        integer, intent(In) :: n
        integer, intent(In) :: n0
        integer, intent(In) :: nmax
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
        real(kind=4), dimension(kp+2) , intent(In) :: z2
! 
! 
      if(n == n0.or.n == nmax.or.mod(n,1000) == 0.) then
       do k = 1,km
       do j = 1,jm
       do i = 1,im
        a1(i,j,k) = real(dxl(i-1)+dx1(i))
        a2(i,j,k) = real(dyl(j-1)+dy1(j))
        a3(i,j,k) = real(z2(k))
       end do
       end do
       end do

       open(unit=22,file=data22,form='unformatted',status='unknown')
        write(22) im,jm,km
        write(22) (((real(a1(i,j,k)),i=1,im),j=1,jm),k=1,km), (((real(a3(i,j,k)),i=1,im),j=1,jm), &
      k=1,km), (((real(a2(i,j,k)),i=1,im),j=1,jm),k=1,km)
       close(unit=22)

       open(unit=23,file=data23,form='unformatted',status='unknown')
        write(23) im,jm,km,4
        write(23) (((real(0.5*(u(i-1,j,k)+u(i,j,k))), i=1,im),j=1,jm),k=1,km), (((real(0.5*(w(i,j, &
      k-1)+w(i,j,k))), i=1,im),j=1,jm),k=1,km), (((real(0.5*(v(i,j-1,k)+v(i,j,k))), i=1,im),j=1, &
      jm),k=1,km), (((real(amask1(i,j,k)),i=1,im),j=1,jm),k=1,km)
       close(unit=23)
      endif      
! 
      return
      end subroutine anime




end module module_anime

