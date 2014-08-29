module module_aveflow

contains

      subroutine aveflow(n,n1,km,jm,im,aveu,avev,avew,avep,avel,aveuu,avevv,aveww,avesm,avesmsm, &
      uwfx,avesu,avesv,avesw,avesuu,avesvv,avesww,u,v,w,p,sm,nmax,uwfxs,data10,time,data11)
      use common_sn ! create_new_include_statements() line 102
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avel
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avep
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avesm
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avesmsm
        real(kind=4), dimension(ip,kp) , intent(Out) :: avesu
        real(kind=4), dimension(ip,kp) , intent(Out) :: avesuu
        real(kind=4), dimension(ip,kp) , intent(Out) :: avesv
        real(kind=4), dimension(ip,kp) , intent(Out) :: avesvv
        real(kind=4), dimension(ip,kp) , intent(Out) :: avesw
        real(kind=4), dimension(ip,kp) , intent(Out) :: avesww
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: aveu
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: aveuu
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avev
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avevv
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: avew
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: aveww
        character(len=70), intent(In) :: data10
        character(len=70), intent(In) :: data11
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
        integer, intent(In) :: n
        integer, intent(In) :: n1
        integer, intent(In) :: nmax
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1) , intent(In) :: p
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1) , intent(In) :: sm
        real(kind=4), intent(In) :: time
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(ip,jp,kp) , intent(Out) :: uwfx
        real(kind=4), dimension(ip,kp) , intent(InOut) :: uwfxs
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
! 
! 
      if(n == n1) then
        do k = 1,km
        do j = 1,jm
        do i = 1,im
          aveu(i,j,k) = 0.0
          avev(i,j,k) = 0.0
          avew(i,j,k) = 0.0
          avep(i,j,k) = 0.0
          avel(i,j,k) = 0.0
          aveuu(i,j,k) = 0.0
          avevv(i,j,k) = 0.0
          aveww(i,j,k) = 0.0
          avesm(i,j,k) = 0.0
          avesmsm(i,j,k) = 0.0
          uwfx(i,j,k) = 0.0
        end do
        end do
        end do
        do k = 1,km
        do i = 1,im
          avesu(i,k) = 0.0
          avesv(i,k) = 0.0
          avesw(i,k) = 0.0
          avesuu(i,k) = 0.0
          avesvv(i,k) = 0.0
          avesww(i,k) = 0.0
        end do
        end do 
      end if
! 
      if(n >= n1) then
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        aveu(i,j,k) = aveu(i,j,k)+u(i,j,k)
        avev(i,j,k) = avev(i,j,k)+v(i,j,k)
        avew(i,j,k) = avew(i,j,k)+w(i,j,k)
        avep(i,j,k) = avep(i,j,k)+p(i,j,k)
        aveuu(i,j,k) = aveuu(i,j,k)+u(i,j,k)**2
        avevv(i,j,k) = avevv(i,j,k)+v(i,j,k)**2
        aveww(i,j,k) = aveww(i,j,k)+w(i,j,k)**2
        avesm(i,j,k) = avesm(i,j,k)+sm(i,j,k)
        avesmsm(i,j,k) = avesmsm(i,j,k)+sm(i,j,k)**2
        uwfx(i,j,k) = uwfx(i,j,k)+0.5*(u(i,j,k-1)+u(i,j,k)) *0.5*(w(i,j,k-1)+w(i+1,j,k-1))
      end do
      end do
      end do  
! 
      endif
! 
      if(n == nmax) then
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        aveu(i,j,k) = aveu(i,j,k)/float(nmax-n1+1)
        avev(i,j,k) = avev(i,j,k)/float(nmax-n1+1)
        avew(i,j,k) = avew(i,j,k)/float(nmax-n1+1)
        avep(i,j,k) = avep(i,j,k)/float(nmax-n1+1)
        avel(i,j,k) = avel(i,j,k)/float(nmax-n1+1)
        aveuu(i,j,k) = aveuu(i,j,k)/float(nmax-n1+1)
        avevv(i,j,k) = avevv(i,j,k)/float(nmax-n1+1)
        aveww(i,j,k) = aveww(i,j,k)/float(nmax-n1+1)
        avesm(i,j,k) = avesm(i,j,k)/float(nmax-n1+1)
        avesmsm(i,j,k) = avesmsm(i,j,k)/float(nmax-n1+1)
      end do
      end do
      end do 
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        uwfx(i,j,k) = uwfx(i,j,k)/float(nmax-n1+1) -0.5*(aveu(i,j,k-1)+aveu(i,j,k)) *0.5*(avew(i,j, &
      k-1)+avew(i+1,j,k-1))
      end do
      end do
      end do
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        avesu(i,k) = avesu(i,k)+aveu(i,j,k)
        avesv(i,k) = avesv(i,k)+avev(i,j,k)
        avesw(i,k) = avesw(i,k)+avew(i,j,k)
        avesuu(i,k) = avesuu(i,k)+aveuu(i,j,k)
        avesvv(i,k) = avesvv(i,k)+avevv(i,j,k)
        avesww(i,k) = avesww(i,k)+aveww(i,j,k)
        uwfxs(i,k) = uwfxs(i,k)+uwfx(i,j,k)
      end do
      end do
      end do        
      do k = 1,km
      do i = 1,im
        avesu(i,k) = avesu(i,k)/float(jm)
        avesv(i,k) = avesv(i,k)/float(jm)
        avesw(i,k) = avesw(i,k)/float(jm)
        avesuu(i,k) = avesuu(i,k)/float(jm)
        avesvv(i,k) = avesvv(i,k)/float(jm)
        avesww(i,k) = avesww(i,k)/float(jm)
        uwfxs(i,k) = uwfxs(i,k)/float(jm)
      end do
      end do
! 
        do k = 1,km
        do j = 1,jm
        do i = 1,im
          aveuu(i,j,k) = sqrt(abs(aveuu(i,j,k)-aveu(i,j,k)**2))
          avevv(i,j,k) = sqrt(abs(avevv(i,j,k)-avev(i,j,k)**2))
          aveww(i,j,k) = sqrt(abs(aveww(i,j,k)-avew(i,j,k)**2))
        end do
        end do
        end do
        do k = 1,km
        do i = 1,im
          avesuu(i,k) = sqrt(abs(avesuu(i,k)-avesu(i,k)**2))
          avesvv(i,k) = sqrt(abs(avesvv(i,k)-avesv(i,k)**2))
          avesww(i,k) = sqrt(abs(avesww(i,k)-avesw(i,k)**2))
        end do
        end do
! 
#ifndef NO_IO
        open(unit=10,file=data10,form='unformatted',status='unknown')
          write(10) n,time
          write(10) (((aveu(i,j,k),i=1,im),j=1,jm),k=1,km)
          write(10) (((avew(i,j,k),i=1,im),j=1,jm),k=1,km)
          write(10) (((avev(i,j,k),i=1,im),j=1,jm),k=1,km)
        close(unit=10)
! 
        open(unit=11,file=data11,form='unformatted',status='unknown')
          write(11) n,time
          write(11) (((aveuu(i,j,k),i=1,im),j=1,jm),k=1,km)
          write(11) (((aveww(i,j,k),i=1,im),j=1,jm),k=1,km)
          write(11) (((avevv(i,j,k),i=1,im),j=1,jm),k=1,km)
          write(11) (((uwfx(i,j,k),i=1,im),j=1,jm),k=1,km)
        close(unit=11)
#endif      
      endif
! 

      return
      end subroutine aveflow

end module module_aveflow

