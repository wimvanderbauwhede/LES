module module_timseris

contains

      subroutine timseris(n,dt,u,v,w)
      use common_sn
        real(kind=4), intent(In) :: dt
        integer, intent(In) :: n
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w

! 
! 
! -------jikeiretsu data-------------
      if(mod(n,2) == 0) then
      write(50,7000) float(n)*dt ,u(170, 65,4),u(170, 65,9),u(170, 65,14) ,u(198,100,4),u(198,100, &
      9),u(198,100,14) ,u(238, 95,4),u(238,100,9),u(238,100,14) ,u(248, 80,4),u(238, 80,9),u(238, &
       80,14) ,v(170, 65,4),v(170, 65,9),v(170, 65,14) ,v(198,100,4),v(198,100,9),v(198,100,14) , &
      v(238, 95,4),v(238,100,9),v(238,100,14) ,v(248, 80,4),v(238, 80,9),v(238, 80,14) ,w(170, 65, &
      4),w(170, 65,9),w(170, 65,14) ,w(198,100,4),w(198,100,9),w(198,100,14) ,w(238, 95,4),w(238, &
      100,9),w(238,100,14) ,w(248, 80,4),w(238, 80,9),w(238, 80,14)
 7000 format(37e20.10)
      endif 
!

      return
      end subroutine timseris

end module module_timseris

