module module_bondFG

contains

      subroutine bondfg(km,jm,f,im,g,h)
      use common_sn ! create_new_include_statements() line 102
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(Out) :: h
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
! 
! --inflow condition
      do k = 1,km
      do j = 1,jm
        f( 0,j,k) = f(1  ,j,k)
      end do
      end do
! --sideflow condition
      do k = 1,km
      do i = 1,im
        g(i, 0,k) = g(i,jm  ,k)
      end do
      end do
! --ground and top condition
      do j = 1,jm
      do i = 1,im
        h(i,j, 0) = 0.0
        h(i,j,km) = 0.0
      end do
      end do
!                                               
      return                                   
      end subroutine bondFG                                    

end module module_bondFG

