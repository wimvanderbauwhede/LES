module module_timdata

contains

      subroutine timdata()
      use common_sn ! create_new_include_statements() line 102
! 
! 
! -----time series of wind and concentration ----------       
      open(50,file='winddata.dat')
!               
      return
      end subroutine timdata

end module module_timdata

