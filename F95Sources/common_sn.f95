module common_sn
      use params_common_sn ! context_free_refactorings() line 275
!      An IMPLICIT statement specifies a type and size for all user-defined names
!      that begin with any letter, either a single letter or in a range of letters,
!      appearing in the specification. 

      implicit real*4(a-h,o-z)
      implicit integer(i-n) 
!! Original line !!       real*4 nou1,nou2,nou3,nou4,nou5,nou6,nou7,nou8,nou9
! 
!       parameter(ip=150,jp=150,kp=90)
! 
!! Original line !!       character*70 data10,data11,data20,data21,data22,data23 ,data24,data25,data26,data27,data30 ,data31,data41 ,data40,data12,data13,data19,data29 ,data50,data51,data52,data53,data54 ,fname
! 
! 
! 
! 
! 
! 
! 
! --les
! --ifdata
! --press
! --vel2,velFG
! --stretch
! 
      real a1(1:ip,1:jp+1,1:kp+1),a2(1:ip,1:jp+1,1:kp+1) ,a3(1:ip,1:jp+1,1:kp+1)

end module common_sn
