module params_common_sn
    integer, parameter  :: ipmax = 150;       integer, parameter  :: jpmax = 150
    !integer, parameter  :: ipmax = 254;       integer, parameter  :: jpmax = 253
#ifndef TEST_SMALL_DOMAIN
      integer, parameter  :: ip = 150;       integer, parameter  :: jp = 150;       integer, parameter  :: kp = 90
      !integer, parameter  :: ip = 254;       integer, parameter  :: jp = 253;       integer, parameter  :: kp = 94
#else      
      integer, parameter  :: ip = 25;       integer, parameter  :: jp = 25;       integer, parameter  :: kp = 90
#endif      
end module params_common_sn
