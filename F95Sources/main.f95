      program main
        use module_init 
        use module_grid 
        use module_set 
        use module_timdata 
#ifdef TIMSERIS_FIXED
        use module_timseris
#endif
        use module_aveflow
        use module_ifdata 
#if IANIME == 1      
        use module_anime 
#endif      
#ifdef _OPENCL_LES_WV
        use module_LES_combined_ocl
#else        
        use module_velnw 
        use module_bondv1 
        use module_velFG 
#if IFBF == 1      
        use module_feedbf 
#endif      
        use module_les 
        use module_press 
        use module_adam 
#endif        
        use common_sn 
        real(kind=4) :: alpha
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: amask1
        real(kind=4), dimension(ip,jp,kp)  :: avel
        real(kind=4), dimension(ip,jp,kp)  :: avep
        real(kind=4), dimension(ip,jp,kp)  :: avesm
        real(kind=4), dimension(ip,jp,kp)  :: avesmsm
        real(kind=4), dimension(ip,kp)  :: avesu
        real(kind=4), dimension(ip,kp)  :: avesuu
        real(kind=4), dimension(ip,kp)  :: avesv
        real(kind=4), dimension(ip,kp)  :: avesvv
        real(kind=4), dimension(ip,kp)  :: avesw
        real(kind=4), dimension(ip,kp)  :: avesww
        real(kind=4), dimension(ip,jp,kp)  :: aveu
        real(kind=4), dimension(ip,jp,kp)  :: aveuu
        real(kind=4), dimension(ip,jp,kp)  :: avev
        real(kind=4), dimension(ip,jp,kp)  :: avevv
        real(kind=4), dimension(ip,jp,kp)  :: avew
        real(kind=4), dimension(ip,jp,kp)  :: aveww
        real(kind=4) :: beta
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1)  :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: cmask1
        real(kind=4), dimension(ip,jp,kp)  :: cn1
        real(kind=4), dimension(ip)  :: cn2l
        real(kind=4), dimension(ip)  :: cn2s
        real(kind=4), dimension(jp)  :: cn3l
        real(kind=4), dimension(jp)  :: cn3s
        real(kind=4), dimension(kp)  :: cn4l
        real(kind=4), dimension(kp)  :: cn4s
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: cov1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov9
        character(len=70) :: data10
        character(len=70) :: data11
        character(len=70) :: data20
        character(len=70) :: data21
        character(len=70) :: data22
        character(len=70) :: data23
        character(len=70) :: data24
        character(len=70) :: data25
        character(len=70) :: data26
        character(len=70) :: data27
        character(len=70) :: data30
        character(len=70) :: data31
        real(kind=4), dimension(kp)  :: delx1
        real(kind=4), dimension(0:ip,jp,kp)  :: dfu1
        real(kind=4), dimension(ip,0:jp,kp)  :: dfv1
        real(kind=4), dimension(ip,jp,kp)  :: dfw1
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu9
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: dmask1
        real(kind=4) :: dt
        real(kind=4), dimension(-1:ip+1)  :: dx1
        real(kind=4), dimension(0:ip)  :: dxl
        real(kind=4), dimension(0:ip)  :: dxs
        real(kind=4), dimension(0:jp+1)  :: dy1
        real(kind=4), dimension(0:jp)  :: dyl
        real(kind=4), dimension(0:jp)  :: dys
        real(kind=4), dimension(-1:kp+2)  :: dzn
        real(kind=4), dimension(-1:kp+2)  :: dzs
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: f
#if ICAL == 1
        real(kind=4), dimension(ip,jp,kp)  :: fghold
#endif
        real(kind=4), dimension(ip,jp,kp)  :: fold
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fz
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: g
        real(kind=4), dimension(ip,jp,kp)  :: gold
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: h
        real(kind=4), dimension(ip,jp,kp)  :: hold
#ifndef _OPENCL_LES_WV
        real(kind=4), dimension(ip,jp,kp)  :: fghold
#endif
        integer :: ianime
        integer :: ical
        integer :: ifbf
        integer :: im
        integer :: jm
        integer :: km
        integer :: n
        integer :: n0
        integer :: n1
        integer :: nmax
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou9
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1)  :: p
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: rhs
        real(kind=4) :: ro
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1)  :: sm
        real(kind=4) :: time
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: usum
        real(kind=4), dimension(ip,jp,kp)  :: uwfx
        real(kind=4), dimension(ip,kp)  :: uwfxs
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: v
        real(kind=4) :: vn
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: wsum
        real(kind=4), dimension(kp+2)  :: z2
        !WV NOTE ipmax, jpmax are needed because the zbm array is read from a file and needs to be large enough
        !WV Otherwise the only supported domain is 150x150x90
        real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1)  :: zbm
#ifdef TIMINGS
        real (kind=4), dimension(0:9) :: timestamp
#endif
! -----------------------------------------------------------------------
!
      call set(data10,data11,data20,data21,data22,data23,data24,data25,data26,data27,data30,data31, &
      im,jm,km,ifbf,ianime,ical,n0,n1,nmax,dt,ro,vn,alpha,beta)
      call grid(dx1,dxl,dy1,dyl,z2,dzn,dzs,dxs,dys)
      call timdata()
      call init(km,jm,im,u,v,w,p,cn2s,dxs,cn2l,cn3s,dys,cn3l,dzs,cn4s,cn4l,cn1,amask1,bmask1, &
      cmask1,dmask1,zbm,z2,dzn)
      n=n0
      call ifdata( &
#if ICAL == 1
      data30,data31, fold,gold,hold,fghold, time &
#endif
      n,u,im,jm,km,v,w,p,usum,vsum,wsum, &
      delx1,dx1,dy1,dzn,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h,z2,dt,dxs,cov1, &
      cov2,cov3,dfu1,vn,cov4,cov5,cov6,dfv1,cov7,cov8,cov9,dfw1,dzs,nou1,nou5,nou9,nou2,nou3,nou4, &
      nou6,nou7,nou8,bmask1,cmask1,dmask1,alpha,beta,fx,fy,fz,amask1,zbm)

#ifdef _OPENCL_LES_WV
      call initialise_LES_kernel( &
            p,u,v,w,usum,vsum,wsum,f,g,h,fold,gold,hold, &
            diu1, diu2, diu3, diu4, diu5, diu6, diu7, diu8, diu9, &
            amask1, bmask1, cmask1, dmask1, &
            cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s, &
            rhs, sm, dxs, dys, dzs, dx1, dy1, dzn, z2, &
            dt, im, jm, km &
              )
#endif


#ifdef VERBOSE
#ifdef _OPENCL_LES_WV
    print *,'MAIN: calling OpenCL run_LES_kernel for ', nmax-n0+1, ' time steps, domain = ',im,'x',jm,'x',km
#else
    print *,'MAIN: running reference LES code for ', nmax-n0+1, ' time steps, domain = ',im,'x',jm,'x',km
#endif
#endif
! --main loop
#ifdef TIMINGS
    nmax=201
    call cpu_time(timestamp(8))
#endif
      do n = n0,nmax
        time = float(n-1)*dt
! -------calculate turbulent flow--------c
#ifdef _OPENCL_LES_WV
      call run_LES_kernel ( &
            n, nmax &
            )
#else
! -------calculate turbulent flow--------c
#ifdef TIMINGS
         print *, 'run_LES_reference: time step = ',n
#endif
        ! ========================================================================================================================================================
        ! ========================================================================================================================================================
#ifdef TIMINGS
        call cpu_time(timestamp(0))
#endif
        call velnw(km,jm,im,p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)
#ifdef TIMINGS
        call cpu_time(timestamp(1))
#endif
        call bondv1(jm,u,z2,dzn,v,w,km,n,im,dt,dxs)
#ifdef TIMINGS
        call cpu_time(timestamp(2))
#endif
        call velfg(km,jm,im,dx1,cov1,cov2,cov3,dfu1,diu1,diu2,dy1,diu3,dzn,vn,f,cov4,cov5,cov6,dfv1, &
      diu4,diu5,diu6,g,cov7,cov8,cov9,dfw1,diu7,diu8,diu9,dzs,h,nou1,u,nou5,v,nou9,w,nou2,nou3, &
      nou4,nou6,nou7,nou8)
#ifdef TIMINGS
        call cpu_time(timestamp(3))
#endif
#if IFBF == 1
        call feedbf(km,jm,im,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f,g, &
      h)
#endif
#ifdef TIMINGS
        call cpu_time(timestamp(4))
#endif
        call les(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h)
#ifdef TIMINGS
        call cpu_time(timestamp(5))
#endif
        call adam(n,nmax,data21,fold,im,jm,km,gold,hold,fghold,f,g,h)
#ifdef TIMINGS
        call cpu_time(timestamp(6))
#endif
        call press(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,f,g,h,dt,cn1,cn2l,p,cn2s,cn3l,cn3s,cn4l,cn4s,n, &
      nmax,data20,usum,vsum,wsum)
#ifdef TIMINGS
        call cpu_time(timestamp(7))
        do i=1, 7
            print '("Time for state ",i2," = ",f6.3," s")',i,timestamp(i)-timestamp(i-1)
        end do
#endif

#endif		
! -------data output ---------------------c
! WV: This is clearly broken, as the dimensions for u/v/w are 150x150x90
#ifdef TIMSERIS_FIXED
        call timseris(n,dt,u,v,w)
#endif
        call aveflow(n,n1,km,jm,im,aveu,avev,avew,avep,avel,aveuu,avevv,aveww,avesm,avesmsm,uwfx, &
      avesu,avesv,avesw,avesuu,avesvv,avesww,u,v,w,p,sm,nmax,uwfxs,data10,time,data11)
#if IANIME == 1
        call anime(n,n0,nmax,km,jm,im,dxl,dx1,dyl,dy1,z2,data22,data23,u,w,v,amask1)
#endif
!
      end do

#ifdef TIMINGS
    call cpu_time(timestamp(9))
    print *,"Total time:" ,timestamp(9)-timestamp(8),"s for ",nmax-n0,"iterations"
#endif
!

      end program



