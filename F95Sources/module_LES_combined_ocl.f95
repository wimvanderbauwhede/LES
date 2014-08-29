! This is a cleaned-up version of the module template, without debug statements.
! For the full version see Revision: 7743
!#define FIXED_NTH
!#define APPROX_PAR_SOR
module module_LES_combined_ocl
    use module_LES_conversions
!    use module_LES_tests
    integer :: init_write_uvw_p_uvwsum_to_file = 0
    integer :: init_initialise_LES_kernel = 0
    integer :: init_write_uvw_p_to_file = 0
    integer :: init_write_fgh_old_to_file = 0
    integer :: init_run_LES_kernel = 0
contains
    subroutine initialise_LES_kernel (         p,u,v,w,usum,vsum,wsum,f,g,h,fold,gold,hold,         diu1, diu2, diu3, diu4, diu5, diu6, diu7, diu8, diu9,         amask1, bmask1, cmask1, dmask1,         cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s,         rhs, sm, dxs, dys, dzs, dx1, dy1, dzn, z2,         dt, im, jm, km         )
        use oclWrapper
        use params_common_sn
        implicit none
        ! Arguments
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1), intent(InOut)  :: p
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: usum
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: v
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(InOut) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: wsum
        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut)  :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut)  :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut)  :: h
        real(kind=4), dimension(ip,jp,kp), intent(InOut)  :: fold
        real(kind=4), dimension(ip,jp,kp), intent(InOut)  :: gold
        real(kind=4), dimension(ip,jp,kp), intent(InOut)  :: hold
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu9
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In)  :: amask1
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1), intent(In)  :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(In)  :: cmask1
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In)  :: dmask1
        real(kind=4), dimension(ip,jp,kp), intent(In)  :: cn1
        real(kind=4), dimension(ip), intent(In)  :: cn2l
        real(kind=4), dimension(ip), intent(In)  :: cn2s
        real(kind=4), dimension(jp), intent(In)  :: cn3l
        real(kind=4), dimension(jp), intent(In)  :: cn3s
        real(kind=4), dimension(kp), intent(In)  :: cn4l
        real(kind=4), dimension(kp), intent(In)  :: cn4s
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In)  :: rhs
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1), intent(In)  :: sm
        real(kind=4), dimension(0:ip), intent(In)  :: dxs
        real(kind=4), dimension(0:jp), intent(In)  :: dys
        real(kind=4), dimension(-1:kp+2), intent(In)  :: dzs
        real(kind=4), dimension(-1:ip+1), intent(In)  :: dx1
        real(kind=4), dimension(0:jp+1), intent(In)  :: dy1
        real(kind=4), dimension(-1:kp+2), intent(In)  :: dzn
        real(kind=4), dimension(kp+2)  :: z2        
        real(kind=4), intent(In) :: dt
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
!        integer, intent(In) :: nmax
        ! -----------------------------------------------------------------------
        ! Combined arrays for OpenCL kernels
        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1)  :: uvw
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp) :: uvwsum
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp)  :: fgh
        real(kind=4), dimension(0:3,ip,jp,kp)  :: fgh_old
!        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: cov
        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: diu
!        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: nou
        real(kind=4), dimension(0:3,-1:ip+1,-1:jp+1,0:kp+1)  :: mask1
        ! Extra arguments for kernel
        real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1) :: p_scratch
        real(kind=4), dimension(256) :: chunks_num, chunks_denom
        real(kind=4), dimension(256) :: val_ptr
        integer, dimension(256) :: n_ptr
        integer, dimension(256) :: state_ptr
        
        ! OpenCL-specific declarations
        integer :: initialise_LES_kernel_nunits, initialise_LES_kernel_nthreads, initialise_LES_kernel_globalrange, initialise_LES_kernel_localrange
        integer :: init_initialise_LES_kernel_local
        ! integer :: srclen, klen
        character(len=*), parameter :: srcstr="../OpenCL/Kernels/LES_combined_kernel_mono.cl"
        character(len=*), parameter :: kstr="LES_combined_kernel_mono"
        ! character(len=*), parameter :: koptsstr=""
        !This is a hook to insert the actual buffer declarations!
        
        ! OpenCL buffer declarations
        integer(8) :: cn1_buf
        integer(8) :: cn2l_buf
        integer(8) :: cn2s_buf
        integer(8) :: cn3l_buf
        integer(8) :: cn3s_buf
        integer(8) :: cn4l_buf
        integer(8) :: cn4s_buf
        integer(8) :: rhs_buf
        integer(8) :: sm_buf
        integer(8) :: dxs_buf
        integer(8) :: dys_buf
        integer(8) :: dzs_buf
        integer(8) :: dx1_buf
        integer(8) :: dy1_buf
        integer(8) :: dzn_buf
        integer(8) :: z2_buf
        integer(8) :: uvw_buf
        integer(8) :: uvwsum_buf
        integer(8) :: fgh_buf
        integer(8) :: fgh_old_buf
        integer(8) :: diu_buf
        integer(8) :: mask1_buf
        integer(8) :: p_scratch_buf
        integer(8) :: chunks_num_buf
        integer(8) :: chunks_denom_buf
        integer(8) :: val_ptr_buf
        integer(8) :: n_ptr_buf
        integer(8) :: state_ptr_buf
        ! OpenCL buffer size declarations
        integer, dimension(3):: cn1_sz
        integer, dimension(1):: cn2l_sz
        integer, dimension(1):: cn2s_sz
        integer, dimension(1):: cn3l_sz
        integer, dimension(1):: cn3s_sz
        integer, dimension(1):: cn4l_sz
        integer, dimension(1):: cn4s_sz
        integer, dimension(3):: rhs_sz
        integer, dimension(3):: sm_sz
        integer, dimension(1):: dxs_sz
        integer, dimension(1):: dys_sz
        integer, dimension(1):: dzs_sz
        integer, dimension(1):: dx1_sz
        integer, dimension(1):: dy1_sz
        integer, dimension(1):: dzn_sz
        integer, dimension(1):: z2_sz
        integer, dimension(4):: uvw_sz
        integer, dimension(4):: uvwsum_sz
        integer, dimension(4):: fgh_sz
        integer, dimension(4):: fgh_old_sz
        integer, dimension(4):: diu_sz
        integer, dimension(4):: mask1_sz
        integer, dimension(4):: p_scratch_sz
        integer, dimension(1):: chunks_num_sz
        integer, dimension(1):: chunks_denom_sz
        integer, dimension(1):: val_ptr_sz
        integer, dimension(1):: n_ptr_sz
        integer, dimension(1):: state_ptr_sz
    integer :: init_ocl_local
    init_ocl_local = 0
        
        ! Convert to new format
        call convert_to_uvw(u,v,w,uvw)
        call convert_to_fgh(usum,vsum,wsum,uvwsum)
        call convert_to_fgh(f,g,h,fgh)
        call convert_to_fgh_old(fold,gold,hold, fgh_old)
        ! The following are all read-only
        call convert_to_9vec(diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,diu)
        call convert_to_mask1(amask1,bmask1,cmask1,dmask1,mask1)
        state_ptr(1) = 0
        p_scratch(0,:,:,:) = p(:,:,:)
        p_scratch(1,:,:,:) = p(:,:,:)
        chunks_num=0.0
        chunks_denom=0.0
        
        if ( init_ocl_local /= 1 ) then 
          init_initialise_LES_kernel_local = 1
        if ( init_initialise_LES_kernel /= 1 ) then 
          init_initialise_LES_kernel = 1
        !  call oclInitOpts(srcstr,kstr,koptsstr)
          call oclInit(srcstr,kstr)
        
          call oclGetMaxComputeUnits(initialise_LES_kernel_nunits)
          call oclGetNThreadsHint(initialise_LES_kernel_nthreads)
        
        !  print *, "Compute units:",initialise_LES_kernel_nunits," Threads:",initialise_LES_kernel_nthreads
        end if
        
        ! OpenCL buffer sizes
        cn1_sz = shape(cn1)
        cn2l_sz = shape(cn2l)
        cn2s_sz = shape(cn2s)
        cn3l_sz = shape(cn3l)
        cn3s_sz = shape(cn3s)
        cn4l_sz = shape(cn4l)
        cn4s_sz = shape(cn4s)
        rhs_sz = shape(rhs)
        sm_sz = shape(sm)
        dxs_sz = shape(dxs)
        dys_sz = shape(dys)
        dzs_sz = shape(dzs)
        dx1_sz = shape(dx1)
        dy1_sz = shape(dy1)
        dzn_sz = shape(dzn)
        z2_sz = shape(z2)
        uvw_sz = shape(uvw)
        uvwsum_sz = shape(uvwsum)
        fgh_sz = shape(fgh)
        fgh_old_sz = shape(fgh_old)
        diu_sz = shape(diu)
        mask1_sz = shape(mask1)
        p_scratch_sz = shape(p_scratch)
        chunks_num_sz = shape(chunks_num)
        chunks_denom_sz = shape(chunks_denom)
        val_ptr_sz = shape(val_ptr)
        n_ptr_sz = shape(n_ptr)
        state_ptr_sz = shape(state_ptr)
        
        ! Create OpenCL buffers
        call oclMake3DFloatArrayReadBuffer(cn1_buf,cn1_sz ,cn1)
        call oclMake1DFloatArrayReadBuffer(cn2l_buf,cn2l_sz ,cn2l)
        call oclMake1DFloatArrayReadBuffer(cn2s_buf,cn2s_sz ,cn2s)
        call oclMake1DFloatArrayReadBuffer(cn3l_buf,cn3l_sz ,cn3l)
        call oclMake1DFloatArrayReadBuffer(cn3s_buf,cn3s_sz ,cn3s)
        call oclMake1DFloatArrayReadBuffer(cn4l_buf,cn4l_sz ,cn4l)
        call oclMake1DFloatArrayReadBuffer(cn4s_buf,cn4s_sz ,cn4s)
        call oclMake3DFloatArrayReadBuffer(rhs_buf,rhs_sz ,rhs)
        call oclMake3DFloatArrayReadBuffer(sm_buf,sm_sz ,sm)
        call oclMake1DFloatArrayReadBuffer(dxs_buf,dxs_sz ,dxs)
        call oclMake1DFloatArrayReadBuffer(dys_buf,dys_sz ,dys)
        call oclMake1DFloatArrayReadBuffer(dzs_buf,dzs_sz ,dzs)
        call oclMake1DFloatArrayReadBuffer(dx1_buf,dx1_sz ,dx1)
        call oclMake1DFloatArrayReadBuffer(dy1_buf,dy1_sz ,dy1)
        call oclMake1DFloatArrayReadBuffer(dzn_buf,dzn_sz ,dzn)
        call oclMake1DFloatArrayReadBuffer(z2_buf,z2_sz ,z2)
        call oclMake4DFloatArrayReadWriteBuffer(uvw_buf,uvw_sz ,uvw)
        call oclMake4DFloatArrayReadWriteBuffer(uvwsum_buf,uvwsum_sz ,uvwsum)
        call oclMake4DFloatArrayReadWriteBuffer(fgh_buf,fgh_sz ,fgh)
        call oclMake4DFloatArrayReadWriteBuffer(fgh_old_buf,fgh_old_sz ,fgh_old)
        call oclMake4DFloatArrayReadWriteBuffer(diu_buf,diu_sz ,diu)
        call oclMake4DFloatArrayReadWriteBuffer(mask1_buf,mask1_sz ,mask1)
        call oclMake4DFloatArrayReadBuffer(p_scratch_buf,p_scratch_sz ,p_scratch)
        call oclMake1DFloatArrayReadWriteBuffer(chunks_num_buf,chunks_num_sz ,chunks_num)
        call oclMake1DFloatArrayReadWriteBuffer(chunks_denom_buf,chunks_denom_sz ,chunks_denom)
        call oclMake1DFloatArrayReadWriteBuffer(val_ptr_buf,val_ptr_sz ,val_ptr)
        call oclMake1DIntArrayReadWriteBuffer(n_ptr_buf,n_ptr_sz ,n_ptr)
        call oclMake1DIntArrayReadWriteBuffer(state_ptr_buf,state_ptr_sz ,state_ptr)
        
        ! Set OpenCL argument order
        call oclSetFloatArrayArg(0, p_scratch_buf )
        call oclSetFloatArrayArg(1, uvw_buf )
        call oclSetFloatArrayArg(2, uvwsum_buf )
        call oclSetFloatArrayArg(3, fgh_buf )
        call oclSetFloatArrayArg(4, fgh_old_buf )
        call oclSetFloatArrayArg(5, rhs_buf )
        call oclSetFloatArrayArg(6, mask1_buf )
        call oclSetFloatArrayArg(7, diu_buf )
        call oclSetFloatArrayArg(8, sm_buf )
        call oclSetFloatArrayArg(9, dxs_buf )
        call oclSetFloatArrayArg(10, dys_buf )
        call oclSetFloatArrayArg(11, dzs_buf )
        call oclSetFloatArrayArg(12, dx1_buf )
        call oclSetFloatArrayArg(13, dy1_buf )
        call oclSetFloatArrayArg(14, dzn_buf )
        call oclSetFloatArrayArg(15, z2_buf )
        call oclSetFloatArrayArg(16, cn1_buf )
        call oclSetFloatArrayArg(17, cn2l_buf )
        call oclSetFloatArrayArg(18, cn2s_buf )
        call oclSetFloatArrayArg(19, cn3l_buf )
        call oclSetFloatArrayArg(20, cn3s_buf )
        call oclSetFloatArrayArg(21, cn4l_buf )
        call oclSetFloatArrayArg(22, cn4s_buf )
        call oclSetFloatArrayArg(23, val_ptr_buf )
        call oclSetFloatArrayArg(24, chunks_num_buf )
        call oclSetFloatArrayArg(25, chunks_denom_buf )
        call oclSetIntArrayArg(26, n_ptr_buf )
        call oclSetIntArrayArg(27, state_ptr_buf )
        call oclSetFloatConstArg(28, dt )
        call oclSetIntConstArg(29, im )
        call oclSetIntConstArg(30, jm )
        call oclSetIntConstArg(31, km )
        
        end if
        
        initialise_LES_kernel_globalrange = 1
        initialise_LES_kernel_localrange = 1
        if (initialise_LES_kernel_localrange == 0) then
          call padRange(initialise_LES_kernel_globalrange,initialise_LES_kernel_nunits*initialise_LES_kernel_nthreads)
        end if
        
        ! Copy all arrays required for the full run
        call oclWrite3DFloatArrayBuffer(cn1_buf,cn1_sz,cn1)
        call oclWrite1DFloatArrayBuffer(cn2l_buf,cn2l_sz,cn2l)
        call oclWrite1DFloatArrayBuffer(cn2s_buf,cn2s_sz,cn2s)
        call oclWrite1DFloatArrayBuffer(cn3l_buf,cn3l_sz,cn3l)
        call oclWrite1DFloatArrayBuffer(cn3s_buf,cn3s_sz,cn3s)
        call oclWrite1DFloatArrayBuffer(cn4l_buf,cn4l_sz,cn4l)
        call oclWrite1DFloatArrayBuffer(cn4s_buf,cn4s_sz,cn4s)
        call oclWrite3DFloatArrayBuffer(rhs_buf,rhs_sz,rhs)
        call oclWrite3DFloatArrayBuffer(sm_buf,sm_sz,sm)
        call oclWrite1DFloatArrayBuffer(dxs_buf,dxs_sz,dxs)
        call oclWrite1DFloatArrayBuffer(dys_buf,dys_sz,dys)
        call oclWrite1DFloatArrayBuffer(dzs_buf,dzs_sz,dzs)
        call oclWrite1DFloatArrayBuffer(dx1_buf,dx1_sz,dx1)
        call oclWrite1DFloatArrayBuffer(dy1_buf,dy1_sz,dy1)
        call oclWrite1DFloatArrayBuffer(dzn_buf,dzn_sz,dzn)
        call oclWrite1DFloatArrayBuffer(z2_buf,z2_sz,z2)
        call oclWrite4DFloatArrayBuffer(uvw_buf,uvw_sz,uvw)
        call oclWrite4DFloatArrayBuffer(uvwsum_buf,uvwsum_sz,uvwsum)
        call oclWrite4DFloatArrayBuffer(fgh_buf,fgh_sz,fgh)
        call oclWrite4DFloatArrayBuffer(fgh_old_buf,fgh_old_sz,fgh_old)
        call oclWrite4DFloatArrayBuffer(diu_buf,diu_sz,diu)
        call oclWrite4DFloatArrayBuffer(mask1_buf,mask1_sz,mask1)
        call oclWrite4DFloatArrayBuffer(p_scratch_buf,p_scratch_sz,p_scratch)
        call oclWrite1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
        call oclWrite1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz,chunks_denom)
        call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz,val_ptr)
        call oclWrite1DIntArrayBuffer(n_ptr_buf,n_ptr_sz,n_ptr)
        call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        
        ! call LES_combined(p_scratch, uvw, uvwsum, fgh, fgh_old, rhs, mask1, diu, sm,            dxs, dys, dzs, dx1, dy1, dzn, z2, cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s,            val_ptr, chunks_num, chunks_denom, n_ptr, state_ptr, dt, im, jm, km             )
        call runOcl(initialise_LES_kernel_globalrange,initialise_LES_kernel_localrange)
        
        ! Read back Read and ReadWrite arrays
        call oclRead4DFloatArrayBuffer(uvw_buf,uvw_sz,uvw)
        call oclRead4DFloatArrayBuffer(uvwsum_buf,uvwsum_sz,uvwsum)
        call oclRead4DFloatArrayBuffer(fgh_buf,fgh_sz,fgh)
        call oclRead4DFloatArrayBuffer(fgh_old_buf,fgh_old_sz,fgh_old)
        call oclRead4DFloatArrayBuffer(diu_buf,diu_sz,diu)
        call oclRead4DFloatArrayBuffer(mask1_buf,mask1_sz,mask1)
        call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
        call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz,chunks_denom)
        call oclRead1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz,val_ptr)
        call oclRead1DIntArrayBuffer(n_ptr_buf,n_ptr_sz,n_ptr)
        call oclRead1DIntArrayBuffer(state_ptr_buf,state_ptr_sz,state_ptr)
        
        ! Following buffers are used in the loop, assign to module-level buffer array for convenience
        oclBuffers(1) = p_scratch_buf ! BOUNDP
        oclBuffers(2) = uvw_buf ! BOUNDP
        oclBuffers(3) = uvwsum_buf ! BOUNDP
        oclBuffers(4) = fgh_buf ! DEBUG
        oclBuffers(5) = fgh_old_buf ! ADAM
        oclBuffers(6) = val_ptr_buf ! SOR, ADJ
        oclBuffers(7) = chunks_num_buf ! RHSAV, SOR, PAV
        oclBuffers(8) = chunks_denom_buf ! RHSAV, PAV
        oclBuffers(9) = n_ptr_buf ! BONDV1, SOR
        oclBuffers(10) = state_ptr_buf ! ALL
        oclNunits = initialise_LES_kernel_nunits
        oclNthreadsHint = initialise_LES_kernel_nthreads
    end subroutine initialise_LES_kernel 
    ! --------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------






    subroutine run_LES_kernel (         n, nmax         )
        use oclWrapper
        use oclWrapper
        use params_common_sn
        implicit none
!        integer, intent(In) :: im
!        integer, intent(In) :: jm
!        integer, intent(In) :: km
        integer, intent(In) :: n, nmax
        
        ! -----------------------------------------------------------------------
        ! arrays for OpenCL kernels
!        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1), intent(In)  :: p
!        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
!        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
!        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In):: f
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: g
!        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In)  :: h
!        real(kind=4), dimension(ip,jp,kp) , intent(In) :: fold
!        real(kind=4), dimension(ip,jp,kp), intent(In) :: gold
!        real(kind=4), dimension(ip,jp,kp) , intent(In) :: hold
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: usumo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: vsumo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: wsumo
        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1)  :: uvw
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp) :: uvwsum
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp)  :: fgh
        real(kind=4), dimension(0:3,ip,jp,kp)  :: fgh_old
        real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1)  :: po
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: uo
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: vo
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: wo
!        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fo
!        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: go
!        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: ho
        real(kind=4), dimension(ip,jp,kp)  :: foldo
        real(kind=4), dimension(ip,jp,kp) :: goldo
        real(kind=4), dimension(ip,jp,kp)  :: holdo
        real(kind=4), dimension(256) :: val_ptr
        real(kind=4), dimension(256) :: chunks_num, chunks_denom
        integer, dimension(256) :: n_ptr
        integer, dimension(256) :: state_ptr
        
        integer(8) :: p_buf
        integer(8) :: uvw_buf
        integer(8) :: uvwsum_buf
        integer(8) :: fgh_old_buf, fgh_buf
        integer(8) :: chunks_num_buf, chunks_denom_buf
        integer(8) :: val_ptr_buf, n_ptr_buf, state_ptr_buf
        integer, dimension(4) :: p_sz, uvw_sz, uvwsum_sz, fgh_old_sz,fgh_sz
        integer, dimension(1) :: state_ptr_sz, n_ptr_sz, val_ptr_sz, chunks_num_sz, chunks_denom_sz
        integer :: state, nn, ii
        real(kind=4) :: aaa, bbb, rhsav, area, pav, pco, sor
        integer, parameter  :: nmaxp = 50
        ! --only used mac method
        real, parameter  :: pjuge = 0.0001
        integer :: ngroups, nrd, iter
        real(kind=4) :: avg_iter
        integer, parameter :: ST_INIT=0, ST_VELNW__BONDV1_INIT_UVW=1, ST_BONDV1_CALC_UOUT=2
        integer, parameter :: ST_BONDV1_CALC_UVW=3, ST_VELFG__FEEDBF__LES_CALC_SM=4, ST_LES_BOUND_SM=5, ST_LES_CALC_VISC__ADAM=6
        integer, parameter :: ST_BONDV1_CALC_UVW__VELFG__FEEDBF__LES_CALC_SM=30, ST_VELFG=31, ST_FEEDBF__LES_CALC_SM=32
        integer, parameter :: ST_PRESS_RHSAV=7, ST_PRESS_SOR=8, ST_PRESS_PAV=9, ST_PRESS_ADJ=10, ST_PRESS_BOUNDP=11, ST_DONE=12
        foldo=0
        goldo=0
        holdo=0
        p_buf = oclBuffers(1) 
        uvw_buf = oclBuffers(2)
        uvwsum_buf = oclBuffers(3)
        fgh_buf = oclBuffers(4)
        fgh_old_buf = oclBuffers(5)
        val_ptr_buf=oclBuffers(6) ! SOR, ADJ
        chunks_num_buf = oclBuffers(7) ! RHSAV, SOR, PAV
        chunks_denom_buf = oclBuffers(8) ! RHSAV, PAV
        n_ptr_buf = oclBuffers(9) ! BONDV1, SOR
        state_ptr_buf = oclBuffers(10) ! ALL
        p_sz = shape(po)
        uvw_sz = shape(uvw)
        uvwsum_sz = shape(uvwsum)
        fgh_sz = shape(fgh)
        fgh_old_sz = shape(fgh_old)
        val_ptr_sz = shape(val_ptr)
        chunks_num_sz = shape(chunks_num)
        chunks_denom_sz = shape(chunks_denom)
        n_ptr_sz = shape(n_ptr)
        state_ptr_sz = shape(state_ptr)
        n_ptr(1)=n
        ! ========================================================================================================================================================
        ! ========================================================================================================================================================
        ! 2. Run the time/state nested loops, copying only time and state
        do state = ST_VELNW__BONDV1_INIT_UVW, ST_PRESS_BOUNDP
            state_ptr(1)=state
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz, state_ptr)
            select case (state)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_VELNW__BONDV1_INIT_UVW)
                    oclGlobalRange=(ip+1)*jp*kp
                    oclLocalRange=0
                    call oclWrite1DIntArrayBuffer(n_ptr_buf,n_ptr_sz, n_ptr)
                    call runOcl(oclGlobalRange,oclLocalRange)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_BONDV1_CALC_UOUT) ! REDUCTION
                    oclGlobalRange = oclNunits*128 ! OclNunits and then loop over j, at least for Mac
                    oclLocalRange = 128
                    ngroups = oclNunits
                    call runOcl(oclGlobalRange,oclLocalRange)
                    call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                    call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz, chunks_denom)
                    aaa=0.0
                    bbb=0.0
                    do ii = 1,ngroups
                        aaa = max(aaa,chunks_num(ii))
                        bbb = min(bbb,chunks_denom(ii))
                    end do
                    val_ptr(1) = (aaa + bbb)*0.5
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_BONDV1_CALC_UVW)
                    oclGlobalRange=(kp*jp)+(kp+2)*(ip+2)+(ip+3)*(jp+3)
                    oclLocalRange=0
                    call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz, val_ptr)
                    call runOcl(oclGlobalRange,oclLocalRange)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_VELFG__FEEDBF__LES_CALC_SM)
!#define NEW_VELFG
                    oclGlobalRange=ip*jp*kp
                    oclLocalRange=0
                    call runOcl(oclGlobalRange,oclLocalRange)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_LES_BOUND_SM)
                    oclGlobalRange = (jp+3)*(kp+2) + (kp+2)*(ip+2) + (jp+3)*(ip+2)
                    oclLocalRange = 0
                    call runOcl(oclGlobalRange,oclLocalRange)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_LES_CALC_VISC__ADAM)
                    oclGlobalRange=ip*jp*kp
                    oclLocalRange=0
                    call runOcl(oclGlobalRange,oclLocalRange)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_PRESS_RHSAV)
                    oclGlobalRange = 128*oclNunits
                    oclLocalRange = 128
                    ngroups = oclNunits
                    call runOcl(oclGlobalRange,oclLocalRange)
                    call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                    call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz, chunks_denom)
                    ! Calc the average over the compute units
                    rhsav = 0.0
                    area = 0.0
                    do ii = 1,ngroups ! number of work groups
                        !                        print *,'RHSAV CHUNK',ii, chunks_num(ii)
                        !                        print *,ii, chunks_denom(ii)
                        rhsav = rhsav + chunks_num(ii)
                        area = area + chunks_denom(ii)
                    end do
                    rhsav = rhsav / area
                    val_ptr(1) = rhsav
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_PRESS_SOR)
                    call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz, val_ptr)
                    iter = 0
                    sor = pjuge*1.1 ! just to be larger than pjuge
                    do while (sor > pjuge .and. iter < nmaxp)
                        iter = iter + 2
                        !  This is the "regular" approach
                        do nrd = 0,1
                                oclGlobalRange = (jp+2)*(kp+2)
                                oclLocalRange = jp+2
                                ngroups = kp+2
                                n_ptr(1)=nrd
                                call oclWrite1DIntArrayBuffer(n_ptr_buf,n_ptr_sz, n_ptr)
                                call runOcl(oclGlobalRange,oclLocalRange)
                                if (nrd == 1) then
                                    !read back chunks_num. Only once!
                                    call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                                    sor = 0.0
                                    do ii = 1,ngroups
                                        sor = sor + chunks_num(ii)
                                    end do
                                    sor = sqrt(sor) ! This is different from the original code but it is correct
                                end if
                            end do ! nrd
                            avg_iter = iter
                        end do ! while iter
! -----------------------------------------------------------------------------------------------------------------------------
                    case (ST_PRESS_PAV)
                        oclGlobalRange = 128*oclNunits ! so must loop over kp/oclNunits
                        oclLocalRange = 128
                        ngroups = oclNunits
                        call runOcl(oclGlobalRange,oclLocalRange)
                        call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                        call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz, chunks_denom)
                        ! Calc the average over the compute units
                        pav = 0.0
                        pco = 0.0
                        do ii = 1,ngroups
                            pav = pav + chunks_num(ii)
                            pco = pco + chunks_denom(ii)
                        end do
                        pav = pav / pco
                        val_ptr(1) = pav
! -----------------------------------------------------------------------------------------------------------------------------
                    case (ST_PRESS_ADJ)
                        oclGlobalRange=ip*jp*kp
                        oclLocalRange=0
                        call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz, val_ptr)
                        call runOcl(oclGlobalRange,oclLocalRange)
! -----------------------------------------------------------------------------------------------------------------------------
                    case (ST_PRESS_BOUNDP)
                        oclGlobalRange = (jp+2)*(kp+2) + (kp+2)*(ip+2) + (jp+2)*(ip+2)
                        oclLocalRange = 0
                        call runOcl(oclGlobalRange,oclLocalRange)
                        !
                        if ((mod(n,1000) == 0 .or. n == nmax)) then
                            ! read back results and write to file
                            call oclRead4DFloatArrayBuffer(uvw_buf,uvw_sz,uvw)
                            call convert_from_uvw(uvw,uo,vo,wo)
                            call oclRead4DFloatArrayBuffer(uvwsum_buf,uvwsum_sz,uvwsum)
                            call convert_from_fgh(uvwsum,usumo,vsumo,wsumo)
                            call oclRead4DFloatArrayBuffer(p_buf,p_sz,po)
                            nn = n/1000
                            print *, 'timestep: ',nn,' pressure at centre: ',po(0,ip/2,jp/2,kp/2),                                 'vel at centre: ',                                 uo(ip/2,jp/2,kp/2),vo(ip/2,jp/2,kp/2),wo(ip/2,jp/2,kp/2)
                        end if
                end select
            end do ! states loop
        end subroutine run_LES_kernel
     
        ! --------------------------------------------------------------------------------
        ! --------------------------------------------------------------------------------
        ! Auxiliary subroutines for file I/O
        ! --------------------------------------------------------------------------------
! TODO: could refactor into separate module
        subroutine write_uvw_p_to_file(u,v,w,p,im,jm,km,data20)
            use oclWrapper
            use params_common_sn
            real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1)  :: p
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: u
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: v
            real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: w
            integer :: im, jm, km
            character(len=70) :: data20
            open(unit=20,file=data20,form='unformatted',status='unknown')
            write(20) (((u(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((v(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((w(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((p(0,i,j,k),i=1,im),j=1,jm),k=1,km)
            close(unit=20)
        end subroutine write_uvw_p_to_file
        subroutine write_uvw_p_uvwsum_to_file(u,v,w,p,usum,vsum,wsum,im,jm,km,data20)
            use oclWrapper
            use params_common_sn
            real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1), intent(In)  :: p
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(In)  :: u
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
            real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(In)  :: w
            real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: usum
            real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: vsum
            real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: wsum
            integer :: im, jm, km
            character(len=70) :: data20
            open(unit=20,file=data20,form='unformatted',status='unknown')
            write(20) (((u(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((v(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((w(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((p(0,i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((usum(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((vsum(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((wsum(i,j,k),i=1,im),j=1,jm),k=1,km)
            close(unit=20)
        end subroutine write_uvw_p_uvwsum_to_file
        subroutine write_fgh_old_to_file(fold,gold,hold,im,jm,km,data21)
            use oclWrapper
            use params_common_sn
            real(kind=4), dimension(ip,jp,kp)  :: fold
            real(kind=4), dimension(ip,jp,kp)  :: gold
            real(kind=4), dimension(ip,jp,kp)  :: hold
            integer :: im, jm, km
            character(len=70) :: data21
            open(unit=21,file=data21,form='unformatted',status='unknown')
            write(21) (((fold(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(21) (((gold(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(21) (((hold(i,j,k),i=1,im),j=1,jm),k=1,km)
            close(unit=21)
        end subroutine write_fgh_old_to_file
    
    end module module_LES_combined_ocl
