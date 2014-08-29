module module_LES_conversions
    contains
    subroutine convert_to_uvw(u,v,w,uvw)
        use params_common_sn
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: u
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: v
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: w
        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1)  :: uvw        
        integer :: ii,jj,kk
        do jj = -1,jp+1
            do ii = 0,ip+1
                do kk = 0,kp+1
                    uvw(0,ii,jj,kk) = u(ii,jj,kk)
                    uvw(1,ii,jj,kk) = v(ii,jj,kk)
                    uvw(2,ii,jj,kk) = w(ii,jj,kk)
                    uvw(3,ii,jj,kk) = 0.0
                end do
               uvw(2,ii,jj,-1)=w(ii,jj,-1)
               uvw(1,ii,jj,-1)=0.0 !undefined!
               uvw(0,ii,jj,-1)=0.0 !undefined!
               uvw(3,ii,jj,-1)=0.0 !undefined!
           end do
       end do
    end subroutine convert_to_uvw

    subroutine convert_to_fgh(f,g,h,fgh)
        use params_common_sn
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: f,g,h
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp)  :: fgh
        integer :: ii,jj,kk
        do ii = 0,ip
            do jj = 0,jp
                do kk = 0,kp
                    fgh(0,ii,jj,kk) = f(ii,jj,kk)
                    fgh(1,ii,jj,kk) = g(ii,jj,kk)
                    fgh(2,ii,jj,kk) = h(ii,jj,kk)
                end do
            end do
        end do
    end subroutine convert_to_fgh

    subroutine convert_to_fgh_old(fold,gold,hold, fgh_old)
        use params_common_sn
        real(kind=4), dimension(ip,jp,kp)  :: fold,gold,hold
        real(kind=4), dimension(0:3,ip,jp,kp)  :: fgh_old
        integer :: ii,jj,kk
        do ii = 1,ip
            do jj = 1,jp
                do kk = 1,kp
                    fgh_old(0,ii,jj,kk) = fold(ii,jj,kk)
                    fgh_old(1,ii,jj,kk) = gold(ii,jj,kk)
                    fgh_old(2,ii,jj,kk) = hold(ii,jj,kk)
                end do
            end do
        end do
    end subroutine convert_to_fgh_old   
  
    subroutine convert_to_9vec(cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9,cov)
        use params_common_sn
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: cov1,cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov2, cov3, cov4, cov6, cov7, cov8, cov9
        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: cov ! We use 16 positions for alignment!
        integer :: ii,jj,kk
        do kk = 0,kp+2
            do jj = 0,jp+2
                do ii = 0,ip+2
                    cov(1,ii,jj,kk) = cov1(ii,jj,kk)
                    cov(2,ii,jj,kk) = cov2(ii,jj,kk)
                    cov(3,ii,jj,kk) = cov3(ii,jj,kk)
                    cov(4,ii,jj,kk) = cov4(ii,jj,kk)
                    cov(5,ii,jj,kk) = cov5(ii,jj,kk)
                    cov(6,ii,jj,kk) = cov6(ii,jj,kk)
                    cov(7,ii,jj,kk) = cov7(ii,jj,kk)
                    cov(8,ii,jj,kk) = cov8(ii,jj,kk)
                    cov(9,ii,jj,kk) = cov9(ii,jj,kk)
                end do
                cov(1,-1,jj,kk) = cov1(-1,jj,kk)
                cov(5,-1,jj,kk) = cov5(-1,jj,kk)
            end do
        end do
    end subroutine convert_to_9vec

 
    subroutine convert_from_9vec(cov,cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9)
        use params_common_sn
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: cov1,cov5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov2, cov3, cov4, cov6, cov7, cov8, cov9
        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: cov ! We use 16 positions for alignment!
        integer :: ii,jj,kk
        do kk = 0,kp+2
            do jj = 0,jp+2
                do ii = 0,ip+2
                     cov1(ii,jj,kk)= cov(1,ii,jj,kk)  
                     cov2(ii,jj,kk)= cov(2,ii,jj,kk)  
                     cov3(ii,jj,kk)= cov(3,ii,jj,kk)  
                     cov4(ii,jj,kk)= cov(4,ii,jj,kk)  
                     cov5(ii,jj,kk)= cov(5,ii,jj,kk)  
                     cov6(ii,jj,kk)= cov(6,ii,jj,kk)  
                     cov7(ii,jj,kk)= cov(7,ii,jj,kk)  
                     cov8(ii,jj,kk)= cov(8,ii,jj,kk)  
                     cov9(ii,jj,kk)= cov(9,ii,jj,kk)  
                end do
                 cov1(-1,jj,kk)= cov(1,-1,jj,kk)  
                 cov5(-1,jj,kk)= cov(5,-1,jj,kk)  
            end do
        end do
    end subroutine convert_from_9vec
    
    subroutine convert_to_bcdmask1(bmask1,cmask1,dmask1,mask1)
        use params_common_sn
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1)  :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: cmask1
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: dmask1
       
        real(kind=4), dimension(0:3,-1:ip+1,-1:jp+1,0:kp+1)  :: mask1
        integer :: ii,jj,kk
        do kk = 0,kp+1
            do jj = 0,jp+1
                do ii = 0,ip+1
                    mask1(0,ii,jj,kk)=bmask1(ii,jj,kk)
                    mask1(1,ii,jj,kk)=cmask1(ii,jj,kk)
                    mask1(2,ii,jj,kk)=dmask1(ii,jj,kk)
                end do
            end do
        end do
        do kk = 0,kp+1
            do jj = 0,jp+1
                mask1(0,-1,jj,kk)=bmask1(-1,jj,kk)
            end do
        end do
        do kk = 0,kp+1
            do ii = 0,ip+1
                mask1(1,ii,-1,kk)=cmask1(ii,-1,kk)
            end do
        end do
    end subroutine convert_to_bcdmask1

    subroutine convert_to_mask1(amask1,bmask1,cmask1,dmask1,mask1)
        use params_common_sn
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: amask1
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1)  :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: cmask1
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: dmask1
       
        real(kind=4), dimension(0:3,-1:ip+1,-1:jp+1,0:kp+1)  :: mask1
        integer :: ii,jj,kk
        do kk = 0,kp+1
            do jj = 0,jp+1
                do ii = 0,ip+1
                    mask1(0,ii,jj,kk)=bmask1(ii,jj,kk)
                    mask1(1,ii,jj,kk)=cmask1(ii,jj,kk)
                    mask1(2,ii,jj,kk)=dmask1(ii,jj,kk)
                    mask1(3,ii,jj,kk)=amask1(ii,jj,kk)
                end do
            end do
        end do
        do kk = 0,kp+1
            do jj = 0,jp+1
                mask1(0,-1,jj,kk)=bmask1(-1,jj,kk)
            end do
        end do
        do kk = 0,kp+1
            do ii = 0,ip+1
                mask1(1,ii,-1,kk)=cmask1(ii,-1,kk)
            end do
        end do
    end subroutine convert_to_mask1

    subroutine convert_to_cn234ls(cn2l,cn2s,cn3l,cn3s,cn4l,cn4s,cn234ls)
        use params_common_sn
        real(kind=4), dimension(ip)  :: cn2l
        real(kind=4), dimension(ip)  :: cn2s
        real(kind=4), dimension(jp)  :: cn3l
        real(kind=4), dimension(jp)  :: cn3s
        real(kind=4), dimension(kp)  :: cn4l
        real(kind=4), dimension(kp)  :: cn4s
        real(kind=4), dimension(1:2*(ip+jp+kp))  :: cn234ls
!        integer :: ii,jj,kk

        cn234ls(1:ip)= cn2l
        cn234ls(ip+1:2*ip) = cn2s
        cn234ls(2*ip+1:2*ip+jp)=cn3l
        cn234ls(2*ip+jp+1:2*ip+2*jp)=cn3s
        cn234ls(2*ip+2*jp+1:2*ip+2*jp+kp)=cn4l
        cn234ls(2*ip+2*jp+kp+1:2*ip+2*jp+2*kp)=cn4s

    end subroutine convert_to_cn234ls

    subroutine convert_from_uvw(uvw,u,v,w)
        use params_common_sn
        real(kind=4), dimension(    0:ip+1,-1:jp+1, 0:kp+1)  :: u
        real(kind=4), dimension(    0:ip+1,-1:jp+1, 0:kp+1)  :: v
        real(kind=4), dimension(    0:ip+1,-1:jp+1,-1:kp+1)  :: w
        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1)  :: uvw        
        integer :: ii,jj,kk
        do jj = -1,jp+1
            do ii = 0,ip+1
                do kk = 0,kp+1
                    u(ii,jj,kk) =    uvw(0,ii,jj,kk-1+1)
                    v(ii,jj,kk) =    uvw(1,ii,jj,kk-1+1)
                    w(ii,jj,kk) =    uvw(2,ii,jj,kk+1)
                end do
                w(ii,jj,-1)= uvw(2,ii,jj,-1+1)
            end do
        end do
    end subroutine convert_from_uvw

    subroutine convert_from_fgh(fgh,f,g,h)
        use params_common_sn
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: f,g,h
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp)  :: fgh
        !real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: f,g,h
        !real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,0:kp+1)  :: fgh
        integer :: ii,jj,kk
        do ii = 0,ip
            do jj = 0,jp
                do kk = 0,kp
                    f(ii,jj,kk) =    fgh(0,ii,jj,kk) 
                    g(ii,jj,kk) =    fgh(1,ii,jj,kk) 
                    h(ii,jj,kk) =    fgh(2,ii,jj,kk) 
                end do
            end do
        end do
    end subroutine convert_from_fgh

    subroutine convert_from_fgh_old(fgh_old,fold,gold,hold)
        use params_common_sn
        real(kind=4), dimension(ip,jp,kp)  :: fold,gold,hold
        real(kind=4), dimension(0:3,ip,jp,kp)  :: fgh_old
        integer :: ii,jj,kk
        do ii = 1,ip
            do jj = 1,jp
                do kk = 1,kp
!                print *,ii,jj,kk
                    fold(ii,jj,kk) = fgh_old(0,ii,jj,kk) 
                    gold(ii,jj,kk) = fgh_old(1,ii,jj,kk) 
                    hold(ii,jj,kk) = fgh_old(2,ii,jj,kk) 
                    !fghold(ii,jj,kk) = fgh_old(3,ii,jj,kk) 
               end do
            end do
        end do
    end subroutine convert_from_fgh_old

end module module_LES_conversions

