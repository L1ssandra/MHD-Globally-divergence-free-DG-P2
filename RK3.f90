    subroutine RK3
    
    include 'com.txt'
    
    t = 0
    count = 0
    
    call set_bc
    
    call TVB_Limiter
    
    call calculate_umax
        
    call calculate_totaldiv
    
    if (totaldiv > 1e-7) then
        call calculate_Az
        call div_free_Balsara
        call calculate_totaldiv
        call calculate_umax
    end if
        
    print *,t,"  ",umax,"  ",totaldiv
    
    do while (t < tend)
        
        call calculate_dt
    
        if (t + dt >= tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        ! Stage 1
        call set_bc
    
        call Lh
        
        call LhEz
        
        uI = uh + dt*du
        
        BxI = Bx + dt*dEz1
        ByI = By + dt*dEz2
        
        uh0 = uh
        Bx0 = Bx
        By0 = By
        
        uh = uI
        Bx = BxI
        By = ByI
        
        call TVB_Limiter
        
        call div_free_Balsara
        
        !call pp_Limiter
            
        call set_bc
        
        if (RKorder == 3) then
            
            !Stage 2
            call Lh
            
            call LhEz
            
            uh = uI + dt*du
            Bx = BxI + dt*dEz1
            By = ByI + dt*dEz2
        
            uII = (3d0/4d0)*uh0 + (1d0/4d0)*uh
            BxII = (3d0/4d0)*Bx0 + (1d0/4d0)*Bx
            ByII = (3d0/4d0)*By0 + (1d0/4d0)*By
        
            uh = uII
            Bx = BxII
            By = ByII
        
            call TVB_Limiter
            
            call div_free_Balsara
            
            !call pp_Limiter
            
            !Stage 3
            call set_bc
        
            call Lh
            
            call LhEz
            
            Bx = BxII + dt*dEz1
            By = ByII + dt*dEz2
            uh = uII + dt*du
        
            uh = (1d0/3d0)*uh0 + (2d0/3d0)*uh
            Bx = (1d0/3d0)*Bx0 + (2d0/3d0)*Bx
            By = (1d0/3d0)*By0 + (2d0/3d0)*By
        
            call TVB_Limiter
            
            call div_free_Balsara
            
            !call pp_Limiter
            
        end if
        
        count = count + 1
        
        call calculate_umax
        
        call calculate_totaldiv
        
        call calculate_pmin
        
        print *,t,"  ",umax,"  ",totaldiv
        
    end do
    
    end subroutine RK3