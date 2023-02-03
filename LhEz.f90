    subroutine LhEz
    
    include 'com.txt'
    
    !include 'init1.txt'
    
    ! Ez on interface
    !EzRL = -Fxhat(:,:,:,6)
    !EzUD = Fyhat(:,:,:,5)
    !EzRL = -0.5*(FR(:,:,:,6) + FL(:,:,:,6))
    !EzUD = 0.5*(FU(:,:,:,5) + FD(:,:,:,5))
    
    ! 用更精确的特征值计算Ez
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                SRmax = UR(i,j,j1,2)/UR(i,j,j1,1)
                SLmax = UL(i + 1,j,j1,2)/UL(i + 1,j,j1,1)
                ! L-F Flux
                !SR = 0.3*min(abs(SRmax),abs(SLmax))
                SR = max(abs(SRmax),abs(SLmax))
                EzRL(i,j,j1) = -0.5*(FR(i,j,j1,7) + FL(i + 1,j,j1,7) - SR*(UL(i + 1,j,j1,7) - UR(i,j,j1,7)))
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                SRmax = UU(i,j,i1,3)/UU(i,j,i1,1)
                SLmax = UD(i,j + 1,i1,3)/UD(i,j + 1,i1,1)
                ! L-F Flux
                !SR = 0.3*min(abs(SRmax),abs(SLmax))
                SR = max(abs(SRmax),abs(SLmax))
                EzUD(i,j,i1) = 0.5*(FU(i,j,i1,6) + FD(i,j + 1,i1,6) - SR*(UD(i,j + 1,i1,6) - UU(i,j,i1,6)))
            end do
        end do
    end do
    
    ! calculate the value of U at vertex
    URU = 0
    ULU = 0
    URD = 0
    ULD = 0
    do i = 0,Nx1
        do j = 0,Ny1
            do d = 1,dimPk
                URU(i,j,:) = URU(i,j,:) + uh(i,j,d,:)*phiRU(d)
                ULU(i,j,:) = ULU(i,j,:) + uh(i,j,d,:)*phiLU(d)
                URD(i,j,:) = URD(i,j,:) + uh(i,j,d,:)*phiRD(d)
                ULD(i,j,:) = ULD(i,j,:) + uh(i,j,d,:)*phiLD(d)
            end do
        end do
    end do
    
    ! calculate Ez at vertex
    do i = 0,Nx
        do j = 0,Ny
            
            URU1 = ULD(i + 1,j + 1,:)
            ULU1 = URD(i,j + 1,:)
            URD1 = ULU(i + 1,j,:)
            ULD1 = URU(i,j,:)
            
            !if (i == 20 .and. j == 19) then
            !if (flux_type == 1) then
            call LF_Flux_2D
            !else if (flux_type == 2) then
            !    call HLL_Flux_2D
            !else if (flux_type == 3) then
                !call HLLC_Flux_2D
            !end if
            !end if
            EzVertex(i,j) = Ezhat
            !EzVertex(i,j) = 0.25*(EzRL(i,j,NumGLP) + EzRL(i,j + 1,1) + EzUD(i,j,NumGLP) + EzUD(i + 1,j,1))
            !if (i == 20 .and. j == 19) then
            !    print *,EzRL(i,j,NumGLP),EzRL(i,j + 1,1),EzUD(i,j,NumGLP),EzUD(i + 1,j,1)
            !    print *,Ezhat,0.25*(EzRL(i,j,NumGLP) + EzRL(i,j + 1,1) + EzUD(i,j,NumGLP) + EzUD(i + 1,j,1))
            !    print *," "
            !end if
        end do
    end do
    
    dEz1 = 0
    dEz2 = 0
    
    ! DG scheme of Ez1
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,k + 1
                do j1 = 1,NumGLP
                    dEz1(i,j,d) = dEz1(i,j,d) + 0.5*weight(j1)*EzRL(i,j,j1)*EzyG(j1,d)
                    !dEz1(i,j,d) = dEz1(i,j,d) + 0.5*weight(j1)*Ezf(xa + i*hx - t + dt,ya + (j - 0.5)*hy + hy1*lambda(j1) - t + dt)*EzyG(j1,d)
                end do
                dEz1(i,j,d) = (dEz1(i,j,d) - EzVertex(i,j)*EzU(d)/hy + EzVertex(i,j - 1)*EzD(d)/hy)/mmE(d)
                !dEz1(i,j,d) = (dEz1(i,j,d) - EzRL(i,j,NumGLP)*EzU(d)/hy + EzRL(i,j,1)*EzD(d)/hy)/mmE(d)
                !dEz1(i,j,d) = (dEz1(i,j,d) - Ezf(xa + i*hx - t + dt,ya + j*hy - t + dt)*EzU(d)/hy + Ezf(xa + i*hx - t + dt,ya + (j - 1)*hy - t + dt)*EzD(d)/hy)/mmE(d)
            end do
        end do
    end do
    
    ! DG scheme of Ez2
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,k + 1
                do i1 = 1,NumGLP
                    dEz2(i,j,d) = dEz2(i,j,d) - 0.5*weight(i1)*EzUD(i,j,i1)*EzxG(i1,d)
                    !dEz2(i,j,d) = dEz2(i,j,d) - 0.5*weight(i1)*Ezf(xa + (i - 0.5)*hx + hx1*lambda(i1) - t + dt,ya + j*hy - t + dt)*EzxG(i1,d)
                end do
                dEz2(i,j,d) = (dEz2(i,j,d) + EzVertex(i,j)*EzR(d)/hx - EzVertex(i - 1,j)*EzL(d)/hx)/mmE(d)
                !dEz2(i,j,d) = (dEz2(i,j,d) + EzUD(i,j,NumGLP)*EzR(d)/hx - EzUD(i,j,1)*EzL(d)/hx)/mmE(d)
                !dEz2(i,j,d) = (dEz2(i,j,d) + Ezf(xa + i*hx - t + dt,ya + j*hy - t + dt)*EzR(d)/hx - Ezf(xa + (i - 1)*hx - t + dt,ya + j*hy - t + dt)*EzL(d)/hx)/mmE(d)
            end do
        end do
    end do
    
    end subroutine LhEz