    subroutine pp_Limiter
    
    include 'com.txt'
    
    real eta,epsilon,eta1,pbar,pq
    
    epsilon = 1e-13
    
    ! Limiting the density
    do i = 1,Nx
        do j = 1,Ny
            eta = 1
            uhGLL = 0
            do d = 1,dimPk
                do n = 1,NumEq
                    uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,d,n)*phiGLL(:,:,d,:)
                end do
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do d = 1,2
                        ! eta1 = (rhobar - epsilon)/(rhobar - rho(xq))
                        eta1 = abs((uh(i,j,1,1) - epsilon)/(uh(i,j,1,1) - uhGLL(i1,j1,1,d)))
                        if (eta1 < 1) then
                            eta = eta1
                        end if
                    end do
                end do
            end do
            uh(i,j,2:dimPk,1) = eta*uh(i,j,2:dimPk,1)
        end do
    end do
    
    ! Limiting the pressure
    do i = 1,Nx
        do j = 1,Ny
            eta = 1
            eta1 = 1
            uhGLL = 0
            do d = 1,dimPk
                do n = 1,NumEq
                    uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,d,n)*phiGLL(:,:,d,:)
                end do
            end do
            pbar = pressure(uh(i,j,1,1),uh(i,j,1,2),uh(i,j,1,3),uh(i,j,1,4),uh(i,j,1,5),uh(i,j,1,6),uh(i,j,1,7),uh(i,j,1,8),gamma)
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do d = 1,2
                        pq = pressure(uhGLL(i1,j1,1,d),uhGLL(i1,j1,2,d),uhGLL(i1,j1,3,d),uhGLL(i1,j1,4,d),uhGLL(i1,j1,5,d),uhGLL(i1,j1,6,d),uhGLL(i1,j1,7,d),uhGLL(i1,j1,8,d),gamma)
                        if (pq < 0) then
                            eta1 = abs(pbar/(pbar - pq))
                        end if
                        if (eta1 < eta) then
                            eta = eta1
                        end if
                    end do
                end do
            end do
            
            if (eta < 1) then
                eta = 0.9*eta
            end if
            
            uh(i,j,2:dimPk,:) = eta*uh(i,j,2:dimPk,:)
            
        end do
    end do
    
    end subroutine pp_Limiter