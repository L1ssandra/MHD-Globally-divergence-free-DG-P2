    subroutine save_solution
    
    include 'com.txt'
    
    include 'init1.txt'
    
    open(unit = 1,file = 'Q1.txt')
    open(unit = 2,file = 'Q2.txt')
    open(unit = 3,file = 'Q3.txt')
    open(unit = 4,file = 'Q4.txt')
    open(unit = 5,file = 'Q5.txt')
    open(unit = 6,file = 'Q6.txt')
    open(unit = 7,file = 'Q7.txt')
    open(unit = 8,file = 'Q8.txt')
    
    uG = 0
    !call div_free_Balsara
    
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,NumEq
                do d = 1,dimPk
                    uG(i,j,:,:,n) = uG(i,j,:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do i1 = 1,NumGLP
            do j = 1,Ny
                do j1 = 1,NumGLP
                    write(1,*) uG(i,j,i1,j1,1)
                    write(2,*) uG(i,j,i1,j1,2)
                    write(3,*) uG(i,j,i1,j1,3)
                    write(4,*) uG(i,j,i1,j1,4)
                    write(5,*) uG(i,j,i1,j1,5)! - ureal(i,j,i1,j1,5)
                    write(6,*) uG(i,j,i1,j1,6)! - ureal(i,j,i1,j1,6)
                    write(7,*) uG(i,j,i1,j1,7)
                    write(8,*) uG(i,j,i1,j1,8)
                end do
            end do
        end do
    end do
    
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)
    
    open(unit = 1,file = 'Bx.txt')
    open(unit = 2,file = 'By.txt')
    open(unit = 3,file = 'Ez.txt')
    
    BxG = 0
    ByG = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,k + 1
                BxG(i,j,:) = BxG(i,j,:) + Bx(i,j,d)*EzG(:,d)
                ByG(i,j,:) = ByG(i,j,:) + By(i,j,d)*EzG(:,d)
            end do
        end do
    end do
    
    do i = 1,Nx
        do i1 = 1,NumGLP
            do j = 1,Ny
                do j1 = 1,NumGLP
                    
                    if (i1 == NumGLP) then
                        write(1,*) BxG(i,j,j1)! - Ez1real(i,j,j1)
                        !write(1,*) -FR(i,j,j1,6)
                        !write(1,*) EzRL(i,j,j1)!UR(i,j,j1,6)!
                        !write(2,*) Ez1real(i,j,j1)!Ezf(xa + i*hx - t + dt,Yc(j) + hy1*lambda(j1) - t + dt)
                    else
                        write(1,*) 0
                        !write(2,*) 0
                    end if
                    
                    if (j1 == NumGLP) then
                        write(2,*) ByG(i,j,i1)! - Ez2real(i,j,i1)
                        !write(2,*) EzUD(i,j,i1) - Ezf(Xc(i) + hx1*lambda(i1) - t + dt,ya + j*hy - t + dt)
                    else
                        write(2,*) 0
                    end if
                    
                    if (i1 == NumGLP .and. j1 /= NumGLP) then
                        write(3,*) EzRL(i,j,j1) !- Ezf(xa + i*hx - t + dt,Yc(j) + hy1*lambda(j1) - t + dt)
                    else if (i1 /= NumGLP .and. j1 == NumGLP) then
                        write(3,*) EzUD(i,j,i1) !- Ezf(Xc(i) + hx1*lambda(i1) - t + dt,ya + j*hy - t + dt)
                    else if (i1 == NumGLP .and. j1 == NumGLP) then
                        write(3,*) EzVertex(i,j) !- Ezf(Xc(i) + hx1*lambda(i1) - t + dt,Yc(j) + hy1*lambda(j1) - t + dt)
                    else
                        write(3,*) 0
                    end if
                    
                end do
            end do
        end do
    end do
    
    end subroutine save_solution