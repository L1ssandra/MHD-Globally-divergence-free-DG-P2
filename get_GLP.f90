    subroutine get_GLP
    
    include 'com.txt'
    
    if (NumGLP == 2) then
        lambda(1) = -0.5773502691896257645091488
        lambda(2) = 0.5773502691896257645091488
        
        weight(1) = 1
        weight(2) = 1
        
        lambdaL(1) = -1
        lambdaL(2) = 1
    else if (NumGLP == 3) then
        lambda(1) = -1
        lambda(2) = 0
        !lambda(3) = 1
        
        weight(1) = 1d0/3d0
        weight(2) = 4d0/3d0
        !weight(3) = 1d0/3d0
    else if (NumGLP == 4) then
        lambda(1) = -0.8611363115940525752239465
        lambda(2) = -0.3399810435848562648026658
        !lambda(3) = 0.3399810435848562648026658
        !lambda(4) = 0.8611363115940525752239465  
        
        weight(1) = 0.3478548451374538573730639     
        weight(2) = 0.6521451548625461426269361
        !weight(3) = 0.6521451548625461426269361     
        !weight(4) = 0.3478548451374538573730639
        
        ! Gauss-Lobatto Points
        lambdaL(1) = -1
        lambdaL(2) = -0.447213595499957939282
        !lambdaL(3) = 0.447213595499957939282
        !lambdaL(4) = 1
    else if (NumGLP == 5) then
        lambda(1) = -1
        lambda(2) = -0.6546536707079771437983
        !lambda(3) = 0
        !lambda(4) = 0.6546536707079771437983
        !lambda(5) = 1
        
        weight(1) = 0.1
        weight(2) = 49d0/90d0
        !weight(3) = 32d0/45d0
        !weight(4) = 49d0/90d0
        !weight(5) = 0.1
        
        
    else if (NumGLP == 6) then
        lambda(1) = -0.9324695142031520278123016     
        lambda(2) = -0.6612093864662645136613996    
        lambda(3) = -0.2386191860831969086305017     
        lambda(4) = 0.2386191860831969086305017     
        lambda(5) = 0.6612093864662645136613996     
        lambda(6) = 0.9324695142031520278123016     
        
        weight(1) = 0.1713244923791703450402961
        weight(2) = 0.3607615730481386075698335
        weight(3) = 0.4679139345726910473898703
        weight(4) = 0.4679139345726910473898703
        weight(5) = 0.3607615730481386075698335
        weight(6) = 0.1713244923791703450402961
        
        !lambda(1) = -1
        !lambda(2) = -0.8611363115940525752239465
        !lambda(3) = -0.3399810435848562648026658
        !lambda(4) = 0.3399810435848562648026658
        !lambda(5) = 0.8611363115940525752239465  
        !lambda(6) = 1
        
        !weight(1) = 0
        !weight(2) = 0.3478548451374538573730639     
        !weight(3) = 0.6521451548625461426269361
        !weight(4) = 0.6521451548625461426269361     
        !weight(5) = 0.3478548451374538573730639
        !weight(6) = 0
    end if
    
    end subroutine get_GLP