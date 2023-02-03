    subroutine MHD_flux(Uh,Fx,Fy)
    
    real Uh(6),Fx(6),Fy(6)
    real rho,u,v,p,E,B1,B2,gamma,gamma1
    real S,T,K
      
    gamma = 5d0/3d0
    gamma1 = gamma - 1
    
    rho = Uh(1)
    u = Uh(2)/rho
    v = Uh(3)/rho
    E = Uh(4)
    B1 = Uh(5)
    B2 = Uh(6)
    
    p = gamma1*(E - 0.5*rho*(u**2 + v**2) - 0.5*(B1**2 + B2**2))
    S = p + 0.5*(B1**2 + B2**2)
    T = E + S
    K = u*B1 + v*B2
    
    Fx(1) = Uh(2)
    Fx(2) = rho*u**2 + SM - B1**2
    Fx(3) = rho*u*v - B1*B2
    Fx(4) = T*u - K*B1
    Fx(5) = 0
    Fx(6) = u*B2 - v*B1
    
    Fy(1) = Uh(3)
    Fy(2) = rho*u*v - B1*B2
    Fy(3) = rho*v**2 + S - B2**2
    Fy(4) = T*v - K*B2
    Fy(5) = v*B1 - u*B2
    Fy(6) = 0
    
    end subroutine MHD_flux