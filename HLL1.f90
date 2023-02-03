    subroutine HLL1(UR,UL,FR,FL,SR,SL,Uhat,Fhat)
    
    real UR(6),UL(6),FR(6),FL(6),SR,SL,Uhat(6),Fhat(6)
    
    if (SR < 0) then
        Fhat = FR
        Uhat = UR
    else if (SL > 0) then
        Fhat = FL
        Uhat = UL
    else
        Fhat = ( SR*FL - SL*FR + SL*SR*(UR - UL) )/(SR - SL)
        Uhat = ( SR*UR - SL*UL - (FR - FL) )/(SR - SL)
    end if
    
    end subroutine HLL1