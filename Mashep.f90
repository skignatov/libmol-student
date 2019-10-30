Program MasEp

    USE, INTRINSIC :: IEEE_ARITHMETIC
    
    Implicit Real(8) (A-H,O-Z)
    
! Mashine epsylon
    eps=1.d0
    Do While (1.d0+eps/2.d0 > 1.d0)
        eps=eps/2.d0
    Enddo

! Mashine zero
    zero=1.d0
    Do While (zero/2.d0 > 0.d0)
        zero=zero/2.d0
    Enddo
    
    Write(*,*)eps,zero

! Max number of significant figures
    Write(*,*)INT(1.d0-DLOG10(eps*10.d0))

! Maximum number
    xmax=1.d0
    Do While(IEEE_IS_FINITE(xmax*2.d0)) 
        xmax=xmax*2.d0
    Enddo
    Write(*,*)xmax
    
    
    Pause
    
    End
    
    
    