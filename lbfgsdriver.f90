Subroutine LBFGSdriver(N,X,F,G)
Use NonLinOpt, Only: CompFG
Implicit Real(8) (A-H,O-Z)

      Integer(4), parameter:: MSAVE=7
      Real(8) X(n),G(n),DIAG(n),W(n*(2*MSAVE +1)+2*MSAVE),Vnorm,Discrep,Discrep1
      
      DOUBLE PRECISION F,EPS,XTOL,GTOL,T1,T2,STPMIN,STPMAX
      INTEGER IPRINT(2),IFLAG,ICALL,N,M,MP,LP,J
      LOGICAL DIAGCO
!C
!C                          JORGE NOCEDAL
!C                        *** July 1990 ***
!C     The driver for LBFGS must always declare LB2 as EXTERNAL
!C
      EXTERNAL LB2
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      
      stpmin=1.d-12
      stpmax=1.d12
      gtol=0.5d0
      
!C
!      N=100
      M=5
      IPRINT(1)= 1
      IPRINT(2)= 0
      lp=6
      mp=6
      
!C
!C     We do not wish to provide the diagonal matrices Hk0, and 
!C     therefore set DIAGCO to FALSE.
!C
      DIAGCO= .FALSE.
      EPS= 1.0D-6       !-5
      XTOL= 1.0D-12  !16
      ICALL=0
      IFLAG=0
      MaxCalls=1000

Do iCall=1,MaxCalls
    
    Call CompFG(0,N,X,F,G)
    gnorm=DSQRT(dot_product(G,G))
        
    CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
   
    IF(IFLAG.LE.0) Exit
    
        Write(6,'('' Cyc:'',i5,2x,'' F:'',g15.8,2x,'' Gnorm:'',g15.8)')icall,f,gnorm
        !!Write(6,'('' X:''<n>f10.5,2x,/'' G:'',<n>g10.3)')X(1:n),G(1:n)
        Call PrintCoord(n,X)
    
Enddo 

5 continue

    Call CompFG(0,N,X,F,G)

10 continue    

    End
