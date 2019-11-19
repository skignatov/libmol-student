Program NLopt

Use NonLinOpt
    
Implicit Real(8) (A-H,O-Z)

Real(8) X0(100),G(100)

!
! Function to be minimized    
!
    iFunc=3 	! Type of function in CompFG
    N=2     	! Number of variables
    If (iFunc==2) n=4

!
! Starting point
!
    X0(1)=2.d0       
    X0(2)=15.d0
    X0(3)=10.d0
    X0(4)=-21.d0

!
! Optimization options
!
    MaxIter=1000
    TolG=1.d-10
    TolX=1.d-7
    iPrint=1
    iu=6

!
! File opening and initial printing
!
    Open(6,file='NonLinOpt.out')
    
    Write(6,'('' Type of function:'',i5)')iFunc
    Write(6,'('' Num of variables:'',i5)')N
    Write(6,'('' Tolerance on X  :'',e12.2)')TolX
    Write(6,'('' Tolerance on G  :'',e12.2)')TolG
    Write(6,'(/'' X0   ='',<N>f10.5)')(X0(j),j=1,N)
    Write(6,*)

!
! Optimization
!
    
    Call DFPMin(N,X0,iter,F,G)
 
!
! Print results
!   
    Write(6,'(//'' Optimization results:'')')
    Write(6,'('' Iter   = '',i10)')iter
    Write(6,'('' X      = '',<N>f10.5)')(X0(j),j=1,N)
    Write(6,'('' F      = '',f10.5)')F
    Write(6,'('' G      = '',<N>f10.5)')(G(j),j=1,N)
    Write(6,'('' F calls= '',i10)')nfcalls
    Write(6,'('' G calls= '',i10)')ngcalls
    
Close(6)

End
!*************************************************************************************************
Subroutine CompF(n,X,F)

Use NonLinOpt, Only: nfcalls,iFunc

Implicit Real(8) (A-H,O-Z)

Real(8) X(n)

nfcalls=nfcalls+1

If (iFunc==1) Then      ! Rosenbrock function: Xmin=(1,1) Fmin=0  Xhard=(-1.5,3)
    x1=X(1)
    x2=X(2)
    F=100.d0*(x2-x1*x1)**2+(1.d0-x1)**2
ElseIf (iFunc==2) Then  ! Powell function:  Xmin=(0,0,0,0)  Fmin=0  
    x1=X(1)
    x2=X(2)
    x3=X(3)
    x4=X(4)
    F=(x1+10.d0*x2)**2 + 5.d0*(x3-x4)**2 + (x2-2.d0*x3)**4 +10.d0*(x1-x4)**4
ElseIf (iFunc==3) Then  ! Test function: Xmin=(1,10) Fmin=0, sensitive to starting point
    x1=X(1)
    x2=X(2)
    da=0.1d0
    F=0.d0
    Do i=1,10
        a=da*Dble(i)
        F=F+((DEXP(-a*x1)-DEXP(-a*x2))-(DEXP(-a)-DEXP(-10.d0*a)))**2
    Enddo
Endif

End Subroutine CompF

