Module Vars
    
Implicit Real(8) (A-H,O-Z)
Integer(4),parameter::MaxAt=100,MaxSubStr=20

Integer(4) Numat,NA(MaxAt),N
Real(8) C(3,MaxAt),Energy,Gnorm

Real(8) eps/5.d0/,sigma/2.d0/

Integer(4) nmax,ifunc

End module
!*********************************************************************************
Program GeoOpt
Use NonLinOpt, Only: DFPmin
Use Vars

Implicit Real(8) (A-H,O-Z)
Real(8) X(3*MaxAt),G(3*MaxAt),xx(3),Cmin(3,MaxAt)

! Program GeoOpt optimizes the structure of LJ clusters using various methods of non-linear optimization

Open(5,file='geoopt.inp')
Open(6,file='geoopt.out')

Call ReadNXYZ(5,Numat,NA,C,'*Geo')

Numat=13
Box=10.d0
Call RANDOM_SEED()

Do i=1,Numat
    NA(i)=3
    Call RANDOM_NUMBER(xx)
    C(1:3,i)=Box*(xx(1:3)-0.5d0)
Enddo


Call PrintMol(6,Numat,NA,C,' Initial coordinates:')

n=0
Do i=1,Numat
    Do k=1,3
        n=n+1
        X(n)=C(k,i)
    Enddo
Enddo

!
! Optimization methods
!
!Call SteepDesc(n,X,F,G)       	! Steepest Descent (simple with GoldenRatio 1D-min)
!Call DFPmin(n,X,iter,F,G) 	! DFP  (with Brent 1D-min)
Call lbfgsdriver(N,X,F,G)	! BFGS 


Call X2C(n,X)

Energy=F
Gnorm=DSQRT(dot_product(G,G))
Write(6,'(/'' Final Energy:'',f15.6)')Energy
Write(6,'( '' Final Gnorm :'',f15.6)')Gnorm
Call PrintMol(6,Numat,NA,C,' Optimized coordinates:')

Call Hessian(N,X)

End
!***********************************************
Subroutine X2C(n,X)
Use Vars, Only:Numat,C
Implicit Real(8) (A-H,O-Z)

Real(8) X(n)

n=0
Do i=1,Numat
    Do k=1,3
        n=n+1
        C(k,i)=X(n)
    Enddo
Enddo

End
!************************************************
Subroutine CompF(N,X,F)
Use Vars, Only: MaxAt,Numat,C,eps,sigma

Implicit Real(8) (A-H,O-Z)
Real(8) X(n)

Call X2C(n,X)

F=0.d0
Do i=2,Numat
    Do j=1,i-1
        rij=Distance(i,j,MaxAt,C)
        xx=(sigma/rij)**6
        F=F+4.d0*eps*(xx*xx-xx)
    Enddo
Enddo

    End
!************************************************
Subroutine PrintCoord(n,X)
Use Vars, Only: Numat,NA,C
Implicit Real(8) (A-H,O-Z)

Real(8) X(n)

Call X2C(n,X)

Call PrintMol(6,Numat,NA,C,' Current coordinates:')
        
    End
!***********************************************
Subroutine SteepDesc(n,X,F,G)
Use NonLinOpt, Only: CompFG
Implicit Real(8) (A-H,O-Z)

Real(8) X(n),X1(n),Xb(n),Xc(n),G(n),P(n)

MaxIter=1000000
Step=0.01d0

Do iter=1,MaxIter
    Call CompFG(0,n,X,F,G)
    gnorm=dot_product(G,G)
    If (gnorm<1.d-8) Exit
    Write(6,'('' Iter ='',i6,'' F ='',g20.5,5x,''Gnorm ='',g20.5)')iter,F,gnorm
    If (iter/100*100==iter) Call PrintCoord(n,X)
    P=-G
    a=0.d0
    b=step
    Call GoldSecMin(a,b,xmin,n,X,P)
    X=X+xmin*P
Enddo

Call CompFG(0,n,X,F,G)

    End
!*******************************************
Subroutine GoldSecMin(a,b,xmin,n,X0,P)
    
Implicit Real(8) (A-H,O-Z)
Real(8),parameter:: gr=0.5d0*(1.d0+DSQRT(5.d0)),eps=1.d-3
Integer(4) MaxIter/1000/
Common /VARS/ Ncall
Real(8) X0(n),P(n)

Ncall=0

!a=1.d0
!b=2.d0
x1=b-(b-a)/gr;f1=F(x1,n,X0,P)
x2=a+(b-a)/gr;f2=F(x2,n,X0,P)

Do iter=1,MaxIter
    If (DABS(b-a)<eps) Goto 10
    If (f1>=f2) Then
        a=x1; fa=f1
        x1=x2; f1=f2
        x2=a+(b-a)/gr;f2=F(x2,n,X0,P)
    Else
        b=x2; fb=f2
        x2=x1;f2=f1
        x1=b-(b-a)/gr;f1=F(x1,n,X0,P)
    Endif
Enddo
Write(*,'('' WARNING! Maximum number of iterations is achieved!'')')

10 xmin=0.5d0*(b+a)
!Write(*,'('' Xmin  :'',f10.5)')Xmin
!Write(*,'('' Fmin  :'',f10.5)')F(xmin)
!Write(*,'('' Iter  :'' i10  )')iter
!Write(*,'('' Fcalls:'',i10  )')Ncall
!pause   

End
!****************************************
Function F(x,n,X0,P)
Use NonLinOpt, Only:CompFG
Implicit Real(8) (A-H,O-Z)
Common/VARS/ Ncall
Real(8) X0(n),P(n),X1(n),G(n)

X1=X0+x*P
Call CompFG(1,n,X1,FF,G)
F=FF

End