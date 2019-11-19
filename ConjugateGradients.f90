Module Vars
    
    Implicit Real(8) (A-H,O-Z)

    PARAMETER (NMAX=50)
    
    Real(8) X0(nmax),P(nmax)
    Integer(4) n,ifunc/0/,igrad/0/,iFuncType
    
End module
!********************************************
program ConjugateGradients

Use Vars, Only: n,nmax,iFuncType,ifunc,igrad

Implicit Real(8) (A-H,O-Z)

Real(8) X(nmax),G(nmax)

Open(6,File='ConGrad.out')

iFuncType=3 ! 1-Quadratic, 2-Rosenbrock, 3-Powell

n=2
If (iFuncType==3) n=4

X(1)= 3.d0
X(2)=-1.d0
X(3)= 0.d0
X(4)= 1.d0

Method=0
Call ConGrad(Method,X,n,iter,F,G,gnorm)            ! Conjugate gradients (Mode: 0 - Polak-Ribiere version, 1 - Fletcher-Reeves version, -1 - Steepest Descent)

Write(6,'(//'' Optimization finished with Niter ='',i6,''  Nfunc ='',i8,''  Ngrad ='',i8)')iter,ifunc,igrad
Write(6,'('' Final F ='',f12.6,''  Final Gnorm ='',f12.6)')F,Gnorm
Write(6,'('' X ='',<n>f10.5)')(X(j),j=1,n)
Write(6,'('' G ='',<n>f10.5)')(G(j),j=1,n)

Close(6)

end program 
!*********************************************************************
Subroutine CompFG(n,XX,F,G)

Use Vars, Only: iFuncType

Implicit Real(8) (A-H,O-Z)

Real(8) XX(n),G(n)


If (iFuncType==1) Then      ! Quadratic F=x**2+y**2

    x=XX(1)
    y=XX(2)
    
    F=x*x+y*y
    
    G(1)=2*x
    G(2)=2*y

Elseif (iFuncType==2) Then    ! Rosenbrock function (Xmin=(a, a**2), Fmin=0), hard starting point (-1.2,1)
    
    x=XX(1)
    y=XX(2)

    a=1.d0
    b=100.d0
    
    F=(a-x)**2 + b*(y-x*x)**2
    
    G(1)=-2.d0*(a-x)-4.d0*b*x*(y-x*x)
    G(2)=2.d0*b*(y-x*x)

ElseIf (iFuncType==3) Then      !Powell function (xmin=(0,0,0,0) F=0, Hard case X0=(3,-1,0,1))

    x1=XX(1)
    x2=XX(2)
    x3=XX(3)
    x4=XX(4)
    
    F= (x1+10.d0*x2)**2 + 5.d0*(x3-x4)**2 + (x2-2.d0*x3)**4 + 10.d0*(x1-x4)**4
    
    G(1)=2.d0*(x1 + 10.d0*x2) + 40.d0*(x1 - x4)**3
    G(2)=20.d0*(x1 + 10.d0*x2) + 4.d0*(x2 - 2*x3)**3
    G(3)=-8.d0*(x2 - 2.d0*x3)**3 + 10.d0*(x3 - x4)
    G(4)=-40.d0*(x1 - x4)**3 - 10.d0*(x3 - x4)
    
Endif

End
!**********************************************************************
SUBROUTINE ConGrad(Mode,X,n,iter,fret,G,gnorm) 

Use Vars, Only: nmax,ifunc

    Implicit Real(8) (A-H,O-Z)
    
    EXTERNAL func 
    
    PARAMETER (ITMAX=10000,EPS=1.d-10) 
    Real(8) X(n),P(n),H(n),G(n) 
    
    ! Conjugate gradients / steepest descent unconstrained optimization program
    ! Mode: 0 - Conjugate gradients (Polak-Ribiere)
    !       1 - Conjugate gradients (Fletcher-Reeves)
    !      -1 - Steepest descent
    
    ! Stopping criteria 
    ftol=1.d-6
    gtol=1.d-6
    
    ! F/G calculation at the initial point
    fp=func(X)                             
    call dfunc(X,G) 

    ! Initialize search direction and gradients
    P(1:n)=-G(1:n)
    H(1:n)=P(1:n)
    G(1:n)=H(1:n)
    
    do iter=1,ITMAX      ! Main loop 
        
        ! Linear search
        call linmin(X,G,n,fret,xmin)            
        
        ! F update and G calculation at the point located during linear search
        fchange=2.d0*Dabs(fret-fx)/(Dabs(fret)+Dabs(fx)+EPS)
        fp=fret 
        call dfunc(X,G) 
        
        Write(6,'('' Iter:'',i6,''  FCN ='',i8,''  F ='',f12.5,''  X ='',<n>f8.4,''  G ='',<n>f8.4)')iter,ifunc,fret,(X(j),j=1,n),(G(j),j=1,n)
        
        ! Check for stopping criterion on F
        if (fchange<=FTOL) Then
            Write(6,'(/'' Stopping criteria on Delta-F satisfied in CG program: Relative-Delta-F ='',4f20.12)')fchange
            return 
        Endif

        ! Prepare gradient corrections
        gg=0.d0
        dgg=0.d0
        do j=1,n 
            If (Mode==-1) Then
                gg=gg+G(j)**2                   ! Steepest descent
                G(j)=-G(j)
                Cycle
            Endif
            gg=gg+P(j)**2 
            If (Mode==0) Then
                dgg=dgg+(G(j)+P(j))*G(j)        ! CG Polak-Ribiere formula
            ElseIf (Mode==1) Then
                dgg=dgg+G(j)**2                 ! CG Fletcher-Reeves formula
            Endif
        enddo
        gnorm=Dsqrt(gg)
        
        ! Check for stopping criterion on G
        if (gg<gtol) Then
            Write(6,'(/'' Stopping criteria on Gnorm satisfied in CG program. Gnorm ='',4f20.12)')gnorm
            Return
        Endif
        
        If (Mode==-1) Cycle
        
        ! Update search direction and store gradients for CG
        gam=dgg/gg                          
        do j=1,n 
            P(j)=-G(j) 
            H(j)=P(j)+gam*H(j) 
            G(j)=H(j)
        enddo 

    enddo
    
    Write(6,'(''CG maximum iterations exceeded: '',i8)')ItMax

End
!**********************************************************
Function Func(X)
    
Use Vars, Only: n,ifunc
    
Implicit Real(8) (A-H,O-Z)
    
Real(8) X(n),G(n)
    
    Call CompFG(n,X,Func,G)
    
    ifunc=ifunc+1
    
    Write(6,'(i5,f12.5,<n>f10.4)')ifunc,Func,(X(j),j=1,n)
End
!**********************************************************************
Subroutine DFunc(X,G)
    Use Vars, only: n,igrad
    Implicit Real(8) (A-H,O-Z)
    Real(8) X(n),G(n)
    Call CompFG(n,X,F,G)
    igrad=igrad+1
End
!**********************************************************************
Function f1dim(x) 

Use Vars

Implicit Real(8) (A-H,O-Z)

Real(8) X1(n) 

X1(1:n)=X0(1:n)+x*P(1:n) 
f1dim=func(X1)

END 
!**********************************************************
Function df1dim(x) 

Use Vars

Implicit Real(8) (A-H,O-Z)

Real(8) dF(n),X1(n) 

X1(1:n)=X0(1:n)+x*P(1:n) 

call dfunc(X1,dF) 

df1dim=0.d0 
do i=1,n
    df1dim=df1dim+dF(i)*P(i)    !/pnorm**2
enddo

End
!**********************************************************
Subroutine linmin(X0,P,n,fret,xmin) 

Use Vars, Only: X0com=>X0,Pcom=>P,nmax

Implicit Real(8) (A-H,O-Z)

Real(8) X0(n),P(n)
PARAMETER (TOL=1.d-5)   !Maximum anticipated n, and TOL passed to brent. 
EXTERNAL f1dim,df1dim


ncom=n              !Set up the common block. 
do i=1,n 
    X0com(i)=X0(i) 
    Pcom(i)=P(i) 
enddo 

ax=0.d0               !Initial guess for brackets. 
xx=0.01d0 

call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)            ! Bracketing the interval with a minimum
fret=dbrent(ax,xx,bx,f1dim,df1dim,tol,xmin)     ! Find 1D-minimum at [ax,bx]

do i=1,n            !Construct the vector results to return. 
    P(i)=xmin*P(i)
    X0(i)=X0(i)+P(i)
enddo

End
!**********************************************************
Subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)

Implicit Real(8) (A-H,O-Z)

EXTERNAL func
PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)

fa=func(ax)
fb=func(bx)

if(fb.gt.fa)then                !Switch roles of a and b so that we can go downhill in the
    dum=ax                      !direction from a to b.
    ax=bx
    bx=dum
    dum=fb
    fb=fa
    fa=dum
endif

cx=bx+GOLD*(bx-ax)                  !First guess for c.
fc=func(cx)

1 if(fb.ge.fc)then                  !\do while": keep returning here until we bracket.
    r=(bx-ax)*(fb-fc)               !Compute u by parabolic extrapolation from a; b; c. TINY
    q=(bx-cx)*(fb-fa)               !is used to prevent any possible division by zero.
    u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*dsign(dmax1(dabs(q-r),TINY),q-r))
    ulim=bx+GLIMIT*(cx-bx)          !We won't go farther than this. Test various possibilities:
    if((bx-u)*(u-cx).gt.0.d0)then     !Parabolic u is between b and c: try it.
        fu=func(u)
        if(fu.lt.fc)then            !Got a minimum between b and c.
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
        else if(fu.gt.fb)then       !Got a minimum between between a and u.
            cx=u
            fc=fu
            return
        endif
        u=cx+GOLD*(cx-bx)           !Parabolic t was no use. Use default magnication.
        fu=func(u)
    else if((cx-u)*(u-ulim).gt.0.d0)then  !Parabolic t is between c and its allowed
        fu=func(u)                      !limit.
        if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
        endif
    else if((u-ulim)*(ulim-cx).ge.0.d0)then   !Limit parabolic u to maximum allowed value.
        u=ulim 
        fu=func(u)
    else                                    !Reject parabolic u, use default magnication.
        u=cx+GOLD*(cx-bx)
        fu=func(u)
    endif
    ax=bx                                   !Eliminate oldest point and continue.
    bx=cx
    cx=u
    fa=fb
    fb=fc
    fc=fu
    goto 1
endif

End
!********************************************************
FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)

Implicit Real(8) (A-H,O-Z)

EXTERNAL df,f
PARAMETER (ITMAX=100,ZEPS=1.0d-6)
LOGICAL ok1,ok2         !Will be used as flags for whether proposed steps are acceptable or not

a=dmin1(ax,cx) 
b=dmax1(ax,cx)
v=bx
w=v
x=v
e=0.d0
fx=f(x)
fv=fx
fw=fx
dx=df(x)        
dv=dx
dw=dx
do iter=1,ITMAX
    xm=0.5d0*(a+b)
    tol1=tol*Dabs(x)+ZEPS
    tol2=2.d0*tol1
    if(dabs(x-xm).le.(tol2-0.5d0*(b-a))) goto 3
    if(dabs(e).gt.tol1) then
        d1=2.d0*(b-a)                       !Initialize these d's to an out-of-bracket value.
        d2=d1
        if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)    !Secant method with one point.
        if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)    !And the other.
        u1=x+d1
        u2=x+d2
        ok1=((a-u1)*(u1-b).gt.0.d0).and.(dx*d1.le.0.d0)
        ok2=((a-u2)*(u2-b).gt.0.d0).and.(dx*d2.le.0.d0)
        olde=e                              !Movement on the step before last.
        e=d
        if(.not.(ok1.or.ok2))then           !Take only an acceptable d, and if both
            goto 1
        else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
                d=d1
            else
                d=d2
            endif
        else if (ok1) then
            d=d1
        else
            d=d2
        endif
        if(dabs(d).gt.dabs(0.5d0*olde))goto 1
        u=x+d
        if( (u-a).lt.tol2 .or. (b-u).lt.tol2 ) d=dsign(tol1,xm-x)
        goto 2
    endif
1   if(dx.ge.0.d0) then                       !Decide which segment by the sign of the derivative.
        e=a-x
    else
        e=b-x
    endif
    d=0.5d0*e                               !Bisect, not golden section.
2   if(dabs(d).ge.tol1) then
        u=x+d
        fu=f(u)
    else
        u=x+dsign(tol1,d)
        fu=f(u)
        if(fu.gt.fx)goto 3                  !If the minimum step in the downhill direction takes us uphill,
    Endif
    du=df(u)                            !Now all the housekeeping, sigh.
    if (fu.le.fx) then
        if(u.ge.x) then
            a=x
        else
            b=x
        endif
        v=w
        fv=fw
        dv=dw
        w=x
        fw=fx
        dw=dx
        x=u
        fx=fu
        dx=du
    else
        if(u.lt.x) then
            a=u
        else
            b=u
        endif
        if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
        else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
        endif
    endif
enddo
Write(6,'(''dbrent exceeded maximum iterations'')')

3 xmin=x
dbrent=fx

End

