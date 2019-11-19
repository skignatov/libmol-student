Module NonLinOpt

Implicit real(8) (A-H,O-Z)
    
Integer(4) nfcalls/0/,ngcalls/0/,iFunc,iu/6/,iPrint/0/,MaxIter/1000/
Real(8) tolx/1.2d-7/,tolg/1.d-7/
    
CONTAINS

SUBROUTINE DFPmin(n,p,iter,fret,g) 

Implicit Real(8) (A-H,O-Z)

Logical check
Real(8) dg(n),g(n),hdg(n),hess(n,n),pnew(n),xi(n)
Real(8) P(n)
Real(8), parameter ::STPMX=100.d0, EPS=3.d-10

!USES df unc, f unc, lnsrch 
!Given a starting point p(l:n) that is a vector of length n, the Broyden-Fletcher-Goldfarb- 
!Shanno variant of Davidon-Fletcher-Powell minimization is performed on a function func, 
!using its gradient as calculated by a routine df unc. The convergence requirement on zeroing 
!the gradient is input as gtol. Returned quantities are p(l:n) (the location of the mini- 
!mum), iter (the number of iterations that were performed), and fret (the minimum value 
!of the function). 
!
!The routine lnsrch is called to perform approximate line minimizations. 
!Parameters: NMAX is the maximum anticipated value of n; ITMAX is the maximum allowed 
!number of iterations; STPMX is the scaled maximum step length allowed in line searches; 
!TOLX is the convergence criterion on x values. 


!fp=func(p)              !Calculate starting function value and gradient,
call CompFG(0,n,p,fp,g)

sum=0.d0
hess=0.d0
do i=1,n             !and initialize the inverse Hessian to the unit matrix.
    hess(i,i)=1.
    xi(i)=-g(i)     !Initial line direction.
    sum=sum+p(i)**2
enddo 

stpmax=STPMX*max(dsqrt(sum),dble(n))

do its=1,MaxIter ! Main loop over the iterations.
    iter=its
    call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check)
    !The new function evaluation occurs in lnsrch; save the function value in fp for the next
    !line search. It is usually safe to ignore the value of check.
    
    fp=fret
    do i=1,n
        xi(i)=pnew(i)-p(i)  !Update the line direction,
        p(i)=pnew(i)        !and the current point.
    enddo
    
    test=0.d0               !Test for convergence on dx.
    do i=1,n
        temp=Dabs(xi(i))/max(Dabs(p(i)),1.d0)
        if(temp.gt.test)test=temp
    enddo

    dg(1:n)=g(1:n)
!    call dfunc(p,g)     !and get the new gradient.
    Call CompFG(2,n,p,fp,g)

    If (iprint>0) Write(iu,'('' Cyc:'',i5,2x,''X:''<n>f10.5,2x,'' F:'',g15.8,2x,''G:'',<n>g10.3)')its,p(1:n),fp,g(1:n)

    if(test.lt.TOLX) Then
        If (iPrint>0) Write(iu,'('' Test on X satisfied'')')
        return
    Endif

    test=0.d0           !Test for convergence on zero gradient.
    den=max(fret,1.d0)
    do i=1,n
        temp=Dabs(g(i))*max(Dabs(p(i)),1.d0)/den
        if(temp.gt.test)test=temp
    enddo

    if(test.lt.tolg) Then
        If (iPrint>0) Write(iu,'('' Test on G satisfied'')')
        return
    Endif

   
    dg(1:n)=g(1:n)-dg(1:n)          !Compute difference of gradients,
    do i=1,n                        !and difference times current matrix.
        hdg(i)=0.d0
        do j=1,n
            hdg(i)=hdg(i)+hess(i,j)*dg(j)
        enddo
    enddo

    fac=0.d0        !Calculate dot products for the denominators.
    fae=0.d0
    sumdg=0.d0
    sumxi=0.d0
    do i=1,n
        fac=fac+dg(i)*xi(i)
        fae=fae+dg(i)*hdg(i)
        sumdg=sumdg+dg(i)**2
        sumxi=sumxi+xi(i)**2
    enddo

    if(fac.gt.Dsqrt(EPS*sumdg*sumxi))then       !Skip update if fac not suciently positive.
        fac=1.d0/fac
        fad=1.d0/fae
        do i=1,n                                !The vector that makes BFGS dierent from DFP:
            dg(i)=fac*xi(i)-fad*hdg(i)
        enddo
        do i=1,n             !The BFGS updating formula:
            do j=i,n
                hess(i,j)=hess(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
                hess(j,i)=hess(i,j)
            enddo
        enddo
    endif

    do i=1,n                !Now calculate the next direction to go,
        xi(i)=0.d0
        do j=1,n
            xi(i)=xi(i)-hess(i,j)*g(j)
        enddo
    enddo
enddo                       !and go back for another iteration.

If (iPrint>0) Write(iu,'('' Max number of iterations is exceeded.'')')

END Subroutine DFPMin
!**********************************************************************************************************
SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check)

Implicit Real(8) (A-H,O-Z)

PARAMETER (ALF=1.d-4)
!,TOLX=1.e-7)
INTEGER(4) n
LOGICAL check
REAL(8) f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF
!EXTERNAL func

!C USES func 
!Given an n-dimensional point xold(1:n), the value of the function and gradient there,
!fold and g(1:n), and a direction p(1:n), nds a new point x(1:n) along the direction
!p from xold where the function func has decreased \suciently." The new function value
!is returned in f. stpmax is an input quantity that limits the length of the steps so that you
!do not try to evaluate the function in regions where it is undened or subject to overflow.
!p is usually the Newton direction. The output quantity check is false on a normal exit.
!It is true when x is too close to xold. In a minimization algorithm, this usually signals
!convergence and can be ignored. However, in a zero-nding algorithm the calling program
!should check whether the convergence is spurious.
!Parameters: ALF ensures sucient decrease in function value; TOLX is the convergence
!criterion on x.

INTEGER i
Real(8) a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam

check=.false.
sum=0.
do i=1,n
    sum=sum+p(i)*p(i)
enddo
sum=Dsqrt(sum)

if(sum.gt.stpmax)then   !Scale if attempted step is too big.
    do i=1,n
        p(i)=p(i)*stpmax/sum
    enddo
endif

slope=0.
do i=1,n
    slope=slope+g(i)*p(i)
enddo
if(slope.ge.0.) Stop 'roundoff problem in lnsrch'

test=0.     !Compute min.
do i=1,n
    temp=abs(p(i))/max(abs(xold(i)),1.d0)
    if(temp.gt.test)test=temp
enddo

alamin=TOLX/test
alam=1.d0       ! Always try full Newton step rst.
1 continue      ! Start of iteration loop.
    do i=1,n
        x(i)=xold(i)+alam*p(i)
    enddo
    !f=func(x)
    Call CompFG(1,n,X,f,G)
    
    if(alam.lt.alamin)then                      !Convergence on x. For zero nding, the calling program should verify the convergence.
        do i=1,n
            x(i)=xold(i)
        enddo
        check=.true.
        return
    else if(f.le.fold+ALF*alam*slope)then       !Sucient function decrease.
        return
    else                                        !Backtrack.
        if(alam.eq.1.)then                      !First time.
            tmplam=-slope/(2.*(f-fold-slope))
        else                                    !Subsequent backtracks.
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
                tmplam=-slope/(2.d0*b)
            else
                disc=b*b-3.d0*a*slope
                if(disc.lt.0.d0)then
                    tmplam=0.5d0*alam
                else if(b.le.0.d0)then
                    tmplam=(-b+Dsqrt(disc))/(3.d0*a)
                else
                    tmplam=-slope/(b+Dsqrt(disc))
                endif
            endif
            if(tmplam.gt.0.5d0*alam)tmplam=0.5*alam
        endif
    endif
    alam2=alam
    f2=f
    alam=max(tmplam,0.1d0*alam)
goto 1  !Try again.

END Subroutine lnsrch
!*******************************************************************
Subroutine CompFG(Mode,n,X,f,G)
Implicit Real(8) (A-H,O-Z)

Real(8) X(n),G(n)

    If (Mode==0.or.Mode==1) Then
        Call CompF(n,x,f)
    Endif
    If (Mode==0.or.Mode==2) Then
        ngcalls=ngcalls+1
!        Do i=1,n
!            x0=x(i)
!            dx=0.01d0*x0
!            x(i)=x0-dx
!            Call CompF(n,X,Fm)
!            x(i)=x0+dx
!            Call CompF(n,X,Fp)
!            x(i)=x0
!            G(i)=(Fp-Fm)/(2.d0*dx)
!        Enddo
        dx=0.01d0
        d2=dx+dx
        Do i=1,N
            x0=X(i)
            X(i)=x0-d2
            Call CompF(N,X,Fm2)
            X(i)=x0-dx
            Call CompF(N,X,Fm1)
            X(i)=x0+dx
            Call CompF(N,X,Fp1)
            X(i)=x0+d2
            Call CompF(N,X,Fp2)
            X(i)=x0
            G(i)=(-Fp2+8.d0*(Fp1-Fm1)+Fm2)/(12.d0*dx)
        Enddo
    Endif

End Subroutine CompFG
!*************************************************************************
End Module NonLinOpt