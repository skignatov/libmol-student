Program GoldSecMin
    
Implicit Real(8) (A-H,O-Z)
Real(8),parameter:: gr=0.5d0*(1.d0+DSQRT(5.d0)),eps=1.d-12
Integer(4) MaxIter/1000/
Common /VARS/ Ncall

Ncall=0

a=1.d0
b=2.d0
x1=b-(b-a)/gr;f1=F(x1)
x2=a+(b-a)/gr;f2=F(x2)

Do iter=1,MaxIter
    If (DABS(b-a)<eps) Goto 10
    If (f1>=f2) Then
        a=x1; fa=f1
        x1=x2; f1=f2
        x2=a+(b-a)/gr;f2=F(x2)
    Else
        b=x2; fb=f2
        x2=x1;f2=f1
        x1=b-(b-a)/gr;f1=F(x1)
    Endif
Enddo
Write(*,'('' WARNING! Maximum number of iterations is achieved!'')')

10 xmin=0.5d0*(b+a)
Write(*,'('' Xmin  :'',f10.5)')Xmin
Write(*,'('' Fmin  :'',f10.5)')F(xmin)
Write(*,'('' Iter  :'' i10  )')iter
Write(*,'('' Fcalls:'',i10  )')Ncall
   
pause   

End
!****************************************
Function F(x)

Implicit Real(8) (A-H,O-Z)
Common/VARS/ Ncall

F=-Sin(x)
Ncall=Ncall+1

End