Subroutine Hessian(N,X)
Use Elements, Only: AMS
Use Vars, Only: Numat,NA,C,Energy,Gnorm
Use NonLinOpt, Only: CompFG

Implicit Real(8) (A-H,O-Z)

Real(8) X(N),W(N),XM2(1,N),Gm2(N),Gm1(N),Gp2(N),Gp1(N),AMass(Numat)
Real(8) H(3*Numat,3*Numat),V(3*Numat,3*Numat)

Integer(4) iCode(1,2),iVec(n)
Character(10) Col0vec,Col0mat(n),Cvec(n),buf,Aname(Numat)

dx=0.01d0
d2=dx+dx
Do i=1,N
    x0=X(i)
    X(i)=x0-d2
    Call CompFG(0,N,X,Fm2,Gm2)
    X(i)=x0-dx
    Call CompFG(0,N,X,Fm1,Gm1)
    X(i)=x0+dx
    Call CompFG(0,N,X,Fp1,Gp1)
    X(i)=x0+d2
    Call CompFG(0,N,X,Fp2,Gp2)
    X(i)=x0
    H(i,1:N)=(-Gp2(1:N)+8.d0*(Gp1(1:N)-Gm1(1:N))+Gm2(1:N))/(12.d0*dx)
Enddo

call X2C(N,X)       ! Restore geometry

Do i=1,Numat
    Amass(i)=AMS(NA(i))
    Call SetAName(NA(i),i,Aname(i))
Enddo

! Diagonal matrix of 1/mass^(1/2)
Do i=1,Numat
    tmp=1.d0/DSQRT(AMS(NA(i)))
    k=3*(i-1)
    XM2(1,k+1:k+3)=tmp
Enddo

Call JacobiSorted(H,N,N,W,V,NRrot,1)

factor=108.516d0
Do i=1,N
    xm2i=XM2(1,i)     !1.d0/DSQRT(AMASS(i/3+1))
    wi=W(i)
    sgn=1.d0
    If (wi<0.d0) sgn=-sgn
    W(i)=sgn*DSQRT(DABS(Wi))*xm2i*factor
Enddo

    ik=0
    Do i=1,Numat
        Do k=1,3
            ik=ik+1
            buf='          '
            If (k==1) buf=Aname(i)
            If (k==1) buf(9:9)='X'
            If (k==2) buf(9:9)='Y'
            If (k==3) buf(9:9)='Z'
            Col0mat(ik)=buf
        Enddo
    Enddo
    
    Write(6,'(/'' Vibrational frequencies (cm-1) and mass-weighted cartesian modes of normal vibrations (3N vectors):'')')
    iCode(1,1)=1
    iCode(1,2)=1
    Ivec=0
    lenCol0=10
    Cvec   ='        '
    Col0vec='        '
!   call PrintVecMat(iUnit,nf,nd,NumCol,ldvec,ldmat,Mvec,N,lenCvec,iCode,Rvec,Ivec,Cvec,Amat,lenCol0,Col0vec,Col0mat)
   call PrintVecMat(    6,10, 3,    10,    1,    n,   1,N,     10,iCode,   W,Ivec,Cvec,   V,lenCol0,Col0vec,Col0mat)
!    call PrintVecMat(    6,10, 3,    10,    1,MaxOpt,  1,N,     10,iCode,   W,Ivec,Cvec,    ,lenCol0,Col0vec,Col0mat)
    Write(6,'(/'' Atomic masses (amu) and vibrationals frequencies (cm-1):'')')
    iend=0
    Do i=1,Numat
        ibeg=iend+1
        iend=ibeg+2
        Write(6,'(i5,2x,a10,f15.5,5x,3f15.5)')i,Aname(i),Amass(i),W(ibeg:iend)
    Enddo
    Write(6,*)

Open(7,file='geoopt.dat')
nFreq=3*Numat
Call PrintXYZFreq(7,Numat,nFreq,Energy,Gnorm,NA,C,W,V,3*Numat,'! Optimized LM structure')
Close(7)

    End
!********************************************************************
Subroutine PrintXYZFreq(iu,Numat,nFreq,Energy,Gnorm,NA,C,Freq,V,LDV,Str)
Implicit Real(8) (A-H,O-Z)

Integer(4) NA(Numat)
Real(8) C(3,Numat),Freq(nFreq),V(LDV,nFreq)
Character(*)Str

ln=Len(Str)
Write(iu,'(a<ln>)')Str
Write(iu,'(2i5,2f20.10)')Numat,nFreq,Energy,Gnorm
Do i=1,Numat
Write(iu,'(i5,3f20.10)')NA(i),C(1:3,i)
Enddo
If (nFreq>0) Write(iu,'(10f12.5)')Freq(1:nFreq)
Do i=1,nFreq
    Write(7,'(10f12.5)')V(1:3*Numat,i)
Enddo
    
End
!*********************************************************************************************************************
Subroutine PrintVecMat(iUnit,nf,nd,NumCol,ldvec,ldmat,Mvec,N,lenCvec,iCode,Rvec,Ivec,Cvec,Amat,lenCol0,Col0vec,Col0mat)

!Use BigArrays, Only: Amat=>V
Implicit Real(8) (A-H,O-Z)

Real(8) Rvec(ldvec,N),Amat(ldmat,N)
Integer(4) iCode(ldvec,2),Ivec(ldvec,N)
Character(lenCol0) Col0vec(1:Mvec),Col0mat(N),str0,str1,buf
Character(nf) Cvec(ldvec,N)

! Subroutine PrintVecMat1 prints square matrix Amat(N,N) in rectangular form and accompanying vectors X(N) (X may be real, integer, or character)
! to unit iUnit, by format f<nf>.<nd> and in NumCol columns 
!
! Input parameters:
!
! iUnit         - printing unit number
! nf, nd        - format f<nf>.<nd> to print the matrix values
! NumCol        - number of columns to print out the matrix and vectors
! LDvec         - leading dimensions for accompanying vectors Rvec,Ivec,Cvec, and code matrix iCode
! LDmat         - leading dimension for matrix Amat
! Mvec          - number of accompanying vectors to print out before Amat
! N             - dimension of matrix A and vectors Rvec,iVec,Cvec
! iCode(ldvec,2)- matrix of print codes to print the accompanying vectors: 
!                       iCode(i,1)=1 - print real vector Rvec() on i-th line
!                       iCode(i,1)=2 - print integer vector Ivec() on i-th line
!                       iCode(i,1)=3 - print character vector Cvec() on i-th line
!                       j=iCode(i,2) - gives  ordering number of vector j in arrays Rvec(j,1:N), Ivec(j,1:N), or Cvec(j,1:N) to be printed on i-th line
! Rvec(ldvec,n) - vector of real values to be printed
! Ivec(ldvec,n) - vector of integer values
! Cvec(ldvec,n) - vector of character values (len=lenCvec)
! Amat(lda,n)   - matrix to be printed out
! lenCol0       - length of character vectors to be printed on the left from vectors and matrix
! Col0vec(1:N)  - vector of character strings (len=LenCol0) to be printed on the left from vectors
! Col0mat(1:N)  - vector of character strings (len=LenCol0) to be printed on the left from matrix


    mc=lenCol0
    str0=repeat(' ',LenCol0)

    nx=nd-1
    ni=nf-nx
    
	jb=0
	20 ja=jb+1
		jb=jb+NumCol
		If(jb>N)jb=N
        
        Write(iUnit,'(/6x,<mc>a1,<NumCol>(i<ni>,<nx>x))')(str0(l:l),l=1,mc),(j,j=ja,jb)    ! Print upper ordering numbers
        
        Do k=1,Mvec
            
            ik1=iCode(k,1)
            ik2=iCode(k,2)
            
            str1=Col0Vec(k)

            If (ik1==1) Then        ! Print real vectors
		        Write(iUnit,'(6x,<mc>a1,<NumCol>(f<nf>.<nd>))')(str1(l:l),l=1,mc),(Rvec(ik2,j),j=ja,jb)
                
            Elseif (ik1==2) Then    ! Print integer vectors
		        Write(iUnit,'(6x,<mc>a1,<NumCol>(i<ni>,<nx>x))')(str1(l:l),l=1,mc),(IVec(ik2,j),j=ja,jb)
                
            Elseif (ik1==3) Then    ! Print character vectors
		        Write(iUnit,'(6x,<mc>a1,<NumCol>(a<nf>))')(str1(l:l),l=1,mc),(CVec(ik2,j),j=ja,jb)
            Endif
        Enddo
        
        Write(iUnit,*)
        Do i=1,N               ! Print matrix
            buf=Col0mat(i)
            Write(iUnit,'(i4,2x,<mc>a1,<NumCol>(f<nf>.<nd>))')i,(buf(l:l),l=1,mc),(Amat(i,j),j=ja,jb)
        Enddo

    If (jb<N) Goto 20
    Write(iUnit,*)

End
