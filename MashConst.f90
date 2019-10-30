Program MashConst
    
    Implicit Real(8) (A-H,O-Z)
    
    ! Calculation of Mashine Epsylon and Mashine Null
    
    
    Open(6,file='MashConst.out')
    
    ! Mashine Null
    
    A=1.d0
    i=0
    Do While (A>0.d0)
        B=A
        A=A/2.d0
        i=i+1
        !write(6,*)i,a,b
    Enddo
    Write(6,'(/''Number of divisons = '',i25   )')i
    Write(6,'( ''Mashine Null - the smallest value (obtained by division 1 by 2 ) distinguishable from 0  : '',1pE25.16)')B

    
    ! Mashine Eps
    
    A=1.d0
    i=0
    Do While (1.d0+A /= 1.d0)
        B=A
        A=A/2.d0
        i=i+1
    Enddo
    Write(6,'(/''Number of divisons = '',i25   )')i
    Write(6,'( ''Mashine Eps - the smallest number A (obtained by divison 1 by 2) which gives (1+A) /= 1  : '',1pE25.16)')B


    
    A=1.d0
    i=0
    Do While (1.d0-A /= 1.d0)
        B=A
        A=A/2.d0
        i=i+1
    Enddo
    Write(6,'(/''Number of divisons = '',i25   )')i
    Write(6,'( ''Mashine Eps - the smallest number A (obtained by divison 1 by 2) which gives (1-A) /= 1  : '',1pE25.16)')B


    A=1.d0
    i=0
    Do While (A /= A*2.d0)
        B=A
        A=A*2.d0
        i=i+1
        !write(6,*)i,a,b
    Enddo
    Write(6,'(/''Number of divisons = '',i25   )')i
    Write(6,'( ''Max number - the largest value A (obtained by multiplication 1 by 2 ) distinguishable from A*2  : '',1pE25.16)')B

    
    
    Call machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,maxexp,eps,epsneg,xmin,xmax)
    
    Write(6,'(//// '' *** Results of Machar suborutine, from Numerical Recipes in F77 ***'')')
    Write(6,'(// ''ibeta - the floating-point radix                                        = '',i25   )')ibeta
    Write(6,'(// ''it -  the number of base-ibeta digits in the floating-point mantissa    = '',i25   )')it
    Write(6,'(// ''iexp - the number of bits in the exponent (including its sign or bias). = '',i25   )')iexp
    Write(6,'(// ''negep - the exponent of the smallest power of ibeta that, subtracted''/ &
                 ''from 1.0, gives something different from 1.0.                           = '',i25   )')negep
    Write(6,'(// ''machep - exponent of the smallest (most negative) power of ibeta that, ''/ &
                 ''added to 1.0, gives something different from 1.0                        = '',i25   )')machep
    Write(6,'(// ''minexp - the smallest (most negative) power of ibeta consistent with'' / &
                 ''there being no leading zeros in the mantissa                            = '',i25   )')minexp
    Write(6,'(// ''maxexp - the smallest (positive) power of ibeta that causes overflow    = '',i25   )')maxexp
    Write(6,'(// ''eps - Smallest positive number that, added to 1.0, is not equal to 1.0  = '',1pE25.16)')eps
    Write(6,'(// ''epsneg - Smallest positive number that, subtracted from 1.0, is not '' / &
                 ''equal to 1.0                                                            = '',1pE25.16)')epsneg
    Write(6,'(// ''xmin - Smallest representable positive number  (ibeta**minexp)          = '',1pE25.16)')xmin
    Write(6,'(// ''xmax - Largest representable positive number (1-epsneg)*(ibeta**maxexp) = '',1pE25.16)')xmax
    
    Write(6,'(// ''irnd - A code in the range 0:5, giving information on what kind of''/ &
                 ''rounding is done in addition, and on how underflow is handled.          = '',i25   )')irnd
    Write(6,'(   ''   irnd = 0-3 - truncating imstead of rounding (not desirable)''/ &
                 ''   irnd = 4   - rounding (IEEE not compliant)''/ &
                 ''   irnd = 5   - IEEE compliant rounding ''/& 
                 ''   irnd = 3-5 also imply the IEEE compliant underflow processing'')')
    Write(6,'(// ''ngrd - the number of "guard digits" used when truncating the product of''/ &
                 ''two mantissas to fit the representation                                 = '',i25   )')ngrd
    
End
!****************************************************************************************************
SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,maxexp,eps,epsneg,xmin,xmax)
Implicit Real(8) (A-H,O-Z)

INTEGER ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
REAL(8) eps,epsneg,xmax,xmin
!Determines and returns machine-specic parameters aecting floating-point arithmetic. Returned
!values include ibeta, the floating-point radix; it, the number of base-ibeta digits
!in the floating-point mantissa; eps, the smallest positive number that, added to 1.0, is not
!equal to 1.0; epsneg, the smallest positive number that, subtracted from 1.0, is not equal to
!1.0; xmin, the smallest representable positive number; and xmax, the largest representable
!positive number. See text for description of other returned parameters.
INTEGER i,itemp,iz,j,k,mx,nxres
REAL(8) a,b,beta,betah,betain,one,t,temp,temp1,tempa,two,y,z ,zero,CONV

!CONV(i)=Dble(i)     !Change to dble(i), and change REAL declaration above to
!                    !DOUBLE PRECISION to one=CONV(1) nd double precision parameters.

one=Conv(1)
two=one+one
zero=one-one
a=one               !Determine ibeta and beta by the method of M. Malcolm.
1 continue
a=a+a
temp=a+one
temp1=temp-a
if ((temp1-one).eq.zero) goto 1
b=one
2 continue
b=b+b
temp=a+b
itemp=int(temp-a)
if (itemp.eq.0) goto 2
ibeta=itemp
beta=CONV(ibeta)
it=0                !Determine it and irnd.
b=one
3 continue
it=it+1
b=b*beta
temp=b+one
temp1=temp-b
if (temp1-one.eq.zero) goto 3
irnd=0
betah=beta/two
temp=a+betah
if (temp-a.ne.zero) irnd=1
tempa=a+beta
temp=tempa+betah
if ((irnd.eq.0).and.(temp-tempa.ne.zero)) irnd=2
negep=it+3          !Determine negep and epsneg.
betain=one/beta
a=one
do i=1, negep
a=a*betain
enddo !11
b=a
4 continue
temp=one-a
if (temp-one.ne.zero) goto 5
a=a*beta
negep=negep-1
goto 4
5 negep=-negep
epsneg=a
machep=-it-3        !Determine machep and eps.
a=b
6 continue
  
  
temp=one+a
if (temp-one.ne.zero) goto 7
a=a*beta
machep=machep+1
goto 6
7 eps=a
ngrd=0              !Determine ngrd.
temp=one+eps
if ((irnd.eq.0).and.(temp*one-one.ne.zero)) ngrd=1
i=0                 !Determine iexp.
k=1
z=betain
t=one+eps
nxres=0
8 continue          ! Loop until an underflow occurs, then exit.
y=z
z=y*y
a=z*one             ! Check here for the underflow.
temp=z*t
if ((a+a.eq.zero).or.(abs(z).ge.y)) goto 9
temp1=temp*betain
if (temp1*beta.eq.z) goto 9
i=i+1
k=k+k
goto 8
9 if (ibeta.ne.10) then
iexp=i+1
mx=k+k
else                !For decimal machines only.
iexp=2
iz=ibeta
10 if (k.ge.iz) then
iz=iz*ibeta
iexp=iexp+1
goto 10
endif
mx=iz+iz-1
endif
20 xmin=y           !To determine minexp and xmin, loop until an underflow oc
y=y*betain            !curs, then exit.
a=y*one             !Check here for the underflow.
temp=y*t
if (((a+a).ne.zero).and.(abs(y).lt.xmin)) then
k=k+1
temp1=temp*betain
if ((temp1*beta.ne.y).or.(temp.eq.y)) then
goto 20
else
nxres=3
xmin=y
endif
endif
minexp=-k           !Determine maxexp, xmax.
if ((mx.le.k+k-3).and.(ibeta.ne.10)) then
mx=mx+mx
iexp=iexp+1
endif
maxexp=mx+minexp
irnd=irnd+nxres         !Adjust irnd to reflect partial underflow.
if (irnd.ge.2) maxexp=maxexp-2  !Adjust for IEEE-style machines.

i=maxexp+minexp
!Adjust for machines with implicit leading bit in binary mantissa, and machines with radix
!point at extreme right of mantissa.
if ((ibeta.eq.2).and.(i.eq.0)) maxexp=maxexp-1


if (i.gt.20) maxexp=maxexp-1
if (a.ne.y) maxexp=maxexp-2
xmax=one-epsneg
if (xmax*one.ne.xmax) xmax=one-beta*epsneg
xmax=xmax/(beta*beta*beta*xmin)
i=maxexp+minexp+3
do  j=1,i
if (ibeta.eq.2) xmax=xmax+xmax
if (ibeta.ne.2) xmax=xmax*beta
enddo !12
return
END

Real(8) Function Conv(i)
CONV=Dble(i) 
End