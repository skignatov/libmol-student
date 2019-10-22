Subroutine RotABG(A,B,G,iDeg,N,MaxAt,C0,C)
Implicit Real(8) (A-H,O-Z)

Real(8) C0(3,MaxAt),C(3,MaxAt),U(3,3)
Real(8) d/1.570796326794897d0/  ! Pi/2

!
! Subroutine RotABG rotates coordinates C0 by Euler angles A,B,G.
! (Radians if iDeg=0, Degrees if iDeg/=0)
!

If (iDeg==0) Then
    x=1.d0
Else
    x=0.01745329251994330d0 ! Conversion factor Deg->Rad: Pi/180 
Endif

! Uncomment +-d for another [Olkhovsky] definition of RotMat
f=a*x   !+d
t=b*x
p=g*x   !-d


cf=DCos(f)
ct=DCos(t)
cp=DCos(p)
sf=DSin(f)
st=DSin(t)
sp=DSin(p)


! Definition by [Korn]: U=A3(alp)*A2(bet)*A3(gam)
! Matrix Z1*Y2*Z3 of active transformation in Wikipedia
U(1,1) =  cf*ct*cp-sf*sp		
U(1,2) = -cf*ct*sp-sf*cp		
U(1,3) =  cf*st				
U(2,1) =  sf*ct*cp+cf*sp		
U(2,2) = -sf*ct*sp+cf*cp		
U(2,3) =  sf*st				
U(3,1) = -st*cp				
U(3,2) =  st*sp				
U(3,3) =  ct					
! For passive rotation, transpose the matrix (same rotation axis, same angles, but now the coordinate system rotates, rather than the vector)

! Transposed definition by [Olkhovsky],p.160, [Bukhgolz]: U=A3(alp)*A1(bet)*A3(gam) 
! Matrix Z1*X2*Z3 of active transformation in Wikipedia
!U(1,1) = cp*cf-ct*sf*sp		
!U(2,1) = cp*sf+ct*cf*sp		
!U(3,1) = sp*st		
!U(1,2) = -sp*cf-ct*sf*cp		
!U(2,2) =-sp*sf+ct*cf*cp 		
!U(3,2) = cp*st 		
!U(1,3) = st*sf		
!U(2,3) = -st*cf		
!U(3,3) = ct 			

Do i=1,N
	Do k=1,3
		C(k,i)=U(k,1)*C0(1,i)+U(k,2)*C0(2,i)+U(k,3)*C0(3,i)
	Enddo
Enddo

End
