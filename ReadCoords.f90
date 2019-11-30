Module Vars
    
    Implicit Real(8) (A-H,O-Z)
    
    Integer(4),parameter::MaxAt=20,MaxSubStr=10
    Integer(4)  NA(MaxAt),      &   ! Atomic numbers of real atoms
                NAX(MaxAt),     &   ! Atomic numbers of all atoms including dummy atoms (XX)
                NABC(3,MaxAt),  &   ! Reference atom numbers in Z-matrix
                Numat,          &   ! Number of real atoms
                Natoms              ! Number of all atoms
        
    Real(8) C(3,MaxAt),         &   ! Coords of real atoms
            CX(3,MaxAt),        &   ! Coords of all atoms including dummies (XX)
            A(3,MaxAt)              ! Z-matrix parameters
    
    Character(255) Str,SubStr(MaxSubStr),UpperCase
    
    External UpperCase
    
End module
!********************************************************************************
Program ReadCoords
Use Vars

Open(5,file='readcoords.dat')
Open(6,file='readcoords.out')

Do While(.not.EOF(5))
    Read(5,'(a255)')Str
    If (Len_Trim(Str)==0) Cycle
    Str=UpperCase(Str)
    If (Index(Str,'*GEO')==1) Then
        Call ReadGZsimple(5)
        Exit
    Endif
Enddo

Write(6,*)
Call PrintMol(6,Natoms,NAX,CX,' Coordinates with dummy atoms')
Write(6,*)
Call PrintMol(6,Numat,NA,C,' Coordinates without dummy atoms')

End
!*******************************************************************************
Subroutine ReadGZsimple(iu)
Use Vars
    
Implicit Real(8)(A-H,O-Z)

! Reading simple Z-matrix (without variables)
Natoms=0    ! Number of atoms including dummies (XX atoms). Numat - number of real atoms
Do While(.not.EOF(5))
    Read(5,'(a255)')Str
    If (Len_Trim(Str)==0) Exit
    Str=AdjustL(Str)
    If (Str(1:1)=='!') Cycle
    Natoms=Natoms+1
    Call SubString(Str,MaxSubstr,nsubstr,SubStr)
    Call ParseAname(SubStr(1),NAX(Natoms),itmp)
    If (Natoms==1) Cycle
    Call ParseAname(SubStr(2),itmp,NABC(1,Natoms))
    Read(SubStr(3),*)A(1,Natoms)
    If (Natoms==2) Cycle
    Call ParseAname(SubStr(4),itmp,NABC(2,Natoms))
    Read(SubStr(5),*)A(2,Natoms)
    If (Natoms==3) Cycle
    Call ParseAname(SubStr(6),itmp,NABC(3,Natoms))
    Read(SubStr(7),*)A(3,Natoms)
Enddo

Write(6,'(/'' Z-matrix'')')
Do i=1,Natoms
    Write(6,'(i5,3f12.5,3i5)')NAX(i),A(1:3,i),NABC(1:3,i)
Enddo
     
Call ZM2XYZ     !It's reasonble to make it as a subroutine because ZM->XYZ transformation is usually doing many times during the optimization 
    
    
End
!********************************************************************
Subroutine ZM2XYZ
Use Vars

Implicit Real(8) (A-H,O-Z)

Real(8) X(3),Y(3),Z(3)


CX=0.d0
If (Natoms==1) Return
CX(1,2)=A(1,2)
If (Natoms==2) Return
CX(1,3)=A(1,3)*DCOSD(A(2,3))
CX(2,3)=A(1,3)*DSIND(A(2,3))
If (NABC(1,3)==2) CX(1,3)=CX(1,2)-CX(1,3)
If (Natoms==3) Return

Do i=4,Natoms
    R=A(1,i)
    Alp=A(2,i)
    Theta=A(3,i)
	n1=NABC(1,i)
    n2=NABC(2,i)
    n3=NABC(3,i)

    IF (DABS(Alp-180.D0)<1.D-6) THEN
        X=C(1:3,n1)-C(1:3,n2)
		C(1:3,i)=C(1:3,n1)+A(1,i)*X/DSQRT(dot_product(X,X))
        Cycle
	EndIf

	X=CX(1:3,n3)-CX(1:3,n2)
	Y=CX(1:3,n1)-CX(1:3,n2)
    Call CrossProduct(X,Y,Z)
    Call CrossProduct(Y,Z,X)
    X=X/DSQRT(dot_product(X,X))
    Y=Y/DSQRT(dot_product(Y,Y))
    Z=Z/DSQRT(dot_product(Z,Z))
	d=-R*DCOSD(Alp)
	e= R*DSIND(Alp)*DCOSD(Theta)
	h= R*DSIND(Alp)*DSIND(Theta)
	CX(1:3,i)=CX(1:3,n3)+e*X+d*Y+h*Z
Enddo

Numat=0
Do i=1,Natoms
    If (NAX(i)==98.or.NAX(i)==99) Cycle
    Numat=Numat+1
    NA(Numat)=NAX(i)
    C(1:3,Numat)=CX(1:3,i)
Enddo
    
    END
!***********************************************************************
Subroutine PrintMol(iu,Numat,NA,C,Str)

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,Numat)
Integer(4) NA(Numat)
Character(10) Aname
Character(*) Str

ls=Len_Trim(Str)
Write(iu,'(/,a<ls>)')Str(1:ls)
Do i=1,Numat
    Call SetAName(NA(i),i,Aname)
    Write(iu,'(a5,3f12.6)')Aname,C(1:3,i)
Enddo

End
