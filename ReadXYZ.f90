Program ReadXYZ

Use Elements, Only: AMS
    
Implicit Real(8) (A-H,O-Z)

!
! Explicit types definition
!
Integer(4), parameter::MaxAt=100,MaxSubStr=10

Character(255) Comment,FilInp,FilOut,Str,SubStr(MaxSubStr)
Character(10) Aname
Integer(4) NA(MaxAt)
Real(8) C(3,MaxAt),CM(3),T(3,3),PMI(3),Paxes(3,3),R(3),R1(3),U(3,3)

!
! File opening 
!
FilInp='c2h6.inp'
FilOut=FilInp
Call FileExtension(FilOut,'.out')
Open(5,File=FilInp)
Open(6,File=FilOut)

!
! Reading geometry
!
Call ReadNXYZ(5,MaxAt,Numat,NA,C,'*Geo')

Write(6,'(30x,''*** Program ReadXYZ ***'')')
Write(6,'('' Input file:'', a255)')FilInp
Write(6,'('' Numat     :'', i5)')Numat
Write(6,'(''*Geo'')')
!
! Calculation of TotMass and CM
!
TotMass=0.d0
CM=0.d0
Do i=1,Numat
    Call SetAName(NA(i),i,Aname)
    Write(6,'(a10,3f15.5)')AName,C(1:3,i)
    amsi=AMS(NA(i))
    CM=CM+amsi*C(1:3,i)
    TotMass=TotMass+amsi
Enddo
CM=CM/TotMass

Write(6,'(/'' Molecular mass  :'',f15.5)')TotMass
Write(6,'( '' Center of masses:'',3f15.10)')CM

!
! Reducing molecule to CM
!
Do i=1,Numat
    C(1:3,i)=C(1:3,i)-CM
Enddo

!
! Printing
!
Call PrintMol(6,Numat,NA,C,' Molecular coordinates reduced to mass center:')

T=0.d0
Do i=1,Numat
    amsi=AMS(NA(i))
    x=C(1,i)
    y=C(2,i)
    z=C(3,i)
    T(1,1)=T(1,1)+amsi*(y*y+z*z)
    T(2,2)=T(2,2)+amsi*(x*x+z*z)
    T(3,3)=T(3,3)+amsi*(x*x+y*y)
    T(2,1)=T(2,1)-amsi*y*x
    T(3,1)=T(3,1)-amsi*z*x
    T(3,2)=T(3,2)-amsi*z*y
Enddo
T(1,2)=T(2,1)
T(1,3)=T(3,1)
T(2,3)=T(3,2)

Write(6,'(/'' Tensor of inertia:'')')
Do i=1,3
    Write(6,'(3f12.6)')T(1:3,i)
Enddo

Call JacobiSorted(T,3,3,PMI,PAxes,nrot,1)

Write(6,'(/'' Principal moments of inertia (aem*A**2):'')')
Write(6,'(3f12.6)')PMI

Write(6,'(/'' Principal axes:'')')
Do i=1,3
    Write(6,'(3f12.6)')PAxes(i,1:3)
Enddo

Do i=1,Numat
    R=C(1:3,i)
    R1=matmul(Transpose(Paxes),R)
    C(1:3,i)=R1
Enddo

Call PrintMol(6,Numat,NA,C,' Molecule reduced to principal axes:')

Alp=5.d0
Do k=1,100
ca=DCOSD(Alp)
sa=DSIND(Alp)
U=0.d0
U(1,1)=ca
U(2,2)=ca
U(1,2)=-sa
U(2,1)=sa
U(3,3)=1.d0
Do i=1,Numat
    R=C(1:3,i)
    R1=matmul(U,R)
    C(1:3,i)=R1
Enddo
Call PrintMol(6,Numat,NA,C,' Rotated molecule:')
Enddo


Close(5)
Close(6)

End
!*****************************************************************
Subroutine ReadNXYZ(iu,MaxAt,Numat,NA,C,SearchStr)

Implicit Real(8) (A-H,O-Z)

Character(*) SearchStr
Character(Len(SearchStr)) SearchStr1
Character(255) Str,SubStr(10),Str1
Integer(4),parameter::MaxSubStr=10
Integer(4) NA(MaxAt)
Real(8) C(3,MaxAt)

SearchStr1=SearchStr
Call UCase(SearchStr1)

Numat=0
Do While(.not.EOF(iu))
    Read(iu,'(a255)')Str
    If (Len_Trim(Str)==0) Cycle
    Str1=Str
    Call UCase(Str1)
    If (INDEX(Str1,Trim(SearchStr1))==1) Then
        Do While(.not.EOF(iu))
            Read(5,'(a255)')Str
            If (Len_Trim(Str)==0) Exit
            Numat=Numat+1
            Call SubString(Str,10,nsubstr,SubStr)
            Call ParseAname(SubStr(1),NA(Numat))
            Read(SubStr(2),'(f255.10)')C(1,Numat)
            Read(SubStr(3),'(f255.10)')C(2,Numat)
            Read(SubStr(4),'(f255.10)')C(3,Numat)
        Enddo
    Endif
Enddo

    End
!******************************************************************
Subroutine PrintMol(iu,Numat,NA,C,Str)

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,Numat)
Integer(4) NA(Numat)
Character(10) Aname
Character(*) Str

ls=Len_Trim(Str)
    
Write(iu,'(/,a<ls>)')Str(1:ls)
Write(iu,'(''*Geo'')')
Do i=1,Numat
    Call SetAName(NA(i),i,Aname)
    Write(iu,'(a5,3f12.6)')Aname,C(1:3,i)
Enddo

End
