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
Real(8) C(3,MaxAt)

!
! Work starts here
!
FilInp='h2o.xyz'
FilOut=FilInp
Call FilExtension(FilOut,'.out')
Open(5,File=FilInp)
Open(6,File=FilOut)

Call ReadNXYZ(5,MaxAt,Numat,C,'*Geo')

Write(6,'(30x,''*** Program ReadXYZ ***'')')
Write(6,'('' Input file:'', a255)')FilInp
Write(6,'('' Numat     :'', i5)')Numat

TotMass=0.d0
Do i=1,Numat
    Call SetAName(NA(i),i,Aname)
    Write(6,'(a10,3f15.5)')AName,C(1:3,i)
    TotMass=TotMass+AMS(NA(i))
Enddo

Write(6,'(/'' Molecular mass = '',f15.5)')TotMass

Close(5)
Close(6)

End
!*****************************************************************
Subroutine ReadNXYZ(iu,MaxAt,Numat,NA,C,SearchStr)

Implicit Real(8) (A-H,O-Z)

Character(*) SearchStr
Integer(4) NA(MaxAt)
Real(8) C(3,MaxAt)


Numat=0
Do While(.not.EOF(iu))
    Read(iu,'(a255)')Str
    If (Len_Trim(Str)==0) Cycle
    If (INDEX(Str,Len_Trim(SearchStr))==1) Then
        Do While(.not.EOF(iu))
            Read(5,'(a255)')Str
            If (Len_Trim(Str)==0) Exit
            Numat=Numat+1
            Call SubString(Str,MaxSubStr,nsubstr,SubStr)
            Call ParseAname(SubStr(1),NA(Numat))
            Read(SubStr(2),'(f255.10)')C(1,Numat)
            Read(SubStr(3),'(f255.10)')C(3,Numat)
            Read(SubStr(4),'(f255.10)')C(3,Numat)
        Enddo
    Endif
Enddo

End

