Program ReadXYZ

Use Elements, Only: AMS
    
Implicit Real(8) (A-H,O-Z)

!
! Explicit types definition
!
Integer(4), parameter::MaxAt=100,MaxSubStr=10

Character(255) Comment,FilInp,FilOut,Str,SubStr(MaxSubStr)
Character(10) Aname,Anm(MaxAt)
Integer(4) NA(MaxAt)
Real(8) C(3,MaxAt)

!
! Work starts here
!
Narg=Command_argument_count()
If (Narg>0) Then
    Call GET_COMMAND_ARGUMENT(1,FilInp)
Else
    Stop ' No arguments provided! Usage: console4 <file.inp>'
Endif

!FilInp='h2o.xyz'
FilOut=FilInp
Call FileExtension(FilOut,'.out')
Open(5,File=FilInp)
Open(6,File=FilOut)

Call ReadNXYZ(5,1,0,0,1,MaxAt,NA,C,NumAt,Anm,'*GEO')

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