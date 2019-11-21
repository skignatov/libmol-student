Module Vars
    
    Implicit Real(8) (A-H,O-Z)

    Integer(4),parameter::MaxAt=10,MaxSubStr=10,MaxMol=500     
    Real(8) Box(3,2),CI(3,MaxAt),TotMassI,CMI(3),CJ(3,MaxAt),rmin
    Integer(4) NAI(MaxAt),NumatI,Nmol
    
End module
!**************************************************
Program RandomMoleculs
    
    Use Elements
    Use Vars
    
    Implicit Real(8) (A-H,O-Z)
    
    Character(255) Comment,Str,SubStr(MaxSubStr)
    Character(10) Aname
    Integer(4) NA(MaxAt)
    Real(8) C(3,MaxAt,MaxMol),U(3,3), X(3)

    Open(6,File='h2o.out')
! 
! Reading molecule
!
Call ReadMol(5)
 
!
! Box definition
!
Nmol=500
Box(1:3,1)=(/-5.d0,-5.d0,-5.d0/)
Box(1:3,2)=(/5.d0,5.d0,5.d0/)
rmin=1.5d0

Vol=1.d-30
Do k=1,3
    Vol=Vol*(Box(k,2)-Box(k,1))
Enddo
TotalMass=Dble(NMol)*TotMassI*1.67d-27
rho=TotalMass/Vol 

!
! Filling the box
!
A: Do k=1,Nmol
    icyc=0
    Do
        Call Place
        Call CheckDistance(k,C,iOK)
        If (iOK==1) Exit
        icyc=icyc+1
        If (icyc>1000000) Then
            Write(6,'('' Cannot place molecule'',i5,'' into the box! Stop.'')')k
            Nmol=k
            Exit A
        Endif
    Enddo
    C(1:3,1:NumatI,k)=CJ(1:3,1:NumatI)
Enddo A

Write(6,'(/'' Density, g/cm3 ='',g20.4)')rho*0.001d0
Write(6,'(/)')
Write(6,'(''*Geo'')')
Call PrintNXYZ(6,C,'*RotGeo')

end 
!**********************************************************
Subroutine CheckDistance(k,C,iOK)
Use Vars
Implicit Real(8) (A-H,O-Z)

Real(8) C(3,MaxAt,MaxMol)

iOK=1
Do imol=1,k-1
    Do ia=1,NumatI
        Do ja=1,NumatI
            rij=DSQRT((C(1,ia,imol)-CJ(1,ja))**2+   &
                      (C(2,ia,imol)-CJ(2,ja))**2+   &
                      (C(3,ia,imol)-CJ(3,ja))**2)
            If (rij<rmin) Then
                iOK=0
                Return
            Endif
        Enddo
    Enddo
Enddo

End
!***********************************************************
Subroutine ReadNXYZ(iu,Numat,NA,C,SearchStr)

Use Vars, Only: MaxAt,MaxSubStr

Implicit Real(8) (A-H,O-Z)

Character(*) SearchStr
Character(255) Str,SubStr(MaxSubStr),UpperCase
Integer(4) NA(MaxAt)
Real(8) C(3,MaxAt)

    Do While(.not.EOF(iu))
        Read(iu,'(a255)')Str
        If (Len_Trim(Str)==0) Cycle
        If (INDEX(UpperCase(Str),Trim(UpperCase(SearchStr)))==1) Then
            Numat=0
            Do While(.not.EOF(iu))
                Read(iu,'(a255)')Str
                If (INDEX(Str,'!')==1) Cycle
                If (Len_Trim(Str)==0) Exit
                Numat=Numat+1
                Call SubString(Str,MaxSubStr,nSubStr,SubStr)
                Call ParseAname(SubStr(1),NA(Numat))
                Read(SubStr(2),'(f255.10)')C(1,Numat)
                Read(SubStr(3),'(f255.10)')C(2,Numat)
                Read(SubStr(4),'(f255.10)')C(3,Numat)
            Enddo
            Exit
        Endif
    Enddo
   
    End
    !**************************************************
Subroutine PrintNXYZ(iu,C,Str)

Use Vars, Only:Nmol,MaxMol,Numat=>NumatI,NA=>NAI,MaxAt

Implicit Real(8) (A-H,O-Z)
Character(*) Str
Real(8) C(3,MaxAt,MaxMol)
Character(10) Aname

ll=Len_Trim(Str)
   
    Do j=1,NMol        
    Do i=1,Numat
        Call SetAName(NA(i),i,Aname)
        Write(iu,'(a10,3f15.6)')Aname,C(1:3,i,j)
    Enddo
    Enddo
    
End
!*****************************************************
Subroutine Place

Use Elements
Use Vars, NA=>NAI,Numat=>NumatI
       
Implicit Real(8) (A-H,O-Z)

Real(8) U(3,3),X(3)

! Random Euler Angles
    Call RANDOM_NUMBER(X)
    X=360.d0*X
    X(2)=X(2)*0.5d0

! Rotation of molecule
    irad=0
    Call RotMat2(X(1),X(2),X(3),irad,U)
    Do i=1,Numat
        CJ(1:3,i)=matmul(U,CI(1:3,i))
    Enddo

! Random coordinates
    Call RANDOM_NUMBER(X)
    X=Box(1:3,1)+(Box(1:3,2)-Box(1:3,1))*X
    Do i=1,NumAt
        CJ(1:3,i)=CJ(1:3,i)+X
    Enddo
End
!***********************************************
Subroutine ReadMol(iu)
Use Vars, NA=>NAI,Numat=>NumatI,CM=>CMI,TotMass=>TotMassI
Use Elements, Only:AMS

Implicit Real(8) (A-H,O-Z)

Character(10) Aname

! Reading
    Open(iu,File='h2o.inp')
    Call ReadNXYZ(iu,Numat,NA,CI,'*Geo')
    Close(iu)
! Coordinates about mass center 
    TotMass=0.d0
    CM=0.d0
    Do i=1,Numat
        Call SetAName(NA(i),i,Aname)
 !       Write(6,'(a10,3f15.5)')Aname,C(1:3,i)
        TotMass=TotMass+AMS(NA(i))
        CM=CM+AMS(NA(i))*CI(1:3,i)
    Enddo
    CM=CM/TotMass
    Write(6,'(/'' TotMass = '',f15.5)')TotMass
    Write(6,'(/'' Mass center: '')')
    Write(6,'(3f15.6)')CM
    
    Do i=1,Numat
        CI(1:3,i)=CI(1:3,i)-CM
    Enddo
    
End
