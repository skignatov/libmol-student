Subroutine ReadNXYZ(iUnit,Mode,iSkip,LeftDel,iForm,MaxAt,NA,C,NumAt,AName,SearchStr)
Implicit Real(8) (A-H,O-Z)

Character(*) SearchStr
Real(8) C(3,MaxAt)
Integer(4) NA(MaxAt)
Character(10) AName(MaxAt)
Character(255) String,ToUpperCase
Character(132) Str132

! Subroutine ReadNXYZ reads from file iUnit the atomic numbers NA and cartesian
! coordinates of atoms C(3,MaxAt). Number of atoms read in is NumAt.
! The read format is determined by iForm:
!	Mode==1		search first occurence of SearchStr in the file iUnit 
!               starting from the beginning of file
!	Mode==<+n>	search n-th  occurence of SearchStr in the file
!	Mode==<-n>	search n-th  occurence of SearchStr from the end of file 
!               (i.e. Mode==-1 corresponds to the last occurrence)
!	Mode==0		start reading coordinates from the (iSkip+1)th line of the file
!               without search of SearchStr in the file
!	if iUnit<0 -reading from the current position
!	iSkip		skip iSkip lines before reading the XYZ coordinates
!	iForm==1	the coordinates are in NXYZ format
!	iForm==2	the coordinates are in NDXYZ format
!	iForm==3	the coordinates are in DNXYZ
!	LeftDel		ignore first LeftDel symbols in the string with coordinates

!
! Find the starting position
!
If (iUnit>=0) Then
	Rewind(iUnit)
	iu=iUnit
Else
	iu=-iUnit
Endif

ifound=0
line=0
If (Mode>0) Then
	Do While (.not.EOF(iu))
		Read(iu,'(a255)') String
		line=line+1
		String=ToUpperCase(String)
		If (INDEX(String,Trim(SearchStr))>0) Then
			ifound=ifound+1
			linefound=line
			If (ifound==Mode) Goto 10
		Endif
	Enddo
ElseIf (Mode<0) Then
	Do While (.not.EOF(iu))
		Read(iu,'(a255)') String
		line=line+1
		String=ToUpperCase(String)
		If (INDEX(String,Trim(SearchStr))>0) Then
			ifound=ifound+1
			linefound=line
		Endif
	Enddo
	Rewind(iu)
	If (Mode==-1) Then
		Do i=1,linefound
			Read(iu,*)
		Enddo
		Goto 10
	Endif
	ifound1=ifound+Mode+1
	line=0
	ifound=0
	Do While (.not.EOF(iu))
		Read(iu,'(a255)') String
		line=line+1
		String=ToUpperCase(String)
		If (INDEX(String,Trim(SearchStr))>0) Then
			ifound=ifound+1
			linefound=line
			If (ifound==ifound1) Goto 10
		Endif
	Enddo
Endif

10 Continue
Do i=1,iSkip
	Read(iu,*)
Enddo

!
! Read in the coordinates
!

NumAt=0
Do While (.not.EOF(iu))
	Read(iu,'(a255)')String
	String=ToUpperCase(String)
	If (Len_Trim(String)==0) Exit
	If (INDEX(String,'------')>0) Exit
	If (INDEX(String,'......')>0) Exit
	If (INDEX(String,'$END')>0) Exit
	NumAt=NumAt+1
	Str132=String(LeftDel+1:132)
	Call ParseXYZ(iForm,NNuc,X,Y,Z,Str132)
	Call ReadWordLeft(Str132,AName(NumAt))
	NA(NumAt)=NNuc
	C(1,NumAt)=X
	C(2,NumAt)=Y
	C(3,NumAt)=Z
Enddo


End
!******************************************************************************
Subroutine ReadWordLeft(String,Str2)

Character(len=*) String,Str2
Character(len=Len(String)) Str1

! Subroutine ReadWordLeft reads a word (sequence of non-blank characters)
! from the beginning of Str1 and puts it to Str2

Str1=AdjustL(String)
l1=Len(Str1)
l2=Len(Str2)

If (l1==0.or.Len_Trim(Str1)==0) Then
	Str2(1:l2)=' '
	Return
Endif	

k=0
Do i=1,l1
	If (Str1(i:i)==' ') Exit
	k=k+1
	If (k>l2) Exit
	Str2(k:k)=Str1(i:i)
Enddo

If (k<l2) Str2(k+1:l2)=' '

Return
End


