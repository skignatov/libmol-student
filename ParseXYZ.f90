Subroutine ParseXYZ(Mode,NNuc,X,Y,Z,line)
Use Elements, Only: ElName
Implicit Real(8) (A-H,O-Z)
Character(132) Line,lin(4),lin1
Character(255) ToUpperCase
Character(5) Anam,Anam1

! Subroutine ParseXYZ parses the input string Line extracting the
! atom number NNuc, and cartesian coords X, Y, Z.
! The line format is determined by Mode
! Mode==1 : AN, X, Y, Z
! Mode==2 : AN, dummy, X, Y, Z
! Mode==3 : dummy, AN, X, Y, Z
! where dummy is a variable to be skipped
! AN - element symbol, element symbol followed by number, or atomic number.

lin=repeat(' ',132)
lin1=repeat(' ',132)
k=0
Do i=1,5
	If (Mode==1.and.i==5) Cycle
	line=AdjustL(line)
	ll=Len_Trim(line)
	If (ll==0) Cycle
	iend=INDEX(line,' ')-1
	If (iend<=0.or.iend>ll) iend=ll
	lin1=line(1:iend)
	line(1:iend)=' '
	If (Mode==2.and.i==2) Cycle
	If (Mode==3.and.i==1) Cycle
	k=k+1
	lin(k)=lin1
	lin1=repeat(' ',132)
Enddo

lin1=lin(2)
Read(lin1,*)x
lin1=lin(3)
Read(lin1,*)y
lin1=lin(4)
Read(lin1,*)z

lin1=AdjustL(Lin(1))
ll=Len_Trim(lin1)
ic=ICHAR(lin1(1:1))
If (ic>=48.and.ic<=57) Then
	iend1=INDEX(lin1,' ')-1
	iend2=INDEX(lin1,'.')-1
	If (iend2<=0) iend2=999999
	iend=MIN0(iend1,iend2)
	lin(1)=lin1(1:iend)
	Read(lin(1),*)NNuc
Else
	Do i=1,ll
		ic=ICHAR(lin1(i:i))
		If (ic>=48.and.ic<=57) Exit
	Enddo
	iend=MIN0(i-1,ll)
	anam=lin1(1:iend)
	Anam=ToUpperCase(Anam)
	Do NNuc=1,99
		Anam1=ToUpperCase(Elname(Nnuc))
		If (Anam==Anam1) Exit
	Enddo
Endif

End