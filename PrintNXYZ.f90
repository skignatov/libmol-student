Subroutine PrintNXYZ(iU,iForm,N,NA,C,iLabel,Str)
Implicit Real(8) (A-H,O-Z)

Real(8) C(3,N)
Integer(4) NA(N)
Character(10) AName
Character(*) Str


!
! Subroutine PrintNXYZ print out the coordinates C(3,N) from array C(3,MaxAt) to the file iU 
! by format iForm = 0(DNXYZ), 1(NXYZ), 2(NDXYZ)
! Title is Str if iLabel>0
!


If (iLabel/=0) ls=Len_Trim(Str)
If (iLabel>0) Then
	Write(iU,'(<ls>a1)')(Str(i:i),i=1,ls)
ElseIf (iLabel<0) Then
	Write(iU,'(<ls>a1,i4)')(Str(i:i),i=1,ls),-iLabel
Endif

Do i=1,N
	nai=na(i)
	Call SetAName(NAi,i,AName)
	If		(iForm==1) Then					! NXYZ
		Write(iU,'(i4,2x,3f15.8)')NA(i),(C(k,i),k=1,3)
	ElseIf	(iForm==2) Then					! NDXYZ
		Write(iU,'(i4,2x,a5,2x,3f15.8)')NA(i),AName,(C(k,i),k=1,3)
	ElseIf	(iForm==3.or.iForm==0) Then		! DNXYZ (Default format)
		Write(iU,'(a5,i4,2x,3f22.14)')AName,NA(i),(C(k,i),k=1,3)
	ElseIf	(iForm==4) Then		! NXYZ with ANAME
		Write(iU,'(a5,2x,3f22.14)')AName,(C(k,i),k=1,3)
	Endif
Enddo
Write(iU,*)

End
