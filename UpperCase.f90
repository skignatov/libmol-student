Character(len=255) Function UpperCase(string)
	Character(len=*) string

    ll=Len_Trim(string)
	l1=Len(UpperCase)
	If (ll>l1) ll=l1
	UpperCase=repeat(' ',ll)
	Do i=1,ll
	    ic=ICHAR(string(i:i))
	    If (ic>=97 .and. ic<=122) Then
			UpperCase(i:i)=CHAR(ic-32)
		Else
		    UpperCase(i:i)=CHAR(ic)
		Endif
    Enddo

End
!************************************************************
Subroutine UCase(Str)

Character(*) Str
Character(1) sym

ls=Len_Trim(Str)
Do i=1,ls
    ic=ICHAR(Str(i:i))
    If (ic>=97.and.ic<=122) Str(i:i)=CHAR(ic-32)
Enddo

End

