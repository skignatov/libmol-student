	 Character(len=255) Function ToUpperCase(string)
	 Character(len=*) string

	 ll=Len_Trim(string)
!	 If (ll>132) ll=132
	 ToUpperCase=''
	 Do i=1,ll
	    ic=ICHAR(string(i:i))
	    If (ic>=97 .and. ic<=122) Then
			ToUpperCase(i:i)=CHAR(ic-32)
		Else
		    ToUpperCase(i:i)=CHAR(ic)
		Endif
     Enddo
	 Return
	 End
