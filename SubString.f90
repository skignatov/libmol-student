Subroutine SubString(String,MaxSubStr,NSubstr,SubStr)

Character(*) string
Character(*) substr(MaxSubStr)
		
ls=Len(string)
lss=Len(substr)
ll=Len_Trim(string)

do i=1,MaxSubStr
	substr(i)=repeat(' ',lss)
enddo

nsubstr=0

     i=0
1	 i=i+1
        If (i>ls) Return
	    If (string(i:i)==' ') goto 1
		k=0
		nsubstr=nsubstr+1
Cyc2:	Do j=i,ll
			If (string(j:j)==' ') exit Cyc2
			k=k+1
		SubStr(nsubstr)(k:k)=string(j:j)
		Enddo cyc2
	    i=i+k
		If (i>=ll) Return
	  goto 1

End
