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
                Call ParseAname(SubStr(1),NA(Numat),nn)
                Read(SubStr(2),'(f255.10)')C(1,Numat)
                Read(SubStr(3),'(f255.10)')C(2,Numat)
                Read(SubStr(4),'(f255.10)')C(3,Numat)
            Enddo
            Exit
        Endif
    Enddo
   
End
