!*****************************************************
Subroutine ParseAname(Aname,NA,n)

Use Elements, Only: ElName

Implicit Real(8) (A-H,O-Z)

Character(*) Aname
Character symb*1,symb1*1,EN*2,buf*10

ll=Len(Aname)
n=0
kbuf=0
buf=repeat(' ',10)
Do i=1,ll
    symb=Aname(i:i)
    ic=ICHAR(symb)      ! переводим символ в ASCII-код
    If (ic>=48.and.ic<=57) Then
        Aname(i:i)=' '   ! Если символ-цифра (ASCII-коды от 48 до 57, то заменяем ее на пробел)
        kbuf=kbuf+1
        buf(kbuf:kbuf)=symb
    Endif
    If (ic>=97.and.ic<=122) Aname(i:i)=CHAR(ic-32)             ! Заменяем строчные буквы на заглавные
Enddo

If (kbuf>0) Then
    buf=AdjustR(buf)
    Read(buf,'(i10)')n
Endif

Do i=1,ll
    EN=ElName(i)
    Call UCase(EN)
    If (INDEX(Aname,EN)>0) Then
        NA=i
        Return
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
