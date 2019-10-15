!*****************************************************
Subroutine ParseAname(Aname,NA)

Use Elements, Only: ElName

Implicit Real(8) (A-H,O-Z)

Character(*) Aname
Character(1) symb,symb1

ll=Len(Aname)
Do i=1,ll
    symb=Aname(i:i)
    ic=ICHAR(symb)      ! переводим символ в ASCII-код
    If (ic>=48.and.ic<=57) Aname(i:i)=' '   ! Если символ-цифра (ASCII-коды от 48 до 57, то заменяем ее на пробел)
    If (ic>=97.and.ic<=122) Aname(i:i)=CHAR(ic-32)             ! Заменяем строчные буквы на заглавные
Enddo

Do i=1,ll
    If (INDEX(Aname,ElName(i))>0) Then
        NA=i
        Return
    Endif
Enddo

End
