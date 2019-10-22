!*****************************************************
Subroutine ParseAname(Aname,NA)

Use Elements, Only: ElName

Implicit Real(8) (A-H,O-Z)

Character(*) Aname
Character(1) symb,symb1
Character(2) EN

ll=Len(Aname)
Do i=1,ll
    symb=Aname(i:i)
    ic=ICHAR(symb)      ! ��������� ������ � ASCII-���
    If (ic>=48.and.ic<=57) Aname(i:i)=' '   ! ���� ������-����� (ASCII-���� �� 48 �� 57, �� �������� �� �� ������)
    If (ic>=97.and.ic<=122) Aname(i:i)=CHAR(ic-32)             ! �������� �������� ����� �� ���������
Enddo

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
