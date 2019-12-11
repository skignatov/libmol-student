Subroutine PrintMol(iu,Numat,NA,C,Str)

Implicit Real(8) (A-H,O-Z)

Real(8) C(3,Numat)
Integer(4) NA(Numat)
Character(10) Aname
Character(*) Str

ls=Len_Trim(Str)
Write(iu,'(/,a<ls>)')Str(1:ls)
Do i=1,Numat
    Call SetAName(NA(i),i,Aname)
    Write(iu,'(a5,3f12.6)')Aname,C(1:3,i)
Enddo

End
