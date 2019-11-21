Module EemData
        
    Real(8) Chi(120)/3.8196d0,4*0.d0,5.7254d0,6.8418d0,8.5000d0,112*0.d0/
    Real(8) Eta(120)/9.8832d0,4*0.d0,6.9235d0,6.3404d0,7.8386d0,112*0.d0/
    
End module
!***********************************************************
Program EEM
    Use LUdecomposition
    Use EEmData
    
    Implicit Real(8) (A-H,O-Z)
    
    Integer(4), parameter::MaxAt=10
    Real(8) C(3,MaxAt)
    Integer(4) NA(MaxAt)
    Real(8),allocatable:: A(:,:),X(:),B(:)
 
    Open(5,file='c2h6.inp')
    Open(6,file='c2h6.out')
    
    Call ReadNXYZ(5,MaxAt,Numat,NA,C,'*Geo')
    
    n1=Numat+1
    Allocate(A(n1,n1),B(n1),X(n1))
    
    Q=0.d0
    
!    Chi(1:8)=(/3.8196d0,0.d0,0.d0,0.d0,0.d0,5.7254d0,6.8418d0,8.5000d0/)
!    Eta(1:8)=(/9.8832d0,0.d0,0.d0,0.d0,0.d0,6.9235d0,6.3404d0,7.8386d0/)

    Do i=1,Numat
        A(i,i)=Eta(NA(i))
        B(i)=-Chi(NA(i))
        Do j=1,i-1
            rij=Distance(i,j,MaxAt,C)
            tmp=1.d0/rij
            A(i,j)=tmp
            A(j,i)=tmp
        Enddo
        A(i,n1)=-1.d0
        A(n1,i)=1.d0
    Enddo
    A(n1,n1)=0.d0
    B(n1)=Q
    
    Call LU_Solve(n1,A,n1,B,X)

    Do i=1,Numat
        Write(6,'('' Q'',i2,f10.4)')i,X(i)
    Enddo
    
End