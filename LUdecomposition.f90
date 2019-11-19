Module LUdecomposition

INTERFACE outerprod
Module Procedure outerprod
END INTERFACE

INTERFACE swap
MODULE PROCEDURE swap 
END INTERFACE

INTERFACE imaxloc
MODULE PROCEDURE imaxloc
END INTERFACE

CONTAINS 

FUNCTION outerprod(a,b)
REAL(8), DIMENSION(:), INTENT(IN) :: a,b
REAL(8), DIMENSION(size(a),size(b)) :: outerprod
outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod

SUBROUTINE swap(a,b)
REAL(8), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(8), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap

FUNCTION imaxloc(arr)
!Index of maxloc on an array.
REAL(8), DIMENSION(:), INTENT(IN) :: arr
INTEGER(4) :: imaxloc_d
INTEGER(4), DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc=imax(1)
END FUNCTION imaxloc

Subroutine LU_Solve(N,A,LDA,B,X)
Implicit Real(8) (A-H,O-Z)

Real(8) A(LDA,N),B(N),X(N),AA(LDA,N),BB(N)
Integer(4) Indx(N)

AA=A
Call LU_Decomposition(N,AA,LDA,Indx,d)

BB=B
Call LU_BackSubstitution(N,AA,LDA,Indx,BB)
X=BB

End Subroutine LU_Solve
!***********************************************************************************************    
Subroutine LU_Decomposition(N,A,LDA,Indx,d)
!Use LUdecomposition
Implicit Real(8) (A-H,O-Z)

Real(8) A(LDA,N)
Integer(4) Indx(N)

!Given an N x N input matrix a, this routine replaces it by the LU decomposition of a
!rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!output vector of length N that records the row permutation effected by the partial pivoting;
!d is output as +1 depending on whether the number of row interchanges was even or odd,
!respectively. This routine is used in combination with lubksb to solve linear equations or
!invert a matrix.

Real(8) vv(N)    ! vv stores the implicit scaling of each row
Real(8), parameter :: TINY=1.0d-20  ! Small number


d=1.0                                                   ! No row interchanges yet.
vv=maxval(DABS(a),dim=2)                                ! Loop over rows to get the implicit scaling information
if (any(vv == 0.0)) Then
    Write(*,'(//'' ERROR! Singular matrix A in LU_Decomposition'')')   
    Stop
Endif

vv=1.d0/vv                                              ! Save the scaling.

do j=1,n
    imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))         ! Find the pivot row.
    if (j /= imax) then                                 ! Do we need to interchange rows?
        call swap(a(imax,:),a(j,:))                     ! Yes, do so...
        d=-d                                            ! ...and change the parity of d.
        vv(imax)=vv(j)                                  ! Also interchange the scale factor.
    end if
    indx(j)=imax
    if (a(j,j) == 0.0) a(j,j)=TINY                      !If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
                                                        !For some applications on singular matrices, it is desirable to substitute TINY for zero.
    
    a(j+1:n,j)=a(j+1:n,j)/a(j,j)                        ! Divide by the pivot element.        
    a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))  ! Reduce remaining submatrix.

enddo

End subroutine LU_Decomposition
!******************************************************************************************************************
SUBROUTINE LU_BackSubstitution(N,A,LDA,Indx,B)

IMPLICIT Real(8) (A-H,O-Z)

Real(8) A(LDA,N),B(N)
Integer(4) Indx(N)

!Solves the set of N linear equations A ú X = B. 
! A     - LU-decomposed equations matrix
! Indx  - permutation vector after LU-decomposition
! B     - On input: right-hand-side vector B. On exit: solution vector X 
!
! A and Indx are not modified and can be left in place for successive calls with different right-hand sides b. 
! B can begin with many zero elements (ii>0), so it is efficient for use in matrix inversion.

! If ii is set to a positive value, it is the index of the first nonvanishing element of b. 
ii=0    

! Forward substitution with permutations
do i=1,n
    ll=indx(i)
    summ=b(ll)
    b(ll)=b(i)
    if (ii > 0) then
        summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
    else if (summ /= 0.0) then
        ii=i        !A nonzero element was encountered, so from now on we will have to do the dot product above.
    Endif
    b(i)=summ
enddo

!Back substitution
do i=n,1,-1       
    b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
end do

End subroutine LU_BackSubstitution


End module LUdecomposition
!***********************************************************************************************    
