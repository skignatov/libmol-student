      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER*4 i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))

14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+ &
     g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END
!***********************************************************************************
Subroutine JacobiSorted(A,N,NP,E,V,NRrot,iSort)
Implicit Real(8) (A-H,O-Z)

Real(8) A(NP,NP),E(NP),V(NP,NP)
Integer(4) iOrder(np),Es(NP),Vs(NP,NP)

! JacobiSorted gives the sorted Eigenvalues and Eigenvectors of real symmetric matrix 
!   iSort>0 - EV are sorted in Ascending order
!   iSort<0 - EV are sorted in Descending order
!   iSort=0 - EV are not sorted

Call Jacobi(A,N,NP,E,V,NRrot)

If (iSort==0) Return

If (iSort>0) Then                   ! Ascending Order
	Call HeapSort(NP,E,V,iOrder)
Else If (iSort<0) Then              ! Descending order
	E=-E
	Call HeapSort(NP,E,V,iOrder)
	E=-E
Endif
	

End	
!*************************************************************************
Subroutine HeapSort(n,ra,rb,ia)
Implicit Real(8) (A-H,O-Z)

Real(8) ra(n),rb(n,n),tmp(n,n)
Integer(4) ia(n)

! HeapSort algorithm from "Numerical Recipes in F90"
! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is
! input; ra is replaced on output by its sorted rearrangement.

! rb the second nxn array with columns sorted as ra. ia - indexes of permutations

! CITED REFERENCES AND FURTHER READING:
! Knuth, D.E. 1973, Sorting and Searching, vol. 3 of The Art of Computer Programming (Reading,
! MA: Addison-Wesley), x5.2.3. [1]
! Sedgewick, R. 1988, Algorithms, 2nd ed. (Reading, MA: Addison-Wesley), Chapter 11. [2]
! 8.4 Indexing and Ranking

if (n.lt.2) return

Do i=1,n
	ia(i)=i
Enddo

!The index l will be decremented from its initial value down to 1 during the \hiring" (heap
!creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
!down to 1 during the \retirement-and-promotion" (heap selection) phase.

l=n/2+1
ir=n
10 continue
if(l.gt.1)then		!Still in hiring phase.
	l=l-1
	rra=ra(l)
	iia=ia(l)			!*
else				!In retirement-and-promotion phase.
	rra=ra(ir)		!Clear a space at end of array.
	iia=ia(ir)			!*
	ra(ir)=ra(1)	!Retire the top of the heap into it.
	ia(ir)=ia(1)		!*
	ir=ir-1			!Decrease the size of the corporation.
	if(ir.eq.1)then	!Done with the last promotion.
		ra(1)=rra	!The least competent worker of all!
		ia(1)=iia		!*
		Goto 30
	endif
endif

i=l					!Whether in the hiring phase or promotion phase, we here
					!set up to sift down element 
j=l+l				!rra to its proper level.
20 if(j.le.ir)then	!Do while j.le.ir:"
if(j.lt.ir)then
	if(ra(j).lt.ra(j+1))j=j+1		!Compare to the better underling.
endif
if(rra.lt.ra(j))then				!Demote rra.
	ra(i)=ra(j)
	ia(i)=ia(j)		!*
	i=j
	j=j+j
else								!This is rra's level. Set j to terminate the sift-down.
	j=ir+1
endif
goto 20
endif
	
ra(i)=rra			!Put rra into its slot.
ia(i)=iia			!*
goto 10


! Permutations of the second aray
30 Continue
Do i=1,n
	j=ia(i)
	tmp(1:n,i)=rb(1:n,j)
Enddo
rb=tmp

END










