

CASE (3)                ! Lennard-Jonnes potential
    Do i=1,Numat
        nai=NA(i)
        Do j=1,i-1
            naj=NA(j)

            rij=Distance(i,j,MaxAt,C)

            rCut=rCutPar(1,ip)                          ! Cutoff distance
            If (rCut>0.d0.and.rij>rCut) Cycle		! Do not consider potential if atoms are far from each another
            nprs=nprs+1

            eps=Xpar(1,ip)	! Take epsylon and sigma from table of atomic parameters
            sgm=Xpar(2,ip)
        
            rr=(sgm/rij)        ! Calculate energy
            rr2=rr*rr
            rr6=rr2*rr2*rr2
            Eattr=-4.d0*eps*rr6
            Erep=-Eattr*rr6
            fij=Erep+Eattr
            F=F+fij
            
				! calculate gradients analytically
            If (iGnumeric==0) Then            ! Analytic gradient components
                Gr=24.d0*eps/(sgm*rij)*rr*(2.d0*rr6*rr6-rr6)      ! G(r)=-(1/r)(dU/dr)
                dCij(1:3)=Gr*(C(1:3,i)-C(1:3,j))
                dFdC(1:3,i)=dFdC(1:3,i)-dCij(1:3)
                dFdC(1:3,j)=dFdC(1:3,j)+dCij(1:3)
            Endif

            
        Enddo
    Enddo

