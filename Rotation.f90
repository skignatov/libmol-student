Module Rotations

Implicit Real(8) (A-H,O-Z)

Real(8), parameter:: Pi=3.141592653589793d0
Real(8), parameter:: deg2rad=Pi/180.d0
Real(8), parameter:: rad2deg=1.d0/deg2rad


CONTAINS
! Subroutines for rotation parameters conversion:
! Subroutine ABG2QN(mode,irad,alp,bet,gam,Q)    - Euler angles (alp,bet,gam)    <-->    Quaternion Q(0:3)
! Subroutine RM2QN(mode,U,Q)                    - Rotation matrix U(3,3)        <-->    Quaternion Q(0:3)
! Subroutine QN2WA(mode,irad,Q,W,angle)         - Quaternion Q(0:3)             <-->    Rotation axis W(3) + angle
! Subroutine RM2WA(mode,irad,U,W,angle)         - Rotation matrix U(3,3)        <-->    Rotation axis W(3) + angle
! Subroutine ABG2RM(mode,irad,alp,bet,gam,U)    - Euler angles (alp,bet,gam)    <-->    Rotation matrix U(3,3)
! Subroutine ABG2WA(mode,irad,alp,bet,gam,W,angle) - Euler angles (alp,bet,gam) <-->    Rotation axis W(3) + angle
!**************************************************************
Subroutine ABG2QN(mode,irad,alp,bet,gam,Q)
Implicit Real(8) (A-H,O-Z)
   
Real(8) Q(0:3)

! Subroutine ABG2QN converts Euler angles alp,bet,gam (active ZXZ definition) to quaternion Q (mode>=0) and back (mode<0)
! Angles in degrees (irad=0) or in radians (irad=1)
! From https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions

q0=Q(0)
q1=Q(1)
q2=Q(2)
q3=Q(3)

conv=1.d0
If (mode>=0) Then
    If (irad==0) conv=deg2rad
    a2=alp*conv*0.5d0
    b2=bet*conv*0.5d0
    a2=gam*conv*0.5d0
    sb=DSIN(b2)
    cb=DCOS(b2)
    Q(1)=DCOS(a2-g2)*sb
    Q(2)=DSIN(a2-g2)*sb
    Q(3)=DSIN(a2+g2)*cb
    Q(0)=DCOS(a2+g2)*cb
Else
    If (irad==0) conv=rad2deg
    alp=DATAN2(2.d0*(q0*q1+q2*q3),1.d0-2.d0*(q1*q1-q2*q2)) *conv
    bet=DASIN(2.d0*(q0*q2-q3*q1))                          *conv
    gam=DATAN2(2.d0*(q0*q3+q1*q2),1.d0-2.d0*(q2*q2+q3*q3)) *conv
Endif
    
End    
!**************************************************************
Subroutine RM2QN(mode,U,Q)

Implicit Real(8) (A-H,O-Z)

Real(8) U(3,3),Q(0:3)
! Subroutine RM2QM converts rotation matrix U to quatrenion Q (mode>=0) and back (mode<0)

If (mode>=0) Then
    trace=U(1,1)+U(2,2)+U(3,3)
    If (trace>1.d-8) Then
        Q(0)=0.5d0*DSQRT(1.d0+U(1,1)+U(2,2)+U(3,3))
        factor=1.d0/(4.d0*Q(0))
        Q(1)=factor*(U(3,2)-U(2,3))
        Q(2)=factor*(U(1,3)-U(3,1))
        Q(3)=factor*(U(2,2)-U(1,2))
    Else
        Q(1)=0.5d0*DSQRT(1.d0+U(1,1)-U(2,2)-U(3,3))
        factor=1.d0/(4.d0*Q(1))
        Q(2)=factor*(U(1,2)+U(2,1))
        Q(3)=factor*(U(1,3)+U(3,1))
        Q(0)=factor*(U(3,2)+U(2,3))
    Endif
Else
    qi=Q(1)
    qj=Q(2)
    qk=Q(3)
    qr=Q(0)
    U(1,1)=0.5d0-qj*qj-qk*qk
    U(2,2)=0.5d0-qi*qi-qk*qk
    U(3,3)=0.5d0-qi*qi-qj*qj
    U(1,2)=qi*qj-qk*qr
    U(1,3)=qi*qk+qj*qr
    U(2,1)=qi*qj+qk*qr
    U(2,3)=qj*qk-qi*qr
    U(3,1)=qi*qk-qj*qr
    U(3,2)=qi*qr+qj*qk
    U=U*2.d0
Endif

End
!**************************************************************
Subroutine QN2WA(mode,irad,Q,W,angle)

Implicit Real(8) (A-H,O-Z)

Real(8) U(3,3),W(3),Q(0:3)
! Subroutine QN2WA converts quaterbnion Q to (rotation vector W and angle) (mode>=0) and back (mode<0)

If (mode>=0) Then
    qn=Sum(Q(1:3)**2)
    qn=1.d0/DSQRT(qn)
    W(1:3)=Q(1:3)*qn
    angle=2.d0*DACOS(Q(0))
    If (irad/=1) angle=angle*rad2deg
Else
    ang=angle
    If (irad/=1) ang=ang*rad2deg
    sa2=DSIN(ang*0.5d0)
    Q(1:3)=W(1:3)*sa2
    Q(0)=DCOS(ang*0.5d0)
Endif

End
!**************************************************************
Subroutine RM2WA(mode,irad,U,W,angle)

Implicit Real(8) (A-H,O-Z)

Real(8) U(3,3),W(3),Q(0:3)
! Subroutine RM2WA converts rotation matrix U to (rotation vector W and angle) (mode>=0) and back (mode<0)
! angle in radians if irad==1 and in degrees otherwise

If (mode>=0) Then
    Call RM2QN(1,U,Q)
    Call QN2WA(1,irad,Q,W,angle)
Else
    e1=W(1)
    e2=W(2)
    e3=W(3)
    ang=angle
    If (irad/=1) ang=ang*deg2rad
    ca=DCOS(ang)
    sa=DSIN(ang)
    c1=1.d0-ca
    U(1,1)=c1*e1*e1+ca
    U(2,2)=c1*e2*e2+ca
    U(3,3)=c1*e3*e3+ca
    U(1,2)=c1*e1*e2-e3*sa
    U(1,3)=c1*e1*e3+e2*sa
    U(2,1)=c1*e2*e1+e3*sa
    U(2,3)=c1*e2*e3-e1*sa
    U(3,1)=c1*e3*e1-e2*sa
    U(3,2)=c1*e3*e2+e1*sa
Endif

End
!**************************************************************
Subroutine ABG2RM(mode,irad,alp,bet,gam,U)

Implicit Real(8) (A-H,O-Z)

Real(8) U(3,3),W(3),Q(0:3)

If (mode>=0) Then
    Call ABG2QN(1,irad,alp,bet,gam,Q)
    Call RM2QN(-1,U,Q)
Else
    Call RM2QN(1,U,Q)
    Call ABG2QN(-1,irad,alp,bet,gam,Q)
Endif

End
!**************************************************************
Subroutine ABG2WA(mode,irad,alp,bet,gam,W,angle)

Implicit Real(8) (A-H,O-Z)

Real(8) U(3,3),W(3),Q(0:3)

If (mode>=0) Then
    Call ABG2QN(1,irad,alp,bet,gam,Q)
    Call QN2WA(1,irad,Q,W,angle)
Else
    Call QN2WA(-1,irad,Q,W,angle)
    Call ABG2QN(-1,irad,alp,bet,gam,Q)
Endif    

End
   

End module

