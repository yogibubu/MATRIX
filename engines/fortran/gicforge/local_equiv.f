C MATRIX GICForge local-equivalence helpers.
C Strict Fortran77-compatible kernels used to align the Fortran path with
C Python NEO local stretch/angle SALC classification.
C
C These routines do not build final GIC rows.  They classify local primitives
C before the existing local SVD/Jacobi machinery is called.

      Subroutine ORCLBND(NBOND,IAT,JAT,C,ZEFF,CLASS,NCLASS)
      Implicit Real*8(A-H,O-Z)
      Integer NBOND,IAT(*),JAT(*),CLASS(*),NCLASS
      Dimension C(3,*),ZEFF(*)
      Dimension K1(200),K2(200),KD(200)
      Data TZ/5.0D-4/,TD/1.0D-3/

      NCLASS=0
      Do 100 I=1,NBOND
       Z1=ZEFF(IAT(I))
       Z2=ZEFF(JAT(I))
       If(Z1.gt.Z2) then
        TMP=Z1
        Z1=Z2
        Z2=TMP
       EndIf
       D=DSQRT((C(1,IAT(I))-C(1,JAT(I)))**2+
     $         (C(2,IAT(I))-C(2,JAT(I)))**2+
     $         (C(3,IAT(I))-C(3,JAT(I)))**2)
       CLASS(I)=0
       Do 20 K=1,NCLASS
        If(DABS(Z1-DBLE(K1(K))*TZ).le.TZ.and.
     $     DABS(Z2-DBLE(K2(K))*TZ).le.TZ.and.
     $     DABS(D-DBLE(KD(K))*TD).le.TD) then
         CLASS(I)=K
         Goto 30
        EndIf
   20  Continue
   30  Continue
       If(CLASS(I).eq.0) then
        NCLASS=NCLASS+1
        CLASS(I)=NCLASS
        K1(NCLASS)=NINT(Z1/TZ)
        K2(NCLASS)=NINT(Z2/TZ)
        KD(NCLASS)=NINT(D/TD)
       EndIf
  100 Continue
      Return
      End

      Subroutine ORCLLIG(CENTER,NNEI,NEI,C,ZEFF,CLASS,NCLASS)
      Implicit Real*8(A-H,O-Z)
      Integer CENTER,NNEI,NEI(*),CLASS(*),NCLASS
      Dimension C(3,*),ZEFF(*)
      Dimension KZ(20),KD(20)
      Data TZ/5.0D-4/,TD/1.0D-3/

      NCLASS=0
      Do 100 I=1,NNEI
       IA=NEI(I)
       D=DSQRT((C(1,IA)-C(1,CENTER))**2+
     $         (C(2,IA)-C(2,CENTER))**2+
     $         (C(3,IA)-C(3,CENTER))**2)
       CLASS(I)=0
       Do 20 K=1,NCLASS
        If(DABS(ZEFF(IA)-DBLE(KZ(K))*TZ).le.TZ.and.
     $     DABS(D-DBLE(KD(K))*TD).le.TD) then
         CLASS(I)=K
         Goto 30
        EndIf
   20  Continue
   30  Continue
       If(CLASS(I).eq.0) then
        NCLASS=NCLASS+1
        CLASS(I)=NCLASS
        KZ(NCLASS)=NINT(ZEFF(IA)/TZ)
        KD(NCLASS)=NINT(D/TD)
       EndIf
  100 Continue
      Return
      End

      Subroutine ORCLTPL(CENTER,NNEI,NEI,C,ITPL,SCORE)
      Implicit Real*8(A-H,O-Z)
      Integer CENTER,NNEI,NEI(*),ITPL
      Dimension C(3,*),A(36),T(36)

      Call ORCLCOS(CENTER,NNEI,NEI,C,A)
      ITPL=0
      SCORE=1.0D30
      Do 100 JTPL=1,2
       Call ORCLTCOS(NNEI,JTPL,T,NT)
       If(NT.ne.NNEI*(NNEI-1)/2) Goto 100
       S=0.0D0
       Do 20 I=1,NT
        S=S+(A(I)-T(I))**2
   20  Continue
       S=DSQRT(S/DBLE(NT))
       If(S.lt.SCORE) then
        SCORE=S
        ITPL=JTPL
       EndIf
  100 Continue
      If(SCORE.gt.1.2D-1) ITPL=0
      Return
      End

      Subroutine ORCLCOS(CENTER,NNEI,NEI,C,A)
      Implicit Real*8(A-H,O-Z)
      Integer CENTER,NNEI,NEI(*)
      Dimension C(3,*),A(*),V(3,20)
      K=0
      Do 20 I=1,NNEI
       IA=NEI(I)
       RN=0.0D0
       Do 10 M=1,3
        V(M,I)=C(M,IA)-C(M,CENTER)
        RN=RN+V(M,I)**2
   10  Continue
       RN=DSQRT(RN)
       If(RN.gt.0.0D0) then
        Do 15 M=1,3
         V(M,I)=V(M,I)/RN
   15   Continue
       EndIf
   20 Continue
      Do 40 I=1,NNEI-1
       Do 30 J=I+1,NNEI
        K=K+1
        A(K)=V(1,I)*V(1,J)+V(2,I)*V(2,J)+V(3,I)*V(3,J)
   30  Continue
   40 Continue
      Call ORCLSORT(K,A)
      Return
      End

      Subroutine ORCLTCOS(NNEI,JTPL,A,N)
      Implicit Real*8(A-H,O-Z)
      Integer NNEI,JTPL,N
      Dimension A(*),V(3,20)
      N=0
      If(NNEI.eq.5.and.JTPL.eq.1) then
C      trigonal bipyramidal
       Call ORCLPOL(3,0.0D0,0.0D0,V)
       V(1,4)=0.0D0
       V(2,4)=0.0D0
       V(3,4)=1.0D0
       V(1,5)=0.0D0
       V(2,5)=0.0D0
       V(3,5)=-1.0D0
      ElseIf(NNEI.eq.6.and.JTPL.eq.1) then
C      octahedral
       Call ORCLSET(V,1, 1.0D0, 0.0D0, 0.0D0)
       Call ORCLSET(V,2,-1.0D0, 0.0D0, 0.0D0)
       Call ORCLSET(V,3, 0.0D0, 1.0D0, 0.0D0)
       Call ORCLSET(V,4, 0.0D0,-1.0D0, 0.0D0)
       Call ORCLSET(V,5, 0.0D0, 0.0D0, 1.0D0)
       Call ORCLSET(V,6, 0.0D0, 0.0D0,-1.0D0)
      ElseIf(NNEI.eq.7.and.JTPL.eq.1) then
C      pentagonal bipyramidal
       Call ORCLPOL(5,0.0D0,0.0D0,V)
       Call ORCLSET(V,6,0.0D0,0.0D0,1.0D0)
       Call ORCLSET(V,7,0.0D0,0.0D0,-1.0D0)
      ElseIf(NNEI.eq.8.and.JTPL.eq.1) then
C      square antiprismatic
       Call ORCLPOL(4,0.45D0,0.0D0,V)
       Call ORCLPOL2(4,-0.45D0,0.7853981633974483D0,V,5)
      ElseIf(NNEI.eq.9.and.JTPL.eq.1) then
C      tricapped trigonal prism
       Call ORCLPOL(3,0.58D0,0.0D0,V)
       Call ORCLPOL2(3,-0.58D0,0.0D0,V,4)
       Call ORCLPOL2(3,0.0D0,1.0471975511965977D0,V,7)
      Else
       Return
      EndIf
      K=0
      Do 40 I=1,NNEI-1
       Do 30 J=I+1,NNEI
        K=K+1
        A(K)=V(1,I)*V(1,J)+V(2,I)*V(2,J)+V(3,I)*V(3,J)
   30  Continue
   40 Continue
      N=K
      Call ORCLSORT(N,A)
      Return
      End

      Subroutine ORCLSET(V,I,X,Y,Z)
      Implicit Real*8(A-H,O-Z)
      Dimension V(3,*)
      V(1,I)=X
      V(2,I)=Y
      V(3,I)=Z
      Return
      End

      Subroutine ORCLPOL(N,Z,PHASE,V)
      Implicit Real*8(A-H,O-Z)
      Dimension V(3,*)
      Call ORCLPOL2(N,Z,PHASE,V,1)
      Return
      End

      Subroutine ORCLPOL2(N,Z,PHASE,V,IOFF)
      Implicit Real*8(A-H,O-Z)
      Dimension V(3,*)
      PI=DACOS(-1.0D0)
      R=DSQRT(DMAX1(0.0D0,1.0D0-Z*Z))
      Do 10 I=1,N
       T=PHASE+2.0D0*PI*DBLE(I-1)/DBLE(N)
       V(1,IOFF+I-1)=R*DCOS(T)
       V(2,IOFF+I-1)=R*DSIN(T)
       V(3,IOFF+I-1)=Z
   10 Continue
      Return
      End

      Subroutine ORCLSORT(N,A)
      Implicit Real*8(A-H,O-Z)
      Dimension A(*)
      Do 20 I=1,N-1
       IM=I
       Do 10 J=I+1,N
        If(A(J).lt.A(IM)) IM=J
   10  Continue
       If(IM.ne.I) then
        T=A(I)
        A(I)=A(IM)
        A(IM)=T
       EndIf
   20 Continue
      Return
      End
