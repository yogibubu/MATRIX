      PROGRAM TEST_KRAITCHMAN
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NMAX=50)
        INTEGER NAtoms
        REAL*8 C(3,NMAX), CS(3,NMAX)
        REAL*8 AtMass0(NMAX), AtMass(NMAX)
        REAL*8 DImat(3,3), CISub(3), XYprod(3)
        REAL*8 MASSFAC, DeltaM
        INTEGER I,J
        REAL*8 dB(3)
C--------------------------------------------------------------
        OPEN(10, FILE='newkra.inp', STATUS='OLD')
        READ(10,*) NAtoms
        READ(10,*) (AtMass0(I), I=1,NAtoms)
        READ(10,*) (AtMass(I), I=1,NAtoms)
        READ(10,*) MASSFAC
        DO I=1,3
          READ(10,*) (C(I,J), J=1,NAtoms)
        END DO
        CLOSE(10)
C Example: assign 1σ uncertainties in MHz (typical)
        dB(1)=0.005D0
        dB(2)=0.005D0
        dB(3)=0.005D0
C--------------------------------------------------------------
C duplicate for isotopologue
        DO I=1,NAtoms
          DO J=1,3
            CS(J,I)=C(J,I)
          END DO
        END DO

        DeltaM = AtMass(2) - AtMass0(2)

        CALL K_DeltaI_EckartAware(NAtoms, MASSFAC,
     &       AtMass0, AtMass, C, CS, DImat)

        CALL K_KraitchmanFromDeltaI(DImat, DeltaM, CISub, XYprod)

        WRITE(*,*) 'ΔI matrix (amu*Å^2):'
        DO I=1,3
          WRITE(*,'(3F15.8)') (DImat(I,J),J=1,3)
        END DO
        WRITE(*,*) 'Substitution coordinates (|a|,|b|,|c|) (Å):'
        WRITE(*,'(3F12.6)') (CISub(I),I=1,3)
        WRITE(*,*) 'Cross products (-ΔIxy/Δm etc):'
        WRITE(*,'(3F12.6)') (XYprod(I),I=1,3)

        CALL Print_PROSPE_Output_EXTENDED(DImat, DeltaM, CISub)

        CALL Print_CostainErrors(DImat, dB, DeltaM, CISub)

C       CALL Compute_rm_Structure(Niso, NAtoms, MassIso, Rs, Rm)

        END

C===============================================================
C Helper: transpose inline data to column-wise (3xN)
C===============================================================
      SUBROUTINE TRANSPOSE3(A)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 A(3,3), TMP(3,3)
        INTEGER I,J
        DO I=1,3
          DO J=1,3
            TMP(I,J)=A(J,I)
          END DO
        END DO
        DO I=1,3
          DO J=1,3
            A(I,J)=TMP(I,J)
          END DO
        END DO
        RETURN
      END

C===============================================================
C  Shift to COM and rotate coordinates to their own principal axes
C===============================================================
      SUBROUTINE ROTATE_TO_PA(N, MASSFAC, MASS, CIN, COUT, Ldiag, R)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER N
        REAL*8 MASSFAC, MASS(*), CIN(3,*), COUT(3,*)
        REAL*8 Ldiag(3,3), R(3,3)
        REAL*8 CTMP(3,1000), TI(3,3)
        INTEGER IAT,I,J
C-- copy
        DO IAT=1,N
          CTMP(1,IAT)=CIN(1,IAT)
          CTMP(2,IAT)=CIN(2,IAT)
          CTMP(3,IAT)=CIN(3,IAT)
        END DO
C-- center to its own COM
        CALL COM_SHIFT(N, MASSFAC, MASS, CTMP)
C-- inertia in current frame
        CALL INERTIA_TENSOR(N, MASSFAC, MASS, CTMP, TI)
C-- diagonalize -> principal moments Ldiag and eigenvectors R
        CALL JACOBI_SYM3(TI, Ldiag, R)
C-- rotate coordinates into their own PA: COUT = R^T * CTMP
        DO IAT=1,N
          DO I=1,3
            COUT(I,IAT)= R(1,I)*CTMP(1,IAT)
     &                  +R(2,I)*CTMP(2,IAT)
     &                  +R(3,I)*CTMP(3,IAT)
          END DO
        END DO
        RETURN
      END

C===============================================================
C  Build ΔI in the ref-PA frame, accounting for different Eckart frames
C===============================================================
      SUBROUTINE K_DeltaI_EckartAware(NAtoms, MASSFAC, AtMass0, AtMass,
     &                               C, CS, DImat)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NAtoms
        REAL*8 MASSFAC, AtMass0(*), AtMass(*)
        REAL*8 C(3,*), CS(3,*), DImat(3,3)
        REAL*8 C0PA(3,1000), CiPA(3,1000)
        REAL*8 L0(3,3), Li(3,3), R0(3,3), Ri(3,3)
        REAL*8 TMP(3,3), Rir0t(3,3)

C-- Put REF (C, AtMass0) in its own PA
        CALL ROTATE_TO_PA(NAtoms, MASSFAC, AtMass0, C,  C0PA, L0, R0)
C-- Put ISO (CS, AtMass) in its own PA (even if CS==C initially)
        CALL ROTATE_TO_PA(NAtoms, MASSFAC, AtMass,  CS, CiPA, Li, Ri)

C-- Rotate ISO principal moments into REF-PA:  Rir0t = Ri * R0^T
        CALL MATMUL3(Ri, R0, Rir0t, .TRUE.)
C   TMP = Rir0t * Li
        CALL MATMUL3(Rir0t, Li, TMP, .FALSE.)
C   Iiso_in_ref = TMP * Rir0t^T  ;  ΔI = Iiso_in_ref - L0
        CALL MATMUL3(TMP, Rir0t, DImat, .TRUE.)
        DImat(1,1)=DImat(1,1)-L0(1,1)
        DImat(2,2)=DImat(2,2)-L0(2,2)
        DImat(3,3)=DImat(3,3)-L0(3,3)
        DImat(1,2)=DImat(1,2)-L0(1,2)
        DImat(2,1)=DImat(1,2)
        DImat(1,3)=DImat(1,3)-L0(1,3)
        DImat(3,1)=DImat(1,3)
        DImat(2,3)=DImat(2,3)-L0(2,3)
        DImat(3,2)=DImat(2,3)
        RETURN
      END
C===============================================================
C  Compute Kraitchman coordinates (|a|,|b|,|c|) and signed (a,b,c)
C  from full DeltaI matrix and atomic mass difference Δm
C===============================================================
      SUBROUTINE K_KraitchmanFromDeltaI(DI, DeltaM, CISub, XYprod)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 DI(3,3), DeltaM, CISub(3), XYprod(3)
        REAL*8 a2, b2, c2, sgn(3), small
        INTEGER I

        small = 1D-10
        IF (DABS(DeltaM) .LT. small) THEN
          WRITE(*,*) '*** ERROR: DeltaM too small in Kraitchman'
          RETURN
        ENDIF

C-- standard Kraitchman absolute values
        a2 =  DI(2,2)/DABS(DeltaM)
        b2 =  DI(1,1)/DABS(DeltaM)
        c2 =  DI(3,3)/DABS(DeltaM)

        IF (a2 .LT. 0D0) a2=0D0
        IF (b2 .LT. 0D0) b2=0D0
        IF (c2 .LT. 0D0) c2=0D0

        CISub(1)=DSQRT(a2)
        CISub(2)=DSQRT(b2)
        CISub(3)=DSQRT(c2)

C-- cross-products (used for sign determination)
        XYprod(1) = -DI(1,2)/DeltaM
        XYprod(2) = -DI(1,3)/DeltaM
        XYprod(3) = -DI(2,3)/DeltaM

C-- determine signs using cross-product relations
        sgn(1)=1D0
        sgn(2)=1D0
        sgn(3)=1D0

        IF (XYprod(1) .LT. 0D0) sgn(1)=-sgn(2)
        IF (XYprod(2) .LT. 0D0) sgn(1)=-sgn(3)
        IF (XYprod(3) .LT. 0D0) sgn(2)=-sgn(3)

C-- apply signs
        DO I=1,3
          CISub(I)=CISub(I)*sgn(I)
        END DO

        RETURN
      END
C===============================================================
C  Shift coordinates to their own center of mass
C===============================================================
      SUBROUTINE COM_SHIFT(N, MASSFAC, M, C)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER N
        REAL*8 MASSFAC, M(*), C(3,*)
        REAL*8 MTOT, COM1, COM2, COM3
        INTEGER I

        MTOT = 0D0
        COM1 = 0D0
        COM2 = 0D0
        COM3 = 0D0
        DO 10 I=1,N
          MTOT = MTOT + MASSFAC*M(I)
          COM1 = COM1 + MASSFAC*M(I)*C(1,I)
          COM2 = COM2 + MASSFAC*M(I)*C(2,I)
          COM3 = COM3 + MASSFAC*M(I)*C(3,I)
10      CONTINUE
        IF (MTOT .EQ. 0D0) RETURN
        COM1 = COM1/MTOT
        COM2 = COM2/MTOT
        COM3 = COM3/MTOT
        DO 20 I=1,N
          C(1,I) = C(1,I) - COM1
          C(2,I) = C(2,I) - COM2
          C(3,I) = C(3,I) - COM3
20      CONTINUE
        RETURN
      END

C===============================================================
C  Full 3x3 inertia tensor (amu*Ang^2) from point masses
C===============================================================
      SUBROUTINE INERTIA_TENSOR(N, MASSFAC, M, C, I)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER N
        REAL*8 MASSFAC, M(*), C(3,*), I(3,3)
        REAL*8 MI, X, Y, Z, R2
        INTEGER K

        I(1,1)=0D0
        I(1,2)=0D0
        I(1,3)=0D0
        I(2,1)=0D0
        I(2,2)=0D0
        I(2,3)=0D0
        I(3,1)=0D0
        I(3,2)=0D0
        I(3,3)=0D0

        DO 30 K=1,N
          MI = MASSFAC*M(K)
          X = C(1,K)
          Y = C(2,K)
          Z = C(3,K)
          R2 = X*X + Y*Y + Z*Z
          I(1,1) = I(1,1) + MI*(R2 - X*X)
          I(2,2) = I(2,2) + MI*(R2 - Y*Y)
          I(3,3) = I(3,3) + MI*(R2 - Z*Z)
          I(1,2) = I(1,2) - MI*(X*Y)
          I(1,3) = I(1,3) - MI*(X*Z)
          I(2,3) = I(2,3) - MI*(Y*Z)
30      CONTINUE
        I(2,1) = I(1,2)
        I(3,1) = I(1,3)
        I(3,2) = I(2,3)
        RETURN
      END

C===============================================================
C  Jacobi diagonalization for real symmetric 3x3 matrix
C  Returns D = diag(eigvals), R = eigenvectors (columns)
C===============================================================
      SUBROUTINE JACOBI_SYM3(A, D, R)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 A(3,3), D(3,3), R(3,3)
        REAL*8 B12, B13, B23
        INTEGER I, J, ITER

C-- initialise R=I, D=A
        DO 110 I=1,3
          DO 100 J=1,3
            R(I,J) = 0D0
            D(I,J) = A(I,J)
100       CONTINUE
          R(I,I) = 1D0
110     CONTINUE

C-- up to 50 sweeps (enough for 3x3)
        DO 200 ITER=1,50
          B12 = DABS(D(1,2))
          B13 = DABS(D(1,3))
          B23 = DABS(D(2,3))
          IF (B12+B13+B23 .LT. 1D-20) GOTO 210
          CALL JACOBI_ROT(D, R, 1, 2)
          CALL JACOBI_ROT(D, R, 1, 3)
          CALL JACOBI_ROT(D, R, 2, 3)
200     CONTINUE
210     CONTINUE

C-- zero off-diagonal explicitly
        DO 220 I=1,3
          DO 215 J=1,3
            IF (I .NE. J) D(I,J)=0D0
215       CONTINUE
220     CONTINUE
        RETURN
      END

C===============================================================
C  One Jacobi rotation on (p,q) for symmetric 3x3
C===============================================================
      SUBROUTINE JACOBI_ROT(D, R, P, Q)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 D(3,3), R(3,3)
        INTEGER P, Q, K
        REAL*8 APP, AQQ, APQ, PHI, C, S, TMP

        IF (DABS(D(P,Q)) .LT. 1D-30) RETURN
        APP = D(P,P)
        AQQ = D(Q,Q)
        APQ = D(P,Q)
        PHI = 0.5D0 * DATAN2( 2D0*APQ, (AQQ-APP) )
        C = DCOS(PHI)
        S = DSIN(PHI)

C-- rotate D = G^T D G (update rows/cols p,q)
        DO 310 K=1,3
          TMP   =  C*D(K,P) + S*D(K,Q)
          D(K,Q)= -S*D(K,P) + C*D(K,Q)
          D(K,P)=  TMP
310     CONTINUE
        DO 320 K=1,3
          TMP   =  C*D(P,K) + S*D(Q,K)
          D(Q,K)= -S*D(P,K) + C*D(Q,K)
          D(P,K)=  TMP
320     CONTINUE
        D(P,Q)=0D0
        D(Q,P)=0D0

C-- accumulate eigenvectors R = R * G
        DO 330 K=1,3
          TMP   =  C*R(K,P) + S*R(K,Q)
          R(K,Q)= -S*R(K,P) + C*R(K,Q)
          R(K,P)=  TMP
330     CONTINUE
        RETURN
      END

C===============================================================
C  3x3 matrix multiply: C = A * (B or B^T)
C  If BTRANS=.TRUE. then second operand is B^T
C===============================================================
      SUBROUTINE MATMUL3(A, B, C, BTRANS)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 A(3,3), B(3,3), C(3,3)
        LOGICAL BTRANS
        INTEGER I, J, K
        REAL*8 S

        DO 430 I=1,3
          DO 420 J=1,3
            S = 0D0
            IF (.NOT. BTRANS) THEN
              DO 410 K=1,3
                S = S + A(I,K)*B(K,J)
410           CONTINUE
            ELSE
              DO 415 K=1,3
                S = S + A(I,K)*B(J,K)
415           CONTINUE
            END IF
            C(I,J) = S
420       CONTINUE
430     CONTINUE
        RETURN
      END

C=======================================================================
C  PROSPE-STYLE OUTPUT COMPATIBLE WITH TEST_KRAITCHMAN
C=======================================================================
      SUBROUTINE Print_PROSPE_Output_EXTENDED(DImat, DeltaM, CISub)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DImat(3,3), CISub(3)
      REAL*8 H, CSPD, CONV, PI, FACT
      REAL*8 Ia, Ib, Ic, A, B, C, Iav, Delt
      INTEGER I

C-- Physical constants
      H     = 6.62607015D-34
      CSPD  = 2.99792458D8
      CONV  = 1.66053906660D-47
      PI    = 4.D0*DATAN(1.D0)
      FACT  = H / (8.D0*PI*PI*CSPD*CONV*1.D6)

C-- Extract diagonal ΔI (in amu·Å²)
      Ia = DImat(1,1)
      Ib = DImat(2,2)
      Ic = DImat(3,3)
      Delt = Ic - Ia - Ib
      Iav  = (Ia + Ib + Ic)/3.D0

C-- Print PROSPE-style block
      WRITE(*,*)
      WRITE(*,*) '--- PROSPE-style results ---'
      WRITE(*,'(A,3F14.6)') 'Moments of inertia (amu·Å²):', Ia,Ib,Ic
      WRITE(*,'(A,F14.6)') 'Inertial defect Δ (amu·Å²) =', Delt
      WRITE(*,'(A,F14.6)') 'Average moment Iav (amu·Å²) =', Iav
      WRITE(*,*)
      WRITE(*,'(A,3F14.6)') 'ΔI (amu·Å²):', Ia,Ib,Ic
C     WRITE(*,'(A,3F14.3)') 'Equivalent ΔB (MHz):', A,B,C
C     WRITE(*,'(A,3F14.3)') 'Rotational constants (MHz):', A,B,C
      WRITE(*,*) '-------------------------------------------'
      RETURN
      END

C=======================================================================
C  COMPUTE RM STRUCTURE FROM RS (Watson linear correction)
C  Fits r_m = r_s + A*(m-1) + B*(m^{-1/2}-1)
C=======================================================================
      SUBROUTINE Compute_rm_Structure(Niso,NAtoms,MassIso,Rs,Rm)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER Niso,NAtoms,J,K,I
      DIMENSION MassIso(Niso,NAtoms), Rs(Niso,NAtoms,3), Rm(NAtoms,3)
      REAL*8 A(3), B(3), M, W1, W2
      REAL*8 S1,S2,S12,S22,SY1,SY2,Det

C-- Fit Watson’s correction for each atom and coordinate
      DO J=1,NAtoms
        DO I=1,3
          S1=0D0
          S2=0D0
          S12=0D0
          S22=0D0
          SY1=0D0
          SY2=0D0
          DO K=1,Niso
            M  = MassIso(K,J)
            W1 = M - 1D0
            W2 = 1D0/DSQRT(M) - 1D0
            S1  = S1  + W1
            S2  = S2  + W2
            S12 = S12 + W1*W2
            S22 = S22 + W2*W2
            SY1 = SY1 + W1*Rs(K,J,I)
            SY2 = SY2 + W2*Rs(K,J,I)
          END DO

C-- Solve 2x2 normal equations for A,B
          Det = S1*S22 - S2*S12
          IF (DABS(Det).GT.1D-12) THEN
             A(I) = (SY1*S22 - SY2*S12)/Det
             B(I) = (SY2*S1  - SY1*S2 )/Det
          ELSE
             A(I)=0D0
             B(I)=0D0
          END IF

C-- Compute rm coordinate for first isotopologue (reference)
          M  = MassIso(1,J)
          Rm(J,I)=Rs(1,J,I)+A(I)*(M-1D0)
     &                         +B(I)*(1D0/DSQRT(M)-1D0)
        END DO
      END DO
      RETURN
      END

C===============================================================
C  Estimate uncertainties on inertial parameters (Costain 1958)
C===============================================================
      SUBROUTINE Print_CostainErrors(DImat, dB, DeltaM, CISub)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 DImat(3,3), dB(3), DeltaM, CISub(3)
        REAL*8 H, Cspd, Conv, FACT
        REAL*8 dI(3), dA2(3), dA(3)
        REAL*8 a2(3)
        INTEGER I

C-- Physical constants
        H    = 6.62607015D-34
        Cspd = 2.99792458D8
        Conv = 1.66053906660D-47
        FACT = H / (8.D0*3.14159265358979D0**2*Cspd*Conv*1.0D6)

C-- Propagate uncertainty: dI = (FACT/B^2)*dB
        DO I=1,3
          IF (DImat(I,I).GT.0D0) THEN
             dI(I) = (FACT / (DImat(I,I)**2)) * dB(I)
          ELSE
             dI(I) = 0D0
          ENDIF
        END DO

C-- Now uncertainties on squared coordinates: d(|a|^2) = dI_corr / |Δm|
        DO I=1,3
          a2(I) = DImat(I,I)/DABS(DeltaM)
          dA2(I) = dI(I)/DABS(DeltaM)
          IF (a2(I).GT.0D0) THEN
             dA(I) = 0.5D0 * dA2(I) / DSQRT(a2(I))
          ELSE
             dA(I) = 0D0
          ENDIF
        END DO

C-- Print results
        WRITE(*,*)
        WRITE(*,*) '--- Costain-style Uncertainty Estimates ---'
        WRITE(*,'(A,3E13.5)') 'σ(Ia), σ(Ib), σ(Ic) (amu·Å²) =',
     &         dI(1), dI(2), dI(3)
        WRITE(*,'(A,3E13.5)') 'σ(|a|), σ(|b|), σ(|c|) (Å) =',
     &         dA(1), dA(2), dA(3)
        WRITE(*,*) '-------------------------------------------'
        RETURN
      END

