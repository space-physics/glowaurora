      SUBROUTINE RPOLY(OP, DEGREE, ZEROR, ZEROI,
     * FAIL)
C FINDS THE ZEROS OF A REAL POLYNOMIAL
C OP  - DOUBLE PRECISION VECTOR OF COEFFICIENTS IN
C       ORDER OF DECREASING POWERS.
C DEGREE   - INTEGER DEGREE OF POLYNOMIAL.
C ZEROR, ZEROI - OUTPUT DOUBLE PRECISION VECTORS OF
C                REAL AND IMAGINARY PARTS OF THE
C                ZEROS.
C FAIL  - OUTPUT LOGICAL PARAMETER, TRUE ONLY IF
C         LEADING COEFFICIENT IS ZERO OR IF RPOLY
C         HAS FOUND FEWER THAN DEGREE ZEROS.
C         IN THE LATTER CASE DEGREE IS RESET TO
C         THE NUMBER OF ZEROS FOUND.
C TO CHANGE THE SIZE OF POLYNOMIALS WHICH CAN BE
C SOLVED, RESET THE DIMENSIONS OF THE ARRAYS IN THE
C COMMON AREA AND IN THE FOLLOWING DECLARATIONS.
C THE SUBROUTINE USES SINGLE PRECISION CALCULATIONS
C FOR SCALING, BOUNDS AND ERROR CALCULATIONS. ALL
C CALCULATIONS FOR THE ITERATIONS ARE DONE IN DOUBLE
C PRECISION.
      Double Precision, Intent(IN) :: OP(101)
      Integer, Intent(INOUT) :: DEGREE
      Double Precision, Intent(OUT) :: ZEROR(100), ZEROI(100)
      Logical, Intent(OUT) :: FAIL

      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      Double Precision ETA, ARE, MRE, L
      INTEGER N, NN
      DOUBLE PRECISION TEMP(101),
     * T, AA, BB, CC, DABS,
     * FACTOR
      Double Precision PT(101), LO, MAX, MIN, XX, YY, COSR,
     * SINR, XXX, X, SC, BND, XM, FF, DF, DX, INFIN,
     * SMALNO, BASE
      INTEGER CNT, NZ, I, J, JJ, NM1
      LOGICAL ZEROK
C THE FOLLOWING STATEMENTS SET MACHINE CONSTANTS USED
C IN VARIOUS PARTS OF THE PROGRAM. THE MEANING OF THE
C FOUR CONSTANTS ARE...
C ETA     THE MAXIMUM RELATIVE REPRESENTATION ERROR
C         WHICH CAN BE DESCRIBED AS THE SMALLEST
C         POSITIVE FLOATING POINT NUMBER SUCH THAT
C         1.D0+ETA IS GREATER THAN 1.
C INFINY  THE LARGEST FLOATING-POINT NUMBER.
C SMALNO  THE SMALLEST POSITIVE FLOATING-POINT NUMBER
C         IF THE EXPONENT RANGE DIFFERS IN SINGLE AND
C         DOUBLE PRECISION THEN SMALNO AND INFIN
C         SHOULD INDICATE THE SMALLER RANGE.
C BASE    THE BASE OF THE FLOATING-POINT NUMBER
C         SYSTEM USED.
      BASE = 8.
      ETA = .5*BASE**(1-26)
      INFIN = 3.402823E38
      SMALNO = 1.402E-45
C ARE AND MRE REFER TO THE UNIT ERROR IN + AND *
C RESPECTIVELY. THEY ARE ASSUMED TO BE THE SAME AS
C ETA.
      ARE = ETA
      MRE = ETA
      LO = SMALNO/ETA
C INITIALIZATION OF CONSTANTS FOR SHIFT ROTATION
      XX = .70710678
      YY = -XX
      COSR = -.069756474
      SINR = .99756405
      FAIL = .FALSE.
      N = DEGREE
      NN = N + 1
C ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
      IF (OP(1).NE.0.D0) GO TO 10
      FAIL = .TRUE.
      DEGREE = 0
      RETURN
C REMOVE THE ZEROS AT THE ORIGIN IF ANY
   10 IF (OP(NN).NE.0.0D0) GO TO 20
      J = DEGREE - N + 1
      ZEROR(J) = 0.D0
      ZEROI(J) = 0.D0
      NN = NN - 1
      N = N - 1
      GO TO 10
C MAKE A COPY OF THE COEFFICIENTS
   20 DO 30 I=1,NN
        P(I) = OP(I)
   30 CONTINUE
C START THE ALGORITHM FOR ONE ZERO
   40 IF (N.GT.2) GO TO 60
      IF (N.LT.1) RETURN
C CALCULATE THE FINAL ZERO OR PAIR OF ZEROS
      IF (N.EQ.2) GO TO 50
      ZEROR(DEGREE) = -P(2)/P(1)
      ZEROI(DEGREE) = 0.0D0
      RETURN
   50 CALL QUAD(P(1), P(2), P(3), ZEROR(DEGREE-1),
     * ZEROI(DEGREE-1), ZEROR(DEGREE), ZEROI(DEGREE))
      RETURN
C FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
   60 MAX = 0.
      MIN = INFIN
      DO 70 I=1,NN
        X = ABS(SNGL(P(I)))
        IF (X.GT.MAX) MAX = X
        IF (X.NE.0. .AND. X.LT.MIN) MIN = X
   70 CONTINUE
C SCALE IF THERE ARE LARGE OR VERY SMALL COEFFICIENTS
C COMPUTES A SCALE FACTOR TO MULTIPLY THE
C COEFFICIENTS OF THE POLYNOMIAL. THE SCALING IS DONE
C TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
C INTERFERING WITH THE CONVERGENCE CRITERION.
C THE FACTOR IS A POWER OF THE BASE
      SC = LO/MIN
      IF (SC.GT.1.0) GO TO 80
      IF (MAX.LT.10.) GO TO 110
      IF (SC.EQ.0.) SC = SMALNO
      GO TO 90
   80 IF (INFIN/SC.LT.MAX) GO TO 110
   90 L = LOG(SC)/LOG(BASE) + .5
      FACTOR = (BASE*1.0D0)**L
      IF (FACTOR.EQ.1.D0) GO TO 110
      DO 100 I=1,NN
        P(I) = FACTOR*P(I)
  100 CONTINUE
C COMPUTE LOWER BOUND ON MODULI OF ZEROS.
  110 DO 120 I=1,NN
        PT(I) = ABS(SNGL(P(I)))
  120 CONTINUE
      PT(NN) = -PT(NN)
C COMPUTE UPPER ESTIMATE OF BOUND
      X = EXP((LOG(-PT(NN))-LOG(PT(1)))/FLOAT(N))
      IF (PT(N).EQ.0.) GO TO 130
C IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
      XM = -PT(NN)/PT(N)
      IF (XM.LT.X) X = XM
C CHOP THE INTERVAL (0,X) UNTIL FF .LE. 0
  130 XM = X*.1
      FF = PT(1)
      DO 140 I=2,NN
        FF = FF*XM + PT(I)
  140 CONTINUE
      IF (FF.LE.0.) GO TO 150
      X = XM
      GO TO 130
  150 DX = X
C DO NEWTON ITERATION UNTIL X CONVERGES TO TWO
C DECIMAL PLACES
  160 IF (ABS(DX/X).LE..005) GO TO 180
      FF = PT(1)
      DF = FF
      DO 170 I=2,N
        FF = FF*X + PT(I)
        DF = DF*X + FF
  170 CONTINUE
      FF = FF*X + PT(NN)
      DX = FF/DF
      X = X - DX
      GO TO 160
  180 BND = X
C COMPUTE THE DERIVATIVE AS THE INTIAL K POLYNOMIAL
C AND DO 5 STEPS WITH NO SHIFT
      NM1 = N - 1
      DO 190 I=2,N
        K(I) = FLOAT(NN-I)*P(I)/FLOAT(N)
  190 CONTINUE
      K(1) = P(1)
      AA = P(NN)
      BB = P(N)
      ZEROK = K(N).EQ.0.D0
      DO 230 JJ=1,5
        CC = K(N)
        IF (ZEROK) GO TO 210
C USE SCALED FORM OF RECURRENCE IF VALUE OF K AT 0 IS
C NONZERO
        T = -AA/CC
        DO 200 I=1,NM1
          J = NN - I
          K(J) = T*K(J-1) + P(J)
  200   CONTINUE
        K(1) = P(1)
        ZEROK = DABS(K(N)).LE.DABS(BB)*ETA*10.
        GO TO 230
C USE UNSCALED FORM OF RECURRENCE
  210   DO 220 I=1,NM1
          J = NN - I
          K(J) = K(J-1)
  220   CONTINUE
        K(1) = 0.D0
        ZEROK = K(N).EQ.0.D0
  230 CONTINUE
C SAVE K FOR RESTARTS WITH NEW SHIFTS
      DO 240 I=1,N
        TEMP(I) = K(I)
  240 CONTINUE
C LOOP TO SELECT THE QUADRATIC  CORRESPONDING TO EACH
C NEW SHIFT
      DO 280 CNT=1,20
C QUADRATIC CORRESPONDS TO A DOUBLE SHIFT TO A
C NON-REAL POINT AND ITS COMPLEX CONJUGATE. THE POINT
C HAS MODULUS BND AND AMPLITUDE ROTATED BY 94 DEGREES
C FROM THE PREVIOUS SHIFT
        XXX = COSR*XX - SINR*YY
        YY = SINR*XX + COSR*YY
        XX = XXX
        SR = BND*XX
        SI = BND*YY
        U = -2.0D0*SR
        V = BND
C SECOND STAGE CALCULATION, FIXED QUADRATIC
        CALL FXSHFR(20*CNT, NZ)
        IF (NZ.EQ.0) GO TO 260
C THE SECOND STAGE JUMPS DIRECTLY TO ONE OF THE THIRD
C STAGE ITERATIONS AND RETURNS HERE IF SUCCESSFUL.
C DEFLATE THE POLYNOMIAL, STORE THE ZERO OR ZEROS AND
C RETURN TO THE MAIN ALGORITHM.
        J = DEGREE - N + 1
        ZEROR(J) = SZR
        ZEROI(J) = SZI
        NN = NN - NZ
        N = NN - 1
        DO 250 I=1,NN
          P(I) = QP(I)
  250   CONTINUE
        IF (NZ.EQ.1) GO TO 40
        ZEROR(J+1) = LZR
        ZEROI(J+1) = LZI
        GO TO 40
C IF THE ITERATION IS UNSUCCESSFUL ANOTHER QUADRATIC
C IS CHOSEN AFTER RESTORING K
  260   DO 270 I=1,N
          K(I) = TEMP(I)
  270   CONTINUE
  280 CONTINUE
C RETURN WITH FAILURE IF NO CONVERGENCE WITH 20
C SHIFTS
      FAIL = .TRUE.
      DEGREE = DEGREE - N
      END SUBROUTINE RPOLY

      SUBROUTINE FXSHFR(L2, NZ)
C COMPUTES UP TO  L2  FIXED SHIFT K-POLYNOMIALS,
C TESTING FOR CONVERGENCE IN THE LINEAR OR QUADRATIC
C CASE. INITIATES ONE OF THE VARIABLE SHIFT
C ITERATIONS AND RETURNS WITH THE NUMBER OF ZEROS
C FOUND.
C L2 - LIMIT OF FIXED SHIFT STEPS
C NZ - NUMBER OF ZEROS FOUND
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      DOUBLE PRECISION ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION SVU, SVV, UI, VI, S
      DOUBLE PRECISION BETAS, BETAV, OSS, OVV, SS, VV, TS, TV,
     * OTS, OTV, TVV, TSS
      INTEGER L2, NZ, TYPE, I, J, IFLAG
      LOGICAL VPASS, SPASS, VTRY, STRY
      NZ = 0
      BETAV = .25
      BETAS = .25
      OSS = SR
      OVV = V
C EVALUATE POLYNOMIAL BY SYNTHETIC DIVISION
      CALL QUADSD(NN, U, V, P, QP, A, B)
      CALL CALCSC(TYPE)
      DO 80 J=1,L2
C CALCULATE NEXT K POLYNOMIAL AND ESTIMATE V
        CALL NEXTK(TYPE)
        CALL CALCSC(TYPE)
        CALL NEWEST(TYPE, UI, VI)
        VV = VI
C ESTIMATE S
        SS = 0.
        IF (K(N).NE.0.D0) SS = -P(NN)/K(N)
        TV = 1.
        TS = 1.
        IF (J.EQ.1 .OR. TYPE.EQ.3) GO TO 70
C COMPUTE RELATIVE MEASURES OF CONVERGENCE OF S AND V
C SEQUENCES
        IF (VV.NE.0.) TV = ABS((VV-OVV)/VV)
        IF (SS.NE.0.) TS = ABS((SS-OSS)/SS)
C IF DECREASING, MULTIPLY TWO MOST RECENT
C CONVERGENCE MEASURES
        TVV = 1.
        IF (TV.LT.OTV) TVV = TV*OTV
        TSS = 1.
        IF (TS.LT.OTS) TSS = TS*OTS
C COMPARE WITH CONVERGENCE CRITERIA
        VPASS = TVV.LT.BETAV
        SPASS = TSS.LT.BETAS
        IF (.NOT.(SPASS .OR. VPASS)) GO TO 70
C AT LEAST ONE SEQUENCE HAS PASSED THE CONVERGENCE
C TEST. STORE VARIABLES BEFORE ITERATING
        SVU = U
        SVV = V
        DO 10 I=1,N
          SVK(I) = K(I)
   10   CONTINUE
        S = SS
C CHOOSE ITERATION ACCORDING TO THE FASTEST
C CONVERGING SEQUENCE
        VTRY = .FALSE.
        STRY = .FALSE.
        IF (SPASS .AND. ((.NOT.VPASS) .OR.
     *   TSS.LT.TVV)) GO TO 40
   20   CALL QUADIT(UI, VI, NZ)
        IF (NZ.GT.0) RETURN
C QUADRATIC ITERATION HAS FAILED. FLAG THAT IT HAS
C BEEN TRIED AND DECREASE THE CONVERGENCE CRITERION.
        VTRY = .TRUE.
        BETAV = BETAV*.25
C TRY LINEAR ITERATION IF IT HAS NOT BEEN TRIED AND
C THE S SEQUENCE IS CONVERGING
        IF (STRY .OR. (.NOT.SPASS)) GO TO 50
        DO 30 I=1,N
          K(I) = SVK(I)
   30   CONTINUE
   40   CALL REALIT(S, NZ, IFLAG)
        IF (NZ.GT.0) RETURN
C LINEAR ITERATION HAS FAILED. FLAG THAT IT HAS BEEN
C TRIED AND DECREASE THE CONVERGENCE CRITERION
        STRY = .TRUE.
        BETAS = BETAS*.25
        IF (IFLAG.EQ.0) GO TO 50
C IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL
C ZERO ATTEMPT QUADRATIC INTERATION
        UI = -(S+S)
        VI = S*S
        GO TO 20
C RESTORE VARIABLES
   50   U = SVU
        V = SVV
        DO 60 I=1,N
          K(I) = SVK(I)
   60   CONTINUE
C TRY QUADRATIC ITERATION IF IT HAS NOT BEEN TRIED
C AND THE V SEQUENCE IS CONVERGING
        IF (VPASS .AND. (.NOT.VTRY)) GO TO 20
C RECOMPUTE QP AND SCALAR VALUES TO CONTINUE THE
C SECOND STAGE
        CALL QUADSD(NN, U, V, P, QP, A, B)
        CALL CALCSC(TYPE)
   70   OVV = VV
        OSS = SS
        OTV = TV
        OTS = TS
   80 CONTINUE
      END SUBROUTINE FXSHFR

      SUBROUTINE QUADIT(UU, VV, NZ)
C VARIABLE-SHIFT K-POLYNOMIAL ITERATION FOR A
C QUADRATIC FACTOR CONVERGES ONLY IF THE ZEROS ARE
C EQUIMODULAR OR NEARLY SO.
C UU,VV - COEFFICIENTS OF STARTING QUADRATIC
C NZ - NUMBER OF ZERO FOUND
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      DOUBLE PRECISION ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION UI, VI, UU, VV, DABS
      DOUBLE PRECISION MP, OMP, EE, RELSTP, T, ZM
      INTEGER NZ, TYPE, I, J
      LOGICAL TRIED
      NZ = 0
      TRIED = .FALSE.
      U = UU
      V = VV
      J = 0
C MAIN LOOP
   10 CALL QUAD(1.D0, U, V, SZR, SZI, LZR, LZI)
C RETURN IF ROOTS OF THE QUADRATIC ARE REAL AND NOT
C CLOSE TO MULTIPLE OR NEARLY EQUAL AND  OF OPPOSITE
C SIGN
      IF (DABS(DABS(SZR)-DABS(LZR)).GT..01D0*
     * DABS(LZR)) RETURN
C EVALUATE POLYNOMIAL BY QUADRATIC SYNTHETIC DIVISION
      CALL QUADSD(NN, U, V, P, QP, A, B)
      MP = DABS(A-SZR*B) + DABS(SZI*B)
C COMPUTE A RIGOROUS  BOUND ON THE ROUNDING ERROR IN
C EVALUTING P
      ZM = SQRT(ABS(SNGL(V)))
      EE = 2.*ABS(SNGL(QP(1)))
      T = -SZR*B
      DO 20 I=2,N
        EE = EE*ZM + ABS(SNGL(QP(I)))
   20 CONTINUE
      EE = EE*ZM + ABS(SNGL(A)+T)
      EE = (5.*MRE+4.*ARE)*EE - (5.*MRE+2.*ARE)*
     * (ABS(SNGL(A)+T)+ABS(SNGL(B))*ZM) +
     * 2.*ARE*ABS(T)
C ITERATION HAS CONVERGED SUFFICIENTLY IF THE
C POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND
      IF (MP.GT.20.*EE) GO TO 30
      NZ = 2
      RETURN
   30 J = J + 1
C STOP ITERATION AFTER 20 STEPS
      IF (J.GT.20) RETURN
      IF (J.LT.2) GO TO 50
      IF (RELSTP.GT..01 .OR. MP.LT.OMP .OR. TRIED)
     * GO TO 50
C A CLUSTER APPEARS TO BE STALLING THE CONVERGENCE.
C FIVE FIXED SHIFT STEPS ARE TAKEN WITH A U,V CLOSE
C TO THE CLUSTER
      IF (RELSTP.LT.ETA) RELSTP = ETA
      RELSTP = SQRT(RELSTP)
      U = U - U*RELSTP
      V = V + V*RELSTP
      CALL QUADSD(NN, U, V, P, QP, A, B)
      DO 40 I=1,5
        CALL CALCSC(TYPE)
        CALL NEXTK(TYPE)
   40 CONTINUE
      TRIED = .TRUE.
      J = 0
   50 OMP = MP
C CALCULATE NEXT K POLYNOMIAL AND NEW U AND V
      CALL CALCSC(TYPE)
      CALL NEXTK(TYPE)
      CALL CALCSC(TYPE)
      CALL NEWEST(TYPE, UI, VI)
C IF VI IS ZERO THE ITERATION IS NOT CONVERGING
      IF (VI.EQ.0.D0) RETURN
      RELSTP = DABS((VI-V)/VI)
      U = UI
      V = VI
      GO TO 10
      END SUBROUTINE QUADIT

      SUBROUTINE REALIT(SSS, NZ, IFLAG)
C VARIABLE-SHIFT H POLYNOMIAL ITERATION FOR A REAL
C ZERO.
C SSS   - STARTING ITERATE
C NZ    - NUMBER OF ZERO FOUND
C IFLAG - FLAG TO INDICATE A PAIR OF ZEROS NEAR REAL
C         AXIS.
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      DOUBLE PRECISION ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION PV, KV, T, S, SSS
      DOUBLE PRECISION MS, MP, OMP, EE
      INTEGER NZ, IFLAG, I, J, NM1
      NM1 = N - 1
      NZ = 0
      S = SSS
      IFLAG = 0
      J = 0
C MAIN LOOP
   10 PV = P(1)
C EVALUATE P AT S
      QP(1) = PV
      DO 20 I=2,NN
        PV = PV*S + P(I)
        QP(I) = PV
   20 CONTINUE
      MP = DABS(PV)
C COMPUTE A RIGOROUS BOUND ON THE ERROR IN EVALUATING
C P
      MS = DABS(S)
      EE = (MRE/(ARE+MRE))*ABS(SNGL(QP(1)))
      DO 30 I=2,NN
        EE = EE*MS + ABS(SNGL(QP(I)))
   30 CONTINUE
C ITERATION HAS CONVERGED SUFFICIENTLY IF THE
C POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND
      IF (MP.GT.20.*((ARE+MRE)*EE-MRE*MP)) GO TO 40
      NZ = 1
      SZR = S
      SZI = 0.D0
      RETURN
   40 J = J + 1
C STOP ITERATION AFTER 10 STEPS
      IF (J.GT.10) RETURN
      IF (J.LT.2) GO TO 50
      IF (DABS(T).GT..001*DABS(S-T) .OR. MP.LE.OMP)
     * GO TO 50
C A CLUSTER OF ZEROS NEAR THE REAL AXIS HAS BEEN
C ENCOUNTERED RETURN WITH IFLAG SET TO INITIATE A
C QUADRATIC ITERATION
      IFLAG = 1
      SSS = S
      RETURN
C RETURN IF THE POLYNOMIAL VALUE HAS INCREASED
C SIGNIFICANTLY
   50 OMP = MP
C COMPUTE T, THE NEXT POLYNOMIAL, AND THE NEW ITERATE
      KV = K(1)
      QK(1) = KV
      DO 60 I=2,N
        KV = KV*S + K(I)
        QK(I) = KV
   60 CONTINUE
      IF (DABS(KV).LE.DABS(K(N))*10.*ETA) GO TO 80
C USE THE SCALED FORM OF THE RECURRENCE IF THE VALUE
C OF K AT S IS NONZERO
      T = -PV/KV
      K(1) = QP(1)
      DO 70 I=2,N
        K(I) = T*QK(I-1) + QP(I)
   70 CONTINUE
      GO TO 100
C USE UNSCALED FORM
   80 K(1) = 0.0D0
      DO 90 I=2,N
        K(I) = QK(I-1)
   90 CONTINUE
  100 KV = K(1)
      DO 110 I=2,N
        KV = KV*S + K(I)
  110 CONTINUE
      T = 0.D0
      IF (DABS(KV).GT.DABS(K(N))*10.*ETA) T = -PV/KV
      S = S + T
      GO TO 10
      END SUBROUTINE REALIT

      SUBROUTINE CALCSC(TYPE)
C THIS ROUTINE CALCULATES SCALAR QUANTITIES USED TO
C COMPUTE THE NEXT K POLYNOMIAL AND NEW ESTIMATES OF
C THE QUADRATIC COEFFICIENTS.
C TYPE - INTEGER VARIABLE SET HERE INDICATING HOW THE
C CALCULATIONS ARE NORMALIZED TO AVOID OVERFLOW
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      DOUBLE PRECISION ETA, ARE, MRE
      INTEGER N, NN
      INTEGER TYPE
C SYNTHETIC DIVISION OF K BY THE QUADRATIC 1,U,V
      CALL QUADSD(N, U, V, K, QK, C, D)
      IF (DABS(C).GT.DABS(K(N))*100.*ETA) GO TO 10
      IF (DABS(D).GT.DABS(K(N-1))*100.*ETA) GO TO 10
      TYPE = 3
C TYPE=3 INDICATES THE QUADRATIC IS ALMOST A FACTOR
C OF K
      RETURN
   10 IF (DABS(D).LT.DABS(C)) GO TO 20
      TYPE = 2
C TYPE=2 INDICATES THAT ALL FORMULAS ARE DIVIDED BY D
      E = A/D
      F = C/D
      G = U*B
      H = V*B
      A3 = (A+G)*E + H*(B/D)
      A1 = B*F - A
      A7 = (F+U)*A + H
      RETURN
   20 TYPE = 1
C TYPE=1 INDICATES THAT ALL FORMULAS ARE DIVIDED BY C
      E = A/C
      F = D/C
      G = U*E
      H = V*B
      A3 = A*E + (H/C+G)*B
      A1 = B - A*(D/C)
      A7 = A + G*D + H*F
      END SUBROUTINE CALCSC

      SUBROUTINE NEXTK(TYPE)
C COMPUTES THE NEXT K POLYNOMIALS USING SCALARS
C COMPUTED IN CALCSC
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      DOUBLE PRECISION ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION TEMP
      INTEGER TYPE
      IF (TYPE.EQ.3) GO TO 40
      TEMP = A
      IF (TYPE.EQ.1) TEMP = B
      IF (DABS(A1).GT.DABS(TEMP)*ETA*10.) GO TO 20
C IF A1 IS NEARLY ZERO THEN USE A SPECIAL FORM OF THE
C RECURRENCE
      K(1) = 0.D0
      K(2) = -A7*QP(1)
      DO 10 I=3,N
        K(I) = A3*QK(I-2) - A7*QP(I-1)
   10 CONTINUE
      RETURN
C USE SCALED FORM OF THE RECURRENCE
   20 A7 = A7/A1
      A3 = A3/A1
      K(1) = QP(1)
      K(2) = QP(2) - A7*QP(1)
      DO 30 I=3,N
        K(I) = A3*QK(I-2) - A7*QP(I-1) + QP(I)
   30 CONTINUE
      RETURN
C USE UNSCALED FORM OF THE RECURRENCE IF TYPE IS 3
   40 K(1) = 0.D0
      K(2) = 0.D0
      DO 50 I=3,N
        K(I) = QK(I-2)
   50 CONTINUE
      END SUBROUTINE NEXTK

      SUBROUTINE NEWEST(TYPE, UU, VV)
C COMPUTE NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS
C USING THE SCALARS COMPUTED IN CALCSC.
      COMMON /GLOBAL/ P, QP, K, QK, SVK, SR, SI, U,
     * V, A, B, C, D, A1, A2, A3, A6, A7, E, F, G,
     * H, SZR, SZI, LZR, LZI, ETA, ARE, MRE, N, NN
      DOUBLE PRECISION P(101), QP(101), K(101),
     * QK(101), SVK(101), SR, SI, U, V, A, B, C, D,
     * A1, A2, A3, A6, A7, E, F, G, H, SZR, SZI,
     * LZR, LZI
      DOUBLE PRECISION ETA, ARE, MRE
      INTEGER N, NN
      DOUBLE PRECISION A4, A5, B1, B2, C1, C2, C3,
     * C4, TEMP, UU, VV
      INTEGER TYPE
C USE FORMULAS APPROPRIATE TO SETTING OF TYPE.
      IF (TYPE.EQ.3) GO TO 30
      IF (TYPE.EQ.2) GO TO 10
      A4 = A + U*B + H*F
      A5 = C + (U+V*F)*D
      GO TO 20
   10 A4 = (A+G)*F + H
      A5 = (F+U)*C + V*D
C EVALUATE NEW QUADRATIC COEFFICIENTS.
   20 B1 = -K(N)/P(NN)
      B2 = -(K(N-1)+B1*P(N))/P(NN)
      C1 = V*B2*A1
      C2 = B1*A7
      C3 = B1*B1*A3
      C4 = C1 - C2 - C3
      TEMP = A5 + B1*A4 - C4
      IF (TEMP.EQ.0.D0) GO TO 30
      UU = U - (U*(C3+C2)+V*(B1*A1+B2*A7))/TEMP
      VV = V*(1.+C4/TEMP)
      RETURN
C IF TYPE=3 THE QUADRATIC IS ZEROED
   30 UU = 0.D0
      VV = 0.D0
      END SUBROUTINE NEWEST

      SUBROUTINE QUADSD(NN, U, V, P, Q, A, B)
C DIVIDES P BY THE QUADRATIC  1,U,V  PLACING THE
C QUOTIENT IN Q AND THE REMAINDER IN A,B
      DOUBLE PRECISION P(NN), Q(NN), U, V, A, B, C
      INTEGER I
      B = P(1)
      Q(1) = B
      A = P(2) - U*B
      Q(2) = A
      DO 10 I=3,NN
        C = P(I) - U*A - V*B
        Q(I) = C
        B = A
        A = C
   10 CONTINUE
      END SUBROUTINE QUADSD

      SUBROUTINE QUAD(A, B1, C, SR, SI, LR, LI)
C CALCULATE THE ZEROS OF THE QUADRATIC A*Z**2+B1*Z+C.
C THE QUADRATIC FORMULA, MODIFIED TO AVOID
C OVERFLOW, IS USED TO FIND THE LARGER ZERO IF THE
C ZEROS ARE REAL AND BOTH ZEROS ARE COMPLEX.
C THE SMALLER REAL ZERO IS FOUND DIRECTLY FROM THE
C PRODUCT OF THE ZEROS C/A.
      DOUBLE PRECISION A, B1, C, SR, SI, LR, LI, B,
     * D, E
      IF (A.NE.0.D0) GO TO 20
      SR = 0.D0
      IF (B1.NE.0.D0) SR = -C/B1
      LR = 0.D0
   10 SI = 0.D0
      LI = 0.D0
      RETURN
   20 IF (C.NE.0.D0) GO TO 30
      SR = 0.D0
      LR = -B1/A
      GO TO 10
C COMPUTE DISCRIMINANT AVOIDING OVERFLOW
   30 B = B1/2.D0
      IF (DABS(B).LT.DABS(C)) GO TO 40
      E = 1.D0 - (A/B)*(C/B)
      D = DSQRT(DABS(E))*DABS(B)
      GO TO 50
   40 E = A
      IF (C.LT.0.D0) E = -A
      E = B*(B/DABS(C)) - E
      D = DSQRT(DABS(E))*DSQRT(DABS(C))
   50 IF (E.LT.0.D0) GO TO 60
C REAL ZEROS
      IF (B.GE.0.D0) D = -D
      LR = (-B+D)/A
      SR = 0.D0
      IF (LR.NE.0.D0) SR = (C/LR)/A
      GO TO 10
C COMPLEX CONJUGATE ZEROS
   60 SR = -B/A
      LR = SR
      SI = DABS(D/A)
      LI = -SI
      END SUBROUTINE QUAD
