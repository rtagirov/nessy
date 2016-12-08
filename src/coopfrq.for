      module MOD_COOPFRQ

      contains

!***  NON-LTE CONTINUOUS OPACITY AT CURRENT DEPTH FOR ALL FREQUENCIES
!***  NOTE: ONLY TRUE OPACITY, WITHOUT THOMSON SCATTERING TERM.

      SUBROUTINE COOPFRQ (NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $                   WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,
     $                   RNEL,TL)
      use MOD_GAUNTFF
      implicit real*8(a-h,o-z)

      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N),NFEDGE(N),EN(N)
      DIMENSION NOM(N)
      DIMENSION OPAC(NF),ETAC(NF),XLAMBDA(NF),EXPFAC(NF)
      DIMENSION SIGMAKI(NF,N)
      DIMENSION GFF(-1:6),GIIIX(-1:6)
      parameter C1 = 1.4388d0  ! = H * C / K    ( CM * ANGSTROEM )
      parameter C2 = 3.972d-16 ! = 2 * H * C    ( G * CM**3 / S**2 )
      parameter C3 = 2.07d-16  ! = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
      parameter CFF= 1.370d-23 ! = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN PAGE 100 )
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL
      TLOG=LOG10(TL)
      OPAC(1:NF)=.0
      ETAC(1:NF)=.0

      !***  BOUND-FREE  ************************************************
      DO NUP=2,N
        !***  UPPER LEVEL MUST BE GROUND STATE OF HIGHER IONIZATION STAGES
        IF ((NOM(NUP) .NE. NOM(NUP-1)) .OR.
     $      (NCHARG(NUP) .NE. NCHARG(NUP-1)+1)) cycle ! i.e. jump to ENDDO
        DO LOW=1,NUP-1
          !***  LOWER LEVEL MUST BELONG TO THE SAME ELEMENT AS UPPER LEVEL
          IF (NOM(LOW) .NE. NOM(NUP)) cycle

          !***  CHARGES MUST DIFFER BY 1
          IF (NCHARG(NUP).NE.NCHARG(LOW)+1 ) cycle

          EDGE=EION(LOW)-ELEVEL(LOW)
          EXPEDGE=EXP(C1*EDGE/TL)
          WE=C3*RNEL*ENTOTL/T32 *WEIGHT(LOW)/WEIGHT(NUP)*EXPEDGE
          NFLOW=NFEDGE(LOW)
          DO K=1,NFLOW
            SIGMA=SIGMAKI(K,LOW)
            !*** RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
            G=WE*EXPFAC(K)
            EMINDU=G*EN(NUP)*SIGMA
            SUM=    EN(LOW)*SIGMA-EMINDU
            OPAC(K)=OPAC(K)+SUM
            ETAC(K)=ETAC(K)+EMINDU
          ENDDO
        ENDDO
      ENDDO
      !***  FREE-FREE  *******************************************************
      !***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
      DO K=1,NF
        W=1.d8/XLAMBDA(K)
        W3=W*W*W
        XLAMLOG=LOG10(XLAMBDA(K))
        DO IG = 0,6
          IF (IG.GT.0) THEN
            CALL GFFLOG (GIIIX(IG),IG,XLAMLOG,TLOG)
          ELSE
            GIIIX(IG)=0.
          ENDIF
            GFF(IG)=GIIIX(IG)*IG*IG
        ENDDO
        !***  PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
        PRESIG =CFF/W3/ROOTTL
        DO I=1,N
          SIGMAFF=PRESIG*GFF(NCHARG(I))
          SUM=RNEL*ENTOTL*EN(I)*SIGMAFF
          EMINDU=SUM*EXPFAC(K)
          SUM=SUM-EMINDU
          OPAC(K)=OPAC(K)+SUM
          ETAC(K)=ETAC(K)+EMINDU
        ENDDO
        ETAC(K)=ETAC(K)*C2*W3*ENTOTL*RSTAR
        OPAC(K)=OPAC(K)*ENTOTL*RSTAR
      ENDDO

      RETURN
      END SUBROUTINE
      end module
