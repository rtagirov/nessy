      module MOD_PRITAU

      contains

      SUBROUTINE PRITAU(MODHEAD,JOBNUM,RSTAR,ND,RADIUS,RNE,ENTOT,T,
     $                  POPNUM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,
     $                  EION,EINST,ALPHA,SEXPO,AGAUNT,NOM,XLAMBDA,
     $                  FWEIGHT,TAUROSS,WAVARR,SIGARR, NF, NFDIM)
C***********************************************************************
C***  PRINTOUT OF THE NLTE OPTICAL DEPTH SCALES (ROSSELAND, THOMSON)
C***********************************************************************
      USE MOD_LIPO
      use MOD_OPAROSS
      IMPLICIT REAL*8(A-H,O-Z)
     
      DIMENSION RADIUS(ND),RNE(ND),ENTOT(ND),T(ND)
      DIMENSION TAUTHOM(ND),TAUROSS(ND)
      DIMENSION POPNUM(ND,N)
      DIMENSION EN(N)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N),ALPHA(N),SEXPO(N)
      DIMENSION NOM(N)
      DIMENSION EINST(N,N)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)

	  DIMENSION WAVARR(N,NFDIM),SIGARR(N,NFDIM)

      CHARACTER*10 LEVEL(N)
      CHARACTER MODHEAD*104
      character*8,dimension(N) :: AGAUNT
     
C***  SIGMATH = THOMSON CROSS-SECTION FOR ELECTRON SCATTERING  (IN CM**2)
      DATA SIGMATH / 6.652E-25 /
     
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,
     $       ///,20X,'O P T I C A L   D E P T H   S C A L E S',
     $       /,20X,39('='),//,
     $       1X,'DEPTH',7X,'R-1',6X,'LOG(R-1)',3X,'EL. TEMPERATURE',2X,
     $       'LOG(PARTICLE DENSITY)',2X,'LOG(EL. DENSITY)',2X,
     $       'TAU-THOMSON',2X,'TAU-ROSSELAND',/,
     $       1X,'INDEX',31X,'(KELVIN)',8X,'(ATOMS PER CM+3)',5X,
     $       '(EL. PER CM+3)',6X,'(NLTE)',8X,'(NLTE)',/)
     
      RSIGTH=RSTAR*SIGMATH
      TAUTHOM(1)=0.0
      TAUROSS(1)=0.0
C***  LOOP OVER ALL DEPTH POINTS  --------------------------------------
      DO 2 L=1,ND
      RL=RADIUS(L)
      RL1=RL-1.
      IF (L .LT. ND) THEN
         RLOG=LOG10(RL1)
      ELSE
         RLOG=-9.99d+99
      ENDIF
      TL=T(L)
      ENTOTL=ENTOT(L)
      DENS=LOG10(ENTOTL)
      RNEL=RNE(L)
      ENEL=RNEL*ENTOTL
      ENELOG=LOG10(ENEL)
      DO 20 I=1,N
   20 EN(I)=POPNUM(L,I)
C***  CHANGES BY MARGIT HABERREITER
      CALL OPAROSS(OPARL,EN,TL,RNEL,ENTOTL,RSTAR,N,
     $             LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $             ALPHA,SEXPO,AGAUNT,NF,NFDIM,XLAMBDA,FWEIGHT,NOM,
     $             WAVARR,SIGARR)

      IF (L .GT. 1) THEN
         DR=RM1-RL
         ENEMEAN=0.5d0*(ENEL+ENELM1)
         TAUTHOM(L)=TAUTHOM(L-1)+RSIGTH*ENEMEAN*DR
         OPAMEAN=0.5d0*(OPARL+OPARM1)
         TAUROSS(L)=TAUROSS(L-1)+OPAMEAN*DR
      ENDIF
      PRINT 3, L,RL1,RLOG,TL,DENS,ENELOG,TAUTHOM(L),TAUROSS(L)
    3 FORMAT (2X,I3,4X,G10.3,2X,G10.3,5X,F9.0,12X,F7.3,13X,F8.3,8X,F7.3,
     $        7X,F7.3)
      RM1=RL
      ENELM1=ENEL
      OPARM1=OPARL
    2 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
     
      PRINT 17
   17 FORMAT (5/,1X,'RADII AND TEMPERATURES FOR DIFFERENT OPTICAL',
     $       ' DEPTHS:      RADIUS',5X,'TEMPERATURE',/,58X,24('-'))
C***  CALCULATE RADII AND CORRESPONDING TEMPERATURES FOR TAU-THOMSON = 1/3,
C***                                                                 = 2/3,
C***                                                                 = 1.0
      TAU13=0.333333333333d0
      IF (TAUTHOM(ND) .LT. TAU13) THEN
         R13=1.d0
         T13=T(ND)
      ELSE
         CALL LIPO (R13,TAU13,RADIUS,TAUTHOM,ND)
         CALL LIPO (T13,R13,T,RADIUS,ND)
      ENDIF
      TAU23=0.666666666666d0
      IF (TAUTHOM(ND) .LT. TAU23) THEN
         R23=1.d0
         T23=T(ND)
      ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUTHOM,ND)
         CALL LIPO (T23,R23,T,RADIUS,ND)
      ENDIF
      TAU1=1.d0
      IF (TAUTHOM(ND) .LT. TAU1 ) THEN
         R1 =1.d0
         T1 =T(ND)
      ELSE
         CALL LIPO (R1 ,TAU1 ,RADIUS,TAUTHOM,ND)
         CALL LIPO (T1 ,R1 ,T,RADIUS,ND)
      ENDIF
      PRINT 18, R13,T13,R23,T23,R1,T1
   18 FORMAT (35X,'TAU-THOMSON = 1/3:',5X,F7.3,6X,F9.0,/,47X,'= 2/3:',
     $        5X,F7.3,6X,F9.0,/,47X,'= 1.0:',5X,F7.3,6X,F9.0,/)
C***  CALCULATE RADII AND CORRESPONDING TEMPERATURES FOR TAU-ROSSELAND = 1/3,
C***                                                                   = 2/3,
C***                                                                   = 1.0
      IF (TAUROSS(ND) .LT. TAU13) THEN
         R13=1.d0
         T13=T(ND)
      ELSE
         CALL LIPO (R13,TAU13,RADIUS,TAUROSS,ND)
         CALL LIPO (T13,R13,T,RADIUS,ND)
      ENDIF
      IF (TAUROSS(ND) .LT. TAU23) THEN
         R23=1.d0
         T23=T(ND)
      ELSE
         CALL LIPO (R23,TAU23,RADIUS,TAUROSS,ND)
         CALL LIPO (T23,R23,T,RADIUS,ND)
      ENDIF
      IF (TAUROSS(ND) .LT. TAU1 ) THEN
         R1 =1.d0
         T1 =T(ND)
      ELSE
         CALL LIPO (R1 ,TAU1 ,RADIUS,TAUROSS,ND)
         CALL LIPO (T1 ,R1 ,T,RADIUS,ND)
      ENDIF
      PRINT 19, R13,T13,R23,T23,R1,T1
   19 FORMAT (33X,'TAU-ROSSELAND = 1/3:',5X,F7.3,6X,F9.0,/,47X,'= 2/3:',
     $        5X,F7.3,6X,F9.0,/,47X,'= 1.0:',5X,F7.3,6X,F9.0)
     
      RETURN
      END subroutine
      end module
