      module MOD_OMEG
      contains
      SUBROUTINE OMEG(N,TL,NCHARG,ELEVEL,EINST, OMEGA, NUP, LOW,
     $                EION,COCO,KEYCOL,WEIGHT,ALTESUM,NATOM,NOM,KODAT)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION EINST(N,N)
      DIMENSION NCHARG(N),ELEVEL(N),WEIGHT(N)
      DIMENSION EION(N),ALTESUM(4,N)
      DIMENSION COCO(N,N,4)
      DIMENSION NOM(N)
      DIMENSION KODAT(30)
      CHARACTER*4 KEYCOL(N,N)
     
C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
      
	if (natom.gt.30) stop 'omeg/natom.gt.30'
      TROOT=SQRT(TL)
      T32=TL*TROOT
     
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      WN2=WAVENUM*WAVENUM
      WN3=WN2*WAVENUM

C***  ONLY FOR HELIUM ...
      IF (NOM(LOW) .EQ. KODAT(1)) THEN
C***  HE I  *******************************************************************
C***  OMEGA(UP-LOW) IS CALCULATED DEPENDING ON KEYWORD KEYCOL  *********
       IF (NCHARG(LOW).EQ.0) THEN
C***  'JEFF': OPTICALLY PERMITTED  TRANSITIONS. JEFFERIES P. 118 (EQ. 6.24)
C             (VAN REGENMORTER FORMULA) 3.24 = 2.16 / 0.667
         IF (KEYCOL(NUP,LOW) .EQ. 'JEFF') THEN
C                                  ====
              OMEGA=3.24*EINST(NUP,LOW)/WN2/T32/(C1*WAVENUM/TL)**1.68
     
C***  'BFK ': OPTICALLY FORBIDDEN TRANSITIONS BETWEEN N=2, N=1 OR WITHIN N=2:
C***          BERRINGTON, FON + KINGSTON ('BFK') 1982, MNRAS 200, 347
         ELSE IF (KEYCOL(NUP,LOW) .EQ. 'BFK1') THEN
C                                       ====
              PBFK=COCO(NUP,LOW,1)/TROOT
     +                          +COCO(NUP,LOW,2)+COCO(NUP,LOW,3)*TROOT
     +                          +COCO(NUP,LOW,4)*T32
              OMEGA=PBFK*WEIGHT(LOW)/WEIGHT(NUP)
     
         ELSE IF (KEYCOL(NUP,LOW) .EQ. 'BFK2') THEN
C                                       ====
              PBFK=COCO(NUP,LOW,1)/TL
     +                          +COCO(NUP,LOW,2)/TROOT+COCO(NUP,LOW,3)
     +                          +COCO(NUP,LOW,4)*TL
              OMEGA=PBFK*WEIGHT(LOW)/WEIGHT(NUP)
     
C***  'BKMS' OR 'BKGR':
C***  FORBIDDEN TRANSITIONS BETWEEN N.GT.2, N=1 OR N.GT.3, N.NE.1:
C***  BENSON + KULANDER ('BK..') 1972, SOLAR PHYSICS 27, 305 (FORMULA 3)
C***  '..MS': MIHALAS + STONE (REF. 11)
C***  '..GR': GREEN (REF. 7)
C***  ATTENTION: COCO(.,.,3) := 1.-ALPHA !!!!!!!!!!!!!!!!!
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'BKMS'
C                                        ====
     $        .OR. KEYCOL(NUP,LOW) .EQ. 'BKGR') THEN
C                                        ====
              OMEGA=COCO(NUP,LOW,1)*TL**COCO(NUP,LOW,2)*
     *               EXP(COCO(NUP,LOW,3)*C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)
     
C***  'UPS.': OPTICALLY FORBIDDEN TRANSITIONS BETWEEN N.GT.2, N.NE.1:
C***        EFFECTIVE COLLISIONAL STRENGTH (UPSILON) FROM WERNER SCHMUTZ
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'UPS0') THEN
C                                        ====
            OMEGA=.0
     
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'UPS1') THEN
C                                        ====
              UPSILON=0.05
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
     
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'UPS2') THEN
C                                        ====
              UPSILON=1.
              OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT

C***  'NONE': TRANSITIONS (FROM N.GT.4) WITH UNKNOWN COLLISIONAL COEFFICIENTS
C***          (COLLISIONAL CROSS SECTION SIGMA(LOW,UP) IS SET TO  PI*A0**2)
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'NONE') THEN
C                                        ====
              OMEGA=5.465E-11*TROOT*(1.+C1*WAVENUM/TL)*WEIGHT(LOW)/
     /               WEIGHT(NUP)
     
C***  'ZERO': NO COLLISIONAL TRANSITION
          ELSE IF (KEYCOL(NUP,LOW) .EQ. 'ZERO') THEN
C                                        ====
              OMEGA=0.
     
          ELSE
C***  NO/UNKNOWN KEYWORD 'KEYCOL' DECODED
              PRINT *,KEYCOL(NUP,LOW),NUP,LOW
              print *,' OMEG: HE I: WRONG KEYWORD KEYCOL '
              STOP 'ERROR'
          ENDIF
     
C***  HE II  ******************************************************************
       ELSE IF (NCHARG(LOW).EQ.1) THEN
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 119 (EQ. 6.25)
          OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT
C***  HE III  *****************************************************************
       ELSE
          print *,' OMEG: HE III - LINE TRANSITION  '
          STOP 'ERROR'
       ENDIF
      ELSE IF (NOM(LOW) .EQ. KODAT(2)) THEN
C***  H I  *************************************************************
         IF (NCHARG(LOW) .EQ. 0) THEN
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 118 (EQ. 6.24)
            OMEGA=3.24*EINST(NUP,LOW)/WN2/T32/(C1*WAVENUM/TL)**1.68

C***  H II  ************************************************************
         ELSE
            print *,' OMEG: H II - LINE TRANSITION'
            STOP 'ERROR'
         ENDIF
C***  CARBON AND ALL METALS ************************************************************
C***  METALS: INDEX GE 3
      ELSE IF (NOM(LOW) .GE. 3) THEN
C***  ALL IONIZATION STAGES
C***  UPSILON TEMPERATURE INDEPENDENT
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 118 (EQ. 6.24)
C                                        
         IF (KEYCOL(NUP,LOW) .EQ. 'CARB') THEN
C                                  ====
            UPSILON=COCO(NUP,LOW,1)
            IF (UPSILON.LE.0.) THEN
               print *,' OMEG: C - UPSILON .LE. 0'
               STOP 'ERROR'
            ENDIF
            OMEGA=8.6287E-6*UPSILON/WEIGHT(NUP)/TROOT
         ELSE
C***  OMEGA (UP-LOW) ACCORDING TO JEFFERIES P. 119 (EQ. 6.25)
            OMEGA=4.06*EINST(NUP,LOW)/WN3/TROOT
c            CALL REMARK ('C - LINE TRANSITION')
c            STOP 'ERROR'
         ENDIF
C***=========================================================
      ELSE
         print *,' OMEG: LINE TRANSITION OF UNKNOWN ELEMENT'
         STOP 'ERROR'
      ENDIF
     
      RETURN
      END subroutine
      end module
