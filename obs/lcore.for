      MODULE MOD_LCORE

      CONTAINS

      SUBROUTINE LCORE(XRED,XBLUE,GAMMAL,LASTIND,INDLOW,INDNUP,
     $                 GAMMAR,IPRILC,MODHEAD,JOBNUM,
     $                 L,ND,VELO,GRADI,RADIUS,SLOLD,XMAX,ERXMIN,
     $                 VDOP,RSTAR,ENTOT,EN,NF,XLAMBDA,ELEVEL,NOM,
     $                 NDIM,N,NCHARG,WEIGHT,EINST,LINE)
     
C*******************************************************************************
C***  DETERMINATION OF THE CMF-FREQUENCIES XRED, XBLUE WHICH CONFINE THE
C***  OPTICALLY THICK LINE CORES
C*******************************************************************************

      USE MOD_LIOP
      use MOD_COFREQ
      use MOD_PRILC

      implicit real*8(a-h,o-z)
      integer,intent(in) :: JOBNUM
      COMMON / COMFUN / DELTAV,XMIN
      LOGICAL LINE(LASTIND)
      DIMENSION TAUMIN(0:2),GAMPRI(0:2)
      LOGICAL THIN
      DIMENSION VELO(ND),GRADI(ND),RADIUS(ND),ENTOT(ND)
      DIMENSION SLOLD(LASTIND),XRED(LASTIND),XBLUE(LASTIND)
      DIMENSION EINST(NDIM,NDIM),WEIGHT(NDIM),EN(NDIM),NCHARG(NDIM)
      DIMENSION XLAMBDA(NF)
      DIMENSION ELEVEL(NDIM)
      DIMENSION NOM(N)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      CHARACTER MODHEAD*104
      REAL*8,dimension(1) :: OPAL,ETAL ! Dimension ND=1 in LIOP

C***  BRANCH FOR GAMMAL=0 : LINE CORES ARE ZERO
      IF (GAMMAL .LE. .0 .AND. GAMMAR .LE. .0 ) THEN
      DO 1 IND=1,LASTIND
      XRED(IND)=.0
      XBLUE(IND)=.0
    1 CONTINUE
      RETURN
      ENDIF
     
C***  BRANCH FOR GAMMAL .NE. .0
C***  VL = VELOCITY IN DOPPLER UNITS ,  GL = GRADIENT IN DOPPLER UNITS
      VL=VELO(L)/VDOP
      GL=GRADI(L)/VDOP
      RL=RADIUS(L)
      DO 3 IND=1,LASTIND
      LOW=INDLOW(IND)
      NUP=INDNUP(IND)
      XRED(IND)=0.d0
      XBLUE(IND)=0.d0
      OPAL(1)=.0d0
C***  TAKE THE APPROPRIATE GAMMA : GAMMAR FOR RESONANCE LINE, GAMMAL ELSE
      IF (LOW .EQ. 1) THEN
c*** ws 4.7.03 try this to dampen the oscillations ...
c         GAMMA=GAMMAR*L
         GAMMA=GAMMAR
      ELSE
         IF ((NOM(LOW) .NE. NOM(LOW-1)) .OR. (NOM(LOW).EQ.NOM(LOW-1)
     $        .AND. NCHARG(LOW).NE.NCHARG(LOW-1))) THEN
c*** ws 4.7.03 try this to dampen the oscillations ...
c            GAMMA=GAMMAR*L
            GAMMA=GAMMAR
         ELSE
c*** ws 4.7.03 try this to dampen the oscillations ...
c            GAMMA=GAMMAL*L
            GAMMA=GAMMAL
         ENDIF
      ENDIF
      IF (GAMMA .LE. .0) GOTO 2
C***  RUDIMENTAL LINES : CORE IS SET TO ZERO
      IF (EINST(LOW,NUP) .EQ. -2.d0) GOTO 2
C***  BOUNDARY POINTS: OPTICALLY THIN (ZERO CORE)
      IF (L .EQ. ND) GOTO 2
      IF (L .EQ. 1) GOTO 2
C***  LINES WHICH WERE NOT TREATED IN THE RADIATION TRANSFER
C***  ARE ASSUMED TO HAVE ZERO CORE
      IF (.NOT.LINE(IND)) GOTO 2

      XLAM=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
C*** CALCULATE LINE OPACITY OPAL AT CURRENT DEPTH POINT AND LINE
      CALL LIOP(EINST(NUP,LOW),WEIGHT(LOW),WEIGHT(NUP),LOW,NUP,
     $1, XLAM, ENTOT(L:L), EN, RSTAR, OPAL, ETAL, VDOP, N)
C***  LASER SECURITY : CORE OF LASER LINES ARE SET TO ZERO
      IF (OPAL(1) .LE. .0) GOTO 2
C***  IF AMBIENT CONTINUUM IS OPTICALLY THICK :
C***  RADIAL DIRECTION -------------------------------------------------
      TAU=OPAL(1)/GL
      GDT=GAMMA /TAU
      DV1=VELO(1)/VDOP - VL
      DVND=VL - VELO(ND)/VDOP
      DELTAV=DMIN1(DV1,DVND)
      CALL COFREQ (XRR,XBR,XMAX,ERXMIN,GDT   ,THIN)
      IF (THIN) GOTO 2
C***  TRANSVERSAL DIRECTION  -------------------------------------------
      R1=RADIUS(1)
      TAUT=OPAL(1)*RL/VL
      GDT=GAMMA /TAUT
      DELTAV=SQRT(1.-RL*RL/R1/R1)*VELO(1)/VDOP
      CALL COFREQ (XRT,XBT,XMAX,ERXMIN,   GDT,THIN)
      IF (THIN) GOTO 2
C***  FIND THE MINIMUM CORE FROM BOTH DIRECTIONS -----------------------
      XRED(IND)=DMAX1(XRR,XRT)
      XBLUE(IND)=DMIN1(XBR,XBT)
C***  STORE OLD LINE SOURCE FUNCTION
    6 SLOLD(IND)=ETAL(1)/OPAL(1)

    2 CONTINUE
C***  STORE OPTICAL DEPTHS FOR PRINTOUT
      DO 4 I=0,2
      IF (IND .EQ. IPRILC+I) THEN
            TAUMIN(I)=OPAL(1)/DMAX1(GL,VL/RL)
            GAMPRI(I)=GAMMA
            ENDIF
    4 CONTINUE
    3 CONTINUE
     
      IF (IPRILC .GT. 0)
     $      CALL PRILC (IPRILC,LASTIND,XRED,XBLUE,TAUMIN,L,ND,ERXMIN,
     $                  MODHEAD,JOBNUM,GAMPRI)
     
      RETURN

      END SUBROUTINE

      END MODULE
