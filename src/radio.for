C**********  MODULNAME: RADIO     ******* 06/08/87  20.14.58.******    99 KARTEN
      module MOD_RADIO
      contains
      SUBROUTINE RADIO (NDIM,N,ENLTE,TL,WEIGHT,NCHARG,EION,ELEVEL,EINST,
     $       RRATE,XLAMBDA,FWEIGHT,XJC,NF,L,XJL,ND,LASTIND,SIGMAKI,NOM)
     
C*******************************************************************************
C***  RADIATIVE RATES ARE CALCULATED AT ONE DEPTH POINT L FROM THE GIVEN
C***  RADIATION FIELD AND ATOMIC CROSS SECTIONS
C***  RADIATIVE TRANSITION RATES STORED IN MATRIX RRATE  ***********
C***  CALLED BY NONLTEPOP
C*******************************************************************************
      use MOD_XRUDI
      use MOD_ISRCHFGT
      use UTILS
      use IFCORE
      implicit real*8(a-h,o-z)
      real*8,intent(in   ) :: XJL(ND,LASTIND)
      DIMENSION RRATE(NDIM,NDIM),EINST(NDIM,NDIM)
      DIMENSION NCHARG(NDIM),ELEVEL(NDIM)
      DIMENSION EION(NDIM),ENLTE(NDIM),WEIGHT(NDIM)
      DIMENSION NOM(N)
      DIMENSION XJC(ND,NF)
      DIMENSION SIGMAKI (NF,N)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1 / 1.4388D0 /
C***  C2 = 2 * H * C    ( H AND C IN CGS UNITS )
      DATA C2 / 3.9724D-16 /
C***  C3 = 4 * PI / H / C ( CGS UNITS )
      DATA C3 / 0.06327D18 /
     
C***  LOOP OVER ALL TRANSITIONS  ---------------------------------------
      IND=0
      DO 1 NUP=2,N
      DO 1 LOW=1,NUP-1
      IF (NOM(LOW) .NE. NOM(NUP)) GOTO 14
      IF (NCHARG(LOW).NE.NCHARG(NUP)) GOTO 8
     
C***  LINE TRANSITION   ************************************************
      IND=IND+1
      WAVENUM=ELEVEL(NUP)-ELEVEL(LOW)
      W3=WAVENUM*WAVENUM*WAVENUM

      !***  CHECK WHETHER THIS TRANSITION IS ONLY RUDIMENTAL
      IF (EINST(LOW,NUP) .EQ. -2.) THEN
        !***  LINE IS RUDIMENTAL
        CALL XRUDI (XJ,WAVENUM,XJC,XLAMBDA,ND,NF,L)
      ELSE
        !*** LINE IS EXPLICIT
        XJ=XJL(L,IND)
      ENDIF
      if(XJ<=0.) then
        print '("radio: XJ = ",e12.3," EINST(",i0,",",i0,") = ",e12.3)',
     &           XJ,LOW,NUP,EINST(LOW,NUP)
        if(EINST(LOW,NUP)==-2.) then
        print '("radio: WV# = ",e12.3," XJC(",i0,"/",i0,",:",i0,")")',
     &           WAVENUM,L,ND,NF
        print '("radio: XJC = ",5(e12.3,X))',XJC
        endif
        !call tracebackqq('',USER_EXIT_CODE=-1_4)
      endif
      call assert(XJ>0., 'XJ==0')
C***  CALCULATION OF LINE RATES
   11 CONTINUE
      EMINDU=EINST(NUP,LOW)*XJ/C2/W3
c	print *,NUP,LOW,'RADIO,NUP,LOW, EMINDU=', EMINDU
c	PRINT *,NUP,LOW,'RADIO,NUP,LOW, EINST(NUP,LOW)', EINST(NUP,LOW)
      RRATE(LOW,NUP)=EMINDU*WEIGHT(NUP)/WEIGHT(LOW)
      RRATE(NUP,LOW)=EINST(NUP,LOW)+EMINDU
c	print *,NUP,LOW,'RADIO: rrate(nup,low)', rrate(nup,low)
      GOTO 1
     
    8 CONTINUE
C***  CHARGE DIFFERENCE MUST BE 1
      IF (NCHARG(NUP) .NE. NCHARG(LOW)+1 ) GOTO 14
C***  UPPER LEVEL MUST BE A GROUND STATE
      IF (NCHARG(NUP) .NE. NCHARG(NUP-1)+1) GOTO 14
     
C***  CONTINUUM TRANSITION    **************************************************
C***  SIGMAKI = PRECALCULATED CROSS SECTION IN CM**2
C***  EDGE = THRESHOLD ENERGY IN KAYSER *******
      EDGE = EION(LOW)-ELEVEL(LOW)
      EDGELAM=1.d8/EDGE
     
C***  RATE INTEGRAL
      SUM=.0D0
      REC=.0D0
C***  FIND EDGE FREQUENCY INDEX
      NFEDGE=ISRCHFGT (NF,XLAMBDA,1,EDGELAM) - 1
      DO 2 K=1,NFEDGE
      WAVENUM=1.E8/XLAMBDA(K)
      W3=WAVENUM*WAVENUM*WAVENUM
      XJCLK=XJC(L,K)
      SIGMA=SIGMAKI(K,LOW)
C***  IONIZATION ***********************
      SUM=SUM+SIGMA*XJCLK*FWEIGHT(K)/WAVENUM
C***  RECOMBINATION   **************
      REC=REC +
     $    SIGMA*(XJCLK+C2*W3)*FWEIGHT(K)*EXP(-C1*WAVENUM/TL)/WAVENUM
    2 CONTINUE
     
C***  CALCULATION OF BOUND-FREE RATES
      RRATE(LOW,NUP)=C3*SUM
      RRATE(NUP,LOW)=C3*REC*ENLTE(LOW)/ENLTE(NUP)
      GOTO 1
     
C***  LEVELS BELONG TO DIFFERENT ELEMENTS,
C***  CHARGE DIFFERENCE NOT 1  OR UPPER LEVEL NO GROUND STATE: ZERO RATE
   14 RRATE(LOW,NUP)=.0D0
      RRATE(NUP,LOW)=.0D0
    1 CONTINUE
	
cmh	OPEN (UNIT=9,FILE='RRATE')
cmh	DO J=1,29
cmh	DO I=1,29
cmh	WRITE (UNIT=9,FMT=*),J,L,RRATE(J,I) 
cmh	END DO
cmh	END DO
cmh	CLOSE (UNIT=9)
C***  ENDLOOP  ---------------------------------------------------------
     
C***  DIAGONAL ELEMENTS ARE SET TO ZERO
      DO 9 J=1,N
    9 RRATE(J,J)=.0
     
      RETURN
      END SUBROUTINE
      END MODULE
