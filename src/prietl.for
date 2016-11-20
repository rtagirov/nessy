      module MOD_PRIETL
      contains
C**********  MODULNAME: PRIETL    ******* 24/03/87  19.46.17.******    78 KARTEN
      SUBROUTINE PRIETL (IND,XLAM,ND,OPA,OPAL,ETA,ETAL,RADIUS,JOBNUM,
     $      LSOPA,XJLMEAN,MODHEAD,ASF,BSF,T,ETLKEY)
      use MOD_BNUE
      use MOD_TRADFUN

      implicit real*8(a-h,o-z)

      DIMENSION OPA(ND),OPAL(ND),ETA(ND),ETAL(ND),RADIUS(ND)
      CHARACTER MODHEAD*104
      DIMENSION XJLMEAN(ND),ASF(ND),BSF(ND),T(ND)
      LOGICAL ETLKEY

C***  WPI = SQRT(PI)
      DATA WPI /1.7724538509055d0/
     
      IF (AHELP.EQ.5HIHELP) GOTO 1
      AHELP=5HIHELP
      PRINT 2,MODHEAD,JOBNUM
     
    2 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,//,10X,
     $ 'LINE OPACITY, EMISSIVITY AND SOURCE FUNCTION',
     $ /,10X,44('-'),//,
     $ ' LINE  L LINE OP.  LINE/CONT. TOT.OPT.DEPTH    R    ',
     $ 'LINE EMISS.  LINE SOURCE  TOT. SOURCE  ',
     $ 'J (FORMAL SOLUTION)  ETLA SOURCE F.   ',/,
     $ ' IND.       /RSTAR               (LINE CENTER) (TAU=1) ',
     $ '              TRAD/K        TRAD/K     ',
     $ '(ERG/CM+2)   TRAD/K   SCATT.   /B(T)  ',//)
     
    1 TAU=.0
      RTAU1=.0
      DO 6 L=1,ND
      OPALC=OPAL(L)/WPI
      ETALC=ETAL(L)/WPI
      OPATOT=OPA(L)+OPALC
      ETATOT=ETA(L)+ETALC
      IF (L.EQ.ND) GOTO 7
      TAUOLD=TAU
      OPATOTP=OPA(L+1)+OPAL(L+1)/WPI
      TAU=TAU+0.5d0*(OPATOT+OPATOTP)*(RADIUS(L)-RADIUS(L+1))
      IF( TAUOLD.GE.1. .OR. TAU.LT.1. )  GOTO 7
      Q=(1.d0-TAUOLD)/(TAU-TAUOLD)
      RTAU1=(1.-Q)*RADIUS(L)+Q*RADIUS(L+1)
    7 IF(((L-1)/LSOPA)*LSOPA.NE.(L-1) .AND. L.NE.ND) GOTO 6
      S=ETALC/OPALC
      TRADLIN=TRADFUN (XLAM,S)
      S=ETATOT/OPATOT
      TRADTOT=TRADFUN (XLAM,S)
     
      IF (ETLKEY) THEN
C***     BRANCH FOR LINES TREATED BY ETLA
         TRADJ=TRADFUN (XLAM,XJLMEAN(L))
         IF ( L .EQ. ND ) THEN
            PRINT 9,IND,L,OPALC,OPALC/OPA(L),TAU,RTAU1,ETALC,TRADLIN,
     $      TRADTOT,XJLMEAN(L),TRADJ,ASF(L),BSF(L)/BNUE(XLAM,T(L))
            ELSE
            PRINT 8,      L,OPALC,OPALC/OPA(L),TAUOLD,   ETALC,TRADLIN,
     $      TRADTOT,XJLMEAN(L),TRADJ,ASF(L),BSF(L)/BNUE(XLAM,T(L))
            ENDIF
         ELSE
C***     BRANCH FOR LINES TREATED ONLY BY FORMAL SOLUTION
         IF ( L .EQ. ND ) THEN
            PRINT 4,IND,L,OPALC,OPALC/OPA(L),TAU,RTAU1,ETALC,TRADLIN,
     $      TRADTOT
            ELSE
            PRINT 3,      L,OPALC,OPALC/OPA(L),TAUOLD,   ETALC,TRADLIN,
     $      TRADTOT
            ENDIF
         ENDIF
    6 CONTINUE
     
      RETURN
     
    9 FORMAT (I5,I3,1P,E9.2,2X,E9.2,6X,E8.2,0P,F9.3,3X,1P,E8.2,0P,
     $    F13.0,F13.0,1P,E13.2,0P,F8.0,F8.3,F8.3,/)
    8 FORMAT (   8X,I3,1P,E9.2,2X,E9.2,6X,E8.2,   12X    ,1P,E8.2,0P,
     $    F13.0,F13.0,1P,E13.2,0P,F8.0,F8.3,F8.3)
    4 FORMAT (I5,   I3,1P,E9.2,2X,E9.2,6X,E8.2,0P,F9.3,3X,1P,E8.2,0P,
     $    F13.0,F13.0,'    ---  ONLY FORMAL SOLUTION  ---',/)
    3 FORMAT (   8X,I3,1P,E9.2,2X,E9.2,6X,E8.2,   12X    ,1P,E8.2,0P,
     $    F13.0,F13.0)
     
      END subroutine
      end module
