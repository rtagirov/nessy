      module MOD_PRIOPA
      contains
C**********  MODULNAME: PRIOPA    ******* 24/03/87  21.31.14.******    54 KARTEN
      SUBROUTINE PRIOPA (XLAM,K,ND,LSTEP,R,
     $ OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,JOBNUM,MODHEAD)
C***  PRINTOUT OF THE CONTINUUM OPACITIES ETC.
C***  CALLED BY WRCONT, COMO, FORMAL, 
      use MOD_TRADFUN
      IMPLICIT REAL*8(A-H,O-Z)
     
      parameter (one=1.d0)
      DIMENSION OPA(ND),ETA(ND),THOMSON(ND),R(ND),IWARN(ND)
      CHARACTER*10 MAINPRO(ND),MAINLEV(ND)
      CHARACTER MODHEAD*104, ihelp*5
     
      IF (IHELP.EQ.'IHELP') then
      IHELP='IHELP'
      PRINT 2,MODHEAD,JOBNUM
    2 FORMAT (10H1$$$$$$$$$,/1X,A,20X,'JOB NO.',I5,//,10X,
     $ 'CONTINUOUS OPACITY, EMISSIVITY AND SOURCE FUNCTION',
     $ /,10X,50('-'),//,
     $ ' FREQUENCY  DEPTH      OPACITY  THOMSON   OPTICAL    R(TAU=1)  '
     $ 'MAIN CONTRIBUTION      LASER    EMISSIVITY   SOURCE FUNCTION',/,
     $ '   INDEX    INDEX   (PER RSTAR) FRACTION   DEPTH               '
     $ 'PROCESS     LEVEL      WARNING    (...)        TRAD/KELVIN')
      else
	print 20
   20 FORMAT (/,
     $ ' FREQUENCY  DEPTH      OPACITY  THOMSON   OPTICAL    R(TAU=1)  '
     $ 'MAIN CONTRIBUTION      LASER    EMISSIVITY   SOURCE FUNCTION',/,
     $ '   INDEX    INDEX   (PER RSTAR) FRACTION   DEPTH               '
     $ 'PROCESS     LEVEL      WARNING    (...)        TRAD/KELVIN')
	endif
c    1 PRINT 3
c    3 FORMAT (1H )
      print '(A,F10.1)',' Wavelength= ',XLAM
      TAU=.0
      RTAU1=.0
      DO 6 L=1,ND
      IF (L.EQ.ND) GOTO 7
      TAUOLD=TAU
      TAU=TAU+0.5d0*(OPA(L)+OPA(L+1))*(R(L)-R(L+1))
      IF( TAUOLD.GE.one .OR. TAU.LT.one )  GOTO 7
      Q=(one-TAUOLD)/(TAU-TAUOLD)
      RTAU1=(one-Q)*R(L)+Q*R(L+1)
    7 IF(((L-1)/LSTEP)*LSTEP.NE.(L-1) .AND. L.NE.ND) GOTO 6
      IF (OPA(L).LE..0) THEN
            TRAD=.0
c		  write (100,*) k, l, TRAD
            ELSE
            S=ETA(L)/OPA(L)/(one-THOMSON(L))
            TRAD=TRADFUN (XLAM,S)
c		  write (100,*) k, l, TRAD
            ENDIF
      IF (L.EQ.ND) THEN
            PRINT 9,K,L,OPA(L),THOMSON(L),TAU,RTAU1
     $                  ,MAINPRO(L),MAINLEV(L),IWARN(L),ETA(L),TRAD
           WRITE(100,*) K,L,OPA(L),THOMSON(L),TAU,RTAU1
     $                  ,ETA(L),TRAD
            ELSE
            PRINT 5,K,L,OPA(L),THOMSON(L),TAUOLD
     $                  ,MAINPRO(L),MAINLEV(L),IWARN(L),ETA(L),TRAD
            
            WRITE(100,*) K,L,OPA(L),THOMSON(L),TAUOLD,RTAU1
     $                  ,ETA(L),TRAD
	ENDIF

c***  open output file to write opacity, emissivity, TRAD
	
c
    6 CONTINUE
	height1=(RTAU1-1.d0)*6.96D10/100000.d0
	PRINT *,' TAU=1 in km: ',height1
	write (33,*) xlam,height1
      RETURN
     
    5 FORMAT (I6,I10,1P,E15.3,0P,F7.3,2X,1P,E10.2, 10X ,3X,
     $                          A10,5X,A10,5X,A1,1P,E11.3,0P,F13.0)
    9 FORMAT (I6,I10,1P,E15.3,0P,F7.3,2X,1P,E10.2,0P,F10.3,3X,
     $                          A10,5X,A10,5X,A1,1P,E11.3,0P,F13.0)
     
      END subroutine
      end module
