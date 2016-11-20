      module MOD_PRILC
      contains
C**********  MODULNAME: PRILC     ******* 24/03/87  21.29.41.******    29 KARTEN
      SUBROUTINE PRILC (IPRILC,INDLAST,XRED,XBLUE,TAUMIN,L,ND,ERXMIN,
     $                  MODHEAD,JOBNUM,GAMPRI)
C***  PRINTOUT OF SCHARMER LINE CORES  *********************************
     
      use MOD_ERF_INF
      implicit real*8(a-h,o-z)

      DIMENSION XRED(2),XBLUE(2),NCHARG(2)
      DIMENSION TAUMIN(0:2),AMP(0:2),GAMPRI(0:2)
      CHARACTER MODHEAD*104
     
	IPR=IPRILC
      IF (IPRILC.gt.INDLAST-2) IPR=INDLAST-2
C***  PRINT THE HEADER OF THE TABLE
      IF (L .EQ. 1) THEN
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,//)
      PRINT 10,(IPR+I,GAMPRI(I),I=0,1)
   10 FORMAT (       10X,'SCHARMER-CORES',//,5X,
     $   2('   -----  LINE',I3,'  ---  GAMMA=',F5.1,'  -----'),//,
     $'  L   ',2('      XRED           XBLUE      TAUMIN      AMP',
     $'          '),/)
c      PRINT 10,(IPR+I,GAMPRI(I),I=0,2)
c   10 FORMAT (       10X,'SCHARMER-CORES',//,5X,
c     $   3('   -----  LINE',I3,'  ---  GAMMA=',F5.1,'  -----'),//,
c     $ '  L   ',3('  XRED   XBLUE  TAUMIN      AMP          '),/)
      ENDIF
     
      DO 2 I=0,2
      FC=(ERF_INF(XBLUE(IPRILC+I))-
     &    ERF_INF(XRED(IPRILC+I)))/(1.-2.*ERXMIN)
      AMP(I)=1./(1.-FC)
    2 CONTINUE
      PRINT 12, L,
     $        (XRED(IPRILC+I),XBLUE(IPRILC+I),TAUMIN(I),AMP(I),I=0,1)
c     $        (XRED(IPRILC+I),XBLUE(IPRILC+I),TAUMIN(I),AMP(I),I=0,2)
   12 FORMAT (I3,2X,2(2E15.6,1P,G11.2,G11.2,0P,6X))
c   12 FORMAT (I3,2X,3(2F7.2,1P,G11.2,G11.2,0P,6X))
      RETURN
      END subroutine
      end module
