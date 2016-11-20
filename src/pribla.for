      module MOD_PRIBLA
      contains
      SUBROUTINE PRIBLA (LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,
     $	                 SCAFAC,ABSFAC)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION XLAMBDA(NF),ENTOT(ND),DENS(21)
      CHARACTER MODHEAD*104
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      DIMENSION SCAFAC(ND,NF),ABSFAC(ND,NF)

      IF (LBLAON.NE.1) THEN
      PRINT 100
  100 FORMAT(1X,10('*'),' LINE BLANKETING NOT ACTIVE ',10('*'))
      ELSE
      PRINT 1,MODHEAD,JOBNUM,NBINW
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,
     $            //,10X,'BIN WIDTH: ',I5,' A',
     $            //,10X,'LINE SCATTERING FACTORS',/,
     $               10X,'LINE ABSORPTION FACTORS',/,10X,22('-'),/,
     $ 10X,' FREQUENCY               LOG ELECTRON DENSITIES',/)
     
      nde=nd-6
      NDP=20
c      nda=23
c      STEP=(NDe-nda)/(ndp-1)
      step=2.
      nda=nde-nint(step*(ndp-1))
      if (nda.le.10) nda=10
      IF (STEP.LE.1.) then
         print *,' nf,nd ',nf,nd
         print *,' nda,nde,step ',nda,nde,step
         STop ' pribla - step.le.1.'
      endif

      DO 2 L=1,NDP
    2 DENS(L)=ENTOT(NINT((L-1)*STEP)+nda)

      PRINT 3,(NINT((L-1)*STEP)+nda,L=1,NDP),(log10(DENS(L)),L=1,NDP)
    3 FORMAT(20x,20I5/20X,20F5.1/20x,100('-'))

      DO 4 K=1,NF
      IF (XLAMBDA(K).GT.ALMAX) GOTO 4
      PRINT 5,K,XLAMBDA(K),
     &      (SCAFAC(NINT((L-1)*STEP)+nda,K),L=1,NDP),
     &      (ABSFAC(NINT((L-1)*STEP)+nda,K),L=1,NDP)
    5 FORMAT (I7,F13.2,20F5.1,/,20X,20F5.1)
    4 CONTINUE
      ENDIF
      RETURN
      END subroutine
      end module
