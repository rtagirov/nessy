      module MOD_PRIGH
      contains
      SUBROUTINE PRIGH (LPRIH,ND,RADIUS,HTOT,GTOT,ETOT,TEFF,ENTOT,RNE,
     $                  RSTAR,T,VELO,GRADI,ATMASS,ABXYZ,NATOM)

      USE CONSTANTS, ONLY: CLIGHT_CGS,PI
      IMPLICIT NONE
      integer,intent(in   ) :: LPRIH,ND,NATOM
      real*8, intent(in   ) :: RADIUS,HTOT,GTOT,ETOT,ENTOT
      real*8, intent(in   ) :: RNE,RSTAR,T,VELO,GRADI,ATMASS,ABXYZ
      real*8, intent(inout) :: TEFF
      integer :: L,N
      real*8  :: CONST,TL,WEI,AMUMEA
        
      real*8,parameter :: four=4.d0
      DIMENSION RADIUS(ND),HTOT(ND),GTOT(ND),ETOT(ND),ENTOT(ND),RNE(ND),
     $          T(ND),VELO(ND),GRADI(ND),
     $          ATMASS(NATOM),ABXYZ(NATOM)
      DIMENSION AMUMEA(70)
      !real*8,parameter :: PI=3.141592654d0
      real*8,parameter :: PISIG=5.5411d+4
      real*8,parameter :: AMU=1.66053111d-24,C=CLIGHT_CGS

C  PURE HELIUM
      CONST=four*PI/C/four/AMU/RSTAR
      IF (TEFF.LE.0..AND.HTOT(1).GT.0.)
     &             TEFF=(four*HTOT(1)*RADIUS(1)*RADIUS(1)*PISIG)**0.25
      PRINT 10, TEFF
   10 FORMAT (///,20X,'ASTROPHYSICAL FLUX AS A FUNCTION OF DEPTH',/,
     $            20X,'=========================================',
     $   //,20X,'TEFF = ',F20.0
     $   //,20X,'DEPTH INDEX    FLUX (ERG/CM2/S)    T (KELVIN)',
     $          '            T/TEFF'/)

      DO 1 L=1,ND
      IF(((L-1)/LPRIH)*LPRIH.NE.(L-1) .AND. L.NE.ND) GOTO 1
      IF (HTOT(L).LE.0.) THEN
      TL=0.
      ELSE
      TL=(four*HTOT(L)*RADIUS(L)*RADIUS(L)*PISIG)**0.25
      ENDIF
      PRINT 13,L,four*HTOT(L),TL,TL/TEFF,GTOT(L)*CONST/ENTOT(L),
     &           ETOT(L)
   13 FORMAT (20X,I10,1PE20.5,0PF15.0,F20.5,2E20.3)
    1 CONTINUE
C
C PRINT OUT FOR TESTS
      PRINT *,' NATOM'
      PRINT *,NATOM
      PRINT *,' ATMASS'
      PRINT *,ATMASS
      PRINT *,' ABXYZ'
      PRINT *,ABXYZ
      PRINT *,' RSTAR'
      PRINT *,RSTAR
      PRINT *,' TEFF'
      PRINT *,TEFF
      PRINT *,' RADIUS'
      PRINT *,RADIUS
      PRINT *,' HTOT'
      PRINT *,HTOT
      PRINT *,' GTOT'
      PRINT *,GTOT
      PRINT *,' ETOT'
      PRINT *,ETOT
      PRINT *,' ENTOT'
      PRINT *,ENTOT
      PRINT *,' RNE'
      PRINT *,RNE
      PRINT *,' T'
      PRINT *,T
      PRINT *,' VELO'
      PRINT *,VELO
      PRINT *,' GRADI'
      PRINT *,GRADI
      IF (ND.GT.70) THEN
         PRINT *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         PRINT *,' ND.GT.70 - REST OF SUBROUTINE SKIPED'
         PRINT *,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         RETURN
         ENDIF
C
C MEAN ATOMIC WEIGHT
      WEI=0.
      DO 20 N=1,NATOM
   20    WEI=WEI+ATMASS(N)*ABXYZ(N)
      PRINT *,' MEAN ATOMIC WEIGHT:',WEI
C      DO 11 L=1,ND
C      AMUMEA(L)=

      RETURN
      END subroutine
      end module
