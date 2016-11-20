      MODULE MOD_BFCROSS
      contains
      SUBROUTINE BFCROSS (SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,NDIM,
     $                    XLAMBDA,ALPHA,SEXPO,AGAUNT,NOM,WAVARR,SIGARR,
     $					NFDIM)
C***********************************************************************
C***  THIS ROUTINE PREPARES AN ARRAY SIGMAKI WITH THE BOUND-FREE CROSS SECTIONS
C***  ( IN CM**2) TO AVOID UNNECCESSARY MULTIPLE CALCULATIONS
C***********************************************************************
      USE MOD_PHOTOCS_M 
      !implicit real*8(a-h,o-z)
      implicit none
      integer,intent(in) :: NF,N,NFDIM,NDIM
      integer,intent(in) :: NOM(N),NCHARG(NDIM)
      real*8,intent(in)  :: EION(NDIM),ELEVEL(NDIM),EINST(NDIM,NDIM)
      real*8,intent(out),dimension(NF,N)       :: SIGMAKI
      real*8,intent(in), dimension(*)          :: ALPHA,SEXPO
      character*8,intent(in) :: agaunt(NDIM)
      real*8,intent(in), dimension(NDIM,NFDIM) :: WAVARR,SIGARR
      real*8,intent(in), dimension(NF)         :: XLAMBDA
      integer :: NUP,K,LOW
      real*8:: SIGMATH,EDGE,WAVENUM
      DO 8 NUP=2,N
      DO 8 LOW=1,NUP-1
      !***  LEVELS MUST BELONG TO THE SAME ELEMENT
      IF (NOM(LOW) .NE. NOM(NUP)) GOTO 8
      !***  CHARGE DIFFERENCE MUST BE 1
      IF (NCHARG(NUP) .NE. NCHARG(LOW)+1 ) GOTO 8
      !***  UPPER LEVEL MUST BE GROUND STATE
      IF (NCHARG(NUP) .EQ. NCHARG(NUP-1)) GOTO 8
      !***  EINST = THRESHOLD CROSS SECTION IN 10**-18 CM**2
      SIGMATH=EINST(LOW,NUP)*1.d-18
      !***  EDGE = THRESHOLD ENERGY IN KAYSER *****
      EDGE=EION(LOW)-ELEVEL(LOW)
      DO 9 K=1,NF
        WAVENUM=1.E8/XLAMBDA(K)
        IF (WAVENUM.LT.EDGE) THEN
          SIGMAKI(K,LOW)=.0
        ELSE
          CALL PHOTOCS_M (SIGMAKI(K,LOW),SIGMATH,EDGE,WAVENUM,ALPHA,
     $                    SEXPO,AGAUNT,LOW,N,WAVARR,SIGARR,NDIM,NFDIM)
        ENDIF
    9 CONTINUE
    8 CONTINUE
      print *,N
      RETURN
      END SUBROUTINE
      END MODULE
