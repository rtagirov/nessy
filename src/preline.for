      module MOD_PRELINE

      contains

      SUBROUTINE PRELINE(NUP,LOW,IND,N,LRUD,XLAM,ND,NFL,LINE,
     $                   XJLMEAN,ELEVEL,NL,EINST,INDNUP,INDLOW,LASTIND)

      use common_block

!     PREPARING SOME QUANTITIES FOR THE CONSIDERED LINE TRANSITION

      implicit real*8(a - h, o - z)
     
	  parameter (zero = 0.0d0)

      DIMENSION EINST(N, N)
      DIMENSION ELEVEL(N)

      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)
!      DIMENSION BMHO(NFL),BMNO(NFL),BMHI(NFL),BMNI(NFL)
      DIMENSION XJLMEAN(ND)
!      DIMENSION XJ(NFL,ND),XH(NFL,ND),XK(NFL,ND),XN(NFL,ND)
      character*7 name, line(nl)

!      logical :: lte_line
     
!     LOOP OVER ALL POSSIBLE LINE INDICES IND

      DO 16 IND = 1, LASTIND

   	  write(name, 1) 'LINE', ind
    1 FORMAT (A4,I3)

      IF (LINE(NL).EQ.NAME) GOTO 6
   16 ENDDO

      PRINT 8, LINE(NL)
    8 FORMAT (10X,'NON-FATAL ERROR: UNRECOGNIZED LINE OPTION CARD=',A8)
      PRINT*, 'The number of lines in DATOM file does not match the number of lines indicated in the CARDS file. See decetl.for'

      LINE(NL)='LINE000'
      NUP=0
      RETURN
     
    6 CONTINUE
C***  FIND THE LEVEL INDICES NUP, LOW
      NUP = INDNUP(IND)
      LOW = INDLOW(IND)

C***  RUDIMENTAL LINES (EINST(LOW,NUP)=-2.) WILL BE MARKED BY LRUD=0
      LRUD=1

      IF (EINST(LOW,NUP) .EQ. -2.d0) THEN
         LINE(NL)='LINE000'
         LRUD=0
         RETURN
      ENDIF

      XLAM=1.d8/(ELEVEL(NUP)-ELEVEL(LOW))
     
C***  INITIALIZE ALL ARRAYS FOR THE INTEGRATION OF MOMENTS
!      DO 11 KL=1,NFL
!      BMHO(KL)=zero
!      BMNO(KL)=zero
!      BMHI(KL)=zero
!   11 BMNI(KL)=zero
      do L = 1, ND; XJLMEAN(L) = zero; enddo
!      DO 12 KL=1,NFL
!      XJ(KL,L)=zero
!      XH(KL,L)=zero
!      XK(KL,L)=zero
!   12 XN(KL,L)=zero
     
      return

      end subroutine

      end module
