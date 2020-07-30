      module MOD_PRIRAT
      character*(*),private,parameter :: FMT_HEADER = 
     $  '(10H1$$$$$$$$$,/,1X,  A  ,20X,"JOB NO.",I5,//)'
      contains
C**********  MODULNAME: PRIRAT    ******* 24/03/87  21.31.18.******   130 KARTEN
      SUBROUTINE PRIRAT (ITNE,N,LEVEL,L,CRATE,RRATE,RATCO,EN,
     $           IFRRA,MODHEAD,JOBNUM,NETTO )
C***  OUTPUT OF RATE MATRIX ON FT09   *********************************

!      use MOD_AMBIPOLAR

      implicit none
      integer, intent(in) :: JOBNUM, ITNE, N, L, IFRRA, NETTO
      real*8,intent(in)    :: EN(N)
      real*8,intent(in)    :: CRATE(N,N)
      real*8,intent(in)    :: RRATE(N,N)
      CHARACTER,intent(in) :: LEVEL(N)*10
      CHARACTER,intent(in) :: MODHEAD*104
      real*8,intent(inout) :: RATCO(N,N)

      CHARACTER*12 NUMBERS(10)
      logical,save :: FIRST_TIME__ = .true.

      integer :: I,J,M,J1,J2  ! loop variables
      real*8 :: SCALE

C***  OPEN OUTPUT FILE AND WRITE HEADLINE FOR THE RATE PRINTOUT
      IF (L.EQ.1.OR.L.EQ.IFRRA) THEN
            OPEN (9,FILE='RATES')
            ENDIF
      WRITE (9,FMT_HEADER)  MODHEAD,JOBNUM
   !14 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,//)
     
C*******************************************************************************
C***  ADD RADIATIVE AND COLLISIONAL RATE COEFF. INTO MATRIX RATCO
      RATCO(:N,:N)=CRATE(:N,:N)+RRATE(:N,:N)
     
      IF (NETTO .NE. 1) THEN
        DO J1=1,N,10
          J2=MIN0(N,J1+9)
          WRITE (9,1) L, (LEVEL(J),J=J1,J2)
    1     FORMAT(//,20X,'TOTAL (COLL.+RAD.) RATE COEFF. AT DEPTH POINT NO.'
     $    ,I3,//,11X,10(2X,A10))
          WRITE (9,4)
    4     FORMAT (1X)
          DO I=1,N
            ENCODE (120,40,NUMBERS) (RATCO(I,J),J=J1,J2)
            where(NUMBERS(:10) == '    0.00E+00') NUMBERS(:10)='    0       '
            WRITE (9,42) LEVEL(I),(NUMBERS(M),M=1,10)
          ENDDO
        ENDDO
      ENDIF
C*******************************************************************************
C***  RATIO COLLISIONAL / RADIATIVE RATE COEFFICIENTS
      IF (NETTO .NE. 1) THEN
        ratio: DO J1=1,N,10
          J2=MIN0(N,J1+9)
          WRITE (9,61) L, (LEVEL(J),J=J1,J2)
  61      FORMAT (//,20X,'RATIO COLLISIONAL / RADIATIVE RATE COEFF. AT L='
     $          ,I3,//,11X,10(2X,A10))
          WRITE (9,4)
          L65: DO I=1,N
            L66: DO J=J1,J2
              IF (RRATE(I,J) .NE. .0) THEN 
                ENCODE (12,40,NUMBERS(1+J-J1)) CRATE(I,J)/RRATE(I,J)
                ELSE
                NUMBERS(1+J-J1)=' '
              ENDIF
            ENDDO L66
            WRITE (9,42) LEVEL(I),(NUMBERS(M),M=1,1+J2-J1)
          ENDDO L65
        ENDDO ratio
      ENDIF
     
C*******************************************************************************
C***  ADDING (ABSOLUTE) RATES INTO MATRIX RATCO
      forall(I=1:N)
        RATCO(I,:N)=EN(I)*(CRATE(I,:N)+RRATE(I,:N))
      endforall

      IF (NETTO .NE. 1) THEN
        total:DO J1=1,N,10
          J2=MIN0(N,J1+9)
          WRITE (9,21) L, (LEVEL(J),J=J1,J2)
  21      FORMAT(//,20X,'  TOTAL RATES                  AT DEPTH POINT NO.'
     $       ,I3,//,11X,10(2X,A10))
          WRITE (9,4)
          DO I=1,N
            ENCODE (120,40,NUMBERS) (RATCO(I,J),J=J1,J2)
            where(NUMBERS(:10) == '    0.00E+00') NUMBERS(:10)='    0       '
            WRITE (9,42) LEVEL(I),(NUMBERS(M),M=1,10)
          ENDDO
        ENDDO total
      ENDIF
C*******************************************************************************
      net: DO J1=1,N,10
        J2=MIN0(N,J1+9)
        WRITE (9,3) L, (LEVEL(J),J=J1,J2)
   3    FORMAT(//,20X,'  NET RATES                    AT DEPTH POINT NO.'
     $           ,I3,//,11X,10(2X,A10))
        WRITE (9,4)
        netPrint: DO I=1,N
          ENCODE (120,40,NUMBERS) (RATCO(I,J)-RATCO(J,I),J=J1,J2)
  40      FORMAT (1P,10E12.2)
          where(NUMBERS(:10) == '    0.00E+00') NUMBERS(:10)='    0       '
          WRITE (9,42) LEVEL(I),(NUMBERS(M),M=1,10)
  42      FORMAT (1X,A10,10A12)
        ENDDO netPrint
      ENDDO net
     
C*******************************************************************************
      L53: DO J1=1,N,10
        J2=MIN0(N,J1+9)
        WRITE (9,9) L, (LEVEL(J),J=J1,J2)
    9   FORMAT(//,20X,'   RELATIVE NET RATES          AT DEPTH POINT NO.'
     $    ,I3,//,11X,10(2X,A10))
        WRITE (9,4)
        DO I=1,N
          SCALE=0.
          DO J=1,N
            SCALE=SCALE+ABS(RATCO(I,J)-RATCO(J,I))
          ENDDO
          SCALE=SCALE/2.
          IF (SCALE.LE.0.) SCALE=1. ! i.e. if all(RATCO(I,:N)-RATCO(:N,I))==0
          ENCODE (120,40,NUMBERS) ((RATCO(I,J)-RATCO(J,I))/SCALE,J=J1,J2)
          where(NUMBERS(:10) == '    0.00E+00') NUMBERS(:10)='    0       '
          WRITE (9,42) LEVEL(I),(NUMBERS(M),M=1,10)
        ENDDO
      ENDDO L53

      !**********************************************************************
      !*** Print the RATCO,CRATE,RRATE as data files for further analysis
      FirstTime: IF (FIRST_TIME__) THEN
        !*** Open Files and print Headers and Element Names
        !*** Start INIT
        !*** Save the Element/Level names
        open (12,FILE='RATES.NAMES',status='replace')
        write (12,FMT_HEADER) MODHEAD,JOBNUM
        write(12,'("     N:",I6)'),N
        write(12,'(A10)')(LEVEL(J),J=1,N)
        close(12);
        !*** the others
        open(12,FILE='RATES.RATCO',status='replace')
        open(13,FILE='RATES.CRATE',status='replace')
        open(14,FILE='RATES.RRATE',status='replace')
        write(12,'("HEADER:7")');  write (12,FMT_HEADER) MODHEAD,JOBNUM
        write(13,'("HEADER:7")');  write (13,FMT_HEADER) MODHEAD,JOBNUM
        write(14,'("HEADER:7")');  write (14,FMT_HEADER) MODHEAD,JOBNUM
        write(12,'("RATCO:",I5,I5," *")'),N,N
        write(13,'("CRATE:",I5,I5," *")'),N,N
        write(14,'("RRATE:",I5,I5," *")'),N,N
      ENDIF FirstTime
      write(12,'("    L:",I5)'),L
      write(13,'("    L:",I5)'),L
      write(14,'("    L:",I5)'),L
      do I=1,N
        write(12,'(X,pe20.12E3$)'),(RATCO(I,J),J=1,N); write(12,'()')
        write(13,'(X,pe20.12E3$)'),(CRATE(I,J),J=1,N); write(13,'()')
        write(14,'(X,pe20.12E3$)'),(RRATE(I,J),J=1,N); write(14,'()')
      enddo
      FIRST_TIME__ = .false.
      RETURN
      END subroutine
      end module
