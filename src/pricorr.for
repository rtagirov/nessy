      MODULE MOD_PRICORR

      CONTAINS

      SUBROUTINE PRICORR(POPNUM, POP1, LEVEL, N, ND, MODHEAD, LSPOP, CORMAX_FIN,
     $                   NCHARG, RNE, JOBNUM, REDUCE, GAMMAL, GAMMAR, T,
     $                   ELEMENT, NATOM, NFIRST, NLAST)

!     PRINTOUT OF CORRECTION FACTORS OF POPULATION NUMBERS
!     RELATIVE TO THE LAST ITERATION

      use mod_tictoc

      use file_operations
      use common_block

      IMPLICIT REAL*8(A - H, O - Z)

      INTEGER, INTENT(IN) :: JOBNUM

      DIMENSION POPNUM(ND, N), POP1(ND, N)
      DIMENSION NCHARG(N), RNE(ND), T(ND)
      DIMENSION NFIRST(NATOM), NLAST(NATOM)

      REAL*8, ALLOCATABLE :: ABXYZ_new(:)

      CHARACTER LEVEL(N)*10
      CHARACTER ELEMENT(NATOM)*10
      CHARACTER MODHEAD*100
      CHARACTER PRILINE*130, NUMBER*12

      REAL*8,  DIMENSION(NATOM) :: CORMAX_ELE
      REAL*8,  DIMENSION(N) ::     CORMAX_LEV

      INTEGER, DIMENSION(NATOM) :: DCM_ELE

      REAL*8,  INTENT(OUT) ::      CORMAX_FIN

      REAL*8 ::                    CORMAX_ALL
      REAL*8 ::                    CORMAX_NHE
      REAL*8 ::                    CORMAX_NHN

      INTEGER ::                   DCM_ALL
      INTEGER ::                   DCM_NHE
      INTEGER ::                   DCM_NHN

      REAL*8 ::                    CM
 
      INTEGER ::                   DCM, DCM_FIN

      CHARACTER(:), ALLOCATABLE :: CONV_FILE

      CHARACTER(LEN = 10) ::       ELE

      CHARACTER(LEN = 8) ::        time

      INTEGER, ALLOCATABLE, DIMENSION(:) :: DEPTHS

      INTEGER :: DEPTHS_DIM, DEPTH_IND

      LOGICAL :: ALLE, NOHE, NOHN

      LOGICAL NEGPOPL
C***  COMMON /GIIIERR/ COUNTS THE ERROR MESSAGES FROM SUBR. GAUNTFF
      COMMON /GIIIERR/  NTUP,NTLOW,NFUP,NFLOW,NHELP
C***  NEGINTL COUNTS THE ERROR MESSAGES  FROM SUBR. LINPOP

!      START = LAMBDA_ITER .EQ. 1! .AND. .NOT. OLDSTART .AND. .NOT. LBKG

      CORMAX_ELE(:) = 0.0D0
      CORMAX_LEV(:) = 0.0D0

      DCM_ELE(:) =    0

      IF (LSPOP .GE. 1) PRINT 30
   30 FORMAT (10H1$$$$$$$$$)
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (              /,1X,  A  ,20X,'JOB NO.',I5,//,10X,
     $ 'RATIO OF POPULATION NUMBERS TO THOSE OF THE LAST ITERATION:')
      IF (LSPOP.LT.1) GOTO 8
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2,(LEVEL(J),J=J1,J2)
    2 FORMAT (//,' DEPTH',10(2X,A10))
      IF (N.LT.J2) PRINT 5
    5 FORMAT (1X)
      PRINT 5
     
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      ENCODE (130,6,PRILINE) L
    6 FORMAT (I6)
     
      DO 12 J=J1,J2
      IF (POP1(L,J) .NE. .0) THEN
            ENCODE (12,11,NUMBER) POPNUM(L,J)/POP1(L,J)
   11       FORMAT (F12.4)
            ELSE
            NUMBER='    INFINITE'
            ENDIF
      I=7+(J-J1)*12
      PRILINE(I:I+11)=NUMBER
   12 CONTINUE
     
      PRINT 13,PRILINE
   13 FORMAT (A)
    3 CONTINUE
     
      IF (J2.EQ.N) GOTO 8
      J1=J1+10
      GOTO 4
    8 CONTINUE

!----------------------------------------------------------------------------------
! CALCULATION BLOCK START
!----------------------------------------------------------------------------------

      DEPTHS_DIM = ND - 3

      ALLOCATE(DEPTHS(DEPTHS_DIM))

      DEPTH_IND = 1

      DO L = 3, ND - 1; DEPTHS(DEPTH_IND) = L; DEPTH_IND = DEPTH_IND + 1; ENDDO

      CORMAX_ALL = 0.0D0
      CORMAX_NHE = 0.0D0
      CORMAX_NHN = 0.0D0

      DCM_ALL = 0
      DCM_NHE = 0
      DCM_NHN = 0

      DO NA = 1, NATOM

         CALL CALC_CONV_IND(NFIRST(NA), NLAST(NA), ELEMENT(NA), DEPTHS, DEPTHS_DIM, POPNUM, POP1, LEVEL, N, ND, CM, DCM)

         CORMAX_ELE(NA) = CM

         DCM_ELE(NA) = DCM

         ALLE = CORMAX_ELE(NA) .GT. CORMAX_ALL
         NOHE = CORMAX_ELE(NA) .GT. CORMAX_NHE .AND. TRIM(ADJUSTL(ELEMENT(NA))) .NE. 'HELIUM'

         IF (ALLE) THEN; CORMAX_ALL = CORMAX_ELE(NA); DCM_ALL = DCM_ELE(NA); ENDIF

         IF (NOHE) THEN; CORMAX_NHE = CORMAX_ELE(NA); DCM_NHE = DCM_ELE(NA); ENDIF

         DO J = NFIRST(NA), NLAST(NA)

            CALL CALC_CONV_IND(J, J, LEVEL(J), DEPTHS, DEPTHS_DIM, POPNUM, POP1, LEVEL, N, ND, CM, DCM)

            CORMAX_LEV(J) = CM

            NOHN = CORMAX_LEV(J) .GT. CORMAX_NHN .AND.
     $             TRIM(ADJUSTL(ELEMENT(NA))) .NE. 'HELIUM' .AND.
     $             LEVEL(J) .NE. 'NeII-----1'

            IF (NOHN) THEN; CORMAX_NHN = CORMAX_LEV(J); DCM_NHN = DCM; ENDIF

         ENDDO

      ENDDO

!      CORMAX_FIN = CORMAX_ELE(1); DCM_FIN = DCM_ELE(1)
      CORMAX_FIN = CORMAX_NHN; DCM_FIN = DCM_NHN
!      CORMAX_FIN = CORMAX_ALL; DCM_FIN = DCM_ALL

!----------------------------------------------------------------------------------
! PRINT OUT BLOCK START
!----------------------------------------------------------------------------------

      WRITE(*, '(/, A)') '---------------------------------------------------------
     $--------------------------------------------------------'

      WRITE(*, 10002) DEPTHS(1), DEPTHS(DEPTHS_DIM), CORMAX_ELEC, LOG10(CORMAX_ELEC)

      WRITE(*, '(/, A)') '---------------------------------------------------------
     $--------------------------------------------------------'

!==========================================================================================

      if (full_conv_print) then

         IF (LAMBDA_ITER .EQ. 1) THEN; CALL MKDIR(CONV_DIR); CALL CLEAN_DIR(CONV_DIR); ENDIF

         DO NA = 1, NATOM

            CONV_FILE = CONV_DIR//ELEMENT(NA)

            CALL OPEN_TO_APPEND(999 + NA, CONV_FILE)

            IF (LAMBDA_ITER .EQ. 1) THEN

                WRITE(999 + NA, '(3x, A, 5x, A, 6x, A, 1x, $)') 'LI', 'MAX', 'DI'

                DO J = NFIRST(NA), NLAST(NA); WRITE(999 + NA, '(2x, A10, $)') LEVEL(J); ENDDO

                WRITE(999 + NA, '(/)')

            ENDIF

            WRITE(999 + NA, '(I5, 2x, ES9.3, 2x, I3, $)') LAMBDA_ITER, CORMAX_ELE(NA), DCM_ELE(NA)

            DO J = NFIRST(NA), NLAST(NA)

               IF (J .NE. NLAST(NA)) WRITE(999 + NA, '(3x, ES9.3, $)') CORMAX_LEV(J)
               IF (J .EQ. NLAST(NA)) WRITE(999 + NA, '(3x, ES9.3)')    CORMAX_LEV(J)

            ENDDO

            CLOSE(UNIT = 999 + NA)

         ENDDO

!==========================================================================================

         CALL OPEN_TO_APPEND(1000 + NATOM, CONV_DIR//'ELEM')

         IF (LAMBDA_ITER .EQ. 1) THEN

            WRITE(1000 + NATOM, '(3x,A,$)') 'LI'

            DO NA = 1, NATOM

               ELE = ELEMENT(NA)

               IF (NA .EQ. 1) WRITE(1000 + NATOM, '(5x,A,$)') ELE(1 : 3)

               IF (NA .NE. 1) WRITE(1000 + NATOM, '(8x,A,$)') ELE(1 : 3)

            ENDDO

            WRITE(1000 + NATOM, '(8x,A,/)') 'ELE'

         ENDIF

         WRITE(1000 + NATOM, '(I5,$)') LAMBDA_ITER

         DO NA = 1, NATOM; WRITE(1000 + NATOM, '(2x,ES9.3,$)') CORMAX_ELE(NA); ENDDO

         WRITE(1000 + NATOM, '(2x,ES9.3)'), CORMAX_ELEC

         CLOSE(UNIT = 1000 + NATOM)

!==========================================================================================

         CALL OPEN_TO_APPEND(1001 + NATOM, CONV_DIR//'ALL')

         IF (LAMBDA_ITER .EQ. 1) THEN

            WRITE(1001 + NATOM, '(3x,A,7x,A,6x,A,3(8x,A),5x,A,/)') 'LI', 'TIME', 'ALL', 'NHE', 'NHN', 'HYD', 'DCMF'

         ENDIF

         time = trim(adjustl(writeTOC()))

         if (len_trim(time) == 5) time = '00:'//time

         WRITE(1001 + NATOM, '(I5,5x,A8,4(2x,ES9.3),2x,I3)') LAMBDA_ITER, time, CORMAX_ALL,
     $                                                       CORMAX_NHE, CORMAX_NHN, CORMAX_ELE(1), DCM_FIN

      else

         CALL OPEN_TO_APPEND(19462, 'conv.out')
        
         time = trim(adjustl(writeTOC()))

         if (len_trim(time) == 5) time = '00:'//time

         WRITE(19462, '(A8,2x,ES9.3)') time, cormax_fin

         if (hpc_run .and. cormax_fin >= 0.1 .and. lambda_iter > 40) then

            call system('rm BROYDEN')
            call system('rm fort.99')
            call system('rm MODFILE')
            call system('rm MODHIST')
            call system('rm molconc.out.inp')
            call system('rm NEWBROYDEN')
            call system('rm RADIOC')
            call system('rm RADIOCL')
            call system('rm RADIOL')

            stop 'MAX. NUMBER OF JOBS EXCEEDED'

         endif

      endif

!----------------------------------------------------------------------------------
!PRINT OUT BLOCK END
!----------------------------------------------------------------------------------
 
      DEALLOCATE(DEPTHS)
     
      IF (GAMMAL.GT..0)
     $   PRINT 14, GAMMAL,GAMMAR
   14    FORMAT (/,10X,'SCHARMER-PARAMETERS:',  10X,
     $          '   GAMMAL=',F5.1,/,
     $        40X,'   GAMMAR=',F7.1)
     
      IF (REDUCE .NE. 1.) PRINT 10,REDUCE
   10 FORMAT (/,10X,'CORRECTIONS REDUCED BY FACTOR',F5.2)
     
C***  INHIBIT NEGATIVE POP. NUMBERS
      NEGWARN=0
      DO 19 L=1,ND
      if (.NOT. allocated(ABXYZ_new)) allocate(ABXYZ_new(NATOM))
      ABXYZ_new(1:NATOM)=ABXYZn(1:NATOM,L)
      NEGPOPL=.FALSE.
      DO 16 J=1,N
      IF (POPNUM(L,J) .LE. .0) THEN
            NEGPOPL=.TRUE.
            NEGWARN=NEGWARN+1
            POPNUM(L,J)=0.5*ABS(POP1(L,J))
      ENDIF
   16 CONTINUE
     
      IF (NEGPOPL) THEN
C***     RENORMALIZATION OF THE POPULATION NUMBERS
         DO 15 NA=1,NATOM
         NFIRNA=NFIRST(NA)
         NLANA=NLAST(NA)
         SUM=0.0
         DO 18 J=NFIRNA,NLANA
   18    SUM=SUM+POPNUM(L,J)
         SUM=SUM/ABXYZ_new(NA)
         DO 15 J=NFIRNA,NLANA
         POPNUM(L,J)=POPNUM(L,J)/SUM
   15    CONTINUE
C***     RE-ADJUST THE ELECTRON DENSITY
         RNEL=0.0
         DO 20 J=1,N
         RNEL=RNEL+NCHARG(J)*POPNUM(L,J)
   20    CONTINUE
         RNE(L)=RNEL
      ENDIF
      if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)
   19 CONTINUE
     
      IF (NEGWARN .GT. 0) PRINT 17,NEGWARN
   17 FORMAT (/,10X,'WARNING: NEGATIVE POP. NUMBERS ARE SET TO 0.5 ',
     $ 'TIMES THEIR OLD VALUE:',I4,' TIMES')
C***  PRINTOUT OF ERROR MESSAGES FROM SUBROUTINE GAUNTFF
      IF (LSPOP .EQ. 1) THEN
      IF (NTUP.GT.0)
     $PRINT  *, NTUP,' CALLS OF GAUNTFF BEYOND UPPER TEMPERATURE BOUND'
      IF (NTLOW.GT.0)
     $PRINT *,NTLOW,' CALLS OF GAUNTFF BEYOND LOWER TEMPERATURE BOUND'
      IF (NFUP.GT.0)
     $PRINT *,NFUP,' CALLS OF GAUNTFF BEYOND UPPER FREQUENCY BOUND'
      IF (NFLOW.GT.0)
     $PRINT *,NFLOW,' CALLS OF GAUNTFF BEYOND LOWER FREQUENCY BOUND'
      ENDIF
     
C***  PRINTOUT OF ERROR MESSAGES FROM SUBROUTINE LINPOP
      IF (NEGINTL .GT. 0)
     $  PRINT *,'           WARNING: NEGATIVE LINE INTENSITIES ',
     $      'REPLACED BY THEIR ABSOLUTE VALUE: ',NEGINTL,' TIMES'

      if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)

      RETURN

10002 FORMAT(/,15x,'DEPTHS: ',I2,' TO ',I3,5X,'CORMAX= ',ES9.3,5X,'LOG=',F5.2,
     $         5X,'ELEMENT: ELECTRONS')

      END SUBROUTINE


      SUBROUTINE CALC_CONV_IND(NFIRNA, NLANA, ELEMENT, DEPTHS, DEPTHS_DIM,
     $                         POP, PREV_POP, LEVEL, N, ND, CORMAX, DCM)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: DEPTHS_DIM, NFIRNA, NLANA

      INTEGER, DIMENSION(DEPTHS_DIM), INTENT(IN) :: DEPTHS

      CHARACTER(LEN = *), INTENT(IN) :: ELEMENT

      INTEGER, INTENT(IN) :: N, ND

      REAL*8, DIMENSION(ND, N), INTENT(IN) :: POP, PREV_POP

      CHARACTER (LEN = 10), DIMENSION(N), INTENT(IN) :: LEVEL

      REAL*8,  INTENT(OUT) :: CORMAX

      INTEGER, INTENT(OUT) :: DCM

      CHARACTER (LEN = 10) :: LEVMAX, LEVMIN, LEVMA2, LEVMI2

      REAL*8 :: Q, QMAX, QMA2, QMIN, QMI2

      INTEGER :: LMAX, LMA2, LMIN, LMI2

      REAL*8 :: CORLOG

      INTEGER :: I, J, DEPTH_IND

C***  CALCULATE MINIMUM AND MAXIMUM OF ALL CORRECTION FACTORS
      LMAX = 0
      LMA2 = 0
      LMIN = 0
      LMI2 = 0
      QMAX = 1.0D0
      QMA2 = 1.0D0
      QMIN = 1.0D0
      QMI2 = 1.0D0
      LEVMAX = 'NONE'
      LEVMA2 = 'NONE'
      LEVMIN = 'NONE'
      LEVMI2 = 'NONE'

      DCM = 0
      CORMAX = 0.0D0

      DEPTH: DO I = 1, DEPTHS_DIM

              DEPTH_IND = DEPTHS(I)

         LEVELS: DO J = NFIRNA, NLANA

                 IF (PREV_POP(DEPTH_IND, J) .EQ. 0.0D0) CYCLE LEVELS

                 Q = POP(DEPTH_IND, J) / PREV_POP(DEPTH_IND, J)

                 IF (Q .GT. QMA2) THEN

                    IF (Q .GT. QMAX) THEN

                       QMA2 = QMAX
                       QMAX = Q
                       LEVMA2 = LEVMAX
                       LEVMAX = LEVEL(J)
                       LMA2 = LMAX
                       LMAX = DEPTH_IND

                    ELSE

                       QMA2 = Q
                       LEVMA2 = LEVEL(J)
                       LMA2 = DEPTH_IND

                    ENDIF

                 ENDIF

                 IF (Q .LT. QMI2) THEN

                    IF (Q .LT. QMIN) THEN

                       QMI2 = QMIN
                       QMIN = Q
                       LEVMI2 = LEVMIN
                       LEVMIN = LEVEL(J)
                       LMI2 = LMIN
                       LMIN = DEPTH_IND

                    ELSE

                       QMI2 = Q
                       LEVMI2 = LEVEL(J)
                       LMI2 = DEPTH_IND

                    ENDIF

                 ENDIF

         ENDDO LEVELS

      ENDDO DEPTH

      CORMAX = DMAX1(QMAX - 1.0D0, 1.0D0 - QMIN)

      IF (ABS(CORMAX - QMAX + 1) .LT. 1D-15) DCM = LMAX

      IF (ABS(CORMAX + QMIN - 1) .LT. 1D-15) DCM = LMIN

      IF (CORMAX .GT. 0.0D0) CORLOG = LOG10(CORMAX)

      IF (NFIRNA .NE. NLANA) THEN

         WRITE(*, '(/,A)') '---------------------------------------------------------
     $--------------------------------------------------------'

         WRITE(*, 9) QMAX, LEVMAX, LMAX, QMA2, LEVMA2, LMA2, QMIN, LEVMIN, LMIN, QMI2, LEVMI2, LMI2,
     $            DEPTHS(1), DEPTHS(DEPTHS_DIM), CORMAX, CORLOG, ELEMENT

      ENDIF

    9 FORMAT(/,15X, 'QMAX1: ', ES9.3,'  (', A10,'  L=', I3,')',
     $         4X,  'QMAX2: ', ES9.3,'  (', A10,'  L=', I3,')',
     $       /,15X, 'QMIN1: ', ES9.3,'  (', A10,'  L=', I3,')',
     $         4X,  'QMIN2: ', ES9.3,'  (', A10,'  L=', I3,')', /,
     $       /,15X, 'DEPTHS: ', I2, ' TO ', I3,
     $         5X,  'CORMAX =',ES9.3,5X,'LOG(CORMAX) = ',F5.2,
     $         5X,  'ELEMENT: ', A)

      END SUBROUTINE

      END MODULE
