      MODULE MOD_SETXJL

      CONTAINS

      SUBROUTINE SETXJL(LASTIND, INDLOW, INDNUP, SLNEW, SLOLD, OPAL, XJLAPP,
     $                  NF, XLAMBDA, SCNEW, OPAC,
     $                  NFL, PHI, PWEIGHT, N, EINST, ELEVEL, EN, WEIGHT,
     $                  ND, XJL, ENTOTL, RSTAR, VDOP, DELTAX, XMAX, L, LOE,
     $                  AccFact, NODM, LEVEL, NFIRST, NLAST, NATOM, ENLTE, ITNEL)

!     CALCULATE LINE RADIATION FIELD WITH APPROXIMATE LAMBDA OPERATOR TERMS

      USE MOD_LIPO
      USE MOD_LIOP
      USE MOD_XRUDI
      USE CONSTANTS
      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)

      PARAMETER (one = 1.d0, two = 2.d0)

      integer, intent(in) :: lastind, N

      DIMENSION EINST(N, N), ELEVEL(N), EN(N + 1), WEIGHT(N)

      DIMENSION XLAMBDA(NF),SCNEW(NF),OPAC(NF)
      DIMENSION PHI(NFL),PWEIGHT(NFL)
      DIMENSION INDNUP(LASTIND),INDLOW(LASTIND)

      INTEGER, INTENT(IN) :: NATOM

      INTEGER, DIMENSION(NATOM), INTENT(IN) :: NFIRST, NLAST

      REAL*8 :: etalind, opalind

      REAL*8, DIMENSION(LASTIND), INTENT(OUT) :: opal

      REAL*8, DIMENSION(LASTIND), INTENT(IN) ::  XJL

      REAL*8, DIMENSION(LASTIND), INTENT(OUT) :: XJLAPP

!     DEPTH INDEX
      INTEGER, INTENT(IN) :: L

!     the Local approximate lambda-Operator Element corresponding to depth L for all lines
      REAL*8, DIMENSION(LASTIND), INTENT(IN) ::  LOE

      REAL*8, DIMENSION(LASTIND), INTENT(IN) ::  SLOLD
      REAL*8, DIMENSION(LASTIND), INTENT(OUT) :: SLNEW, AccFact

      CHARACTER*10, DIMENSION(N), INTENT(IN) :: LEVEL

      REAL*8, DIMENSION(N), INTENT(IN) :: ENLTE

      LOGICAL, INTENT(IN) :: NODM

      INTEGER, INTENT(IN) :: ITNEL

      REAL*8 :: UpperLevelDep, LowerLevelDep

      REAL*8 :: AccFactCutOff

      LOGICAL :: PRINT_COND, DEPTH_ACC_COND

      opal(1 : lastind) = 0.0d0

      DEPTH_ACC_COND = L .NE. ND

      IF (ITNEL .EQ. 1) WRITE(*, '(/,2x,A,5x,A,5x,A,3x,A,7x,A,11x,A,12x,A,13x,A,19x,A,17x,A,14x,A,11x,A,14x,A,13x,A,/)'),
     $                           'LI', 'L', 'BI', 'IND', 'LL', 'UL', 'SLO', 'SLN', 'DE', 'JO', 'DJ', 'DJ/JO', 'ULD', 'LLD'

      XJLAPP(1 : LASTIND) = XJL(1 : LASTIND)

      SLNEW(1 : LASTIND) = 0.0D0

      AccFact(1 : LASTIND) = 0.0D0

!     LOOP OVER ALL LINES
      DO IND = 1, LASTIND

      AccFactCutOff = 10.0D+0

      PRINT_COND = IND .EQ. 1

      LOW = INDLOW(IND)
      NUP = INDNUP(IND)

      IF (LOW .GE. NFIRST(2) .AND. NUP .LE. NLAST(2)) AccFactCutOff = 0.0D0 ! Helium is not accelerated

!     FOR RUDIMENTAL LINES, ZERO CORE IS ASSUMED
      IF (EINST(LOW, NUP) .EQ. -2.0D0) CYCLE

      UpperLevelDep = EN(NUP) / ENLTE(NUP)
      LowerLevelDep = EN(LOW) / ENLTE(LOW)

      XLAM = 1.0D8 / (ELEVEL(NUP) - ELEVEL(LOW))

      if (ind .eq. 1) print*, 'setxjl:', en(low), en(nup)

!      stop

      CALL LIOP_SBE(EINST(NUP, LOW), WEIGHT(LOW), WEIGHT(NUP), LOW, NUP,
     $              XLAM, ENTOTL, EN(1 : N), RSTAR, opalind, etalind, VDOP, N)

C***  LASER SECURITY

      IF (opalind .LE. 0.0D0) THEN

         IF (PRINT_COND) WRITE(*, 10000), LAMBDA_ITER, L, ITNEL, IND,
     $                                    LEVEL(LOW), LEVEL(NUP),
     $                                    UpperLevelDep, LowerLevelDep,
     $                                    '<---------------- OPAL  IS NEGATIVE'

         CYCLE

      ENDIF

      opal(ind) = opalind

      slnew(ind) = etalind / opalind

!      if (ind .eq. 1) then

!         write(*, '(A,2x,i4,4(2x,e15.7))'), 'setxjl:', ind, slnew(ind), etalind, opalind, opal(ind)

!         stop

!      endif

      IF (SLOLD(IND) .EQ. 0.0D0) THEN

         IF (PRINT_COND) WRITE(*, 10000), LAMBDA_ITER, L, ITNEL, IND, 
     $                                    LEVEL(LOW), LEVEL(NUP),
     $                                    UpperLevelDep, LowerLevelDep,
     $                                    '<---------------- SLOLD IS ZERO'

         CALL SET_ACC_FACT(DEPTH_ACC_COND, AccFactCutOff, LOE(IND), AccFact(IND))

         CYCLE

      ENDIF

      CALL SET_ACC_FACT(DEPTH_ACC_COND, AccFactCutOff, LOE(IND), AccFact(IND))

      DELTASL = SLNEW(IND) - SLOLD(IND)

      XJLAPP(IND) = XJLAPP(IND) + AccFact(IND) * DELTASL

      IF (PRINT_COND) WRITE(*, 10001), LAMBDA_ITER, L, ITNEL, IND,
     $                                 LEVEL(LOW), LEVEL(NUP),
     $                                 SLOLD(IND), SLNEW(IND),
     $                                 LOE(IND),
     $                                 XJL(IND), AccFact(IND) * DELTASL,
     $                                 AccFact(IND) * DELTASL / XJL(IND),
     $                                 UpperLevelDep, LowerLevelDep


      IF (.NOT. NODM .AND. PRINT_COND) WRITE(*, '(3x,A)'), 'NEWTON'

      IF (NODM .AND. PRINT_COND)       WRITE(*, '(3x,A)'), 'BROYDEN'

      ENDDO

10000 FORMAT(I5,3x,    ! LI
     $       I3,3x,    ! L
     $       I3,3x,    ! BI
     $       I3,3x,    ! IND
     $       A10,3x,   ! LL
     $       A10,113x, ! UL
     $       E15.7,1x, ! ULD
     $       E15.7,3x, ! LLD
     $       A)

10001 FORMAT(I5,3x,      ! LI
     $       I3,3x,      ! L
     $       I3,3x,      ! BI
     $       I3,3x,      ! IND
     $       A10,3x,     ! LL
     $       A10,1x,     ! UL
     $       E15.7,1x,   ! SLO
     $       E15.7,1x,   ! SLN
     $       E23.15,1x,  ! DE
     $       E15.7,1x,   ! JO
     $       E15.7,1x,   ! DJ
     $       E15.7,1x,   ! R
     $       E15.7,1x,   ! ULD
     $       E15.7,3x,$) ! LLD

      END SUBROUTINE


      SUBROUTINE SET_ACC_FACT(AccCond, AccFactCutOff, LOE, AccFact)

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: AccCond

      REAL*8, INTENT(IN) ::  AccFactCutOff, LOE

      REAL*8, INTENT(OUT) :: AccFact

      AccFact = 0.0D0

      IF (AccCond) AccFact = LOE

      IF (AccFact .GT. AccFactCutOff) AccFact = AccFactCutOff

      END SUBROUTINE


      END MODULE
