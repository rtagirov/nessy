      MODULE MOD_COMA

      CONTAINS

      SUBROUTINE COMA(CRATE,RRATE,RATCO,DM,N,NRANK,V1,ABXYZ,
     $                ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,EION,WEIGHT,ALTESUM,
     $                XLAMBDA,FWEIGHT,XJC,NF,L,XJL,ND,XJLAPP,SLOLD,LASTIND,INDLOW,
     $                INDNUP,NOM,NATOM,KODAT,NFIRST,NLAST,PHI,PWEIGHT,DELTAX,XMAX,
     $                NFL,OPAC,SCNEW,DOPA,DETA,OPAL,SLNEW,DOPAL,DETAL,SIGMAKI,
     $                ETAC,NFEDGE,EXPFAC,NOTEMP,NODM,
     $                WCHARM,EN,RSTAR,SCOLD,VDOP,COCO,KEYCOL,
     $                POPHIIL,POPHML,POPHIL,LOE,ITNEL,LEVEL,JOBNUM,IRESTA)

      use MOD_DERIV
      use MOD_COOPFRQ
      use MOD_SETXJL
      use MOD_RADNET
      use MOD_COLLI
      use MOD_DLIOP
      use MOD_DCOOP
      USE FILE_OPERATIONS
      USE COMMON_BLOCK

C*******************************************************************************
C***  THIS ROUTINE SETS UP THE RATE COEFFICIENT MATRIX RATCO
C***  AND ITS VECTOR DERIVATIVE DM (LINEARIZED MATRIX)
C***  IF NODM IS TRUE, THE ROUTINE DOES NOT SET UP THE DM MATRIX
C*** CALLED BY LINPOP
C*** ENE:
C*** RNE: RELATIVE ELECTRON DENSITY
C***	ENTOTL: PARTICLE DENSITY AT DEPTH POINT L
C*******************************************************************************
C fflow ------------------ Flow chart ------------------------
C 1     | -> STEAL, line 128 (steal.for, line 0)
C :     :\
C 2     | | -> LINPOP, line 212 (linpop.for, line 0)
C :     : :\
C 3     | | | -> COMA, line 208 (coma.for, line 0)
C 4     | | | | -> COOPFRQ, line 65 (coopfrq.for, line 1)
C 5     | | | | | -> GFFLOG, line 71 (<unknown>)
C 4     | | | | -> SETXJL, line 83 (setxjl.for, line 1)
C 5     | | | | | -> LIOP, line 30 (liop.for, line 1)
C 5     | | | | | -> LIPO, line 43 (lipo.for, line 1)
C 5     | | | | | -> XRUDI, line 45 (xrudi.for, line 1)
C 4     | | | | -> COLLI, line 92 (colli.for, line 1)
C 4     | | | | -> RADNET, line 98 (radnet.for, line 1)
C 5     | | | | | -> XRUDI, line 49 (xrudi.for, line 1)
C 4     | | | | -> DCOOP, line 152 (dcoop.for, line 0)
C 5     | | | | | -> GFFLOG, line 97 (<unknown>)
C 5     | | | | | -> GFFLOG, line 147 (<unknown>)
C 4     | | | | -> DLIOP, line 156 (dliop.for, line 1)
C 4     | | | | -> DERIV, line 163 (deriv.for, line 1)
C :
C :
C 1     | -> STEAL, line 217 (steal.for, line 0)
C :      :
C 2     | | -> LINPOP, line 212 (linpop.for, line 0)
C :	: :\
C 3     | | | -> COMA, line 208 (coma.for, line 0)
C       : : :

      IMPLICIT NONE

      integer,intent(in) ::  LASTIND, NATOM, NRANK,  N
      integer,intent(in) ::  KODAT(*), NOM(N)
      integer,dimension(*):: INDNUP, INDLOW
      integer,intent(in) ::  ND, NF,  NFL
      integer, dimension(N)  ::  NCHARG
      integer, dimension(NATOM) ::  NFIRST,NLAST
      real*8, dimension(NRANK),intent(inout) ::  V1
      real*8             :: TL, RSTAR
      real*8             :: COCO(N,N,4), VDOP
      real*8             :: EINST(N,N), ENE,DELTAX
      real*8             :: ALTESUM(4,*), XMAX

      real*8, dimension(lastind) :: opal

      integer, dimension(*)      :: NFEDGE
      real*8,  dimension(*)      :: PHI
      real*8,  dimension(*)      :: ELEVEL,EION,WEIGHT
      real*8,  dimension(*)      :: XLAMBDA,PWEIGHT,DETAL
      real*8,  dimension(*)      :: EXPFAC,DOPAL,FWEIGHT
      real*8,  dimension(NF, N)  :: SIGMAKI
      real*8,  dimension(NF)     :: ETAC,DETA
      real*8,  dimension(NF)     :: OPAC,DOPA
      real*8,  allocatable       :: ABXYZ(:)
      real*8,  dimension(ND, NF) :: XJC, WCHARM
      real*8,  dimension(NF, ND) :: SCOLD

      real*8, dimension(N) :: ENLTE
      real*8, dimension(N + 1) :: EN

      REAL*8, DIMENSION(LASTIND), INTENT(OUT) :: XJLAPP

      REAL*8, DIMENSION(N,  N),        INTENT(INOUT) :: CRATE, RRATE
      REAL*8, DIMENSION(NRANK, NRANK), INTENT(INOUT) :: RATCO

      INTEGER, INTENT(IN) :: L

      REAL*8, DIMENSION(LASTIND), INTENT(IN) :: SLOLD

!     the Local approximate lambda-Operator Element corresponding to depth L for all lines
      REAL*8, DIMENSION(LASTIND) :: LOE

      REAL*8, DIMENSION(LASTIND) :: XJL
      REAL*8, DIMENSION(LASTIND) :: AccFact
      REAL*8, DIMENSION(LASTIND) :: SLNEW

      LOGICAL NOTEMP, NODM, DEPTH_PR_COND, JOBNUM_PR_COND, ITNEL_COND, PR_COND

      integer :: NFIRNA, NA, I, J, K, NPLUS1, NPLUS2, NUP, IND, LOW
      integer :: NLANA
      REAL*8  :: ENTOTL, RNEL, DLOWUP

      real*8, dimension(NF) :: SCNEW, XJCAPP

      REAL*8, DIMENSION(NRANK, NRANK), INTENT(INOUT) :: DM

      REAL*8 :: CP

      integer, parameter :: IONE = 1

      CHARACTER*4 KEYCOL(N, N)

      REAL*8, INTENT(IN) :: POPHIIL, POPHML, POPHIL

      CHARACTER*10, DIMENSION(N), INTENT(IN) :: LEVEL

      INTEGER, INTENT(IN) :: ITNEL, JOBNUM, IRESTA

      REAL*8, ALLOCATABLE, DIMENSION(:, :) :: TNETRATE, RNETRATE, CNETRATE

      REAL*8 ::                               TRUPDOWN, TRDOWNUP

      INTEGER :: N1, N2

      NPLUS1 = N + 1
      NPLUS2 = N + 2

C***  CALCULATE CONTINUUM OPACITIES AND EMISSIVITIES FROM CURRENT POPNUMBERS.
C***  ONLY TRUE OPACITIES ARE ACCOUNTED FOR.

      RNEL = EN(NPLUS1)

      ENTOTL = ENE / RNEL

      CALL COOPFRQ(NF,OPAC,ETAC,XLAMBDA,EXPFAC,SIGMAKI,N,NCHARG,
     $             WEIGHT,ELEVEL,EION,NFEDGE,EN,NOM,RSTAR,ENTOTL,RNEL,TL)
C***  CALCULATE NEW SOURCE FUNCTION AND SCHARMER'S RADIATION FIELD
C***  LOOP OVER ALL CONT. FRQUENCIES  ----------------------------------

      DO K = 1, NF

!     PREVENT DIVIDE CHECK ERRORS

         IF (OPAC(K) .GT. 0.0d0) THEN

!       OPAC can not be < 0 because of laser security check in ccore

            SCNEW(K) = ETAC(K) / OPAC(K)

         ELSE

            SCNEW(K) = 0.0d0

         ENDIF

!         write(*, '(A,2x,2(i4,2x),4(e15.7,2x))'), 'coma xjc check:',
!     $        k, l, xjc(k, l), wcharm(l, k), scnew(k), scold(l, k)

         XJCAPP(K) = XJC(L, K) + WCHARM(L, K) * (SCNEW(K) - SCOLD(K, L))

      ENDDO

C***  ENDLOOP  ---------------------------------------------------------

      ITNEL_COND = ITNEL .GT. 1

      DEPTH_PR_COND = L .LE. 61 .AND. L .GE. 52

      JOBNUM_PR_COND = .TRUE.

      PR_COND = ITNEL_COND .AND. DEPTH_PR_COND .AND. JOBNUM_PR_COND

! Rinat Tagirov: printing the rates
!**********************************************************************************************************************

      IF (PR_COND) THEN

      N1 = NFIRST(1)
      N2 = NLAST(1)

      ALLOCATE(TNETRATE(N2 - N1 + 1, N2 - N1 + 1))
      ALLOCATE(CNETRATE(N2 - N1 + 1, N2 - N1 + 1))
      ALLOCATE(RNETRATE(N2 - N1 + 1, N2 - N1 + 1))
!      ALLOCATE(RBRACKET(N2 - N1 + 1, N2 - N1 + 1))

      TNETRATE(:, :) = 0.0D0
      CNETRATE(:, :) = 0.0D0
      RNETRATE(:, :) = 0.0D0
!      RBRACKET(:, :) = 0.0D0

      DO I = N1, N2

         DO J = N1, N2

           IF (J .GT. I) THEN

               TRUPDOWN = EN(J) * (CRATE(J, I) + RRATE(J, I))
               TRDOWNUP = EN(I) * (CRATE(I, J) + RRATE(I, J))

               TNETRATE(I, J) = TRUPDOWN - TRDOWNUP

               CNETRATE(I, J) = EN(J) * CRATE(J, I) - EN(I) * CRATE(I, J)

               RNETRATE(I, J) = EN(J) * RRATE(J, I) - EN(I) * RRATE(I, J)

!               RBRACKET(I, J) = (EN(J) * RRATE(J, I) - EN(I) * RRATE(I, J))
!     $                          / (EN(J) * EINST(J, I))

           ENDIF

         ENDDO

      ENDDO

!      NRRM_FILE = TRIM(ADJUSTL(NLTE_DIR_3//NRRM_FILE_NAME))
!      NCRM_FILE = TRIM(ADJUSTL(NLTE_DIR_3//NCRM_FILE_NAME))
!      NTRM_FILE = TRIM(ADJUSTL(NLTE_DIR_3//NTRM_FILE_NAME))

!      CALL PRINT_RATE_MATRIX(648, NTRM_FILE, TNETRATE(N1 : N2, N1 : N2),
!     $                       N1, N2, NODM, JOBNUM, L, ITNEL)

!      CALL PRINT_RATE_MATRIX(295, NCRM_FILE, CNETRATE(N1 : N2, N1 : N2),
!     $                       N1, N2, NODM, JOBNUM, L, ITNEL)

!      CALL PRINT_RATE_MATRIX(257, NRRM_FILE, RNETRATE(N1 : N2, N1 : N2),
!     $                       N1, N2, NODM, JOBNUM, L, ITNEL)

      DEALLOCATE(TNETRATE)
      DEALLOCATE(CNETRATE)
      DEALLOCATE(RNETRATE)

      ENDIF

!**********************************************************************************************************************

 
C***  SETUP THE COLLISIONAL AND RADIATIVE RATE COEFFICIENTS
CMH - new: POPHIIL: population numbers of HII at depthpoint L
CMH - new: needed to calculate new collision cross sections for Hminus

      CALL COLLI(N,ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,CRATE,
     $           EION,COCO,KEYCOL,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     $           POPHIIL, POPHML, POPHIL, LEVEL, JOBNUM, L)

!     CALCULATE LINE RADIATION FIELD WITH APPROXIMATE LAMBDA OPERATOR TERMS

      CALL SETXJL(LASTIND, INDLOW, INDNUP, SLNEW(1 : LASTIND),
     $            SLOLD(1 : LASTIND), OPAL, XJLAPP(1 : LASTIND),
     $            NF, XLAMBDA, SCNEW, OPAC,
     $            NFL, PHI, PWEIGHT, N, EINST, ELEVEL, EN, WEIGHT,
     $            ND, XJL(1 : LASTIND), ENTOTL, RSTAR, VDOP, DELTAX, XMAX, L,
     $            LOE(1 : LASTIND), AccFact(1 : LASTIND),
     $            NODM, LEVEL, NFIRST, NLAST, NATOM, ENLTE, ITNEL)

!     RADIATIVE RATES ARE CALCULATED WITH THE MODIFIED RADIATION FIELD

      CALL RADNET(N, ENLTE, TL, WEIGHT, NCHARG ,EION, ELEVEL, EINST,
     $            SLNEW(1 : LASTIND), EN, NOM, RRATE, XLAMBDA, FWEIGHT,
     $            XJCAPP(1 : NF), NF, XJLAPP(1 : LASTIND), SIGMAKI, LASTIND,
     $            LEVEL, L, JOBNUM, ITNEL)

!     ADD RADIATIVE AND COLLISIONAL TERMS INTO RATE COEFFICIENT MATRIX RATCO

      RATCO(1 : N, 1 : N) = 0.0D0

      RATCO(1 : N, 1 : N) = - RRATE(1 : N, 1 : N) - CRATE(1 : N, 1 : N)

C***  DIAGONAL ELEMENTS: -SUM OF THE ROW (I.E. OVER COLUMN INDEX)

      FORALL(I = 1 : N)

        RATCO(I, I) = -SUM(RATCO(I, 1 : N))

      END FORALL

C***  COLUMN NLAST(NA): NUMBER CONSERVATION FOR EACH ELEMENT (NA)
C***  REMARK: TOTAL NUMBER CONSERVATION IS IMPLICITLY ENSURED

      DO NA = 1, NATOM

         NFIRNA = NFIRST(NA)
         NLANA = NLAST(NA)

         RATCO(NFIRNA : NLANA, NLANA) = 1.0D0

      ENDDO

C***  COLUMN N+1 : CHARGE CONSERVATION
      DO 4 I=1,N
    4 RATCO(I,NPLUS1)=NCHARG(I)
      RATCO(NPLUS1,NPLUS1)=-1.
C***  ROW N+1 : ZERO
      DO 5 J=1,N
    5 RATCO(NPLUS1,J)=.0d0
     
C***  IF TEMPERATURE CORRECTIONS ARE ALLOWED (ENERGY EQUATION):
      IF (.NOT. NOTEMP) THEN
C***  ADDITIONAL ROW AND COLUMN = 0
      DO 11 J=1,NPLUS1
      RATCO(J,NPLUS2)=.0d0
      RATCO(NPLUS2,J)=.0d0
   11 CONTINUE
C***  LOWER RIGHT CORNER ELEMENT:
      RATCO(NPLUS2,NPLUS2)= 1.d0
      ENDIF

      DM(1 : NRANK, 1 : NRANK) = 0.0D0

C***  IF NODM IS TRUE, THE DM MATRIX DOES NOT HAVE TO BE SET UP
      IF (NODM) GOTO 101

C***  DERIVATIVE MATRIX DM
C***  FIRST TERMS : THE ORIGINAL MATRIX RATCO

      DM(1 : NRANK, 1 : NRANK) = RATCO(1 : NRANK, 1 : NRANK)

      DO 10 I = 1, NPLUS1

C***  CONSTRUCT DERIVATIVE VECTORS DOPA, DETA WITH RESPECT TO EN(I)

      CALL DCOOP(I,DOPA,DETA,XLAMBDA,NF,TL,RNEL,ENTOTL,EN,RSTAR,
     $           WCHARM,ND,L,NFEDGE,EXPFAC,N,NCHARG,WEIGHT,
     $           ELEVEL,EION,NOM,EINST,SIGMAKI)

C***  CONSTRUCT DERIVATIVE VECTORS DOPAL, DETAL (LINES) WITH RESPECT TO EN(I)
      CALL DLIOP(I,ENTOTL,DOPAL,DETAL,VDOP,RSTAR,N,
     $           EINST,WEIGHT,ELEVEL,LASTIND,INDLOW,INDNUP)
     
      IND = 0

      DO 8 NUP = 2, N

      DO 7 LOW = 1, NUP - 1

C***  COMPUTE ELEMENT (LOW,NUP) OF MATRIX DM (DERIVATIVE WITH RESPECT TO EN(I))
      CALL DERIV(DLOWUP,I,NUP,LOW,IND,
     $           NPLUS1,EN,CRATE,RRATE,EXPFAC,NFEDGE,
     $           WCHARM,ND,L,TL,ENLTE,PHI,PWEIGHT,NFL,DELTAX,XMAX,
     $           DETAL,DOPAL,SLNEW(1 : LASTIND),OPAL,XJLAPP,XJCAPP,
     $           FWEIGHT,DOPA,DETA,OPAC,SCNEW,XLAMBDA,NF,
     $           N,NCHARG,WEIGHT,ELEVEL,NOM,EINST,SIGMAKI,LASTIND,
     $           XJL(1 : LASTIND), AccFact(1 : LASTIND), SLOLD(1 : LASTIND))

      DM(I, NUP) = DM(I, NUP) + DLOWUP

C***  NOTE THAT DUPLOW = - DLOWUP
      DM(I, LOW) = DM(I, LOW) - DLOWUP

    7 CONTINUE

    8 CONTINUE

   10 CONTINUE

C***  COLUMNS NLAST(NA)  (I.E. COLUMNS CONTAINING THE EQUATIONS OF NUMBER
C***  CONSERVATION FOR ELEMENT NA)  ARE NOT CHANGED

      DO NA = 1, NATOM

         NLANA = NLAST(NA)

         DM(1 : NPLUS1, NLANA) = RATCO(1 : NPLUS1, NLANA)

      ENDDO

C***  CONTINUE HERE IF CALCULATION OF DM NOT NEEDED (NODM = TRUE)
  101 CONTINUE

C***  V1 = RIGHT-HAND SIDE VECTOR (INHOMOGENEITY)

      V1(1 : NPLUS1) = 0.0D0

C***          NLAST(NA)-TH ELEMENT = ABXYZ(NA)  (NUMBER CONSERVATION)

      DO NA = 1, NATOM

        NLANA = NLAST(NA)

        V1(NLANA) = ABXYZ(NA)

      ENDDO

C***  ADDITIONAL ELEMENT (N+2) FOR THE ENERGY EQUATION
      IF (.NOT. NOTEMP) THEN

         V1(NPLUS2) = TL

      ENDIF

      RETURN

10000 FORMAT(A,1x,I3,1x,I3,1x,I3,1x,E15.7,1x,E15.7)
!11000 FORMAT(A,1x,I3,1x,I3,1x,I3,1x,E15.7,1x,E15.7)

      END SUBROUTINE


      SUBROUTINE PRINT_RATE_MATRIX(FileUnit,
     $                             FileName,
     $                             NetRateMatrix,
     $                             NStart,
     $                             NStop,
     $                             NODM,
     $                             JOBNUM,
     $                             DepthInd,
     $                             IterNum)

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      CHARACTER (LEN = *), INTENT(IN) :: FileName

      INTEGER, INTENT(IN) :: FileUnit

      INTEGER, INTENT(IN) :: NStart, NStop ! The indices constraining the square (corresponding to any certain chemical element) of the original rate matrix

      INTEGER, INTENT(IN) :: JOBNUM, DepthInd, IterNum

      REAL*8, DIMENSION(:,:), INTENT(IN) :: NetRateMatrix

      LOGICAL, INTENT(IN) :: NODM

      INTEGER :: I, J

      CALL OPEN_TO_APPEND(FileUnit, FileName)

      WRITE(FileUnit, '(A,I3,A,I3,A,I3,A,L,/)') 'JOBNUM = ', JOBNUM,
     $                                          ', L = ', DepthInd,
     $                                          ', ITNEL = ', IterNum - 1,
     $                                          ', NODM = ', NODM


      WRITE(FileUnit, '(A,/)') '----------------------------------------
     $------------------------------------------------------------------
     $-----------------------------------------------------------------
     $-----------'

      DO I = NStart, NStop

         WRITE(FileUnit, '($, 13x, I2)') I

      ENDDO
      
      WRITE(FileUnit, '(/)')

      DO I = NStart, NStop

         WRITE(FileUnit, '($, I2)') I

         DO J = NStart, NStop

            WRITE(FileUnit, '($, E15.7)') NetRateMatrix(I, J)

         ENDDO

         WRITE(FileUnit, '(/)')

      ENDDO

      WRITE(FileUnit, '(A,/)') '----------------------------------------
     $------------------------------------------------------------------
     $-----------------------------------------------------------------
     $-----------'

      IF (FileName .EQ. NLTE_DIR_3//NTRM_FILE_NAME) THEN

      DO J = NStart, NStop

         IF (J .EQ. NStart) THEN
            
         WRITE(FileUnit, '(2x, $, E15.7)') SUM(NetRateMatrix(NStart,
     $                                         NStart + 1 : NStop))

         ELSEIF (J .EQ. NStop) THEN

         WRITE(FileUnit, '($, E15.7, /)') SUM(NetRateMatrix(NStart : 
     $                                        NStop - 1, NStop))

         ELSE

         WRITE(FileUnit,'($, E15.7)')
     $   SUM(NetRateMatrix(NStart : J - 1, J))
     $ - SUM(NetRateMatrix(J, J + 1 : NStop))

         ENDIF

      ENDDO

      WRITE(FileUnit, '(/,A,/)') '--------------------------------------
     $------------------------------------------------------------------
     $------------------------------------------------------------------
     $------------'

      ENDIF

!      IF (DepthInd .EQ. DPN) CLOSE(FileUnit)

      END SUBROUTINE


!      SUBROUTINE PRINT_BROYDEN_STEP(FileUnit,
!     $                              FileName,
!     $                              FileNameLength,
!     $                              TNETRATE,
!     $                              CNETRATE,
!     $                              RNETRATE,
!     $                              RBRACKET,
!     $                              LowerLevel,
!     $                              UpperLevel,
!     $                              JOBNUM,
!     $                              DepthInd,
!     $                              BroydenIterNum,
!     $                              NODM)
!
!      USE FILE_OPERATIONS
!      USE COMMON_BLOCK
!
!      IMPLICIT NONE
!
!!      COMMON /DepthPointsNum/ DPN
!
!!      INTEGER :: DPN
!
!      INTEGER, INTENT(IN) :: FileNameLength
!
!      CHARACTER (LEN = FileNameLength), INTENT(IN) :: FileName
!
!      INTEGER, INTENT(IN) :: FileUnit
!
!      CHARACTER (LEN = 10), INTENT(IN) :: LowerLevel,
!     $                                    UpperLevel
!
!      REAL*8, INTENT(IN) :: TNETRATE,
!     $                      CNETRATE,
!     $                      RNETRATE,
!     $                      RBRACKET
!
!      INTEGER, INTENT(IN) :: JOBNUM,
!     $                       BroydenIterNum
!
!      INTEGER, INTENT(IN) :: DepthInd
!
!      LOGICAL, INTENT(IN) :: NODM
!
!      CALL OPEN_TO_APPEND(FileUnit, FileName)
!
!      IF (LAMBDA_ITER .EQ. 1 .AND. DepthInd .EQ. 1 .AND. BroydenIterNum .EQ. 2) THEN
!
!         WRITE(FileUnit, '(2x,A,10x,A,10x,A,4x,A,4x,A,
!     $                     7x,A,4x,A,7x,A,12x,A,12x,A,12x,A,/)') 'LowLev',
!     $                                                           'UpLev',
!     $                                                           'LI',
!     $                                                           'JOBNUM',
!     $                                                           'L',
!     $                                                           'BI',
!     $                                                           'NODM',
!     $                                                           'TNETRATE',
!     $                                                           'CNETRATE',
!     $                                                           'RNETRATE',
!     $                                                           'RBRACKET'
!
!      ENDIF
!
!      WRITE(FileUnit, '(A10,2x,A2,2x,A10,5x,I3,5x,I3,5x,I3,5x,I3,5x,L,5x,
!     $                  E15.7,5x,E15.7,5x,E15.7,5x,E15.7)') LowerLevel, '->', UpperLevel,
!     $                                                      LAMBDA_ITER,
!     $                                                      JOBNUM,
!     $                                                      DepthInd,
!     $                                                      BroydenIterNum,
!     $                                                      NODM,
!     $                                                      TNETRATE,
!     $                                                      CNETRATE,
!     $                                                      RNETRATE,
!     $                                                      RBRACKET
!
!      IF (DepthInd .EQ. DPN) CLOSE(FileUnit)
!
!      END SUBROUTINE


      END MODULE
