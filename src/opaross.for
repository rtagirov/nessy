      module MOD_OPAROSS

      contains

      SUBROUTINE OPAROSS(OPARL,EN,TL,RNEL,ENTOTL,RSTAR,N,
     $                   LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $                   ALPHA,SEXPO,AGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
     $                   WAVARR,SIGARR)

C***  CALLED BY GREYM
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT DEPTH POINT L
c***  AND FREQUENCY POINT K=1,NF
C***  FOR GIVEN POPNUMBERS

      use MOD_COOP

      IMPLICIT REAL*8(A - H, O - Z)

      real*8, dimension(NF, N) :: DUMMY2

      DIMENSION EINST(N, N)
      DIMENSION EN(N)
      DIMENSION NCHARG(N),WEIGHT(N),ELEVEL(N),EION(N),ALPHA(N),SEXPO(N)
      DIMENSION NOM(N)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)
      CHARACTER*10 LEVEL(N),MAINPRO,MAINLEV

      character*8, dimension(N) :: AGAUNT

      real*8,dimension(1) :: TL_,RNEL_,ENTOTL_,OPA_,ETA_,THOMSON_
      INTEGER,DIMENSION(1) :: IWARN_
      character*10,dimension(1) :: MAINPRO_,MAINLEV_
C******************************************************************************
CMH  CHANGES BY MARGIT HABERREITER, 20 MAY, 2002
CMH  CALLED BY GREYM
	  DIMENSION WAVARR(N, NF), SIGARR(N, NF)
C******************************************************************************
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
      DATA C1,C2 / 1.4388E8, 3.9724E+8 /
C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI  (ERG/CM**2/S/STERAD/KELVIN**4)
      DATA STEBOL / 1.8046E-5 /
     
C***  FREQUENCY INTEGRATION ***
      SUM=.0

      DO 1 K = 1, NF

      XLAM =     XLAMBDA(K)

      TL_      = TL
      RNEL_    = RNEL
      ENTOTL_  = ENTOTL

      CALL COOP(XLAM,1,TL_,RNEL_,EN,ENTOTL_,RSTAR,OPA_,ETA_,
     $          THOMSON_,IWARN_,MAINPRO_,MAINLEV_,NOM,N,LEVEL,NCHARG,WEIGHT,
     $          ELEVEL,EION,EINST,ALPHA,SEXPO,AGAUNT,0,DUMMY2,
     $          WAVARR(1 : N, 1 : NF),SIGARR(1 : N, 1 : NF),NF)

      TL       = TL_(1)
      RNEL     = RNEL_(1)
      ENTOTL   = ENTOTL_(1)
      OPA      = OPA_(1)
      ETA      = ETA_(1)
      THOMSON  = THOMSON_(1)
      IWARN    = IWARN_(1)
      MAINPRO  = MAINPRO_(1)
      MAINLEV  = MAINLEV_(1)

C***  DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
C***  FOR SMALL VALUES OF TL (I.E. FOR LARGE VALUES OF RMAX IN CASE OF A
C***  SPHERICAL TEMPERATURE STRUCTURE): PREVENTION OF AR004 (BAD SCALAR
C***  ARGUMENT TO ARLIB MATH ROUTINE)
      HNUEKT=C1/XLAM/TL
      IF (HNUEKT .GT. 5000.) THEN
          DBDT=C1*C2/XLAM/XLAM/XLAM/XLAM/TL/TL
      ELSE
          EXFAC=EXP(HNUEKT)
          DBDT=C1*C2/XLAM/XLAM/XLAM/XLAM/TL/TL/(EXFAC+1./EXFAC-2.)
      ENDIF
      SUM=SUM+DBDT*FWEIGHT(K)/OPA
    1 CONTINUE
      OPARL=4.*STEBOL*TL*TL*TL/SUM
      RETURN
      END subroutine
      end module
