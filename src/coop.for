      module MOD_COOP

      type ROW
        integer :: LBKG_IDX
        real*8,allocatable :: V(:)
      endtype
      type(ROW),private,allocatable:: ROWS(:)

      contains

      SUBROUTINE COOP(XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $                OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $                N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $                ALPHA,SEXPO,AGAUNT,K,SIGMAKI,WAVARR,SIGARR,NF,NFDIM)

!***  NON-LTE CONTINUOUS OPACITY AT GIVEN FREQUENCY POINT K (XLAM)
!***  OPACITY AT DEPTH POINT L: OPAL
!***  SUMMATION OVER ALL DEPTH POINTS LEADS TO THE OPACITY AT GIVEN
!***  FREQUENCY POINT L
!**** N:      NUMBER OF LEVELS
!MH   ND:     NUMBER OF DEPTH POINTS
!MH   ND  = 1 WHEN COOP IS CALLED FROM OPAROSS
!MH   XLAM    : FREQUENCY (LOOP OVER NF IN OPAROSS)
!MH   T(ND)   : TEMPERATURE AT EACH DEPTH POINT ND
!MH   RNE(ND)  : RELATIVE ELECTRON DENSITY AT EACH DEPTH POINT ND
!MH   POPNUM(ND,N): ARRAY WITH POPULATION NUMBERS FOR EACH DEPTH POINT ND AND LEVEL N
!MH   ENTOT(ND) : TOTAL ELECTRON DENSITY FOR EACH DEPTH POINT ND
!MH   RSTAR   :
!MH   OPA     :
!MH   ETA     :
!MH   ak_bol   = 1.380622590000000E-016; Boltzmann constant in erg/K
!MH   EXPFAC   =  exp(h*nu/k/T)
!MH   W3      =W*W*W, W is wavenumber
!MH   C2      = 2 * H * C    ( G * CM**3 / S**2 )
!MH   C2*W3   = 2 * H * NU^3/ C^2
!MH   CALLED BY OPAROSS, CCORE, COMO, ETL, WRCONT, FIOSS
!MH   ND = 1 WHEN COOP IS CALLED FROM OPAROSS
!MH   CALLED BY OPAROSS:  LBKG = TRUE
!MH            ABLIN,EMLIN IS THEN WRITTEN TO *.LBKG FILE
!MH                      OPACITY.FOR HAS TO BE RUN IN ORDER TO BIN
!MH   CALLED BY FIOSS:   LBKG = FALSE
!MH                      SECOND TIME FIOSS IS RUN:
!MH                      LINE OPACITY, EMISSIVITY MUST NOT BE ADDED TO TOTAL OPACITY, EMISSIVITY
!***********************************************************************
!MH  VERY IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MH  ND = 1 WHEN COOP IS CALLED FROM OPAROSS!!!!!!!!!!!!!!
!***********************************************************************

      use MOD_PHOTOCS
      use MOD_GAUNTFF
      use MOD_HMINUSFF
      use MOD_LINSCA

      use common_block
      use file_operations
      use phys

      implicit none

      !global in/out
      INTEGER:: N, ND ! N = Number of levels, ND - number of depth points
      INTEGER, dimension(ND) :: IWARN

      INTEGER :: K, NF, NFDIM
      
      CHARACTER*10                        :: MAINPRO(ND), MAINLEV(ND)

      INTEGER, intent(in),  dimension(N)  :: NCHARG, NOM
      REAL*8,  intent(out), dimension(ND) :: OPA, ETA, THOMSON
      
      !global in
      
      REAL*8 :: SIGMAKI(NF, N), EINST(N, N)
      REAL*8 :: POPNUM(ND, N)

      REAL*8,dimension(N):: WEIGHT,ELEVEL,EION,ALPHA,SEXPO
      REAL*8,dimension(ND):: T, RNE, ENTOT

      real*8, dimension(N, NFDIM) :: SIGARR, WAVARR

      REAL*8       :: XLAM
      character*8  :: agaunt(N) 
      REAL*8       :: RSTAR
      CHARACTER*10 :: LEVEL(N)

      !functions
      REAL*8  :: SQRT ! FUNCTION

      !local
      REAL*8  :: G, WE, EDGE, T32, ROOTTL, TL, W, W3
      REAL*8  :: OPAL, OPAMAX, ETAL, EMINDU  ! induced emmission
      REAL*8  :: C1,C2,C3,CFF,AK_BOL, ELDEN, EXPFAC, PRESIG, SUM
      REAL*8  :: GFF(-1:6),GIIIX(-1:6)
      REAL*8  :: SCAFAC,ABSFAC,ALEMIS, ALCROSS
      real*8 :: LINOP(ND)
      INTEGER :: L, J, I ! loop variables
      INTEGER :: JM  ! JM = J-1
      INTEGER :: L_H0, IG
      REAL*8 :: SIGMAE, SIGMA, SIGMATH, SIGMAFF, SIGHMFF

!     C1 = H * C / K    (CM * ANGSTROEM)
      DATA AK_BOL, C1 /1.38062259d-16, 1.4388d0/

!     C2 = 2 * H * C (G * CM**3 / S**2)
!     C2 * W3 = 2 * H * NU^3 / C^2 (pre-factor of the Planck-function)
      DATA C2 /3.972d-16/

!     SIGMAE = ELECTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE /0.6652d-24/

!     C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON = 0.5*(h^2/(2*pi*m_e*k_bol))^(3/2)
      DATA C3 /2.07d-16/

!     CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION (ALLEN PAGE 100)
      DATA CFF /1.370d-23/

!     W = wavenumber in cm^(-1)
      W = 1.d8 / XLAM; W3 = W * W * W

!     Read in all LINOP at once
      CALL RDOPAC(XLAM, LINOP, NDPMIN, XLBKG1, XLBKG2)

!**********************************************************************
!***  LOOP OVER ALL DEPTH POINTS
      DO L = 1, ND

      OPAMAX=.0d0
      OPAL=.0d0
      ETAL=.0d0
      IWARN(L)=1H 
      TL=T(L)
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL

!***  BOUND-FREE  ******************************************************
!***  I = LOW      J = UP
      DO J=2,N
!***  UPPER LEVEL MUST BE GROUND STATE OF HIGHER IONIZATION STAGES
        IF ((NOM(J) .NE. NOM(J-1)) .OR. (NCHARG(J) .NE. NCHARG(J-1)+1)) GOTO 5
        JM=J-1
        DO I=1,JM
!***      LOWER LEVEL I MUST BELONG TO THE SAME ELEMENT AS UPPER LEVEL J
          IF (NOM(I) .NE. NOM(J)) GOTO 2
!***      CHARGES MUST DIFFER BY 1
          IF (NCHARG(J).NE.NCHARG(I)+1 ) GOTO 2
          EDGE=EION(I)-ELEVEL(I)
          IF (W.LT.EDGE) GOTO 2
!***      CALCULATE SIGMA, THE FREQUENCY-DEPENDENT CROSS SECTION
!***      IF ( K .GT. 0 ) IT IS ASSUMED THAT THE BOUND-FREE
!         CROSS SECTIONS SIGMA
!***      HAVE BEEN ALREADY CALCULATED ( ARRAY SIGMAKI )
          IF (K .GT. 0) THEN
              SIGMA=SIGMAKI(K,I)
          ELSE
              SIGMATH=EINST(I,J)*1.d-18
!****     **************************************************************
!***      Changes by Margit Haberreiter
!MH       SIGMA IN CM^2

              CALL PHOTOCS(SIGMA,SIGMATH,EDGE,W,ALPHA,SEXPO,AGAUNT,I,WAVARR,SIGARR,N,NFDIM)

          ENDIF
!***      RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
          
          WE=C3*RNE(L)*ENTOT(L)/T32

          !  EDGE-W = E_ion = \chi_ion (1/cm) -- Saha Boltzmann Equation
          ! (EDGE-W)*C1 = h*frequency/k_bolz  (K)
          ! TL = Temperature at Depth point L

          G=WEIGHT(I)/WEIGHT(J)*WE*EXP(C1*(EDGE-W)/TL)
          EMINDU=G*POPNUM(L,J)*SIGMA

          SUM=POPNUM(L,I)*SIGMA-EMINDU
          OPAL=OPAL+SUM
          ETAL=ETAL+EMINDU
!         ************************************
!***  LASER WARNING IF STIMULATED EMISSION EXCEEDS ABSORPTION IN ONE TRANSITION
          IF (SUM.LT.0d0) IWARN(L)=1H* 
          IF(SUM.LT.OPAMAX) GOTO 2
          OPAMAX=SUM
          MAINPRO(L)='BOUND-FREE'
          MAINLEV(L)=LEVEL(I)
    2   ENDDO
    5 ENDDO

!***  FREE-FREE  *******************************************************
!***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
      DO IG = 0,6
        IF (IG.GT.0) THEN
          CALL GAUNTFF (GIIIX(IG),IG,XLAM,TL)
        ELSE
          GIIIX(IG)=0d0
        ENDIF
        GFF(IG)=GIIIX(IG)*IG*IG
      ENDDO
!***  PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
      PRESIG =CFF/W3/ROOTTL
      EXPFAC=EXP(-C1*W/TL)
      l_h0 = 0
      DO I=1,N
!-hm-addition
!MH   coop HYDROGEN FF OPACITY
         if (level(i).eq.'H I......1') l_h0 = i
         SIGMAFF=PRESIG*GFF(NCHARG(I))

         SUM=RNE(L)*ENTOT(L)*POPNUM(L,I)*SIGMAFF
         EMINDU=SUM*EXPFAC
         SUM=SUM-EMINDU
         OPAL=OPAL+SUM
         ETAL=ETAL+EMINDU
         IF (SUM.GE.OPAMAX) THEN
        OPAMAX=SUM
        MAINPRO(L)='FREE-FREE'
        MAINLEV(L)=LEVEL(I)
         ENDIF
      ENDDO

      ETAL = ETAL * C2 * W3
!-hm-addition
!*** hier kommt H^- ff rein...
      if (l_h0.eq.0) then
         do i=1,N
           print *, 'Level(',i,')  = ', level(i)
         enddo
         print *,' HI ground level not found'
         stop ' no HI'
      else
         i=l_h0
!     opacity fit is per El-pressure and per Hydrogen atom
         call hminusff (sighmff,xlam,TL)
         sighmff=sighmff*ak_bol*Tl
!MH   NEW ENTRY FOR SUM

         SUM= RNE(L)*ENTOT(L)*POPNUM(L,I)*SIGhmff 
!**************************************************
          EMINDU=SUM*EXPFAC
!     fit includes stim correction:---  SUM=SUM-EMINDU
          OPAL=OPAL+SUM
!MH   MARGIT HABERREITER KORRECTION BY 1./(1.-EXPFAC)
!MH   W3=W*W*W, W is wavenumber
!MH   C2 = 2 * H * C    ( G * CM**3 / S**2 )
!MH   C2*W3 = 2 * H * NU^3/ C^2
          ETAL=ETAL+EMINDU*C2*W3/(1.d0-EXPFAC)
          if (sum.gt.opamax) then
            MAINPRO(L)='Hminus-ff'
            MAINLEV(L)=LEVEL(I)
          endif
      endif

!***  THOMSON SCATTERING ***********************************************
!MH   NEW ENTRY FOR SUM
!MH   SUM MUST NOT BE ADDED TO FORMER SUM!
!MH   as hminus ff is already added to ETAL
      SUM=RNE(L)*SIGMAE
!     WRITE (990,*) L,'RNE(L)*SIGMAE = ',RNE(L)*SIGMAE
!     add line-blanketing enhancement
      IF (XLAM.GT.227.83d0 .AND. XLAM.LE.3100d0 ) THEN

         ELDEN=RNE(L)*ENTOT(L)

         CALL LINSCA(XLAM, ELDEN, SCAFAC, ABSFAC)

        ! MULTIPLY THE TRUE ABSORPTION TO ACCOUNT FOR LINE ABSORPTION
         OPAL = OPAL * ABSFAC
        ! MULTIPLY SUM TO ACOUNT FOR LINE-SCATTERING
         SUM = SUM * SCAFAC
        ENDIF
        IF (SUM.GE.OPAMAX) THEN
          MAINPRO(L)='THOMSON'
          MAINLEV(L)='ELECTRON'
        ENDIF

        OPAL = OPAL + SUM

!MH    CHANGES BY MARGIT HABERREITER
!MH    RDOPAC: READS BINNED LINE OPACITY DATA FROM *.LBKG FILES
!MH    LINOP, LINEM AT DEPTHPOINT L AND FREQUENCY XLAM
!MH    IF coop IS CALLED FROM WRCONT
!        IF (LBKG) THEN

!MH    OUTWARD OF TEMPERATURE MINIMUM EMISSION AND ABSORPTION SET TO ZERO
          if (l .le. NDPMIN) then
            !alemis =0.
            alcross=linop(L)/ENTOT(L)
            alemis = BNUE(XLAM,T(NDPMIN))*alcross
          else
          !MH**  INWARD OF TEMPERATURE MINIMUM

            alcross = linop(L)/ENTOT(L)
            alemis  = BNUE(XLAM,T(L))*alcross
          endif
          OPAL=OPAL+alcross
          ETAL=ETAL+alemis

        ENDIF

        ! THOMSON = RELATIVE FRACTION FROM THE TOTAL OPACITY
        THOMSON(L) = SUM / OPAL
        OPA(L) = OPAL * ENTOT(L) * RSTAR
        ETA(L) = ETAL * ENTOT(L) * RSTAR
      ENDDO

      return

      end subroutine

      SUBROUTINE COOP_OPAROSS(XLAM,T,RNE,POPNUM,ENTOT,RSTAR,
     $                        OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $                        N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $                        ALPHA,SEXPO,AGAUNT,K,SIGMAKI,WAVARR,SIGARR,NF,NFDIM)

!***  NON-LTE CONTINUOUS OPACITY AT GIVEN FREQUENCY POINT K (XLAM)
!***  OPACITY AT DEPTH POINT L: OPAL
!***  SUMMATION OVER ALL DEPTH POINTS LEADS TO THE OPACITY AT GIVEN
!***  FREQUENCY POINT L
!**** N:      NUMBER OF LEVELS
!MH   ND:     NUMBER OF DEPTH POINTS
!MH   ND  = 1 WHEN COOP IS CALLED FROM OPAROSS
!MH   XLAM    : FREQUENCY (LOOP OVER NF IN OPAROSS)
!MH   T(ND)   : TEMPERATURE AT EACH DEPTH POINT ND
!MH   RNE(ND)  : RELATIVE ELECTRON DENSITY AT EACH DEPTH POINT ND
!MH   POPNUM(ND,N): ARRAY WITH POPULATION NUMBERS FOR EACH DEPTH POINT ND AND LEVEL N
!MH   ENTOT(ND) : TOTAL ELECTRON DENSITY FOR EACH DEPTH POINT ND
!MH   RSTAR   :
!MH   OPA     :
!MH   ETA     :
!MH   ak_bol   = 1.380622590000000E-016; Boltzmann constant in erg/K
!MH   EXPFAC   =  exp(h*nu/k/T)
!MH   W3      =W*W*W, W is wavenumber
!MH   C2      = 2 * H * C    ( G * CM**3 / S**2 )
!MH   C2*W3   = 2 * H * NU^3/ C^2
!MH   CALLED BY OPAROSS, CCORE, COMO, ETL, WRCONT, FIOSS
!MH   ND = 1 WHEN COOP IS CALLED FROM OPAROSS
!MH   CALLED BY OPAROSS:  LBKG = TRUE
!MH            ABLIN,EMLIN IS THEN WRITTEN TO *.LBKG FILE
!MH                      OPACITY.FOR HAS TO BE RUN IN ORDER TO BIN
!MH   CALLED BY FIOSS:   LBKG = FALSE
!MH                      SECOND TIME FIOSS IS RUN:
!MH                      LINE OPACITY, EMISSIVITY MUST NOT BE ADDED TO TOTAL OPACITY, EMISSIVITY
!***********************************************************************
!MH  VERY IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MH  ND = 1 WHEN COOP IS CALLED FROM OPAROSS!!!!!!!!!!!!!!
!***********************************************************************

      use MOD_PHOTOCS
      use MOD_GAUNTFF
      use MOD_HMINUSFF
      use MOD_LINSCA

      use common_block
      use file_operations
      use phys

      implicit none

      !global in / out
      INTEGER :: N ! Number of levels
      INTEGER :: IWARN

      INTEGER :: K, NF, NFDIM
      
      CHARACTER*10                        :: MAINPRO, MAINLEV

      INTEGER, intent(in),  dimension(N)  :: NCHARG, NOM

      REAL*8,  intent(out) :: OPA, ETA, THOMSON
      
      ! global in
      
      REAL*8 :: SIGMAKI(NF, N), EINST(N, N)
      REAL*8 :: POPNUM(N)

      REAL*8, dimension(N) :: WEIGHT, ELEVEL, EION, ALPHA, SEXPO

      REAL*8 :: T, RNE, ENTOT

      real*8, dimension(N, NFDIM) :: SIGARR, WAVARR

      REAL*8       :: XLAM
      character*8  :: agaunt(N) 
      REAL*8       :: RSTAR
      CHARACTER*10 :: LEVEL(N)

      !functions
      REAL*8  :: SQRT ! FUNCTION

      !local
      REAL*8  :: G, WE, EDGE, T32, ROOTTL, TL, W, W3
      REAL*8  :: OPAL, OPAMAX, ETAL, EMINDU  ! induced emmission
      REAL*8  :: C1,C2,C3,CFF,AK_BOL, ELDEN, EXPFAC, PRESIG, SUM
      REAL*8  :: GFF(-1:6),GIIIX(-1:6)
      REAL*8  :: SCAFAC,ABSFAC,ALEMIS, ALCROSS

      INTEGER :: L, J, I ! loop variables
      INTEGER :: JM  ! JM = J-1
      INTEGER :: L_H0, IG
      REAL*8 :: SIGMAE, SIGMA, SIGMATH, SIGMAFF, SIGHMFF

!     C1 = H * C / K    (CM * ANGSTROEM)
      DATA AK_BOL, C1 /1.38062259d-16, 1.4388d0/

!     C2 = 2 * H * C (G * CM**3 / S**2)
!     C2 * W3 = 2 * H * NU^3 / C^2 (pre-factor of the Planck-function)
      DATA C2 /3.972d-16/

!     SIGMAE = ELECTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE /0.6652d-24/

!     C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON = 0.5*(h^2/(2*pi*m_e*k_bol))^(3/2)
      DATA C3 /2.07d-16/

!     CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION (ALLEN PAGE 100)
      DATA CFF /1.370d-23/

!     W = wavenumber in cm^(-1)
      W = 1.d8 / XLAM; W3 = W * W * W

      OPAMAX = 0.0d0
      OPAL = 0.0d0
      ETAL = 0.0d0
      IWARN=1H 
      TL = T
      ROOTTL=SQRT(TL)
      T32=TL*ROOTTL

!***  BOUND-FREE ******************************************************
!***  I = LOW J = UP
      DO J=2,N
!***  UPPER LEVEL MUST BE GROUND STATE OF HIGHER IONIZATION STAGES
        IF ((NOM(J) .NE. NOM(J-1)) .OR. (NCHARG(J) .NE. NCHARG(J-1)+1)) GOTO 5
        JM=J-1
        DO I=1,JM
!***      LOWER LEVEL I MUST BELONG TO THE SAME ELEMENT AS UPPER LEVEL J
          IF (NOM(I) .NE. NOM(J)) GOTO 2
!***      CHARGES MUST DIFFER BY 1
          IF (NCHARG(J).NE.NCHARG(I)+1 ) GOTO 2
          EDGE=EION(I)-ELEVEL(I)
          IF (W.LT.EDGE) GOTO 2
!***      CALCULATE SIGMA, THE FREQUENCY-DEPENDENT CROSS SECTION
!***      IF ( K .GT. 0 ) IT IS ASSUMED THAT THE BOUND-FREE
!         CROSS SECTIONS SIGMA
!***      HAVE BEEN ALREADY CALCULATED ( ARRAY SIGMAKI )
          IF (K .GT. 0) THEN
              SIGMA=SIGMAKI(K,I)
          ELSE
              SIGMATH=EINST(I,J)*1.d-18
!****     **************************************************************
!***      Changes by Margit Haberreiter
!MH       SIGMA IN CM^2

              CALL PHOTOCS(SIGMA,SIGMATH,EDGE,W,ALPHA,SEXPO,AGAUNT,I,WAVARR,SIGARR,N,NFDIM)

          ENDIF
!***      RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
          
          WE = C3 * RNE * ENTOT / T32

          !  EDGE-W = E_ion = \chi_ion (1/cm) -- Saha Boltzmann Equation
          ! (EDGE-W)*C1 = h*frequency/k_bolz  (K)
          ! TL = Temperature at Depth point L

          G=WEIGHT(I)/WEIGHT(J)*WE*EXP(C1*(EDGE-W)/TL)
          EMINDU=G*POPNUM(J)*SIGMA

          SUM=POPNUM(I)*SIGMA-EMINDU
          OPAL=OPAL+SUM
          ETAL=ETAL+EMINDU
!         ************************************
!***  LASER WARNING IF STIMULATED EMISSION EXCEEDS ABSORPTION IN ONE TRANSITION
          IF (SUM.LT.0d0) IWARN=1H* 
          IF(SUM.LT.OPAMAX) GOTO 2
          OPAMAX=SUM
          MAINPRO = 'BOUND-FREE'
          MAINLEV = LEVEL(I)
    2   ENDDO
    5 ENDDO

!***  FREE-FREE  *******************************************************
!***  PRECALCULATE FREE-FREE GAUNT FACTORS FOR THE DIFFERENT ION CHARGES
      DO IG = 0,6
        IF (IG.GT.0) THEN
          CALL GAUNTFF (GIIIX(IG),IG,XLAM,TL)
        ELSE
          GIIIX(IG)=0d0
        ENDIF
        GFF(IG)=GIIIX(IG)*IG*IG
      ENDDO
!***  PRECALCULATE SIGMAFF, LEAVING OUT THE FACTOR NCHARGE*NCHARGE
      PRESIG =CFF/W3/ROOTTL
      EXPFAC=EXP(-C1*W/TL)
      l_h0 = 0
      DO I=1,N
!-hm-addition
!MH   coop HYDROGEN FF OPACITY
         if (level(i).eq.'H I......1') l_h0 = i
         SIGMAFF=PRESIG*GFF(NCHARG(I))

         SUM = RNE * ENTOT * POPNUM(I) * SIGMAFF
         EMINDU=SUM*EXPFAC
         SUM=SUM-EMINDU
         OPAL=OPAL+SUM
         ETAL=ETAL+EMINDU
         IF (SUM.GE.OPAMAX) THEN
         OPAMAX = SUM
         MAINPRO = 'FREE-FREE'
         MAINLEV = LEVEL(I)
         ENDIF
      ENDDO

      ETAL = ETAL * C2 * W3
!-hm-addition
!*** hier kommt H^- ff rein...
      if (l_h0.eq.0) then
         do i=1,N
           print *, 'Level(',i,')  = ', level(i)
         enddo
         print *,' HI ground level not found'
         stop ' no HI'
      else
         i=l_h0
!     opacity fit is per El-pressure and per Hydrogen atom
         call hminusff (sighmff,xlam,TL)
         sighmff=sighmff*ak_bol*Tl
!MH   NEW ENTRY FOR SUM

         SUM = RNE * ENTOT * POPNUM(I)*SIGhmff 
!**************************************************
          EMINDU=SUM*EXPFAC
!     fit includes stim correction:---  SUM=SUM-EMINDU
          OPAL=OPAL+SUM
!MH   MARGIT HABERREITER KORRECTION BY 1./(1.-EXPFAC)
!MH   W3=W*W*W, W is wavenumber
!MH   C2 = 2 * H * C    ( G * CM**3 / S**2 )
!MH   C2*W3 = 2 * H * NU^3/ C^2
          ETAL=ETAL+EMINDU*C2*W3/(1.d0-EXPFAC)
          if (sum.gt.opamax) then
            MAINPRO = 'Hminus-ff'
            MAINLEV = LEVEL(I)
          endif
      endif

!***  THOMSON SCATTERING ***********************************************
!MH   NEW ENTRY FOR SUM
!MH   SUM MUST NOT BE ADDED TO FORMER SUM!
!MH   as hminus ff is already added to ETAL
      SUM = RNE * SIGMAE

!     add line-blanketing enhancement
      IF (XLAM.GT.227.83d0 .AND. XLAM.LE.3100d0 ) THEN

         ELDEN = RNE * ENTOT

         CALL LINSCA(XLAM, ELDEN, SCAFAC, ABSFAC)

        ! MULTIPLY THE TRUE ABSORPTION TO ACCOUNT FOR LINE ABSORPTION
         OPAL = OPAL * ABSFAC
        ! MULTIPLY SUM TO ACOUNT FOR LINE-SCATTERING
         SUM = SUM * SCAFAC
        ENDIF
        IF (SUM.GE.OPAMAX) THEN
          MAINPRO = 'THOMSON'
          MAINLEV = 'ELECTRON'
        ENDIF

        OPAL = OPAL + SUM

        ! THOMSON = RELATIVE FRACTION FROM THE TOTAL OPACITY
        THOMSON = SUM / OPAL
        OPA = OPAL * ENTOT * RSTAR
        ETA = ETAL * ENTOT * RSTAR

      return

      end subroutine

      subroutine RDOPAC(XLAM,LINOP,NDPMIN,XLBKG1,XLBKG2)

!     PROGRAM READS THE BINNED LINE OPACITIES AND INCLUDES THEM IN THE HMINUS CODE TO CALCULATE THE CONTINUUMS OPACITIES

      use utils

      implicit none

      real*8, intent(out) :: LINOP(:)
      real*8, intent(in)  :: XLAM
      integer,intent(in)  :: XLBKG1,XLBKG2,NDPMIN
      integer             :: LAM
      character*50        :: flnam
      integer             :: i,IDX, nLINOP

      if(.not.allocated(ROWS)) allocate(ROWS(0))

!     ABBIN: BINNED LINE OPACITY -  OPACITY DISTRIBUTION FUNCTION
!     XLAM: WAVELENGTH GIVEN IN FGRID INCLUDING EDGE WAVELENGHTS

      LINOP = 0.0
      IDX = -1

!      if (xlbkg1 .ne. 0) print*, 'xlbkg1 = ', xlbkg1
!      if (xlbkg2 .ne. 0) print*, 'xlbkg2 = ', xlbkg2

      if (XLAM < XLBKG1) return
      if (XLAM > XLBKG2) return

      nLINOP = size(LINOP)
      LAM = XLAM
      LAM = ((XLAM)/10);   LAM = LAM*10+5

      !*** Search for frequency index
      SEARCH_IDX: do i=1,size(ROWS)
        if(ROWS(i)%LBKG_IDX == LAM) then
          IDX=i
          if(size(ROWS(i)%V) <= nLINOP) exit SEARCH_IDX !** reread the LINOP
          LINOP = ROWS(i)%V(1:nLINOP)
          return !** if found and correct size, return
        endif
      enddo SEARCH_IDX

      write(flnam, '(i6)') LAM

      flnam = './lbkg/'//trim(adjustl(flnam))//'.lbkg'

      write(*, '(A5,1x,F9.3,1x,I4,1x,A16)'), 'check', xlam, lam, flnam

      open(300, file = trim(adjustl(flnam)), status = 'old', action = 'read', err = 999)

      read(300, '(A)')           ! read the header
      !***********************************************************
      !*  L = 1 : MOST OUTWARD DEPTH POINT
      !*  L = ND: MOST INWARD DEPTH POINT
      !*  If depth point is outward of temperature minimum,
      !*  use LINOP of temperature minimum
      read(300,'(1pe12.5)') LINOP
      LINOP(1:NDPMIN) = LINOP(NDPMIN)
      close (300)
      !************************************************************
      if(IDX>0) then 
        !*** Change the size of ROWS(IDX)%V and replace with new value
        deallocate(ROWS(IDX)%V)
        allocate(ROWS(IDX)%V(nLINOP))
        ROWS(IDX)%V = LINOP
      else
        call append_ROW(LAM,LINOP)  ! ROW(end+1)%(LBKG_IDX, V) = (LAM, LINOP)
      endif

      if  (any(linop  <  0.)) then
        print *,xlam,'linop= ',linop
        call ERROR('rdopac: negative linop in ' // flnam)
      endif
      RETURN

      !***************************************************************
      !* Error handling
  999 print *, 'RDOPAC: IO Error opening file ' // flnam
      print '("XLBKG1,2 = ", i7," ",i7," xlam =",e10.4," lam=",i6)',
     $    XLBKG1,XLBKG2,xlam,lam
      call ERROR('RDOPAC: IO Error opening file ' // flnam)

      END subroutine

      !****************************************
      !* Append a row to the array of rows.
      subroutine append_ROW(LAM,LINOP)
      integer,intent(in) :: LAM
      real*8,intent(in) :: LINOP(:)
      type(ROW),allocatable :: old_ROWS(:)
      integer :: nrow
      nrow = size(ROWS)
      !* make a copy of the original, reallocate by n+1 and copy back the old one to the first n rows.
      allocate(old_ROWS(nrow))
      old_ROWS = ROWS
      deallocate(ROWS)
      allocate(ROWS(nrow+1))
      ROWS(1:nrow) = old_ROWS

      !* Set Data for the last row (LBKG_IDX, V)
      ROWS(nrow+1)%LBKG_IDX = LAM
      allocate(ROWS(nrow+1)%V(size(LINOP)))
      ROWS(nrow+1)%V = LINOP

      end subroutine append_ROW

      end module
