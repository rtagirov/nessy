      module MOD_COOP_M

      integer :: NDPMIN_WARN=0 !* set to true if an NDPMIN sanity check fails

      contains

      SUBROUTINE COOP_M(XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $                  OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $                  N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $                  ALPHA,SEXPO,AGAUNT,K,SIGMAKI,WAVARR,SIGARR,
     $                  LBKG,XLBKG1,XLBKG2,NF)
      use MOD_BNUE

!***  NON-LTE CONTINUOUS OPACITY AT GIVEN FREQUENCY POINT K (XLAM)
!***  OPACITY AT DEPTH POINT L: OPAL
!***  SUMMATION OVER ALL DEPTH POINTS LEADS TO THE OPACITY AT GIVEN
!***  FREQUENCY POINT L
!**** N:      NUMBER OF LEVELS
!MH   ND:     NUMBER OF DEPTH POINTS
!MH   ND  = 1 WHEN COOP_M IS CALLED FROM OPAROSS_M
!MH   XLAM    : FREQUENCY (LOOP OVER NF IN OPAROSS_M)
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
!MH   CALLED BY OPAROSS_M, CCORE, COMO, ETL, WRCONT, FIOSS8
!MH   ND = 1 WHEN COOP_M IS CALLED FROM OPAROSS_M
!MH   CALLED BY OPAROSS:  NLTELBKG=1
!MH            ABLIN,EMLIN IS THEN WRITTEN TO *.LBKG FILE
!MH                      OPACITY.FOR HAS TO BE RUN IN ORDER TO BIN
!MH   CALLED BY FIOSS8:   NLTELBKG=0
!MH                      SECOND TIME FIOSS8 IS RUN:
!MH                      LINE OPACITY, EMISSIVITY MUST NOT BE ADDED TO TOTAL OPACITY, EMISSIVITY
!***********************************************************************
!MH  VERY IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MH  ND = 1 WHEN COOP_M IS CALLED FROM OPAROSS_M!!!!!!!!!!!!!!
!***********************************************************************

      use MOD_PHOTOCS_M
      use MOD_GAUNTFF
      use MOD_HMINUSFF
      use MOD_LINSCA
      use MOD_RDOPAC

      IMPLICIT NONE
!MH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
!MH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
      !global in/out
      INTEGER:: N, ND ! max #depthpoints, #depthpoints ! N = Number of levels
      INTEGER, dimension(ND):: IWARN
      INTEGER :: K, NF
      
      INTEGER,dimension(N)  :: NCHARG, NOM
      CHARACTER*10          :: MAINPRO(ND),MAINLEV(ND)
      REAL*8, dimension(ND) :: OPA,ETA,THOMSON
      
      !global in
      
      REAL*8 :: SIGMAKI(NF, N), EINST(N, N)
      REAL*8 :: POPNUM(ND, N)

      REAL*8,dimension(N):: WEIGHT,ELEVEL,EION,ALPHA,SEXPO
      REAL*8,dimension(ND):: T, RNE, ENTOT
      REAL*8,dimension(N,NF) ::SIGARR,WAVARR
      REAL*8  :: XLAM
      character*8 :: agaunt(N) 
      integer XLBKG1,XLBKG2
      REAL*8  :: RSTAR
      LOGICAL :: LBKG
      CHARACTER*10 ::LEVEL(N)
      !functions
      REAL*8  :: SQRT  ! FUNCTION

      !local
      REAL*8  :: G, WE, EDGE, T32, ROOTTL, TL, W, W3
      REAL*8  :: OPAL, OPAMAX, ETAL, EMINDU  ! induced emmission
      REAL*8  :: C1,C2,C3,CFF,AK_BOL, ELDEN, EXPFAC, PRESIG, SUM
      REAL*8  :: GFF(-1:6),GIIIX(-1:6)
      REAL*8  :: SCAFAC,ABSFAC,ALEMIS, ALCROSS
      real*8 :: LINOP(ND)
      INTEGER :: L, J, I ! loop variables
      INTEGER :: JM  ! JM = J-1
      INTEGER :: L_H0, IG, NLTELBKG
      REAL*8  :: SIGMAE, SIGMA, SIGMATH, SIGMAFF,SIGHMFF
!***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA AK_BOL, C1 / 1.38062259d-16, 1.4388d0 /
!***  C2 = 2 * H * C    ( G * CM**3 / S**2 )
!MH   C2*W3 = 2 * H * NU^3 / C^2 (pre-factor of the Planck-function)
      DATA C2 / 3.972d-16 /
!***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652d-24 /
!***  C3 = RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
!        = 0.5 ( h^2/(2*pi*m_e*k_bol) )^(3/2)

      DATA C3 / 2.07d-16 /
!***  CFF = COEFFICIENT FOR FREE-FREE CROSS SECTION ( ALLEN PAGE 100 )
      DATA CFF / 1.370d-23 /
!***  W= wavenumber in cm^(-1)
!MH   NDPMIN IS THE DEPTH POINT OF THE TEMPERATURE MINIMUM
!MH   Model C,P: NDPMIN = 55
!MH   Model A:   NDPMIN = 53
!MH   TABASUN:   NDPMIN = 1
      integer:: NDPMIN 

       ! Read in all LINOP at once.
      W=1.d8/XLAM; W3=W*W*W
      IF ( ND.EQ.1 )THEN
        NDPMIN = 1
        NLTELBKG=0
      ELSE
!**********************************************************************
!***  LINEBLANKETING=0
!MH   DETERMINE DEPTH INDEX OF TEMPERATURE MINIMUM
        NDPMIN = 0       ! Sanity check 1
        DO L=2,ND-1       ! Find minimum
          if ((T(L) .LT. T(L-1)) .AND. (T(L) .LT. T(L+1))) THEN
            NDPMIN = L
          endif
        ENDDO
        if (LBKG) then
          NLTELBKG=1
        else
          NLTELBKG=0
        endif
        if (NDPMIN .eq. 0) then     ! Sanity check 1 - finish
          if(NDPMIN_WARN==0)
     $      print *,'coop_m: something wrong! NDPMIN = ',NDPMIN, 
     $              ' setting NDPMIN to 1'
          !pause
          NDPMIN = 1      ! continue, dont abort
          NDPMIN_WARN=NDPMIN_WARN+1
        endif
        if ((T(NDPMIN) .gt. 5000d0) .and. (ND .gt. 1)) then !Sanity check 2
          print *,'something wrong with Temperature minimum!',
     $                 NDPMIN, T(NDPMIN)
          write (6,*)'something wrong with Temperature minimum!',
     $                   NDPMIN, T(NDPMIN)
          pause
        endif
      ENDIF

      CALL RDOPAC(XLAM,LINOP,NDPMIN,XLBKG1,XLBKG2)

!**********************************************************************
!***  LOOP OVER ALL DEPTH POINTS
      DO L=1,ND
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
        IF ((NOM(J) .NE. NOM(J-1)) .OR. 
     $      (NCHARG(J) .NE. NCHARG(J-1)+1)) GOTO 5
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

              CALL PHOTOCS_M(SIGMA,SIGMATH,EDGE,W,ALPHA,SEXPO,AGAUNT,I,
     $                       WAVARR(1 : N, 1 : NF),SIGARR(1 : N, 1 : NF), N, NF)

          ENDIF
!***      RECIPROCAL STATISTICAL WEIGHT OF FREE ELECTRON
          
          WE=C3*RNE(L)*ENTOT(L)/T32

          !  EDGE-W = E_ion = \chi_ion (1/cm) -- Saha Boltzmann Equation
          ! (EDGE-W)*C1 = h*frequency/k_bolz  (K)
          ! TL = Temperature at Depth point L

          G=WEIGHT(I)/WEIGHT(J)*WE*EXP(C1*(EDGE-W)/TL)
          EMINDU=G*POPNUM(L,J)*SIGMA

!          if (L .eq. 64 .and. I .eq. 2) print*, 'coop_m 1: ', popnum(64, i)

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
!MH   COOP_M HYDROGEN FF OPACITY
         if (level(i).eq.'H I......1') l_h0 = i
         SIGMAFF=PRESIG*GFF(NCHARG(I))

!         if (L .eq. 64 .and. I .eq. 2) print*, 'coop_m 1: ', popnum(64, i)

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

      ETAL=ETAL*C2*W3
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

        !***  MULTIPLY THE TRUE ABSORPTION TO ACCOUNT FOR LINE ABSORPTION
         OPAL=OPAL*ABSFAC
        !***  MULTIPLY SUM TO ACOUNT FOR LINE-SCATTERING
         SUM=SUM*SCAFAC
        ENDIF
        IF (SUM.GE.OPAMAX) THEN
          MAINPRO(L)='THOMSON'
          MAINLEV(L)='ELECTRON'
        ENDIF
        OPAL=OPAL+SUM
!**********************************************************************
!MH**  CHANGES BY MARGIT HABERREITER
!MH**  RDOPAC: READS BINNED LINE OPACITY DATA FROM *.LBKG FILES
!MH**  LINOP, LINEM AT DEPTHPOINT L AND FREQUENCY XLAM
!MH**  IF COOP_M IS CALLED FROM WRCONT
        IF (NLTELBKG .EQ. 1) THEN
          ! PRINT *,l,xlam, 'READING NON-LTE LINEBLANKETING DATA'
          !MH**  LINEMMIN: EMISSIVITY AT DEPTH POINT 84 (MODELC)
          !MH** LINEMMIN: NOT YET DEVIDED BY ENTOT(L)
          !MH** check temperature minimum
          if (NDPMIN .NE. 1 .AND. (  T(NDPMIN) .gt. 5000d0)) then
            print *, 'coop_m: NDPMIN wrong! ', NDPMIN, T(NDPMIN)
            pause
          endif

          !MH**  OUTWARD OF TEMPERATURE MINIMUM EMISSION AND ABSORPTION SET TO ZERO
          if (l .le. NDPMIN) then
            !alemis =0.
            alcross=linop(L)/ENTOT(L)
            alemis = BNUE(XLAM,T(NDPMIN))*alcross
          else
          !MH**  INWARD OF TEMPERATURE MINIMUM
          ! if (l .gt. NDPMIN) then
            alcross = linop(L)/ENTOT(L)
            alemis  = BNUE(XLAM,T(L))*alcross
          endif
          OPAL=OPAL+alcross
          ETAL=ETAL+alemis

        ENDIF
        !***  THOMSON = RELATIVE FRACTION FROM THE TOTAL OPACITY
        THOMSON(L)=SUM/OPAL
        OPA(L)=OPAL*ENTOT(L)*RSTAR
        ETA(L)=ETAL*ENTOT(L)*RSTAR
      ENDDO

      RETURN

      END subroutine

      end module
