      PROGRAM FIOSS

!      use utils

      use MOD_COOP
      use MOD_DATOM
      use MOD_READMOD
      use MOD_READPOP
      use MOD_READRAD
      use MOD_COMPRO
      use MOD_CWGAUSS9
      use MOD_DECF_SYN
      use MOD_DIFFUS
      use MOD_ELIMIN
      use MOD_INTRFC
      use MOD_PRIDAT
      use MOD_WMOM0_F
      use MOD_TRAPLO
      use MOD_SABOLT
      use MOD_PREPR_F
      use MOD_OBSINT
      use MOD_LPPLOT
      use MOD_PRIPRO
      use MOD_PREF_SYN
      use MOD_TICTOC
      use OPINT
      use CONSTANTS,only:CLIGHT_SI

      use vardatom_full
      use vardatom_nlte
      use varhminus
      use common_block
      use file_operations
      use auxfioss
      use mod_synopa
      use mod_quant
      use geo_mesh

!             intrfc passes new a AGAUNT !!
!*******************************************************************************
!***  Formal Integral in the Observers frame for Sperically Symmetric geometry
!*******************************************************************************
!  FIOSS PC version: new MODFILE and POPNUM read
!  FIOSS7.01            routine INTRFC with fix of quantum numbers of HeI
!  FIOSS7-alpha/SYNSPEC calling obsint
!
!  FIOSS6-alpha/SYNSPEC new feature Tion_pot and dil for non-LTE populations
!
!  FIOSS5-alpha/SYNSPEC alpha version of FIOSS5
!
!  FIOSS5/SYNSPEC Includes treating of clumping (via use of p-grid)
!                     and electron-scattering iterations
!         FIOSS5/SYNSPEC is the SYNSPEC-opacity version of test-fioss9 in crlib
!
!  FIOSS4 is the version used for the SS-paper
!  difference between fioss3 and fioss4 is mainly in the interpolation of between
!  2 grid points:
!                fioss3: interpolate between the nearest grid points
!                fioss4: take nearest grid point 
!                (to find the location in the code search for "nearest")
!                also in fioss4: modified TRAPLO: prints nr of entries
!
!MH*	XLAM: WAVELENGTH IN VACUUM
!*******************************************************************************

      implicit real*8(a - h, o - z)

      !***  Constant for thermal Gauss-Profile (= m(e)/(4k)) (cgs?)
      PARAMETER (GAUKONST = 1.649538d-12)

      !*** in principle NFODIM should be equal to NFMAX in OPINT.FOR
      !        however since the routine XY likes to cut the number of
      !        points in order to be able to add other lines, it is a
      !        good idea to make NFODIM quite a bit larger

      real*8, allocatable, dimension(:) ::       R

      real*8, allocatable, dimension(:) ::       VDU

      real*8, allocatable, dimension(:) ::       eta, opa, thomson, dil

      real*8, allocatable, dimension(:) ::       dummy1_nf, dummy2_nf, dummy3_nf, dummy4_nf

      real*8, allocatable, dimension(:) ::       ZRAY, RRAY

      real*8, allocatable, dimension(:) ::       XCMF

      integer, allocatable, dimension(:) ::      levelpl, nfedge, itne, iwarn, irind

      integer, allocatable, dimension(:, :) ::   iback

      real*8, allocatable, dimension(:, :) ::    WLK

      real*8, allocatable, dimension(:, :) ::    Z2D

      real*8, allocatable, dimension(:) ::       Z1D

      real*8, allocatable, dimension(:, :) ::    Tion_pot

      real*8, allocatable, dimension(:) ::       VERTVELO, VELOVAR

      real*8 :: mu

      integer :: NVD

      LOGICAL :: LOPA

      CHARACTER KARTE*80, MODHEAD*104, PHEAD*28

      CHARACTER flnam*70

      CHARACTER CLVFLNAM*70

      LOGICAL NORM,PLOT,TRANS,FIN,PROLIB,COROT,VAR,ADDVELO
      INTEGER LMAX
      LOGICAL EXLOOP

      LOGICAL :: NTF, NFE

      INTEGER :: IND, J, MLC

      INTEGER :: MP, SP, EP, IS, FP, N_CLV, LP

      REAL*8 :: XLMIN, XLMAX

      !***  CLIGHT = SPEED OF LIGHT IN KM/SECOND
      real*8,parameter :: CLIGHT =CLIGHT_SI/1d3
      DATA PI / 3.141592654 /
      !***  SIGMAE = ELCTRON SCATTERING CROSS SECTION  ( CM**2 )
      DATA SIGMAE / 0.6652d-24 /
      !***  DEFAULT VALUES
      DATA LSOPA,LSDWL,LSPRO,IFIRST,NFOBS,NPHI/-1,-1,-1,1,65,1/
      DATA NPHIP,LPSTI,LPSTA,LPENI,LPEND/45,0,1,0,1/
      DATA JFIRSI,JLASI/0,0/

      real*8,allocatable :: dummy2(:, :) ! do not allocate - make sure noone uses it

      REAL*8, ALLOCATABLE :: WAV_CLV(:), FLUX_CLV(:, :)
       
      real*8 :: dummy0

      real*8, allocatable :: PROFILE(:),PROFN(:),DLAM(:),EMINT(:)

      logical :: CORE,FILE_EXIST
      real*8, allocatable, dimension(:,:) :: XJK,CWK,DINT,XJ
      real*8, allocatable, dimension(:) :: XNU

      character (len = 8) :: rwlae_str

      !MH   VARIABILITY OF VELOCITY ACTIV/NONACTIVE
      !MH var = true then variation of Doppler velocity considered
      ! VAR = .TRUE.
      VAR	= .FALSE. 
      !MH ADDVELO TRUE then ASPLUND 2000, A&A 359, 729 is considered
      ADDVELO = .FALSE.
      ! ADDVELO = .TRUE.
      LOPA	= .TRUE.
!***  TAUMAX = MAXIMUM OPTICAL DEPTH WHERE INTEGRATION IS TRUNCATED
!     taumax is not used anymore since all points are needed for an
!     evaluation of the electron scattering integral
!***  XMAX = width in doppler units added to the line interval
!     only for line line index selections
!***  XN = RESOLUTION ELEMENTS WITHIN (MINIMUM) TURBULENCE VELOCITY
!     so far, in all tests there was no significant improvement by 
!     using more than 3 points per V-DOP
      XN = 2.

      XMAX = 3.

      CALL TIC()

      call datom('full',
     $           N,
     $           LEVEL,
     $           NCHARG,
     $           WEIGHT,
     $           ELEVEL,
     $           EION,
     $           MAINQN,
     $           EINST,
     $           ALPHA,
     $           SEXPO,
     $           AGAUNT,
     $           COCO,
     $           KEYCOL,
     $           ALTESUM,
     $           INDNUP,
     $           INDLOW,
     $           LASTIND,
     $           NATOM,
     $           ELEMENT,
     $           SYMBOL,
     $           NOM,
     $           KODAT,
     $           ATMASS,
     $           STAGE,
     $           NFIRST,
     $           NLAST,
     $           WAVARR,
     $           SIGARR,
     $           eleatnum,
     $           levatnum,
     $           NFDIM)

      call datom('nlte',
     $           N_nlte,
     $           level_nlte,
     $           ncharg_nlte,
     $           weight_nlte,
     $           elevel_nlte,
     $           eion_nlte,
     $           mainqn_nlte,
     $           einst_nlte,
     $           alpha_nlte,
     $           sexpo_nlte,
     $           agaunt_nlte,
     $           coco_nlte,
     $           keycol_nlte,
     $           altesum_nlte,
     $           indnup_nlte,
     $           indlow_nlte,
     $           lastind_nlte,
     $           natom_nlte,
     $           element_nlte,
     $           symbol_nlte,
     $           nom_nlte,
     $           kodat_nlte,
     $           atmass_nlte,
     $           stage_nlte,
     $           nfirst_nlte,
     $           nlast_nlte,
     $           wavarr_nlte,
     $           sigarr_nlte,
     $           eleatnum_nlte,
     $           levatnum_nlte,
     $           nfdim) ! NFDIM is known from varhminus module

!     READING OF THE MODEL FILE ------------------------------------

      ND = read_param('ND')
      NP = read_param('NP')
      NF = read_param('NF')

      ndaddim = 2 * ND + 12

      allocate(entot(ND))
      allocate(T(ND), RNE(ND))
      allocate(XJC(ND), XJCARR(ND, NF), WCHARM(ND, NF), XJL(ND, lastind_nlte))
      allocate(R(ND), TAUROSS(ND))
      allocate(VELO(ND), GRADI(ND), VDU(ND))

      allocate(HTOT(ND), GTOT(ND), XTOT(ND), ETOT(ND))

      allocate(eddi(3, ND), eddarr(3, ND, NF))

      allocate(dummy1_nf(NF), dummy2_nf(NF), dummy3_nf(NF), dummy4_nf(NF))

      allocate(P(NP), Z1D(ND * NP), Z2D(ND, NP))

      allocate(popnum(ND, N), pop1(ND, N), pop2(ND, N), pop3(ND, N))

      allocate(eta(ND), opa(ND))

      allocate(thomson(ND))

      allocate(U(ND, NP), WLK(ND, NP))

      allocate(ZRAY(ndaddim), xcmf(ndaddim), rray(ndaddim))

      allocate(IRIND(NDADDIM), IBACK(NDADDIM, NP))

      allocate(Tion_pot(ND, 3), dil(ND))

      allocate(ITNE(ND), IWARN(ND))

      allocate(NFEDGE(N), LEVELPL(N))

      allocate(mainpro(ND), mainlev(ND))

      allocate(abxyz(NATOM))

      allocate(vertvelo(ND))

      allocate(velovar(ND))

      IFL = 3; open(IFL, file = 'MODFILE', STATUS='OLD'); rewind IFL

      lblank=0

      CALL READMOD(IFL,N,ND,TEFF,R,NP,P,Z1D,ENTOT,VELO,GRADI,RSTAR,VDOP,NF,
     $             dummy1_nf(1 : NF),dummy2_nf(1 : NF),dummy3_nf(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      close(IFL)

      open(unit = 3, file = '../mu')

      read(3, *) mu

      close(3)

      R = R - 1.0

      R = R / mu

      R = R + 1.0

      Z1D = Z1D - 1.0

      Z1D = Z1D / mu

      Z1D = Z1D + 1.0

      GRADI = mu * GRADI

      IFL = 3; open(IFL, file = 'POPNUM', STATUS='OLD')
!     also reads number of depth points, temperature, rne
      call readpop(ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,jobnum)
      close (ifl)

      NDPMIN = tempmin(T, ND)

!     read the radiation field from files RADIOC and RADIOL (pop1 is used as dummy storage)
      CALL READRAD(NF,ND,POP1,XJCARR,XJC,XJL,
     $             HTOT,GTOT,XTOT,ETOT,dummy4_nf,TOTIN,TOTOUT,
     $             NCHARG_nlte,EDDARR,EDDI,NOM_nlte,WCHARM,N_nlte,lastind_nlte,
     $             EINST_nlte,MODHEAD,JOBNUM)

      !***  READING VERTICAL VELOCITY ASPLUND 2000, A&A 359, 729
      if (ADDVELO) then
        open (IFL,file='VELO',STATUS='OLD')
        DO L=1,ND
          read (UNIT=IFL, fmt=*), XX,VERTVELO(L),YY
        ENDDO
        close (IFL)
      endif
      !*** READ VARIATION OF VELOCITY BY RICHARD WACHTER
      IF (VAR) THEN
        open (IFL,file='VELOVAR',STATUS='OLD')
        DO L=1,ND
          read (UNIT=IFL, fmt=*), VELOVAR(L)
        ENDDO
        close (IFL)
      ENDIF

      call ZGRID(R, P, Z1D, ND, NP) ! calculate the Z1D grid

      Z2D = RESHAPE(Z1D, (/ND, NP/))

!     PRINTOUT OF THE ATOMIC DATA
      call PRIDAT(N,LEVEL,NCHARG, WEIGHT,ELEVEL,EION,EINST,
     $            KODAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $            NATOM,ELEMENT,NOM,ABXYZ,ATMASS)

      !***  ADDITIONAL CALCULATION FOR SUBROUTINE PRIPRO
      RSTAR2 = RSTAR * RSTAR

      !***  DECODING INPUT CARDS
      JFIRST=1
      JLAST=NP-1
      NORM=.FALSE.
      PLOT=.FALSE.
      TRANS=.FALSE.
      PROLIB=.FALSE.
      FIN=.FALSE.
      COROT=.FALSE.
      PHI=PI/2.
      iTionsel=-1
      !***  wavelength-shift of the profile wavelength scale
      shift=0.
      !***  adding spectral range for electron scattering wings
      FELSCA=0.
      !***  MAXITER = max iterations for the el'scattering contribution
      MAXITER = 1

      CALL READ_NLTE_TRA(NFE, NUMTRA, ELID, CHARGE, SWL, SWU, AUL, WAVTRA)

      !***  LOOP FOR EVERY DETECTED LINE - OPTION CARD   ---------------
      open (1,file='CARDS',STATUS='OLD')
      vdopp=vdop
      VSINI=0.
      MLC = 0
      MAIN_LOOP: DO  ! FIN exits the loop
      vdop =vdopp
      FMAX = 0.
      FMIN = 0.
      CALL DECF_SYN (KARTE,PLOT,NFOBS,LSOPA,FMAX,FMIN,
     $        MAXITER,felsca,RWLAE,PHEAD,PROLIB,VSINI,SHIFT,
     $        LSPRO,LSDWL,NORM,TRANS,FIN,VDOP,
     $        NPHIP,LPSTI,LPENI,JFIRSI,JLASI,COROT,iTionsel,XLMIN,XLMAX)

!RT:  BLOCK FOR CALCULATION OF NLTE LINES LISTED IN THE DATOM FILE
!RT:  FOR FURTHER DETAILS SEE LINOP.FOR
!===========================================================================

      NTC = 0

      IF (MLC .EQ. 0 .AND. NFE .EQ. .TRUE.) THEN

          NTF = .FALSE.

          DO IND = 1, NUMTRA

             IF (WAVTRA(IND) .GE. XLMIN .AND. WAVTRA(IND) .LE. XLMAX) THEN 

                 NTC = NTC + 1

                 NTF = .TRUE.

             ENDIF

          ENDDO

          IF (NTF) THEN 

!              ALLOCATE(NLTE_ABS(NTC))
!              ALLOCATE(NLTE_EMI(NTC))
              ALLOCATE(NTI(NTC))

              J = 1

              DO IND = 1, NUMTRA

                 IF (WAVTRA(IND) .GE. XLMIN .AND. WAVTRA(IND) .LE. XLMAX) THEN

                     NTI(J) = IND; J = J + 1

                 ENDIF

              ENDDO

              CALL READ_NLTETRAPOP(NTI, NTPL, NTPU, NTDL, NTDU)

!              DEALLOCATE(NTI)

          ENDIF

      ENDIF

      MLC = MLC + 1

!===============================================================================

      if(allocated(PROFILE )) deallocate(PROFILE)
      if(allocated(PROFN   )) deallocate(PROFN  )
      if(allocated(DLAM    )) deallocate(DLAM   )
      if(allocated(EMINT   )) deallocate(EMINT  )
      allocate(PROFILE(NFOBS))
      allocate(PROFN(NFOBS))
      allocate(DLAM(NFOBS))
      allocate(EMINT(NFOBS))
      if (FIN) exit MAIN_LOOP ! that is GOTO END

!***  Calculate the ionisation temperatures
      print *,' Tion-sel ', itionsel
      call SABOLT(ENTOT,RNE,POPNUM,T,ND,N,NCHARG,WEIGHT,ELEVEL,EION,
     $            KODAT,NOM,NATOM,XTOT,HTOT,GTOT,R,MODHEAD,JOBNUM,
     $            Npot,Tion_pot,dil,teff,iTionsel)








      vdopp=vdop
      print *
      print *,' Macro turbulence used to calculate the profile:'
      !***  VDOP IS USED AS REFERENCE TO CALCULATE FREQUENCY GRID
      !***  VDOPP IS PHYSICAL VELOCITY 


      vdop=vdopp
      print *,'                   V-turb = ',vdopp,' km/s'
      ! VDOP is used to calculate the integration X-grid
      print *,'                   V-grid = ',vdop ,' km/s'

      !***  ADDING VERTICAL VELOCITY TO VELO(L)
      amaxvdu=0.

      !**********************************************************
      !***  CHANGES by Margit Haberreiter 22.8.03:
      !***  replace expansion velocity by tabulated
      !***  horizontally averaged vertical velocity VERTVELO(L) IN KM/S
      !***  VELOVAR: VARIATION OF VERTICAL VELOCITY IN M/S
      !***  Reference: Asplund et al. 2000, A&A 359, 729
      IF(VAR) THEN
        VERTVELO(1:ND)=VERTVELO(1:ND)+(VELOVAR(1:ND)/1000.)

        PRINT *,'FIOSS: VELOVAR(:)',VELOVAR(1:ND),VERTVELO(1:ND)
        VELO(1:ND)= VERTVELO(1:ND)
      ENDIF
      ! VELO(L)= VERTVELO(L)*5.
      if (ADDVELO) then
        VELO(1:ND)= VERTVELO(1:ND)
        print *, 'FIOSS: ASPLUND VERTICAL VELO ACTIVE'
      else
        !**********************************************************
        !***  INTRODUCING DIMENSIONLESS VELOCITY UNITS
        print *, 'FIOSS: ASPLUND VERTICAL VELO NOT ACTIVE'
        VDU (1:ND)=VELO(1:ND)/VDOP
      endif
      amaxvdu=maxval(abs(VDU))
      print *,' Max-VDU= ',amaxvdu

      !*** FMAX input is in km/s
      if (FMAX.gt.0.) FMAX=FMAX/vdop
      if (FMIN.gt.0.) FMIN=FMIN/vdop
      !*** V-sini is in km/s
      VSIDU=VSINI/VDOP
      IF (VSINI.NE.0.) THEN








         NPHI=NPHIP
         LPEND=NPHI
         IF (LPSTI.NE.0) LPSTA=LPSTI
         IF (LPSTI.NE.0) LPEND=LPSTI
         IF (LPENI.NE.0) LPEND=LPENI
      ENDIF
      IF (JFIRSI.NE.0) JFIRST=JFIRSI
      IF (JFIRSI.NE.0) JLAST=JFIRSI
      IF (JLASI .NE.0) JLAST=JLASI
      !*** width of electron scattering (in doppler units)
      esca_wd = sqrt(T(ND)/GAUKONST)/1.e5/vdop
      ! esca_wd = esca_wd * 2.3
      esca_wd = esca_wd * felsca
      print *,' e-folding width of e-scattering ',esca_wd,' [V-Dop] = ',
     &         esca_wd*vdop,' [km/s]'

      xlam=rwlae
      print*, 'RWLAE = ', rwlae

!***  PREPARATION OF LINE QUANTITIES (ALSO FOR BLENDING LINES)
      CALL PREF_SYN(KARTE,N_nlte,ELEVEL_nlte,LINE,INDLOW_nlte,INDNUP_nlte,LASTIND_nlte,
     $              VDOP,FMAX,FMIN,XMAX,VDU(1),VSIDU,esca_wd,
     $              DXOBS,NFOBS,XLAM,FREMAX,
     $              NF,FNUEC)

      IF (LINE .EQ. 0) cycle MAIN_LOOP       ! go back to DECF_SYN
!***  replace the wavelength XLAM by the reference RWLAE
!MH*  XLAM: WAVELENGTH IN VACUUM
      IF (RWLAE.GT.0.) THEN
         PRINT *,' ',XLAM,' replaced by ',RWLAE
         XLAM=RWLAE
      ENDIF
    
      !***  DEFINING ZERO-POINT OF THE OBSERVER'S FRAME FREQUENCY
      xobs0 = FREMAX-DXOBS
      CALL DIFFUS (XLAM,T,R,ND,BCORE,DBDR)   !BCORE=Plank (XLAM, T) at R(ND), DBDR=d(BCORE)/dR at R=ND
      ncoop=n
      CALL COOP(XLAM,ND,T,RNE,POPNUM,ENTOT,RSTAR,
     $          OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $          N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $          ALPHA,SEXPO,AGAUNT,0,DUMMY2,
     $          WAVARR,SIGARR,NF,NFDIM)

!     CALCULATION OF THE CONTINUUM RADIATION FIELD XJC AT THE LINE FREQUENCY
      CALL ELIMIN(XLAM,FNUCONT,DUMMY0,U,Z2D,XJC,R,P,BCORE,DBDR,OPA,ETA,THOMSON,EDDI,ND,NP)

      print *,' Continuum Flux interpolated from the model: ',FNUEC
      print *,'      "      "  from ELIMIN ',FNUCONT,
     $        ' at l=',XLAM
      refcon=fnucont




      PRINT *,' VERSION 8 / SYNSPEC'
      print *,' xobs0, dxobs ',xobs0,dxobs
      do K=1,NFOBS
        XO=XOBS0+K*DXOBS
        DLAM(K)=-XO*XLAM*VDOP/CLIGHT
        PROFILE(K)=0.0
      enddo

!*** create the wavelength array for the opacity routine
!       note: the dlam array starts with the zero element
!             the vopa array with the first
!***  CHANGES by Margit Haberreiter 22.8.03
!***  Fehler bei neg vel:      !     vopa0 = xobs0+dxobs       + vdu(1) + 1./xn/1000. + vsidu 
!***  negative velocities: vdu(1) not max. value
!***  instead of vdu(1) use amaxvdu
!***  groesster abs Wert von VDU
      vopa0 = xobs0+dxobs       + vdu(1) + 1./xn/1000. + vsidu
      nvopa = vopa0*xn 
      nvopa = nvopa + 1
      print *, '1. FIOSS NVOPA=', NVOPA
      vopa0 = real(nvopa)/xn 
!***  changed by MH 22.8.03
!      vopam = xobs0+dxobs*NFOBS - vdu(1) - vsidu
      vopam = xobs0+dxobs*NFOBS -amaxvdu - vsidu
      print *,' xobs0, xobsm, NFOBS'
      print *,  xobs0,xobs0+dxobs*NFOBS, NFOBS
      dvopa = -1./xn
      nvopa = (vopam - vopa0)/dvopa + 1 + 1
      PRINT*, 'NVOPA 2 = ', NVOPA
      !call setNVDIM(NVOPA)
      call OPINT_INIT(NVOPA,NFOBS)
      print *, '2. FIOSS NVOPA=', NVOPA
      print *,' vopa0,   vopam,   dvopa,    nvopa,   nvdim'
      print *,vopa0,vopam,dvopa,nvopa,nvdim()
      
      print *, '2. FIOSS NVOPA=', NVOPA

      NVD = NVDIM()

      allocate(XJK(ND, NVD))
      allocate(CWK(ND, NVD))
      allocate(DINT(NVD, 2 * ND))
      allocate(XJ(ND, NVD))
      allocate(XNU(NVD))

      if (nvopa.gt.nvdim()) then
         print *,'nvopa,nvdim',nvopa,nvdim()
         print *,' opa/eta - dimension insufficient'
         stop ' change dimensions'
      endif

      do kopa=1,nvopa
         vopa(kopa)=vopa0+(kopa-1)*dvopa
         XNU(Kopa)=vopa(kopa)*VDOP*1.E5
      enddo

!***  Compute Integration Weights for thermal Gauss Profile
      dstep=abs(dvopa*vdop*1.e5)
      print *, ' vdop ',vdop,' dstep ',dstep

      if (maxiter.gt.1) CALL CWGAUSS9(CWK,XNU,ND,NVOPA,T,xjc,dstep)

      if(any(isnan(XJC))) stop 'fioss: NAN XJC'
      if(any(isnan(CWK))) stop 'fioss: NAN XJC'
      if(any(isnan(XNU))) stop 'fioss: NAN XJC'
      if(any(isnan(T))) stop 'fioss: NAN XJC'

!     store quantities for Ivan's routines
      call quant(nd,n,t,rne,entot,popnum,vdop,vdopp,nvopa,rwlae,Npot,Tion_pot,dil,teff,xjc)

      !*** call ivan's routines
      rewind 555
      rewind 19
      rewind 24
      rewind 25
      rewind 26
      rewind 27
      rewind 28
      rewind 29
      rewind 55
      rewind 56

      call intrfc(ncharg,weight,elevel,eion,
     *            einst,alpha,sexpo,agaunt,natom,
     *            symbol,nfirst,nlast, WAVARR, SIGARR, N, NFDIM)

      PRINT *,'FIOSS: Time elapsed after INTRFC: ',TOC()

      !$$$ calculate opacities
      !*****************************************************************
      !***  MARGIT HABERREITER
      !***  FROM SYNOPA OPAC IS CALLED THE LAST TIME IN FIOSS RUN
      !***  OPENING FILE FOR EMLIN, ABLIN OUTPUT
      !***  WRITTEN IN OPAC, CALLED BY SYNOPA
      !***  FREQUENCY GRID WRITTEN IN INTRFC
      write (flnam,'(F14.0)') XLAM
      flnam = adjustl(trim(flnam)//'lopa')
      open (200,file=flnam,status='unknown')
      !*** micha: Open the ABEMLIN File according to the setting in CARDS
      !***  1) get the full filename
      write(flnam,'(F14.0)') XLAM
      flnam=adjustl(trim(flnam)//'abemlin')
      flnam=adjustl(adjustr(CARDS%ABEMLIN_PATH)//flnam)
      !*** block for ABEMLIN write/read
      !***  2) decide what to do - auto: read, write
      if(CARDS.ABEMLIN==CARD_PARAMS.ABEMLIN_AUTO) then !** automode -
          inquire(file=flnam,exist=FILE_EXIST)
          if(FILE_EXIST) then             !**   if file exist read,
            CARDS.ABEMLIN=CARD_PARAMS.ABEMLIN_READ

          else                            !** else write
            CARDS.ABEMLIN=CARD_PARAMS.ABEMLIN_WRITE

          endif
      endif
      !*** 3) open the file
      select case(CARDS.ABEMLIN)
        case(CARD_PARAMS.ABEMLIN_READ)  !** read values from file
          open(201,file=flnam,status='OLD')
          close(200)
        case(CARD_PARAMS.ABEMLIN_WRITE) !** write values to file
          open(201,file=flnam,status='REPLACE',action='write')
      end select

      !*****************************************************************
   
      call synopa(WAVARR, SIGARR, N, NFDIM)
  
   
      PRINT *,'FIOSS: Time elapsed after SYNOPA: ',TOC()
      close (unit=200)
      close (unit=201)
      if(CARDS%PRINT_TAU) then
        write(flnam,'(F14.0)') XLAM
        flnam=adjustl(trim(flnam)//'tau')
        open  (unit=9998,name=flnam)
      endif




      DO L=1,ND









        opamax=0.

        !MH   WAVELENGTH FOR EACH POINT FROM 1 TO NVOPA
        !MH   WLRANGE: WAVELENGTH RANGE OVER HALF THE LINE PROFILE
        !MH   DLAM: DELTA LAMBDA BETWEEN EACH OF NVOPA POINTS   
        WL =0.
        WLRANGE = xlam*fmax*vdop/clight
        xlam0  = xlam-WLRANGE
        dlmbd =    2.*WLRANGE/(nvopa-1)
        do k=1,nvopa 
          WL = xlam0+(k-1)*dlmbd     
          ! write (9996,*) WL,l,opatot(k,l)
          opatot(k,l)=opatot(k,l)*rstar
          if (opatot(k,l).gt.opamax) then
            !if (l.eq.39) write (91,692) k,opatot(k,l)
            opamax=opatot(k,l)
            kmax=k
          endif
        enddo
!test
!         write(91,692) l,opamax
 692     format(i5,1p7e10.2)
        etamax=0.
        do k=1,nvopa
          ETATOT(k,L)=ETATOT(k,L)*rstar
          ! if (etatot(k,l).gt.etamax) etamax=etatot(k,l)
          if (kmax.eq.k) etamax=etatot(k,l)
          !if (l.eq.39) write (92,692) k,etatot(k,l)
        enddo
        !test
        ! write(92,692) l,etamax/opamax
      ENDDO

!***  Compute Integration Weights for Moment0
      CALL WMOM0_F (ND,NP,R,p,WLK)

      !***  For the first loop, the line dependent radiation field
      !***  must be set to the continuum radiation field

      XJ(L, :) = XJC(L)

      !***  Reset the (old) Profile which is used as Loop - Terminator
      !*    to be sure, that the second Loop is executed
      !*    (look for exloop condition)
      PROFN(1:NFOBS)=-99.0

      !***  HERE COMES THE BIG LOOP, WHERE THE LINE RADIATION FIELD WILL BE
      !***  ITERATIVELY COMPUTED

      DO_ITER:DO M=1,MAXITER
        print *, 'ITERATION No. ', m,' of ', maxiter
        !***  Reset frequency dependet mean intensity (moment 0)
        XJK(:,:)=0.0

        !*** Reset the Profile
        PROFILE(:)=0.

        !***  LOOP FOR EACH IMPACT PARAMETER ===========================

         call system('mkdir -p contr')

         write(rwlae_str, '(i8)') int(rwlae)

!         open(250, file='../contr.txt',access='append')
         open(250, file='./contr/'//trim(adjustl(rwlae_str))//'.contr')
         write(250,*), rwlae

!         N_CLV = 100
!         N_CLV = 20
         N_CLV = 10
!         N_CLV = 1
!         N_CLV = NFOBS
         LP = NP - ND

         ALLOCATE(WAV_CLV(N_CLV))
         ALLOCATE(FLUX_CLV(LP, N_CLV))

         DO JP = JFIRST, JLAST

          if (Jfirst.eq.Jlast) print *,'fioss.for: JP= ',jp
          !***  LOOP FOR EACH ANGLE TO THE ROTATION AXIS
          !***  RESET EMINT
          DO LPHI=LPSTA,LPEND

       

            EMINT(1:NFOBS)=0.
            !***  Reset DINT (Stored EMINT for Computation of Moment 0 (COMPXJ))
            DINT(:,:)=0.

            IRAY=ND*(JP-1)+1
            IF (NPHI.GT.1) PHI=PI*(LPHI-1)/(NPHI-1)

            CALL PREPR_F(Z2D(1 : ND, JP),P,ND,NP,JP,LTOT,LMAX,WE,CORE,VDU,R,
     $                   IRIND,IBACK,RRAY,ZRAY,XCMF,NDADDIM)

            CALL OBSINT(LTOT,CORE,BCORE,DBDR,P,
     $                  IRIND,RRAY,ZRAY,XCMF,
     $                  ND,NP,NVD,JP,NVOPA,VOPA0,DVOPA,
     $                  EMINT,XOBS0,DXOBS,NFOBS,XN,
     $                  ENTOT,RNE,SIGMAE,RSTAR,
     $                  XJ,DINT,RWLAE,DLAM)

            !***  Compute the Line Intensity Field (Moment 0)
            !* add into array XJK
            if (m.lt.maxiter) CALL COMPXJ9(ND,NP,JP,NVopa,DINT,XJK,WLK,IBACK)

            !***  Compute the Profile by adding the new Intensities per JP
            !     add into array PROFILE

            IF (NPHI .GT. 1) THEN

                WPHI = 1. / (NPHI - 1)
                IF (LPHI.EQ.1 .OR. LPHI.EQ.NPHI) WPHI=WPHI/2.
                WPHI=WPHI*WE
                CALL COMPRO (PROFILE,EMINT,NFOBS,WPHI,JFIRST,JLAST)

            ELSE

                CALL COMPRO (PROFILE,EMINT,NFOBS,WE,JFIRST,JLAST)

            ENDIF

          ENDDO ! LPHI

!         AVERAGING THE SPECTRUM FOR CLV CALCULATIONS

          IF (JP .LE. LP) THEN

             IS = INT(NFOBS / N_CLV)
             FP = INT(IS / 2)

             DO K = 1, N_CLV

                IF (N_CLV .GT. 1 .AND. N_CLV .LT. NFOBS) THEN

                   MP = FP + IS * (K - 1)
                   SP = IS * (K - 1)
                   EP = IS * K

                ELSE IF (N_CLV .EQ. 1) THEN

                   MP = FP
                   SP = 1
                   EP = NFOBS
                
                ELSE IF (N_CLV .EQ. NFOBS) THEN

                   MP = K
                   SP = K
                   EP = K

                ELSE

                   PRINT*, 'N_CLV = ', N_CLV, ' IS NOT RECOGNIZED. ABORT.'; STOP

                ENDIF

                WAV_CLV(K) = DLAM(MP) + RWLAE

                if (wav_clv(k) <= 290.0d0) FLUX_CLV(JP, K) = 0.0d0

                if (wav_clv(k) > 290.0d0)  FLUX_CLV(JP, K) = SUM(emint(SP : EP)) / IS

             ENDDO

          ENDIF

        ENDDO ! LOOP OVER JP (IMPACT PARAMETERS)

        close(250)

        !***  THE EXTREME FREQUENCY POINTS DEFINE THE REFERENCE CONTINUUM
        print *
        PRINT *,"Loop-Nr: ",M
        print '(A,1p5e12.3)','profile(1) ... (3)',(profile(k),k=1,3)
        REFCON=0.5*(PROFILE(1)+PROFILE(NFOBS))

        prof1=PROFILE(1)       ;   prof2=PROFILE(NFOBS)
        do K=1,NFOBS
          if (norm) then
              CON=((NFOBS-K)*PROF1+(K-1)*PROF2)/real(NFOBS-1)
              print *,' K=',K,PROFILE(K)
              PROFILE(K)=PROFILE(K)/CON






           endif
        enddo
        print '(A,1pe12.4,A,e12.4)',
     &      ' Profile Reference continuum: ',refcon,' Elimin:',fnucont

        !***  Compare the old profile with the new One
        !***    end the DO LOOP if they agree
        EXLOOP = all(PROFILE(1:NFOBS)-PROFN(1:NFOBS) <= 5e-3 )
        if (EXLOOP) exit DO_ITER

        !***  Copy the new profile for the next comparison
        PROFN(1:NFOBS)=PROFILE(1:NFOBS)

        !***  Compute the new mean intensity XJ for use in Obsint
        !***  Electron-Scattering Computation and fold with thermal
        !***  Gauss-Profile
        if (m.lt.maxiter) CALL COMPGAU9(XJ,XJK,CWK,XNU,T,DSTEP,ND,NVOPA)

      ENDDO DO_ITER
!***  200: END OF el-sca MAXITER LOOP ----------------------------------------------

      IF (EXLOOP) THEN
        PRINT *,'Abbruch bei Konvergenz'
        PRINT *,'**********************'
      ELSE
        PRINT *,'Abbruch nach Max. Anzahl Iterationen'
        PRINT *,'************************************'
      ENDIF
!***  OUTPUT FOR THE DETECTED LINE:

      LOW=INDLOW(LINE)
      NUP=INDNUP(LINE)

      IF (PLOT)
     $   CALL LPPLOT (6,XOBS0,DXOBS,NFOBS,PROFILE,KARTE,MODHEAD,JOBNUM)

      IF ((TRANS.OR.PROLIB).AND.LSDWL.LE.1) then
         IF (shift.gt.0) then
            RWLAE=XLAM*(1.+shift/3.d5)
            print '(A,F10.0,A)',' Wavelength shift by ',shift,' [km/s]'
         endif

         ! write the <wl>.mdisp & <wl>.title file
         
         CALL TRAPLO(PROFILE,DLAM,NFOBS,KARTE,MODHEAD,JOBNUM,RWLAE,PHEAD,PROLIB)

      endif
 

    !   stop

      IF (LSDWL.LE.1)
     $CALL       PRIPRO (XLAM,VDOP,NFOBS,PROFILE,XOBS0,DXOBS,JOBNUM,
     $     VSINI,REFCON,MODHEAD,DLAM,LSPRO,IFIRST,NPHI,LPSTA,LPEND,
     $          XN,XMAX,JFIRST,JLAST,P,WE,PHEAD,PROLIB,
     $          KARTE,WEIGHT(LOW),WEIGHT(NUP),LEVEL(LOW),LEVEL(NUP),
     $          EINST(NUP,LOW),FNUEC,RSTAR2,VELO(1)
     $         ,EQWI)
      ENDDO MAIN_LOOP
      !GOTO 1
!***  ENDLOOP    -------------------------------------------------------
!   20 CONTINUE
      PRINT *,' END OF INPUT REACHED'
      ! CLOSE (UNIT=880) 
      ! CLOSE (UNIT=890) 
      ! CLOSE (UNIT=870)
      ! close (unit=9996)
      ! close (unit=9997)
      if(CARDS%PRINT_TAU) close (unit=9998)
      IF (TRANS.OR.PROLIB) THEN
        CLOSE (1)
        ! CALL JSYMSET (2LG1,'TRANSFER')
        ! CALL REMARK ('PLOT DATA TO BE ROUTED')
      ENDIF
 
      WRITE(CLVFLNAM, '(F14.0)') XLAM; CLVFLNAM = ADJUSTL(TRIM(CLVFLNAM)//'clv')

      CALL OPEN_TO_APPEND(100, CLVFLNAM)

      DO K = 1, N_CLV

          WRITE(100, '(E15.7,1x,$)'), WAV_CLV(K)

          DO JP = JFIRST, LP

             IF (JP .LT. LP) WRITE(100, '(E15.7,1x,$)'), FLUX_CLV(JP, K)
             IF (JP .EQ. LP) WRITE(100, '(E15.7)'),      FLUX_CLV(JP, K)

          ENDDO

      ENDDO

      CLOSE(100)
 
      STOP 'O.K.'

      END
