      module MOD_INIBL0
      character*(*),private,parameter :: FMT_LOPA='(1pe12.5)'   ! format 310
      real*8, allocatable ::wav_oi(:), opac_oi(:), opac_f(:,:)
      contains

       SUBROUTINE inibl0(WAVARR, SIGARR, N, NF)
       use MOD_SYNSUBM
       use SYNTHP_CONT,only: FREQC,ABSOC,NFCONT
       use UTILS,only:assert
       use MOD_DECF_SYN
       use CONSTANTS,only: CLIGHT_SI
       use MOD_BALINI
       use MOD_HYDTAB
       use MOD_chemeq
C       use MOD_LINOP_MS
C
C     AUXILIARY INITIALIZATION PROCEDURE
C
      implicit none
      integer,intent(in) :: N,      NF
      real*8, intent(in) :: WAVARR, SIGARR

      integer :: IFHE2,IHE1,IHE144,IHE2UV,IHE2RE,IHE2VI,IHE2L,ILWHE2
      integer :: INLIST,NLTOFF,IEMOFF,I,I1,I2,ID,ILVCS,IBVCS,MHE10,MHE20
      real*8  :: ALAST,CUTOF0,CUTOFS,RELOP,SPACE,VELMAX
      real*8  :: FRLAST,XX,ACOR,WPH,A1,A2,A3
      real*8  :: ALAM0,ALAM1,FRMIN,FRLI0,FRLIM,TSTD,DSTD

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'

      real*8,parameter :: un=1.
      logical lwph
      real*8,dimension(MFREQ) :: ABSO,EMIS !,SCAT(2)

      DIMENSION WAVARR(N, NF),SIGARR(N, NF)
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE,CUTOF0,CUTOFS,TSTD,DSTD
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
      COMMON/HE2PAR/IFHE2,IHE2L,ILWHE2,MHE10,MHE20
      common/wprob/wph(mdepth,40),acor(mdepth),lwph
      common/linoff/velmax,nltoff,iemoff
      !*** DEBUG
      integer :: NWLHYD(MLINH)
      real*8  :: WLHYD(MLINH,MHWL)
      real*8  :: PRFHYD(MLINH,MDEPTH,MHWL)
      integer :: j
      COMMON/HYDVCS/PRFHYD,WLHYD,NWLHYD
      !*** DEBUG END
C
C --------------------------------------------------------------
C Parameters controlling an evaluation of the synthetic spectrum
C --------------------------------------------------------------
C
C     ALAM0, ALAM1 - synthetic spectrum is evaluated between wavelengths
C                    ALAM0 (initial) and ALAM1 (final), given in Anstroms
C     CUTOF0       - cutoff parameter for normal lines (given in Angstroms)
C                    ie the maximum distance from the line center, in
C                    which the opacity in the line is allowd to contribute
C                    to the total opacity (recommended 5 - 10)
C     RELOP        - the minimum value of the ratio (opacity in the line
C                    center)/(opacity in continuum), for which is the line
C                    taken into account (usually 1d-4 to 1d-3)
C     SPACE        - the maximum distance of two neighbouring frequency
C                    points fro evaluating the spectrum; in Angstroms
C
C     INLTE        = 0  -  pure LTE (no line in NLTE)
C                  ne.0 -  NLTE option, ie one or more lines treated
C                          in the exact or approximate NLTE approach
C                          in this case, other input parameters have to
C                          be specified - unit 10 (see procedure NLTE)
C     ICONTL       = 1  -  Lyman and Balmer lines are considered as an
C                          continuum opacity source
C                  ne.1 -  Lyman and Balmer lines are not considered as
C                          an continuum opacity source
C     INLIST       = 1  -  line list is read in the "original" format
C                          (see procedure INISET)
C                  ne.1 -  line list is read in the "new" format
C                          (see procedure INISET)
C     IFHE2        gt.0 -  He II line opacity in the first four series
C                          (Lyman, Balmer, Paschen, Brackett)
C                          for lines with lambda < 3900 A
C                          is taken into account even if line list
C                          does not contain any He II lines (i.e.
C                          He II lines are treated as the hydrogen lines)
C
C     and, finally, parameters that were not present in the previous
C     version, namely ILVCS,IBVCS,IHE1,IHE447,IHE2UV,IHE2VI,IHE2RE.
C
C     IBVCS       = 0  - means that Balmer lines are calculated by
C                        approxiamte formulae
C                 > 0  - means that H-alpha to H-delta are calculated
C                        in detail, using the Vidal, Cooper, Smith tables;
C                        the tables are stored in file FOR0xx.dat,
C                        where xx=IBVCS;
C                        higher Balmer lines are calculated as before
C
C     the meaning of other parameters is quite analogous, for the
C     following lines
C
C     IHE1      - He I lines at 4026, 4387, and 4922 Angstroms
C                 (tables calculated by L. Shamey)
C     IHE144    - He I line at 4471 angstroms; tables calculated by
C                 Barnard, Cooper, and Shamey
C
C     IHE2UV    - for the He II lines calculated by Schoening and Butler,
C                 whose wavelengths are in the UV part of the spectrum,
C                 ie. the lines 2-3, 3-5, 3-6, 3-7, 3-8, 3-9, and 3-10.
C     IHE2VI    - the same for the "visible" lines, ie. 3-4, 4-8, 4-9,
C                 4-10, 4-11, 4-12, 4-13, 4-14, and 4-15.
C     IHE2RE    - the same for IR and "red" lines, ie. 4-5, 4-6, and 4-7.
C
      velw(:) = 0d0   !*** set velw to zero. It is not set anywhere else, 
                      !*** but used! (velw < velmax is always true) --micha
      READ(55,*) INLTE,INLIST,IFHE2  ! Read in the options from the (cryptic) file fort.55
      READ(55,*) ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
      READ(55,*) ALAM0,ALAST,CUTOF0,RELOP,SPACE
cc    READ(81,*) INLTE,INLIST,IFHE2
cc    READ(81,*) ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
cc    READ(81,*) ALAM0,ALAST,CUTOF0,RELOP,SPACE
      VELMAX=3.e10
      NLTOFF=0
      IEMOFF=0
cc    READ(81,*,END=5,ERR=5) VELMAX,NLTOFF,IEMOFF
      READ(55,*,END=5,ERR=5) VELMAX,NLTOFF,IEMOFF
    5 write(6,602) velmax*1.e-5,nltoff,iemoff
  602 format(//' velmax(velocity for line rejection),nltoff,iemoff',
     *  f10.1,2i3)
cccc  CLOSE (UNIT=81)
C

      IF(IMODE.EQ.-1) THEN
         INLTE=0
         CUTOF0=0.
      END IF
!alex ALAMC=(ALAM0+ALAST)*0.5
!      ALAMC=(WLAM(1)+WLAM(NFREQ))*0.5
!      SPACF=2.997925E18/ALAMC/ALAMC*SPACE
!      SPACL=SPACE
!alex use predefined wavelength grid
      ALAST = WLAM(1)
      ALAM0 = WLAM(NFREQ)
      !WRITE(6,601) ALAM0,ALAST,CUTOF0,RELOP !,SPACF,SPACL
      CUTOF0=0.1*CUTOF0
      SPACE=SPACE*0.1
!alex use predefined wavelength grid
      ALAM0=1.D-1*ALAM0
      ALAST=1.D-1*ALAST
!      =2.997925D17/ALAM0
      FRLAST=CLIGHT_SI*1d9/ALAST
c     if(nfreq.gt.mfreq) nfreq=mfreq
c
c     frequency points
c
!     do 10 ij=1,nfreq
!calex    freq(ij)=frfirs-(ij-1)*spacf
!        FREQ(IJ) = FOBS(IJ)
!alex    WRITE (6,*) IJ,WOBS(IJ),FOBS(IJ)
!  10 continue
C
C     frequency points for evaluating the contiuum opacity
C
      !nfcont=2
cdis  freqc(1)=freq(1)
c     freqc(1)=2.997925D18/wmin
c     wobs(nobs+1) = wmin
cdis  freqc(2)=freq(nfreq)
C     freqc(2)=2.997925D18/wmax
c     wobs(nobs+2) = wmax
C
C     INTERPOLATION COEFFICIENTS (FOR FREQUENCY INTERPOLATION)
C
      !*** calculate the weights for integration
      i1=1
      i2=1
!       if(allocated(FRXIDX)) deallocate(FRXIDX)
!       allocate(FRXIDX(NFCONT))
      FRXIDX(1)=1
      DO i=1,NFCONT()-1
        XX=FREQC(i+1)-FREQC(i)
        call assert(XX > 0.,'XX<=0')
        i1=i2
        do while(i2<NFREQ .and. FREQ(i2) < FREQC(i+1))
          i2=i2+1
        enddo
        FRXIDX(i+1)=i2
        FRX1(i1:i2)=(FREQ(i1:i2)-FREQC(i))/XX
        FRX2(i1:i2)=(FREQC(i+1)-FREQ(i1:i2))/XX
      ENDDO
      WLAM(1:NFREQ)=CLIGHT_SI*1d10/FREQ(1:NFREQ)
      ! XX=FREQC(2)-FREQC(1)
      ! FRX1(1:NFREQ)=(FREQ(1:NFREQ)-FREQC(1))/XX
      ! FRX2(1:NFREQ)=(FREQC(2)-FREQ(1:NFREQ))/XX

C  the following 4 statements commeneted out for fioss:
c      NFREQ=2
c      FREQ(1)=2.997925D17/ALAM0
c      FREQ(2)=FRLAST
c      CALL SIGAVS
      IF(IBVCS.GT.0)  THEN
        IF(IBVCS/=24) THEN
          print *,'call HYDINI'
          CALL HYDINI(IBVCS)
        ELSE
          print *,'call BALINI'
          CALL BALINI(IBVCS)
        ENDIF
      ENDIF
      IF(IHE144.GT.0) CALL HE1INI(1,IHE144)
      IF(IHE1.GT.0)   CALL HE1INI(2,IHE1)
      IF(IHE2UV.GT.0) CALL HE2INI(1,IHE2UV)
      IF(IHE2VI.GT.0) CALL HE2INI(2,IHE2VI)
      IF(IHE2RE.GT.0) CALL HE2INI(3,IHE2RE)
c
c   set up also occupation probabilities for the first 40 levels of
c   hydrogen
c
      do id=1,nd
        if(lwph) then
          acor(id)=0.09*elec(id)**.166666667/sqrt(temp(id))
          do i=1,40
              wph(id,i)=wn(dble(i),acor(id),elec(id),un)
          enddo
        else
          wph(id,1:40)=1.
        endif
      enddo
C
c     pretabulate expansion coefficients for the Voigt function
c
      CALL PRETAB
c
c     calculate the characteristic standard opacity
c
      IF(IMODE.LT.2) THEN
cmh      CALL CROSET
C***  CHANGES BY MARGIT HABERREITER ***
         CALL CROSET(WAVARR,SIGARR, N, NF)
         print*, 'CROSET'
 
          call readmollines
          print*, 'readmollines'
 
          call readMolconc(ND)
          print*, 'readMolconc'

       
    

         DO ID=1,ND

       
            CALL OPAC(ID,0,ABSO,EMIS,WAVARR,SIGARR,N,NF)

        

            ABSTD(ID)=MINVAL(ABSOC(:NFCONT()))
         ENDDO

   
       

   

         IF(INLTE.GT.0) CALL NLTE(0,1,1,1,A1,A2,A3)
         IF(cards.ABEMLIN/=card_params.ABEMLIN_READ)
     *    CALL INILIN(INLIST)

      END IF


  502 FORMAT(16I5)
  503 FORMAT(4F10.3,2E10.3)
  601 FORMAT('1----------------------------------------------'/
     *       ' BASIC INPUT PARAMETERS FOR SYNTHETIC SPECTRA'/
     *       ' ---------------------------------------------'/
     *       ' INITIAL LAMBDA',28X,'=',F10.3,' ANGSTROMS'/
     *       ' FINAL   LAMBDA',28X,'=',F10.3,' ANGSTROMS'/
     *       ' CUTOFF PARAMETER',26X,'=',F10.3,' ANGSTROMS'/
     *       ' MINIMUM VALUE OF (LINE OPAC.)/(CONT.OPAC) =',1PE10.1/
     *       ' MAXIMUM FREQUENCY SPACING',17X,'=',1PE10.3,'  I.E.',
     *          0PF6.3,'  ANGSTROMS'/
     *       ' ---------------------------------------------'/)
      RETURN
      END subroutine
C
C ********************************************************************
C ********************************************************************
C


      SUBROUTINE OPAC(ID, MODE, ABSO, EMIS, WAVARR, SIGARR, N, NF)
C     ========================================
C
C     Absorption, emission, and scattering coefficients
C     at depth ID and for several frequencies (some or all)
C
C     Input: ID    - depth index
C            CROSS - two dimensional array of photoionization
C                    cross-sections
C     Output: ABSO - array of absorption coefficient
C             EMIS - array of emission coefficient
C             SCAT - array of scattering coefficient (all scattering
C                    mechanisms except electron scattering)
C
C             mode = 0 - only continuum opacities
C             mode = 1 - continuum + line opacities
C
CMH   SIGEL: electron scattering in cm2/electron
CMH   k           = Boltzmann-Konstante
CMH   CBF         = (1/2) *  (h2/(2 pi m_e k^(3/2) ;
CMH                 a constant in a Saha-Boltzmann formula (see e.g. Mihalas 1978, p.113; called C_I there)
CMH   CFF         = (4 e6/3 c h) * (2 pi/3 k m_e3)^(1/2)   (see e.g. Mihalas 1978, pp. 101-102)
CMH   EBF:  = emission bound free?
CMH   HK          = 4.79928144D-11 (Planck/Boltzmann = h/k in K*s)
CMH   HKT         = h/k/T
CMH   HINV  = 1/h
CMH   HKF         = h*nu/k/T
CMH   KT          = 1/k/T
CMH   X           = EXP(-h*nu/k/T)
CMH   X1          = 1-EXP(-h*nu/k/T)
CMH   BN          = 2*h/c^2 ?????(units are not correct)
CMH   BN      = C2 = 2 * H * C = 8.2779372e-006  (DIMENSION ANGSTROEM**3  * ERG/SEC/HZ/CM**2)
CMH   BN          = 1.4388e-2 m*K (Unsoeld und Baschek, Planck-Strahlungskonstante c2)
CMH   BNU         = BN*(nu^3)
CMH   CSB         = (1/2) *  (h2/(2 pi m_e k))^(3/2)
CMH   CFF         = (4 e6/3 c h) * (2 pi/3 k m_e3)^(1/2)
C
      use MOD_HMINUSFF
      use MOD_LINOP_MS
      use MOD_SYNSUBM
      use MOD_DECF_SYN
      use UTILS
      use SYNTHP_CONT, only: ABSOC,EMISC,FREQC,SCATC,NFCONT
      use constants, only:   CLIGHT_SI, CLIGHT_CGS, BoltzmannConstantEV,
     $                       PlanckConstantEV
      use MOD_chemeq
      use MOD_BNUE

!*****************************************************************************************************************
!Module introduced by Rinat Tagirov for proper calculation of Gaunt factors, see the module itself for the details
      USE MOD_GFF_TEMP

!Module introduced by Rinat Tagirov for simulation of hydrogen ionization step
      USE MOD_HYD_ION_STEP
!*****************************************************************************************************************

      implicit none

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/LINDAT.FOR'
      
      integer,intent(in) :: id, NF, N
      integer,intent(in) :: MODE
      real*8, intent(in)   ,dimension(N, NF)  :: WAVARR,SIGARR
      real*8, intent(inout),dimension(MFREQ) :: ABSO,EMIS

      !*** DEBUG
      integer :: ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
      !*** DEBUG END
      COMMON/PHOPAR/CROSS(MCROSS,MFCONT)
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
     $ IOPFE1
      COMMON/HYLPAR/IHYL,ILOWH,M10,M20
      COMMON/HE2PAR/IFHE2,IHE2L,ILWHE2,MHE10,MHE20
      COMMON/BLAPAR/RELOP,SPACE,CUTOF0,CUTOFS,TSTD,DSTD
      common/wprob/wph(mdepth,40),acor(mdepth),lwph
      real*8,parameter :: frh = 3.28805d+15
      real*8,parameter :: UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.68E8
      real*8,parameter :: HINV=1.5092973D26
      integer :: N0I,N1I,NKE, indo
      integer :: I,II,I1,I2,IE,IL,IJ,IQ,IT, ILOWH
      integer :: IHYL,IHE2L,IFHE2,ILWHE2, IRSCT
      integer :: M10,M20,MHE10,MHE20
      integer :: IOPADD,IOPHMI,IOPH2P,IOPHLI,IOPHE1,IOPHE2,IOPFE1
      integer :: ID_TMP, NFREQ_TMP

      real*8  :: ABAD,EMAD,SCAD,X1,X
      real*8  :: T,T1,TK,ANE,SRT,SGFF,CON,CONTS,FR,FR15,BNU,HKF,HKT
      real*8  :: ABF,EBF,AFF,ABLY,EMLY,SCLY,XN,SG,XW,WNSTAR,XNSTAR
      real*8  :: ACOR,DW,WPH,CROSS,XX,XLAM,SIGHMFF,FACTOR,CH
      real*8  :: HKFM,SFF,SF1,SF2,BNE,RELOP,CUTOF0,CUTOFS
      real*8  :: SPACE,TSTD,DSTD
      logical :: lwph
      real*8,dimension(MFREQ) :: ABLIN,EMLIN
      real*8 molopac(NFREQ), molemiss(NFREQ)
      real*8 :: freqt, lambdat, contf, totFe, totFeI, totFeII, totH
      integer :: ind(1)

!******************************************************
!RINAT TAGIROV
      REAL*8  :: AbsCoefH, AbsCoefHM, AbsCoefRest

      REAL*8  :: HydIonDeg, HeIonDeg

!      REAL*8  :: SG_WRONG

      INTEGER :: UnitRatios

      REAL*8 ::  RI, f_0
!******************************************************

!******************************************************************************************************************************************
!RINAT TAGIROV
      UnitRatios = 10506
      OPEN(UNIT=UnitRatios, FILE="absorption_coefratio.out",
     $ ACTION="write", ACCESS="append") ! file for printing out the ratios of absorption coefficients
!2001  FORMAT('wavelength',10x,'height_index',10x,'temperature',10x,
!     $ 'n_e',10x,'HM/(H+HM)',10x,'rest/(H+HM)')
2000  FORMAT(E15.7,8x,I2,8x,E15.7,8x,E15.7,8x,E15.7,8x,
     $ E15.7,8x,E15.7,8x,E15.7,8x,E15.7,8x,E15.7)
!      WRITE(UnitRatios, 2001)
!******************************************************************************************************************************************

CMH   IELHM SET TO 1, FIRST ELEMENT IS HMINUS
CMH   PRINT *,'SUBROUTINE OPAC : IELHM SET TO ZERO'
      IELHM = 1
C     print *,'opad......'
      T=TEMP(ID)
      ANE=ELEC(ID)
      T1=UN/T
      HKT=HK*T1
      TK=HKT*HINV
      SRT=UN/SQRT(T)
      SGFF=CFF*SRT
      CON=CSB*T1*SRT
      conts=1.e-36/con
C
C     Opacity and emissivity in continuum
C     **** calculated only in the first and the last frequency *****

      FCONT_LOOP: DO IJ=1,NFCONT()
        FR=FREQC(IJ)
        FR15=FR*TEN15
        BNU=BN*FR15*FR15*FR15
        HKF=HKT*FR
        ABF=0.
        EBF=0.
        AFF=0.
        ABLY=0.
        EMLY=0.
        SCLY=0.

!******************************************************************************************************************************************
!RINAT TAGIROV
        AbsCoefRest = 0.

!        CALL HYD_ION_STEP(ID, 65, POPUL(1:12, ID), ANE)
!******************************************************************************************************************************************

        NION_LOOP: DO IL=1,NION
          N0I=NFIRST(IL)
          N1I=NLAST(IL)
          NKE=NNEXT(IL)
          XN=POPUL(NKE,ID)

          ! Bound-free contribution
          ! + pseudo-continuum for H (accounting for disolved fraction)

          N0I_LOOP: DO II=N0I,N1I
            IF(IBF(II).EQ.0) cycle N0I_LOOP
            sg=0.
            IF(IFWOP(II).LT.0) THEN
              sg=0.
              ! SG=SGMERG(II,ID,FR,DSG)
            ELSE
              If(lwph .and. il.eq.ielh) Then
                iq=nquant(ii)
                xw=1./(iq*iq)-fr/frh
                if(xw.gt.0.) then
                    xnstar=1./sqrt(xw)
                    if(xnstar.gt.iq+1) then
                      wnstar=wn(xnstar,acor(id),ane,un)
                      dw=(wph(id,iq)-wnstar)/wph(id,iq)
                    else
                      dw=0.
                    end if
                else
                    dw=1.
                endif
                if(dw.gt.0.) then
                   sg=sigk(fr, ii, 1, WAVARR, SIGARR, N, NF)*dw
                endif
              Else
                SG=CROSS(II,IJ)
              Endif
            END IF
            ABF=ABF+SG*POPUL(II,ID)
            XX=SG*XN*EXP(ENION(II)*TK)
            IF(XX.ge.conts) EBF=EBF+XX*CON*G(II)/G(NKE)
          ENDDO N0I_LOOP
          IT=IFREE(IL)
          IF(IT.EQ.0) cycle NION_LOOP

          !*** Free-free contribution

          IE = IL

!******************************************************************************************************************************************
!RINAT TAGIROV
          xlam = CLIGHT_SI*1d10/fr ! This variable has been pulled out of the following IF condition because it was needed for output below 
!******************************************************************************************************************************************

          IF(IE.EQ.IELHM) THEN
            !*********************************************************
            !***  MARGIT HABERREITER
            !MH   USE SUBROUTINE HMINUSFF TO CALCULATE HMINUS OPACITY
            !MH   units sighmff: cm^4/dyne
            !***  xlam in Angstom
            !***  cl speed of light in cm/s
            !***  factor = k_Boltzmann*T*electron density
            !***  xlam: wavelength in Angstrom for hminusff subroutine
            !MH  ORIGINALLY:
            !   65        SFF=SFFHMI(XN,FR,T)
!            xlam = CLIGHT_SI*1d10/fr ! This variable has been pulled out of this IF condition because it was needed for output below
            call hminusff(sighmff,xlam,T)
            factor = 1.38062259d-16*T*XN
            sff = sighmff*factor
!***********************************************************************************
!RINAT TAGIROV
            IF (IJ .EQ. 1) THEN
               AbsCoefHM = ANE * SFF
            ENDIF
!***********************************************************************************

          ELSE  !if IE /= IELHM
            CH=IZ(IL)*IZ(IL)
            SF1=CH*XN*SGFF/(FR*FR*FR)

            ! The following expression is the so-called modified free-free
            ! opacity, ie. allowing for the photoionization from higher,
            ! non-explicit, LTE energy levels of the ion IL

            HKFM=HKT*MIN(FF(IL),FR)
            SF2=EXP(HKFM)
            if(IT.EQ.2) then
!              SG=GFREE(T,FR/CH)

!RINAT TAGIROV:
!GFF_TEMP is the gaunt factor for any Rydberg system as a function of temperature and frequency
!See the subroutine for more details

!              PRINT*, 'ACHTUNG: FR =' , FR

              CALL GFF_TEMP(1, 1.0D0, PlanckConstantEV * FR,
     $                      BoltzmannConstantEV * T, SG)

!              SG = SG * 15.0D0

              SF2=SF2+SG-UN
            endif

            X=EXP(-HKF)
            X1=UN-X
!***********************************************************************************
!RINAT TAGIROV, ALEXANDER SHAPIRO
            SFF=X1*SF1*SF2 !!! X1 has been added to account for the induced emission for every element except H^- (the induced emission for H^- was taken into account in hminusff)
!***********************************************************************************
          ENDIF ! IE==IELHM

          AFF=AFF+SFF

!***********************************************************************************
!RINAT TAGIROV
          IF ((IL .EQ. 2) .AND. (IJ .EQ. 1)) THEN
             AbsCoefH = ANE * SFF
          ENDIF

          IF ((IL .NE. 1) .AND. (IL .NE. 2) .AND. (IJ .EQ. 1)) THEN ! Overall absorption coefficient due to all elements except H and H^-
!          IF ((IL .EQ. 4) .AND. (IJ .EQ. 1)) THEN                  ! HeI absorbtion coefficient
             AbsCoefRest = AbsCoefRest + ANE * SFF
          ENDIF
!***********************************************************************************

        ENDDO NION_LOOP

!**************************************************************************************************
!Introduction of the plasma refractive index by Rinat Tagirov (see Tapping & DeTracey, 1990, SoPh)

!        f_0 = 9.0D3 * DSQRT(ANE)

!        IF (f_0 / FR .GE. 1.0D0) THEN

!           AFF = 1.0D40

!        ELSE

!        PRINT*, 'AAAAAAAAAAAAAAAA', f_0, FR, ANE

!           RI = DSQRT(1.0D0 - (f_0 / FR)**2.0D0)

!        PRINT*, 'AAAAAAAAAAAAAAA2', RI

!           AFF = AFF / RI

!        ENDIF

!**************************************************************************************************

        !* Additional opacities

        CALL OPADD(0,IJ,ID,FR,ABAD,EMAD,SCAD)
        ! IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY)

        !* Total continuum opacity and emissivity
        !*
        BNE=BNU*X*ANE
CMH       ABSOC(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ANE*SIGEL+ABAD+ABLY
CMH         HMINUSFF fit includes stimulated emission
cmh         correction by X1 = 1. - exp(-h*nu/k*T) obsolete for Hminus
        ABSOC(IJ)=ABF+ANE*(AFF-X*EBF)+ANE*SIGEL+ABAD+ABLY
        EMISC(IJ)=BNE*(AFF/X1+EBF)+EMAD+EMLY
        SCATC(IJ)=SCAD+SCLY+ANE*SIGEL
        if (ij .eq. 5) then 
   !     write(*, 615) ID, ANE*SIGEL,SCAD, ABSOC(IJ)+SCATC(IJ)
        endif

! 615    FORMAT('239  ', I2, 3E12.4)

!***********************************************************************************
!RINAT TAGIROV
        IF ((IJ .EQ. 1) .AND. (MODE .EQ. 0)) THEN

           HydIonDeg = POPUL(12, ID) / SUM(POPUL(1:12, ID))

           HeIonDeg =  POPUL(23, ID) / SUM(POPUL(13:23, ID))

!GFF_TEMP is the gaunt factor for any Rydberg system as a function of temperature and frequency
!See the subroutine for more details

           CALL GFF_TEMP(1, 1.0D0, PlanckConstantEV * FR,
     $                  BoltzmannConstantEV * T, SG)

           WRITE(UnitRatios, 2000), xlam, ID, T, ANE,
     $     AbsCoefHM * 1.0D+2 / (AbsCoefH + AbsCoefHM),
     $     AbsCoefRest * 1.0D+2 / (AbsCoefH + AbsCoefHM),
     $     GFREE(T, FR / CH), SG, HydIonDeg * 1.0D+2, HeIonDeg * 1.0D+2

        ENDIF
!***********************************************************************************

      ENDDO FCONT_LOOP

!***********************************************************************************
!RINAT TAGIROV
      CLOSE(UnitRatios)
!***********************************************************************************

      if(mode.eq.0) return

      !AVAB=(ABSOC(1)+ABSOC(NFCONT))*0.5*RELOP
      !***************************************************************
      !***  interpolated continuum and hydrogen line opacity and
      !***  emissivity for all frequencies
      i2=0
      do i=1,NFCONT()-1
        i1=i2+1
        i2=FRXIDX(i+1)
        ABSO(i1:i2)=FRX1(i1:i2)*ABSOC(i+1)+FRX2(i1:i2)*ABSOC(i)
        EMIS(i1:i2)=FRX1(i1:i2)*EMISC(i+1)+FRX2(i1:i2)*EMISC(i)
      enddo

!       ABSO(1:NFREQ)=FRX1(1:NFREQ)*ABSOC(2)+FRX2(1:NFREQ)*ABSOC(1)
!       EMIS(1:NFREQ)=FRX1(1:NFREQ)*EMISC(2)+FRX2(1:NFREQ)*EMISC(1)
      !
      !**** Opacity and emissivity in lines ****
      !
      !*** If the ABEMLIN key is set to READ then read the precalculated
      !*** ABLIN/EMLIN from files, otherwise calculate and write out the
      !*** lopa and if the ABEMLIN key is set to WRITE then write out
      !*** the ABLIN/EMLIN variables.
      IF(cards.ABEMLIN==card_params.ABEMLIN_READ) THEN
        read(201,*) NFREQ_TMP,ID_TMP
        call assert(ID_TMP==ID,'opac:ABEMLIN R:ID /= FILE ID')
        call assert(NFREQ==NFREQ_TMP,'opac:ABEMLIN R:NFREQ/=FILE NFREQ')
        read(201,*) ABLIN(1:NFREQ)
        read(201,*) EMLIN(1:NFREQ)
      ELSE

!        PRINT*, 'WE ARE HERE'

        CALL LINOP_MS(ID,ABLIN,EMLIN)
 !       print*, emlin

!          PRINT*, 'LA 1'
!          write(*, *) NFREQ, ID
!          PRINT*, 'LA 2'
!          write(*, *) ABLIN(1 : NFREQ)
!          PRINT*, 'LA 3'
!          write(*, *) EMLIN(1 : NFREQ)



        if(cards.ABEMLIN==card_params.ABEMLIN_WRITE) then
          write(201,*) NFREQ,ID
          write(201,*) ABLIN(1:NFREQ)
          write(201,*) EMLIN(1:NFREQ)
        endif
      ENDIF  ! EMABLIN /= EMABLIN_READ
      ! CALL LINOP(ID,ABLIN,EMLIN,AVAB)
      ! print *,(abso(ij),emlin(ipj),ij=1,nfreq)
      !***********************************************************
      !***  MARGIT HABERREITER
      !***  WRITING OUTPUT FOR NON-LTE BLANKETING
      !***  WHICH WILL BE THE INBUT FOR HMINUS-ROUTINE
      !***  ID: depth point
      !***  NFREQ: number of frequency points
      write (200,*) nfreq, id
    

! *****************************
       freqt=(freq(1)+freq(NFREQ))/2.

       lambdat=(clight_cgs/freqt)*1.d8

       contf=1.

!====================================================================
!FUDGE REGULAR

       if ((lambdat .lt. 3200.) .and. (lambdat .gt. 1600.)) then
        
       ind=minloc(abs(wav_f(1:Nfudge)-lambdat))  

       contf=ffactor(ind(1))
  
       endif 
!====================================================================
!FUDGE BIG

!       if ((lambdat .lt. 3600.) .and. (lambdat .gt. 1600.)) then
        
!       ind=minloc(abs(wav_f(1:Nfudge)-lambdat))  

!       contf = 1.0D6
  
!       endif 

!===================================================================

       totFe=popul(101,id)+popul(102,id)+popul(103,id)+popul(104,id)+popul(105,id)+popul(106,id)   
       totFeII=popul(106, id)
       totFeI=totFe-totFeII

       totH=sum(popul(1:12,id))

!        open(unit=400, file='conc.txt', access='append')
!       write(400, '(i3, 4e12.5)'), id, molconc(4, id)/totH, molconc(3, id)/totH, totFeI/totH, totFeII/totH 
!        close (unit=400)
!       print*, "CO test", id, T, 20161.*molconc(4, id)/totH

!     print*, "CN test", id, T, 4.72251d+07*molconc(3, id)/totH




!       contf=1.+(contf-1.)*20161.*molconc(4, id)/totH   ! CO case

!        contf=1.+(contf-1.)*4.72251d+07*molconc(3, id)/totH   ! CN case       


!       contf=1.+(contf-1.)*8.986*totFeI/totFe   ! FeI case   

!       contf=1.+(contf-1.)*1.125*totFeII/totFe ! FeII case 





!       print*, "Fe POPNUM I", id, totFeII/totFe ! totFeII/totFe, (totFeI+totFeII)/totFe


!       print*, "Fe POPNUM I", id, totFeI/totFe, totFeII/ANE ! totFeII/totFe, (totFeI+totFeII)/totFe
!       print*, "Fe POPNUM I", id, popul(99,id)/popul(100,id), totFe
!        print*, "Fe POPNUM I", id, popul(107,id)/popul(108,id), totFe

!      ATTENTION!!

!       contf=1.

!      ATTENTION!!

!************* opacities FALC ******************


      ABLIN(1:NFREQ)=ABLIN(1:NFREQ)+ABSO(1:NFREQ)*(contf-1.)
       
      
      do i=1, NFREQ 

    
      EMLIN(i)=EMLIN(i)+ABSO(i)*(contf-1.)*PLAN(max(NDPMIN,id))
    


      enddo







       do i=1,nfreq
!!!       write (200,310) ablin(i),emlin(i)
        write (200,FMT_LOPA) ablin(i)
        IF ((ABLIN(I) .LT. 0.) .OR. (EMLIN(I) .LT. 0.)) THEN
          PRINT '(i0,X,i0," ",$)',I,ID
          PRINT '("inibl0: SYNSUBM, OPAC: NEGATIVE OPACITY,'//
     &             ' EMISSIVITY : ",$)' ! WARNING!!!
          PRINT *,I, ABLIN(I),EMLIN(I)
        ENDIF
      enddo

  300 CONTINUE




!************* opacities FALC ******************

!      if (id .eq. 66) then
!      open (unit=100, file="../opac1_nf.txt",access="APPEND")
!      write(100, '(f10.4, e12.5, e12.5, e12.5)'), lambdat, 
!     * sum(ABSO(1:NFREQ))/NFREQ, sum(ABLIN(1:NFREQ))/NFREQ,
!     *  sum(ABSO(1:NFREQ))/sum(ABLIN(1:NFREQ))
!      close (100)
!      endif


!      if (id .eq. 66) then
!      open (unit=100, file="../total_nf.txt",access="APPEND")
!      write(100, '(f10.4, e12.5)'), lambdat, 
!     * sum(ABSO(1:NFREQ) + ABLIN(1:NFREQ))/NFREQ
!      close (100)
!      endif


!************* opacities FALC ******************


! *****************************
!     Test for continuum calculations


          


!      ABSO(1:NFREQ)=ABSO(1:NFREQ)*contf
      
!      do i=1, NFREQ 

!      if (id .le. 55) then
!      EMIS(i)=min(EMIS(i)*contf, ABSO(i)*BNUE(lambdat, 4500.))
!      else

!      EMIS(i)=EMIS(i)*contf
!      endif


!      enddo

!     Rinat, use it to calculate continuum

      ABSO(1:NFREQ)=ABSO(1:NFREQ)+ABLIN(1:NFREQ)*0.
      EMIS(1:NFREQ)=EMIS(1:NFREQ)+EMLIN(1:NFREQ)*0.


!      ABLIN(1:NFREQ)= ABLIN(1:NFREQ)*0.
!      EMLIN(1:NFREQ)= EMLIN(1:NFREQ)*0.

!      print*, 'contf', contf

!     Test for continuum calculations

!      print*, NFREQ
    
!      do i=1, NFREQ 
!      print*, i, freq(i)
!      enddo

 !     call readmollines
 !     call tester
 !    call calcmolopac(NFREQ, freq, T, ID, molopac, molemiss)
      
!      do i=1, NFREQ

!      print*, 'ratio test', ID, molopac(i)/ABSO(i), molemiss(i)/EMIS(i)

!      enddo


  !    print*, '************************************************'
  !    print*, 'scale test', ID, ABSO(200), molopac(200)
  !    print*, '************************************************' 

  !    ABSO(1:NFREQ)=ABSO(1:NFREQ)+molopac(1:NFREQ)
  !    EMIS(1:NFREQ)=EMIS(1:NFREQ)+molemiss(1:NFREQ)
      
 
  

 
c     write (200,310) ABSO(IJ),EMIS(IJ)

  225 CONTINUE
C
C     **** Detailed opacity and emissivity in hydrogen lines ****
C          (for IHYL=1)
C
C
      IF(IHYL.GE.0) THEN
        IF(IBVCS==24) THEN
          CALL HYDLIN(ID,T,ANE,ABLIN,EMLIN)
        ELSE
          CALL HYDLIN_IVANY(ID,T,ANE,ABLIN,EMLIN)
        ENDIF
        ABSO(1:NFREQ)=ABSO(1:NFREQ)+ABLIN(1:NFREQ)
        EMIS(1:NFREQ)=EMIS(1:NFREQ)+EMLIN(1:NFREQ)
        IF(ANY(ABLIN<0)) PRINT *,'ABLIN HYD < 0 @ ',ID
        IF(ANY(EMLIN<0)) PRINT *,'ABLIN HYD < 0 @ ',ID
      END IF
C
C     **** Detailed opacity and emissivity in HE II lines ****
C          (for IHE2L=1)
C
      IF(IHE2L.GE.0) THEN
        CALL HE2LIN(ID,T,ANE,ABLIN,EMLIN)
        ABSO(1:NFREQ)=ABSO(1:NFREQ)+ABLIN(1:NFREQ)
        EMIS(1:NFREQ)=EMIS(1:NFREQ)+EMLIN(1:NFREQ)
c     write (200,310) ABSO(IJ),EMIS(IJ)
      END IF
C
C     opacity due to detailed photoinization cross-section
C     (from tables; including resonance features)
C     The two routines may be called and correspond to different formats
C     as well as difference in INPUT!
C
C    COMMENTED OUT HERE !!!
c      CALL PHTION(ID,ABSO,EMIS)
c      CALL PHTX(ID,ABSO,EMIS)
C


        

!         do i=nfreq,1, -1
!         write(100, '(i3, f10.4, e12.5)'), id, 
!     * clight_cgs/freq(i)*1.d8, ABSO(i)
!         enddo


    !    print*, nfreq

    !    stop
     
!*      indo=nint(nfreq/40.)


   

!*      do i=1, 20

!*      wav_oi(i)=clight_cgs/freq(indo+2*indo*(i-1))*1.d8
!*      opac_oi(i)=sum(ABSO(2*indo*(i-1):2*indo*i))/(2.*indo) 
      
!*      enddo
   
   
   


      RETURN
      END SUBROUTINE
      
      function getContIdx(idxFreq)
      use MOD_ERROR
      use UTILS, only: ASSERT
      use SYNTHP_CONT, only: NFCONT
      implicit none
      include '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      integer i,getContIdx,idxFreq
      call assert(LOC(FRXIDX)/=0,'FRXIDX is not yet available')
      call assert(idxFreq>=1,'Frequency Index can not be negative')
      do i=1,NFCONT()-1
        if(idxFreq <= FRXIDX(i)) then
          getContIdx=i
          return
        endif
      enddo
      if(idxFreq<=FRXIDX(NFCONT())) then
        getContIdx=i
        return
      endif
      print *,'inibl0:getContIdx:ERROR:',
     &                    idxFreq,NFCONT(),FRXIDX(NFCONT())
      call error('Could not find FrxIdx')
      end function getContIdx
      
      
      
      SUBROUTINE HYDLIN_IVANY(ID,T,ANE,ABSOH,EMISH)
      use constants
      use MOD_HYDTAB,only: WLINE, MLINH,MHWL,NWLH,ILIN0
      use MOD_SYNSUBM, only: STARK0,DIVSTR,STARKA,FEAUTR
C
C     opacity and emissivity of hydrogen lines
C
      implicit real*8(a-h,o-z)
      INCLUDE '../inc/PARAMS.FOR'
      integer,intent(in   ) :: ID
      real*8, intent(in   ) :: T,ANE
      real*8, intent(inout) :: ABSOH(MFREQ),EMISH(MFREQ)
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      real*8,parameter:: FRH1=3.28805E15,FRH2=FRH1/4.,UN=1.,SIXTH=1./6.
      real*8,parameter:: CPP=4.1412E-16,CCOR=0.09,CPJ=157803.,CLST=1.1E3
      real*8,parameter:: CDOP=1.284523E12,TWO=2.
      real*8,parameter:: CID=0.02654 ! = pi*e^2/mc in cgs
      real*8,parameter:: C00=1.25E-9
      real*8,parameter:: CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/(CLIGHT_SI*1e10)
      real*8,parameter:: CID1=0.01497
      logical :: lwph,lquasi
      real*8  :: PJ(40),PRF0(99),OSCB(8),
     *          ABSO(MFREQ),EMIS(MFREQ)
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2IR
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
     $ IOPFE1
      COMMON/HYLPAR/IHYL,ILOWH,M10,M20
      COMMON/AUXVCS/XK,FXK,BETAD,DBETA
      integer :: NWLHYD(MLINH)
      real*8  :: WLHYD(MLINH,MHWL)
      real*8  :: PRFHYD(MLINH,MDEPTH,MHWL)
      real*8  :: PRF
      COMMON/HYDVCS/PRFHYD,WLHYD,NWLHYD
      !COMMON/BALVCS/PRFBAL(8,MDEPTH,36),WLBAL(8,36),NWLH(8)
      common/wprob/wph(mdepth,40),acor(mdepth),lwph
      real*8 WLINE2(8) !* Only used for checks
      DATA WLINE2  /6562.80, 4861.32, 4340.46, 4101.73,
     *             3970.07, 3889.05, 3835.38, 3797.90/
      DATA OSCB   /0.6407,  0.1193,   0.04467,  0.02209,
     *             1.27D-2, 8.036D-3, 5.429D-3, 3.851D-3/ !* Only used for checks
      DATA FRH    /3.289017E15/ 
      PARAMETER (HINV=1.5092973D26)
      logical :: explicitHy
      exp10(X)=exp(2.30258509299404568402d0*X)
C
      if(iath.le.0) return
      izz=1
      i0=1
      i1=nfreq
      ABSO (I0:I1)=0.
      EMIS (I0:I1)=0.
      ABSOH(I0:I1)=0.
      EMISH(I0:I1)=0.
      T1=UN/T
      SQT=SQRT(T)
      ANES=EXP(SIXTH*LOG(ANE)) !** ANES=ANE^(1/6)
C
C     populations of the first 40 levels of hydrogen
C
      ANP=POPUL(NKH,ID)
      !** PP = CPP * ANE*ANP/T^(3/2)
      PP=CPP*ANE*ANP*T1/SQT
      NLH=N1H-N0HN+1
      if(ifwop(n1h).lt.0) nlh=nlh-1
      DO IL=1,40
        X=IL*IL
        IF(IL.LE.NLH) THEN
          PJ(IL)=POPUL(N0HN+IL-1,ID)/X
        ELSE
          !*    =4.1412E-16*ANE*ANP/T^3/2 *exp(157803./(IL^2*T)
          PJ(IL)=PP*EXP(CPJ/X*T1)
        ENDIF
      ENDDO
      p2=pp*exp(cpj4*t1)
C
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      F00=C00*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID))
C
C     determination of the last observable line
C
      MLST=CLST*EXP(-LOG(ANE)/7.5)
      IF(MLST.GT.40) MLST=40
C
C -------------------------------------------------------------------
C     overall loop over spectral series (only in the infrared region)
C -------------------------------------------------------------------
C
      ISERL=ILOWH
      ISERU=ILOWH
      IF(WLAM(I0).GT.17000..AND.WLAM(I1).LT.21000.) THEN
         ISERL=3
         ISERU=4
       ELSE IF(WLAM(I0).GT.22700.) THEN
         ISERL=4
         ISERU=5
         IF(WLAM(I0).GT.32800.) ISERU=6
         IF(WLAM(I0).GT.44660.) ISERU=7
      END IF
C
      DO 200 I=ISERL,ISERU
      II=I*I
      XII=UN/dble(II)
      PLTEI=PP*EXP(CPJ*T1*XII)*II
      POPI=PJ(I)*II
      IF(I.EQ.1) FRH=3.28805E15
      FEDGE=FRH*XII
      FLST=FEDGE-FRH/MLST**2
      IF(I.LE.NLH) FEDGE=ENION(I+N0HN-1)*HINV
C
C     determination of which hydrogen lines contribute in a current
C     frequency region
C
      M1=M10
      IF(I.LT.ILOWH) M1=ILOWH-1
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.3..and.M1.LE.6.AND.I.EQ.2) GO TO 10
      IF(grav.lt.3..and.M1.LE.4.AND.I.EQ.1) GO TO 10
      IF(M1.GT.MLST .and. .not.lwph) GO TO 120
      M1=M1-1
      M2=M20+3
      IF(M1.LT.I+1) M1=I+1
   10 CONTINUE
      if(grav.gt.3.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3
      end if
      if(grav.gt.6.) then
         m2=m2+2
         m1=m1-1
         if(m1.gt.i+6) m1=m1-1
      end if
      IF(M1.LT.I+1) M1=I+1
      IF(M2.GT.40) M2=40
      if(id.eq.1) write(6,666) i,m1,m2
  666 format(/' hydrogen lines contribute - ilow=',i2,', iup from ',i3,
     *       ' to',i3/)
C
      A=0.
      E=0.
      ABSO(I0:I1)=0.
      EMIS(I0:I1)=0.
C
C     loop over lines which contribute at given wavelength region
C
      DO 100 J=M1,M2
         !IF(I.EQ.2.AND.J.LE.10.AND.IBVCS.GT.0) THEN
        explicitHy = IBVCS /= 0 .and. I<=4 .and. J<=22
        !explicitHy = IBVCS /= 0 .and. I==2 .and. J<=8
        IF(explicitHy) explicitHy = (ILIN0(I,J) /=0 )
        !explicitHy=explicitHy.or.(I==2.and.J<=8)
        EXPLICIT_HY: IF(explicitHy) THEN
          ILINE=ILIN0(I,J)
          !IF(I==2.and.J<=8.and.ILINE==0) ILINE=J-2
          NWL=NWLHYD(ILINE)
          PRF0(1:NWL)=PRFHYD(ILINE,ID,1:NWL) !-10.
          CALL STARK0(I,J,1,XKIJ,WL0,FIJ,FIJ0)
          FID=(FIJ+FIJ0)*CID * WL0**2/2.99792458d18/F00
          DO IJ=I0,I1
            AL=ABS(WLAM(IJ)-WLINE(I,J))/F00
            IF(AL.LT.1.E-4) AL=1.E-4
            AL=LOG10(AL)
            !*** Interpolate the PRF (phi_v)
            !** Find Index of the appropiate hydrogen broadening
            HYD_MIN_IW:DO IWL=1,NWL-1
              IW0=IWL
              IF(AL.LE.WLHYD(ILINE,IWL+1)) exit HYD_MIN_IW
            ENDDO HYD_MIN_IW
            IW1=IW0+1
            !** linear interpolation in log10
            PRF=(PRF0(IW0)*(WLHYD(ILINE,IW1)-AL)+PRF0(IW1)*
     *          (AL-WLHYD(ILINE,IW0)))/
     *          (WLHYD(ILINE,IW1)-WLHYD(ILINE,IW0))
            SG=EXP10(PRF)*FID*wph(id,j)
            ABSO(IJ)=ABSO(IJ)+SG
            EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
              if(ID==50) then  !*** DEBUG
                write(9913,'(7(p1e12.4,X),$)')
     &           ! 1           2     3 4  5   6   7
     &          WLAM(IJ),WLINE(I,J),AL,SG,AL, ANE,FID!*** DEBUG
                 !                            8    9
                write(9913,'(100(i4,X),$)'), IW0,ILINE,I,J
              endif
          ENDDO
        ELSE EXPLICIT_HY
          CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
          FXK=F00*XKIJ
          FXK1=UN/FXK
          DOP=DOP0/WL0
          DBETA=WL0*WL0*CINV*FXK1
          BETAD=DOP*DBETA
          FID=CID*FIJ*DBETA*wph(id,j)
          FID0=CID1*FIJ0/DOP*wph(id,j)
          CALL DIVSTR(BETAD,AD,DIV)
          lquasi=i.eq.1.and.j.eq.2.and.iophli.ge.3
          DO 60 IJ=I0,I1
            fr=freq(ij)
            IF(FREQ(IJ).GT.FLST .and. .not.lwph) cycle
            BETA=ABS(WLAM(IJ)-WL0)*FXK1
c              if(beta.gt.5.e4) go to 60
            SG=STARKA(BETA,BETAD,AD,DIV,TWO)*FID
            if(iophli.eq.2.and.i.eq.1.and.j.eq.2)
     *        sg=sg*feautr(fr,id)
            if(fid0.gt.0.) then
              xd=beta/betad
              if(xd.lt.5.) sg=sg+exp(-xd*xd)*fid0
            end if
            ABSO(IJ)=ABSO(IJ)+SG
            EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
   60     ENDDO
        END IF EXPLICIT_HY
  100 ENDDO
C
C     resulting absorption and emission coefficients
C
      EMIS(I0:I1)=EMIS(I0:I1)*II
      ABSO(I0:I1)=ABSO(I0:I1)*POPI-EMIS(I0:I1)
C
C     pseudocontinuum for wavelengths smaller than the last observable
C     line
C
  120 CONTINUE
      if(lwph) GO TO 150
      selectcase(I)
        case(1); SGF=6.313E-18
        case(2); SGF=1.387E-17
        case(3); SGF=2.156E-17
        case(4); SGF=2.929E-17
        case(5); SGF=3.705E-17
        case(6); SGF=4.483E-17
      endselect
!       IF(I.EQ.1) THEN
!          SGF=6.313E-18
!        ELSE IF(I.EQ.2) THEN
!          SGF=1.387E-17
!        ELSE IF(I.EQ.3) THEN
!          SGF=2.156E-17
!        ELSE IF(I.EQ.4) THEN
!          SGF=2.929E-17
!        ELSE IF(I.EQ.5) THEN
!          SGF=3.705E-17
!        ELSE IF(I.EQ.6) THEN
!          SGF=4.483E-17
!       END IF
      E=PLTEI*SGF
      A=POPI*SGF
      DO 130 IJ=I0,I1
         F=FREQ(IJ)
         IF(F.GE.FEDGE.OR.F.LT.FLST) cycle
         EMIS(IJ)=E*EXP(-4.79928E-11*F*T1)
         ABSO(IJ)=A-EMIS(IJ)
  130 ENDDO
  150 CONTINUE
C
C     add contributions for different spectral series
C
      ABSOH(I0:I1)=ABSOH(I0:I1)+ABSO(I0:I1)
      EMISH(I0:I1)=EMISH(I0:I1)+EMIS(I0:I1)
  200 CONTINUE
C
C     finally, multiply EMISH by 2h nu^3/c^2
C
      DO 210 IJ=I0,I1
         F=FREQ(IJ)
         F15=F*1.E-15
         EMISH(IJ)=1.4743E-2*F15*F15*F15*EMISH(IJ)
  210 CONTINUE
      RETURN
      END SUBROUTINE
      
     

      FUNCTION airlambda(vaclambda)
!    translate vacuum to the airlambda lambda in A                             
                                                                                
      IMPLICIT NONE
      real*8 vaclambda, airlambda, sig, n

      sig=1.d4/(vaclambda)
      n=1.d0+6.4328d-5+(2.94981d-2)/(146.-sig**2.)+(2.554d-4)/
     $ (41-sig**2.)
      airlambda=vaclambda/n
      return
      END FUNCTION





 
      end module
