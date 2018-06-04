      module mod_opac

      character*(*), private, parameter :: FMT_LOPA = '(1pe12.5)'

      contains

      SUBROUTINE OPAC(ID, MODE, ABSO, EMIS, WAVARR, SIGARR, N, NFDIM)

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

      use MOD_HMINUSFF
      use MOD_LINOP
      use MOD_SYNSUBM
      use MOD_DECF_SYN
      use UTILS
      use SYNTHP_CONT, only: ABSOC,EMISC,FREQC,NFCONT, absoc_rayleigh
      use constants, only:   CLIGHT_SI, CLIGHT_CGS, BoltzmannConstantEV,
     $                       PlanckConstantEV
      use MOD_chemeq

      use mod_gff_temp

      use common_block
      use phys

      implicit none

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/LINDAT.FOR'
      
      integer,intent(in) :: id, NFDIM, N
      integer,intent(in) :: MODE

      real*8, intent(in), dimension(N, NFDIM) :: WAVARR, SIGARR

      real*8, intent(inout),dimension(MFREQ) :: ABSO, EMIS

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

      integer :: UnitRatios

      integer :: idx

      real*8 :: wvl
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

        AbsCoefRest = 0.

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
                   sg=sigk(fr, ii, 1, WAVARR, SIGARR, N, NFDIM)*dw
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

              CALL GFF_TEMP(1, 1.0D0, PlanckConstantEV * FR, BoltzmannConstantEV * T, SG)

              SF2=SF2+SG-UN
            endif

            X=EXP(-HKF)
            X1=UN-X
!***********************************************************************************
!RINAT TAGIROV, ALEXANDER SHAPIRO
!X1 has been added to account for the induced emission for every element except H^-
!the induced emission for H^- had been taken into account in hminusff
            SFF = X1 * SF1 * SF2 
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

!       Additional opacities

        CALL OPADD(0,IJ,ID,FR,ABAD,EMAD,SCAD)

!       Total continuum opacity and emissivity

        BNE = BNU * X * ANE

CMH     HMINUSFF fit includes stimulated emission
cmh     correction by X1 = 1. - exp(-h*nu/k*T) obsolete for Hminus

        ABSOC(IJ)=ABF+ANE*(AFF-X*EBF)+ANE*SIGEL+ABAD+ABLY
        EMISC(IJ)=BNE*(AFF/X1+EBF)+EMAD+EMLY

        if (rayleigh) then

            wvl = light_speed * 1.0d+8 / fr

            if (wvl .le. 1600.0d0) then

                absoc_rayleigh(ij) = popul(2, id) * sigma_rayleigh(wvl)

                absoc(ij) = absoc(ij) + absoc_rayleigh(ij)
                emisc(ij) = emisc(ij) + absoc_rayleigh(ij) * plan(id)

            endif

        endif

!***********************************************************************************
!RINAT TAGIROV
        IF ((IJ .EQ. 1) .AND. (MODE .EQ. 0)) THEN

           HydIonDeg = POPUL(12, ID) / SUM(POPUL(1:12, ID))

           HeIonDeg =  POPUL(23, ID) / SUM(POPUL(13:23, ID))

!GFF_TEMP is the gaunt factor for any Rydberg system as a function of temperature and frequency
!See the subroutine for more details

           CALL GFF_TEMP(1, 1.0D0, PlanckConstantEV * FR, BoltzmannConstantEV * T, SG)

           WRITE(UnitRatios, 2000), xlam, ID, T, ANE,
     $     AbsCoefHM * 1.0D+2 / (AbsCoefH + AbsCoefHM),
     $     AbsCoefRest * 1.0D+2 / (AbsCoefH + AbsCoefHM),
     $     GFREE(T, FR / CH), SG, HydIonDeg * 1.0D+2, HeIonDeg * 1.0D+2

        ENDIF
!***********************************************************************************

      ENDDO FCONT_LOOP

      CLOSE(UnitRatios)

      if(mode.eq.0) return

      !***************************************************************
      !interpolated continuum and hydrogen line opacity and
      !emissivity for all frequencies
      i2=0
      do i=1,NFCONT()-1
        i1=i2+1
        i2=FRXIDX(i+1)
        ABSO(i1:i2)=FRX1(i1:i2)*ABSOC(i+1)+FRX2(i1:i2)*ABSOC(i)
        EMIS(i1:i2)=FRX1(i1:i2)*EMISC(i+1)+FRX2(i1:i2)*EMISC(i)
      enddo

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

        CALL LINOP(ID, ABLIN, EMLIN)

        if(cards.ABEMLIN==card_params.ABEMLIN_WRITE) then
          write(201,*) NFREQ,ID
          write(201,*) ABLIN(1:NFREQ)
          write(201,*) EMLIN(1:NFREQ)
        endif
      ENDIF  ! EMABLIN /= EMABLIN_READ


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

       totFe=popul(101,id)+popul(102,id)+popul(103,id)+popul(104,id)+popul(105,id)+popul(106,id)   
       totFeII=popul(106, id)
       totFeI=totFe-totFeII

       totH=sum(popul(1:12,id))

       ABLIN(1 : NFREQ) = ABLIN(1 : NFREQ) + ABSO(1 : NFREQ) * (contf - 1.0d0)

       EMLIN(1 : NFREQ) = EMLIN(1 : NFREQ) + ABSO(1 : NFREQ) * (contf - 1.0d0) * PLAN(max(NDPMIN, id))

       do i = 1, nfreq

        write (200,FMT_LOPA) ablin(i)
        IF ((ABLIN(I) .LT. 0.) .OR. (EMLIN(I) .LT. 0.)) THEN
          PRINT '(i0,X,i0," ",$)',I,ID
          PRINT '("opac: NEGATIVE OPACITY,'//
     &             ' EMISSIVITY : ",$)' ! WARNING!!!
          PRINT *,I, ABLIN(I),EMLIN(I)
        ENDIF
      enddo

  300 CONTINUE

!     Rinat, use it to calculate continuum

      ABSO(1:NFREQ)=ABSO(1:NFREQ)+ABLIN(1:NFREQ)!*0.
      EMIS(1:NFREQ)=EMIS(1:NFREQ)+EMLIN(1:NFREQ)!*0.

  225 CONTINUE

!     Detailed opacity and emissivity in hydrogen lines (for IHYL = 1)

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

      END IF

C     opacity due to detailed photoinization cross-section
C     (from tables; including resonance features)
C     The two routines may be called and correspond to different formats
C     as well as difference in INPUT!
C
C    COMMENTED OUT HERE !!!
c      CALL PHTION(ID,ABSO,EMIS)
c      CALL PHTX(ID,ABSO,EMIS)

      return

      end subroutine

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

      end module
