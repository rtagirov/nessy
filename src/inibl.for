       module mod_inibl

       contains

       subroutine inibl(wavarr, sigarr, N, NFDIM)

       use MOD_SYNSUBM
       use SYNTHP_CONT,only: FREQC,ABSOC,NFCONT
       use UTILS,only:assert
       use MOD_DECF_SYN
       use CONSTANTS,only: CLIGHT_SI
       use MOD_BALINI
       use MOD_HYDTAB
       use MOD_chemeq

       use mod_opac

!     auxiliary initialization procedure

      implicit none

      integer,intent(in) :: N,      NFDIM
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
      real*8,dimension(MFREQ) :: ABSO, EMIS

      DIMENSION WAVARR(N, NFDIM), SIGARR(N, NFDIM)

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


      real*8 :: opac_start, opac_finish


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

         CALL CROSET(WAVARR, SIGARR, N, NFDIM)

         call readmollines

         call readMolconc(ND)

         call cpu_time(opac_start)

         DO ID = 1, ND

            CALL OPAC(ID, 0, ABSO, EMIS, WAVARR, SIGARR, N, NFDIM)

            ABSTD(ID)=MINVAL(ABSOC(:NFCONT()))

         ENDDO

         call cpu_time(opac_finish)

         print*, 'inibl: opac execution time = ', (opac_finish - opac_start)

         IF(INLTE.GT.0) CALL NLTE(0,1,1,1,A1,A2,A3)

         IF(cards.ABEMLIN/=card_params.ABEMLIN_READ) CALL INILIN(INLIST)

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
       
      end module
