      MODULE MOD_BALINI
      use MOD_HYDTAB, only: MLINH,MHWL,MHE,MHT,ILIN0,WLINE,NWLH,NTH,NEH

      contains

C ********************************************************************
C ********************************************************************
C
      !*** replaced by HYDINI
      SUBROUTINE BALINI(IBVCS)
      use constants
      use MOD_HYDTAB
C
C     Initializes necessary arrays for evaluating Balmer-alpha thru
C     Balmer-delta absorption profiles from the modified
C     Vidal, Cooper and Smith (Ap.J.Suppl. 25, 37, 1973) tables
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: IBVCS
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
!      COMMON/VCSDAT/WL(36,8),XT(7,8),XNE(11,8),PRF(36,7,11)
!     *              ,NWLH(8),NTH(8),NEH(8)

      COMMON/AUXVCS/XK,FXK,BETAD,DBETA
      integer :: NWLBAL(MLINH)
      real*8  :: WLBAL(MLINH,MHWL)
      real*8  :: PRFBAL(MLINH,MDEPTH,MHWL)
      COMMON/HYDVCS/PRFBAL,WLBAL,NWLBAL
      DATA NLINE  /8/
!       WLINE=0
!       ILIN0=0
      ! OPEN(UNIT=IBVCS,STATUS='OLD',READONLY)
      WLINE(2,3:10) =  (/6562.80, 4861.32, 4340.46, 4101.73,
     *             3970.07, 3889.05, 3835.38, 3797.90/)
      BALINI_LINES: DO ILINE=1,NLINE
        WL0=WLINE(2,ILINE+2)
        ILIN0(2,ILINE+2)=ILINE
        !**
        !* read the modified VCS tables, which have to be stored in file
        !* unit IBVCS (which is the input parameter in the progarm)
        !**
        READ(IBVCS,501) JLINE,ILINE1
        READ(IBVCS,502) NWL,(WL(I,ILINE),I=1,NWL)
        READ(IBVCS,503) NT,(XT(I,ILINE),I=1,NT)
        READ(IBVCS,504) NE,(XNE(I,ILINE),I=1,NE)
        READ(IBVCS,500)
        NWLH(ILINE)=NWL
        NWLBAL(ILINE)=NWL
        NTH(ILINE)=NT
        NEH(ILINE)=NE
        DO I=1,NWL
          IF(WL(I,ILINE).LT.1.E-4) WL(I,ILINE)=1.E-4
          WLBAL(ILINE,I)=LOG10(WL(I,ILINE))
        ENDDO
        !*
        DO IE=1,NE
          DO IT=1,NT
            READ(IBVCS,500)
            READ(IBVCS,505) (PRF(IWL,IT,IE),IWL=1,NWL)
          ENDDO
        ENDDO
        !*****
        !* coefficient for the asymptotic profile is determined from
        !* the input data
        !*
        XCLOG=PRF(NWL,1,1)+2.5*LOG10(WL(NWL,ILINE))+31.5304-
     *         XNE(1,ILINE)-2.*LOG10(WL0)
        XKLOG=0.6666667*(XCLOG-0.176)
        XK=EXP(XKLOG*2.3025851)
        !*
        DO ID=1,ND
          !**********
          !* temperature is modified in order to account for the
          !* effect of turbulent velocity on the Doppler width
          !*
          T=TEMP(ID)+6.06E-9*VTURB(ID)
          ANE=ELEC(ID)
          TL=LOG10(T)
          ANEL=LOG10(ANE)
          F00=1.25E-9*ANE**0.666666667
          FXK=F00*XK
          DOP=1.E8/WL0*SQRT(1.65E8*T)
          DBETA=WL0*WL0/(CLIGHT_SI*1e10)/FXK
          BETAD=DBETA*DOP
          !*******
          !* interpolation to the actual values of temperature and
          !* electron density. The result is stored at array PRFBAL,
          !* having indices
          !* ILINE (line number: 1 for H-alpha,...,4 for H-delta, etc.);
          !* ID - depth index
          !* IWL - wavelength index
          !*
          DO IWL=1,NWL
              CALL INTVCS(PROF,TL,ANEL,IWL,ILINE)
              PRFBAL(ILINE,ID,IWL)=PROF
          ENDDO
        ENDDO
      ENDDO BALINI_LINES
C
  500 FORMAT(1X)
  501 FORMAT(//11X,I1,12X,I1/)
  502 FORMAT(1X,I4,6E10.3,5(/5X,6F10.4))
  503 FORMAT(1X,I4,7F10.3)
  504 FORMAT(1X,I4,6F10.2/5X,5F10.2)
  505 FORMAT(10F8.3)
      RETURN
      END SUBROUTINE

      !**** replaced by INTHYD
      SUBROUTINE INTVCS(W0,X0,Z0,IWL,ILINE)
      use MOD_SYNSUBM
      use MOD_HYDTAB

C
C     Interpolation in temperature and electron density from the
C     modified VCS tables for Balmer lines to the actual valus of
C     temperature and electron density
C
      USE MOD_ERROR
c      INCLUDE 'IMPLIC.FOR'
      INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IWL,ILINE
      real*8, intent(in   ) :: X0,Z0
      real*8, intent(  out) :: W0
      PARAMETER (TWO=2.D0)
      COMMON/AUXVCS/XK,FXK,BETAD,DBETA
!      COMMON/VCSDAT/WL(36,8),XT(7,8),XNE(11,8),PRF(36,7,11)
!     *              ,NWLH(8),NTH(8),NEH(8)
      DIMENSION ZZ(3),XX(3),WX(3),WZ(3)
C
      NX=3
      NZ=3
      NT=NTH(ILINE)
      NE=NEH(ILINE)
      COUNTER=COUNTER+1
C
C     for values lower than the lowest grid value of electron density
C     the profiles are determined by the approximate expression
C     (see STARKA); not by an extrapolation in the VCS tables which may
C     be very inaccurate
C
      IF(Z0.LT.XNE(1,ILINE)*0.99) THEN
         CALL DIVSTR(BETAD,A,DIV)
         W0=STARKA(WL(IWL,ILINE)/FXK,BETAD,A,DIV,TWO)*DBETA
         W0=LOG10(W0)
         GO TO 500
      END IF
C
C     Otherwise, one interpolates (or extrapolates for higher than the
C     highes grid value of electron density) in the VCS tables
C
      DO 10 IZZ=1,NE-1
         IPZ=IZZ
         IF(Z0.LE.XNE(IZZ+1,ILINE)) GO TO 20
   10 CONTINUE
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1  ! N0Z = min(N0Z, NE-NZ+1)
      IF(N0Z.LT.1) THEN
        PRINT '("N0Z=",I7," NE=",I7," NZ=",I7)',N0Z, NE, NZ
        call ERROR('ERROR: synsubm: Index to small: N0Z')
      ENDIF
      N1Z=N0Z+NZ-1
C
      DO 300 IZZ=N0Z,N1Z
         I0Z=IZZ-N0Z+1
         ZZ(I0Z)=XNE(IZZ,ILINE)
C
C     Likewise, the approximate expression instead of extrapolation
C     is used for higher that the highest grid value of temperature,
C     if the Doppler width expressed in beta units (BETAD) is
C     sufficiently large (> 10)
C
         IF(X0.GT.1.1*XT(NT,ILINE).AND.BETAD.GT.10.) THEN
            CALL DIVSTR(BETAD,A,DIV)
            W0=STARKA(WL(IWL,ILINE)/FXK,BETAD,A,DIV,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500
         END IF
C
C     Otherwise, normal inter- or extrapolation
C
C     Both interpolations (in T as well as in electron density) are
C     by default the quadratic interpolations in logarithms
C
         DO 30 IX=1,NT-1
            IPX=IX
            IF(X0.LE.XT(IX+1,ILINE)) GO TO 40
   30    CONTINUE
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
         N1X=N0X+NX-1
         DO 200 IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XT(IX,ILINE)
            WX(I0)=PRF(IWL,IX,IZZ)
  200       CONTINUE
         WZ(I0Z)=YINT(XX,WX,X0)
  300 CONTINUE
      W0=YINT(ZZ,WZ,Z0)
  500 CONTINUE
      RETURN
      END SUBROUTINE
      END MODULE
