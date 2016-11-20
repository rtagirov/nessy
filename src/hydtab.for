      MODULE MOD_HYDTAB
      !** for HYDTAB
      real*8  :: WLINE(4,22)
      integer,parameter :: MLINH=4*22 ! maxLine
      integer,parameter :: MHWL=80   ! maxWavelength
      integer,parameter :: MHE=20   ! max#electrons
      integer,parameter :: MHT=7    ! max#Temp
      !common block /VCSDAT/
      real*8 :: WL(MHWL,MLINH),XT(MHT,MLINH)
      real*8 :: XNE(MHE,MLINH),PRF(MHWL,MHT,MHE)
      integer :: NWLH(MLINH),NTH(MLINH),NEH(MLINH)
      integer :: ILIN0(4,22)
      !end common block VCSDAT
!       contains
!       SUBROUTINE HYDINI(IHYDPR)
!       use UTILS,only:assert
! C
! C     Initializes necessary arrays for evaluating hydrogen line profiles
! C     from the Schoening and Butler or Lemke tables
! C
! !      use MOD_SYNSUBM,only:STARK0
!       use PARAMS_ARRAY,only:NDDIM
!       implicit none
!       INCLUDE '../INCLUDE2/PARAMS.FOR'
!       !*** MODELP: TEMP,VTURB,ELEC
!       INCLUDE '../INCLUDE2/MODELP.FOR'
! !       COMMON/VCSDAT/WL(MHWL,MLINH),XT(MHT,MLINH),
! !      *              XNE(MHE,MLINH),PRF(MHWL,MHT,MHE),
! !      *              NWLH(MLINH),NTH(MLINH),NEH(MLINH)
!       real*8  ::    XK,FXK,BETAD,DBETA
!       COMMON/AUXVCS/XK,FXK,BETAD,DBETA
!       integer :: IHYDP(4)
!       !*** added for debuging
!       !COMMON/BALVCS/PRFBAL(8,MDEPTH,36),WLBAL(8,36),NWLH(8)
!       integer :: NWLHYD(MLINH)
!       real*8  :: WLHYD(MLINH,MHWL)
!       real*8  :: PRFHYD(MLINH,MDEPTH,MHWL)
! 
!       COMMON/HYDVCS/PRFHYD,WLHYD,NWLHYD
!       !*** end debugging
!       CHARACTER*1 :: CHAR
!       integer,save :: INIT = 0
!       integer :: I, J,IZZ,IHYDPR,ILEMKE,NLINE,ILINE,NWL,NT,NE
!       integer :: IE,IT,IWL,ID,IHYDP0,NTAB,ITAB,ILINEB,NLLY
!       integer :: ILI,INE,ILNE
!       real*8  :: TMIN,DLA,DLE,DLT,QLT
!       real*8  :: T,ANE,ANEL,F00,DOP,PROF,TL,ALMIN,ANEMIN
!       real*8  :: FIJ,WL0,FIJ0,XCLOG,XKLOG
! C
!       IF(INIT.EQ.0) THEN
!         DO I=1,4
!           DO J=I+1,22
!             CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0)
!             WLINE(I,J)=WL0
!             ! OSCH(I,J)=FIJ+FIJ0
!           END DO
!           print *
!         END DO
!         INIT=1
!       END IF
!       ILIN0(:4,:22)=0
!       !
!       ! --------------------------------------------
!       !     Schoening-Butler tables - for IHYDPR < 0
!       ! --------------------------------------------
!       !
!       IF(IHYDPR.LT.0) THEN
!       IHYDPR=-IHYDPR
!       ILEMKE=0
!       NLINE=12
!       !
!       ! OPEN(UNIT=IHYDPR,STATUS='OLD')
!       !
!       DO I=1,12
!          READ(IHYDPR,500)
!       END DO
!       DO 100 ILINE=1,NLINE
!         !
!         !     read the tables, which have to be stored in file
!         !     unit IHYDPR (which is the input parameter in the progarm)
!         !
!          READ(IHYDPR,501) I,J
!          IF(ILINE.EQ.12) J=10
!          WL0=WLINE(I,J)
!          ILIN0(I,J)=ILINE
!          READ(IHYDPR,*) CHAR,NWL,(WL(I,ILINE),I=1,NWL)
!          READ(IHYDPR,*) CHAR,NT,(XT(I,ILINE),I=1,NT)
!          READ(IHYDPR,*) CHAR,NE,(XNE(I,ILINE),I=1,NE)
!          READ(IHYDPR,500)
!          NWLH(ILINE)=NWL
!          NWLHYD(ILINE)=NWL
!          NTH(ILINE)=NT
!          NEH(ILINE)=NE
! C
!          DO 10 I=1,NWL
!             IF(WL(I,ILINE).LT.1.E-4) WL(I,ILINE)=1.E-4
!             WLHYD(ILINE,I)=LOG10(WL(I,ILINE))
!    10    ENDDO
! C
!          DO IE=1,NE
!             DO IT=1,NT
!                READ(IHYDPR,500)
!                READ(IHYDPR,*) (PRF(IWL,IT,IE),IWL=1,NWL)
!             ENDDO
!          ENDDO
! C
! C        coefficient for the asymptotic profile is determined from
! C        the input data
! C
!          XCLOG=PRF(NWL,1,1)+2.5*LOG10(WL(NWL,ILINE))+31.5304-
!      *         XNE(1,ILINE)-2.*LOG10(WL0)
!          XKLOG=0.6666667*(XCLOG-0.176)
!          XK=EXP(XKLOG*2.3025851)
! C
!          DO ID=1,ND
! C
! C           temperature is modified in order to account for the
! C           effect of turbulent velocity on the Doppler width
! C
!             T=TEMP(ID)+6.06E-9*VTURB(ID)
!             ANE=ELEC(ID)
!             TL=LOG10(T)
!             ANEL=LOG10(ANE)
!             F00=1.25E-9*ANE**0.666666667
!             FXK=F00*XK
!             DOP=1.E8/WL0*SQRT(1.65E8*T)
!             DBETA=WL0*WL0/2.997925E18/FXK
!             BETAD=DBETA*DOP
! C
! C       interpolation to the actual values of temperature and electron
! C       density. The result is stored at array PRFHYD, having indices
! C       ILINE (line number: 1 for L-alpha,..., 4 for H-delta, etc.);
! C                           5 for H-alpha,..., 8 for H-delta, etc.)
! C       ID - depth index
! C       IWL - wavelength index 
! C
!             DO IWL=1,NWL
!                CALL INTHYD(PROF,TL,ANEL,IWL,ILINE,ILEMKE)
!                PRFHYD(ILINE,ID,IWL)=PROF
!             END DO
!          END DO
!   100 ENDDO
! C
!   500 FORMAT(1X)
!   501 FORMAT(12X,I1,9X,I1)
!   502 FORMAT(2X,I4,2x,6E12.4/8X,6E12.4,4(/4X,6F12.4))
!   503 FORMAT(2X,I4,F10.3,5F12.3,(/4X,6F12.3))
!   504 FORMAT(2X,I4,F10.3,5F12.2,2(/4X,6F12.2))
!   505 FORMAT(10F8.3)
! C
!       IHYDPR=-IHYDPR
!       RETURN
!       END IF
! C
! C ---------------------------------
! C     read Lemke tables
! C ---------------------------------
! C
!       ILEMKE=1
!       ILINE=0
!       IHYDP0=IHYDPR
!       IF(IHYDPR.LT.100) THEN
!          NTAB=1
!          IHYDP(1)=IHYDPR
!        ELSE IF(IHYDP0.LT.10000) THEN
!          NTAB=2
!          IHYDP(1)=IHYDP0/100
!          IHYDP(2)=IHYDP0-100*IHYDP(1)
!        ELSE IF(IHYDP0.LT.1000000) THEN
!          NTAB=3
!          IHYDP(1)=IHYDP0/10000
!          IHYDP(2)=(IHYDP0-10000*IHYDP(1))/100
!          IHYDP(3)=IHYDP0-10000*IHYDP(1)-100*IHYDP(2)
!        ELSE
!          NTAB=4
!          IHYDP(1)=IHYDP0/1000000
!          IHYDP(2)=(IHYDP0-1000000*IHYDP(1))/10000
!          IHYDP(3)=(IHYDP0-1000000*IHYDP(1)-10000*IHYDP(2))/100
!          IHYDP(4)=IHYDP0-1000000*IHYDP(1)-10000*IHYDP(2)-
!      *            100*IHYDP(3)
!       END IF
! c      
!       ILINE=0
!       DO ITAB=1,NTAB
!       ILINEB=ILINE
!       IHYDPR=IHYDP(ITAB)
!       READ(IHYDPR,*) NLLY
!       HEADER: DO ILI=1,NLLY
!          ILINE=ILINE+1
!          READ(IHYDPR,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,
!      *                  NWL,NE,NT
!          WL0=WLINE(I,J)
!          ILIN0(I,J)=ILINE
!          NWLH(ILINE)=NWL
!          NWLHYD(ILINE)=NWL
!          NTH(ILINE)=NT
!          NEH(ILINE)=NE
!          DO IWL=1,NWL
!             WL(IWL,ILINE)=ALMIN+(IWL-1)*DLA
!             WLHYD(ILINE,IWL)=WL(IWL,ILINE)
!             WL(IWL,ILINE)=EXP(2.3025851*WL(IWL,ILINE))
!          ENDDO
!          DO INE=1,NE
!             XNE(INE,ILINE)=ANEMIN+(INE-1)*DLE
!          ENDDO
!          DO IT=1,NT
!             XT(IT,ILINE)=TMIN+(IT-1)*DLT
!          ENDDO
!       ENDDO HEADER
! c
!       DO ILI=1,NLLY         
!          ILNE=ILINEB+ILI
!          call assert(ILINE<=MLINH)
!          NWL=NWLH(ILNE)
!          READ(IHYDPR,500)
!          DO INE=1,NEH(ILNE)
!             DO IT=1,NTH(ILNE)
! c              READ(IHYDPR,506) QLT,(PRF(IWL,IT,INE),IWL=1,NWL)
!                READ(IHYDPR,*) QLT,(PRF(IWL,IT,INE),IWL=1,NWL)
!             END DO
!          END DO
!   506    FORMAT((f8.4,8f9.4))
! C
! C        coefficient for the asymptotic profile is determined from
! C        the input data
! C
!          XCLOG=PRF(NWL,1,1)+2.5*WLHYD(ILNE,NWL)-0.477121
!          XKLOG=0.6666667*XCLOG
!          XK=EXP(XKLOG*2.3025851)
! C
!          DO ID=1,ND
! C
! C           temperature is modified in order to account for the
! C           effect of turbulent velocity on the Doppler width
! C
!             T=TEMP(ID)+6.06E-9*VTURB(ID)
!             ANE=ELEC(ID)
!             TL=LOG10(T)
!             ANEL=LOG10(ANE)
!             F00=1.25E-9*ANE**0.666666667
!             FXK=F00*XK
!             DOP=1.E8/WL0*SQRT(1.65E8*T)
!             DBETA=WL0*WL0/2.997925E18/FXK
!             BETAD=DBETA*DOP
! C
! C       interpolation to the actual values of temperature and electron
! C       density. The result is stored at array PRFHYD, having indices
! C       ILINE - line number
! C       ID    - depth index
! C       IWL   - wavelength index 
! C
!             DO IWL=1,NWL
!                CALL INTHYD(PROF,TL,ANEL,IWL,ILNE,ILEMKE)
!                PRFHYD(ILNE,ID,IWL)=PROF
!             END DO
!         END DO
!       END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C ********************************************************************
! C
! !*** micha: add key ILEMKE to the interface
!       SUBROUTINE INTHYD(W0,X0,Z0,IWL,ILINE,ILEMKE)
!       implicit none
!       real*8,intent(out) :: W0
!       real*8,intent(in)  :: X0,Z0
!       integer,intent(in) :: IWL,ILINE,ILEMKE
! !       real*8 :: WL,XT,XNE,PRF
! !       integer:: NWLH,NTH,NEH
! !       COMMON/VCSDAT/WL(36,8),XT(7,8),XNE(11,8),PRF(36,7,11),
! !      *              NWLH(8),NTH(8),NEH(8)
!       real*8 ::     XK,FXK,BETAD,DBETA
!       COMMON/AUXVCS/XK,FXK,BETAD,DBETA
! 
! C
! C     Interpolation in temperature and electron density from the
! C     Schoening and Butler tables for hydrogen lines to the actual valus of
! C     temperature and electron density
! C
!       integer :: NX,NZ,NT,NE,IZZ,IPZ,N0X,N1X,N0Z,N1Z,I0Z,IX,IPX,I0
!       real*8  :: BETA,A,DIV
!       INCLUDE '../INCLUDE2/PARAMS.FOR'
!       real*8,parameter :: TWO=2d0
!       real*8,dimension(3) :: ZZ,XX,WX,WZ
! C
!       NX=3
!       NZ=3
!       NT=NTH(ILINE)
!       NE=NEH(ILINE)
!       BETA=WL(IWL,ILINE)/FXK
!       IF(ILEMKE.EQ.1) THEN
!          BETA=WL(IWL,ILINE)/XK
!          NX=2
!          NZ=2
!       END IF
! C
! C     for values lower than the lowest grid value of electron density
! C     the profiles are determined by the approximate expression
! C     (see STARKA); not by an extrapolation in the HYD tables which may
! C     be very inaccurate
! C
!       IF(Z0.LT.XNE(1,ILINE)*0.99.OR.Z0.GT.XNE(NE,ILINE)*1.01) THEN
!          !*** changed: add BETAD
!          CALL DIVSTR(BETAD,A,DIV)
!          W0=STARKA(BETA,BETAD,A,DIV,TWO)*DBETA
!          W0=LOG10(W0)
!          GO TO 500 ! RETURN
!       END IF
! C
! C     Otherwise, one interpolates (or extrapolates for higher than the
! C     highes grid value of electron density) in the HYD tables
! C
!       DO 10 IZZ=1,NE-1
!          IPZ=IZZ
!          IF(Z0.LE.XNE(IZZ+1,ILINE)) GO TO 20
!    10 ENDDO
!    20 N0Z=IPZ-NZ/2+1
!       IF(N0Z.LT.1) N0Z=1
!       IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
!       N1Z=N0Z+NZ-1
! C
!       DO 300 IZZ=N0Z,N1Z
!          I0Z=IZZ-N0Z+1
!          ZZ(I0Z)=XNE(IZZ,ILINE)
! C
! C     Likewise, the approximate expression instead of extrapolation
! C     is used for higher that the highest grid value of temperature,
! C     if the Doppler width expressed in beta units (BETAD) is
! C     sufficiently large (> 10)
! C
!          IF(X0.GT.1.01*XT(NT,ILINE).AND.BETAD.GT.10.) THEN
!             CALL DIVSTR(BETAD,A,DIV)
!             W0=STARKA(BETA,BETAD,A,DIV,TWO)*DBETA
!             W0=LOG10(W0)
!             GO TO 500 ! RETURN
!          END IF
! C
! C     Otherwise, normal inter- or extrapolation
! C
! C     Both interpolations (in T as well as in electron density) are
! C     by default the quadratic interpolations in logarithms
! C
!          DO 30 IX=1,NT-1
!             IPX=IX
!             IF(X0.LE.XT(IX+1,ILINE)) GO TO 40
!    30    ENDDO
!    40    N0X=IPX-NX/2+1
!          IF(N0X.LT.1) N0X=1
!          IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
!          N1X=N0X+NX-1
!          DO 200 IX=N0X,N1X
!             I0=IX-N0X+1
!             XX(I0)=XT(IX,ILINE)
!             WX(I0)=PRF(IWL,IX,IZZ)
!   200    ENDDO
!          IF(WX(1).LT.-99..OR.WX(2).LT.-99..OR.WX(3).LT.-99.) THEN
!             CALL DIVSTR(BETAD,A,DIV)
!             W0=STARKA(BETA,BETAD,A,DIV,TWO)*DBETA
!             W0=LOG10(W0)
!             GO TO 500 ! RETURN
!           ELSE
!             WZ(I0Z)=YINT(XX,WX,X0)
!          END IF
!   300 ENDDO
!       W0=YINT(ZZ,WZ,Z0)
!   500 CONTINUE
!       RETURN
!       END SUBROUTINE
! C
! C ********************************************************************
! C
      END MODULE