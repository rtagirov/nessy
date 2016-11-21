      module MOD_SYNSUBM
      contains
      FUNCTION PHE1(ID,FREQ,ILINE)
      use constants
C     ============================
C
C     Absorption profile for four lines of He I, given by
C     Barnard, Cooper, Smith (1974) JQSRT 14, 1025 for the 4471 line;
C     Shamey (1969) PhD thesis, for other lines
C
C     Input: ID    - depth index
C            FREQ  - frequency
C            ILINE - index of the line ( = 1  for 4471,
C                                        = 2  for 4387,
C                                        = 3  for 4026,
C                                        = 4  for 4922)
C
C     Output: PHE1 - profile coefficient in frequency units,
C                    normalized to sqrt(pi) [not unity]
C
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID,ILINE
      real*8, intent(in   ) :: FREQ
      real*8 :: PHE1
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      COMMON/PROHE1/PRFHE1(50,4,8,3),DLMHE1(50,8,3),XNEHE1(8),
     *              NWLAM(8,4)
      COMMON/PRO447/PRF447(80,4,7),DLM447(80,7),XNE447(7)
      DIMENSION WLAM0(4),XT(4),XX(3),WX(3),YY(2),PP(2),ZZ(3),WZ(3)
      DATA WLAM0 / 4471.50, 4387.93, 4026.20, 4921.93/
      DATA XT / 3.699, 4.000, 4.301, 4.602/
C
      integer,parameter :: NX=3
      integer,parameter :: NZ=3
      integer,parameter :: NY=2
      integer,parameter :: NT=4
      T=TEMP(ID)
      TL=LOG10(T)
      ANE=ELEC(ID)
      ANEL=LOG10(ANE)
      !ALAM=2.99792458E18/FREQ
      ALAM=CLIGHT_SI*1e10/FREQ
      DLAM=ALAM-WLAM0(ILINE)
      DOPL=SQRT(4.125E7*T+VTURB(ID))*WLAM0(ILINE)/CLIGHT_CGS
C
      IF(ILINE.EQ.1.AND.ANEL.GE.XNE447(1)) GO TO 10
      IF(ILINE.NE.1.AND.ANEL.GE.XNEHE1(1)) GO TO 10
      !***
      !***  isolated line approximation for low electron densities
      !***
      A=WTOT(T,ANE,ID,ILINE)/DOPL
      V=ABS(DLAM)/DOPL
      V1=ABS(ALAM-4471.682)/DOPL
      PHE1=VOIGTK(A,V)
      IF(ILINE.EQ.1) PHE1=(8.*PHE1+VOIGTK(A,V1))/9.
      RETURN
      !***
      !***  otherwise, interpolation (or extrapolation) in tables
      !***
   10 continue 
      !NX=3
      !NZ=3
      !NY=2
      !NT=4
      NE=8
      ILNE=ILINE-1
      IF(ILINE.EQ.1) NE=7
      !***
      !***  Interpolation in electron density
      !***
      !** Find the table that fits the #electrons 
      !**   (assumes a monotonic increasing table in #electrons)
      DO JZ=1,NE-1
        IPZ=JZ
	  IF(ILINE.EQ.1.AND.ANEL.LE.XNE447(JZ+1)) GO TO 30
	  IF(ILINE.NE.1.AND.ANEL.LE.XNEHE1(JZ+1)) GO TO 30
      ENDDO
   30 N0Z=IPZ-NZ/2+1  !N0Z is the Table index-0.5
      !** Force N0Z into [1, NE-2]
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
      N1Z=N0Z+NZ-1

      DO 300 JZ=N0Z,N1Z
        I0Z=JZ-N0Z+1
        IF(ILINE.EQ.1) ZZ(I0Z)=XNE447(JZ) !XNE447 = #electrons
        IF(ILINE.NE.1) ZZ(I0Z)=XNEHE1(JZ)
        !***
        !***  Interpolation in temperature
        !***
        DO IX=1,NT-1
          IPX=IX
          IF(TL.LE.XT(IX+1)) GO TO 50
        ENDDO
   50   N0X=IPX-NX/2+1     ! Temperature Index (in [1,2])
        IF(N0X.LT.1) N0X=1
        IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
        N1X=N0X+NX-1
        DO 200 IX=N0X,N1X
          I0X=IX-N0X+1
          XX(I0X)=XT(IX)
          !***
          !***  Interpolation in wavelength
          !***
          !***  1. For delta lambda beyond tabulated values - special
          !***   extrapolation (Cooper's suggestion)
          !***
          NLST=NWLAM(JZ,ILINE)
          IF(ILINE.EQ.1) THEN
            D1=DLM447(1,JZ)
            D2=DLM447(NLST,JZ)
            IF(DLAM.LT.D1) THEN
              PRF=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D1,PRF447(1,IX,JZ))
              GO TO 150
            ELSE IF(DLAM.GT.D2) THEN
              PRF=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D2,
     *                   PRF447(NLST,IX,JZ))
              GO TO 150
            END IF
          ELSE
            D1=DLMHE1(1,JZ,ILNE)
            D2=DLMHE1(NLST,JZ,ILNE)
            IF(DLAM.LT.D1) THEN
              PRF=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D1,
     *                        PRFHE1(1,IX,JZ,ILNE))
              GO TO 150
            ELSE IF(DLAM.GT.D2) THEN
              PRF=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D2,
     *                   PRFHE1(NLST,IX,JZ,ILNE))
              GO TO 150
            END IF
          END IF
          !***
          !***  normal linear interpolation in wavelength
          !***  (for 4471, linear interpolation in logarithms)
          !***
          MAX_LAM_LT:DO IY=1,NLST-1
            IPY=IY
            IF(ILINE.EQ.1) THEN
              IF(DLAM.LE.DLM447(IY+1,JZ))  EXIT MAX_LAM_LT
            ELSE
              IF(DLAM.LE.DLMHE1(IY+1,JZ,ILNE)) EXIT MAX_LAM_LT
            ENDIF
          ENDDO MAX_LAM_LT
          N0Y=IPY-NY/2+1
          IF(N0Y.LT.1) N0Y=1
          IF(N0Y.GT.NLST-NY+1) N0Y=NLST-NY+1
          N1Y=N0Y+NY-1
          DO IY=N0Y,N1Y
            I0=IY-N0Y+1
            IF(ILINE.EQ.1) THEN
              YY(I0)=DLM447(IY,JZ)
              PP(I0)=LOG(PRF447(IY,IX,JZ))
            ELSE
              YY(I0)=DLMHE1(IY,JZ,ILNE)
              PP(I0)=PRFHE1(IY,IX,JZ,ILNE)
            ENDIF
          ENDDO
          IF(ILINE.NE.1) THEN
            WX(I0X)=(PP(2)*(DLAM-YY(1))+PP(1)*(YY(2)-DLAM))/
     *                (YY(2)-YY(1))
          ELSE
            ! IF(YY(1).NE.0.) THEN
            !   YY(1)=LOG(ABS(YY(1)))
            ! ELSE
            !   YY(1)=-8.
            ! END IF
            ! IF(YY(2).NE.0.) THEN
            !   YY(2)=LOG(ABS(YY(2)))
            ! ELSE
            !   YY(2)=-8.
            ! END IF
            ! IF(DLAM.NE.0.) THEN
            !    DLAM=LOG(ABS(DLAM))
            ! ELSE
            !   DLAM=-8.
            ! END IF
             WX(I0X)=(PP(2)*(DLAM-YY(1))+PP(1)*(YY(2)-DLAM))/
     *                (YY(2)-YY(1))
             WX(I0X)=EXP(WX(I0X))
	   END IF
	   GO TO 200
  150      WX(I0X)=PRF
  200   ENDDO
	WZ(I0Z)=YINT(XX,WX,TL)
  300 ENDDO
      W0=YINT(ZZ,WZ,ANEL)
      PHE1=W0*DOPL*1.772454
      RETURN
      END function
C
C ********************************************************************
C ********************************************************************
C
C

      SUBROUTINE PHE2(ISPEC,ID,ABLIN,EMLIN)
C
C     Evaluation of the opacity and emissivity in a given He II line,
C     using profile coefficients calculated by Schoening and Butler.
C
C     Input: ISPEC - line index, defined in HE2INI
C            ID    - depth index
C     Output: ABLIN - absorption coefficient
C             EMLIN - emission coefficient
C
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ISPEC,ID
      real*8, intent(  out) :: ABLIN,EMLIN
        
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      
      DIMENSION ABLIN(1),EMLIN(1),OSCHE2(19),PRF0(99),WLL(40)
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,
     *        1.604E-2, 1.023E-2, 6.980E-3,
     *        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,
     *        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,
     *        1.038,    1.793E-1, 6.549E-2/
C
C     ILINE - line index
C
      ILINE=ISPEC-5
C
      DO 10 IWL=1,NWLHE2(ILINE)
         PRF0(IWL)=PRFHE2(ILINE,ID,IWL)
         WLL(IWL)=WLHE2(ILINE,IWL)
   10 CONTINUE
C
      I=ILHE2(ILINE)
      J=IUHE2(ILINE)
      II=I*I
      JJ=J*J
      IF(I.LE.2) THEN
         WLINE=227.838/(1./II-1./JJ)
       ELSE
         WLINE=227.7776/(1./II-1./JJ)
      END IF
      T=TEMP(ID)
C
C     He III population (either LTE or NLTE, depending on input model)
C
      IF(IELHE2.GT.0.and.inlte.gt.0) THEN
         PP=POPUL(NNEXT(IELHE2),ID)
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1
       ELSE
         PP=RRR(ID,3,2)
         NLHE2=0
      END IF
C
C     population of the lower level of the given transition
C     (again either LTE or NLTE)
C
      PP=PP*ELEC(ID)*4.1412E-16/T/SQRT(T)*II
      IF(I.LE.NLHE2.and.inlte.gt.0) THEN
         POPI=POPUL(NFIRST(IELHE2)+I-1,ID)
       ELSE
         POPI=PP*EXP(631479./T/II)
      END IF
C
C     population of the upper level of the given transition
C     (again either LTE or NLTE)
C
      IF(J.LE.NLHE2) THEN
         POPJ=POPUL(NFIRST(IELHE2)+J-1,ID)*II/JJ
       ELSE
         POPJ=PP*EXP(631479./T/JJ)
      END IF

C
C     loop over frequency points - opacity and emissivity in the given line
C     absorption coefficent is found by interpolating in previously
C     calculated tables, based on calculations of Schoening and Butler
C     (see procedure HE2INI)
C
      FID=0.02654*OSCHE2(ILINE)
      DO 50 IJ=1,NFREQ
         AL=ABS(WLAM(IJ)-WLINE)
         IF(AL.LT.1.E-4) AL=1.E-4
         AL=LOG10(AL)
         DO 20 IWL=1,NWLHE2(ILINE)-1
            IW0=IWL
            IF(AL.LE.WLL(IWL+1)) GO TO 30
   20    CONTINUE
   30    IW1=IW0+1
         PRF=(PRF0(IW0)*(WLL(IW1)-AL)+PRF0(IW1)*(AL-WLL(IW0)))/
     *       (WLL(IW1)-WLL(IW0))
         SG=EXP(PRF*2.3025851)*FID
         ABLIN(IJ)=ABLIN(IJ)+SG*(POPI-POPJ)
         EMLIN(IJ)=EMLIN(IJ)+SG*POPJ*1.4747E-2*(FREQ(IJ)*1.E-15)**3
   50 CONTINUE
      RETURN
      END subroutine
      
C
C
C
C     ****************************************************************
C
C

      SUBROUTINE CARBON(IB,FR,SG)
C     ===========================
C
C     Photoionization cross-section for neutral carbon 2p1D and 2p1S
C     levels (G.B.Taylor - private communication)
C
      INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IB
      real*8, intent(in   ) :: FR
      real*8, intent(inout) :: SG
      
        
      DIMENSION FR2(34),SG2(34),FR3(45),SG3(45)
      DATA FR2/ 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82,
     *                0.83,       0.85, 0.86, 0.87, 0.88, 0.89, 0.90,
     *          0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
     *          1.00, 1.10, 1.20, 1.30, 1.45, 1.50, 1.60, 1.80, 2./
      DATA SG2/ 12.04, 12.03, 12.09, 12.26, 12.60, 13.24, 14.36, 16.24,
     *          19.28, 23.94, 37.41, 42.88, 44.76, 43.41, 40.46, 37.19,
     *          34.26, 31.82, 29.96, 28.57, 27.68, 27.37, 27.84, 29.69,
     *          34.45, 46.35, 13.80, 11.54, 10.40,  8.96,  8.54,  7.47,
     *           6.53,  5.66/
      DATA FR3/ 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82,
     *          0.84, 0.86, 0.864,0.866,0.868,0.87, 0.874,0.876,0.88,
     *          0.882,0.884,0.886,0.888,0.89 ,0.894,0.896,0.898,0.90,
     *          0.904,0.908,0.910,0.920,0.94, 0.98, 1.00, 1.10, 1.20,
     *          1.26, 1.34, 1.36, 1.40, 1.46, 1.60, 1.70, 1.80, 2./
      DATA SG3/ 13.94, 13.29, 12.56, 11.73, 10.82, 10.18,  8.62,  7.27,
     *           5.74,  4.14,  4.61,  5.92,  6.94,  8.34, 10.21, 16.12,
     *          20.64, 34.56, 44.82, 57.71, 73.09, 89.99,106.38,127.08,
     *         128.38,124.44,117.17, 99.32, 82.95, 76.05, 52.65, 33.23,
     *          21.29, 18.69, 12.62, 11.44,  9.77,  7.53, 10.47,  9.65,
     *          10.19,  7.28,  6.70,  6.11,  4.96/
      DATA NC2,NC3/34,45/
      DATA FR0/3.28805E15/
      F=FR/FR0
      IF(IB.NE.-602) GO TO 25
      J=2
      IF(F.LE.FR2(1)) GO TO 20
      DO 10 I=2,NC2
         J=I
         IF(F.GT.FR2(I-1).AND.F.LE.FR2(I)) GO TO 20
   10 CONTINUE
   20 SG=(F-FR2(J-1))/(FR2(J)-FR2(J-1))*(SG2(J)-SG2(J-1))+SG2(J-1)
      SG=SG*1.E-18
   25 IF(IB.NE.-603) GO TO 50
      J=2
      IF(F.LE.FR3(1)) GO TO 40
      DO 30 I=2,NC3
         J=I
         IF(F.GT.FR3(I-1).AND.F.LE.FR3(I)) GO TO 40
   30 CONTINUE
   40 SG=(F-FR3(J-1))/(FR3(J)-FR3(J-1))*(SG3(J)-SG3(J-1))+SG3(J-1)
      SG=SG*1.E-18
   50 CONTINUE
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
cmh      SUBROUTINE CROSET
      SUBROUTINE CROSET(WAVARR,SIGARR,NDIM,NFDIM)
      use SYNTHP_CONT,only:FREQC,NFCONT
C     called by INIBL0(WAVARR,SIGARR)
C     driver for evaluating photoinization cross-sections
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'SYNTHP.FOR'
      use SYNTHP_CONT
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: NDIM,NFDIM
      real*8, intent(in   ) :: WAVARR,SIGARR
      
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)
      common/phopar/cross(mcross,mfcont)
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
     $ IOPFE1
C
!      print*, 'HM test', IOPHMI
!      stop

      CROSS(:,:)=0.
      DO 20 IT=1,NLEVE0
         DO 20 IJ=1,NFCONT()
	         FR=FREQC(IJ)
c	CROSS(IT,IJ)=SIGK(FR,IT,0)
C***  CHANGED BY MARGIT HABERREITER *****************************
cmh	CROSS(IT,IJ)=SIGK(FR,IT,0,WAVARR,SIGARR)
CMH	mode changed >0
	CROSS(IT,IJ)=SIGK(FR,IT,0,WAVARR,SIGARR,NDIM,NFDIM)
C****************************************************************
c***	PRINT *,'CROSET AFTER SIGK',IT,IJ,CROSS(IT,IJ)
   20 CONTINUE
C
C   PARAMETERS FOR CALCULATING ADDITIONAL OPACITIES
C
      IF(IOPADD.EQ.0) GO TO 40
      print *,'mode=-1 for OPADD'
      DO 30 IJ=1,NFCONT()
         FR=FREQC(IJ)
C***	Margit Haberreiter
CMH   MODE FOR OPADD EITHER 0 OR -1
C            CALL OPADD(0,IJ,1,FR,X1,X2,X3)
	       CALL OPADD(-1,IJ,1,FR,X1,X2,X3)
   30 CONTINUE
   40 CONTINUE
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C

      SUBROUTINE DIVHE2(BETAD,A,DIV)
C
C     Auxiliary procedure for evaluating approximate Stark profile
C     for He II lines
C     This procedure is quite analogous to DIVSTR for hydrogen;
C     the only difference is a somewhat different definition
C     of the parameter A ,ie. A for He II is equal to A for hydrogen
C     minus ln(2)
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: BETAD
      real*8, intent(inout) :: A,DIV
      PARAMETER (UN=1.,TWO=2.,UNQ=1.25,UNH=1.5,TWH=2.5,FO=4.,FI=5.)
      PARAMETER (CA=0.978,BL=5.821,AL=1.26,CX=0.28,DX=0.0001)
C
      A=UNH*LOG(BETAD)-CA
      IF(BETAD.LT.BL) RETURN
      IF(A.GE.AL) THEN
         X=SQRT(A)*(UN+UNQ*LOG(A)/(FO*A-FI))
      ELSE
         X=SQRT(CX+A)
      ENDIF
      DO 10 I=1,5
         XN=X*(UN-(X*X-TWH*LOG(X)-A)/(TWO*X*X-TWH))
         IF(ABS(XN-X).LE.DX) GO TO 20
         X=XN
   10 CONTINUE
   20 DIV=X
      RETURN
      END SUBROUTINE
C
C *******************************************************************
C *******************************************************************
C
      PURE SUBROUTINE DIVSTR(BETAD,A,DIV)
C     ==============================
C
C     Auxiliary procedure for STARKA - determination of the division
C     point between Doppler and asymptotic Stark profiles
C
C     Input:  BETAD - Doppler width in beta units
C     Output: A     - auxiliary parameter
C                     A=1.5*LOG(BETAD)-1.671
C             DIV   - only for A > 1; division point between Doppler
C                     and asymptotic Stark wing, expressed in units
C                     of betad.
C                     DIV = solution of equation
C                     exp(-(beta/betad)**2)/betad/sqrt(pi)=3*beta**-5/2
C
c      INCLUDE 'IMPLIC.FOR'
      INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: BETAD
      real*8, intent(  out) :: A,DIV
      
        
      PARAMETER (UN=1.,TWO=2.,UNQ=1.25,UNH=1.5,TWH=2.5,FO=4.,FI=5.)
      PARAMETER (CA=1.671,BL=5.821,AL=1.26,CX=0.28,DX=0.0001)
C
      A=UNH*LOG(BETAD)-CA
      IF(BETAD.LT.BL) RETURN
      IF(A.GE.AL) THEN
         X=SQRT(A)*(UN+UNQ*LOG(A)/(FO*A-FI))
      ELSE
         X=SQRT(CX+A)
      ENDIF
      !* Newton Rhapson Method to solve the equation
      DO I=1,5
         XN=X*(UN-(X*X-TWH*LOG(X)-A)/(TWO*X*X-TWH))
         IF(ABS(XN-X).LE.DX) GO TO 20
         X=XN
      ENDDO
   20 DIV=X
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
      PURE FUNCTION EPS(T,ANE,ALAM,ION,N)
C
C   NLTE PARAMETER EPSILON (COLLISIONAL/SPONTANEOUS DEEXCITATION)
C   AFTER  KASTNER, 1981, J.Q.S.R.T. 26, 377
C
c     INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: ION,N
      real*8, intent(in   ) :: T,ANE,ALAM
      real*8 :: EPS
      real*8,parameter :: CK0= 7.75E-8,CK1 = 2.58E-8
      X=1.438E8/ALAM/T
      XKT=12390./ALAM
      TT=0.75*X
      T1=TT+1.
      A=4.36E7*XKT*XKT/(1.-EXP(-X))
      IF(ION.EQ.1) GO TO 10
      B=1.1+LOG(T1/TT)-0.4/T1/T1
      C=X*B*SQRT(T)/XKT/XKT*ANE
      IF(N.EQ.0) C=CK0*C
      IF(N.NE.0) C=CK1*C
      GO TO 20
   10 C=2.16/T/SQRT(T)/X**1.68*ANE
   20 EPS=C/(C+A)
      RETURN
      END FUNCTION

C ********************************************************************
C ********************************************************************
C

      PURE FUNCTION EXTPRF(DLAM,IT,ILINE,ANEL,DLAST,PLAST)
C     ===============================================
C
C     Extrapolation in wavelengths in Shamey, or Barnard, Cooper,
C     Smith tables
C     Special formula suggested by Cooper
C
c     INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IT,ILINE
      real*8, intent(in   ) :: DLAM,ANEL,DLAST,PLAST
      real*8 :: EXTPRF
      real*8,parameter::W0(4,4)=(/ 1.460, 1.269, 1.079, 0.898,
     *             6.130, 5.150, 4.240, 3.450,
     *             4.040, 3.490, 2.960, 2.470,
     *             2.312, 1.963, 1.624, 1.315/)
C
      WE=W0(IT,ILINE)*EXP(ANEL*2.3025851)*1.E-16
      DLASTA=ABS(DLAST)
      D52=DLASTA*DLASTA*SQRT(DLASTA)
      F=D52*(PLAST-WE/3.14159/DLAST/DLAST)
      EXTPRF=(WE/3.14159+F/SQRT(ABS(DLAM)))/DLAM/DLAM
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION FEAUTR(FREQ,ID)
      use constants
C
C     LYMAN-ALPHA STARK BROADENING AFTER N.FEAUTRIER
C
c     INCLUDE 'PARAMS.FOR'
c     INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID
      real*8, intent(in   ) :: FREQ
      real*8  :: FEAUTR
        
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      DIMENSION DL(20),F05(20),F10(20),F20(20),F40(20),X(4)
      DATA F05 / 0.0537, 0.0964, 0.1330, 0.3105, 0.4585, 0.6772, 0.8229,
     *           0.8556, 0.9250, 0.9618, 0.9733, 1.1076, 1.0644, 1.0525,
     *           0.8841, 0.8282, 0.7541, 0.7091, 0.7164, 0.7672/
      DATA F10 / 0.1986, 0.2764, 0.3959, 0.5740, 0.7385, 0.9448, 1.0292,
     *           1.0317, 0.9947, 0.8679, 0.8648, 0.9815, 1.0660, 1.0793,
     *           1.0699, 1.0357, 0.9245, 0.8603, 0.8195, 0.7928/
      DATA F20 / 0.4843, 0.5821, 0.7003, 0.8411, 0.9405, 1.0300, 1.0029,
     *           0.9753, 0.8478, 0.6851, 0.6861, 0.8554, 0.9916, 1.0264,
     *           1.0592, 1.0817, 1.0575, 1.0152, 0.9761, 0.9451/
      DATA F40 / 0.7862, 0.8566, 0.9290, 0.9915, 1.0066, 0.9878, 0.8983,
     *           0.8513, 0.6881, 0.5277, 0.5302, 0.6920, 0.8607, 0.9111,
     *           0.9651, 1.0793, 1.1108, 1.1156, 1.1003, 1.0839/
      DATA DL / -150., -120., -90., -60., -40., -20., -10., -8., -4.,
     *          -2., 2., 4., 8., 10., 20., 40., 60., 90., 120., 150./
      DLAM=CLIGHT_SI*1e10/FREQ-1215.685
      DO 10 I=2,20
         IF(DLAM.LE.DL(I)) GO TO 20
   10 CONTINUE
      I=20
   20 J=I-1
      C=DL(J)-DL(I)
      A=(DLAM-DL(I))/C
      B=(DL(J)-DLAM)/C
      X(1)=F05(J)*A+F05(I)*B
      X(2)=F10(J)*A+F10(I)*B
      X(3)=F20(J)*A+F20(I)*B
      X(4)=F40(J)*A+F40(I)*B
      J=JT(ID)
      Y=TI0(ID)*X(J)+TI1(ID)*X(J-1)+TI2(ID)*X(J-2)
      FEAUTR=0.5*(Y+1.)
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE GAMHE(IND,T,ANE,ANP,ID,GAM)
C
C     NEUTRAL HELIUM STARK BROADENING PARAMETERS
C     AFTER DIMITRIJEVIC AND SAHAL-BRECHOT, 1984, J.Q.S.R.T. 31, 301
C     OR FREUDENSTEIN AND COOPER, 1978, AP.J. 224, 1079  (FOR C(IND).GT.0)
C
c     INCLUDE 'PARAMS.FOR'
c     INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: IND,ID
      real*8, intent(in   ) :: T,ANE,ANP
      real*8, intent(  out) :: GAM
        
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      DIMENSION W(5,20),V(4,20),C(20)
C
C   ELECTRONS T= 5000   10000   20000   40000     LAMBDA
C
      DATA W /  5.990,  6.650,  6.610,  6.210,    3819.60,
     *          2.950,  3.130,  3.230,  3.300,    3867.50,
     *          0.000,  0.000,  0.000,  0.000,    3871.79,
     *          0.142,  0.166,  0.182,  0.190,    3888.65,
     *          0.000,  0.000,  0.000,  0.000,    3926.53,
     *          1.540,  1.480,  1.400,  1.290,    3964.73,
     *         41.600, 50.500, 57.400, 65.800,    4009.27,
     *          1.320,  1.350,  1.380,  1.460,    4120.80,
     *          7.830,  8.750,  8.690,  8.040,    4143.76,
     *          5.830,  6.370,  6.820,  6.990,    4168.97,
     *          0.000,  0.000,  0.000,  0.000,    4437.55,
     *          1.630,  1.610,  1.490,  1.350,    4471.50,
     *          0.588,  0.620,  0.641,  0.659,    4713.20,
     *          2.600,  2.480,  2.240,  1.960,    4921.93,
     *          0.627,  0.597,  0.568,  0.532,    5015.68,
     *          1.050,  1.090,  1.110,  1.140,    5047.74,
     *          0.277,  0.298,  0.296,  0.293,    5875.70,
     *          0.714,  0.666,  0.602,  0.538,    6678.15,
     *          3.490,  3.630,  3.470,  3.190,    4026.20,
     *          4.970,  5.100,  4.810,  4.310,    4387.93/
C
C   PROTONS   T= 5000   10000   20000   40000
C
      DATA V /  1.520,  4.540,  9.140, 10.200,
     *          0.607,  0.710,  0.802,  0.901,
     *          0.000,  0.000,  0.000,  0.000,
     *          0.0396, 0.0434, 0.0476, 0.0526,
     *          0.000,  0.000,  0.000,  0.000,
     *          0.507,  0.585,  0.665,  0.762,
     *          0.930,  1.710, 13.600, 27.200,
     *          0.288,  0.325,  0.365,  0.410,
     *          1.330,  6.800, 12.900, 14.300,
     *          1.100,  1.370,  1.560,  1.760,
     *          0.000,  0.000,  0.000,  0.000,
     *          1.340,  1.690,  1.820,  1.630,
     *          0.128,  0.143,  0.161,  0.181,
     *          2.040,  2.740,  2.950,  2.740,
     *          0.187,  0.210,  0.237,  0.270,
     *          0.231,  0.260,  0.291,  0.327,
     *          0.0591, 0.0650, 0.0719, 0.0799,
     *          0.231,  0.260,  0.295,  0.339,
     *          2.180,  3.760,  4.790,  4.560,
     *          1.860,  5.320,  7.070,  7.150/
      DATA C /2*0.,1.83E-4,0.,1.13E-4,5*0.,1.6E-4,9*0./
C
      IF(W(1,IND).EQ.0.) GO TO 10
      J=JT(ID)
      GAM=((TI0(ID)*W(J,IND)+TI1(ID)*W(J-1,IND)+TI2(ID)*W(J-2,IND))
     *     *ANE
     *    +(TI0(ID)*V(J,IND)+TI1(ID)*V(J-1,IND)+TI2(ID)*V(J-2,IND))
     *     *ANP)*1.884E3/W(5,IND)**2
      IF(GAM.LT.0.) GAM=0.
      RETURN
   10 GAM=C(IND)*T**0.16667*ANE
      RETURN
      END SUBROUTINE
C
C
C     ****************************************************************
C
C




      FUNCTION GAUNT_W(I,FR)
      use constants
C     ====================                                                                                                                                                  C
        INCLUDE '../inc/IMPLIC.FOR'

      integer :: I
      real*8  :: FR
      real*8 :: a1(6), a2(6), b1(6), b2(6), c1(6), c2(6)
      real*8 :: GAUNT_W
 
      a1=(/1.0780, 1.0926, 1.0983, 1.0954, 1.0912, 1.0925/)
      b1=(/-8.754d14,-2.019d14, -9.450d13,-5.188d13,-3.200d13, -2.331d13/)
      c1=(/-1.791d29, 1.836d28, 9.177d27, 3.552d27, 1.576d27, 9.325d26/)

      a2=(/0.798,0.768,0.793, 0.831, 0.758, 0.790/)
      b2=(/5.358d15, 6.242d15,5.480d15,4.094d15,6.633d15,5.808d15/)
      c2=(/-3.484d31,-3.208d31,-2.318d31,-1.430d31,-3.320d31,-2.844d31/)
   !   print*, 'test1'
      if (i .gt. 6) then 
         GAUNT_W=GAUNT(I, FR)
 
         return
      else
         
      GAUNT_W=a1(i)+b1(i)/FR+c1(i)/(FR*FR)   
      return
      endif
     
      END FUNCTION


      FUNCTION GAUNT(I,FR)
      use constants
C     ====================
C
C     Hydrogenic bound-free Gaunt factor for the principal quantum
C     number I and frequency FR
C
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: I
      real*8, intent(in   ) :: FR
      real*8  ::  GAUNT
        
      X=FR/(CLIGHT_SI*1e6)
      GAUNT=1.
      IF(I.GT.10) GO TO 16
      GO TO (1,2,3,4,5,6,7,8,9,10),I
    1 GAUNT=1.2302628+X*(-2.9094219E-3+X*(7.3993579E-6-8.7356966E-9*X))
     *+(12.803223/X-5.5759888)/X
      GO TO 16
    2 GAUNT=1.1595421+X*(-2.0735860E-3+2.7033384E-6*X)+(-1.2709045+
     *(-2.0244141/X+2.1325684)/X)/X
      GO TO 16
    3 GAUNT=1.1450949+X*(-1.9366592E-3+2.3572356E-6*X)+(-0.55936432+
     *(-0.23387146/X+0.52471924)/X)/X
      GO TO 16
    4 GAUNT=1.1306695+X*(-1.3482273E-3+X*(-4.6949424E-6+2.3548636E-8*X))
     *+(-0.31190730+(0.19683564-5.4418565E-2/X)/X)/X
      GO TO 16
    5 GAUNT=1.1190904+X*(-1.0401085E-3+X*(-6.9943488E-6+2.8496742E-8*X))
     *+(-0.16051018+(5.5545091E-2-8.9182854E-3/X)/X)/X
      GO TO 16
    6 GAUNT=1.1168376+X*(-8.9466573E-4+X*(-8.8393133E-6+3.4696768E-8*X))
     *+(-0.13075417+(4.1921183E-2-5.5303574E-3/X)/X)/X
      GO TO 16
    7 GAUNT=1.1128632+X*(-7.4833260E-4+X*(-1.0244504E-5+3.8595771E-8*X))
     *+(-9.5441161E-2+(2.3350812E-2-2.2752881E-3/X)/X)/X
      GO TO 16
    8 GAUNT=1.1093137+X*(-6.2619148E-4+X*(-1.1342068E-5+4.1477731E-8*X))
     *+(-7.1010560E-2+(1.3298411E-2 -9.7200274E-4/X)/X)/X
      GO TO 16
    9 GAUNT=1.1078717+X*(-5.4837392E-4+X*(-1.2157943E-5+4.3796716E-8*X))
     *+(-5.6046560E-2+(8.5139736E-3-4.9576163E-4/X)/X)/X
      GO TO 16
   10 GAUNT=1.1052734+X*(-4.4341570E-4+X*(-1.3235905E-5+4.7003140E-8*X))
     *+(-4.7326370E-2+(6.1516856E-3-2.9467046E-4/X)/X)/X
   16 RETURN
      END FUNCTION
C
C
C     ****************************************************************
C
      FUNCTION GFREE(T,FR)
      use constants
C     ====================
C
C     Hydrogenic free-free Gaunt factor, for temperature T and
C     frequency FR
C
c     INCLUDE 'IMPLIC.FOR'
      INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: T,FR
      real*8  :: GFREE
        
      THET=5040.4/T
      IF(THET.LT.4.E-2) THET=4.E-2
      X=FR/(CLIGHT_SI*1e6)
      IF(X.GT.1) GO TO 10
      IF(X.LT.0.2) X=0.2
      GFREE=(1.0823+2.98E-2/THET)+(6.7E-3+1.12E-2/THET)/X
      RETURN
   10 C1=(3.9999187E-3-7.8622889E-5/THET)/THET+1.070192
      C2=(6.4628601E-2-6.1953813E-4/THET)/THET+2.6061249E-1
      C3=(1.3983474E-5/THET+3.7542343E-2)/THET+5.7917786E-1
      C4=3.4169006E-1+1.1852264E-2/THET
      GFREE=((C4/X-C3)/X+C2)/X+C1
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE GRIEM(ID,ANE,ION,FR,WGR,GAM)
C
C     STARK DAMPING PARAMETER (GAM) CALCULATED FROM INPUT VALUES
C     OF STARK WIDTHS FOR  T=5000, 10000, 20000, 40000 K,
C     AND FOR  NE=1.E16 (FOR NEUTRALS)  OR  NE = 1.E17 (FOR IONS)
C

      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID,ION
      real*8, intent(in   ) :: ANE,FR,WGR
      real*8, intent(  out) :: GAM
        
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      DIMENSION WGR(4)
      J=JT(ID)
      GAM=(TI0(ID)*WGR(J)+TI1(ID)*WGR(J-1)+TI2(ID)*WGR(J-2))
     *    *ANE*1.E-10*FR*1.E-10*FR*4.2E-14
      IF(ION.GT.1) GAM=GAM*0.1
      IF(GAM.LT.0.) GAM=0.
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C

      SUBROUTINE HE1INI(MODE,IH)
      !==============================
      !***
      !***  Initializes necessary arrays for evaluating the He I line
      !***  absorption profiles using data calculated by Barnard, Cooper
      !***  and Smith JQSRT 14, 1025, 1974 (for 4471)
      !***  or Shamey, unpublished PhD thesis, 1969 (for other lines)
      !***
      !***  This procedure is quite analogous to BALINI for Balmer lines
      !***
      !implicit real*8(a-h,o-z)
      implicit none
      integer,intent(in   ) :: MODE,IH
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      !*** PRFHE1(Wavelength,Temperature,Electrondensity,Line Number)
      real*8 :: PRFHE1,DLMHE1,XNEHE1
      integer:: NWLAM
      COMMON/PROHE1/PRFHE1(50,4,8,3),DLMHE1(50,8,3),XNEHE1(8),
     *              NWLAM(8,4)
      !*** PRF447(Wavelength, Temperature, Electrondensity)
      !*** DLM447(Wavelength, Electrondensity)
      !*** XNE447(Electrondensity)
      real*8 :: PRF447,DLM447,XNE447
      COMMON/PRO447/PRF447(80,4,7),DLM447(80,7),XNE447(7)
      integer,parameter :: NT = 4
      integer :: NE, IE, IL, IT,ILN,IE1,I
      integer :: NWL ! #wavelength lines
      real*8 :: WL0,XNE
      !
      ! OPEN(UNIT=IH,STATUS='OLD',READONLY)
      !
      IF(MODE.EQ.1) THEN
      !***
      !***     read the Barnard, Cooper, Smith tables for He I 4471 line,
      !***     which have to be stored in file unit IH
      !***
        NE=7
        MODE1_readtable: DO IE=1,NE
          ! read header line
          ! *** Line #, Freqpoint,  Table#,  #electrons, #wavelength
          READ(IH,501) IL,WL0,IE1,XNE,NWL
          NWLAM(IE,1)=NWL
          XNE447(IE)=LOG10(XNE)
          !*** Read in the
          DO I=1,NWL
            READ(IH,502) DLM447(I,IE), (PRF447(I,IT,IE),IT=1,NT)
          ENDDO
        ENDDO MODE1_readtable
      ELSE
      !***
      !***     read Shamey's tables for He I 4387, 4026, and 4922 lines
      !***     which have to be stored in file unit IH
      !***
        NE=8
        DO ILN=1,3
          DO IE=1,NE
            READ(IH,501) IL,WL0,IE1,XNE,NWL
            NWLAM(IE,ILN+1)=NWL
            XNEHE1(IE)=LOG10(XNE)
            DO I=1,NWL
              READ(IH,*) DLMHE1(I,IE,ILN), (PRFHE1(I,IT,IE,ILN),IT=1,NT)
            ENDDO
          ENDDO
        ENDDO
      END IF
C
  501 FORMAT(/9X,I2,7X,F10.3,13X,I2,6X,E8.1,7X,I3/)
  502 FORMAT(5E10.2)
      RETURN
      END SUBROUTINE

C
C ********************************************************************
C ********************************************************************
C

      SUBROUTINE HE2INI(MODE,IH)
      use constants
C
C     Initializes necessary arrays for evaluating the He II line
C     absorption profiles using data calculated by Schoening and
C     Butler
C
C     This procedure is quite analogous to BALINI for Balmer lines
C
c     INCLUDE 'PARAMS.FOR'
c     INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: MODE,IH
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      DIMENSION NLINE0(3),NLINE1(3)
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
      COMMON/HE2DAT/WLHE(36,19),XTHE(6),XNEHE(11,19),PRFHE(36,6,11),
     *              NWLHE,NTHE,NEHE
      COMMON/AUXHE2/XKE,FXKE,BETADE,DBETAE
      DATA NLINE0 /1,8,17/,
     *     NLINE1 /7,16,19/
C
c      OPEN(UNIT=IH,STATUS='OLD',READONLY)
C
      DO 100 ILINE=NLINE0(MODE),NLINE1(MODE)
C
C     read the Schoening and Butler tables, which have to be stored
C     in file unit IH
C
         READ(IH,501) ILHE2(ILINE),IUHE2(ILINE)
         IF(ILHE2(ILINE).LE.2) THEN
            WL00=227.838
          ELSE
            WL00=227.7776
         END IF
         WL0=WL00/(1./ILHE2(ILINE)**2-1./IUHE2(ILINE)**2)
         READ(IH,*) NWLHE,(WLHE(I,ILINE),I=1,NWLHE)
         READ(IH,503) NTHE,(XTHE(I),I=1,NTHE)
         READ(IH,504) NEHE,(XNEHE(I,ILINE),I=1,NEHE)
         READ(IH,500)
         NWLHE2(ILINE)=NWLHE
C
         DO 10 I=1,NWLHE
            IF(WLHE(I,ILINE).LT.1.E-4) WLHE(I,ILINE)=1.E-4
            WLHE2(ILINE,I)=LOG10(WLHE(I,ILINE))
   10    CONTINUE
C
         DO 20 IE=1,NEHE
            DO 20 IT=1,NTHE
               READ(IH,500)
               READ(IH,505) (PRFHE(IWL,IT,IE),IWL=1,NWLHE)
   20    CONTINUE
C
C        coefficient for the asymptotic profile is determined from
C        the input data
C
         XCLOG=PRFHE(NWLHE,1,1)+2.5*LOG10(WLHE(NWLHE,ILINE))+31.831-
     *         XNEHE(1,ILINE)-2.*LOG10(WL0)
         XKLOG=0.6666667*(XCLOG-0.176)
         XKE=EXP(XKLOG*2.3025851)
         DO 50 ID=1,ND

        !    print*, 'test vturb', ID, VTURB(ID)

            T=TEMP(ID)+2.42E-8*VTURB(ID)
            ANE=ELEC(ID)
            TL=LOG10(T)
            ANEL=LOG10(ANE)
            F00=1.25E-9*ANE**0.666666667
            FXKE=F00*XKE
            DOP=1.E8/WL0*SQRT(4.12E7*T)
            DBETAE=WL0*WL0/(CLIGHT_SI*1e10)/FXKE
            BETADE=DBETAE*DOP
C
C     interpolation to the actual values of temperature and electron
C     density. The result is stored at array PRFHE2, which has indices
C     ILINE  - index of line
C     ID     - depth index
C     IWL    - wavelength index (notice that the wavelength grid may
C              generally be different for different lines
C
            DO 30 IWL=1,NWLHE
               CALL INTHE2(PROF,TL,ANEL,IWL,ILINE)
               PRFHE2(ILINE,ID,IWL)=PROF
   30      CONTINUE
   50   CONTINUE
  100 CONTINUE
C
  500 FORMAT(1X)
  501 FORMAT(//14X,I2,9X,I2/)
  502 FORMAT(2X,I4,1P6E10.3,4(/5X,0P6F10.4)/5X,5F10.4)
  503 FORMAT(2X,I4,F10.3,5F12.3)
  504 FORMAT(2X,I4,F10.2,5F12.2/4X,5F12.2)
  505 FORMAT(10F8.3)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE HE2LIN(ID,T,ANE,ABSOH,EMISH)
      use constants
C
C     opacity and emissivity of He II lines  (these which are not considered
C     explicitly)
C
c     INCLUDE 'PARAMS.FOR'
c     INCLUDE 'MODELP.FOR'
c     INCLUDE 'SYNTHP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID
      real*8, intent(in   ) :: T,ANE
      real*8, intent(inout) :: ABSOH,EMISH

	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
	INCLUDE '../inc/SYNTHP.FOR'

      PARAMETER (SIXTH=1./6., UN=1.)
      DIMENSION PJ(80),SGF(12),FRHE(12),OSCHE2(19),PRF0(99),
     *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
     $ IOPFE1
      COMMON/HE2PAR/IFHE2,IHE2L,ILWHE2,MHE10,MHE20
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
C      DATA FRHE /1.3158143E16,3.2895358E15,1.4624517E15,8.225929E14/
      DATA FRHE /1.3158153D+16, 3.2895381D+15, 1.4624854D+15,
     *           8.2261878D+14, 5.2647201D+14, 3.6560459D+14,
     *           2.6860713D+14, 2.0565220D+14, 1.6249055D+14,
     *           1.3161730D+14, 1.0877460D+14, 9.1400851D+13/
      DATA SGF  /1.578e-18, 3.468e-18, 5.392e-18, 7.322e-18,
     *           9.262e-18, 1.121e-17, 1.316e-17, 1.511e-17,
     *           1.708e-17, 1.904e-17, 2.178e-17, 2.376e-17/
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,
     *        1.604E-2, 1.023E-2, 6.980E-3,
     *        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,
     *        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,
     *        1.038,    1.793E-1, 6.549E-2/
C
      I=ILWHE2
      izz=2
      I0=1
      I1=NFREQ
      DO 1 IJ=I0,I1
         ABSO(IJ)=0.
         EMIS(IJ)=0.
         ABSOH(IJ)=0.
         EMISH(IJ)=0.
    1 CONTINUE
C
C     He III population (either LTE or NLTE, depending on input model)
C
      IF(IELHE2.GT.0) THEN
         ANP=POPUL(NNEXT(IELHE2),ID)
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1
       ELSE
         ANP=RRR(ID,3,2)
         NLHE2=0
      END IF
C
C     populations of the first 60 levels of He II
C
      PP=ANE*4.1412E-16*ANP/T/SQRT(T)
      DO 5 IL=1,60
         X=IL*IL
         IIL=NFIRST(IELHE2)+IL-1
         IF(IL.LE.NLHE2) PJ(IL)=POPUL(IIL,ID)/X
         IF(IL.GT.NLHE2) PJ(IL)=PP*EXP(631479./X/T)
    5 CONTINUE
C
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      ANES=EXP(SIXTH*LOG(ANE))
      F00=3.906e-11*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(4.12E7*T+VTURB(ID))
C
C     determination of the last observable line
C
      MLST=1.82E3/ANE**0.13333333
      IF(MLST.GT.60) MLST=60
C
C -------------------------------------------------------------------
C     overall loop over spectral series (only in the infrared region)
C -------------------------------------------------------------------
C
      ISERU=ILWHE2
      IF(ILWHE2.LE.3) THEN
         ISERL=ILWHE2
       ELSE IF(ILWHE2.LE.5) THEN
         ISERL=ILWHE2-1
       ELSE IF(ILWHE2.LE.7) THEN
         ISERL=ILWHE2-2
       ELSE IF(ILWHE2.LE.9) THEN
         ISERL=ILWHE2-3
       ELSE
         ISERL=ILWHE2-4
      END IF
C
      DO 200 I=ISERL,ISERU
      II=I*I
      PLTEI=PP*EXP(631479./T/II)*II
      POPI=PJ(I)*II
      FLST=FRHE(I)-FRHE(I)*II/MLST**2
C
C     determination of which He II lines contribute in a current
C     frequency region
C
      M1=MHE10
      IF(I.LT.ILWHE2.AND.FRHE(I).GT.FREQ(2)) THEN
         M1=SQRT(FRHE(I)*II/(FRHE(I)-FREQ(2)))
      END IF
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.6..and.M1.LE.6.AND.I.EQ.2) GO TO 10
      IF(grav.lt.6..and.M1.LE.4.AND.I.EQ.1) GO TO 10
      IF(M1.GT.MLST) GO TO 120
      M1=M1-1
      M2=MHE20+3
      IF(M2.GT.60) M2=60
   10 CONTINUE
      if(grav.gt.6.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3
      end if
      IF(M1.LT.I+1) M1=I+1
      IF(M2.GT.60) M2=60
      A=0.
      E=0.
      DO 15 IJ=I0,I1
         ABSO(IJ)=0.
         EMIS(IJ)=0.
   15 CONTINUE
C
C     loop over lines which contribute at given wavelength region
C
      DO 100 J=M1,M2
         JJ=J*J
         IF(I.LE.2) THEN
            WLINE=227.838/(1./II-1./JJ)
          ELSE
            WLINE=227.7776/(1./II-1./JJ)
         END IF
         ILINE=0
         IF(I.EQ.2) THEN
            IF(J.EQ.3.AND.IHE2UV.GT.0) ILINE=1
          ELSE IF(I.EQ.3) THEN
            IF(J.EQ.4.AND.IHE2VI.GT.0) ILINE=8
            IF(J.GT.5.AND.J.LE.10.AND.IHE2UV.GT.0) ILINE=J-3
          ELSE IF(I.EQ.4) THEN
            IF(J.LE.7.AND.IHE2RE.GT.0) ILINE=J+12
            IF(J.GE.8.AND.J.LE.15.AND.IHE2VI.GT.0) ILINE=J+1
         END IF
         IF(ILINE.GT.0) THEN
	    NWL=NWLHE2(ILINE)
            DO 20 IWL=1,NWL
   20          PRF0(IWL)=PRFHE2(ILINE,ID,IWL)
            FID=0.02654*OSCHE2(ILINE)
            DO 50 IJ=I0,I1
               AL=ABS(WLAM(IJ)-WLINE)
               IF(AL.LT.1.E-4) AL=1.E-4
               AL=LOG10(AL)
               DO 30 IWL=1,NWL-1
                  IW0=IWL
                  IF(AL.LE.WLHE2(ILINE,IWL+1)) GO TO 40
   30          CONTINUE
   40          IW1=IW0+1
               PRF=(PRF0(IW0)*(WLHE2(ILINE,IW1)-AL)+PRF0(IW1)*
     *             (AL-WLHE2(ILINE,IW0)))/
     *             (WLHE2(ILINE,IW1)-WLHE2(ILINE,IW0))
               SG=EXP(PRF*2.3025851)*FID
               ABSO(IJ)=ABSO(IJ)+SG
               EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
   50       CONTINUE
          ELSE
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            FXK=F00*XKIJ
c            DOP=DOP0/WLINE
c            DBETA=WLINE*WLINE/2.997925E18/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0/(CLIGHT_SI*1e10)/FXK
            BETAD=DOP*DBETA
            FID=0.02654*FIJ*DBETA
            FID0=0.01497*FIJ0/DOP
            CALL DIVHE2(BETAD,AD,DIV)
            DO 60 IJ=I0,I1
               IF(FREQ(IJ).GT.FLST) GO TO 60
c               BETA=ABS(WLAM(IJ)-WLINE)/FXK
               BETA=ABS(WLAM(IJ)-WL0)/FXK
               SG=STARKA(BETA,BETAD,AD,DIV,UN)*FID
               if(fid0.gt.0.) then
                  xd=beta/betad
                  if(xd.lt.5.) sg=sg+exp(-xd*xd)*fid0
               end if
               ABSO(IJ)=ABSO(IJ)+SG
               EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
   60       CONTINUE
         END IF
  100 CONTINUE
C
C     resulting absorption and emission coefficients
C
      DO 110 IJ=I0,I1
  601    format(1p2e12.4)
         EMIS(IJ)=EMIS(IJ)*II
         ABSO(IJ)=ABSO(IJ)*POPI-EMIS(IJ)
c        write(62,601) abso(ij),emis(ij)
  110 CONTINUE
C
C     pseudocontinuum for wavelengths smaller than the last observable
C     line
C
  120 CONTINUE
      E=PLTEI*SGF(I)
      A=POPI*SGF(I)
      DO 130 IJ=I0,I1
         F=FREQ(IJ)
         IF(F.GT.FRHE(I).OR.F.LT.FLST) GO TO 130
         EMIS(IJ)=E*EXP(-4.79928E-11*F/T)
         ABSO(IJ)=A-EMIS(IJ)
  130 CONTINUE
  150 CONTINUE
C
C     add contributions for different spectral series
C
      DO 160 IJ=I0,I1
         ABSOH(IJ)=ABSOH(IJ)+ABSO(IJ)
         EMISH(IJ)=EMISH(IJ)+EMIS(IJ)
  160 CONTINUE
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
C
C ********************************************************************
C ********************************************************************
C
C
      SUBROUTINE HE2SET
      use SYNTHP_CONT,only:FREQC,NFCONT
      use constants
C
C     Initialization procedure for treating the He II line opacity
C
      implicit real*8(a-h,o-z)

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      COMMON/HE2PAR/IFHE2,IHE2L,ILWHE2,MHE10,MHE20
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2IR
      dimension frhe(12)
      DATA FRHE /1.3158153D+16, 3.2895381D+15, 1.4624854D+15,
     *           8.2261878D+14, 5.2647201D+14, 3.6560459D+14,
     *           2.6860713D+14, 2.0565220D+14, 1.6249055D+14,
     *           1.3161730D+14, 1.0877460D+14, 9.1400851D+13/
C
C     IHE2L=-1  -  He II lines are excluded a priori
C
      IHE2L=-1
      IF(IFHE2.LE.0) RETURN
      IF(FREQC(NFCONT()).GE.1.315812E16) RETURN
      AL0=CLIGHT_SI*1e9/FREQC(1)
      AL1=CLIGHT_SI*1e9/FREQC(NFCONT())
c      IF(AL0.GT.390.) RETURN
      if(grav.lt.6.) then
         IF(AL0.GT.31..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.26.1.AND.AL1.LT.29.8) RETURN
         IF(AL0.GT.24.8.AND.AL1.LT.25.1) RETURN
         IF(AL0.GT.122.1.AND.AL1.LT.162.9) RETURN
         IF(AL0.GT.165.1.AND.AL1.LT.204.9) RETURN
         IF(AL0.GT.109..AND.AL1.LT.120.9) RETURN
         IF(AL0.GT.103..AND.AL1.LT.107.9) RETURN
         IF(AL0.GT.99.7.AND.AL1.LT.102.) RETURN
         IF(AL0.GT.320.8.AND.AL1.LT.364.4) RETURN
         IF(AL0.GT.273.8.AND.AL1.LT.319.8) RETURN
         IF(AL0.GT.251.6.AND.AL1.LT.272.8) RETURN
         IF(AL0.GT.239.0.AND.AL1.LT.250.6) RETURN
         IF(AL0.GT.231.1.AND.AL1.LT.238.0) RETURN
         IF(AL0.GT.225.8.AND.AL1.LT.230.1) RETURN
       else
         IF(AL0.GT.33..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.124.1.AND.AL1.LT.160.9) RETURN
         IF(AL0.GT.167.1.AND.AL1.LT.202.9) RETURN
         IF(AL0.GT.111..AND.AL1.LT.118.9) RETURN
         IF(AL0.GT.322.8.AND.AL1.LT.364.4) RETURN
         IF(AL0.GT.275.8.AND.AL1.LT.317.8) RETURN
         IF(AL0.GT.253.6.AND.AL1.LT.270.8) RETURN
         IF(AL0.GT.241.0.AND.AL1.LT.248.6) RETURN
         IF(AL0.GT.233.1.AND.AL1.LT.236.0) RETURN
      end if
C
C     otherwise, He II lines are included
C
      IHE2L=1
      MHE10=60
      MHE20=60
      IF(AL1.LT.91.) THEN
         ILWHE2=1
       ELSE IF(AL0.LT.204.) THEN
         ILWHE2=2
       ELSE IF(AL0.LT.364.) THEN
         ILWHE2=3
       ELSE IF(AL0.LT.569.) THEN
         ILWHE2=4
       ELSE IF(AL0.LT.819.) THEN
         ILWHE2=5
       ELSE IF(AL0.LT.1116.) THEN
         ILWHE2=6
       ELSE IF(AL0.LT.1457.) THEN
         ILWHE2=7
       ELSE IF(AL0.LT.1844.) THEN
         ILWHE2=8
       ELSE IF(AL0.LT.2277.) THEN
         ILWHE2=9
       ELSE IF(AL0.LT.2756.) THEN
         ILWHE2=10
       ELSE IF(AL0.LT.3279.) THEN
         ILWHE2=11
       ELSE
         ILWHE2=12
      END IF
      FRION=FRHE(ILWHE2)
      FR1=FRION*ILWHE2*ILWHE2
      IF(FRION.GT.FREQ(2)) MHE10=SQRT(FR1/(FRION-FREQC(NFCONT())))
      IF(FRION.GT.FREQ(1)) MHE20=SQRT(FR1/(FRION-FREQC(1)))
      WRITE(6,601) ILWHE2,MHE20+1
  601 FORMAT(1H0/ ' *** HE II LINES CONTRIBUTE'/
     * '     THE NEAREST LINE ON THE SHORT-WAVELENGTH SIDE IS',
     * I3,'  TO ',I3/)
      RETURN
      END subroutine
C
C
C     ****************************************************************
C
C

      FUNCTION HEPHOT(S,L,N,FREQ)

        
C     ===========================
C
C           EVALUATES HE I PHOTOIONIZATION CROSS SECTION USING SEATON
C           FERNLEY'S CUBIC FITS TO THE OPACITY PROJECT CROSS SECTIONS
C           UP TO SOME ENERGY "EFITM" IN THE RESONANCE-FREE ZONE.  BEYOND
C           THIS ENERGY LINEAR FITS TO LOG SIGMA IN LOG (E/E0) ARE USED.
C           THIS EXTRAPOLATION SHOULD BE USED UP TO THE BEGINNING OF THE
C           RESONANCE ZONE "XMAX", BUT AT PRESENT IT IS USED THROUGH IT.
C           BY CHANGING A FEW LINES THAT ARE PRESENTLY COMMENTED OUT,
C           FOR ENERGIES IN THE RESONANCE ZONE A VALUE OF 1/100 OF THE
C           THRESHOLD CROSS SECTION IS USED -- THIS IS PURELY AD HOC AND
C           ONLY A TEMPORARY MEASURE.  OBVIOUSLY ANY OTHER VALUE OR FUNCTIONAL
C           FORM CAN BE INSERTED HERE.
C
C           CALLING SEQUENCE INCLUDES:
C                S = MULTIPLICITY, EITHER 1 OR 3
C                L = ANGULAR MOMENTUM, 0, 1, OR 2;
C                    for L > 2 - hydrogenic expresion
C                FREQ = FREQUENCY
C
C           DGH JUNE 1988 JILA, slightly modified by I.H.
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: S,L,N
      real*8, intent(in   ) :: FREQ
!       integer,intent(inout) :: 
!       real*8, intent(inout) :: 
!       integer,intent(  out) :: 
!       real*8, intent(  out) :: 
!       integer :: 
!       real*8  ::      
      INTEGER SS,LL
      DIMENSION COEF(4,53),IST(3,2),N0(3,2),
     *          FL0(53),A(53),B(53),XFITM(53)
      !DIMENSION XMAX(53)
C
      DATA IST/1,36,20,11,45,28/
      DATA N0/1,2,3,2,2,3/
C
      DATA FL0/
     . 2.521D-01,-5.381D-01,-9.139D-01,-1.175D+00,-1.375D+00,-1.537D+00,
     .-1.674D+00,-1.792D+00,-1.896D+00,-1.989D+00,-4.555D-01,-8.622D-01,
     .-1.137D+00,-1.345D+00,-1.512D+00,-1.653D+00,-1.774D+00,-1.880D+00,
     .-1.974D+00,-9.538D-01,-1.204D+00,-1.398D+00,-1.556D+00,-1.690D+00,
     .-1.806D+00,-1.909D+00,-2.000D+00,-9.537D-01,-1.204D+00,-1.398D+00,
     .-1.556D+00,-1.690D+00,-1.806D+00,-1.909D+00,-2.000D+00,-6.065D-01,
     .-9.578D-01,-1.207D+00,-1.400D+00,-1.558D+00,-1.692D+00,-1.808D+00,
     .-1.910D+00,-2.002D+00,-5.749D-01,-9.352D-01,-1.190D+00,-1.386D+00,
     .-1.547D+00,-1.682D+00,-1.799D+00,-1.902D+00,-1.995D+00/
C
      DATA XFITM/
     . 3.262D-01, 6.135D-01, 9.233D-01, 8.438D-01, 1.020D+00, 1.169D+00,
     . 1.298D+00, 1.411D+00, 1.512D+00, 1.602D+00, 7.228D-01, 1.076D+00,
     . 1.206D+00, 1.404D+00, 1.481D+00, 1.464D+00, 1.581D+00, 1.685D+00,
     . 1.777D+00, 9.586D-01, 1.187D+00, 1.371D+00, 1.524D+00, 1.740D+00,
     . 1.854D+00, 1.955D+00, 2.046D+00, 9.585D-01, 1.041D+00, 1.371D+00,
     . 1.608D+00, 1.739D+00, 1.768D+00, 1.869D+00, 1.803D+00, 7.360D-01,
     . 1.041D+00, 1.272D+00, 1.457D+00, 1.611D+00, 1.741D+00, 1.855D+00,
     . 1.870D+00, 1.804D+00, 9.302D-01, 1.144D+00, 1.028D+00, 1.210D+00,
     . 1.362D+00, 1.646D+00, 1.761D+00, 1.863D+00, 1.954D+00/
C
      DATA A/
     . 6.95319D-01, 1.13101D+00, 1.36313D+00, 1.51684D+00, 1.64767D+00,
     . 1.75643D+00, 1.84458D+00, 1.87243D+00, 1.85628D+00, 1.90889D+00,
     . 9.01802D-01, 1.25389D+00, 1.39033D+00, 1.55226D+00, 1.60658D+00,
     . 1.65930D+00, 1.68855D+00, 1.62477D+00, 1.66726D+00, 1.83599D+00,
     . 2.50403D+00, 3.08564D+00, 3.56545D+00, 4.25922D+00, 4.61346D+00,
     . 4.91417D+00, 5.19211D+00, 1.74181D+00, 2.25756D+00, 2.95625D+00,
     . 3.65899D+00, 4.04397D+00, 4.13410D+00, 4.43538D+00, 4.19583D+00,
     . 1.79027D+00, 2.23543D+00, 2.63942D+00, 3.02461D+00, 3.35018D+00,
     . 3.62067D+00, 3.85218D+00, 3.76689D+00, 3.49318D+00, 1.16294D+00,
     . 1.86467D+00, 2.02110D+00, 2.24231D+00, 2.44240D+00, 2.76594D+00,
     . 2.93230D+00, 3.08109D+00, 3.21069D+00/
C
      DATA B/
     .-1.29000D+00,-2.15771D+00,-2.13263D+00,-2.10272D+00,-2.10861D+00,
     .-2.11507D+00,-2.11710D+00,-2.08531D+00,-2.03296D+00,-2.03441D+00,
     .-1.85905D+00,-2.04057D+00,-2.02189D+00,-2.05930D+00,-2.03403D+00,
     .-2.02071D+00,-1.99956D+00,-1.92851D+00,-1.92905D+00,-4.58608D+00,
     .-4.40022D+00,-4.39154D+00,-4.39676D+00,-4.57631D+00,-4.57120D+00,
     .-4.56188D+00,-4.55915D+00,-4.41218D+00,-4.12940D+00,-4.24401D+00,
     .-4.40783D+00,-4.39930D+00,-4.25981D+00,-4.26804D+00,-4.00419D+00,
     .-4.47251D+00,-3.87960D+00,-3.71668D+00,-3.68461D+00,-3.67173D+00,
     .-3.65991D+00,-3.64968D+00,-3.48666D+00,-3.23985D+00,-2.95758D+00,
     .-3.07110D+00,-2.87157D+00,-2.83137D+00,-2.82132D+00,-2.91084D+00,
     .-2.91159D+00,-2.91336D+00,-2.91296D+00/
C
!       DATA XMAX/
!      . 3.84251D-01, 9.89081D-01, 1.34758D+00, 1.60093D+00, 1.79766D+00,
!      . 1.95867D+00, 2.09411D+00, 2.21217D+00, 2.31555D+00, 2.40796D+00,
!      . 8.54216D-01, 1.25914D+00, 1.52406D+00, 1.74001D+00, 1.91064D+00,
!      . 2.05087D+00, 2.17074D+00, 2.27614D+00, 2.37017D+00, 1.38737D+00,
!      . 1.63108D+00, 1.82209D+00, 1.97849D+00, 2.11121D+00, 2.26666D+00,
!      . 2.38139D+00, 2.47258D+00, 1.36716D+00, 1.60889D+00, 1.86129D+00,
!      . 2.01803D+00, 2.15086D+00, 2.26617D+00, 2.36801D+00, 2.45919D+00,
!      . 9.86915D-01, 1.35996D+00, 1.60525D+00, 1.81878D+00, 1.97653D+00,
!      . 2.10895D+00, 2.22480D+00, 2.32641D+00, 2.47969D+00, 1.06043D+00,
!      . 1.39913D+00, 1.64584D+00, 1.85192D+00, 2.01031D+00, 2.15231D+00,
!      . 2.26861D+00, 2.37124D+00, 2.46305D+00/
C
      DATA ((COEF(I,J),I=1,4),J=1,10)/
     . 8.734D-01,-1.545D+00,-1.093D+00, 5.918D-01, 9.771D-01,-1.567D+00,
     .-4.739D-01,-1.302D-01, 1.174D+00,-1.638D+00,-2.831D-01,-3.281D-02,
     . 1.324D+00,-1.692D+00,-2.916D-01, 9.027D-02, 1.445D+00,-1.761D+00,
     .-1.902D-01, 4.401D-02, 1.546D+00,-1.817D+00,-1.278D-01, 2.293D-02,
     . 1.635D+00,-1.864D+00,-8.252D-02, 9.854D-03, 1.712D+00,-1.903D+00,
     .-5.206D-02, 2.892D-03, 1.782D+00,-1.936D+00,-2.952D-02,-1.405D-03,
     . 1.845D+00,-1.964D+00,-1.152D-02,-4.487D-03/
      DATA ((COEF(I,J),I=1,4),J=11,19)/
     . 7.377D-01,-9.327D-01,-1.466D+00, 6.891D-01, 9.031D-01,-1.157D+00,
     .-7.151D-01, 1.832D-01, 1.031D+00,-1.313D+00,-4.517D-01, 9.207D-02,
     . 1.135D+00,-1.441D+00,-2.724D-01, 3.105D-02, 1.225D+00,-1.536D+00,
     .-1.725D-01, 7.191D-03, 1.302D+00,-1.602D+00,-1.300D-01, 7.345D-03,
     . 1.372D+00,-1.664D+00,-8.204D-02,-1.643D-03, 1.434D+00,-1.715D+00,
     .-4.646D-02,-7.456D-03, 1.491D+00,-1.760D+00,-1.838D-02,-1.152D-02/
      DATA ((COEF(I,J),I=1,4),J=20,27)/
     . 1.258D+00,-3.442D+00,-4.731D-01,-9.522D-02, 1.553D+00,-2.781D+00,
     .-6.841D-01,-4.083D-03, 1.727D+00,-2.494D+00,-5.785D-01,-6.015D-02,
     . 1.853D+00,-2.347D+00,-4.611D-01,-9.615D-02, 1.955D+00,-2.273D+00,
     .-3.457D-01,-1.245D-01, 2.041D+00,-2.226D+00,-2.669D-01,-1.344D-01,
     . 2.115D+00,-2.200D+00,-1.999D-01,-1.410D-01, 2.182D+00,-2.188D+00,
     .-1.405D-01,-1.460D-01/
      DATA ((COEF(I,J),I=1,4),J=28,35)/
     . 1.267D+00,-3.417D+00,-5.038D-01,-1.797D-02, 1.565D+00,-2.781D+00,
     .-6.497D-01,-5.979D-03, 1.741D+00,-2.479D+00,-6.099D-01,-2.227D-02,
     . 1.870D+00,-2.336D+00,-4.899D-01,-6.616D-02, 1.973D+00,-2.253D+00,
     .-3.972D-01,-8.729D-02, 2.061D+00,-2.212D+00,-3.072D-01,-1.060D-01,
     . 2.137D+00,-2.189D+00,-2.352D-01,-1.171D-01, 2.205D+00,-2.186D+00,
     .-1.621D-01,-1.296D-01/
      DATA ((COEF(I,J),I=1,4),J=36,44)/
     . 1.129D+00,-3.149D+00,-1.910D-01,-5.244D-01, 1.431D+00,-2.511D+00,
     .-3.710D-01,-1.933D-01, 1.620D+00,-2.303D+00,-3.045D-01,-1.391D-01,
     . 1.763D+00,-2.235D+00,-1.829D-01,-1.491D-01, 1.879D+00,-2.215D+00,
     .-9.003D-02,-1.537D-01, 1.978D+00,-2.213D+00,-2.066D-02,-1.541D-01,
     . 2.064D+00,-2.220D+00, 3.258D-02,-1.527D-01, 2.140D+00,-2.225D+00,
     . 6.311D-02,-1.455D-01, 2.208D+00,-2.229D+00, 7.977D-02,-1.357D-01/
      DATA ((COEF(I,J),I=1,4),J=45,53)/
     . 1.204D+00,-2.809D+00,-3.094D-01, 1.100D-01, 1.455D+00,-2.254D+00,
     .-4.795D-01, 6.872D-02, 1.619D+00,-2.109D+00,-3.357D-01,-2.532D-02,
     . 1.747D+00,-2.065D+00,-2.317D-01,-5.224D-02, 1.853D+00,-2.058D+00,
     .-1.517D-01,-6.647D-02, 1.943D+00,-2.055D+00,-1.158D-01,-6.081D-02,
     . 2.023D+00,-2.070D+00,-6.470D-02,-6.800D-02, 2.095D+00,-2.088D+00,
     .-2.357D-02,-7.250D-02, 2.160D+00,-2.107D+00, 1.065D-02,-7.542D-02/
C
      IF(L.GT.2) GO TO 20
C
C          SELECT BEGINNING AND END OF COEFFICIENTS
C
      SS=(S+1)/2
      LL=L+1
      NSL0=N0(LL,SS)
      I=IST(LL,SS)+N-NSL0
C
C          EVALUATE CROSS SECTION
C
      FL=LOG10(FREQ/3.28805E15)
      X=FL-FL0(I)
      IF(X.GE.-0.001D0) THEN
         IF(X.LT.XFITM(I)) THEN
            P=COEF(4,I)
            DO 10 K=1,3
               P=X*P+COEF(4-K,I)
   10       CONTINUE
            HEPHOT=1.D-18*1.D1**P
          ELSE
C           OTHERWISE REMOVE INSTRUCTION AND 3 FOLLOWING "C"
C         ELSE IF(X.LT.XMAX(I)) THEN
            HEPHOT=1.D-18*1.D1**(A(I)+B(I)*X)
C         ELSE
C           HEPHOT=1.D-18*1.D1**(COEF(1,I)-2.0D0)
          END IF
       ELSE
          HEPHOT=0.
      END IF
      RETURN
C
C     Hydrogenic expression for L > 2
C      [multiplied by relative population of state (s,l,n), ie.
C       by  stat.weight(s,l)/stat.weight(n)]
C
   20 GN=2.D0*N*N
      HEPHOT=2.815D29/FREQ/FREQ/FREQ/N**5*(2*L+1)*S/GN
      RETURN
      END FUNCTION
C
C
C     ******************************************************************
C
C

      SUBROUTINE HESET(IL,ALM,EXCL,EXCU,ION,IPRF0,ILWN,IUPN)
C     ======================================================
C
C     Auxiliary procedure for INISET - set up quantities:
C     IPRF0      - index for the procedure evaluating standard absorption
C                  profile coefficient for He I lines - see GAMHE
C     ILWN,IUPN  - only in NLTE option is switched on;
C                  indices of the lower and upper level associated with
C                  the given line
C
C     Input: IL - line index
C            ALM - line wavelength in nm
C            EXCL - excitation potential of the lower level (in cm**-1)
C            EXCU - excitation potential of the upper level (in cm**-1)
C            ION  - ionisation degree (1=neutrals, 2=once ionized, etc.)
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: IL,ION
      real*8, intent(in   ) :: ALM,EXCL,EXCU
      integer,intent(inout) :: IPRF0,ILWN,IUPN
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      PARAMETER (HC=5.0344731D15)
      DIMENSION JU(23),NU(23),IT(23)
      DATA IT/1,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0/
      DATA NU/6,6,9,3,8,4,7,5,6,6,5,4,4,4,3,4,3,3,5,5,7,8,10/
      DATA JU/15,3,5,9,5,3,5,3,5,1,1,15,3,5,3,1,15,5,15,5,1,1,1/
C
C ******* He I  ***********
C
      IF(ION.NE.1) GO TO 20
C
C     switch IPRF0 - see GAMHE
C
      ALAM=ALM*10.
      IPRF=0
      IF(ABS(ALAM-3819.60).LT.1.) IPRF=1
      IF(ABS(ALAM-3867.50).LT.1.) IPRF=2
      IF(ABS(ALAM-3871.79).LT.1.) IPRF=3
      IF(ABS(ALAM-3888.65).LT.1.) IPRF=4
      IF(ABS(ALAM-3926.53).LT.1.) IPRF=5
      IF(ABS(ALAM-3964.73).LT.1.) IPRF=6
      IF(ABS(ALAM-4009.27).LT.1.) IPRF=7
      IF(ABS(ALAM-4120.80).LT.1.) IPRF=8
      IF(ABS(ALAM-4143.76).LT.1.) IPRF=9
      IF(ABS(ALAM-4168.97).LT.1.) IPRF=10
      IF(ABS(ALAM-4437.55).LT.1.) IPRF=11
      IF(ABS(ALAM-4471.50).LT.1.) IPRF=12
      IF(ABS(ALAM-4713.20).LT.1.) IPRF=13
      IF(ABS(ALAM-4921.93).LT.1.) IPRF=14
      IF(ABS(ALAM-5015.68).LT.1.) IPRF=15
      IF(ABS(ALAM-5047.74).LT.1.) IPRF=16
      IF(ABS(ALAM-5875.70).LT.1.) IPRF=17
      IF(ABS(ALAM-6678.15).LT.1.) IPRF=18
      IF(ABS(ALAM-4026.20).LT.1.) IPRF=19
      IF(ABS(ALAM-4387.93).LT.1.) IPRF=20
      IF(ABS(ALAM-4023.97).LT.1.) IPRF=21
      IF(ABS(ALAM-3935.91).LT.1.) IPRF=22
      IF(ABS(ALAM-3833.55).LT.1.) IPRF=23
      IF(IPRF.GT.0.AND.IPRF.LE.20) IPRF0=IPRF
C
C     Indices of NLTE levels associated with the given line
C
      IF(INLTE.gt.5.OR.IELHE1.EQ.0) RETURN
      N0I=NFIRST(IELHE1)
      N1I=NLAST(IELHE1)
      EION=ENION(N0I)*HC
      ILW=0
      IUN=0
      NQL=0
      IF(IPRF.GT.0) NQL=NU(IPRF)
c-pr
      print *,' routine HESET - HeI lines'
      print '(A,2F18.2,I7)',' EXCL, EXCU, NQL: ',EXCL,EXCU,NQL
      DO 10 I=N0I,N1I
	 NQ=NQUANT(I)
	 EX=EION-ENION(I)*HC
	 IF(ABS(EXCL-EX).LT.100.) THEN
	    ILW=I
	    IGL=INT(G(I)+0.001)
	 END IF
	 IF(NQ.EQ.NQL) THEN
	    IG=INT(G(I)+0.001)
	    IF(IT(IPRF).EQ.0) THEN
	       IF(NQ.EQ.2.AND.IG.EQ.JU(IPRF)) IUN=I
	       IF(NQ.EQ.3) THEN
	          IF(IG.EQ.JU(IPRF)) THEN
	             IF(IG.EQ.1.OR.IG.EQ.5) IUN=I
	             IF(IG.EQ.3.AND.IGL.EQ.1) IUN=I
		   ELSE
		     IF(IG.EQ.9) IUN=I
                  END IF
	       END IF
	       IF(NQ.EQ.4) THEN
	          IF(IG.EQ.JU(IPRF)) THEN
	             IF(IG.EQ.1.OR.IG.EQ.5.OR.IG.EQ.7) IUN=I
	             IF(IG.EQ.3.AND.IGL.EQ.1) IUN=I
		   ELSE
		     IF(IG.EQ.16) IUN=I
                  END IF
	       END IF
	       IF(IG.EQ.25.OR.IG.EQ.36) IUN=I
	       IF(IG.EQ.49.OR.IG.EQ.64.OR.IG.EQ.81) IUN=I
	       IF(IG.EQ.100.OR.IG.EQ.121.OR.IG.EQ.144) IUN=I
             ELSE
	       IF(NQ.EQ.3) THEN
	          IF(IG.EQ.JU(IPRF)) THEN
	             IF(IG.EQ.9.OR.IG.EQ.15) IUN=I
	             IF(IG.EQ.3.AND.IGL.EQ.9) IUN=I
		   ELSE
		     IF(IG.EQ.27) IUN=I
                  END IF
	       END IF
	       IF(NQ.EQ.4) THEN
	          IF(IG.EQ.JU(IPRF)) THEN
	             IF(IG.EQ.9.OR.IG.EQ.15.OR.IG.EQ.21) IUN=I
	             IF(IG.EQ.3.AND.IGL.EQ.9) IUN=I
		   ELSE
		     IF(IG.EQ.48) IUN=I
                  END IF
	       END IF
	       IF(IG.EQ.75) IUN=I
	       IF(IG.EQ.108.OR.IG.EQ.147.OR.IG.EQ.192) IUN=I
	       IF(IG.EQ.243.OR.IG.EQ.300.OR.IG.EQ.363) IUN=I
	    END IF
	    IF(NQ.EQ.2.AND.IG.EQ.16) IUN=I
	    IF(NQ.EQ.3.AND.IG.EQ.36) IUN=I
	    IF(NQ.EQ.4.AND.IG.EQ.64) IUN=I
	    IF(NQ.EQ.5.AND.IG.EQ.100) IUN=I
	    IF(NQ.EQ.6.AND.IG.EQ.144) IUN=I
	    IF(NQ.EQ.7.AND.IG.EQ.196) IUN=I
	    IF(NQ.EQ.8.AND.IG.EQ.256) IUN=I
	    IF(NQ.EQ.9.AND.IG.EQ.324) IUN=I
	    IF(NQ.EQ.10.AND.IG.EQ.400) IUN=I
	 END IF
   10 CONTINUE
      print *, 'il,iprof,ilw,iupn',il,iprf,ilw,iun

      ILWN=ILW
      IUPN=IUN
C
C ******* He II ***********
C
   20 IF(ION.NE.2.OR.IELHE2.LE.0) RETURN
      N0I=NFIRST(IELHE2)
      NLHE2=NLAST(IELHE2)-N0I+1
      XL=SQRT(1./(1.-EXCL/438916.146))
      ILW=INT(XL)
      IF((FLOAT(ILW)-XL).LT.0.) ILW=ILW+1
      XU=SQRT(1./(1.-EXCU/438916.146))
      IUN=INT(XU)
      IF((FLOAT(IUN)-XU).LT.0.) IUN=IUN+1
      IF(ILW.LE.NLHE2) ILWN=ILW+N0I-1
      IF(IUN.LE.NLHE2) IUPN=IUN+N0I-1
c-pr
      print '(A,2i7)',
     $     ' routine HESET - HeII lines:        ilw,iupn=',ilwn,iupn
      RETURN
      END SUBROUTINE
C
C
C     ****************************************************************
C
C
      FUNCTION HIDALG(IB,FR)
      use constants
C     ======================
C
C     Read table of wavelengths and photo-ionization cross-sections
C     from Hidalgo (1968, Ap. J., 153, 981) for the species indicated by IB
C     (Hidalgo's number = INDEX = -IB-100).
C     Compute linearly interpolated value of the cross-section
C     at the frequency FR.
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IB
      real*8, intent(in   ) :: FR
      real*8 :: HIDALG
      DIMENSION WL1(20),WL2(20),WL(20),SIG0(20,24),SIGS(20)
C
      DATA WL1 /
     *  39.1, 80.9, 97.6,100.1,104.3,107.2,108.7,111.9,113.6,115.4,
     * 117.1,119.0,124.8,126.9,129.1,131.3,133.6,136.0,138.5,141.1/
      DATA WL2 /
     *  68.5, 80.9,100.1,120.9,158.8,165.7,177.3,190.6,200.7,206.2,
     * 211.9,218.0,224.5,231.3,246.3,5*0./
      DATA SIG0 /
     *120*0.,
     *.0460,.2400,.3500,.3700,.4000,.4300,.4400,.4600,.4700,.4900,
     *.5000,.5200,.5700,.6200, 6*0.,
     * 80*0.,
     *.0092,.1000,.1900,.2100,.2300,.2500,.2600,.2900,.3000,.3200,
     *.3400,.3500,.4100,.4300,.4500,.4800,.5000,.5300,.5600,.5900,
     * 20*0.,
     *.3400,.4600,.6300,.7700,.9100,1.080, 14*0.,
     * 20*0.,
     *.0064,.1100,.2200,.4100,.9400,1.000,1.300,1.600, 12*0.,
     * 80*0.,
     *.0370,.0650,.1300,.2400,.5500,.6300,.7700,.9500,1.100,1.250,
     * 10*0.,
     * 40*0.,
     *.0220,.0390,.0800,.1500,.3500,.4000,.4900,.6200,.7200,.7800,
     *.8500,.9300,1.020,
     * 7*0./
C
      INDEX=-IB-100
      NUM=20
      IF(INDEX.GE.13.AND.INDEX.LE.27) NUM=15
      DO 10 I=1,NUM
         IF(INDEX.LT.13) WL(I)=WL1(I)
         IF(INDEX.GE.13) WL(I)=WL2(I)
         SIGS(I)=SIG0(I,INDEX)
   10 CONTINUE
C
      WLAM=CLIGHT_SI*1e10/FR
      IL=1
      IR=NUM
      DO 50 I=1,NUM-1
        IF(WLAM.GE.WL(I).AND.WLAM.LE.WL(I+1)) THEN
          IL=I
          IR=I+1
          GO TO 60
        ENDIF
 50   CONTINUE
C
C     LINEAR INTERPOLATION:
C
 60   SIGM=(SIGS(IR)-SIGS(IL))*(WLAM-WL(IL))/(WL(IR)-WL(IL))
     *      + SIGS(IL)
C
C     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE:
C
       IF(WLAM.LE.WL(1)) SIGM=SIGS(1)
       IF(WLAM.GE.WL(NUM)) SIGM=SIGS(NUM)
C
C     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION:
C
c       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL)
C
      HIDALG=SIGM*1.E-18
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
!       SUBROUTINE HYDLIN(ID,T,ANE,ABSOH,EMISH)
!       use constants
!       use MOD_HYDTAB
! C
! C     opacity and emissivity of hydrogen lines
! C
!       implicit real*8(a-h,o-z)
!       integer,intent(in   ) :: ID
!       real*8, intent(in   ) :: T,ANE
!       real*8, intent(inout) :: ABSOH,EMISH
!       INCLUDE '../inc/PARAMS.FOR'
!       INCLUDE '../inc/MODELP.FOR'
!       INCLUDE '../inc/SYNTHP.FOR'
!       PARAMETER (FRH1=3.28805E15,FRH2=FRH1/4.,UN=1.,SIXTH=1./6.)
!       PARAMETER (CPP=4.1412E-16,CCOR=0.09,CPJ=157803.,CLST=1.1E3)
!       PARAMETER (C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.)
!       PARAMETER (CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/(CLIGHT_SI*1e10))
!       PARAMETER (CID1=0.01497)
!       logical lwph,lquasi
!       DIMENSION PJ(40),PRF0(36),WLINE(8),OSCB(8),
!      *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
!       COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2IR
!       COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
!      $ IOPFE1
!       COMMON/HYLPAR/IHYL,ILOWH,M10,M20
!       COMMON/AUXVCS/XK,FXK,BETAD,DBETA
!       COMMON/BALVCS/PRFBAL(8,MDEPTH,36),WLBAL(8,36),NWLH(8)
!       common/wprob/wph(mdepth,40),acor(mdepth),lwph
!       DATA WLINE  /6562.80, 4861.32, 4340.46, 4101.73,
!      *             3970.07, 3889.05, 3835.38, 3797.90/
!       DATA OSCB   /0.6407,  0.1193,   0.04467,  0.02209,
!      *             1.27D-2, 8.036D-3, 5.429D-3, 3.851D-3/
!       DATA FRH    /3.289017E15/
!       PARAMETER (HINV=1.5092973D26)
! C
!       if(iath.le.0) return
!       izz=1
!       i0=1
!       i1=nfreq
!       ABSO (I0:I1)=0.
!       EMIS (I0:I1)=0.
!       ABSOH(I0:I1)=0.
!       EMISH(I0:I1)=0.
!       T1=UN/T
!       SQT=SQRT(T)
!       ANES=EXP(SIXTH*LOG(ANE))
! C
! C     populations of the first 40 levels of hydrogen
! C
!       ANP=POPUL(NKH,ID)
!       PP=CPP*ANE*ANP*T1/SQT
!       NLH=N1H-N0HN+1
!       if(ifwop(n1h).lt.0) nlh=nlh-1
!       DO IL=1,40
!          X=IL*IL
!          IF(IL.LE.NLH) PJ(IL)=POPUL(N0HN+IL-1,ID)/X
!          IF(IL.GT.NLH) PJ(IL)=PP*EXP(CPJ/X*T1)
!       ENDDO
!       p2=pp*exp(cpj4*t1)
! C
! C     Frequency- and line-independent parameters for evaluating the
! C     asymptotic Stark profile
! C
!       F00=C00*ANES*ANES*ANES*ANES
!       DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID))
! C
! C     determination of the last observable line
! C
!       MLST=CLST*EXP(-LOG(ANE)/7.5)
!       IF(MLST.GT.40) MLST=40
! C
! C -------------------------------------------------------------------
! C     overall loop over spectral series (only in the infrared region)
! C -------------------------------------------------------------------
! C
!       ISERL=ILOWH
!       ISERU=ILOWH
!       IF(WLAM(I0).GT.17000..AND.WLAM(I1).LT.21000.) THEN
!          ISERL=3
!          ISERU=4
!        ELSE IF(WLAM(I0).GT.22700.) THEN
!          ISERL=4
!          ISERU=5
!          IF(WLAM(I0).GT.32800.) ISERU=6
!          IF(WLAM(I0).GT.44660.) ISERU=7
!       END IF
! C
!       DO 200 I=ISERL,ISERU
!       II=I*I
!       XII=UN/II
!       PLTEI=PP*EXP(CPJ*T1*XII)*II
!       POPI=PJ(I)*II
!       IF(I.EQ.1) FRH=3.28805E15
!       FEDGE=FRH*XII
!       FLST=FEDGE-FRH/MLST**2
!       IF(I.LE.NLH) FEDGE=ENION(I+N0HN-1)*HINV
! C
! C     determination of which hydrogen lines contribute in a current
! C     frequency region
! C
!       M1=M10
!       IF(I.LT.ILOWH) M1=ILOWH-1
!       M2=M1+1
!       IF(M1.LT.I+1) M1=I+1
!       IF(grav.lt.3..and.M1.LE.6.AND.I.EQ.2) GO TO 10
!       IF(grav.lt.3..and.M1.LE.4.AND.I.EQ.1) GO TO 10
!       IF(M1.GT.MLST .and. .not.lwph) GO TO 120
!       M1=M1-1
!       M2=M20+3
!       IF(M1.LT.I+1) M1=I+1
!    10 CONTINUE
!       if(grav.gt.3.) then
!          m2=m2+5
!          m1=m1-3
!          if(m1.gt.i+6) m1=m1-3
!       end if
!       if(grav.gt.6.) then
!          m2=m2+2
!          m1=m1-1
!          if(m1.gt.i+6) m1=m1-1
!       end if
!       IF(M1.LT.I+1) M1=I+1
!       IF(M2.GT.40) M2=40
!       if(id.eq.1) write(6,666) i,m1,m2
!   666 format(/' hydrogen lines contribute - ilow=',i2,', iup from ',i3,
!      *       ' to',i3/)
! C
!       A=0.
!       E=0.
!       ABSO(I0:I1)=0.
!       EMIS(I0:I1)=0.
! C
! C     loop over lines which contribute at given wavelength region
! C
!       DO 100 J=M1,M2
! c        IF(I.EQ.1.AND.IOPHLI.GT.0.AND.J.LE.5) GO TO 100
!          IF(I.EQ.2.AND.J.LE.10.AND.IBVCS.GT.0) THEN
! 	    ILINE=J-2
! 	    NWL=NWLH(ILINE)
!             DO 20 IWL=1,NWL
!    20          PRF0(IWL)=PRFBAL(ILINE,ID,IWL)
!             FID=CID*OSCB(ILINE)
!             DO 50 IJ=I0,I1
!                AL=ABS(WLAM(IJ)-WLINE(ILINE))
!                IF(AL.LT.1.E-4) AL=1.E-4
!                AL=LOG10(AL)
!                DO 30 IWL=1,NWL-1
!                   IW0=IWL
!                   IF(AL.LE.WLBAL(ILINE,IWL+1)) GO TO 40
!    30          CONTINUE
!    40          IW1=IW0+1
!                PRF=(PRF0(IW0)*(WLBAL(ILINE,IW1)-AL)+PRF0(IW1)*
!      *             (AL-WLBAL(ILINE,IW0)))/
!      *             (WLBAL(ILINE,IW1)-WLBAL(ILINE,IW0))
!                SG=EXP(PRF*AL10)*FID*wph(id,j)
!                ABSO(IJ)=ABSO(IJ)+SG
!                EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
!    50       CONTINUE
!           ELSE
!             CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
!             FXK=F00*XKIJ
!             FXK1=UN/FXK
!             DOP=DOP0/WL0
!             DBETA=WL0*WL0*CINV*FXK1
!             BETAD=DOP*DBETA
!             FID=CID*FIJ*DBETA*wph(id,j)
!             FID0=CID1*FIJ0/DOP*wph(id,j)
!             CALL DIVSTR(BETAD,AD,DIV)
!             lquasi=i.eq.1.and.j.eq.2.and.iophli.ge.3
!             DO 60 IJ=I0,I1
!                fr=freq(ij)
!                IF(FREQ(IJ).GT.FLST .and. .not.lwph) GO TO 60
!                BETA=ABS(WLAM(IJ)-WL0)*FXK1
! c              if(beta.gt.5.e4) go to 60
!                SG=STARKA(BETA,BETAD,AD,DIV,TWO)*FID
!                if(iophli.eq.2.and.i.eq.1.and.j.eq.2)
!      *            sg=sg*feautr(fr,id)
!                if(fid0.gt.0.) then
!                   xd=beta/betad
!                   if(xd.lt.5.) sg=sg+exp(-xd*xd)*fid0
!                end if
!                ABSO(IJ)=ABSO(IJ)+SG
!                EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
!    60       CONTINUE
!          END IF
!   100 CONTINUE
! C
! C     resulting absorption and emission coefficients
! C
!       EMIS(I0:I1)=EMIS(I0:I1)*II
!       ABSO(I0:I1)=ABSO(I0:I1)*POPI-EMIS(I0:I1)
! C
! C     pseudocontinuum for wavelengths smaller than the last observable
! C     line
! C
!   120 CONTINUE
!       if(lwph) GO TO 150
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
!       E=PLTEI*SGF
!       A=POPI*SGF
!       DO 130 IJ=I0,I1
!          F=FREQ(IJ)
!          IF(F.GE.FEDGE.OR.F.LT.FLST) GO TO 130
!          EMIS(IJ)=E*EXP(-4.79928E-11*F*T1)
!          ABSO(IJ)=A-EMIS(IJ)
!   130 CONTINUE
!   150 CONTINUE
! C
! C     add contributions for different spectral series
! C
!       ABSOH(I0:I1)=ABSOH(I0:I1)+ABSO(I0:I1)
!       EMISH(I0:I1)=EMISH(I0:I1)+EMIS(I0:I1)
!   200 CONTINUE
! C
! C     finally, multiply EMISH by 2h nu^3/c^2
! C
!       DO 210 IJ=I0,I1
!          F=FREQ(IJ)
!          F15=F*1.E-15
!          EMISH(IJ)=1.4743E-2*F15*F15*F15*EMISH(IJ)
!   210 CONTINUE
!       RETURN
!       END SUBROUTINE
      SUBROUTINE HYDLIN(ID,T,ANE,ABSOH,EMISH)
      use constants
      use MOD_HYDTAB,only: WLINE, MLINH,MHWL,NWLH,ILIN0
C
C     opacity and emissivity of hydrogen lines
C
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID
      real*8, intent(in   ) :: T,ANE
      real*8, intent(inout) :: ABSOH,EMISH
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      PARAMETER (FRH1=3.28805E15,FRH2=FRH1/4.,UN=1.,SIXTH=1./6.)
      PARAMETER (CPP=4.1412E-16,CCOR=0.09,CPJ=157803.,CLST=1.1E3)
      PARAMETER (C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.)
      PARAMETER (CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/(CLIGHT_SI*1e10))
      PARAMETER (CID1=0.01497)
      logical lwph,lquasi
      DIMENSION PJ(40),PRF0(36),OSCB(8),
     *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2IR
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
     $ IOPFE1
      COMMON/HYLPAR/IHYL,ILOWH,M10,M20
      COMMON/AUXVCS/XK,FXK,BETAD,DBETA
      !COMMON/BALVCS/PRFBAL(8,MDEPTH,36),WLBAL(8,36),NWLH(8)
      
      integer :: NWLBAL(MLINH)
      real*8  :: WLBAL(MLINH,MHWL)
      real*8  :: PRFBAL(MLINH,MDEPTH,MHWL)
      COMMON/HYDVCS/PRFBAL,WLBAL,NWLBAL
      
      common/wprob/wph(mdepth,40),acor(mdepth),lwph
      
      DATA OSCB   /0.6407,  0.1193,   0.04467,  0.02209,
     *             1.27D-2, 8.036D-3, 5.429D-3, 3.851D-3/
      DATA FRH    /3.289017E15/
      PARAMETER (HINV=1.5092973D26)
C
      WLINE(2,3:10)=(/6562.80, 4861.32, 4340.46, 4101.73,
     *             3970.07, 3889.05, 3835.38, 3797.90/)
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
      ANES=EXP(SIXTH*LOG(ANE))
C
C     populations of the first 40 levels of hydrogen
C
      ANP=POPUL(NKH,ID)
      PP=CPP*ANE*ANP*T1/SQT
      NLH=N1H-N0HN+1
      if(ifwop(n1h).lt.0) nlh=nlh-1
      DO IL=1,40
         X=IL*IL
         IF(IL.LE.NLH) PJ(IL)=POPUL(N0HN+IL-1,ID)/X
         IF(IL.GT.NLH) PJ(IL)=PP*EXP(CPJ/X*T1)
      ENDDO
      p2=pp*exp(cpj4*t1)
C
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      F00=C00*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID))
C

!      print*, 'Test Laiman', vturb(id)


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
      XII=UN/II
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
c        IF(I.EQ.1.AND.IOPHLI.GT.0.AND.J.LE.5) GO TO 100
         IF(I.EQ.2.AND.J.LE.10.AND.IBVCS.GT.0) THEN
          ILINE=ILIN0(I,J)
          NWL=NWLH(ILINE)
            DO 20 IWL=1,NWL
   20          PRF0(IWL)=PRFBAL(ILINE,ID,IWL)
            FID=CID*OSCB(ILINE)
            DO 50 IJ=I0,I1
               AL=ABS(WLAM(IJ)-WLINE(I,J))
               IF(AL.LT.1.E-4) AL=1.E-4
               AL=LOG10(AL)
               DO 30 IWL=1,NWL-1
                  IW0=IWL
                  IF(AL.LE.WLBAL(ILINE,IWL+1)) GO TO 40
   30          CONTINUE
   40          IW1=IW0+1
               PRF=(PRF0(IW0)*(WLBAL(ILINE,IW1)-AL)+PRF0(IW1)*
     *             (AL-WLBAL(ILINE,IW0)))/
     *             (WLBAL(ILINE,IW1)-WLBAL(ILINE,IW0))
               SG=EXP(PRF*AL10)*FID*wph(id,j)
               ABSO(IJ)=ABSO(IJ)+SG
               EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
   50       CONTINUE
          ELSE
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
               IF(FREQ(IJ).GT.FLST .and. .not.lwph) GO TO 60
               BETA=ABS(WLAM(IJ)-WL0)*FXK1
c              if(beta.gt.5.e4) go to 60
               SG=STARKA(BETA,BETAD,AD,DIV,TWO)*FID
               if(iophli.eq.2.and.i.eq.1.and.j.eq.2)
     *            sg=sg*feautr(fr,id)
               if(fid0.gt.0.) then
                  xd=beta/betad
                  if(xd.lt.5.) sg=sg+exp(-xd*xd)*fid0
               end if
               ABSO(IJ)=ABSO(IJ)+SG
               EMIS(IJ)=EMIS(IJ)+SG*PJ(J)
   60       CONTINUE
         END IF
  100 CONTINUE
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
      IF(I.EQ.1) THEN
         SGF=6.313E-18
       ELSE IF(I.EQ.2) THEN
         SGF=1.387E-17
       ELSE IF(I.EQ.3) THEN
         SGF=2.156E-17
       ELSE IF(I.EQ.4) THEN
         SGF=2.929E-17
       ELSE IF(I.EQ.5) THEN
         SGF=3.705E-17
       ELSE IF(I.EQ.6) THEN
         SGF=4.483E-17
      END IF
      E=PLTEI*SGF
      A=POPI*SGF
      DO 130 IJ=I0,I1
         F=FREQ(IJ)
         IF(F.GE.FEDGE.OR.F.LT.FLST) GO TO 130
         EMIS(IJ)=E*EXP(-4.79928E-11*F*T1)
         ABSO(IJ)=A-EMIS(IJ)
  130 CONTINUE
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

C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE HYLSET
      use SYNTHP_CONT,only:FREQC,NFCONT
      use constants
C
C     Initialization procedure for treating the hydrogen line opacity
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'SYNTHP.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/SYNTHP.FOR'
      DIMENSION ALB(15)
      COMMON/HYLPAR/IHYL,ILOWH,M10,M20
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2IR
      DATA ALB /656.28,486.13,434.05,410.17,397.01,
     *          388.91,383.54,379.79,377.06,375.02,
     *          373.44,372.19,371.20,370.39,369.72/
C
C     IHYL=-1  -  hydrogen lines are excluded a priori
C
      IHYL=-1
      if(iath.le.0) return
      IF(FREQC(NFCONT()).GE.3.28805E15) RETURN
      AL0=CLIGHT_SI*1e9/FREQC(1)
      AL1=CLIGHT_SI*1e9/FREQC(NFCONT())
      IF(grav.lt.6.) then
         IF(AL0.GT.160..AND.AL1.LT.364.6) RETURN
         IF(AL0.GT.506..AND.AL1.LT.630.) RETURN
         IF(AL0.GT.680..AND.AL1.LT.820.3) RETURN
C         IF(AL0.GT.1500.) RETURN
       else
         IF(AL0.GT.200..AND.AL1.LT.364.6) RETURN
         IF(AL0.GT.540..AND.AL1.LT.600.) RETURN
         IF(AL0.GT.720..AND.AL1.LT.820.3) RETURN
C         IF(AL0.GT.1500.) RETURN
      end if
C
C     otherwise, hydrogen lines are included
C
      IHYL=0
      M20=40
      IF(AL1.LT.364.) THEN
         ILOWH=1
         FRION=3.28805E15
         M10=SQRT(3.28805E15/(FRION-FREQC(NFCONT())))
         IF(FRION.GT.FREQ(1)) M20=SQRT(3.28805E15/
     *                        (FRION-FREQC(1)))
         IHYL=1
         IF(AL0.GT.123.) IHYL=0
         IF(AL0.GT.104..AND.AL1.LT.120.) IHYL=0
         IF(AL0.GT.98.5.AND.AL1.LT.102.) IHYL=0
         IF(IMODE.EQ.2.OR.ILVCS.GT.0.OR.GRAV.GE.6.) IHYL=1
       ELSE IF(AL1.LT.820.) THEN
         ILOWH=2
         FRION=8.2225E14
         M10=SQRT(3.289017E15/(FRION-FREQC(NFCONT())))
         IF(FRION.GT.FREQ(1)) M20=SQRT(3.289017E15/
     *                        (FRION-FREQC(1)))
         DO 10 I=1,15
            AL=ALB(I)
            IF(AL.LT.AL0-1..OR.AL.GT.AL1+1.) GO TO 10
            IHYL=1
            GO TO 20
   10    CONTINUE
   20    CONTINUE
         IF(IMODE.EQ.2.OR.IBVCS.GT.0.OR.GRAV.GE.6.) IHYL=1
CMH    ELSE IF(AL1.LT.1458.) THEN
CMH	 MORE EXACT wavelength!!! This could also be the case below
       ELSE IF(AL1.LT.1458.4299863152) THEN
         ILOWH=3
         FRION=3.6544142E14
	   if ((FRION-FREQC(NFCONT())) .lt. 0.) then
		print *,'hylset: give more precise wavelength!'
		print *, AL1, 1458.4299863152
		print *, 1.e8/FRION, ' must not be greater than ',1.e8/FREQC(NFCONT())
		pause
	   endif
         M10=SQRT(3.289017E15/(FRION-FREQC(NFCONT())))
         IF(FRION.GT.FREQ(1)) M20=SQRT(3.289017E15/
     *                        (FRION-FREQC(1)))
         IHYL=1
         IF(AL0.GT.1310.) IHYL=0
         IF(AL0.GT.1124..AND.AL1.LT.1250.) IHYL=0
         IF(AL0.GT.1035..AND.AL1.LT.1060.) IHYL=0
         IF(IMODE.EQ.2.OR.GRAV.GE.6.) IHYL=1
       ELSE IF(AL1.LT.2278.) THEN
         ILOWH=4
         FRION=2.0555837E14
	   if ((FRION-FREQC(NFCONT())) .lt. 0.) then
		print *,'hylset: give more precise wavelength!'
		print *, AL1, 2278.
		print *, 1.e8/FRION, ' must not be greater than ',1.e8/FREQC(NFCONT())
		pause
	   endif
         M10=SQRT(3.289017E15/(FRION-FREQC(NFCONT())))
         IHYL=1
       ELSE IF(AL1.LT.3281.) THEN
         ILOWH=5
         FRION=1.315589E14
	   if ((FRION-FREQC(NFCONT())) .lt. 0.) then
		print *,'hylset: give more precise wavelength!'
		print *, AL1, 3281.
		print *, 1.e8/FRION, ' must not be greater than ',1.e8/FREQC(NFCONT())
		pause
	   endif
         M10=SQRT(3.289017E15/(FRION-FREQC(NFCONT())))
         IHYL=1
       ELSE IF(AL1.LT.4466.) THEN
         ILOWH=6
         FRION=9.136394E13
	   if ((FRION-FREQC(NFCONT())) .lt. 0.) then
		print *,'hylset: give more precise wavelength!'
		print *, AL1, 4466.
		print *, 1.e8/FRION, ' must not be greater than ',1.e8/FREQC(NFCONT())
		pause
	   endif
         M10=SQRT(3.289017E15/(FRION-FREQC(NFCONT())))
         IHYL=1
       ELSE
         ILOWH=7
         FRION=6.7120228E13
	   if ((FRION-FREQC(NFCONT())) .lt. 0.) then
		print *,'hylset: FRION lt FREQ(2)'
		print *, AL1, 1458.4299863152
		print *, 1.e8/FRION, ' must not be greater than ',1.e8/FREQC(NFCONT())
		pause
	   endif
         M10=SQRT(3.289017E15/(FRION-FREQC(NFCONT())))
         IHYL=1
      END IF
c      WRITE(6,601) ILOWH,M20+1
c  601 FORMAT(1H0/ ' *** HYDROGEN LINES CONTRIBUTE'/
c     * '     THE NEAREST LINE ON THE SHORT-WAVELENGTH SIDE IS',
c     * I3,'  TO ',I3/)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C
C

      SUBROUTINE INIBLA
      use constants

C     =================
C
C     outprint of the selected lines
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
c      INCLUDE 'SYNTHP.FOR'
c      INCLUDE 'LINDAT.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
	INCLUDE '../inc/SYNTHP.FOR'
	INCLUDE '../inc/LINDAT.FOR'
      REAL*4 TYPION(9)
      logical test
      COMMON/LBLANK/IBLANK,NBLANK
      COMMON/IPRNTR/IPRIN
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
C
      DATA TYPION /4H I  ,4H II ,4H III,4H IV ,4H V  ,
     *             4H VI ,4H VII,4HVIII,4H IX /
      DATA APB,AP0,AP1,AP2,AP3,AP4 /4H    ,4H   .,4H   *,4H  **,4H ***,
     *                              4H****/
CMH	C1:
CMH	DP1: 2.* Boltzmann / ATOMIC MASS UNIT = 2.*1.38e-27/1.66054e-24 g
      PARAMETER (C1=2.3025851, C2=4.2014672, C3=1.4387886)
      PARAMETER (DP0=3.33564E-11, DP1=1.651E8,                   ! DP1 - only 1/2 as we are interested only in projection
     *           VW1=0.42, VW2=0.3, TENM4=1.E-4)
      PARAMETER (UN=1.)
C
c      IF(NLIN.EQ.0) RETURN
      XX=0.5*(FREQ(1)+FREQ(2))  !FIXME!!!
      BNU=BN*(XX*1.E-15)**3
      HKF=HK*XX
      DO 10 ID=1,ND
         T=TEMP(ID)
         ANE=ELEC(ID)
         EXH=EXP(HKF/T)
         EXHK(ID)=UN/EXH
         PLAN(ID)=BNU/(EXH-UN)
         STIM(ID)=UN-EXHK(ID)
         if(iath.gt.0) then
            ANP=POPUL(NKH,ID)
            AH=DENS(ID)/WMM/YTOT-ANP
          else
            ah=rrr(id,1,1)
         end if
         AHE=RRR(ID,1,2)
         VDWC(ID)=(AH+VW1*AHE)*(T*TENM4)**VW2
         DO 5 IAT=1,NATOMS
            IF(AMAS(IAT).GT.0.)
     *      DOPA1(IAT,ID)=UN/(XX*DP0*SQRT(DP1*T/AMAS(IAT)+VTURB(ID)))  ! 1/(Delta nu_D)
      
  !         print*, 'tets vdop', ID, VTURB(ID)

           if(IAT>30)write(991,*)  IAT,ID,DOPA1(IAT,ID)

  !          if(IAT>30) print*, 'mytest', IAT,ID,DOPA1(IAT,ID), 1.d-5*sqrt(DP1*T/AMAS(IAT)), 1.d-5*sqrt(VTURB(ID)) 


CMH	DOPA1: All velcities are calculated in units of dopa1
CMH	DP1*T/AMAS(IAT): Thermal broadening
    5 CONTINUE
   10 CONTINUE
C
C     for ID=IDSTD - output of selected line parameters
C
      IF(IPRIN.LT.-1) RETURN
c-test special print
c-opacity
c         test=.true.
c
         test=.false.
c***      ID=IDSTD
      do id = 2,ND-1,1
      ALM0=CLIGHT_SI*1e10/FREQ(1)
      ALM1=CLIGHT_SI*1e10/FREQ(2)
      if (id.eq.idstd) then
      IF(IMODE.GE.0) WRITE(6,601) IBLANK,ALM0,ALM1
      IF(IMODE.GE.0.OR.(IMODE.EQ.-1.AND.IBLANK.EQ.1)) WRITE(6,602)
      endif
      write (12,*) id
C
      DO 20 IL=1,NLIN0
         ALAM=CLIGHT_SI*1e10/FREQ0(IL)
         IAT=INDAT(IL)/100
         ION=MOD(INDAT(IL),100)
         AGAM = profil(IL,IAT,ID)
         ABCNT=EXP(GF0(IL)-EXCL0(IL)/TEMP(ID))*RRR(ID,ION,IAT)*
     *          DOPA1(IAT,ID)*STIM(ID)
         STR0=ABCNT/ABSTD(ID)
         GF=(GF0(IL)+C2)/C1
         EXCL=EXCL0(IL)/C3
         IF(STR0.LE.1.2) THEN
            WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577))
          ELSE
            WW1=SQRT(LOG(STR0))
         END IF
         IF(STR0.GT.55.) THEN
            WW2=0.5*SQRT(3.14*AGAM*STR0)
            IF(WW2.GT.WW1) WW1=WW2
         END IF
         EQW=ALAM/FREQ0(IL)*1.E3/DOPA1(IAT,ID)*WW1

         STR=EQW*10.
         APR=APB
         IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0
         IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1
         IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2
         IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3
         IF(STR.GE.1.E4) APR=AP4
         il0=0
CMH no line output
!SWITCH: Line Output
!         if (id.eq.idstd)
!      *   WRITE(6,603) IL0,IL,ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
!      *                STR0,EQW,APR,ilown(il),iupn(il)

        if ((id.eq.idstd) .and. (eqw .gt. 1)) then 
         WRITE(6,603) IL0,IL,ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
     *                STR0,EQW,APR,ilown(il),iupn(il)
        endif
c-test
c-opacity
         if (test) then
c         if (alam.gt.583.95 .and. alam.lt.584.75)
         if (alam.gt.303.45 .and. alam.lt.304.2)
c         if (alam.gt.1639. .and. alam.lt.1641.)
     *   WRITE(12,603) IL0,IL,ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
     *                STR0,EQW,APR,ilown(il),iupn(il)
         endif
   20 CONTINUE
c-test special print
      enddo
C
  601 FORMAT(1H1,' DATA FOR SELECTED LINES',I4,'. SET:',
     * ' INTERVAL  ',F9.3,' -',F9.3,' ANGSTROMS'/
     *          '  -----------------------')
  602 FORMAT(1H0/1H ,13X,
     * 'LAMBDA    ATOM    LOG GF       ELO    LINE/CONT',2X,
     * 'EQ.WIDTH'/)
  603 FORMAT(1H ,I3,I7,F10.3,2X,2A4,F7.2,F12.3,1PE11.2,0PF8.1,1X,A4,
     *2i3)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C
      SUBROUTINE INILIN(INLIST)
C     =========================
C
C     read in the input line list,
C     selection of lines that may contribute,
C     set up auxiliary fields containing line parameters,
C
C     Input of line data - unit 19:
C
C     For each line, one (or two) records, containing:
C
C    ALAM    - wavelength (in nm)
C    ANUM    - code of the element and ion (as in Kurucz-Peytremann)
C              (eg. 2.00 = HeI; 26.00 = FeI; 26.01 = FeII; 6.03 = C IV)
C    GF      - log gf
C    EXCL    - excitation potential of the lower level (in cm*-1)
C    QL      - the J quantum number of the lower level
C    EXCU    - excitation potential of the upper level (in cm*-1)
C    QU      - the J quantum number of the upper level
C    AGAM    = 0. - radiation damping taken classical
C            > 0. - the value of Gamma(rad)
C
C     There are now two possibilities, called NEW and OLD, of the next
C     parameters:
C     a) NEW, next parameters are:
C    GS      = 0. - Stark broadening taken classical
C            > 0. - value of log gamma(Stark)
C    GW      = 0. - Van der Waals broadening taken classical
C            > 0. - value of log gamma(VdW)
C    INEXT   = 0  - no other record necessary for a given line
C            > 0  - next record is read, which contains:
C    WGR1,WGR2,WGR3,WGR4 - Stark broadening values from Griem (in Angst)
C                   for T=5000,10000,20000,40000 K, respectively;
C                   and n(el)=1e16 for neutrals, =1e17 for ions.
C    ILWN    = 0  - line taken in LTE (default)
C            > 0  - line taken in NLTE, ILWN is then index of the
C                   lower level
C            =-1  - line taken in approx. NLTE, with Doppler K2 function
C            =-2  - line taken in approx. NLTE, with Lorentz K2 function
C    IUN     = 0  - population of the upper level in LTE (default)
C            > 0  - index of the lower level
C    IPRF    = 0  - Stark broadening determined by GS
C            < 0  - Stark broadening determined by WGR1 - WGR4
C            > 0  - index for a special evaluation of the Stark
C                   broadening (in the present version inly for He I -
C                   see procedure GAMHE)
C      b) OLD, next parameters are
C     IPRF,ILWN,IUN - the same meaning as above
C     next record with WGR1-WGR4 - again the same meaning as above
C     (this record is automatically read if IPRF<0
C
C     The only differences between NEW and OLD is the occurence of
C     GS and GW in NEW, and slightly different format of reading.
C
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
c      INCLUDE 'SYNTHP.FOR'
c      INCLUDE 'LINDAT.FOR'
      use constants

      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: INLIST
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
	INCLUDE '../inc/SYNTHP.FOR'
	INCLUDE '../inc/LINDAT.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE,CUTOF0,CUTOFS,TSTD,DSTD
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2IR
      COMMON/LBLANK/IBLANK,NBLANK
      COMMON/HE2PAR/IFHE2,IHE2L,ILWHE2,MHE10,MHE20
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
c      COMMON/NXTINI/ALM00,ALST00,NXTSET,INLIST
      DATA INLSET /0/
C
      PARAMETER (C1     = 2.3025851,
     *           C2     = 4.2014672,
     *           C3     = 1.4387886,
     *           CNM    = CLIGHT_CGS*1d7,
     *           ANUMIN = 1.9,
     *           ANUMAX = 90.41,
     *           AHE2   = 2.01,
     *           EXT0   = 3.17,
     *           UN     = 1.0,
     *           TEN    = 10.,
     *           HUND   = 1.D2,
     *           TENM4  = 1.D-4,
     *           TENM8  = 1.D-8,
     *           OP4    = 0.4,
     *           AGR0=2.4734E-22,
     *           XEH=13.595, XET=8065.9, XNF=25.,
     *           R02=2.5, R12=45., VW0=4.5E-9)
      PARAMETER (ENHE1=198310.76, ENHE2=438908.85)
C
      IL=0
      INNLT0=0
      IGRIE0=0
      DOPSTD=1.E7/ALAM0*DSTD
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD
      AVAB=ABSTD(IDSTD)*RELOP
      ASTD=0.01
      IF(GRAV.GT.6.) ASTD=0.1
      CUTOFF=CUTOF0
      ALAST=CNM/FRLAST
      IF(INLTE.GE.1.AND.INLSET.EQ.0) THEN
         CALL NLTSET(0,IL,IAT,ION,EXCL,EXCU,IEVEN,INNLT0)
         INLSET=1
      END IF
C
C     first part of reading line list - read only lambda, and
C     skip all lines with wavelength below ALAM0-CUTOFF
C
      ALAM=0.
    7 READ(19,*,END=100) ALAM
      IF(ALAM.LT.ALAM0-CUTOFF) GO TO 7
      BACKSPACE(19)
      GO TO 10
c
    8 write(6,688) alam
  688 format(' error in line list, lambda= ',f12.4/)
   10 ILWN=0
      IUN=0
      IPRF=0
      GS=0.
      GW=0.
      IF(INLIST.EQ.0) THEN
         READ(19,*,END=100,err=8) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
     *                        GS,GW,INEXT
cmh
cmh	if (ANUM.gt.100) then
c	print * ,ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
c     *                        GS,GW,INEXT
c	pause
cmh	endif

         IF(INEXT.NE.0) READ(19,*) ILWN,IUN,IPRF,WGR1,WGR2,WGR3,WGR4
       ELSE IF(INLIST.EQ.1) THEN
         READ(19,511,END=100) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
     *                        IPRF,ILWN,IUN
         IF(IPRF.LT.0) READ(19,*) WGR1,WGR2,WGR3,WGR4
  511    FORMAT(F10.4,2F6.2,F11.3,F4.1,F12.3,F4.1,E10.3,3I3)
       ELSE IF(INLIST.EQ.2) THEN
         READ(19,521,END=100) ALAM,GF,ANUM,EXCL,QL,XLBL,EXCU,QU
  521    FORMAT(F10.4,F7.3,F6.2,F11.3,F5.1,A11,F11.3,F5.1)
         AGAM=0.
      END IF
C
C     first selection : for a given interval a atomic number
C
!       db410=.false.
!       IF(ALAM.EQ.410.9447d0.and.GF==3d-1) THEN
!       PRINT *,ALAM,GF
!       pause
!       db410=.true.
!       ENDIF

      IF(ALAM.GT.ALAST+CUTOFF) GO TO 100
CMH	SET LINES ABOVE ATOMIC NUMBER 30 TO LITHIUM
!       IF(ANUM.GT.ANUMAX) THEN
! 	intanum = 0
! 	diff	= 0.00
! 	intanum = int(anum)
! 	diff = anum - intanum
! 	ANUM = 3.00 + diff
! 	endif
      IF(ANUM.LT.ANUMIN.OR.ANUM.GT.ANUMAX) GO TO 10
      IF(ABS(ANUM-AHE2).LT.TENM4.AND.IFHE2.GT.0) GO TO 10
C
C     second selection : for line strenghts
C
      FR0=CNM/ALAM   ! frequency
      IAT=ANUM  
      FRA=(ANUM-FLOAT(IAT)+TENM4)*HUND   ! hereafter help variable
      ION=INT(FRA)+1
      IF(ION.GT.IONIZ(IAT)) GO TO 10
      IEVEN=1
      IF(EXCL.GT.EXCU) THEN
         FRA=EXCL
         EXCL=EXCU
         EXCU=FRA
         FRA=QL
         QL=QU
         QU=FRA
         IEVEN=0
      END IF
      GFP=C1*GF-C2  
      EPP=C3*EXCL
      if (tstd .eq. 0.) then
        print *,'ATTENTION:'
        print *,'Standard formation height not properly defined'
        pause
      endif
      gx=gfp-epp/tstd
      ab0=0.
CMH	CRITERION TO EVALUATE LINES STILL ACTIVE
      if(gx.gt.-41)
     * AB0=EXP(GFP-EPP/TSTD)*RRR(IDSTD,ION,IAT)/DOPSTD/AVAB
      IF(AB0.LT.UN) GO TO 10
C
C     truncate line list if there are more lines than maximum allowable
C     (given by MLIN0 - see include file LINDAT.FOR)
C
      IL=IL+1
      IF(IL.GT.MLIN0) THEN
         WRITE(6,601) ALAM
	   PRINT *,'SYNSUB LINE-ARRAY not big enough',IL, MLIN0
	   pause
         IL=MLIN0
         GO TO 100
      END IF
C
C     =============================================
C     line is selected, set up necessary parameters
C     =============================================
C
C     evaluation of EXTIN0 - the distance (in delta frequency) where
C     the line is supposed to contribute to the total opacity
C
      IF(IAT.LE.2) THEN
         EXT=SQRT(10.*AB0)
       ELSE IF(IAT.LE.14) THEN
         EX0=AB0*ASTD*10.
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
       ELSE
         EX0=AB0*ASTD
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
      END IF
      EXTIN0=EXT*DOPSTD
C
C     store parameters for selected lines
C
      FREQ0(IL)=FR0
      EXCL0(IL)=EPP
      GF0(IL)=GFP
      EXTIN(IL)=EXTIN0
      INDAT(IL)=100*IAT+ION
C
C     indices for corresponding excitation temperatures of the lower
C     and upper levels
C
      if(excl.ge.enhe2) then
         ipotl(il)=3
       else if(excl.ge.enhe1) then
         ipotl(il)=2
       else
         ipotl(il)=1
      end if
c
      if(excu.ge.enhe2) then
         ipotu(il)=3
       else if(excu.ge.enhe1) then
         ipotu(il)=2
       else
         ipotu(il)=1
      end if
C
C     ****** line broadening parameters *****
C
C     1) natural broadening
C
      IF(AGAM.GT.0.) THEN
         GAMR0(IL)=EXP(C1*AGAM)
       ELSE
         GAMR0(IL)=AGR0*FR0*FR0
      END IF
C
C     if Stark or Van der Waals broadenig assumed classical,
C     evaluate the effective quantum number
C
      IF(GS.EQ.0..OR.GW.EQ.0) THEN
         Z=FLOAT(ION)
         XNEFF2=Z**2*(XEH/(ENEV(IAT,ION)-EXCU/XET))
         IF(XNEFF2.LE.0..OR.XNEFF2.GT.XNF) XNEFF2=XNF
      END IF
C
C     2) Stark broadening
C
      IF(GS.NE.0.) THEN
         GS0(IL)=EXP(C1*GS)
       ELSE
         GS0(IL)=TENM8*XNEFF2*XNEFF2*SQRT(XNEFF2)
      END IF
C
C     3) Van der Waals broadening
C
      IF(GW.NE.0.) THEN
         GW0(IL)=EXP(C1*GW)
       ELSE
         IF(IAT.LT.21) THEN
            R2=R02*(XNEFF2/Z)**2
          ELSE
            R2=(R12-FLOAT(IAT))/Z
                if (R2 .lt. 0.) then
                R2=0.
                print *,'Class. Van der Waals braoadening zero for IAT=',IAT
                endif
         END IF
         GW0(IL)=VW0*R2**OP4
      END IF
C
C     4) parameters for a special profile evaluation:
C
C     a) special He I and He II line broadening parameters
C
      ISPRFF=0
      IF(IAT.LE.2) ISPRFF=ISPEC(IAT,ION,ALAM)
      IF(IAT.EQ.2) CALL HESET(IL,ALAM,EXCL,EXCU,ION,IPRF,ILWN,IUN)
      ISPRF(IL)=ISPRFF
      IF(ISPRFF.GT.0) THEN
         IF(IL.GT.1.AND.ISPRFF.EQ.ISPRF(IL-1).AND.
     *      IHE2UV.GT.0) THEN
         CALL STOP('LINE LIST INCOMPATIBLE WITH SPEC. HE PROFILES')
         END IF
      END IF
      IPRF0(IL)=IPRF
C
C     b) parameters for Griem values of Stark broadening
C
      IF(IPRF.LT.0) THEN
         IGRIE0=IGRIE0+1
         IGRIEM(IL)=IGRIE0
         IF(IGRIE0.GT.MGRIEM) THEN
            WRITE(6,603) ALAM
            GO TO 20
         END IF
         WGR0(1,IGRIE0)=WGR1
         WGR0(2,IGRIE0)=WGR2
         WGR0(3,IGRIE0)=WGR3
         WGR0(4,IGRIE0)=WGR4
      END IF
   20 CONTINUE
C
C     implied NLTE option
C
      if(inlte.eq.-2.or.inlte.eq.12) then
          if(iat.le.20.and.excl.le.1000.) qu=-abs(qu)
       else if(inlte.eq.-3) then
          if(excl.le.1000.) qu=-abs(qu)
       else if(inlte.eq.-4) then
          qu=-abs(qu)
      end if
C
C     NLTE lines initialization
C
      INDNLT(IL)=0
      IF(QU.LT.0..OR.QL.LT.0.) THEN
         ILWN=-1
         QU=ABS(QU)
         QL=ABS(QL)
      END IF
C***  CONDITION FOR NON-LTE

      IF(ILWN.LT.0.AND.INLTE.NE.0) THEN
         INNLT0=INNLT0+1
         INDNLT(IL)=INNLT0
         IF(INNLT0.GT.MNLT) THEN
            WRITE(6,604) ALAM
            GO TO 100
         END IF
         GI=2.*QL+UN
         GJ=2.*QU+UN
         CALL NLTE(1,IL,ILWN,IUN,GI,GJ,EXCU)
         ILOWN(IL)=ILWN
         IUPN(IL)=IUN
      END IF
      IF(ILWN.GT.0.AND.INLTE.NE.0) THEN
         INNLT0=INNLT0+1
         INDNLT(IL)=INNLT0
         IF(INNLT0.GT.MNLT) THEN
            WRITE(6,604) ALAM
            GO TO 100
         END IF
         GI=2.*QL+UN
         GJ=2.*QU+UN
         CALL NLTE(1,IL,ILWN,IUN,GI,GJ,EXCU)
         ILOWN(IL)=ILWN
         IUPN(IL)=IUN
      END IF
      IF(INLTE.GE.1) THEN
	   CALL NLTSET(1,IL,IAT,ION,EXCL,EXCU,IEVEN,INNLT0)
         IF(INDNLT(IL).GT.0) THEN
            IF(INDNLT(IL).GT.MNLT) THEN
               WRITE(6,604) ALAM
               GO TO 100
            END IF
            GI=2.*QL+UN
            GJ=2.*QU+UN
            ILWN=ILOWN(IL)
            IUN=IUPN(IL)
CMH	if upper level LTE and lower level NLTE = ILWN = 0
		IF((ILWN.EQ.IUN) .OR.((ILWN.NE. 0.).and.(iun .eq.0.))) THEN
CMH-original     IF(ILWN.EQ.IUN) THEN
               INDNLT(IL)=0
               ILOWN(IL)=0
               IUPN(IL)=0
             ELSE
               CALL NLTE(1,IL,ILWN,IUN,GI,GJ,EXCU)
            END IF
         END IF
      END IF
      GO TO 10
  100 continue
C
      NLIN0=IL
C
      NSP=0
      DO 260 IL=1,NLIN0
         ISP=ISPRF(IL)
         IF(ISP.GT.5) THEN
            NSP=NSP+1
            ISP0(NSP)=ISP
         END IF
  260 CONTINUE
C
      NNLT=INNLT0
      NGRIEM=IGRIE0
C
      WRITE(6,611) NLIN0,NNLT,NGRIEM
  611 FORMAT(/' LINES - TOTAL        :',I10
     *       /' LINES - NLTE         :',I10
     *       /' LINES - GRIEM        :',I10/)
  601 FORMAT('0 **** MORE LINES THAN MLIN0, LINE LIST TRUNCATED !!'/
     *'       AT LAMBDA',F15.4,'  NM'/)
  602 FORMAT('0 **** MORE LINES WITH SPECIAL PROFILES THAN MPRF'/
     *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
  603 FORMAT('0 **** MORE LINES WITH GRIEM PROFILES THAN MGRIEM'/
     *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
  604 FORMAT('0 **** MORE LINES IN NLTE OPTION THAN MNLT'/
     *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE INIMOD
            USE MOD_STATEQ

C
C   SET UP COMMON/RRRVAL/  - VALUES OF  N(ION)/U(ION) FOR ALL THE ATOMS
C   AND IONS CONSIDERED
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      COMMON/BLAPAR/RELOP,SPACE,CUTOF0,CUTOFS,TSTD,DSTD
      COMMON/HPOPST/HPOP
C
C
c     1. "low-temperature" ionization fractions
c         (using Hamburg partition functions)
c
      if(teff0.ge.0.) then
!     open (unit=500,file='relab.out')
      DO 50 ID=1,ND
         do ii=1,nleve0
! 		if (id .eq. ND) then
! 			WRITE (500,*)		 II,POPUL(II,ID),RELAB(II)
! 			WRITE (500,*)		 II,POPUL(II,ID)*RELAB(II)
! c			PRINT *, 'INIMOD: ', II,POPUL(II,ID),RELAB(II)
! c			pause
! 		endif
            popul(ii,id)=popul(ii,id)*relab(ii)
         end do
         sum=0.
         do 1 ii=1,nleve0
    1       sum=sum+popul(ii,id)
         dens(id)=sum*wmm
         sum=0.
c         do 2 ii=n0h,nkh
c    2       sum=sum+popul(ii,id)
c         hpop=sum
C
         CALL STATE(ID)
         HPOP=DENS(ID)/WMM/YTOT
c         AH=HPOP-POPUL(NKH,ID)
         DO J=1,MION
           DO I=1,MATOM
             RRR(ID,J,I)=RR(I,J)*HPOP
           ENDDO
         ENDDO
         if(TEMP(ID)==0.) then
          print *,'synsubm: inimod: TSTD = 0, ID =',ID
         endif
         IF(ID.NE.IDSTD) GO TO 50
         TSTD=TEMP(ID)
         VTS=VTURB(ID)
         DSTD=SQRT(1.4E7*TSTD+VTS)
         WRITE(6,601) ID,TEMP(ID),ELEC(ID)
c-pr
c         write (6,*) ' INIMOD in SYNSUB: print skipped'
         DO 20 I=2,MATOM
cc ivan used MI1   20       WRITE(6,602) TYPAT(I),(RR(I,J),J=1,MI1)
cc ws 4-6-96 9 is highest or?
   20       WRITE(6,602) TYPAT(I),(RR(I,J),J=1,9)
         WRITE(6,603)
c-pr
c         write (6,*) ' INIMOD in SYNSUB: print skipped'
         DO 30 I=2,MATOM
cc ws   30       WRITE(6,602) TYPAT(I),(PFSTD(J,I),J=1,MI1)
   30       WRITE(6,602) TYPAT(I),(PFSTD(J,I),J=1,9)
   50 CONTINUE
!	close (500)
c
c     2. "high-temperature" ionization fractions
c         (using the Opacity Project ionization fractions)
c
      else
      do id=1,nd
         do ion=1,mion
            do iat=1,matom
               rrr(id,ion,iat)=0.
            end do
         end do
      end do
      CALL FRAC1
      ID=IDSTD
      HPOP=DENS(ID)/WMM/YTOT
      WRITE(6,604) ID,TEMP(ID),ELEC(ID)
      DO 60 I=1,MATOM
         WRITE(6,605) TYPAT(I),(RRR(ID,J,I)/hpop,J=1,MION)
         ioniz(i)=i+1
   60 continue
      end if
C
  601 FORMAT(1H1,' STANDARD DEPTH =',I3,'   T,NE =',F8.1,1PE12.3
     *         /'  --------------'
     *       //' N/U  AT THE STANDARD DEPTH'/
     *         ' --------------------------'//)
  602 FORMAT(1H ,A4,(1P9E9.2))
  603 FORMAT(1H0/' PARTITION FUNCTIONS AT THE STANDARD DEPTH'/
     *          ' ------------------------------------------'//)
  604 FORMAT(/' N/U  AT THE STANDARD DEPTH  - OP DATA',
     * '  (ID =',I3,' ; T,Ne = ',F8.1,1PE12.3,' )'//)
  605 FORMAT(1H ,A4,(1P8E9.2))
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C

      SUBROUTINE INTHE2(W0,X0,Z0,IWL,ILINE)
C
C     Interpolation in temperature and electron density from the
C     Schoening and Butler tables for He II lines to the actual
C     actual values of temperature and electron density
C
C     This procedure is quite analogous to INTVCS for Balmer lines
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IWL,ILINE
      real*8, intent(in   ) :: X0,Z0
      real*8, intent(  out) :: W0
      PARAMETER (UN=1.D0)
      COMMON/HE2DAT/WLHE(36,19),XTHE(6),XNEHE(11,19),PRFHE(36,6,11),
     *              NWLHE,NTHE,NEHE
      COMMON/AUXHE2/XKE,FXKE,BETADE,DBETAE
      DIMENSION ZZ(3),XX(3),WX(3),WZ(3)
C
      NX=3
      NZ=3
C
C     for values lower than the lowest grid value of electron density,
C     and higher values of temperature,
C     the profiles are determined by the approximate expression
C     (see STARKA); not by an extrapolation in the tables which may
C     be very inaccurate
C
c     IF(Z0.LT.XNEHE(1,ILINE)*0.99) THEN
      IF(Z0.LT.XNEHE(1,ILINE)*0.99.OR.X0.GT.1.01*XTHE(NTHE)) THEN
         CALL DIVHE2(BETADE,AE,DIVE)
         W0=STARKA(WLHE(IWL,ILINE)/FXKE,BETADE,AE,DIVE,UN)*DBETAE
         W0=LOG10(W0)
         return
      END IF
C
C     Otherwise, one interpolates (or extrapolates for higher than the
C     highes grid value of electron density) in the Schoening and
C     Butler tables
C
      DO 10 IZZ=1,NEHE-1
         IPZ=IZZ
         IF(Z0.LE.XNEHE(IZZ+1,ILINE)) GO TO 20
   10 CONTINUE
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NEHE-NZ+1) N0Z=NEHE-NZ+1
      N1Z=N0Z+NZ-1
C
      DO 300 IZZ=N0Z,N1Z
         I0Z=IZZ-N0Z+1
         ZZ(I0Z)=XNEHE(IZZ,iline)
C
C     Both interpolations (in T as well as in electron density) are
C     by default the quadratic interpolations in logarithms
C
         DO 30 IX=1,NTHE-1
            IPX=IX
            IF(X0.LE.XTHE(IX+1)) GO TO 40
   30    CONTINUE
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NTHE-NX+1) N0X=NTHE-NX+1
         N1X=N0X+NX-1
         DO 200 IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XTHE(IX)
            WX(I0)=PRFHE(IWL,IX,IZZ)
  200       CONTINUE
         WZ(I0Z)=YINT(XX,WX,X0)
  300 CONTINUE
      W0=YINT(ZZ,WZ,Z0)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C
C ********************************************************************
C ********************************************************************
C
C

      FUNCTION ISPEC(IAT,ION,ALAM)
C
C     Auxiliary procedure for INISET
C
C     Input:  IAT  - atomic number
C             ION  - ion (=1 for neutrals, =2 for once ionized, etc.)
C             ALAM - wavelength in nanometers
C     Output: ISPEC - parameter specifying whether the given line
C                     is taken with a special (pretabulated) absorption
C                     profile - only for hydrogen and helium
C                   = 0  - profile is taken as an ordinary Voigt profile
C                   > 0  - special profile
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IAT,ION
      real*8, intent(in   ) :: ALAM
      integer*8  :: ISPEC
        
      COMMON/DETLIN/ILVCS,IBVCS,IHE1,IHE144,IHE2UV,IHE2VI,IHE2RE
C
      ISPEC=0
      IF(IAT.GT.2) RETURN
C
      IF(IAT.EQ.1) THEN
         ISPEC=1
         RETURN
       ELSE
         IF(ION.EQ.1) THEN
            IF(ABS(ALAM-447.1).LT.0.5.AND.IHE144.GT.0) THEN
              ISPEC=2
              PRINT *, IAT, ION, ALAM, IHE144
            ENDIF
            IF(ABS(ALAM-438.8).LT.0.2.AND.IHE1.GT.0) ISPEC=3
            IF(ABS(ALAM-402.6).LT.0.2.AND.IHE1.GT.0) ISPEC=4
            IF(ABS(ALAM-492.2).LT.0.2.AND.IHE1.GT.0) ISPEC=5
          ELSE
C
            IF(ALAM.LT.163..OR.ALAM.GT.1012.7) RETURN
            IF(ALAM.LT.321.) THEN
               IF(ABS(ALAM-164.0).LT.0.2.AND.IHE2UV.GT.0) ISPEC=6
               IF(ABS(ALAM-320.3).LT.0.2.AND.IHE2UV.GT.0) ISPEC=7
               IF(ABS(ALAM-273.3).LT.0.2.AND.IHE2UV.GT.0) ISPEC=8
               IF(ABS(ALAM-251.1).LT.0.2.AND.IHE2UV.GT.0) ISPEC=9
               IF(ABS(ALAM-238.5).LT.0.2.AND.IHE2UV.GT.0) ISPEC=10
               IF(ABS(ALAM-230.6).LT.0.2.AND.IHE2UV.GT.0) ISPEC=11
               IF(ABS(ALAM-225.3).LT.0.2.AND.IHE2UV.GT.0) ISPEC=12
             ELSE IF(ALAM.LT.541.) THEN
               IF(ALAM.LT.392.3) RETURN
               IF(ABS(ALAM-468.6).LT.0.2.AND.IHE2VI.GT.0) ISPEC=13
               IF(ABS(ALAM-485.9).LT.0.2.AND.IHE2VI.GT.0) ISPEC=14
               IF(ABS(ALAM-454.2).LT.0.2.AND.IHE2VI.GT.0) ISPEC=15
               IF(ABS(ALAM-433.9).LT.0.2.AND.IHE2VI.GT.0) ISPEC=16
               IF(ABS(ALAM-420.0).LT.0.2.AND.IHE2VI.GT.0) ISPEC=17
               IF(ABS(ALAM-410.0).LT.0.2.AND.IHE2VI.GT.0) ISPEC=18
               IF(ABS(ALAM-402.6).LT.0.2.AND.IHE2VI.GT.0) ISPEC=19
               IF(ABS(ALAM-396.8).LT.0.2.AND.IHE2VI.GT.0) ISPEC=20
               IF(ABS(ALAM-392.3).LT.0.2.AND.IHE2VI.GT.0) ISPEC=21
             ELSE
               IF(ABS(ALAM-1012.4).LT.0.2.AND.IHE2RE.GT.0) ISPEC=22
               IF(ABS(ALAM-656.0).LT.0.2.AND.IHE2RE.GT.0) ISPEC=23
               IF(ABS(ALAM-541.2).LT.0.2.AND.IHE2RE.GT.0) ISPEC=24
            END IF
         END IF
      END IF
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION XK2DOP(TAU)
C
C     KERNEL FUNCTION K2  (AUXILIARY PROCEDURE TO NLTE)
C     AFTER  HUMMER,  1981, J.Q.S.R.T. 26, 187
C
c      INCLUDE 'PARAMS.FOR'
      implicit real*8(a-h,o-z)
      real*8, intent(in   ) :: TAU
      real*8 :: XK2DOP
	INCLUDE '../inc/PARAMS.FOR'
      DATA PI2SQ,PISQ /2.506628275D0,  1.772453851D0/
      DATA A0,A1,A2,A3,A4 /
     *  1.D0,  -1.117897000D-1,  -1.249099917D-1,  -9.136358767D-3,
     *         -3.370280896D-4/
      DATA B0,B1,B2,B3,B4,B5 /
     *  1.D0,   1.566124168D-1,   9.013261660D-3,   1.908481163D-4,
     *         -1.547417750D-7,  -6.657439727D-9/
      DATA C0,C1,C2,C3,C4 /
     *  1.0D0,   1.915049608D01,   1.007986843D02,   1.295307533D02,
     *         -3.143372468D01/
      DATA D0,D1,D2,D3,D4,D5/
     *  1.D0,   1.968910391D01,   1.102576321D02,   1.694911399D02,
     *         -1.669969409D01,  -3.666448000D01/
      XK2DOP=1.D0
      IF(TAU.LE.0.) RETURN
      IF(TAU.GT.11.) GO TO 10
      P=A0+TAU*(A1+TAU*(A2+TAU*(A3+TAU*A4)))
      Q=B0+TAU*(B1+TAU*(B2+TAU*(B3+TAU*(B4+TAU*B5))))
      XK2DOP=TAU/PI2SQ*LOG(TAU/PISQ)+P/Q
      RETURN
   10 X=1.D0/LOG(TAU/PISQ)
      P=C0+X*(C1+X*(C2+X*(C3+X*C4)))
      Q=D0+X*(D1+X*(D2+X*(D3+X*(D4+X*D5))))
      XK2DOP=P/Q/2.D0/TAU/SQRT(LOG(TAU/PISQ))
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
C
!       SUBROUTINE LINOP(ID,ABLIN,EMLIN)
! C     =====================================
! C
! C     TOTAL LINE OPACITY (ABLIN)  AND EMISSIVITY (EMLIN)
! C
!       implicit real*8 (a-h,o-z)
!       integer,intent(in   ) :: ID
!       !real*8, intent(in   ) :: AVAB
!       real*8, intent(inout) :: ABLIN,EMLIN
!         
! 	INCLUDE '../inc/PARAMS.FOR'
! 	INCLUDE '../inc/MODELP.FOR'
! 	INCLUDE '../inc/SYNTHP.FOR'
! 	INCLUDE '../inc/LINDAT.FOR'
!       LOGICAL LPR,lvi,lne,ltrad
!       PARAMETER (UN=1., EXT0 = 3.17,  TEN = 10.)
!       PARAMETER (FRHE1=5.945E15, FRHE2=1.316D16)
!       DIMENSION ABLIN(MFREQ),EMLIN(MFREQ),ABLINN(MFREQ),
!      *          XFR(MFREQ)
!       COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
!       common/linoff/velmax,nltoff,iemoff
! C
! CMH	DETERMINE DEPTH INDEX OF TEMPERATURE MINIMUM
! 	DO L=2,MDEPTH-1
! 	IF ((TEMP(L) .LT. TEMP(L-1)) .AND. (TEMP(L) .LT. TEMP(L+1))
!      $	.and. (temp(l) .ne. 0)) THEN
! 				 NDPMIN = L
! c			PRINT *,'LINOP: TEMPERATURE MINIMUM AT DEPTH POINT=',
! c    $				L,TEMP(NDPMIN)
! 		ENDIF
! 	ENDDO
! 	if (NDPMIN .eq. 0) THEN
! 		NDPMIN = 1
! 	endif
! CMH	check temperature minimum - less than 5000.
! 	if (TEMP(NDPMIN) .gt. 5000.) then
! 		print *,'something wrong with Temperature minimum!',
!      $                NDPMIN, TEMP(NDPMIN)
! 		write (6,*)'something wrong with Temperature minimum!',
!      $                NDPMIN, TEMP(NDPMIN)
! 		pause
! 	endif
! c  	PRINT *,'LINOP: TEMPERATURE MINIMUM AT DEPTH POINT=',
! c     $				NDPMIN,TEMP(NDPMIN)
!       DO 10 IJ=1,NFREQ
!          ABLIN(IJ)=0.
!          ABLINN(IJ)=0.
!    10    EMLIN(IJ)=0.
! C
!       IF(NLIN0.EQ.0) RETURN
! C
! C     overall loop over contributing lines
! C
!       TEM=TEMP(ID)
!       TEM1=UN/TEMP(ID)
! C ATTENTION !!!!!!
!       ! wdil(id)=1.
! c	PRINT *,'WDIL(ID) = ',WDIL(ID)
! C	PRINT *, 'ATTENTION: IN SYNSUBM.FOR WDIL SET TO ONE!!'
! c
!       lvi=Velw(id).gt.velmax.and.iemoff.eq.0
!       lne=velw(id).gt.velmax.and.nltoff.gt.0.and.iemoff.gt.0
!       ltrad=trad(1,id).eq.trad(2,id).and.trad(1,id).eq.trad(3,id)
!       ltrad=ltrad.and.trad(1,id).eq.temp(id)
! C ATTENTION !!!!!!!
! c     ltrad=.true.
! C ATTENTION !!!!!!!
!       XX=(FREQ(1)+FREQ(NFREQ))*5.e-16
!       BNU=BN*XX*XX*XX
! c     if(xx.ge.frhe2) then
! c        itrad=3
! c      else if(xx.ge.frhe1) then
! c        itrad=2
! c      else
! c        itrad=1
! c     end if
! c     planw=bnu/(exp(hk*xx*1.e15/trad(itrad,id))-1.)*wdil(id)
!       planw=xjcon(id)
!       DO 100 IL=1,NLIN0
!          INNLT=INDNLT(IL)
!          trl=trad(ipotl(il),id)
! c        tru=trad(ipotu(il),id)
!          tem1=un/trl
!          wdid=wdil(id)
! c	   print*,id,'in linop: wdid = ' ,wdid
! c	   pause
! c
! c     preparation for rejecting lines for v > velmax
! c
!          if(lvi) then
!             if(innlt.eq.0) then
!                go to 100
!              else
!                if(nltoff.ne.0) go to 100
!             end if
!          end if
! c
!          IAT=INDAT(IL)/100
!          ION=MOD(INDAT(IL),100)
!          LPR=.TRUE.
!          ISP=ISPRF(IL)
!          IF(ISP.GT.1.AND.ISP.LE.5) LPR=.FALSE.
!          IF (ISP.GE.6) GO TO 100
!          AGAM = profil(IL,IAT,ID)
!          DOP1=DOPA1(IAT,ID)
!          FR0=FREQ0(IL)
!          IF(INNLT.EQ.0) THEN
!             if(excl0(il).lt.2000.) wdid=1.
!             AB0=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*
!      *          DOP1*STIM(ID)*wdid
!             if(.not.ltrad) then
!                sl0=planw
! c              x=(hk*fr0+excl0(il)*tem1*(trl-tru))/tru
! c              sl0=bnu/(exp(x)-un)
!             end if
!           ELSE IF(INNLT.GT.0) THEN
!             AB0=ABCENT(INNLT,ID)
!             SL0=SLIN(INNLT,ID)
!           ELSE
!             BI=UN
!             BJ=UN
!             II=ILOWN(IL)
!             JJ=IUPN(IL)
!             IF(II.GT.0) BI=BFAC(II,ID)
!             IF(JJ.GT.0) BJ=BFAC(JJ,ID)
!             STIMB=BI*(1.-EXHK(ID)*BJ/BI)
!             AB0=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*
!      *          DOP1*STIMB
!             SL0=PLAN(ID)*(UN-EXHK(ID))/(BI/BJ-EXHK(ID))
!          END IF
! C
! C        set up limiting frequencies where the line I is supposed to
! C        contribute to the opacity
! C
! c        EX0=AB0/AVAB*AGAM
! c        EXT=EXT0
! c        IF(EX0.GT.TEN) EXT=SQRT(EX0)
! c        EXT=EXT/DOP1
! c        FR1=FR0+EXT
! c        FR2=FR0-EXT
! C
! C        set up array of dimensionless frequencies - x
! C        and evaluate limiting frequency indices where the line
! C        had a non-negligible opacity
! C
!          IJ1=1
!          IJ2=NFREQ
!          DO 20 IJ=1,NFREQ
!             FR=FREQ(IJ)
!             XFR(IJ)=ABS(FR-FR0)*DOP1
!    20    CONTINUE
! 
!    
! C
!          IF(INNLT.EQ.0.and.ltrad) THEN
! C
! C        *********
! C        LTE lines
! C        *********
! C
! 	 IF(LPR) THEN
! C
!             DO 40 IJ=IJ1,IJ2
!                ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XFR(IJ))
!    40       CONTINUE
! C
! C        special expressions for 4 selected He I lines
! C
!          ELSE
!             DO 60 IJ=1,NFREQ
!                FR=FREQ(IJ)
!                 ABL=AB0*PHE1(ID,FR,ISP-1)
!                ABLIN(IJ)=ABLIN(IJ)+ABL
!    60       CONTINUE
!          END IF
! C
! C        **********
! C        NLTE LINES
! C        **********
! C
!        ELSE
! 	 IF(LPR) THEN
! C
!             DO 80 IJ=IJ1,IJ2
! !                write (990,'(i3, "ABL=",e12.5," AGAM=",e12.5," AB0="'//
! !      $               ',e12.5,"voigtk=",e12.5)'),
! !      $                ij,ab0,voigtk(agam,xfr(ij))
! 
!                ABL=AB0*VOIGTK(AGAM,XFR(IJ))
!                ABLINN(IJ)=ABLINN(IJ)+ABL
!                if(lne) go to 80
!                EMLIN(IJ)=EMLIN(IJ)+ABL*SL0
!    80       CONTINUE
! C
! C        again, special expressions for 4 selected He I lines
! C
!          ELSE
!             DO 90 IJ=1,NFREQ
!                FR=FREQ(IJ)
!                ABL=AB0*PHE1(ID,FR,ISP-1)
!                ABLINN(IJ)=ABLINN(IJ)+ABL
!                if(lne) go to 90
!                EMLIN(IJ)=EMLIN(IJ)+ABL*SL0
!    90       CONTINUE
!          END IF
!       END IF
! 
!   100 CONTINUE
! C
! C     LTE line emissivity
! C
! CMH	IF THE DEPTPOINT IS OUTWARD OF THE TEMPERATURE MIMIMUM
! CMH	THE PLANCK-FUNCION (OR TEMPERATURE) IS SET TO THE
! CMH	VALUE AT THE TEMPERATURE MINIMUM (ID = 55)
! C	print *,ID,'LINOP: USE THE PLANCK FUNKTION OF ID=55 (FOR MODEL C)'
! 
!       if(velw(id).le.velmax) then
!       DO 110 IJ=1,NFREQ
! 		if (id .lt. NDPMIN) then
! CMH TO DO: automatic detection of the temperature minimum
! CMH TO DO: when FAL is read in!!!!
!           EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(NDPMIN)*WDIL(ID)
! c		print *, 'ID TMIN is set to NDPMIN!!!!'
! c		print *, 'PLAN(NDPMIN)',NDPMIN, PLAN(NDPMIN)
! c		print *, 'PLAN(ID)',NDPMIN, PLAN(ID)
! c		pause
! c		print *,'Automatic detection of tmin necessary!!'
! 		else
!          EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID)*WDIL(ID)
! 		endif
!   110 CONTINUE
!       end if
! 
! C
! 
! C     add opacity of NLTE lines
! C
!       DO 115 IJ=1,NFREQ
!          ABLIN(IJ)=ABLIN(IJ)+ABLINN(IJ)
!   115 CONTINUE
! C
! C     special routine for selected He II lines
! C
!       IF(NSP.EQ.0) RETURN
!       DO 120 IS=1,NSP
!          ISP=ISP0(IS)
!          if(velw(id).gt.velmax.and.nltoff.gt.0) go to 120
!          IF(ISP.GE.6.AND.ISP.LE.24) CALL PHE2(ISP,ID,ABLIN,EMLIN)
!   120 CONTINUE
! C
! 
!       RETURN
!       END SUBROUTINE

C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE NLTE(MODE,IL,ILW,IUN,GI,GJ,EXCU)
C     ===========================================
C
C     Control procedure for the NLTE option
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
c      INCLUDE 'LINDAT.FOR'
      use constants,only:CLIGHT_CGS,CLIGHT_SI
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: MODE,IL,ILW,IUN
      real*8, intent(in   ) :: GI,GJ,EXCU

	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
	INCLUDE '../inc/LINDAT.FOR'
      DIMENSION POPNLT(MLEV,MDEPTH),GNLT(MLEV),
     *          NPLUS(MLEV),IONT(MLEV),
     *          GION(MIOEX),EION(MIOEX)
      real*8,parameter :: C3 = 1.4387886
      
C
C *** For  MODE = 1, and only for INLTE > 1, input f the additional
C        level parameters  and NLTE level populations
C
	PRINT *,'CALLED NLTE CONTROL PROCEDURE!!!!!'
      IF(MODE.GT.0) GO TO 50
      IF(INLTE.LE.5) RETURN
      INLADD=INLTE
!      print*, inladd
!      stop
      READ(INLADD,*,end=40,err=40) NPOP,NIONL,REL
      IF(REL.LE.0.) REL=1.
      WRITE(6,600) NPOP,NIONL,REL
      IF(NPOP.LE.0.OR.NIONL.LE.0) RETURN
      DO 10 ION=1,NIONL
         READ(INLADD,*) GION(ION),EION(ION)
         WRITE(6,601) ION,GION(ION),EION(ION)
   10 CONTINUE
      DO 20 I=1,NPOP
         READ(INLADD,*) GNLT(I),NPLUS(I),IONT(I)
         WRITE(6,602) I,GNLT(I),NPLUS(I),IONT(I)
   20 CONTINUE
      DO 30 J=1,ND
         READ(INLADD,*) (POPNLT(I,J),I=1,NPOP)
         DO 25 I=1,NPOP
   25       POPNLT(I,J)=POPNLT(I,J)*REL
         WRITE(6,603) J,(POPNLT(I,J),I=1,NPOP)
   30 CONTINUE
   40 RETURN
   50 CONTINUE
  600 FORMAT(1H1,'INPUT PARAMETERS FOR THE NLTE OPTION'/
     *           1H ,'------------------------------------'//
     *      ' NPOP =',I3/' NIONL =',I3/' REL  =',1PE10.3//)
  601 FORMAT(1H ,'ION =',I3,'    GION =',F10.2,'   EION =',F15.3)
  602 FORMAT(1H ,'LEVEL =',I3,'    G =',F10.2,'    NPLUS =',I3,'  IONT
     *=',I3)
  603 FORMAT(1H0,'POPULATIONS FOR ID =',I3/(1H ,1P10E10.3))
  605 FORMAT(1H0/)
C
C   -----------------------------------------------------------------
C
C *** FOR  MODE.GT.1  -  CALCULATION OF THE
C     CENTRAL OPACITY (ABCENT) AND THE LINE SOURCE FUNCTION  (SLIN)
C
      ILNLT=INDNLT(IL)

      IF(ILNLT.LE.0) RETURN
      IAT=INDAT(IL)/100
      ION=MOD(INDAT(IL),100)
      EGF=EXP(GF0(IL))
c-pr
c      print *,' egf =',egf
      BNU=BN*(FREQ0(IL)*1.E-15)**3
      DP0=3.33564E-11*FREQ0(IL)
      DP1=1.651E8/AMAS(IAT)
      IF(ILW.LE.0) GO TO 100
      IF(INLTE.LE.5) THEN
C
C     INLTE <=5 - line is a transition between explicit levels of the
C                 input model
C
	 N0I=NFIRST(IEL(ILW))
         NKI=NNEXT(IEL(ILW))
	 E1=ENION(N0I)
	 EL=(E1-ENION(ILW))/BOLK
C	 EU=(E1-ENION(IUN))/BOLK
         E1=E1/BOLK
         EXCU3=EXCU*C3
c-pr
         print *,' NLTE in SYNSUB:'
c	   print *, iun, g(iun)
c         print *,' ilw, g_ilw',ilw,g(ilw)
c         print *,' iun, g_iun',iun,g(iun)
         DO 60 ID=1,ND
            T=TEMP(ID)
C	    PI=POPUL(ILW,ID)/G(ILW)*EXP(-(EXCL0(IL)-EL)/T)
	    PI=POPUL(ILW,ID)/G(ILW)
            IF(IUN.GT.0) THEN
C              PJ=POPUL(IUN,ID)/G(IUN)*EXP(-(EXCU3-EU)/T)
               PJ=POPUL(IUN,ID)/G(IUN)
             ELSE
               PP=POPUL(NKI,ID)*ELEC(ID)/T/SQRT(T)*2.0706E-16/G(NKI)
               PJ=PP*EXP((E1-EXCU3)/T)
C              PJ=PP
            END IF
            X=PI/PJ
            DOP=DP0*SQRT(DP1*T+VTURB(ID))
            SLIN(ILNLT,ID)=BNU/(X-1.)
            ABCENT(ILNLT,ID)=PI*(1.-1./X)*EGF/DOP
   60    CONTINUE
       ELSE IF(INLTE.gt.5) THEN
C
C     INLTE > 5 - line is not a transition between explicit levels of the
C                 input model
C                 corresponding NLTE population were read and stored
C                 at array POPNLT
C
         IONN=IONT(ILW)
         DO 70 ID=1,ND
            T=TEMP(ID)
            PLOW=POPNLT(ILW,ID)*GI/GNLT(ILW)
            PION=POPNLT(NPLUS(ILW),ID)*GION(IONN)/GNLT(NPLUS(ILW))
            PP=PION*ELEC(ID)/T/SQRT(T)*2.0706E-16/GION(IONN)
            PLTE=PP/EXP((EXCL0(IL)-EION(IONN)*C3)/T)*GI
            BI=PLOW/PLTE
            BJ=1.
            IF(IUN.GT.0) THEN
               PUPP=POPNLT(IUN,ID)*GJ/GNLT(IUN)
               PLTE=PP/EXP((EXCU-EION(IONN))*C3/T)*GJ
               BJ=PUPP/PLTE
            END IF
            EXPB=FREQ0(IL)/T*4.79928E-11
            EXPB=EXP(-EXPB)*BJ/BI
            DOP=DP0*SQRT(DP1*T+VTURB(ID))
            ABCENT(ILNLT,ID)=PLOW*(1.-EXPB)*EGF/DOP/GI
            SLIN(ILNLT,ID)=BNU*EXPB/(1.-EXPB)
   70    CONTINUE
      END IF
      RETURN
C
C     Approximate NLTE for resonance lines - second order escape
C     probablity theory form of the source function
C
C     Optical depth scale
C
  100 CONTINUE
      ALMIL=CLIGHT_CGS*1d7/FREQ0(IL)
      HKF=HK*FREQ0(IL)
      DO 110 ID=1,ND
         T=TEMP(ID)
         DOP=DP0*SQRT(DP1*T+VTURB(ID))
         X=EXP(HKF/T)
         ABCENT(ILNLT,ID)=EGF*EXP(-EXCL0(IL)/T)*RRR(ID,ION,IAT)/
     *                    DOP*(1.-1./X)
         AB=ABSTD(ID)+ABCENT(ILNLT,ID)*1.77245
         IF(ID.EQ.1) THEN
            ABM=AB/DENS(1)
            TAU=0.5*DM(1)*ABM
          ELSE
            AB0=AB/DENS(ID)
            TAU=TAU+0.5*(DM(ID)-DM(ID-1))*(AB0+ABM)
            ABM=AB0
         END IF
C
C        approximate epsilon after Kastner
C
         E=EPS(T,ELEC(ID),ALMIL,ION,IUN)
         XK2=XK2DOP(TAU)
         SLIN(ILNLT,ID)=SQRT(E/(E+(1.-E)*XK2))*BNU/(X-1.)
  110 CONTINUE
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C

      SUBROUTINE NLTSET(MODE,IL,IAT,ION,EXCL,EXCU,IEVEN,INNLT0)
C     =========================================================
C
C     NLTE option -  automatic assignement of level indices
C
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: MODE,IL,IAT,ION
      real*8, intent(in   ) :: EXCL,EXCU
      integer,intent(inout) :: IEVEN,INNLT0
!       real*8, intent(inout) :: 
!       integer,intent(  out) :: 
!       real*8, intent(  out) :: 
!       integer :: 
!       real*8  :: 
        
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/LINDAT.FOR'
      PARAMETER (MNION = 20,
     *           MNLEV = 50,
     *           ECONST=  5.03411142E15)
      PARAMETER (INLION = 14,
     *           INLLEV = 13)
      COMMON/NL2PAR/ELIMEV(MNION,MNLEV),ELIMOD(MNION,MNLEV),
     *              ELIML(MNION,MNLEV),
     *              ENREV(MNION,MNLEV),ENROD(MNION,MNLEV),
     *              INDEV(MNION,MNLEV),INDOD(MNION,MNLEV),
     *              INDLV(MNION,MNLEV),
     *              NEVEN(MNION),NODD(MNION),NODD0,NLEVS(MNION),
     *              IATN(MNION),IONN(MNION),NNION
C
C     +++++++++++++++++++++++++++
C     MODE = 0  -  initialization
C     +++++++++++++++++++++++++++
C
      IF(MODE.EQ.0) THEN
         NNION=0
         READ(INLION,*,END=55,ERR=55) NNION
         IF(NNION.LE.0) GO TO 55
         DO 50 I=1,NNION
            READ(INLION,*) IATN(I),IONN(I)
            READ(INLLEV,*) NEVEN(I)
            IF(NEVEN(I).GT.0) THEN
               DO 10 J=1,NEVEN(I)
                  READ(INLLEV,*) II,ELIMEV(I,J),ENREV(I,J)
   10          CONTINUE
            END IF
            READ(INLLEV,*) NODD(I)
            NODD0=NODD(I)
            IF(NODD(I).GT.0) THEN
               DO 20 J=1,NODD(I)
                  READ(INLLEV,*) II,ELIMOD(I,J),ENROD(I,J)
   20          CONTINUE
             ELSE
               NODD(I)=NEVEN(I)
               DO 30 J=1,NODD(I)
                  ELIMOD(I,J)=ELIMEV(I,J)
   30          CONTINUE
            END IF
            INDION=0
            DO 40 IONEX=1,NION
               N0I=NFIRST(IONEX)
               IA=NUMAT(IATM(N0I))
               IF(IA.EQ.IATN(I).AND.IZ(IONEX)-1.EQ.IONN(I)) INDION=IONEX
   40       CONTINUE
            IF(INDION.LE.0) THEN
               CALL STOP(' INCONSISTENCY IN UNIT 14 INPUT - NLTE')
            END IF
            NOFF=NFIRST(INDION)-1
c
            INT=0
            INE=1
            INO=1
   22       INT=INT+1
            IF(INT.GT.NEVEN(I)+NODD(I)) GO TO 24
            IF(ENREV(I,INE).LT.ENROD(I,INO)) THEN
                INDEV(I,INE)=INT+NOFF
                INE=INE+1
                IF(INE.GT.NEVEN(I)) GO TO 28
             ELSE
                INDOD(I,INO)=INT+NOFF
                INO=INO+1
                IF(INO.GT.NODD(I)) GO TO 26
            END IF
            GO TO 22
   24       CONTINUE
            GO TO 50
C
   26       IF(INE.LE.NEVEN(I)) THEN
               DO 27 INE1=INE,NEVEN(I)
                  INT=INT+1
                  INDEV(I,INE1)=INT+NOFF
   27          CONTINUE
            END IF
            GO TO 50
   28       IF(INO.LE.NODD(I)) THEN
               DO 29 INO1=INO,NODD(I)
                  INT=INT+1
                  INDOD(I,INO1)=INT+NOFF
   29          CONTINUE
            END IF
C
   50    CONTINUE
   55    CONTINUE
C
         INDION=NNION
         DO 90 IONEX=1,NION
            N0I=NFIRST(IONEX)
            IA=NUMAT(IATM(N0I))
            IONM1=IZ(IONEX)-1
            IF(IA.EQ.1.OR.IA.EQ.2) GO TO 90
            DO 60 I=1,NNION
               IF(IA.EQ.IATN(I).AND.IONM1.EQ.IONN(I)) GO TO 90
   60       CONTINUE
            IF(NFIRST(IONEX).EQ.NLAST(IONEX)) GO TO 90
            INDION=INDION+1
            EION=ENION(NFIRST(IONEX))
            NLEVS(INDION)=NLAST(IONEX)-NFIRST(IONEX)+1
            NEVEN(INDION)=0
            IATN(INDION)=IA
            IONN(INDION)=IONM1
            DELE=0.
            DO 70 II=NFIRST(IONEX),NLAST(IONEX)
               I=II-NFIRST(IONEX)+1
               E=(EION-ENION(II))*ECONST
               IF(II.LT.NLAST(IONEX)) THEN
                  E1=(EION-ENION(II+1))*ECONST
                  DELE=0.5*(E1-E)
                  ELIML(INDION,I)=E+DELE
                ELSE
                  IF(INLTE.EQ.1) THEN
CCC-------------------------------------------------------------
CCC***  MARGIT HABERREITER
CMH   DELE ONLY ONE TENTH OF DIFFERENCE OF LOWER LEVELS
CCC-------------------------------------------------------------
				   DELE=0.1*(ENION(II-1)-ENION(II))*ECONST
CCC-------------------------------------------------------------
CMH   originally: DELE NOT CALCULATED
CCC-------------------------------------------------------------
                     ELIML(INDION,I)=E+DELE
                   ELSE
                     ELIML(INDION,I)=EION*ECONST
                  END IF
               END IF
               INDLV(INDION,I)=II
   70       CONTINUE
   90    CONTINUE
         NNION=INDION
  100    RETURN
      END IF
C
C
C     ++++++++++++++++++++++++++++++++++++++++++
C     MODE > 0  -  level indices for the line IL
C     ++++++++++++++++++++++++++++++++++++++++++
C
	IF(NNION.LE.0) RETURN
	INION=0
	IONM1=ION-1
	DO 102 I=1,NNION
		IF(IAT.EQ.IATN(I).AND.IONM1.EQ.IONN(I)) INION=I
  102	CONTINUE

	IF(INION.LE.0) RETURN
      IF(NEVEN(INION).EQ.0) IEVEN=2
C
	IF(IEVEN.EQ.1) THEN
		IND=0
		DO 110 J=1,NEVEN(INION)
			IF(EXCL.LE.ELIMEV(INION,J)) THEN
				IND=J
				GO TO 120
			END IF
  110		CONTINUE
	ILWN=0

CMH   OTHERWISE STEP INTO DIFFERENT ELSE IF STATEMENT
CMH	GO TO 200
      GO TO 145
  120	CONTINUE
	ILWN=INDEV(INION,IND)
C
	IND=0
	DO 130 J=1,NODD(INION)
		IF(EXCU.LE.ELIMOD(INION,J)) THEN
			IND=J
			GO TO 140
		END IF
  130	CONTINUE
	IUN=0
CMH   OTHERWISE STEP INTO DIFFERENT ELSE IF STATEMENT
CMH	GO TO 200
	GO TO 145
  140	CONTINUE
	IUN=INDOD(INION,IND)
  145	CONTINUE
	ELSE IF(IEVEN.EQ.0) THEN
         IND=0
         DO 150 J=1,NODD(INION)
            IF(EXCL.LE.ELIMOD(INION,J)) THEN
               IND=J
               GO TO 160
            END IF
  150    CONTINUE
         ILWN=0
         GO TO 200
  160	   CONTINUE
         ILWN=INDOD(INION,IND)
C
         IND=0
         DO 170 J=1,NEVEN(INION)
            IF(EXCU.LE.ELIMEV(INION,J)) THEN
               IND=J
               GO TO 180
            END IF
  170		CONTINUE
	IUN=0
	GO TO 200
  180	CONTINUE
	IUN=INDEV(INION,IND)
  200	CONTINUE
c
c        transition between levels without a distinction in parity
c
      ELSE
         IND=0
         DO 210 J=1,NLEVS(INION)
            IF(EXCL.LE.ELIML(INION,J)) THEN
               IND=J
               GO TO 220
            END IF
  210    CONTINUE
         ILWN=0
         IUN=0
         GO TO 300
  220		CONTINUE
		ILWN=INDLV(INION,IND)
C
		IND=0
		DO 230 J=1,NLEVS(INION)
            IF(EXCU.LE.ELIML(INION,J)) THEN
               IND=J
               GO TO 240
            END IF
  230    CONTINUE
         IUN=0
         GO TO 300
  240    CONTINUE
         IUN=INDLV(INION,IND)
  300    CONTINUE
      END IF
      IF(INLTE.NE.3) THEN
         INNLT0=INNLT0+1
         INDNLT(IL)=INNLT0
        ELSE
         INDNLT(IL)=-1
      END IF
      ILOWN(IL)=ILWN
      IUPN(IL)=IUN
      RETURN
      END SUBROUTINE
C
C
C     ****************************************************************
C
C


      SUBROUTINE OPADD(MODE,IJ,ID,FR,ABAD,EMAD,SCAD)
C     ====================================================
C
C     CALLED BY CROSET,
C	Additional opacities
C     This is basically user-supplied procedure; here are some more
C     important non-standard opacity sources, namely
C     Rayleigh scattering, H- opacity, H2+ opacity, and additional
C     opacity of He I and He II.
C     Inclusion of these opacities is contolled by switches transmitted
C     by COMMON/OPCPAR - see description in START.
C
C     Input parameters:
C     MODE  - controls the nature and the amount of calculations
C           = -1 - (OPADD called from START) evaluation of relevant
C                  depth-dependent quantities (usually photoionization
C                  cross-sections, but also possibly other), which are
C                  stored in array CROS
C           = 0  - evaluation of an additional opacity, emissivity, and
C                  scattering - for procedure OPAC0
C     IJ    - frequency index
C     IT    - number of levels
C     ID    - depth index
C     FR    - frequency
C     CROSS - array of photionization cross-section or other auxiliary
C             parameters for evaluating an additional opacity
C             Note: CROSS is output for MODE=-1, and input for other MODE
C
C     Output:
C
C     ABAD  - absorption coefficient (at frequency point IJ and depth ID)
C     EMAD  - emission coefficient (at frequency point IJ and depth ID)
C     SCAD  - scattering coefficient (at frequency point IJ and depth ID)
C	NDIM  - NUMBER OF LEVELS
C**********************************************************
CMH      INCLUDE 'PARAMS.FOR'
CMH      INCLUDE 'MODELP.FOR'
      use constants
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: MODE,IJ,ID
      real*8, intent(in   ) :: FR
      real*8, intent(  out) :: ABAD,EMAD,SCAD
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      COMMON/PHOPAR/CROSS(MCROSS,MFCONT)
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,
     $ IOPFE1
      PARAMETER (FRRAY =  2.463E15)
      REAL*8 :: AB0,AB1,X1
C	real*8::SIGINT

C	print *,'opadd.....'
      AB1=0.
      AB0=0.
      ABAD=0.
      EMAD=0.
      SCAD=0.
C
      if(iath.gt.0) then
      N0HN=NFIRST(IELH)
      NKH=NKA(IATH)
C
      IF(MODE.GE.0) THEN
         T=TEMP(ID)
         ANE=ELEC(ID)
         HKT=HK/T
         T32=1./T/SQRT(T)
      END IF
C
CMH   IT = NLEVE0?? WARUM?
      IT=NLEVE0
C
C   -----------------------
C   HI  Rayleigh scattering
C   -----------------------
C
      IF(IRSCT.NE.0.AND.IOPHLI.EQ.0) THEN
      IT=IT+1
      IF(MODE.LT.0) THEN
         X=(CLIGHT_SI*1d10/MIN(FR,FRRAY))**2
         CROSS(IT,IJ)=(5.799E-13+(1.422E-6+2.784/X)/X)/X/X
       ELSE
         ABAD=POPUL(N0HN,ID)*CROSS(IT,IJ)
         SCAD=ABAD
      END IF
      END IF
      IF(IOPHMI.NE.0) THEN
C
C   ----------------------------
C   H-  bound-free and free-free
C   ----------------------------
C     Note: IOPHMI must not by taken non-zero if H- is considered
C           explicitly, because H- opacity would be taken twice
C   ----------------------------------------------------------------
CMH   BUT IOPHMI IS NECESSARY FOR THE CALCULATION IN IOP
CMH	MH: H- is considered explicitely
CMH   therefore IOPHMI MUST BE 0!!!
C-------------------------------------------------------------------
      IT=IT+1
      IF(MODE.LT.0) THEN
         CROSS(IT,IJ)=SBFHMI(FR)
       ELSE
          XHM=8762.9/T
          SB=1.0353E-16*T32*EXP(XHM)*POPUL(N0HN,ID)*ANE*CROSS(IT,IJ)
          SF=SFFHMI(POPUL(N0HN,ID),FR,T)*ANE
          AB0=SB+SF
      END IF
      END IF

C
      IF(IOPH2P.GT.0) THEN
C
C   -----------------------------
C   H2+  bound-free and free-free
C   -----------------------------
C
      IT=IT+1
      IF(MODE.LT.0) THEN
         X=FR*1.E-15
         CROSS(IT,IJ)=(-7.342E-3+(-2.409+(1.028+(-4.23E-1+
     *        (1.224E-1-1.351E-2*X)*X)*X)*X)*X)*1.602E-12/BOLK
         IT=IT+1
         X=LOG(FR)
         CROSS(IT,IJ)=-3.0233E3+(3.7797E2+(-1.82496E1+(3.9207E-1-
     *                 3.1672E-3*X)*X)*X)*X
       ELSE
         X2=-CROSS(IT,IJ)/T+CROSS(IT+1,IJ)
         IT=IT+1
         SB=0.
         IF(X2.GT.-150.) SB=POPUL(N0HN,ID)*POPUL(NKH,ID)*EXP(X2)
         AB1=SB
      END IF
      END IF
      end if
C   ----------------------
C   neutral helium opacity
C   ----------------------
C     approximate, hydrogenic, opacity of neutral helium
C     given as a sum of bound-free transitions from averaged
C     levels with principal quantum numbers between that next
C     to the highest level considered explicitly and IOPHE1
C
      IF(IOPHE1.GT.0) THEN
      NHE2=NNEXT(IELHE1)
      I1=NQUANT(NLAST(IELHE1))+1
      IF(MODE.LT.0) THEN
         SG0=2.815E-16/(FR*1.E-15)**3
         DO 20 I=I1,IOPHE1
            IT=IT+1
            FRI=3.29E15/I/I
            SG=0.
            IF(FR.GE.FRI) SG=SG0*GAUNT(I,FR)/I**5
            CROSS(IT,IJ)=SG
   20    CONTINUE
       ELSE
         SB=4.1412E-16*T32
         DO 30 I=I1,IOPHE1
            IT=IT+1
            FRI=3.29E15/I/I
            IF(FR.LT.FRI) GO TO 30
            XI=HKT*FRI
            SG=SB*EXP(XI)*I*I*POPUL(NHE2,ID)*ANE*CROSS(IT,IJ)
            AB0=AB0+SG
   30    CONTINUE
      END IF
      END IF
C
C   --------------
C   ionized helium (analogously as for neutral helium)
C   --------------
C
      IF(IOPHE2.GT.0) THEN
      NHE3=NNEXT(IELHE2)
      I1=NQUANT(NLAST(IELHE2))+1
      IF(MODE.LT.0) THEN
         SG0=4.504E-15/(FR*1.E-15)**3
         DO 40 I=I1,IOPHE2
            IT=IT+1
            FRI=1.316E16/I/I
            SG=0.
            IF(FR.GE.FRI) SG=SG0*GAUNT(I,FR/4.E0)/I**5
            CROSS(IT,IJ)=SG
   40    CONTINUE
       ELSE
         SB=4.1412E-16*T32
         DO 50 I=I1,IOPHE2
            IT=IT+1
            FRI=1.316E16/I/I
            IF(FR.LT.FRI) GO TO 50
            XI=HKT*FRI
            SG=SB*EXP(XI)*I*I*POPUL(NHE3,ID)*ANE*CROSS(IT,IJ)
            AB0=AB0+SG
   50    CONTINUE
      END IF
      END IF
C
C   --------------
C     Finally, actual absorption and emission coefficients
C	"einstein coefficients" mihalas p. 102
C   --------------
      IF(MODE.LT.0) RETURN
      X=EXP(-HKT*FR)
      X1=1.-X
c correction for stimmulated emission (1.-X)
      BNX=BN*(FR*1.E-15)**3*X
      ABAD=ABAD+X1*(AB0+AB1)
      EMAD=EMAD+BNX*(AB0+AB1)
      RETURN
      END SUBROUTINE

C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE OPDATA
C     =================
C
C     Procedure reads photo-ionization cross sections fit coefficients
C     based on Opacity-Project (OP) data from file RBF.DAT
C     Data, as stored, requires linear interpolation.
C
C     Meaning of global variables:
C        NTOTOP    = total number of levels in Opacity Project data
C        IDLVOP() = level identifyer of current level
C        NOP()     = number of fit points for current level
C        XOP(,)    = x     = alog10(nu/nu0)       of fit point
C        SOP(,)    = sigma = alog10(sigma/10^-18) of fit point
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      PARAMETER    (MMAXOP = 200,! maximum number of levels in OP data
     +              MOP    =  15 )! maximum number of fit points per level
      CHARACTER*10  IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
      LOGICAL*2     LOPREA
      COMMON /TOPB/ SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
     +              XOP(MOP,MMAXOP) ,! x = alog10(nu/nu0) of fit point
     +              NOP(MMAXOP)     ,! number of fit points for current level
     +              NTOTOP          ,! total number of levels in OP data
     +              IDLVOP         ,! level identifyer Opacity-Project data
     +              LOPREA   ! .T. OP data read in; .F. OP data not yer read in
      CHARACTER*4 IONID
C
      OPEN (UNIT=40,FILE='RBF.DAT',STATUS='OLD',READONLY)
C     Skip header
      DO IREAD = 1, 21
         READ (40,*)
      END DO
      IOP = 0
C         = initialize sequential level index op Opacity Project data
C     Read number of elements in file
      READ (40,*) NEOP
      DO IEOP = 1, NEOP
C        Skip element name header
         DO IREAD = 1, 3
            READ (40,*)
         END DO
C        Read number of ionization stages of current element in  file
         READ (40,*) NIOP
         DO IIOP = 1, NIOP
C           Read ion identifyer, atomic & electron number, # of levels
C           for current ion
            READ (40,*) IONID, IATOM_OP, IELEC_OP, NLEVEL_OP
            DO ILOP = 1, NLEVEL_OP
C              Increase sequential level index of Opacity Project data
               IOP = IOP+1
C              Read level identifyer and number of sigma fit points
               READ (40,*) IDLVOP(IOP), NOP(IOP)
C              Read normalized log10 frequency and log10 cross section values
               DO IS = 1, NOP(IOP)
                  READ (40,*) INDEX, XOP(IS,IOP), SOP(IS,IOP)
               END DO
            END DO
         END DO
      END DO
      NTOTOP  = IOP
C             = total number of levels in Opacity Project data
      LOPREA  = .TRUE.
C             = set flag as data has been read in
      CLOSE (UNIT=40)
C
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
C
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION profil(IL,IAT,ID)
C     =================================
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
c      INCLUDE 'SYNTHP.FOR'
c      INCLUDE 'LINDAT.FOR'
      implicit none
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
	INCLUDE '../inc/SYNTHP.FOR'
	INCLUDE '../inc/LINDAT.FOR'
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
      DIMENSION WGR(4)
      real*8 :: wgr,pi4
      real*8 :: AGAM,profil
      integer,intent(in) :: IL,IAT,ID
      integer :: IPRF,I,ION
      real*8 :: T,ANE,ANP,GAM,FR,VDWC,DOPA1
      PARAMETER (PI4=7.95774715E-2)
C
      IPRF=IPRF0(IL)
      T=TEMP(ID)
      ANE=ELEC(ID)
C
C     radiative broadening (classical)
C
      AGAM=GAMR0(IL)
C
C     Stark broadening - standard (given in the line list or classical)
C
      IF(IPRF.EQ.0) THEN
         AGAM=AGAM+GS0(IL)*ANE
C
C     Stark broadening - special expressions for He I
C
       ELSE IF(IPRF.GT.0) THEN
         ANP=POPUL(NKH,ID)
         CALL GAMHE(IPRF,T,ANE,ANP,ID,GAM)
         AGAM=AGAM+GAM
C
C     Stark broadening - Griem
C
       ELSE
         DO 10 I=1,4
   10       WGR(I)=WGR0(I,IGRIEM(IL))
         FR=FREQ0(IL)
         ION=MOD(INDAT(IL),100)
         CALL GRIEM(ID,ANE,ION,FR,WGR,GAM)
         AGAM=AGAM+GAM
      END IF
C
C     Van Der Waals broadening
C
      AGAM=AGAM+GW0(IL)*VDWC(ID)
C
C     final Voigt parameter a
C
      AGAM=AGAM*DOPA1(IAT,ID)*PI4
C
      profil=AGAM
      RETURN
      END FUNCTION
C
C
C     ****************************************************************
C

      FUNCTION REIMAN(IB,FR)
C     ======================
C
C     Read table of photon energies and photo-ionization cross-sections
C     from Reilman & Manson (1979, Ap. J. Suppl., 40, 815) for the species
C     indicated by IB
C
C     Compute linearly interpolated value of the cross-section
C     at the frequency FR.
C
C     (At the moment, only a few transitions are considered)
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IB
      real*8, intent(in   ) :: FR
      DIMENSION HEV(30),F0(30),SIG0(30,2),SIGS(30)
C
      DATA HEV /
     * 130.,160.,190.,210.,240.,270.,300.,330.,360.,390.,
     * 420.,450.,480.,510.,540.,570.,600.,630.,660.,690.,
     * 720.,750.,780.,810.,840.,870.,900.,930.,960.,990./
      DATA SIG0 /
     * 3*0.,     4.422E-1, 3.478E-1,
     * 2.794E-1, 2.286E-1, 1.899E-1, 1.598E-1, 1.360E-1,
     * 1.169E-1, 1.013E-1, 8.845E-2, 7.776E-2, 6.877E-2,
     * 6.114E-2, 5.463E-2, 4.904E-2, 4.419E-2, 3.998E-2,
     * 3.629E-2, 3.305E-2, 3.019E-2, 2.766E-2, 2.540E-2,
     * 2.339E-2, 2.158E-2, 1.996E-2, 1.850E-2, 1.718E-2,
     * 4*0.,     1.981E-1, 1.584E-1,
     * 1.290E-1, 1.066E-1, 8.932E-2, 7.567E-2, 6.475E-2,
     * 5.589E-2, 4.862E-2, 4.259E-2, 3.754E-2, 3.329E-2,
     * 2.966E-2, 2.656E-2, 2.388E-2, 2.157E-2, 1.954E-2,
     * 1.777E-2, 1.621E-2, 1.484E-2, 1.362E-2, 1.253E-2,
     * 1.155E-2, 1.067E-2, 9.888E-3, 9.179E-3/
C
      INDEX=-IB-300
      NUM=30
      DO 10 I=1,NUM
         F0(I)=HEV(I)*2.418573E14
         SIGS(I)=SIG0(I,INDEX)
   10 CONTINUE
C
      IL=1
      IR=NUM
      DO 50 I=1,NUM-1
        IF(FR.GE.F0(I).AND.FR.LE.F0(I+1)) THEN
          IL=I
          IR=I+1
          GO TO 60
        ENDIF
 50   CONTINUE
C
C     LINEAR INTERPOLATION:
C
 60   SIGM=(SIGS(IR)-SIGS(IL))*(FR-F0(IL))/(F0(IR)-F0(IL))
     *      + SIGS(IL)
C
C     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE:
C
       IF(FR.LE.F0(1)) SIGM=SIGS(1)
       IF(FR.GE.F0(NUM)) SIGM=SIGS(NUM)
C
C     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION:
C
c       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL)
C
      REIMAN=SIGM*1.E-18
      RETURN
      END FUNCTION
C
C
C     ****************************************************************
C
C

      FUNCTION SBFHE1(II,IB,FR)
C     =========================
C
C     Calculates photoionization cross sections of neutral helium
C     from states with n = 1, 2, 3, 4.
C
C     The levels are either non-averaged (l,s) states, or some
C     averaged levels.
C     The program allows only two standard possibilities of
C     constructing averaged levels:
C     i)  all states within given principal quantum number n (>1) are
C         lumped together
C     ii) all siglet states for given n, and all triplet states for
C         given n are lumped together separately (there are thus two
C         explicit levels for a given n)
C
C     The cross sections are calculated using appropriate averages
C     of the Opacity Project cross sections, calculated by procedure
C     HEPHOT
C
C     Input parameters:
C      II    - index of the lower level (in the numbering of explicit
C              levels)
C      IB    - photoionization switch IBF for the given transition
C            = 10  -  means that the given transition is from an
C                     averaged level
C            = 11  -  the given transition is from non-averaged
C                     singlet state
C            = 13  -  the given transition is from non-averaged
C                     triplet state
C      FR    - frequency
C
c      INCLUDE 'PARAMS.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: II,IB
      real*8, intent(in   ) :: FR
      
	INCLUDE '../inc/PARAMS.FOR'
C
      NI=NQUANT(II)
      IGI=INT(G(II)+0.01)
      IS=IB-10
C
C     ----------------------------------------------------------------
C     IB=11 or 13  - photoionization from an non-averaged (l,s) level
C     ----------------------------------------------------------------
C
      IF(IS.EQ.1.OR.IS.EQ.3) THEN
         IL=(IGI/IS-1)/2
         SBFHE1=HEPHOT(IS,IL,NI,FR)
      END IF
C
C     ----------------------------------------------------------------
C     IS=10 - photoionization from an averaged level
C     ----------------------------------------------------------------
C
      IF(IS.EQ.0) THEN
         IF(NI.EQ.2) THEN
C
C ********    photoionization from an averaged level with n=2
C
            IF(IGI.EQ.4) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,2,FR)+3.D0*HEPHOT(1,1,2,FR))/9.D0
            ELSE IF(IGI.EQ.12) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,2,FR)+3.D0*HEPHOT(3,1,2,FR))/9.D0
            ELSE IF(IGI.EQ.16) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,2,FR)+3.D0*(HEPHOT(1,1,2,FR)+
     *                HEPHOT(3,0,2,FR))+9.D0*HEPHOT(3,1,2,FR))/1.6D1
            ELSE
               GO TO 10
            END IF
C
C
C ********    photoionization from an averaged level with n=3
C
         ELSE IF(NI.EQ.3) THEN
            IF(IGI.EQ.9) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,3,FR)+3.D0*HEPHOT(1,1,3,FR)+
     *                5.D0*HEPHOT(1,2,3,FR))/9.D0
            ELSE IF(IGI.EQ.27) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,3,FR)+3.D0*HEPHOT(3,1,3,FR)+
     *                5.D0*HEPHOT(3,2,3,FR))/9.D0
            ELSE IF(IGI.EQ.36) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,3,FR)+3.D0*HEPHOT(1,1,3,FR)+
     *                5.D0*HEPHOT(1,2,3,FR)+
     *                3.D0*HEPHOT(3,0,3,FR)+9.D0*HEPHOT(3,1,3,FR)+
     *                15.D0*HEPHOT(3,2,3,FR))/3.6D0
            ELSE
               GO TO 10
            END IF
         ELSE IF(NI.EQ.4) THEN
C
C ********    photoionization from an averaged level with n=4
C
            IF(IGI.EQ.16) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+
     *                 5.D0*HEPHOT(1,2,4,FR)+
     *                 7.D0*HEPHOT(1,3,4,FR))/1.6D1
            ELSE IF(IGI.EQ.48) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,4,FR)+3.D0*HEPHOT(3,1,4,FR)+
     *                 5.D0*HEPHOT(3,2,4,FR)+
     *                 7.D0*HEPHOT(3,3,4,FR))/1.6D1
            ELSE IF(IGI.EQ.64) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+
     *                 5.D0*HEPHOT(1,2,4,FR)+
     *                 7.D0*HEPHOT(1,3,4,FR)+
     *                 3.D0*HEPHOT(3,0,4,FR)+
     *                 9.D0*HEPHOT(3,1,4,FR)+
     *                 15.D0*HEPHOT(3,2,4,FR)+
     *                 21.D0*HEPHOT(3,3,4,FR))/6.4D1
            ELSE
               GO TO 10
            END IF
         ELSE
            GO TO 10
         END IF
      END IF
      RETURN
   10 WRITE(6,601) NI,IGI,IS
  601 FORMAT(1H0/' INCONSISTENT INPUT TO PROCEDURE SBFHE1'/
     * ' QUANTUM NUMBER =',I3,'  STATISTICAL WEIGHT',I4,'  S=',I3)
      STOP
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
C

      FUNCTION SBFHMI(FR)
      use constants,only:CLIGHT_CGS

C     ===================
C
C     Bound-free cross-section for H- (negative hydrogen ion)
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: FR
      SBFHMI=0.
      FR0=1.8259E14
      IF(FR.LT.FR0) RETURN
      IF(FR.LT.2.111E14) GO TO 10
      X=CLIGHT_CGS*1d5/FR
      SBFHMI=(6.80133E-3+X*(1.78708E-1+X*(1.6479E-1+X*(-2.04842E-2+X*
     1        5.95244E-4))))*1.E-17
      RETURN
   10 X=CLIGHT_CGS*1d5*(1./FR0-1./FR)
      SBFHMI=(2.69818E-1+X*(2.2019E-1+X*(-4.11288E-2+X*2.73236E-3)))
     1       *X*1.E-17
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION SFFHMI(POPI,FR,T)
C     ==========================
C
C     Free-free cross section for H- (After Kurucz,1970,SAO 309, P.80)
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: POPI,FR,T
CMH      SFFHMI=(1.3727E-25+(4.3748E-10-2.5993E-7/T)/FR)*POPI/FR
	SFFHMI=(1.3727E-25+(4.3748E-10-2.5993E-7/T)/FR)*POPI/FR
      RETURN
      END FUNCTION
C
C
C     ****************************************************************
C

      FUNCTION SGHE12(FR)
C     ===================
C
C     Special formula for the photoionization cross-section from the
C     averaged <n=2> level of He I
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: FR
      
      DATA C1/3.E0/,C2/9.E0/,C3/1.6E1/,
     * A1/6.45105E-18/,A2/3.02E-19/,A3/9.9847E-18/,A4/1.1763673E-17/,
     * A5/3.63662E-19/,A6/-2.783E2/,A7/1.488E1/,A8/-2.311E-1/,
     * E1/3.5E0/,E2/3.6E0/,E3/1.91E0/,E4/2.9E0/,E5/3.3E0/
      X=FR*1.E-15
      XX=LOG(FR)
      SGHE12=(C1*(A1/X**E1+A2/X**E2)+A3/X**E3+C2*(A4/X**E4+A5/X**E5)+
     *       C1*EXP(A6+XX*(A7+XX*A8)))/C3
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C CHANGED BY MARGIT HABERREITER
cmh      FUNCTION SIGK(FR,ITR,MODE)
      FUNCTION SIGK(FR,ITR,MODE,WAVARR,SIGARR,NDIM,NFDIM)
C     ==========================
C     CALLED BY CROSET (AND OTHERS?)
C     driver for evaluating the photoionization cross-sections
C
C     Input: FR  -  frequency
C            ITR -  index of the transition ? mh: rather index of level
C            MODE-  =0 - cross-section equal to zero longward of edge
C                   >0 - cross-section non-zero (extrapolated) longward
C                        of edge
C
c      INCLUDE 'PARAMS.FOR'
      use MOD_INTPLCS
      use constants,only: CLIGHT_CGS,CLIGHT_SI
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ITR,MODE,NDIM,NFDIM
      real*8, intent(in   ) :: FR,WAVARR,SIGARR
      INCLUDE '../inc/PARAMS.FOR'
      real*8,parameter :: SIH0=2.815D29
      DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)
      logical,save :: first=.true.
C	REAL*8 LIONVAC,LVAC
C      PARAMETER (CL    = 2.997925D10)
CMH  INPUT FOR HMINUS BF CALCULATION TAKEN FORM THE PHOTOCS_M - ROUTINE
CMH  H- FROM JOHN A&A193,189 (1988)
      real*8, parameter :: CN(6) =  (/
     &   152.519d0,49.534d0,-118.858d0,92.536d0,-34.194d0,4.982d0 /)
C
      PEACH(X,S,A,B)  =A*X**S*(B+X*(1.-B))*1.E-18
      HENRY(X,S,A,B,C)=A*X**S*(C+X*(B-2.*C+X*(1.+C-B)))*1.E-18
C
      SIGK=0.
     

      IF(IBF(ITR).EQ.0) RETURN
CMH     Energy of Groundstate != 0, now taken into acount in DATOM
CMH	ENION = EION/H/C
CMH	H*C=1.9857e-16
CMH	PLANCK CONSTANT: H=6.6256E-27 ERG*S
CMH     ENION - in luft ausrechnen!
CMH	LIONVAC: IONIZATION WAVELENGTH IN VACCUUM
CMH	LVAC		: VACCUUM WAVELENGTH
CMH	EVALUATE WAVELENGTH
CMH	ABOVE 2000 A REFORM TO AIR WAVELENGTH
c-vac	ENERGY = ENION(ITR)/1.9857e-16
c-vac	LIONVAC = 1.E8/ENERGY
c-vac	LVAC	= 2.997925D18/FR					! WAVELENGTH IN VACUUM, AND CM
c-vac	IF (LVAC .GT. 2000.) THEN
C***  LANG P.204
C          ANU2=(ELEVEL(I)-ELEVEL(J))*(ELEVEL(I)-ELEVEL(J))
c-vac		ANU2 = (1.d4/LIONVAC)*(1.d4/LIONVAC)			! anu in 10-6 m-1 (Edlen 1953)
c-vac          ANAIR=6432.8D-8+2949810./(146.D+8-ANU2)+25540./(41.D+8-ANU2)
C          WLENG=WLENG/(1.+ANAIR)
c-vac		FR0 = ENION(ITR)/6.6256E-27*(1.+ANAIR)	!WAVELENGTH IN AIR
c		PRINT*,ITR,ANAIR, LIONVAC,LVAC,FR,FR0
c-vac	ELSE
	    FR0=ENION(ITR)/6.6256E-27
c-vac	ENDIF
      IF(mode.eq.0  .and. FR.LT.FR0) THEN
	RETURN
	END IF

      
      IF(FR0.LE.0.) RETURN
C
C     IBF(ITR) is the switch controlling the mode of evaluation of the
C        cross-section:
C      = 0  hydrogenic cross-section, with Gaunt factor set to 1
C      = 1  hydrogenic cross-section with exact Gaunt factor
C      = 2  Peach-type expression (see function PEACH)
C      = 3  Henry-type expression (see function HENRY)
C      = 4  Butler new calculations
C      = 9  Opacity project fits (interpolations)
C      < 0  non-standard, user supplied expression (user should update
C           subroutine SPSIGK)
C
C      for H- : for any IBF > 0  - standard expression
C      for He I:
C      for IBF = 11 or = 13  -  Opacity Project cross section
C                Seaton-Ferney's cubic fits, Hummer's procedure (HEPHOT)
C           IBF = 11  means that the multiplicity S=1 (singlet)
C           IBF = 13  means that the multiplicity S=3 (triplet)
C       for IBF = 10  - cross section, based on Opacity Project, but
C                       appropriately averaged for an averaged level
CMH	II, ITR  = INDEX OF THE LEVEL
CMH   IELHM SET TO 1, FIRST ELEMENT IS HMINUS
c	PRINT *,'FUNCTION SIKG : IELHM SET TO ONE'
	IELHM = 1
      IB=IBF(ITR)
     

      II=ITR
      IQ=NQUANT(II)
      IE=IEL(II)
      IF(IB.LT.0) GO TO 60
      IF(IE.EQ.IELHM) GO TO 40
      IF(IE.EQ.IELHE1.AND.IB.GE.10.AND.IB.LE.13) GO TO 50
      IF(IE.EQ.IELHE1.AND.IB.EQ.21) GO TO 55
      CH=IZ(IE)*IZ(IE)
      IQ5=IQ*IQ*IQ*IQ*IQ
      FREL=FR0/FR
      IF(IB.EQ.0) THEN
C
C        hydrogenic expression (for IBF = 0)
C
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5
C
C        exact hydrogenic - with Gaunt factor (for IBF=1)
C
       ELSE IF(IB.EQ.1) THEN
 !        print*, 'BFtest', iq, GAUNT(IQ,FR), GAUNT_W(IQ,FR) 
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5
         SIGK=SIGK*GAUNT(IQ,FR/CH)
!         print*, 'BFtest', iq, sigk
       ELSE IF(IB.EQ.2) THEN
C
C        Peach-type formula (for IBF=2)
C
         SIGK=PEACH(FREL,S0BF(II),ALFBF(II),BETBF(II))
       ELSE IF(IB.EQ.3) THEN
C
C        Henry-type formula (for IBF=3)
C
         SIGK=HENRY(FREL,S0BF(II),ALFBF(II),BETBF(II),GAMBF(II))
C
C        Butler expression
C
       ELSE IF(IB.EQ.4) THEN
         XL=LOG(FREL)
         SL=S0BF(II)+XL*(ALFBF(II)+XL*BETBF(II))
         SIGK=EXP(SL)
C
C     selected Opacity Project data (for IBF=9)
C
       ELSE IF(IB.EQ.9) THEN
C*****CHANGES BY MARGIT HABERREITER
CMH   INTERPOLATE CROSS SECTIONS WITH INTPLCS
	    CALL intplcs(SIGK,FR,WAVARR,SIGARR,ITR,NDIM,NFDIM)
c	    CALL intplcs(SIGK,FR,FR0,ITR)
C         SIGK=TOPBAS(FR,FR0,TYPLEV(II))
C	PRINT *,'SIGK AFTER INTPLCS: SIGK=', SIGK
      END IF
      RETURN
C
C     special expression for H-
C
C  40	   XLAM: WAVELENGHT IN MYCRONS
CMH   CL: SPEED OF LIGHT
   40	XLAM=1.d4*CLIGHT_CGS/FR
	if(first) then
        PRINT *,'sigk ATTENTION: HMINUS EION INSETED BY HAND!!!'
        first=.false.
      endif
	WAVENUM = FR/CLIGHT_CGS
CMH   EION, EDGE, SIGMATH TAKEN FROM COOP_M
	EION = 6090.50000000000
      EDGE=EION
	SIGMATH = 1.e-18
c	   X=EDGE/WAVENUM

      X3=XLAM*XLAM*XLAM
CMH	avoid wavenumber - edge .lt. 0
	if (WAVENUM .lt. EDGE) then
		print *,'warning synsub: wavenumber lt edge'
		print *,WAVENUM ,edge
		SIGK = 0.
	else
	    X=SQRT((WAVENUM-EDGE)*1.d-4)
		if (wavenum.gt.50000.) then
				FLAMDA=127.84486d0*X
			ELSE
				FLAMDA = 0.d0
			DO I=6,1,-1
				FLAMDA=(FLAMDA+CN(I))*X
			ENDDO
		ENDIF
      SIGK=FLAMDA*X*X*X3*SIGMATH
	endif
CMH ORIGINALLY
CMH      X=SQRT((WAVENUM-EDGE)*1.d-4)
CMH      if (wavenum.gt.50000.) then
CMH		FLAMDA=127.84486d0*X
CMH      ELSE
CMH          FLAMDA = 0.d0
CMH          DO I=6,1,-1
CMH			FLAMDA=(FLAMDA+CN(I))*X
CMH          ENDDO
CMH      ENDIF
CMH      SIGK=FLAMDA*X*X*X3*SIGMATH
C***  ORIGINALLY
c   40 SIGK=SBFHMI(FR)
      RETURN
C
C     He I cross-sections
C
   50 SIGK=SBFHE1(II,IB,FR)
      RETURN
   55 X=LOG(CLIGHT_SI*1e10/FR)
      SIGK=EXP(-58.229+X*(4.3965-X*0.22134))/G(II)
      RETURN
C
C     non-standard, user supplied form of cross-section (for IBF < 0)
C
   60 CALL SPSIGK(IB,FR,SIGSP)
      SIGK=SIGSP
      RETURN
      END FUNCTION
C ********************************************************************
C

      SUBROUTINE SPSIGK(IB,FR,SIGSP)
C     ==================================
C
C     Non-standard evaluation of the photoionization cross-sections
C     Basically user-suppled procedure; here are some examples
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IB
      real*8, intent(in   ) :: FR
      real*8, intent(  out) :: SIGSP
      SIGSP=0.
C
C     Special formula for the He I ground state
C
      IF(IB.EQ.-201) SIGSP=7.3E-18*EXP(1.373-2.311E-16*FR)
C
C     Special formula for the averaged <n=2> level of He I
C
      IF(IB.EQ.-202) SIGSP=SGHE12(FR)
C
C     Carbon ground configuration levels 2p2 1D and 1S
C
      IF(IB.EQ.-602.OR.IB.EQ.-603) THEN
         CALL CARBON(IB,FR,SG)
         SIGSP=SG
      END IF
C
C     Hidalgo (Ap.J. 153, 981, 1968) photoionization data
C
      IF(IB.LE.-101.AND.IB.GE.-137) SIGSP=HIDALG(IB,FR)
C
C     Reilman and Manson (Ap.J. Suppl. 40, 815, 1979) photoionization data
C
      IF(IB.LE.-301.AND.IB.GE.-337) SIGSP=REIMAN(IB,FR)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE STARK0(I,J,IZZ,XKIJ,WL0,FIJ,FIJ0)
C
C     Auxiliary procedure for evaluating the approximate Stark profile
C     of hydrogen lines - sets up necessary frequency independent
C     parameters
C
C     Input:  I     - principal quantum number of the lower level
C             J     - principal quantum number of the upper level
C             IZZ   - ionic charge (IZZ=1 for hydrogen, etc.)
C     Output: XKIJ  - coefficients K(i,j) for the Hotzmark profile;
C                     exact up to j=6, asymptotic for higher j
C             WL0   - wavelength of the line i-j
C             FIJ   - Stark f-value for the line i-j
C             FIJ0  - f-value for the undisplaced component of the line
C
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: I,J,IZZ
      real*8, intent(  out) :: XKIJ,WL0,FIJ,FIJ0

        
      PARAMETER (RYD1=911.763811,RYD2=911.495745,CXKIJ=5.5E-5)
      PARAMETER (WI1=911.753578, WI2=227.837832)
      PARAMETER (UN=1.,TEN=10.,TWEN=20.,HUND=100.)
      DIMENSION FSTARK(10,4),XKIJT(5,4),FOSC0(10,4),FADD(5,5)
      DATA XKIJT/3.56E-4,5.23E-4,1.09E-3,1.49E-3,2.25E-3,.0125,.0177,
     * .028,.0348,.0493,.124,.171,.223,.261,.342,.683,.866,1.02,1.19,
     * 1.46/
      DATA FSTARK/  .1387,    .0791,   .02126,   .01394,   .00642,
     *           4.814E-3, 2.779E-3, 2.216E-3, 1.443E-3, 1.201E-3,
     *              .3921,    .1193,   .03766,   .02209,   .01139,
     *           8.036E-3, 5.007E-3,  3.85E-3, 2.658E-3, 2.151E-3,
     *              .6103,    .1506,   .04931,   .02768,   .01485,
     *             .01023, 6.588E-3, 4.996E-3, 3.524E-3, 2.838E-3,
     *              .8163,    .1788,   .05985,   .03189,   .01762,
     *             .01196, 7.825E-3, 5.882E-3, 4.233E-3, 3.375E-3/
      DATA FOSC0 / 0.27746,  0., 0.00773,  0., 0.00134, 0.,
     *             0.000404, 0., 0.000162, 0.,
     *             0.24869,  0., 0.00701,  0., 0.00131, 0.,
     *             0.000422, 0., 0.000177, 0.,
     *             0.23175,  0., 0.00653,  0., 0.00118, 0.,
     *             0.000392, 0., 0.000169, 0.,
     *             0.22148,  0.0005, 0.00563, 0.0004, 0.00108, 0.,
     *             0.000362, 0., 0.000159, 0./
      DATA FADD /  1.231, 0.2069, 7.448E-2, 3.645E-2, 2.104E-2,
     *             1.424, 0.2340, 8.315E-2, 4.038E-2, 2.320E-2,
     *             1.616, 0.2609, 9.163E-2, 4.416E-2, 2.525E-2,
     *             1.807, 0.2876, 1.000E-1, 4.787E-2, 2.724E-2,
     *             1.999, 0.3143, 1.083E-1, 5.152E-2, 2.918E-2/
C
      II=I*I
      JJ=J*J
      JMIN=J-I
      IF(JMIN.LE.5.and.i.le.4) THEN
         XKIJ=XKIJT(JMIN,I)
       ELSE
         XKIJ=CXKIJ*(II*JJ)*(II*JJ)/(JJ-II)
      END IF
      IF(I.LE.4) THEN
         IF(JMIN.LE.10) THEN
            FIJ=FSTARK(JMIN,I)
            FIJ0=FOSC0(JMIN,I)
          ELSE
            CFIJ=((TWEN*I+HUND)*J/(I+TEN)/(JJ-II))
            FIJ=FSTARK(10,I)*CFIJ*CFIJ*CFIJ
            FIJ0=0.
         END IF
       ELSE IF(I.LE.9) THEN
         IF(JMIN.LE.5) THEN
            FIJ=FADD(JMIN,I-4)
            FIJ0=0.
          ELSE
            CFIJ=((TEN*I+25.)*J/(I+5.)/(JJ-II))
            FIJ=FADD(5,I-4)*CFIJ*CFIJ*CFIJ
            FIJ0=0.
         END IF
       ELSE
         CFIJ=UN*J/(JJ-II)
         FIJ=1.96*I*CFIJ*CFIJ*CFIJ
         FIJ0=0.
      END IF
C      if(i.eq.1) then
C         wl0=ryd1/(un/ii-un/jj)
C       else
C         WL0=RYD2/(UN/II-UN/JJ)
C      end if
C     wavelength with an explicit correction to the air wavalength
C
      w0=wi1
      if(izz.eq.2) w0=wi2
      WL0=W0/(UN/II-UN/JJ)
      IF(WL0.GT.2000.) THEN
         ALM=1.E8/(WL0*WL0)
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0/(XN1*1.D-6+UN)
      END IF
C      fij=0.
C      fij0=0.

      RETURN
      END SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION STARKA(BETA,BETAD,A,DIV,FAC)
C
C     Approximate expressions for the hydrogen Stark profile
C
C     Input: BETA  - delta lambda in beta units,
C            BETAD - Doppler width in beta units
C            A     - auxiliary parameter
C                    A=1.5*LOG(BETAD)-1.671
C            DIV   - only for A > 1; division point between Doppler
C                    and asymptotic Stark wing, expressed in units
C                    of betad.
C                    DIV = solution of equation
C                    exp(-(beta/betad)**2)/betad/sqrt(pi)=
C                     = 1.5*FAC*beta**-5/2
C                    (ie. the point where Doppler profile is equal to
C                     the asymptotic Holtsmark)
C                    In order to save computer time, the division point
C                    DIV is calculated in advance by routine DIVSTR.
C            FAC   - factor by which the Holtsmark profile is to be
C                    multiplied to get total Stark Profile
C                    FAC should be taken to 2 for hydrogen, (and =1
C                    for He II)
C
c      INCLUDE 'IMPLIC.FOR'
      INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: BETA,BETAD,A,DIV,FAC
      real*8  :: STARKA
        
      PARAMETER (F0=-0.5758228,F1=0.4796232,F2=0.07209481,AL=1.26)
      PARAMETER (SD=0.5641895,SLO=-2.5,TRHA=1.5,BL1=1.14,BL2=11.4)
      PARAMETER (SAC=0.08)
      XD=BETA/BETAD
C
C     for a > 1 Doppler core + asymptotic Holtzmark wing with division
C               point DIV
C
      IF(A.GT.AL) THEN
         IF(XD.LE.DIV) THEN
            STARKA=SD*EXP(-XD*XD)/BETAD
          ELSE
            STARKA=TRHA*FAC*EXP(SLO*LOG(BETA))
         END IF
       ELSE
C
C     empirical formula for a < 1
C
         IF(BETA.LE.BL1) THEN
            STARKA=SAC
          ELSE IF(BETA.LT.BL2) THEN
            XL=LOG(BETA)
            FL=(F0*XL+F1)*XL
            STARKA=F2*EXP(FL)
          ELSE
            STARKA=TRHA*FAC*EXP(SLO*LOG(BETA))
         END IF
      END IF
      RETURN
      END FUNCTION
c
c
c
      subroutine stop(text)
c
c     stops the program and writes a text
c
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      character*(*) text
      write(6,10) text
   10 format(1x,a)
      stop
      end SUBROUTINE
C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE TINT
C
C     LOGARITHMIC INTERPOLATION COEFFICIENTS FOR INTERPOLATION OF
C     TEMP(ID) TO THE VALUES  5000,10000,20000,40000
C
c     INCLUDE 'PARAMS.FOR'
c     INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      DIMENSION TT(4)
      DATA TT /3.699, 4.000, 4.301, 4.602/
C
      DO 10 ID=1,ND
         T=LOG10(TEMP(ID))
         J=3
         IF(T.GT.TT(3)) J=4
         JT(ID)=J
         X=(TT(J)-TT(J-1))*(TT(J)-TT(J-2))*(TT(J-1)-TT(J-2))
         TI0(ID)=(T-TT(J-2))*(T-TT(J-1))*(TT(J-1)-TT(J-2))/X
         TI1(ID)=(T-TT(J-2))*(TT(J)-T)*(TT(J)-TT(J-2))/X
         TI2(ID)=(T-TT(J-1))*(T-TT(J))*(TT(J)-TT(J-1))/X
   10 CONTINUE
      RETURN
      END SUBROUTINE
C
C
C     ******************************************************************
C
C***  CHANGES BY MARGIT HABERREITER
CMH   FUNCTION TOPBAS(FREQ,FREQ0,TYPLV)
! 	FUNCTION TOPBAS(FREQ,FREQ0,TYPLV)
! C     ==================================
! C
! C     Procedure calculates the photo-ionisation cross section SIGMA in
! C     [cm^2] at frequency FREQ. FREQ0 is the threshold frequency from
! C     level I of ion KI. Threshold cross-sections will be of the order
! C     of the numerical value of 10^-18.
! C     Opacity-Project (OP) interpolation fit formula
! C
! c      INCLUDE 'IMPLIC.FOR'
! 	INCLUDE '../inc/IMPLIC.FOR'
!       real*8, intent(in   ) :: FREQ,FREQ0
!       character*10,intent(in) ::  TYPLV
!       PARAMETER    (E10=2.3025851)
!       PARAMETER    (MMAXOP = 200,! maximum number of levels in OP data
!      +              MOP    =  15 )! maximum number of fit points per level
!       CHARACTER*10  IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
!       LOGICAL*2     LOPREA
!       COMMON /TOPB/ SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
!      +              XOP(MOP,MMAXOP) ,! x = alog10(nu/nu0) of fit point
!      +              NOP(MMAXOP)     ,! number of fit points for current level
!      +              NTOTOP          ,! total number of levels in OP data
!      +              IDLVOP         ,! level identifyer Opacity-Project data
!      +              LOPREA   ! .T. OP data read in; .F. OP data not yer read in
!       DIMENSION XFIT(MOP) ,! local array containing x     for OP data
!      +          SFIT(MOP)  ! local array containing sigma for OP data
! 
! C     Read OP data if not yet done
! C
!       TOPBAS=0.D0
! C      IF (.NOT.LOPREA) CALL OPDATA
!       X = LOG10(FREQ/FREQ0)
!       DO IOP = 1,NTOTOP
!          IF (IDLVOP(IOP).EQ.TYPLV) THEN
! C           level has been detected in OP-data file
!             IF (NOP(IOP).LE.0) GO TO 20
!             DO IFIT = 1,NOP(IOP)
!                XFIT(IFIT) = XOP(IFIT,IOP)
!                SFIT(IFIT) = SOP(IFIT,IOP)
!             END DO
!             SIGM  = YLINTP (X,XFIT,SFIT,NOP(IOP),MOP)
!             SIGM  = 1.D-18*EXP(E10*SIGM)
!             TOPBAS=SIGM
!             GO TO 10
!          END IF
!       END DO
!    10 RETURN
! C     Level is not found ,or no data for this level, in RBF.DAT
!    20 WRITE (61,100) TYPLV
!   100 FORMAT ('SIGMA.......: OP DATA NOT AVAILABLE FOR LEVEL ',A10)
!       RETURN
!       END FUNCTION
c
c

C
C *******************************************************************
C


      function voigte(a,vs)
c     =====================
c
c  computes a voigt function  h = h(a,v)
c  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
c  traving (landolt-b\rnstein, p. 449).
c
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: a,vs
      real*8  :: voigte
      
      dimension ak(19),a1(5)
      data ak      /-1.12470432, -0.15516677,  3.28867591, -2.34357915,
     ,  0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,
     , -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,
     ,-0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278/
      data sqp/1.772453851/,sq2/1.414213562/
c
      v = abs(vs)
      u = a + v
      v2 = v*v
      if (a.eq.0.0) go to 140
      if (a.gt.0.2) go to 120
      if (v.ge.5.0) go to 121
c
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 1
c
  100 quo = 1.
      if (v.lt.2.4) go to 101
      quo = 1./(v2 - 1.5)
      m = 11
      go to 102
c
  101 m = 6
      if (v.lt.1.3) m = 1
  102 do 103 i=1,5
         a1(i) = ak(m)
         m = m + 1
  103 continue
      h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))
      if (k.gt.1) go to 110
c
c a le 0.2  and v lt 5.
c
      h = h1*a + ex*(1. + a*a*(1. - 2.*v2))
      voigte=h
      return
c
  110 pqs = 2./sqp
      h1p = h1 + pqs*ex
      h2p = pqs*h1p - 2.*v2*ex
      h3p = (pqs*(1. - ex*(1. - 2.*v2)) - 2.*v2*h1p)/3. + pqs*h2p
      h4p = (2.*v2*v2*ex - pqs*h1p)/3. + pqs*h3p
      psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
c
c 0.2 lt a le 1.4  and  a + v le 3.2
c
      h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
      voigte=h
      return
c
  120 if (a.gt.1.4.or.u.gt.3.2) go to 130
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 2
      go to 100
c
c a le 0.2  and  v ge 5.
c
  121 h = a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp)
      voigte=h
      return
c
  130 a2 = a*a
      u = sq2*(a2 + v2)
      u2 = 1./(u*u)
c
c a gt 1.4  or  a + v gt 3.2
c
      h = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) +
     ,        u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2))
      voigte=h
      return
c
c a eq 0.
c
  140 h=0.
      if(v2.lt.100.)h=exp(-v2)
      voigte=h
      return
      end function
C
C
C     ****************************************************************
C
C
      pure function wn(xn,a,ane,z)
c
c     evaluation of the occupation probablities for a hydrogenic ion
c     using eqs (4.26), and (4.39) of Hummer,Mihalas Ap.J. 331, 794, 1988.
c     approximate evaluation of Q(beta) - Hummer
c
c     Input: xn  - real number corresponding to quantum number n
c            a   - correlation parameter
c            ane - electron density
c            z   - ionic charge
c
c     Limitations: 
      INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: xn,a,ane,z
      real*8  :: wn
        
      parameter (p1=0.1402,p2=0.1285,p3=1.,p4=3.15,p5=4.,un=1.)
      parameter (tkn=3.01,ckn=5.33333333,cb=8.3e14)
      parameter (f23=-2./3.)
c
c     evaluation of k(n)  (4.24)
c
      if(xn.ge.3.01) then
         xkn=un
       else
         xkn=ckn*xn/(xn+un)/(xn+un)
      end if
c
c     evaluation of beta 
c        = K_n(1/16n^4)(4pi/3)^(-2/3)a_0^-2 C^(1/3)Z^3N_e^(-2/3)  (4.39)
c
      !  cb=1/16 * (4pi/3)^{-2/3} * a_0^{-2} * C^{1/3}
      !  8.3e14*Z^3*K_n / n^4        N_e^(-2/3)
      beta=cb*z*z*z*xkn/xn/xn/xn/xn*exp(f23*log(ane))
c
c     approximate expression for Q(beta)
c
      x=exp(p4*log(un+p3*a))
      c1=p1*(x+p5*z*a*a*a)
      c2=p2*x
      f=(c1*beta*beta*beta)/(un+c2*beta*sqrt(beta))
      wn=f/(un+f)
      return
      end function
C
C ********************************************************************
C ********************************************************************
C

      FUNCTION WTOT(T,ANE,ID,ILINE)
C     =============================
C
C     Evaluates the total (electron + ion) impact Stark width
C     for four HeI lines
C     After Griem (1974); and Barnard, Cooper, Smith (1974) JQSRT 14,
C     1025 for the 4471 line
C
C     Input: T     - temperature
C            ANE   - electron density
C            ID    - depth index
C            ILINE - index of the line ( = 1  for 4471,
C                                        = 2  for 4387,
C                                        = 3  for 4026,
C                                        = 4  for 4922)
C     Output: WTOT - Stark width in Angstroms

      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID,ILINE
      real*8, intent(in   ) :: T,ANE
      real*8  :: WTOT
        
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      DIMENSION ALPH0(4,4),W0(4,4),ALAM0(4)
      DATA ALPH0 / 0.107, 0.119, 0.134, 0.154,
     *             0.206, 0.235, 0.272, 0.317,
     *             0.172, 0.193, 0.218, 0.249,
     *             0.121, 0.136, 0.157, 0.184/
      DATA W0    / 1.460, 1.269, 1.079, 0.898,
     *             6.130, 5.150, 4.240, 3.450,
     *             4.040, 3.490, 2.960, 2.470,
     *             2.312, 1.963, 1.624, 1.315/
      DATA ALAM0 / 4471.50, 4387.93, 4026.20, 4921.93/
C
      I=JT(ID)
      ALPHA=(TI0(ID)*ALPH0(I,ILINE)+TI1(ID)*ALPH0(I-1,ILINE)+
     *      TI2(ID)*ALPH0(I-2,ILINE))*(ANE*1.E-13)**0.25
      WE=   (TI0(ID)*W0(I,ILINE)+TI1(ID)*W0(I-1,ILINE)+
     *      TI2(ID)*W0(I-2,ILINE))*ANE*1.E-16
      F0=1.884E19/ALAM0(ILINE)/ALAM0(ILINE)
      SIG=(4.32E-5*WE/SQRT(T)*F0/ANE**0.3333)**0.3333
      WTOT=WE*(1.+1.36/SIG*ALPHA**0.8889)
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION YINT(XL,YL,XL0)
C
C     Quadratic interpolation routine
C
C     Input:  XL - array of x
C             YL - array of f(x)
C             XL0 - the point x(0) to which one interpolates
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      real*8, intent(in   ) :: XL,YL,XL0
      DIMENSION XL(3),YL(3)
      A0=(XL(2)-XL(1))*(XL(3)-XL(2))*(XL(3)-XL(1))
      A1=(XL0-XL(2))*(XL0-XL(3))*(XL(3)-XL(2))
      A2=(XL0-XL(1))*(XL(3)-XL0)*(XL(3)-XL(1))
      A3=(XL0-XL(1))*(XL0-XL(2))*(XL(2)-XL(1))
      YINT=(YL(1)*A1+YL(2)*A2+YL(3)*A3)/A0
      RETURN
      END FUNCTION
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION YLINTP (XINT,X,Y,N,NTOT)
C     =================================
C
C     linear interpolation routine. Determines YINT = Y(XINT) from
C     grid Y(X) with N points and dimension NTOT.
C
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: N,NTOT
      real*8, intent(in   ) :: XINT,X,Y
      real*8  :: YLINTP
        
      DIMENSION X(NTOT),Y(NTOT)
C
C     bisection (see Numerical Recipes par 3.4 page 90)
      JL = 0
      JU = N+1
10    IF (JU-JL.GT.1) THEN
         JM = (JU+JL)/2
         IF ((X(N).GT.X(1)).EQV.(XINT.GT.X(JM))) THEN
            JL = JM
         ELSE
            JU = JM
         END IF
         GO TO 10
      END IF
      J = JL
      IF (J.EQ.N) J = J-1
      IF (J.EQ.0) J = J+1
CM42 version: enable extrapolation
cccc  IF (J.EQ.0.OR.J.EQ.N) CALL ERRPRI
cccc +                      ('YLININTERP: EXTRAPOLATION')
C
      RC         = (Y(J+1)-Y(J))/(X(J+1)-X(J))
      YLINTP = RC*(XINT-X(J))+Y(J)
C
      RETURN
      END FUNCTION
C
C
C     *******************************************************************
C
C

      SUBROUTINE PRETAB
C     =================
C
C     pretabulate expansion coefficients for the Voigt function
C     200 steps per doppler width - up to 10 Doppler widths
C
c      INCLUDE 'PARAMS.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
      PARAMETER (VSTEPS=200.,MVOI=2001)
      real*8 :: H0TAB,H1TAB,H2TAB
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      DIMENSION TABVI(81),TABH1(81)
      DATA TABVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,
     11.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,
     2 3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,
     3 5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,
     4 9.0,9.2,9.4,9.6,9.8,10.0,10.2,10.4,10.6,10.8,11.0,11.2,11.4,11.6,
     5 11.8,12.0/
      DATA TABH1/-1.12838,-1.10596,-1.04048,-.93703,-.80346,-.64945,
     1-.48552,-.32192,-.16772,-.03012,.08594,.17789,.24537,.28981,
     2.31394,.32130,.31573,.30094,.28027,.25648,.231726,.207528,.184882,
     3.164341,.146128,.130236,.116515,.104739,.094653,.086005,.078565,
     4 .072129,.066526,.061615,.057281,.053430,.049988,.046894,.044098,
     5 .041561,.039250,.035195,.031762,.028824,.026288,.024081,.022146,
     6 .020441,.018929,.017582,.016375,.015291,.014312,.013426,.012620,
     7 .0118860,.0112145,.0105990,.0100332,.0095119,.0090306,.0085852,
     8 .0081722,.0077885,.0074314,.0070985,.0067875,.0064967,.0062243,
     9 .0059688,.0057287,.0055030,.0052903,.0050898,.0049006,.0047217,
     T .0045526,.0043924,.0042405,.0040964,.0039595/
C
      N=MVOI
      DO 10 I=1,N
   10    H0TAB(I)=FLOAT(I-1)/VSTEPS
      CALL INTERP(TABVI,TABH1,H0TAB,H1TAB,81,N,2,0,0)
      DO 20 I=1,N
         VV=(FLOAT(I-1)/VSTEPS)**2
         H0TAB(I)=EXP(-VV)
         H2TAB(I)=H0TAB(I)-(VV+VV)*H0TAB(I)
   20 CONTINUE
      RETURN
      END SUBROUTINE
C
C
C     *******************************************************************
C
C

      FUNCTION VOIGTK(A,V)
C     ====================
C
C     Voigt function after Kurucz (in Computational Astrophysics)
C
c      INCLUDE 'PARAMS.FOR'
      implicit real*8(a-h,o-z)
      real*8, intent(in   ) :: A,V
      real*8  :: VOIGTK
	INCLUDE '../inc/PARAMS.FOR'
      PARAMETER (MVOI=2001)
      PARAMETER (ONE=1., THREE=3., TEN=10., FIFTN=15., TWOH=200.,
     *           C14142=1.4142, C11283=1.12838, C15=1.5,C32=3.2,
     *           C05642=0.5642,C79788=0.79788,C02=0.2,C14=1.4,
     *           C37613=0.37613,C23=2./3.,
     *           CV1=-.122727278,CV2=.532770573,CV3=-.96284325,
     *           CV4=.979895032)
      real*8 :: H0TAB,H1TAB,H2TAB
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      IV=V*TWOH+C15
      IF(A.LT.C02) THEN
         IF(V.LE.TEN) THEN
            VOIGTK=(H2TAB(IV)*A+H1TAB(IV))*A+H0TAB(IV)
          ELSE
            VOIGTK=C05642*A/(V*V)
         END IF
         RETURN
      END IF
      IF(A.GT.C14) GO TO 10
      IF(A+V.GT.C32) GO TO 10
      VV=V*V
      HH1=H1TAB(IV)+H0TAB(IV)*C11283
      HH2=H2TAB(IV)+HH1*C11283-H0TAB(IV)
      HH3=(ONE-H2TAB(IV))*C37613-HH1*C23*VV+HH2*C11283
      HH4=(THREE*HH3-HH1)*C37613+H0TAB(IV)*C23*VV*VV
      VOIGTK=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0TAB(IV))*
     *       (((CV1*A+CV2)*A+CV3)*A+CV4)
      RETURN
   10 AA=A*A
      VV=V*V
      U=(AA+VV)*C14142
      UU=U*U
      VOIGTK=((((AA-TEN*VV)*AA*THREE+FIFTN*VV*VV)/UU+THREE*VV-AA)/UU+
     *       1d0)*A*C79788/U
      RETURN
      END FUNCTION
C
C
C     ****************************************************************
C
C

      SUBROUTINE INTERP(X,Y,XX,YY,NX,NXX,NPOL,ILOGX,ILOGY)
C     ====================================================
C
C     General interpolation procedure of the (NPOL-1)-th order
C
C     for  ILOGX = 1  logarithmic interpolation in X
C     for  ILOGY = 1  logarithmic interpolation in Y
C
C     Input:
C      X    - array of original x-coordinates
C      Y    - array of corresponding functional values Y=y(X)
C      NX   - number of elements in arrays X or Y
C      XX   - array of new x-coordinates (to which is to be
C             interpolated
C      NXX  - number of elements in array XX
C     Output:
C      YY   - interpolated functional values YY=y(XX)
C
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: NX,NXX,NPOL,ILOGX,ILOGY
      real*8, intent(inout) :: X,Y,XX,YY
	INCLUDE '../inc/PARAMS.FOR'
      DIMENSION X(1),Y(1),XX(1),YY(1)
      EXP10(X0)=EXP(X0*2.30258509299405D0)
      IF(NPOL.LE.0.OR.NX.LE.0) GO TO 200
      IF(ILOGX.EQ.0) GO TO 30
      DO 10 I=1,NX
   10    X(I)=LOG10(X(I))
      DO 20 I=1,NXX
   20    XX(I)=LOG10(XX(I))
   30 IF(ILOGY.EQ.0) GO TO 50
      DO 40 I=1,NX
   40    Y(I)=LOG10(Y(I))
   50 NM=(NPOL+1)/2
      NM1=NM+1
      NUP=NX+NM1-NPOL
      DO 100 ID=1,NXX
         XXX=XX(ID)
         DO 60 I=NM1,NUP
            IF(XXX.LE.X(I)) GO TO 70
   60    CONTINUE
         I=NUP
   70    J=I-NM
         JJ=J+NPOL-1
         YYY=0.
         DO 90 K=J,JJ
            T=1.
            DO 80 M=J,JJ
               IF(K.EQ.M) GO TO 80
               T=T*(XXX-X(M))/(X(K)-X(M))
   80       CONTINUE
   90    YYY=Y(K)*T+YYY
         YY(ID)=YYY
  100 CONTINUE
      IF(ILOGX.EQ.0) GO TO 130
      DO 110 I=1,NX
  110    X(I)=EXP10(X(I))
      DO 120 I=1,NXX
  120    XX(I)=EXP10(XX(I))
  130 IF(ILOGY.EQ.0) RETURN
      DO 140 I=1,NX
  140    Y(I)=EXP10(Y(I))
      DO 150 I=1,NXX
  150    YY(I)=EXP10(YY(I))
      RETURN
  200    N=NX
         IF(NXX.GE.NX) N=NXX
      DO 210 I=1,N
         XX(I)=X(I)
  210    YY(I)=Y(I)
      RETURN
      END SUBROUTINE
      
      
      SUBROUTINE HYDINI(IHYDPR)
      use UTILS,only:assert
C
C     Initializes necessary arrays for evaluating hydrogen line profiles
C     from the Schoening and Butler or Lemke tables
C
      use MOD_HYDTAB, only: PRF,WLINE,ILIN0,WL,XT,XNE,NWLH,NTH,NEH
      use MOD_HYDTAB, only: MLINH,MHWL,MHE,MHT
      use PARAMS_ARRAY,only:NDDIM
      implicit none
      INCLUDE '../inc/PARAMS.FOR'
      !*** MODELP: TEMP,VTURB,ELEC
      INCLUDE '../inc/MODELP.FOR'
!       COMMON/VCSDAT/WL(MHWL,MLINH),XT(MHT,MLINH),
!      *              XNE(MHE,MLINH),PRF(MHWL,MHT,MHE),
!      *              NWLH(MLINH),NTH(MLINH),NEH(MLINH)
      real*8  ::    XK,FXK,BETAD,DBETA
      COMMON/AUXVCS/XK,FXK,BETAD,DBETA
      integer :: IHYDP(4)
      !COMMON/BALVCS/PRFBAL(8,MDEPTH,36),WLBAL(8,36),NWLH(8)
      integer :: NWLHYD(MLINH)
      real*8  :: WLHYD(MLINH,MHWL)
      real*8  :: PRFHYD(MLINH,MDEPTH,MHWL)

      COMMON/HYDVCS/PRFHYD,WLHYD,NWLHYD
      !*** end debugging
      CHARACTER*1 :: CHAR
      integer,save :: INIT = 0
      integer :: I, J,IZZ,IHYDPR,ILEMKE,NLINE,ILINE,NWL,NT,NE
      integer :: IE,IT,IWL,ID,IHYDP0,NTAB,ITAB,ILINEB,NLLY
      integer :: ILI,INE,ILNE
      real*8  :: TMIN,DLA,DLE,DLT,QLT
      real*8  :: T,ANE,ANEL,F00,DOP,PROF,TL,ALMIN,ANEMIN
      real*8  :: FIJ,WL0,FIJ0,XCLOG,XKLOG
C
      IF(INIT.EQ.0) THEN
        DO I=1,4
          DO J=I+1,22
            CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0)
            WLINE(I,J)=WL0
            ! OSCH(I,J)=FIJ+FIJ0
          END DO
        END DO
        INIT=1
      END IF
      ILIN0(:4,:22)=0
      !
      ! --------------------------------------------
      !     Schoening-Butler tables - for IHYDPR < 0
      ! --------------------------------------------
      !
      IF(IHYDPR.LT.0) THEN
      IHYDPR=-IHYDPR
      ILEMKE=0
      NLINE=12
      !
      ! OPEN(UNIT=IHYDPR,STATUS='OLD')
      !
      DO I=1,12
         READ(IHYDPR,500)
      END DO
      DO 100 ILINE=1,NLINE
        !
        !     read the tables, which have to be stored in file
        !     unit IHYDPR (which is the input parameter in the progarm)
        !
         READ(IHYDPR,501) I,J
         IF(ILINE.EQ.12) J=10
         WL0=WLINE(I,J)
         ILIN0(I,J)=ILINE
         READ(IHYDPR,*) CHAR,NWL,(WL(I,ILINE),I=1,NWL)
         READ(IHYDPR,*) CHAR,NT,(XT(I,ILINE),I=1,NT)
         READ(IHYDPR,*) CHAR,NE,(XNE(I,ILINE),I=1,NE)
         READ(IHYDPR,500)
         NWLH(ILINE)=NWL
         NWLHYD(ILINE)=NWL
         NTH(ILINE)=NT
         NEH(ILINE)=NE
C
         DO 10 I=1,NWL
            IF(WL(I,ILINE).LT.1.E-4) WL(I,ILINE)=1.E-4
            WLHYD(ILINE,I)=LOG10(WL(I,ILINE))
   10    ENDDO
C
         DO IE=1,NE
            DO IT=1,NT
               READ(IHYDPR,500)
               READ(IHYDPR,*) (PRF(IWL,IT,IE),IWL=1,NWL)
            ENDDO
         ENDDO
C
C        coefficient for the asymptotic profile is determined from
C        the input data
C
         XCLOG=PRF(NWL,1,1)+2.5*LOG10(WL(NWL,ILINE))+31.5304-
     *         XNE(1,ILINE)-2.*LOG10(WL0)
         XKLOG=0.6666667*(XCLOG-0.176)
         XK=EXP(XKLOG*2.3025851)
C
         DO ID=1,ND
C
C           temperature is modified in order to account for the
C           effect of turbulent velocity on the Doppler width
C
            T=TEMP(ID)+6.06E-9*VTURB(ID)
            ANE=ELEC(ID)
            TL=LOG10(T)
            ANEL=LOG10(ANE)
            F00=1.25E-9*ANE**0.666666667
            FXK=F00*XK
            DOP=1.E8/WL0*SQRT(1.65E8*T)
            DBETA=WL0*WL0/2.997925E18/FXK
            BETAD=DBETA*DOP
C
C       interpolation to the actual values of temperature and electron
C       density. The result is stored at array PRFHYD, having indices
C       ILINE (line number: 1 for L-alpha,..., 4 for H-delta, etc.);
C                           5 for H-alpha,..., 8 for H-delta, etc.)
C       ID - depth index
C       IWL - wavelength index 
C
            DO IWL=1,NWL
               CALL INTHYD(PROF,TL,ANEL,IWL,ILINE,ILEMKE)
               PRFHYD(ILINE,ID,IWL)=PROF
            END DO
         END DO
  100 ENDDO
C
  500 FORMAT(1X)
  501 FORMAT(12X,I1,9X,I1)
  502 FORMAT(2X,I4,2x,6E12.4/8X,6E12.4,4(/4X,6F12.4))
  503 FORMAT(2X,I4,F10.3,5F12.3,(/4X,6F12.3))
  504 FORMAT(2X,I4,F10.3,5F12.2,2(/4X,6F12.2))
  505 FORMAT(10F8.3)
C
      IHYDPR=-IHYDPR
      RETURN
      END IF
C
C ---------------------------------
C     read Lemke tables
C ---------------------------------
C
      ILEMKE=1
      ILINE=0
      IHYDP0=IHYDPR
      IF(IHYDPR.LT.100) THEN
         NTAB=1
         IHYDP(1)=IHYDPR
       ELSE IF(IHYDP0.LT.10000) THEN
         NTAB=2
         IHYDP(1)=IHYDP0/100
         IHYDP(2)=IHYDP0-100*IHYDP(1)
       ELSE IF(IHYDP0.LT.1000000) THEN
         NTAB=3
         IHYDP(1)=IHYDP0/10000
         IHYDP(2)=(IHYDP0-10000*IHYDP(1))/100
         IHYDP(3)=IHYDP0-10000*IHYDP(1)-100*IHYDP(2)
       ELSE
         NTAB=4
         IHYDP(1)=IHYDP0/1000000
         IHYDP(2)=(IHYDP0-1000000*IHYDP(1))/10000
         IHYDP(3)=(IHYDP0-1000000*IHYDP(1)-10000*IHYDP(2))/100
         IHYDP(4)=IHYDP0-1000000*IHYDP(1)-10000*IHYDP(2)-
     *            100*IHYDP(3)
      END IF
c      
      ILINE=0
      DO ITAB=1,NTAB
      ILINEB=ILINE
      IHYDPR=IHYDP(ITAB)
      READ(IHYDPR,*) NLLY
      HEADER: DO ILI=1,NLLY
         ILINE=ILINE+1
         READ(IHYDPR,*) I,    J,   ALMIN, ANEMIN,
     *                  TMIN, DLA, DLE,
     *                  DLT,  NWL, NE,    NT
         WL0=WLINE(I,J)
         ILIN0(I,J)=ILINE
         NWLH(ILINE)=NWL
         NWLHYD(ILINE)=NWL
         NTH(ILINE)=NT
         NEH(ILINE)=NE
         DO INE=1,NE
            XNE(INE,ILINE)=ANEMIN+(INE-1)*DLE
         ENDDO
         DO IWL=1,NWL
            WL(IWL,ILINE)=ALMIN+(IWL-1)*DLA
            WLHYD(ILINE,IWL)=WL(IWL,ILINE)
            WL(IWL,ILINE)=EXP(2.3025851*WL(IWL,ILINE))
         ENDDO

         DO IT=1,NT
            XT(IT,ILINE)=TMIN+(IT-1)*DLT
         ENDDO
      ENDDO HEADER

c
      DO ILI=1,NLLY         
         ILNE=ILINEB+ILI
         call assert(ILINE<=MLINH)
         NWL=NWLH(ILNE)
         READ(IHYDPR,500)
         DO INE=1,NEH(ILNE)
            DO IT=1,NTH(ILNE)
c              READ(IHYDPR,506) QLT,(PRF(IWL,IT,INE),IWL=1,NWL)
               READ(IHYDPR,*) QLT,(PRF(IWL,IT,INE),IWL=1,NWL)
            END DO
         END DO
  506    FORMAT((f8.4,8f9.4))
C
C        coefficient for the asymptotic profile is determined from
C        the input data
C
         XCLOG=PRF(NWL,1,1)+2.5*WLHYD(ILNE,NWL)-0.477121
         XKLOG=0.6666667*XCLOG
         XK=EXP(XKLOG*2.3025851)
C
         DO ID=1,ND
C
C           temperature is modified in order to account for the
C           effect of turbulent velocity on the Doppler width
C
            T=TEMP(ID)+6.06E-9*VTURB(ID)




            ANE=ELEC(ID)
            TL=LOG10(T)
            ANEL=LOG10(ANE)
            F00=1.25E-9*ANE**0.666666667
            FXK=F00*XK
            DOP=1.E8/WL0*SQRT(1.65E8*T)
            DBETA=WL0*WL0/2.997925E18/FXK
            BETAD=DBETA*DOP
C
C       interpolation to the actual values of temperature and electron
C       density. The result is stored at array PRFHYD, having indices
C       ILINE - line number
C       ID    - depth index
C       IWL   - wavelength index 
C

            DO IWL=1,NWL
               CALL INTHYD(PROF,TL,ANEL,IWL,ILNE,ILEMKE)
               PRFHYD(ILNE,ID,IWL)=PROF
            END DO
        END DO


      END DO
      END DO
!       CALL INTHYD(PROF,3.398d0,10.5d0,1,ILIN0(2,3),1)
!       print *,PROF,WL(1,ILIN0(2,3))
C
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C
!*** micha: add key ILEMKE to the interface
      SUBROUTINE INTHYD(W0,X0,Z0,IWL,ILINE,ILEMKE)
      use MOD_HYDTAB
      implicit none
      real*8,intent(out) :: W0
      real*8,intent(in)  :: X0,Z0
      integer,intent(in) :: IWL,ILINE,ILEMKE
!       real*8 :: WL,XT,XNE,PRF
!       integer:: NWLH,NTH,NEH
!       COMMON/VCSDAT/WL(36,8),XT(7,8),XNE(11,8),PRF(36,7,11),
!      *              NWLH(8),NTH(8),NEH(8)
      real*8 ::     XK,FXK,BETAD,DBETA
      COMMON/AUXVCS/XK,FXK,BETAD,DBETA

C
C     Interpolation in temperature and electron density from the
C     Schoening and Butler tables for hydrogen lines to the actual valus of
C     temperature and electron density
C
      integer :: NX,NZ,NT,NE,IZZ,IPZ,N0X,N1X,N0Z,N1Z,I0Z,IX,IPX,I0
      real*8  :: BETA,A,DIV
      INCLUDE '../inc/PARAMS.FOR'
      real*8,parameter :: TWO=2d0
      real*8,dimension(3) :: ZZ,XX,WX,WZ
C
      NX=3
      NZ=3
      NT=NTH(ILINE)
      NE=NEH(ILINE)
      BETA=WL(IWL,ILINE)/FXK
      IF(ILEMKE.EQ.1) THEN
         BETA=WL(IWL,ILINE)/XK
         NX=2
         NZ=2
      END IF
C
C     for values lower than the lowest grid value of electron density
C     the profiles are determined by the approximate expression
C     (see STARKA); not by an extrapolation in the HYD tables which may
C     be very inaccurate
C
      IF(Z0.LT.XNE(1,ILINE)*0.99.OR.Z0.GT.XNE(NE,ILINE)*1.01) THEN
         !*** changed: add BETAD
         CALL DIVSTR(BETAD,A,DIV)
         W0=STARKA(BETA,BETAD,A,DIV,TWO)*DBETA
         W0=LOG10(W0)
         GO TO 500 ! RETURN
      END IF
C
C     Otherwise, one interpolates (or extrapolates for higher than the
C     highes grid value of electron density) in the HYD tables
C
      DO 10 IZZ=1,NE-1
         IPZ=IZZ
         IF(Z0.LE.XNE(IZZ+1,ILINE)) GO TO 20
   10 ENDDO
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
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
         IF(X0.GT.1.01*XT(NT,ILINE).AND.BETAD.GT.10.) THEN
            CALL DIVSTR(BETAD,A,DIV)
            W0=STARKA(BETA,BETAD,A,DIV,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500 ! RETURN
         END IF
C
C     Otherwise, normal inter- or extrapolation
C
C     Both interpolations (in T as well as in electron density) are
C     by default the quadratic interpolations in logarithms
C
         !*** Find Temperature Index
         DO 30 IX=1,NT-1
            IPX=IX
            IF(X0.LE.XT(IX+1,ILINE)) GO TO 40
   30    ENDDO
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
         N1X=N0X+NX-1
         !*** Find 
         DO 200 IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XT(IX,ILINE)
            WX(I0)=PRF(IWL,IX,IZZ)
  200    ENDDO
         IF(WX(1).LT.-99..OR.WX(2).LT.-99..OR.WX(3).LT.-99.) THEN
            CALL DIVSTR(BETAD,A,DIV)
            W0=STARKA(BETA,BETAD,A,DIV,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500 ! RETURN
          ELSE
            WZ(I0Z)=YINT(XX,WX,X0)
         END IF
  300 ENDDO
      W0=YINT(ZZ,WZ,Z0)
  500 CONTINUE
      RETURN
      END SUBROUTINE
!************************
 !     SUBROUTINE GBF_AVE(N, Zi, hv, gbf_all, gbf_av) 
 !     IMPLICIT DOUBLE PRECISION (A -H,O-Z)
 !     DIMENSION gbf_all(*) 
 !     gbL_av =0.D0  
 !     DO L=0, N-1 
 !     CALL GAUNT_BF(N, L, Zi, hv, gbf) 
 !     gbf_all(L + 1) =gbf 
 !     gbf_av = gbf_av +(2*L + 1)*gbf
 !     ENDDO 
 !     gbf_av =gbf_av/DFLOAT(N**2) 
 !     RETURN
 !     END
      


  !    SUBROUTINE GAUNT_BF(N, L, Zi, hv, gbf) 
  !    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  !    DATA RY, PI/13.60580436D0, 3.14l592653D0/ 
  !    Ee=hv-Zi**2.D0*RY/DFLOAT(N**2) 
  !    gbf=0.D0 
  !    IF(Ee.LE.1.D-18) RETURN
  !     ETA=DSQRT(Zi**2.D0*RY/Ee) 
  !     RHO=ETA/DFLOAT(N) 
  !     SUM=0.D0
  !     DO K=1,L+1
  !     SUM =SUM + DLOG(K**2.D0 + ETA**2.D0) 
  !     ENDDO
  !     R12=1.D0 + RHO**2.D0 
  !     EL2=ETA**2.D0+(L+11**2.D0
  !     A1=DEXP(-2.D0*ETA*(PI/2.D0-DARTAN(RHO)) + 0.5D0*(FACLN(N + L)
  !     -FACLN(2*L + 1)— FACLN(2*L + 2)—FACLN(N —L —I) + SUM ±



!***********************




      end module
