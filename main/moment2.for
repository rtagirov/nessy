      module MOD_MOMENT2
      contains
      SUBROUTINE MOMENT2 (R,JMAX,P,U,XK)
C***  INTEGRATION OF THE 2. MOMENT XK OF THE RADIATION FIELD U
C***  FEAUTRIER-INTENSITY U(J), IMPACT PARAMETER MESH P(J)
C***  AND RADIUS POINT R ARE GIVEN.
C***  WEIGHTS ARE ACCORDING TO TRAPEZOIDAL RULE IN Z*Z*DZ, Z=SQRT(R*R-P*P)
      implicit real*8(a-h,o-z)

      DIMENSION P(JMAX),U(JMAX)
      RR=R*R
C***  FIRST STEP, IMPLYING P(1)=0
      Z=R
      ZQ=RR
      PJ=P(2)
      ZNQ=RR-PJ*PJ
      ZNEXT=SQRT(ZNQ)
      W=Z*(3*ZQ-ZNQ)-ZNEXT*(ZQ+ZNQ)
      XK=W*U(1)
C***  MIDDLE STEPS
      DO 1 J=3,JMAX
      ZLAST=Z
      ZLQ=ZQ
      Z=ZNEXT
      ZQ=ZNQ
      PJ=P(J)
      ZNQ=RR-PJ*PJ
      ZNEXT=SQRT(ZNQ)
      W=Z*(ZLQ-ZNQ)+ZLAST*(ZLQ+ZQ)-ZNEXT*(ZQ+ZNQ)
    1 XK=XK+W*U(J-1)
C***  LAST STEP, IMPLYING P(JMAX)=R
      W=Z*ZQ
      XK=XK+W*U(JMAX)
      XK=XK/R/RR/12.
      RETURN
      end subroutine
      end module