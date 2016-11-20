      module MOD_WMOM0_F
      contains
      SUBROUTINE WMOM0_F (ND,NP,R,p,WLK)
C***  The Following Modul computes the Integration Weights for the
C***  calculation of the 0-order Moment in COMPXJ
     
C***  INTEGRATION OF THE ZERO-MOMENT OF THE RADIATION FIELD (MEAN INTENSITY)
C***  THE INTEGRATION WEIGHTS ARE GENERATED (ARRAY WLK)
C***  RADIUS-MESH R and P-mesh ARE GIVEN
C***  WEIGHTS ARE ACCORDING TO TRAPEZOIDAL RULE IN Z=SQRT(R*R-P*P)
      implicit real*8(a-h,o-z)
     
c-old      DIMENSION R(ND),Z(ND,NP),WLK(ND,NP)
      DIMENSION R(ND),p(NP),WLK(ND,NP)
      DO L=1,ND
        DO JP=1,NP
          WLK(L,JP)=0.0
        ENDDO
      ENDDO
      DO L=1,ND
        JMAX=NP-(L-1)
        RL2=2.*R(L)
        rr =r(l)*r(l)
C***  FIRST STEP
c-old        ZJ=Z(L,1)
c-old        ZNEXT=Z(L,2)
        zj    = sqrt(rr - p(1)*p(1) )
        znext = sqrt(rr - p(2)*p(2) )
        WLK(L,1)=(ZJ-ZNEXT)/RL2
C***  MIDDLE STEPS
        DO 1 J=3,JMAX
        ZLAST=ZJ
        ZJ=ZNEXT
c-old        ZNEXT=Z(L,J)
        znext = sqrt(rr-p(j)*p(j))
        WLK(L,J-1)=(ZLAST-ZNEXT)/RL2
    1   CONTINUE
C***  LAST STEP, IMPLYING Z(L,JMAX)=.0
    2   WLK(L,JMAX)=ZJ/RL2
    3   CONTINUE
      ENDDO

      RETURN
      END subroutine
      end module