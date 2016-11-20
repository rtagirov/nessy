      module MOD_CUBIC
      contains
      SUBROUTINE CUBIC (L,LTOT,ZRAY,XCMF,P1,P2,P3,P4)
C***  CUBIC INTERPOLATION WITH FIXED DERIVATIVES AT THE MESH POINTS
C***  THESE DERIVATIVES ARE CALCULATED BY LINEAR I5TERPOLATION BETWEEN
C***  NEIGHBOURING POINTS
C***  XCMF(L) = GRID OF MESH POINTS
C***  ZRAY(L) = FUNCTION TO BE INTERPOLATED
C***  P1,P2,P3,P4 = RESULTING COEFFICIENTS
C***  EVALUATION OF THE INTERPOLATED VALUE Z(X) :
C***  Z(X)=P1*(X-XCMF(L-1))**3 + P2*(X-XCMF(L-1))
C***                + P3*(XCMF(L)-X)**3 + P4*(XCMF(L)-X)
C***  THE COEFFICIENT MATRIX H HAS BEEN INVERTED ANALYTICALLY
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER ( ONE = 1.D+0, two = 2.d0, three = 3.D0 )
     
      DIMENSION ZRAY(LTOT),XCMF(LTOT)
     
      IF (L .LT. 2 .OR. L .GT. LTOT) STOP 'ERROR'
     
C***  SET UP THE COEFFICIENT MATRIX
      D1=one/(XCMF(L)-XCMF(L-1))
      D2=D1*D1
      D3=D1*D2
      D23=D2/three
      H11=D3
      H12=-D3
      H13=D23
      H14=two*D23
      H21=-D1
      H22=two*D1
      H23=-0.333333333333333d0
      H24=-0.666666666666666d0
      H31=-D3
      H32=D3
      H33=-two*D23
      H34=-D23
      H41=two*D1
      H42=-D1
      H43=0.666666666666666d0
      H44=0.333333333333333d0
C***  FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
      LA=MAX0(L-2,1)
      LB=MIN0(L+1,LTOT)
C***  FUNCTION TO BE INTERPOLATED: ZRAY
c     the two Y-values at the edges of the interval: F1, F2
      F1=ZRAY(L-1)
      F2=ZRAY(L)
      if (L-2 .eq. LA) then
c        derivatives given by the parabola through L-2, L-1, and L...
         him1 = (XCMF(L-1)-XCMF(LA))
         hi   = (XCMF(L)-XCMF(L-1))
         sim1 = (ZRAY(L-1)-ZRAY(LA))/him1
         si = (ZRAY(L)-ZRAY(L-1))/hi
c        the formula of m. steffen 1990 A&A 239, 443
         pi = (sim1*hi+si*him1)/(him1+hi)
         f3 = (sign(one,sim1)+sign(one,si))*
     &         min(min(abs(sim1),abs(si)),0.5d0*abs(pi))
c         print *,L-1,him1,sim1,hi,si,pi,f3
      else
         F3=(ZRAY(L)-ZRAY(LA))/(XCMF(L)-XCMF(LA))
c         print *,f3
      endif
c     ...  and L-1, L, and L+1
      if (L+1 .eq. LB) then
c        derivatives given by the parabola through L-2, L-1, and L...
         him1 = (XCMF(L)-XCMF(L-1))
         hi   = (XCMF(LB)-XCMF(L))
         sim1 = (ZRAY(L)-ZRAY(L-1))/him1
         si = (ZRAY(LB)-ZRAY(L))/hi
c        the formula of m. steffen 1990 A&A 239, 443
         pi = (sim1*hi+si*him1)/(him1+hi)
         f4 = (sign(one,sim1)+sign(one,si))*
     &         min(min(abs(sim1),abs(si)),0.5d0*abs(pi))
c         print *,L,him1,sim1,hi,si,pi,f4
      else
         F4=(ZRAY(LB)-ZRAY(L-1))/(XCMF(LB)-XCMF(L-1))
c         print *,f4
      endif
C***  CALCULATE POLYNOMIAL COEFFICIENTS: P(VECTOR) = H(MATRIX) * F(VECTOR)
      P1=H11*F1+H12*F2+H13*F3+H14*F4
      P2=H21*F1+H22*F2+H23*F3+H24*F4
      P3=H31*F1+H32*F2+H33*F3+H34*F4
      P4=H41*F1+H42*F2+H43*F3+H44*F4
     
      RETURN
      END subroutine
      end module
