      module MOD_MOMENT1
      contains
      SUBROUTINE MOMENT1 (R,NP,P,U,H)
C***  CALCULATES AN ANGLE INTEGRAL H OF THE RADIATION FIELD U
C***  BESIDES OF THE OUTER BOUNDARY, THIS IS NOT THE 1. MOMENT H,
C***  BUT RATHER AN INTENSITY-LIKE QUANTITY
C***  INTEGRATION WITH TRAPEZOIDAL RULE, WEIGHTS P * DP
      implicit real*8(a-h,o-z)
     
      DIMENSION U(NP),P(NP)
      NPM=NP-1
      A=.0
      B=P(2)
      W=(B-A)*(B+2.*A)
      H=W*U(1)
    1 DO 2 J=2,NPM
      C=P(J+1)
      W=(A+B+C)*(C-A)
      H=H+W*U(J)
    3 A=B
    2 B=C
      W=(B-A)*(2.*B+A)
      H=H+W*U(NP)
    7 H=H/R/R/6.
      RETURN
      end subroutine
      end module