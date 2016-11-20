      module MOD_ERF_INF
      contains
C**********  MODULNAME: ERF_INF       ******* 24/03/87  21.09.24.******    16 KARTEN
      FUNCTION ERF_INF(X)
C***  INTEGRAL OVER NORMALIZED GAUSS FUNCTION FROM -INFINITY TO X

      implicit real*8(a-h,o-z)

      DIMENSION A(4)
      DATA A / .25482 9592 , -.28449 6736,  1.42141 3741,-1.453152027 /

      Z=ABS(X)
      T=1.D0/(1.D0+.3275911*Z)
      ERF_INF=1.06140 5429D0
      DO 1 I=1,4
      K=5-I
    1 ERF_INF=A(K)+T*ERF_INF
      ERF_INF=T*ERF_INF*EXP(-Z*Z)*.5D0
      IF (X .GT. .0) ERF_INF=1.D0-ERF_INF

      RETURN
      END function
      end module