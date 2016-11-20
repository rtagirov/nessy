      module MOD_PLOTTAB
      contains
      SUBROUTINE PLOTTAB (KANAL,X,Y,N,ISYMBOL)
      USE MOD_PLOTCON
      implicit real*8(a-h,o-z)
      real*8 x(n), y(n)
      CALL PLOTCON (KANAL,X,Y,N,ISYMBOL)

      RETURN
      END subroutine
      end module
