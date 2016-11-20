      module MOD_BNUE
      contains
      PURE FUNCTION BNUE(XLAMBDA,T)
      !use MOD_ERROR
C*** PLANCK FUNCTION, LAMBDA IN ANGSTROEM, T IN KELVIN          ****************
C***  BNUE IN CGS UNITS: ERG PER CM**2, PER SEC AND PER HERTZ       ************
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C   (DIMENSION ANGSTROEM**3  * ERG/SEC/HZ/CM**2)
      IMPLICIT NONE
      real*8 :: BNUE, cxt
      real*8,intent(in) :: XLAMBDA,T
      real*8,parameter  :: C1=1.438831D8
      real*8,parameter  :: C2=3.972967D+8
     
      IF (T .LE. .0) THEN
        BNUE=.0
        !print '(A,f0.2,X,f0.2)','T,XLAMBDA = ',T,XLAMBDA
      ELSE
        CXT=C1/XLAMBDA/T
        IF (CXT.GT.500.) THEN
          BNUE=0.0
        ELSE
          BNUE=C2/(EXP(CXT)-1d0)/XLAMBDA**3 !/XLAMBDA/XLAMBDA
        ENDIF
      ENDIF
      END function
      end module