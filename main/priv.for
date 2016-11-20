      module MOD_PRIV
      contains
      SUBROUTINE PRIV (K,XL,LPRIV,ND,NP,RADIUS,Z,VJL,RSTAR)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION Z(ND,NP),VJL(NP,ND) 
      DIMENSION RADIUS(ND) 
 
      IF (ND.LE.2) STOP 'ND.LE.2' 
      IF (K.EQ.1) THEN 
      PRINT 10 
   10 FORMAT (/, 
     $10X,'VARIABLE V AS A FUNCTION OF MU, DEPTH AND WAVELENGTH',/, 
     $40X,'===================================================='/) 
      PRINT *,RSTAR,' RC (CORE RADIUS IN CM)' 
      NPRPT=(ND-2)/LPRIV 
      PRINT *,NPRPT,' PRINTED DEPTH POINTS' 
      DO 111 L=2,ND-1 
      IF(((L-1)/LPRIV)*LPRIV.NE.(L-1) .AND. L.NE.ND) GOTO 111 
      JMAX=NP+1-L 
      PRINT *,L,JMAX,RADIUS(L)*RSTAR,RADIUS(L),' L, JMAX, R, R/RC; MU=' 
   11 FORMAT (10X,' L=',I7,'  RADIUS=',1PE15.5) 
      PRINT 12,( (Z(L,J)/RADIUS(L)),J=1,JMAX) 
   12 FORMAT (8(F10.7)) 
C  12 FORMAT (10X,' MU (J) =',/, (10(1PE13.6)) ) 
  111 CONTINUE 
      ENDIF 
 
      PRINT *,XL,' WAVELENGTH IN A' 
      DO 1 L=2,ND-1 
      IF(((L-1)/LPRIV)*LPRIV.NE.(L-1) .AND. L.NE.ND) GOTO 1 
C***  ALL NON-BOUNDARY POINTS  L= 2 ... ND-1 
      JMAX=NP+1-L 
      PRINT  *,L,JMAX,' L,JMAX; V-NU(J,L)=' 
   14 FORMAT ('   L=',I5,'   J=    1  TO ',I5) 
      PRINT 13,(VJL(J,L),J=1,JMAX) 
   13 FORMAT (8(1PE10.3)) 
C  13 FORMAT (10X,' V-NU (J,L) =',/,(10(1PE13.6))) 
    1 CONTINUE 
 
      RETURN 
      END subroutine

      end module