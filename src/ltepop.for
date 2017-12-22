      module MOD_LTEPOP

      contains

      subroutine LTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,ABXYZ,NFIRST,NLAST,NATOM)

!     POPULATION NUMBERS IN THERMODYNAMIC EQUILIBRIUM FOR ALL ELEMENTS
!     CALLED BY GREYM, POPZERO, LINPOP

      use mod_error

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION ENLTE(N), WEIGHT(N), NCHARG(N), EION(N), ELEVEL(N), NOM(N)
      DIMENSION NFIRST(NATOM), NLAST(NATOM)

      real*8, intent(in), dimension(NATOM) :: ABXYZ

C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1/1.4388d0/
C***  C2 = FACTOR IN SAHA EQ.   ( C.F. MIHALAS P.113 )
      DATA C2/2.07d-16/
      T32=TL*SQRT(TL)
C***  LOOP FOR EACH ELEMENT  -------------------------------------------

      DO 9 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      ENLTE(NFIRNA)=1.
      DO 1 J=NFIRNA+1,NLANA
C***  BOLTZMANN FACTOR
      ENLTE(J)=EXP(C1*(ELEVEL(J-1)-ELEVEL(J))/TL)*WEIGHT(J)/WEIGHT(J-1)
     $      * ENLTE(J-1)

!      print*, tl, na, j, enlte(j)
	IF(NOM(J) .NE. NOM(NFIRNA)) STOP 'LTEPOP:WRONG ELEMENT MEMBERSHIP'
      IF (NCHARG(J) .EQ. NCHARG(J-1) ) GOTO 1
      IF (NCHARG(J) .NE. NCHARG(J-1)+1 )
     $            STOP "LTEPOP: INVALID CHARGE DIFFERENCE"
C***  SAHA FACTOR
      ENLTE(J)=ENLTE(J)*EXP(-C1*EION(J-1)/TL)*T32/ENE/C2
    1 CONTINUE
     
C***  NORMALIZATION
      SUM=.0
      DO 2 J=NFIRNA,NLANA
    2 SUM=SUM+ENLTE(J)
      SUM=SUM/ABXYZ(NA)
	
      DO 3 J=NFIRNA,NLANA

    3 ENLTE(J)=ENLTE(J)/SUM     

    9 CONTINUE
    
      RETURN

      END subroutine

      end module
