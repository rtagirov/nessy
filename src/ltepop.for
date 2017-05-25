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

 
C***  ENDLOOP  ---------------------------------------------------------
    
      RETURN
      END subroutine


      SUBROUTINE LTEPOP_NEW (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
C***********************************************************************
C***  POPULATION NUMBERS IN THERMODYNAMIC EQUILIBRIUM FOR ALL ELEMENTS
C***  CALLED BY GREYM,POPZERO
C***********************************************************************
      use mod_error
      use MOD_PARTF2
      IMPLICIT REAL*8(A-H,O-Z) 
      REAL*8 prt1, prt2, prt3
      DIMENSION ENLTE(N), WEIGHT(N),NCHARG(N),EION(N),ELEVEL(N),NOM(N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
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
      
      SUMPC=0.
      DO J=NFIRNA,NLANA
      SUMPC=SUMPC+WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL)
      ENDDO

      I=NOM(NFIRNA)
 
      
      CALL PARTF2(I,1,TL,PRT1)
      CALL PARTF2(I,2,TL,PRT2)
      CALL PARTF2(I,3,TL,PRT3)
  
      SUMPTH=PRT1+PRT2+PRT3
 !    print *, 'TEST'
 !    print *, J, PRT1, PRT2, PRT3

      corr=SUMPC/SUMPTH

      DO J=NFIRNA,NLANA
      ENLTE(J)=ENLTE(J)*corr
      ENDDO


    9 CONTINUE

      
C***  ENDLOOP  ---------------------------------------------------------
     
      RETURN
      END subroutine

      SUBROUTINE LTEPOP_NEW1 (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                  ABXYZ,NFIRST,NLAST,NATOM)
C***********************************************************************
C***  POPULATION NUMBERS IN THERMODYNAMIC EQUILIBRIUM FOR ALL ELEMENTS
C***  CALLED BY GREYM,POPZERO
C***********************************************************************
      use mod_error
      use MOD_PARTF2
      IMPLICIT REAL*8(A-H,O-Z) 
      REAL*8 prt1, prt2, prt3
      DIMENSION Eip(3,30)
      DIMENSION ENLTE(N), WEIGHT(N),NCHARG(N),EION(N),ELEVEL(N),NOM(N)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
C***  C1 = H * C / K    ( CM * KELVIN )
      DATA C1/1.4388d0/
C***  C2 = FACTOR IN SAHA EQ.   ( C.F. MIHALAS P.113 )
      DATA C2/2.07d-16/
      DATA evtocm/8055.2/
      DATA Eip/'  H ',  0.756095, 13.6157,
     * '  HE',24.6190, 54.4876, 
     * '  LI', 5.39864, 75.619,
     * '  BE', 9.33467, 18.206,
     * '  B ', 8.30868, 25.149,
     * '  C ',11.2748, 24.376,
     * '  N ',14.5638, 29.593,
     * '  O ',13.6355, 35.108,
     * '  F ',17.4452, 34.980,
     * '  NE',21.5922, 41.070,
     * '  NA', 5.1457, 47.290,
     * '  MG', 7.65605, 15.030, 
     * ' AL ', 5.99342, 18.823, 
     * ' SI ', 8.16215, 16.350, 
     * ' P  ',10.5002, 19.720,
     * ' S  ',10.3733, 23.400,
     * ' CL ',12.9843, 23.800,
     * ' AR ',15.7798, 27.620,
     * ' K  ', 4.34624, 31.810,
     * ' CA ', 6.12101, 11.8870,
     * ' SC ', 6.56992, 12.800,
     * ' TI ', 6.82913, 13.630,
     * ' V  ', 6.74844, 14.200,
     * ' CR ', 6.77520, 16.490,
     * ' MN ', 7.44356, 15.640,
     * ' FE ', 7.91253, 16.183,
     * ' CO ', 7.87442, 17.060,
     * ' NI ', 7.64723, 18.168, 
     * ' CU ', 7.7363, 20.292,
     * ' ZN ', 9.40626, 17.964/
  
  


    
      T32=TL*SQRT(TL)
C***  LOOP FOR EACH ELEMENT  -------------------------------------------
      DO 9 NA=1,NATOM
      
      if (NA .eq. 1) then
      CALL PARTF2(NA,1,TL,PRT2)
      CALL PARTF2(NA,2,TL,PRT3)
      CALL PARTF2(NA,3,TL,PRT1)
      goto 10
      endif


      CALL PARTF2(NA,1,TL,PRT1)
      CALL PARTF2(NA,2,TL,PRT2)
      CALL PARTF2(NA,3,TL,PRT3)

 10   F1=EXP(-C1*evtocm*Eip(2,NA)/TL)*T32/(C2*ENE)
      F2=EXP(-C1*evtocm*Eip(3,NA)/TL)*T32/(C2*ENE)
      
      denom=PRT1+PRT2*F1+PRT3*F1*F2
      coeff1=1./denom
      coeff2=coeff1*F1
      coeff3=coeff2*F2


      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)

      if (NA .eq. 1) then

      DO 1 J=NFIRNA, NLANA
    
      IF (NCHARG(J) .EQ. -1) THEN
      ENLTE(J)=coeff1*WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL) 
      ELSEIF (NCHARG(J) .EQ. 0) THEN
      ENLTE(J)=coeff2*WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL)
      ELSEIF (NCHARG(J) .EQ. 1) THEN 
      ENLTE(J)=coeff3*WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL)
      ENDIF
      ENLTE(J)=ENLTE(J)*ABXYZ(NA)

   1  CONTINUE
      
      else

      DO 2 J=NFIRNA, NLANA
C***  BOLTZMANN FACTOR
      IF (NCHARG(J) .EQ. 0) THEN
      ENLTE(J)=coeff1*WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL) 
      ELSEIF (NCHARG(J) .EQ. 1) THEN
      ENLTE(J)=coeff2*WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL)
      ELSEIF (NCHARG(J) .EQ. 2) THEN 
      ENLTE(J)=coeff3*WEIGHT(J)*EXP(-C1*ELEVEL(J)/TL)
      ENDIF
      ENLTE(J)=ENLTE(J)*ABXYZ(NA)

      

 2    CONTINUE    

      endif


    9 CONTINUE

      
C***  ENDLOOP  ---------------------------------------------------------
     
      RETURN
      END subroutine





      end module
