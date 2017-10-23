      MODULE MOD_PHOTOCS

      contains

      SUBROUTINE PHOTOCS(SIGMA,SIGMATH,EDGE,WAVENUM,ALPHA,SEXPO,AGAUNT,LOW,WAVARR,SIGARR,N,NFDIM)

C***  CALCULATES SIGMA(NUE), THE FREQUENCY DEPENDENT PHOTO CROSS SECTION
C***  Returns SIGMA

      USE MOD_CSTABREAD

      implicit none

      REAL*8, intent(out) :: SIGMA
      INTEGER,intent(in) :: N, NFDIM, LOW
      REAL*8,intent(in)  :: SIGMATH,EDGE,WAVENUM
      REAL*8,dimension(N),intent(in) :: ALPHA,SEXPO
      character*8,intent(in) :: agaunt(N)

      REAL*8,dimension(N, NFDIM),intent(in) :: WAVARR, SIGARR

      !local variables
      REAL*8 :: CN, DEN, U, GAUNT, FLAMDA
      REAL*8 :: X, X3, XLAM, X0LN2, X0LN, XLN, XLN2
      INTEGER :: I
      
C******************************************************************************
C***  CHANGES BY MARGIT HABERREITER, 20 MAY, 2002
CMH  LEVLOW NEEDS TO BE DEFINED, AS IT IS USED AS A KEYWORD TO SELECT THE 
CMH  ELEMENT AND LEVEL TO READ THE CONTINUUM OPACITIES FROM AN INPUT TABLE
C      CHARACTER*10 LEVLOW
C******************************************************************************

      REAL*8, PARAMETER :: ONE = 1.D+0
      
      DIMENSION CN(6)

C***  H- FROM JOHN A&A193,189 (1988)
      DATA CN/152.519d0,49.534d0,-118.858d0,92.536d0,-34.194d0,4.982d0/
C***  THE FOLLOWING DATA ARE FOR MIHALAS' GAUNT FACTOR FIT ( HI AND HE II, N=1)
      REAL*8 :: A0,A1,A2,A3, AM1, AM2
      DATA A0,A1,A2,A3,AM1,AM2 /
     $ 1.2302628d0,  -3.1927214d-2, 8.9105122d-4, -1.1544111d-5,
     $ -0.50812150d0, 0.10631895d0 /   
      
      X = EDGE / WAVENUM

C***  VARIABLE AGAUNT IS MISUSED TO CARRY THE KEYWORD 'KOESTER' FOR
C***  NEW PHOTOIONIZATION CROSS SECTIONS (KOESTER ET AL. 1985, A+A 149, 423)
C***  FIT COEFFICIENTS FOR THE MODIFIED (!) FORMULA:

      IF (AGAUNT(LOW) .EQ. 'KOESTER') THEN

          XLN=LOG(1.d8/WAVENUM)
          XLN2=XLN*XLN
          X0LN=LOG(1.d8/EDGE)
          X0LN2=X0LN*X0LN
          SIGMA=SIGMATH*X**ALPHA(LOW)*EXP(SEXPO(LOW)*(XLN2-X0LN2))

C*****************************************************************
C***  changes by Margit Haberreiter
CMH   IN CASE OF TABLE READ CROSS SECTIONS FROM TABLE
CMH   THE VARAIBLE AGAUNT IS MISUSED TO CARRY THE KEYWORD 'TABLE' FOR
CMH   READING THE CONTINUUM CROSS SECTIONS FOR EACH ELEMENT, LEVEL AND WAVENUMBER 
CMH  FROM A TABLE
CMH  else if (AGAUNT(LOW) .EQ. 7HTABLE01) THEN
      else if (AGAUNT(LOW) .EQ. 'TABLE') THEN

        call cstabread(SIGMA, WAVENUM, LOW, WAVARR, SIGARR, N, NFDIM)

C*****************************************************************
CMH  HIER IST DAS PROBLEM, AGAUNT(LOW) WIRD NICHT ERKANNT ??
      else if (AGAUNT(LOW) .EQ. 'H-MINUS') THEN
c***     source: John 1988 A&A 193,189
c        control table: Wishart 1979, MNRAS 187, 59p
         XLAM=1.d4/WAVENUM
         X3=XLAM*XLAM*XLAM
         X=SQRT((WAVENUM-EDGE)*1.d-4)
         if (wavenum.gt.50000.) then
            FLAMDA=127.84486d0*X
         ELSE
            FLAMDA = 0.d0
            DO I=6,1,-1
               FLAMDA=(FLAMDA+CN(I))*X
            ENDDO
         ENDIF

         SIGMA=FLAMDA*X*X*X3*SIGMATH

      ELSE
C***  OLD VERSION: SEATON / HYDROGENIC

          IF (ALPHA(LOW) .NE. .0) THEN
C***       ALPHA(LOW) DEFINED: SEATON FORMULA
           SIGMA=SIGMATH * X**SEXPO(LOW)*(ALPHA(LOW)+(one-ALPHA(LOW))*X)
           ELSE
C***  ALPHA(LOW) NOT DEFINED: HYDROGENIC EXPONENT NUE**(-3)
           SIGMA=SIGMATH*X*X*X
          ENDIF

      ENDIF 
C***  BOUND-FREE GAUNT FACTORS ARE CALCULATED DEPENDING ON KEYWORD AGAUNT
     
C***  POLYNOMIAL FIT FOR GII(N=1) FROM MIHALAS 1967 APJ 149,P187
C***  THIS MIHALAS GAUNT FACTOR IS ONLY VALID FOR GROUND STATE N=1 !
      IF (AGAUNT(LOW) .EQ. 'MIHALAS' ) THEN

            GAUNT=A0+X*(AM1+X*AM2)+(A1+(A2+A3/X)/X)/X
            SIGMA=SIGMA*GAUNT

      ENDIF
C***  GII FIT FROM SEATON (1960), REP. PROG. PHYS. 23, P. 313
C***  THIS FORMULA IS VALID FOR ALL HYDROGENIC LEVELS ( H OR HE II)
C***  VARIABLE SEXPO IS MISUSED TO CARRY THE MAIN QUANTUM NUMBER
      IF (AGAUNT(LOW) .EQ. 'SEATON' ) THEN

            IF (SEXPO(LOW) .LE. .0d0) stop 'ERROR in PHOTOCS: MAIN QUANTUM NUMBER UNDEFINED'

            U=one/X - one
            DEN=(X/SEXPO(LOW))**0.666666666666d0
            GAUNT=one + 0.1728d0 * (U-one) * DEN -
     -            0.0496d0 * (U*(U+1.333333333333d0)+one) * DEN * DEN
            SIGMA=SIGMA*GAUNT

      ENDIF

      RETURN

      END SUBROUTINE

      END MODULE
