      MODULE MOD_PRIMOD

      CONTAINS

      SUBROUTINE PRIMOD(ND,RADIUS,RSTAR,ENTOT,T,VELO,GRADI,NP,
     $                  MODHEAD,JOBNUM,TTABLE,TAUROSS,R23)

!     PRINTOUT OF THE DEPTH-DEPENDENT MODEL SPECIFICATIONS

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)
     
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE
      DIMENSION RADIUS(ND),ENTOT(ND),T(ND),VELO(ND),GRADI(ND)
      DIMENSION TAUROSS(ND)
      LOGICAL TTABLE, SPHERIC
      CHARACTER MODHEAD*104
      
      REAL*8 :: RSTAR

      REAL*8, DIMENSION(ND) :: VEL_DIFF
      REAL*8, DIMENSION(ND) :: GRAD_KMSKM
     
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (10H1$$$$$$$$$,/,1X,  A  ,20X,'JOB NO.',I5,
     $       ///,20X,'M O D E L   S P E C I F I C A T I O N',/,20X,
     $ 38('='),/,20X,'   DEPTH-DEPENDENT QUANTITIES',//,
c     $ ' DEPTH     R-1     RMAX-R    CRITERION        VELOCITY        '
     $ ' DEPTH     R-1     km        HEIGHT km        VELOCITY        '
     $ 'GRADIENT      LOG NUMBER DENSITY    EL. TEMPERATURE  ',
     $ ' LOG TAU-ROSS',/,
     $ ' INDEX                                         ( M/S)     ( M/'
     $ 'S PER RSTAR)     (ATOMS PER CM+3)         (KELVIN)   ',
     $ '   (IN LTE)  ',/)
	R23km=(R23-1.d0)*RSTAR/100000.d0
      DO 20 L=1,ND
      RL1=RADIUS(L)-1.
C***  RL2 IN KM
      RL2=((RADIUS(L)-RADIUS(ND))*RSTAR)/100000.d0
      DENS=LOG10(ENTOT(L))
C      Print *,l,'Tauross(l)=', tauross(l)	
	IF (TAUROSS(L).LE.0.d0) THEN
      ALTAU=-99.
      ELSE
      ALTAU=LOG10(TAUROSS(L))
      ENDIF
   20 PRINT 3,L,RL1,RL2,RL2-R23km,VELO(L)*1000.,GRADI(L)*1000.,
     $	DENS,T(L),ALTAU
c    3 FORMAT(I3, 1X,G14.7,G10.3,3X,A8,F15.4,F15.4,F22.3,F20.0,F13.4)
    3 FORMAT(I3, 1X,G14.7,G11.4,G12.4,F13.4,1x,F15.4,F22.3,F20.0,F13.4)
      PRINT 4,NP
    4 FORMAT (/,' NUMBER OF IMPACT-PARAMETER POINTS  NP =',I3)
      PRINT 9, R23,R23km
    9 FORMAT (/,' RADIUS WHERE TAU-ROSSELAND = 2/3: R23 =',F10.6,
     $' corresponding to km =',F10.3)
      PRINT 11, VPAR2,BETA,VPAR1,RCON,HSCALE,VMIN
   11 FORMAT(/,' VELOCITY LAW:       REF-RADIUS=',F11.4,3H R* ,
     +'   BETA=',F5.1,'   V-INF=',F8.1,' KM/S' ,/,
     +         '              CONNECTION RADIUS=',F11.4,3H R*,/,
     +         ' HYDROSTATIC LAW:  SCALE HEIGHT=',E11.4,3H R*,
     +'   VMIN=',F11.5,' KM/S')
      IF (TTABLE) THEN
         PRINT 6
    6    FORMAT (/,' TEMPERATURE STRUCTURE FROM TABLE INPUT',/)
         ELSE IF (SPHERIC) THEN
            PRINT 7,TEFF
    7       FORMAT (/,' TEMPERATURE STRUCTURE AS FOR A SPHERICAL, ',
     $         'GREY LTE ATMOSPHERE (APPROXIMATELY) WITH TEFF=',F7.0)
            ELSE
            PRINT 8,TEFF
    8       FORMAT (/,' TEMPERATURE STRUCTURE AS FOR A ',
     $         'PLANE-PARALLEL, GREY LTE ATMOSPHERE WITH TEFF=',F7.0)
            ENDIF
     
      IF (TMODIFY .NE. .0) PRINT 5,TMODIFY
    5 FORMAT (/,' TEMPERATURE STRATIFICATION MODIFIED: TMODIFY=',F6.3)
     
      IF (TMIN .GT. .0) PRINT 10, TMIN
   10 FORMAT (' MINIMUM TEMPERATUR SPECIFIED: TMIN=',F7.0)

      DO L = 2, ND - 1

         VEL_DIFF(L) = (VELO(L + 1) - VELO(L - 1)) / 2.0D0

         GRAD_KMSKM(L) = VEL_DIFF(L) * 2.0D0 / (HEIGHT(L + 1) - HEIGHT(L - 1))

      ENDDO

      VEL_DIFF(1) = VELO(2) - VELO(1)

      GRAD_KMSKM(1) = VEL_DIFF(1) / (HEIGHT(2) - HEIGHT(1))

      VEL_DIFF(ND) = VELO(ND) - VELO(ND - 1)

      GRAD_KMSKM(ND) = VEL_DIFF(ND) / (HEIGHT(ND) - HEIGHT(ND - 1))

      OPEN(UNIT = 7, FILE = 'vel_field.out', ACTION = 'WRITE', STATUS = 'REPLACE')

      WRITE(7, 30) 'DI', 'HEIGHT, [km]', 'VEL, [km/s]', 'VEL_DIFF, [km/s]', 'GRAD, [km/s/km]'

      DO L = 1, ND; WRITE(7, 40) L, HEIGHT(L), VELO(L), -VEL_DIFF(L), GRAD_KMSKM(L); ENDDO

      CLOSE(7)

30    FORMAT(1x,A,8x,A,8x,A,7x,A,4x,A,/)
40    FORMAT(I3,5x,ES15.7,5x,ES15.7,5x,ES15.7,5x,ES15.7)
     
      RETURN

      END SUBROUTINE

      END MODULE
