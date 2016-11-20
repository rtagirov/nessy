      module MOD_CCORE
      contains
C**********  MODULNAME: CCORE     ******* 06/08/87  20.27.56.******    87 KARTEN
      SUBROUTINE CCORE (WCHARM,NF,DELTAC,IPRICC,MODHEAD,JOBNUM,
     $      SCOLD,RADIUS,XLAMBDA,ND,T,RNE,POP1,ENTOT,RSTAR,
     $      OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $      NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,SIGMAKI,
     $      WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2,NFDIM)
C*******************************************************************************
C***  DETERMINE THE FREQUENCY WEIGHTS WCHARM REPRESENTING THE APPROXIMATE
C***  LAMBDA-OPERATORS AT EACH CONTINUUM FREQUENCY POINT
C***  THE CONT. SOURCE FUNCTION IS CALCULATED FROM OLD POP.NUMBERS (=SCOLD)
C***  INITIALISED/CHANGED: WCHARM, SCOLD
C***			in coop_m: OPA, ETA, THOMSON, IWARN, MAINPRO,MAINLEV
C***	14.03.2006: removed calculation of WCHARM (micha)
C*******************************************************************************
!--- HMINUS --------
!1	 STEAL			     128, 217
!2     | LINPOP 			     212
!:	 : :					:
!3     | | CCORE 		      	93
!4     | | | COOP_M 		      51
!5     | | | | PHOTOCS_M 	     186
!6     | | | | | cstabread 	      60
!7     | | | | | | intpl 	      28
!5     | | | | GAUNTFF 		     218
!5     | | | | hminusff            253
!5     | | | | LINSCA 		     281
!5     | | | | RDOPAC 		     308
!4     | | | PRICC 		     127
      use MOD_COOP_M
      use MOD_PRICC
	implicit NONE
	! GLOBAL Variables intent(out)
	integer,intent(in) :: JOBNUM
      integer :: NFDIM, NDIM,N
	integer :: NOM(N),  IPRICC
	integer :: ND, NF	! Number of depth- & frequency-points
	
	integer :: XLBKG1,XLBKG2

	real*8,dimension(ND,NF):: WCHARM!the Lambda Operator
	! gobal variables, intent(inout)  (inout because of subprocedures)
	real*8,dimension(ND) :: OPA,ETA,THOMSON
	real*8,dimension(ND) :: RADIUS  ! = [1+delta, ... ,1] given in relative solar units
	
	character,dimension(*) :: MAINPRO*10, MAINLEV*10, LEVEL*10
	real*8,dimension(NF,ND) :: SCOLD
	!global variables, intent(in)
	character::MODHEAD*104
	real*8,dimension(NDIM,NFDIM) :: WAVARR,SIGARR
	real*8,dimension(NF)    :: XLAMBDA	     ! Array of Wavelength
	logical :: LBKG	
	real*8 :: DELTAC, RSTAR
	
	! local variables
	real*8 :: OPAL	! OPAL =(OPA*(1-THOMSON))(L)
	real*8 :: XLAM   			! wavelength at point k(@see XLAMBDA)
      real*8,allocatable :: DUMMY1(:)
      character*8,allocatable :: CDUMMY1(:)
	integer :: K, L			! loop variables			
!	integer iostat
! micha: added dimensions for modules
      integer, dimension(*) :: NCHARG, IWARN
      real*8, dimension(*)  :: T, RNE, POP1, ENTOT,WEIGHT,ELEVEL,EION
      real*8, dimension(*)  :: EINST,SIGMAKI
C***  MARGIT HABERREITER
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
     
C***  LOOP OVER ALL CONTINUUM FREQUENCIES  *************************************
      DO K=1,NF
        XLAM=XLAMBDA(K)
        CALL  COOP_M (XLAM,ND,T,RNE,POP1  ,ENTOT,RSTAR,
     $			  OPA,ETA,THOMSON,IWARN,MAINPRO,MAINLEV,NOM,
     $			  NDIM,N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $			  DUMMY1,DUMMY1,CDUMMY1,K,NF,SIGMAKI,WAVARR,SIGARR,
     $			  LBKG,XLBKG1,XLBKG2,NFDIM)
C***  SOURCE FUNCTION WITH OLD POP.NUMBERS (WITHOUT THOMSON OPACITY)
	  DO L=1,ND
C***  LASER SECURITY - if medium acts as a laser reset it
	    OPAL=OPA(L)*(1d0-THOMSON(L))
	    IF (OPAL .GT. 0d0) THEN
	      SCOLD(K,L)=ETA(L)/OPAL
	    ELSE
	      SCOLD(K,L)=.0d0
	    ENDIF
	  ENDDO
      ENDDO
      
C***  ENDLOOP  *********************************************************
      IF (IPRICC .EQ. 1)  ! IPRICC == "PRINT CCORE" in Cards
     $	CALL PRICC (ND,NF,WCHARM,DELTAC,MODHEAD,JOBNUM)
	
      RETURN
      END subroutine
      end module
