      module MOD_INTPLCS
      contains
	subroutine intplcs(SIGINT,FREQ,WAVARR,SIGARR,ITR,NDIM,NFDIM)
c	SUBROUTINE intplcs(SIGINT,FREQ,FREQ0,ITR) 
c           CALL intplcs(SIGK,FR,FR0,WAVARR,SIGARR,ITR) 
C*** CALLED BY SIGK
C*********************************************************
C***	written by Margit Haberreiter 
C***  TO INTERPOLATE FROM THE COARSE FREQUENCY GRID DATA
C***  WAVARR	:WAVENUMBER IN CM^-1
C***  SIGARR	:CROSS SECTION IN MBARN(?)
C***           INTERPOLATED FROM THE OPACITY AND IRON PROJECT DATA
C***  FR   			: FREQUENCY CONSIDERED IN ANSTROM
C***  FR0 			: THRESHOLD FREQUENCY
C***  ITR				: ATOMIC LEVEL
C***  SIGINT          : LINEARILY INTERPOLATED OPACITY
C***                    CROSS SECTION AT PARTICULAR LEVEL AND FREQUENCY POINT
C**************************
      USE MOD_INTPL
c	INCLUDE '../TLUSTY/PARAMS.FOR'
c     INCLUDE '../TLUSTY/MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ITR,NDIM,NFDIM
      real*8, intent(in   ) :: FREQ,WAVARR,SIGARR
      real*8, intent(  out) :: SIGINT
!       integer :: 
!       real*8  :: 
        
	INCLUDE '../INCLUDE2/PARAMS.FOR'
	INCLUDE '../INCLUDE2/MODELP.FOR'
	DIMENSION WAVARR(NDIM,NFDIM),SIGARR(NDIM,NFDIM)
	real*8::WVNMBR,WAV1,WAV2,SIG1,SIG2
	LOW = ITR
C     CL: C IN CGS UNITS
C     WVNMBR IN CM-1
	WVNMBR = FREQ/CL
C***********************************************************************	
		SIGINT=0.
c		do  k=1, FGARR
		do  k=1, NFDIM-1
			if ((wavarr(low,k) .ge. WVNMBR) .AND. 
     $			(wavarr(low,k+1) .le. WVNMBR)) then
				WAV1=wavarr(low,k+1)
				WAV2=wavarr(low,k)
				SIG1=SIGARR(low,K+1)
				SIG2=SIGARR(low,K)
				call intpl(SIGINT,WVNMBR,WAV1,WAV2,SIG1,SIG2)
C					 INTPL(Y,     X,     X1,  X2,  Y1,  Y2)
c				PRINT *,'!!!FOUND INTERPOLATION ROUTINE FOR',ITR
c				print *,WAV1,WAV2,WVNMBR
c				PRINT *,SIG1,SIG2,SIGINT
				SIGINT=SIGINT*1.d-18
				GOTO 1
			ENDIF
		end do
    1 CONTINUE
	END subroutine
      end module
