      module MOD_SEQUIN
      contains
C**********  MODULNAME: SEQUIN    ******* 24/03/87  22.08.25.******    70 KARTEN
      SUBROUTINE SEQUIN (NDIM,N,X,K,XNEW)
C***  X(K) MUST BE A STRICKTLY MONOTONIC SEQUENCE
C***   THIS SUBROUTINE INSERTS THE POINT XNEW AT SUITABLE POSITION K
C***  NO INSERTION IN CASE OF EXACT COINCIDENCE (K:= NEGATIVE INDEX )

      IMPLICIT REAL*8(A-H,O-Z)
      integer,intent(in) :: NDIM
      integer,intent(inout) ::N,K
      real*8,intent(in) :: XNEW
      real*8,intent(inout) :: X
      DIMENSION X(NDIM)
     
C***  INDEX ESTIMATE BY BISECTION
      NA=1
      A=X(1)
      NB=N
      B=X(N)
      IF (XNEW.EQ.A) GOTO 9
      IF (XNEW.EQ.B) GOTO 10
      IF((XNEW-A)*(XNEW-B).GE..0) GOTO 4
    1 IF(NB-NA .EQ. 1) GOTO 2
      NH=(NA+NB)/2
      H=X(NH)
C      IF (XNEW.EQ.H) GOTO 11
C***  ONLY INSERT A POINT IF IT IS NOT TOO CLOSE TO ANOTHER ONE
      IF (ABS(XNEW-H)/H.LT.0.00005) GOTO 11
      IF((XNEW-A)*(XNEW-H).GT..0) GOTO 3
      NB=NH
      B=H
      GOTO 1
    3 NA=NH
      A=H
      GOTO 1
     
C***  CASES OF EXACT COINCIDENCE
    9 K=-1
      RETURN
   10 K=-N
      RETURN
   11 K=-NH
      IF (XNEW.LT.H) X(NH)=XNEW
      RETURN
     
      ENTRY SEQUINE  (NDIM,N,X,K,XNEW)
C***  ENTRY POINT IF THE INSERTION INDEX K IS GIVEN *********************
      IF (K.EQ.N+1) GOTO 8
      ! CALL ASSERT(K .gt. 0 .and. K .le. N,'ENTRY SEQUINE: WRONG INDEX GIVEN')
      IF (K.LE.0 .OR. K.GT.N) THEN
            WRITE(6,*) 'ENTRY SEQUINE: WRONG INDEX GIVEN'
            STOP 'ERROR'
            ENDIF
      NB=K
     
C***  INSERTION OF XNEW AT INDEX NB
    2 IF(N+1 .GT. NDIM) THEN
            WRITE(6,*) 'ARRAY DIMENSION TOO SMALL'
            STOP 'ERROR'
            ENDIF
      K=NB
    7 DO 5 J=K,N
      I=N+K-J
    5 X(I+1)=X(I)
    8 X(K)=XNEW
      N=N+1
      RETURN
     
C***  XNEW LIES OUTSIDE X(1) ... X(N)
    4 IF(N+1 .GT. NDIM) THEN
            WRITE(6,*) 'ARRAY DIMENSION TOO SMALL'
            STOP 'ERROR'
            ENDIF
      IF((XNEW-A)*(B-A).GT..0) GOTO 6
     
      GOTO 7
    6 N=N+1
      K=N
      X(N)=XNEW
      RETURN
      END subroutine
      end module
