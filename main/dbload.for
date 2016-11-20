**********  MODULNAME: DBLOAD    ******* 3/07/92  22.14.04.******    11 KARTEN
      MODULE MOD_DBLOAD
      contains
      SUBROUTINE DBLOAD (A1,L,NDIM)
C***  LOAD ARRAY IN A1 FORT.19
      USE MOD_ERROR
      use MOD_FORMATS

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION A1(NDIM,NDIM)
      character namer*8, name*8

      write (name,FMT_KEY) 'BROY',L

      read (19,'(A8)') namer
	if (name .ne.namer) then
	   write (6,*) ' Broyden matrix read error'
	   call ERROR('dbload: error: Broyden Matrix Read ')
	endif

c      CALL READMS (19,A1,NDIM*NDIM,NAMER,IERR)
	read (19,*,err=666) A1
      if (any(isnan(A1)))
     $   call error('dbload.for: NaNs encounterd in BROYDEN')
	return

666   continue
c      IF (IERR .LT. 0 .AND. IERR .NE. -10) THEN
      write (6,*) 'ERROR WHEN READING RATCO FROM FILE 19'
      STOP 'ERROR'

      RETURN
      END SUBROUTINE
      END MODULE
