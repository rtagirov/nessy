      MODULE BROYDEN

      contains

      SUBROUTINE DBCLOSE(ND, nrank, A1)
C***  copy file 18 to file.19

      IMPLICIT REAL*8(A-H,O-Z)

      dimension A1(nrank,nrank)
      CHARACTER MODHEAD*104, name*8
      rewind 18
	rewind 19

	read (18,'(a104)',end=1) modhead
	write (19,'(a104)') modhead

      do L=1,nd

         read (18,'(A8)') name
         write (19,'(A8)') name
         read (18,*) A1
         write (19,*) A1
	enddo

1     continue

      close (18)
      close (19)

      RETURN

      END SUBROUTINE

      SUBROUTINE DBOPEN(MODHEAD)
!     OPEN FORT.19

      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)

      CHARACTER MODHEAD*104, MODEL*104

      !***  open the file for the updated Broyden matrix
      open (18,file='NEWBROYDEN',status='unknown',err=888)
      write (18,'(A104)') modhead

      open (19,file='BROYDEN',status='old',err=666)
      NOFILE(1 : DPN) = .FALSE.

!     CHECK WHETHER THESE DATA EXIST AND BELONG TO THE PRESENT MODEL

      read(19, '(A104)', err=777) model

      IF (MODEL .NE. MODHEAD) THEN
!        FILE IS NOT APPROPRIATE - CLEAR THE INDEX ARRAY
         write (6, *) ' BROYDEN file does not belong to this model'
         NOFILE(1 : DPN)= .TRUE.
      ENDIF

      RETURN

666   continue
      NOFILE(1 : DPN) = .TRUE.
      open (19, file='BROYDEN', status='NEW')

      RETURN

777   continue
      write (6,*) 'ERROR WHEN READING MODHEAD FROM TAPE 19'
      STOP 'ERROR ft19'

888   continue
      write (6,*) ' new Broyden file is already existing'
      stop 'error ft18'

      END SUBROUTINE

      SUBROUTINE DBSAVE(A1, L, N)
!     SAVE ARRAY IN A1 FORT.19

      use utils
      USE MOD_FORMATS

      IMPLICIT REAL*8(A - H, O - Z)

      DIMENSION A1(N, N)

      CHARACTER NAME*8

      WRITE(NAME, FMT_KEY) 'BROY', L

      WRITE(18, '(A8)') name
      WRITE(18, *)      A1

      IF (ANY(ISNAN(A1))) CALL ERROR('DBSAVE: NaNs ENCOUNTERED IN BROYDEN')

      RETURN

      END SUBROUTINE

      SUBROUTINE DBLOAD(A1, L, N)
C***  LOAD ARRAY IN A1 FORT.19
      use utils
      use MOD_FORMATS

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION A1(N, N)
      character namer*8, name*8

      write (name,FMT_KEY) 'BROY',L

      read (19,'(A8)') namer
	if (name .ne.namer) then
	   write (6,*) ' Broyden matrix read error'
	   call ERROR('dbload: error: Broyden Matrix Read ')
	endif

	read (19,*,err=666) A1
      if (any(isnan(A1)))
     $   call error('dbload.for: NaNs encounterd in BROYDEN')
	return

666   continue

      write (6,*) 'ERROR WHEN READING RATCO FROM FILE 19'
      STOP 'ERROR'

      RETURN

      END SUBROUTINE

      END MODULE
