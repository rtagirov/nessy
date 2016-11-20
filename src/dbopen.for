      MODULE MOD_DBOPEN

      CONTAINS

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

      END MODULE
