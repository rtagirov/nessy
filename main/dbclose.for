**********  MODULNAME: DBCLOSE    ******* 3/07/92  22.14.04.******    11 KARTEN
      MODULE MOD_DBCLOSE
      contains
      SUBROUTINE DBCLOSE (ND,nrank,A1)
C***  copy file 18 to file.19

      IMPLICIT REAL*8(A-H,O-Z)

      dimension A1(nrank,nrank)
      CHARACTER MODHEAD*104, name*8
      rewind 18
	rewind 19

	read (18,'(a104)',end=1) modhead
	write (19,'(a104)') modhead

      do L=1,nd

c      PRINT*, L, 'DEPTH CYCLE IN DBCLOSE.FOR IS GOING ON'

c	   print *,' DBCLOSE L= ',L,NRANK
c	   print *,'A1',A1
         read (18,'(A8)') name
         write (19,'(A8)') name
         read (18,*) A1
         write (19,*) A1
	enddo

c      CALL CLOSMS (19,IERR)
1     continue
      close (18)
	close (19)

      RETURN
      END SUBROUTINE
      END MODULE
