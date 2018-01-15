      module MOD_READMOD

      contains

      SUBROUTINE READMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                   GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                   ABXYZ,NATOM,MODHEAD,JOBNUM,LBLANK)

      USE MOD_READMS
      USE MOD_READMSI

      use utils

      IMPLICIT REAL*8(A-H,O-Z)

      REAL*8 TEFF,RSTAR, VDOP, ALMIN, ALMAX
      integer jobnum, NSAVE, ND, NF, NP, IPMAX, NBMAX,NBINW
      integer,parameter :: IPDIM=25,NBDIM=99
      integer,intent(in) :: NATOM

      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)

      dimension radius(nd),p(np),z(nd*np),entot(nd),
     $          velo(nd),gradi(nd),xlambda(nf),fweight(nf),
     $          akey(nf),abxyz(natom)

      CHARACTER CREAD*7, MODHEAD*104

      READ(ifl,'(A7)') cread

      if (trim(adjustl(cread)) .ne. trim(adjustl('MODHEAD'))) then

          write (6, *) 'READMOD: KEYWORD MISMATCH'
          write (6, '(A7,2x,A7)') trim(adjustl(cread)), trim(adjustl('MODHEAD'))

          call error('READMOD: KEYWORD MISMATCH')

      endif

      READ (ifl, '(A104)') MODHEAD

      CALL READMSI1(IFL,jobnum,'JOBNUM',IERR)
      CALL READMS1 (IFL,TEFF,'TEFF',IERR)
      CALL READMS1 (IFL,RSTAR,'RSTAR',IERR)
      CALL READMSI1(IFL,NSAVE,'N',IERR)

	if (nsave.ne.N) then
	   write (6,*) ' atomic data file mismatch'
	   stop 'file-read error'
	endif

      CALL READMS  (IFL, ABXYZ,   NATOM, 'ABXYZ',   IERR)
      CALL READMSI1(IFL, ND,             'ND',      IERR)
      CALL READMS  (IFL, RADIUS,  ND,    'R',       IERR)
      CALL READMS  (IFL, ENTOT,   ND,    'ENTOT',   IERR)
      CALL READMS  (IFL, VELO,    ND,    'VELO',    IERR)
      CALL READMS  (IFL, GRADI,   ND,    'GRADI',   IERR)
      CALL READMS1 (IFL, VDOP,           'VDOP',    IERR)
      CALL READMSI1(IFL, NF,             'NF',      IERR)
      CALL READMS  (IFL, XLAMBDA, NF,    'XLAMBDA', IERR)
      CALL READMS  (IFL, FWEIGHT, NF,    'FWEIGHT', IERR)
      CALL READMS  (IFL, AKEY,    NF,    'KEY',     IERR)
      CALL READMSI1(IFL, NP,             'NP',      IERR)
      CALL READMS  (IFL, P,       NP,    'P',       IERR)
      CALL READMS  (IFL, Z,       ND*NP, 'Z',       IERR)
 
C***  LINE BLANKETING TABLE 
      IF (LBLANK.gt.0) THEN
	   IF (IPMAX.GT.IPDIM .OR. NBMAX.GT.NBDIM) THEN
	      WRITE (6,*) 'IP OR NB DIMENSION ERROR  '
	      PAUSE
	      STOP
	      ENDIF
	   print *,' LB-table read from MODFILE' 
         IERR=0
         CALL READMSI1(IFL,IPMAX,'IPMAX',IERR)
         CALL READMSI1(IFL,NBMAX,'NBMAX',IERR)
         CALL READMSI1(IFL,NBINW,'NBINW',IERR)
         CALL READMS1 (IFL,ALMIN,'ALMIN',IERR)
         CALL READMS1 (IFL,ALMAX,'ALMAX',IERR)
         CALL READMS  (IFL,SCAGRI,IPDIM*NBDIM,'SCAGRI',IERR)
         CALL READMS  (IFL,SCAEVT,IPDIM*NBDIM,'SCAEVT',IERR)
         CALL READMS  (IFL,ABSEVT,IPDIM*NBDIM,'ABSEVT',IERR)
         LBLAON=1
      ENDIF
 
      return

      end subroutine

      end module
