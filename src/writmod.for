      module MOD_WRITMOD

      contains

      SUBROUTINE WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,
     $                   GRADI,RSTAR,VDOP,NF,XLAMBDA,FWEIGHT,AKEY,
     $                   ABXYZ,NATOM,MODHEAD,JOBNUM)

      use MOD_WRITMS
      use MOD_WRITMSI

      IMPLICIT REAL*8(A-H,O-Z)
      parameter (IPDIM=25,NBDIM=99)

      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)

      CHARACTER CNAME*7, MODHEAD*104

      dimension radius(nd),p(np),z(nd*np),entot(nd),
     $          velo(nd),gradi(nd),xlambda(nf),fweight(nf),
     $          akey(nf),abxyz(natom)

      CNAME = 'MODHEAD'

      write(ifl, '(A7)')   cname
      write(ifl, '(A104)') MODHEAD 

      CALL WRITMSI1(IFL,JOBNUM,       'JOBNUM', -1,IERR)
      CALL WRITMS1 (IFL,TEFF,         'TEFF',   -1,IERR)
      CALL WRITMS1 (IFL,RSTAR,        'RSTAR',  -1,IERR)
      CALL WRITMSI1(IFL,N,            'N',      -1,IERR)
      CALL WRITMS  (IFL,ABXYZ,  NATOM,'ABXYZ',  -1,IERR)
      CALL WRITMSI1(IFL,ND,           'ND',     -1,IERR)
      CALL WRITMS  (IFL,RADIUS, ND,   'R',      -1,IERR)
      CALL WRITMS  (IFL,ENTOT,  ND,   'ENTOT',  -1,IERR)
      CALL WRITMS  (IFL,VELO,   ND,   'VELO',   -1,IERR)
      CALL WRITMS  (IFL,GRADI,  ND,   'GRADI',  -1,IERR)
      CALL WRITMS1 (IFL,VDOP,         'VDOP',   -1,IERR)
      CALL WRITMSI1(IFL,NF,           'NF',     -1,IERR)
      CALL WRITMS  (IFL,XLAMBDA,NF,   'XLAMBDA',-1,IERR)
      CALL WRITMS  (IFL,FWEIGHT,NF,   'FWEIGHT',-1,IERR)
      CALL WRITMS  (IFL,AKEY,   NF,   'KEY',    -1,IERR)

      CALL WRITMSI1(IFL,NP,       'NP', -1,IERR)
      CALL WRITMS  (IFL,P, NP,    'P',  -1,IERR)
      CALL WRITMS  (IFL,Z, ND*NP, 'Z',  -1,IERR)
 
C***  LINE BLANKETING TABLE 
      IF (LBLAON.GE.1) THEN
	   print *,' WRITMOD:  LB table written to MODFILE'
	   IF (IPMAX.GT.IPDIM .OR. NBMAX.GT.NBDIM) THEN
	      WRITE (6,*) 'IP OR NB DIMENSION ERROR'
	      PAUSE
	      STOP
	   ENDIF
         IERR=0
         CALL WRITMSI1(IFL,IPMAX,             'IPMAX', -1,IERR)
         CALL WRITMSI1(IFL,NBMAX,             'NBMAX', -1,IERR)
         CALL WRITMSI1(IFL,NBINW,             'NBINW', -1,IERR)
         CALL WRITMS1 (IFL,ALMIN,             'ALMIN', -1,IERR)
         CALL WRITMS1 (IFL,ALMAX,             'ALMAX', -1,IERR)
         CALL WRITMS  (IFL,SCAGRI,IPDIM*NBDIM,'SCAGRI',-1,IERR)
         CALL WRITMS  (IFL,SCAEVT,IPDIM*NBDIM,'SCAEVT',-1,IERR)
         CALL WRITMS  (IFL,ABSEVT,IPDIM*NBDIM,'ABSEVT',-1,IERR)
      ENDIF
 
      RETURN
	END subroutine
      end module
