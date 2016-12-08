      module MOD_READPOP

      contains

      subroutine readpop(ifl, T, popnum, pop1, pop2, pop3, rne, n, nd, modhead, jobnum)

      USE MOD_READMS
      USE MOD_READMSI
      USE MOD_ERROR

      IMPLICIT REAL*8(A-H,O-Z)

      real*8    T(ND), RNE(ND), POPNUM(ND*N), POP1(ND*N), POP2(ND*N), POP3(ND*N)
      CHARACTER CNAME*7, CREAD*7, MODHEAD*104, MODREAD*104
      INTEGER   JOBREAD, NSAVE, NDREAD

      CNAME='MODHEAD'

      READ (ifl,'(A7)') cread

      if (cread.ne.cname) then
	   write (6,*) 'READMOD: KEYWORD MISMATCH'
	   write (6,*) cread,cname
         pause
         call error('READPOP: KEYWORD MISMATCH')
	endif
	READ (ifl,'(A104)') MODREAD

      CNAME = 'JOBNUM'

      CALL READMSI1(IFL,JOBread,CNAME,IERR)

      if (jobread.gt.jobnum) jobnum = jobread

      CNAME='N'

      CALL READMSI1(IFL,NSAVE,CNAME,IERR)
      if (nsave .ne. N) then
	   write (6,*) 'READPOP: atomic data file mismatch'
	   stop 'file-read error'
      endif

      CNAME='ND'

      CALL READMSI1(IFL,NDREAD,CNAME,IERR)

      if (nd .ne. NDREAD) then
	   write (6,*) 'READPOP dimension of ND larger than NDREAD'
	   pause
	   stop
      endif

      CNAME='T'
      CALL READMS (IFL,T,ND,CNAME,IERR)

      CNAME='RNE'
      CALL READMS (IFL,RNE,ND,CNAME,IERR)

 
      CNAME='POPNUM'
      CALL READMS(IFL, POPNUM, ND*N, CNAME, IERR)

      CNAME='POP1'
      CALL READMS(IFL, POP1,   ND*N, CNAME, IERR)

      CNAME='POP2'
      CALL READMS(IFL, POP2,   ND*N, CNAME, IERR)

      CNAME='POP3'
      CALL READMS(IFL, POP3,   ND*N, CNAME, IERR)

      return

      end subroutine

      end module
