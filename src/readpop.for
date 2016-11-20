      module MOD_READPOP
      contains
      subroutine readpop (ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,
     $                    modhead,jobnum)
c     subroutine readpop (ifl,T,popnum,pop1,pop2,pop3,rne,xneclc,n,nd,
c234567890 234567890 234567890 234567890 234567890 234567890 234567890 2
      USE MOD_READMS
      USE MOD_READMSI
      USE MOD_ERROR
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 T(ND), RNE(ND), POPNUM(ND*N),POP1(ND*N), POP2(ND*N), POP3(ND*N)
	CHARACTER CNAME*10, CREAD*10, MODHEAD*104, MODREAD*104
      INTEGER JOBREAD,NSAVE,NDREAD
c	PRINT *,'1. IN READPOP: N=',N, 'ifl=',ifl
c	PAUSE
      CNAME='MODHEAD'
c      CALL READMSC(IFL,MODREAD,104,CNAME,IERR)
      READ (ifl,'(A10)') cread
c	PRINT *,'2. IN READPOP, AFTER CREAD: N=',N, 'ifl=',ifl
c	PAUSE
      if (cread.ne.cname) then
	   write (6,*) 'READMOD: KEYWORD MISMATCH'
	   write (6,*) cread,cname
         pause
         call error('READPOP: KEYWORD MISMATCH')
	endif
	READ (ifl,'(A104)') MODREAD
c	PRINT *,'2. IN READPOP AFTER MODREAD: N=',N, 'ifl=',ifl
c	PAUSE
      IF (MODHEAD .NE. MODREAD) then
!	   WRITE (6,*) 'mismatch between POPNUM and MODEL files'
!	   stop
	endif
      CNAME='JOBNUM'
      CALL READMSI1(IFL,JOBread,CNAME,IERR)
      if (jobread.gt.jobnum) then
c	   WRITE (6,*) 'readpop: mismatch between job-numbers'
c	   pause
         jobnum=jobread
c	   stop
	endif
      CNAME='N'
      CALL READMSI1(IFL,NSAVE,CNAME,IERR)
	if (nsave.ne.N) then
	   write (6,*) 'READPOP: atomic data file mismatch'
	   stop 'file-read error'
	endif
      CNAME='ND'
      CALL READMSI1(IFL,NDREAD,CNAME,IERR)
	if (nd.ne.NDREAD) then
	   write (6,*) 'READPOP dimension of ND larger than NDREAD'
	   pause
	   stop
	endif
      CNAME='T'
      CALL READMS (IFL,T,ND,CNAME,IERR)
      CNAME='RNE'
      CALL READMS (IFL,RNE,ND,CNAME,IERR)
c	PRINT *,'3. IN READPOP BEFORE XNECLC: N=',N, 'ifl=',ifl
c	PAUSE
c	CNAME='XNECLC'
c	PRINT *,'3. IN READPOP AFTER XNECLC: N=',N, 'ifl=',ifl
c	PAUSE
c      CALL READMS (IFL,XNECLC,ND,CNAME,IERR)
c	print *,'N=',N, 'ifl=',ifl
 
      CNAME='POPNUM'
c	print *,'popnum N=',N, 'ifl=',ifl
c	pause
      CALL READMS (IFL,POPNUM,ND*N,CNAME,IERR)
      CNAME='POP1'
c	print *,'POPNUM N=',N, 'ifl=',ifl
c	pause
      CALL READMS (IFL,POP1,ND*N,CNAME,IERR)
      CNAME='POP2'
c      print *,'POP2 N=',N, 'ifl=',ifl
c	pause
	CALL READMS (IFL,POP2,ND*N,CNAME,IERR)
      CNAME='POP3'
c      print *,'POP3 N=',N, 'ifl=',ifl
c	pause
	CALL READMS (IFL,POP3,ND*N,CNAME,IERR)

      return
	end subroutine
      end module
