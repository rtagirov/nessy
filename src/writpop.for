      module MOD_WRITPOP
      contains
      subroutine writpop (ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,
     $  modhead,jobnum)
      use MOD_WRITMS
      use MOD_WRITMSI
c      subroutine writpop (ifl,T,popnum,pop1,pop2,pop3,rne,xneclc,n,nd,
c234567890 234567890 234567890 234567890 234567890 234567890 234567890 2
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 T(:),RNE(:),POP1(*),POP2(*),POP3(*)
      real*8 POPNUM(*)
      CHARACTER CNAME*10, MODHEAD*104

      CNAME='MODHEAD'
      write (ifl,'(A10)')  cname
	write (ifl,'(A104)') MODHEAD 
      CNAME='JOBNUM'
      CALL WRITMSI1(IFL,JOBNUM,CNAME,-1,IERR)
      CNAME='N'
      CALL WRITMSI1(IFL,N,CNAME,-1,IERR)
      CNAME='ND'
      CALL WRITMSI1(IFL,ND,CNAME,-1,IERR)
      CNAME='T'
      CALL WRITMS (IFL,T,ND,CNAME,-1,IERR)
      CNAME='RNE'
	CALL WRITMS (IFL,RNE,ND,CNAME,-1,IERR)
c      CNAME='XNECLC'
c      CALL WRITMS (IFL,XNECLC,ND,CNAME,-1,IERR)
      CNAME='POPNUM'
      CALL WRITMS (IFL,POPNUM,ND*N,CNAME,-1,IERR)
      CNAME='POP1'
      CALL WRITMS (IFL,POP1,ND*N,CNAME,-1,IERR)
      CNAME='POP2'
      CALL WRITMS (IFL,POP2,ND*N,CNAME,-1,IERR)
      CNAME='POP3'
      CALL WRITMS (IFL,POP3,ND*N,CNAME,-1,IERR)
      return
      end subroutine
      end module