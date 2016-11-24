      module MOD_WRITPOP

      contains

      subroutine writpop(ifl,T,popnum,pop1,pop2,pop3,rne,n,nd,modhead,jobnum)

      use MOD_WRITMS
      use MOD_WRITMSI

      IMPLICIT REAL*8(A-H,O-Z)

      real*8 T(:),RNE(:),POP1(*),POP2(*),POP3(*)
      real*8 POPNUM(*)

      CHARACTER MODHEAD*104

      write (ifl,'(A7)') 'MODHEAD'
      write (ifl,'(A104)') MODHEAD 

      CALL WRITMSI1(IFL,JOBNUM,'JOBNUM',-1,IERR)

      CALL WRITMSI1(IFL,N,'N',-1,IERR)

      CALL WRITMSI1(IFL,ND,'ND',-1,IERR)

      CALL WRITMS (IFL,T,ND,'T',-1,IERR)

      CALL WRITMS (IFL,RNE,ND,'RNE',-1,IERR)

      CALL WRITMS (IFL,POPNUM,ND*N,'POPNUM',-1,IERR)

      CALL WRITMS (IFL,POP1,ND*N,'POP1',-1,IERR)

      CALL WRITMS (IFL,POP2,ND*N,'POP2',-1,IERR)

      CALL WRITMS (IFL,POP3,ND*N,'POP3',-1,IERR)

      return

      end subroutine

      end module
