      module MOD_WRITRADL

      contains

      subroutine writradl(XJL,XJLMEAN,EINST,NCHARG,NOM,ND,N,LASTIND,MODHEAD,JOBNUM)

      use MOD_WRITMS
      use MOD_WRITMSI
      use MOD_FORMATS

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER (ONE = 1.D+0, TWO = 2.D+0)

      DIMENSION XJL(ND,LASTIND),XJLMEAN(ND),EINST(N,N)
      DIMENSION NCHARG(N),NOM(N)

      CHARACTER MODHEAD*104

      character (len = 8) :: cname8
      character (len = 6) :: cname6

!     LINE RADIATION FIELD XJL
!     (DEPTH VEKTOR FOR EACH LINE TRANSITION LABELLED WITH IND)
      IFL = 4; open (IFL, file='RADIOL', STATUS='UNKNOWN')

      write (ifl,'(A7)')  'MODHEAD'
      write (ifl,'(A104)') MODHEAD 

      cname6='JOBNUM'; CALL WRITMSI1(IFL,jobnum,cname6,-1,IERR)

      IND = 0

      DO J = 2, N

         JM = J - 1

         DO I = 1, JM

            IF ((NOM(I) .NE. NOM(J)) .OR. (NCHARG(I) .NE. NCHARG(J))) cycle

            IND = IND + 1

            IF (EINST(I, J) .EQ. -2.0d0) cycle

            XJLMEAN(1 : ND) = XJL(1 : ND, IND)

            WRITE(cname8,FMT_KEY) 'XJL ',IND

            CALL WRITMS (IFL,XJLMEAN,ND,cname8,-1,IERR)

        ENDDO

      ENDDO

      close(IFL)

      return

      end subroutine

      end module
