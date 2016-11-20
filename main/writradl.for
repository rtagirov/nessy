      module MOD_WRITRADL
      contains
      subroutine writradl (XJL,XJLMEAN,EINST,NCHARG,NOM,nd,N,NDIM,
     $                     LASTIND,MODHEAD,JOBNUM)
      use MOD_WRITMS
      use MOD_WRITMSI
      use MOD_FORMATS
      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER ( ONE = 1.D+0, TWO = 2.D+0 )

      DIMENSION XJL(ND,LASTIND),XJLMEAN(ND),EINST(NDIM,NDIM)
      DIMENSION NCHARG(NDIM),NOM(N)

c      dimension xjcARR(nd,nf),xjc(nd),eddARR(3,nd,nf),eddi(3,nd)
      CHARACTER CNAME*10, MODHEAD*104

C***  LINE RADIATION FIELD XJL  ***********************************************
C***  ( DEPTH VEKTOR FOR EACH LINE TRANSITION LABELLED WITH IND )
      IFL=4
      open (IFL,file='RADIOL',STATUS='UNKNOWN')

      CNAME='MODHEAD'
      write (ifl,'(A10)')  cname
	write (ifl,'(A104)') MODHEAD 
      CNAME='JOBNUM'
      CALL WRITMSI1(IFL,jobnum,CNAME,-1,IERR)

      IND=0
      DO J=2,N
        JM=J-1
        DO I=1,JM
          IF ((NOM(I) .NE. NOM(J)) .OR. (NCHARG(I) .NE. NCHARG(J)))cycle
          IND=IND+1
          IF (EINST(I,J) .EQ. -2.0d0) cycle
          XJLMEAN(1:nd)=XJL(1:nd,IND)
          WRITE(CNAME,FMT_KEY) 'XJL ',IND
          CALL WRITMS (IFL,XJLMEAN,ND,CNAME,-1,IERR)
        ENDDO
      ENDDO
      close (IFL)

	return
	end subroutine
      end module
