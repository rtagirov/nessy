      module MOD_WRITRADC

      contains

      subroutine writradc(xjc2,xjc,eddarr,eddi,emflux,totin,totout,HTOT,GTOT,XTOT,ETOT,wcharm,nd,nf,MODHEAD,JOBNUM)

      use MOD_WRITMS
      use MOD_WRITMSI
      use MOD_FORMATS

      IMPLICIT REAL*8(A - H, O - Z)

      real*8, dimension(ND, NF) :: wcharm

      real*8 :: xjc2(nd, nf), xjc(nd), eddARR(3, nd, nf), eddi(3, nd)

      real*8 :: emflux(NF), HTOT(ND), GTOT(ND), XTOT(ND), ETOT(ND)

      CHARACTER CNAME*10, MODHEAD*104

!     CONTINUUM RADIATION FIELD XJC
!     LOOP OVER ALL CONTINUUM FREQUENCY POINTS
!     (DEPTH VEKTOR AT EACH FREQUENCY POINT)
!     CALLED BY COMO

      IFL = 2; open(IFL, file = 'RADIOC', STATUS = 'UNKNOWN')

      write(ifl, '(A)')   'MODHEAD'
      write(ifl, '(A104)') MODHEAD

      CALL WRITMSI1(IFL, JOBNUM, 'JOBNUM', -1, IERR)
      CALL WRITMSI1(IFL, NF,     'NF',     -1, IERR)

      CALL WRITMS1(IFL, TOTIN,  'TOTIN',  -1, IERR)
      CALL WRITMS1(IFL, TOTOUT, 'TOTOUT', -1, IERR)

      CALL WRITMS(IFL, EMFLUX, NF, 'EMFLUX', -1, IERR)
      CALL WRITMS(IFL, HTOT,   ND, 'HTOT',   -1, IERR)
      CALL WRITMS(IFL, GTOT,   ND, 'GTOT',   -1, IERR)
      CALL WRITMS(IFL, XTOT,   ND, 'XTOT',   -1, IERR)
      CALL WRITMS(IFL, ETOT,   ND, 'ETOT',   -1, IERR)

      DO 6 K = 1, NF

      WRITE(CNAME, FMT_KEY) 'XJC ', K

      xjc(1 : ND) = xjc2(1 : ND, K)

      eddi(1,1:ND)=eddARR(1,1:ND,K)
      eddi(2,1:ND)=eddARR(2,1:ND,K)
      eddi(3,1:ND)=eddARR(3,1:ND,K)

      CALL WRITMS (IFL,XJC,ND,CNAME,-1,IERR)

!     EDDI IS A DUMMY-WRITE
      WRITE(CNAME, FMT_KEY) 'EDDI', K
      CALL WRITMS(IFL,EDDI,3*ND,CNAME,-1,IERR)

      WRITE(CNAME, FMT_KEY) 'WCHA', K
      CALL WRITMS(IFL,WCHARM(:,K),ND,CNAME,-1,IERR)

    6 CONTINUE

      close (IFL)

      return

      end subroutine

      end module
