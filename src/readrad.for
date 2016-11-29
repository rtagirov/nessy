      module MOD_READRAD

      contains

      SUBROUTINE READRAD(NF,ND,DUMMY,XJCREA,XJC,
     $                   XJL,HTOT,GTOT,XTOT,ETOT,EMFLUX,TOTIN,TOTOUT,
     $                   NCHARG,EDDREA,EDDI,NOM,WCHARM,N,
     $                   lastind,EINST,MODHEAD,JOBNUM)

C******************************************************************************
C***  READ THE RADIATION FIELD
C***  XJC: CONTINUUM RADIATION FIELD
C***  XJL: LINE RADIATION FIELD
C***  WCHARM: The Approximated Lambda Operator, Added 13-Mar-2006(Micha)
C*** ATTENTION: Watch Out For The 72charcters/Line Limit!
C******************************************************************************

      USE MOD_ERROR
      USE MOD_READMS
      USE MOD_READMSI
      USE MOD_FORMATS
      USE UTILS

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER ( ONE = 1.D+0, TWO = 2.D+0 )
      Real*8,Dimension(ND,NF)  :: WCHARM
      DIMENSION XJCREA(ND,NF), XJC(ND), EDDREA(3,ND,NF), EDDI(3,ND), DUMMY(ND)
      Real*8  TOTIN,TOTOUT
      Real*8,Dimension(NF)::EMFLUX
      Real*8,Dimension(ND):: HTOT,GTOT,XTOT,ETOT
      Integer JOBRead, NFREA
      DIMENSION EINST(N,N),XJL(ND,lastind)
      DIMENSION NCHARG(N),NOM(N)
      CHARACTER CREAD*7, MODHEAD*104, MODREAD*104

      character (len = 8) :: cname8
      character (len = 7) :: cname7
      character (len = 6) :: cname6

C***  CONTINUUM RADIATION FIELD XJC  *****************************************
C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS

      IFL = 2; open (IFL,File='RADIOC',STATUS='OLD', ACTION='READ')

      cname7='MODHEAD'

      READ (Ifl,'(A7)') cread

      If (cread .ne. cname7) Then
        Write (6,*) 'READRAD: KEYWORD MISMATCH: MODHEAD'
        Write (6,'("/",A,"/!=/",A,"/")') cread, cname7
        Pause
        Call Error('READRAD: KEYWORD MISMATCH: MODHEAD')
      Endif

      READ (Ifl,'(A104)') MODREAD

      CALL READMSI1(IFL,JOBread,'JOBNUM',IERR)

      If (Jobread.Lt.Jobnum-5) Then

        WRITE(6,'(A,I20,A,I20)') 
     $     'Readrad: Mismatch Between Job-Numbers',
     $     Jobread,'=Jobread <Jobnum-5=-5+', Jobnum
        Pause
        Call Error('Readrad: Mismatch Between Job-Numbers')
      Endif
      If (Jobread.Gt.Jobnum) Then
C         Print *,' Jobnum Updated From File RADIOC'
         Jobnum=Jobread
      Endif
      CALL READMSI1(IFL,NFREA,'NF',IERR)
      IF (NFREA.NE.NF) THEN
        WRITE (6,*) 'NF .Ne. NF On File RADIOC'
        PAUSE
        STOP
      ENDIF
      CALL READMS1 (IFL,TOTIN ,'TOTIN' ,IERR)
      CALL READMS1 (IFL,TOTOUT,'TOTOUT',IERR)
      CALL READMS (IFL,EMFLUX,NF,'EMFLUX',IERR)
      CALL READMS (IFL,HTOT,ND,'HTOT',IERR)
      CALL READMS (IFL,GTOT,ND,'GTOT',IERR)
      CALL READMS (IFL,XTOT,ND,'XTOT',IERR)
      CALL READMS (IFL,ETOT,ND,'ETOT',IERR)

      DO K = 1, NF

         write(cname8, FMT_KEY) 'XJC ', K

         CALL READMS(IFL,XJC,ND,cname8,IERR)

         DO L = 1, ND; XJCREA(L,K) = XJC(L); ENDDO

         write(cname8, FMT_KEY) 'EDDI',K

         CALL READMS(IFL,EDDI,3*ND,cname8,IERR)

         Do Ie=1,3
            Do L=1,Nd
               Eddrea(Ie,L,K)=Eddi(Ie,L)
            Enddo
         Enddo

         Write(cname8, FMT_KEY) 'WCHA',K

         CALL READMS(IFL,WCHARM(:,K),ND,cname8,IERR)

      ENDDO

      Close (IFL)

C     LINE RADIATION FIELD XJL
C     (DEPTH VEKTOR FOR EACH LINE TRANSITION LABELLED WITH IND)

      IFL = 4; open(IFL,File='RADIOL',STATUS='UNKNOWN', ACTION='READ')

      cname7='MODHEAD'

      READ (Ifl,'(A7)') cread

      If (cread .ne. cname7) Then
          Write (6,*) 'READMOD: KEYWORD MISMATCH'
          Write (6,'("/",A,"/ != /",A,"/")') cread, cname7
          Call Error('READRAD: KEYWORD MISMATCH',P=.true.)
      Endif

      READ (Ifl,'(A104)') MODREAD
      call assert(MODHEAD == MODREAD,
     &                       'Mismatch Between RADIO And MODEL Files')
      cname6='JOBNUM'
      CALL READMSI1(IFL,JOBread,cname6,IERR)
      call assert(JOBread >= Jobnum-5,'Mismatch Between Job-Numbers')
      if(JOBRead>Jobnum) Jobnum=Jobread

      IND = 0

      do 9 J = 2, N

         JM = J - 1

         do I = 1, JM

            IF ((NOM(I) .NE. NOM(J)) .OR. (NCHARG(I) .NE. NCHARG(J))) GOTO 99

            IND = IND + 1

            IF (EINST(I, J) .EQ. -TWO) GOTO 99

            write(cname8, FMT_KEY) 'XJL ', IND

            CALL READMS(IFL,DUMMY,ND,cname8,IERR)

            do L = 1, ND; XJL(L,IND) = DUMMY(L); enddo

   99    enddo

   9  enddo

      if (ind .gt. lastind) stop 'lastind Smaller Than IND'

      close(IFL)

      return

      end subroutine

      end module
