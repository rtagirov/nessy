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
!      IMPLICIT NONE
      PARAMETER ( ONE = 1.D+0, TWO = 2.D+0 )
      Real*8,Dimension(ND,NF)  :: WCHARM
      DIMENSION XJCREA(ND,NF),XJC(ND)
     $          ,EDDREA(3,ND,NF),EDDI(3,ND),DUMMY(ND)
      Real*8  TOTIN,TOTOUT
      Real*8,Dimension(NF)::EMFLUX
      Real*8,Dimension(ND):: HTOT,GTOT,XTOT,ETOT
      Integer JOBRead, NFREA
      DIMENSION EINST(N,N),XJL(ND,lastind)
      DIMENSION NCHARG(N),NOM(N)
      CHARACTER CREAD*10, MODHEAD*104, MODREAD*104
      Character CNAME*10
C***  CONTINUUM RADIATION FIELD XJC  *****************************************
C***  LOOP OVER ALL CONTINUUM FREQUENCY POINTS
      IFL=2
      Open (IFL,File='RADIOC',STATUS='OLD', ACTION='READ')

      CNAME='MODHEAD'
C      CALL READMSC(IFL,MODREAD,104,CNAME,IERR)
      READ (Ifl,'(A10)') Cread
      If (Cread.Ne.Cname) Then
        Write (6,*) 'READRAD: KEYWORD MISMATCH: MODHEAD'
        Write (6,'("/",A,"/!=/",A,"/")') Cread,Cname
        Pause
        Call Error('READRAD: KEYWORD MISMATCH: MODHEAD')
      Endif
      READ (Ifl,'(A104)') MODREAD
 !     IF (MODHEAD .NE. MODREAD) Then
 !       WRITE (6,*) 'Mismatch Between RADIO And MODEL Files'
 !       PRINT *, MODHEAD, ' .NE. ', MODREAD
 !       CALL ERROR('Missmatch Between RADIO And MODEL Files')
 !       Stop
 !     Endif
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

      DO K=1,NF
         Write(CNAME,FMT_KEY) 'XJC ',K
         CALL READMS (IFL,XJC,ND,CNAME,IERR)
         DO L=1,ND
            XJCREA(L,K)=XJC(L)
         ENDDO
C***  
         Write(CNAME,FMT_KEY) 'EDDI',K
         CALL READMS (IFL,EDDI,3*ND,CNAME,IERR)
         Do Ie=1,3
            Do L=1,Nd
               Eddrea(Ie,L,K)=Eddi(Ie,L)
            Enddo
         Enddo
         Write(CNAME,FMT_KEY) 'WCHA',K
         CALL READMS(IFL,WCHARM(:,K),ND,CNAME,IERR)
      ENDDO
      Close (IFL)
C******************************************************************************
C***  LINE RADIATION FIELD XJL  ***********************************************
C***  ( DEPTH VEKTOR FOR EACH LINE TRANSITION LABELLED WITH IND )
      IFL=4
      Open (IFL,File='RADIOL',STATUS='UNKNOWN', ACTION='READ')
      CNAME='MODHEAD'
C      CALL READMSC(IFL,MODREAD,104,CNAME,IERR)
      READ (Ifl,'(A10)') Cread
      If (Cread.Ne.Cname) Then
         Write (6,*) 'READMOD: KEYWORD MISMATCH'
         Write (6,'("/",A,"/ != /",A,"/")') Cread,Cname
         Call Error('READRAD: KEYWORD MISMATCH',P=.true.)
      Endif
      READ (Ifl,'(A104)') MODREAD
      call assert(MODHEAD == MODREAD,
     &                       'Mismatch Between RADIO And MODEL Files')
      CNAME='JOBNUM'
      CALL READMSI1(IFL,JOBread,CNAME,IERR)
      call assert(JOBread >= Jobnum-5,'Mismatch Between Job-Numbers')
      if(JOBRead>Jobnum) Jobnum=Jobread
      IND=0
      DO 9 J=2,N
      JM=J-1
      DO I=1,JM
      IF ((NOM(I) .NE. NOM(J)) .OR. (NCHARG(I) .NE. NCHARG(J))) GOTO 99
      IND=IND+1
      IF (EINST(I,J) .EQ. - TWO) GOTO 99
      Write(CNAME,FMT_KEY) 'XJL ',IND
      CALL READMS (IFL,DUMMY,ND,CNAME,IERR)
      Do L=1,Nd
         XJL(L,IND)=DUMMY(L)
      Enddo
   99 Enddo
   9  Enddo
      If (Ind.Gt.lastind) Then
         Write (6,*) 'lastind Smaller Than IND'
         Stop
      Endif
      Close (IFL)
C*****************************************************************************
      RETURN
      END Subroutine
      End Module
