      module MOD_REBLANK
      contains
      SUBROUTINE REBLANK (LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)
C  THIS SUBROUTINE READS THE DATA FOR SUBROUTINE LINSCA
C     ==> IP AND DIMENSIONS HAVE TO BE UPDATED, IF THE TABLE IS CHANGED
      use MOD_LINSCA
      IMPLICIT REAL*8(A-H,O-Z)
 
      parameter (IPDIM=25,NBDIM=99)
      DIMENSION XLAMBDA(NF),ENTOT(ND),RNE(ND)
      COMMON /LIBLPAR/ ALMIN, ALMAX, LBLAON, IPMAX, NBMAX, NBINW
      COMMON /LIBLDAT/ SCAGRI(IPDIM), SCAEVT(IPDIM,NBDIM), 
     $                                ABSEVT(IPDIM,NBDIM)
      DIMENSION SCAFAC(ND,NF),ABSFAC(ND,NF)
      LBLAON=0
      IF (LBLANK.LE.-1) THEN
	   print 998
  998    FORMAT(9X,'REBLANK: LB TABLE READ FROM FILE LIBLANK')
         open (7,file='LIBLANK',status='OLD')
         READ (7,*,ERR=100) IPMAX
         IF (IPMAX.GT.IPDIM) STOP ' IPDIM-REBLANK'
         READ (7,*,ERR=101) (SCAGRI(I),I=1,IPMAX)
         READ (7,*,ERR=102) BMAX
         NBMAX=nint(BMAX)
c         IF (NBMAX.GT.NBDIM) STOP ' NBDIM-REBLANK'
         IF (NBMAX.GT.NBDIM) then
            print '(A,A,A)',' Blanketing Table input: ',
     &           'There are more bins than allowed by the dimension ',
     &           'of the routine'
            NBMAX=NBDIM
         endif
         READ (7,*,ERR=102) NBINW
         PRINT *,'    BIN WIDTH: ',NBINW,' A'
         IF (NBINW.ne.50 .and. NBINW.NE.20) then
            print *,' NBINW=',NBINW,' STOP'
            STOP ' not new format'
         endif
         DO 10 NB=1,NBMAX
         READ (7,*,ERR=103) (SCAEVT(I,NB),I=1,IPMAX)
         READ (7,*,ERR=104) (ABSEVT(I,NB),I=1,IPMAX)
c* new files have only 2-entries   READ (7,*,ERR=105) (DUMMY    ,I=1,IPMAX)
C*** For test only - no absorptions 
C         DO I=1,IPMAX
C            ABSEVT(I,NB)=1.
C         ENDDO
   10    ENDDO
c*** use only the first 20 bins ............ 8-tung: manipulation
         IPMAX=IPMAX-1
         ALMIN=227.84
         ALMAX=NBINW*NBMAX - 0.000001
c         NBMAX=21
C         NBMAX=24
         PRINT 999,IPMAX,NBMAX,ALMIN,ALMAX,NBINW
c  999    FORMAT(/9X,'LINE-BLANKETING PARAMETERS'/9X,
c     $     26('=')/,2X,'IP=',I5,5X,'NB=',I5,9X,'LAMBDA-MIN',F10.2,5X,
c     $     'LAMBDA-MAX',F10.2,5x,' BIN-WIDTH',I5,/)
999      FORMAT(9X,'LB PARAMETERS:  IP=',I5,3X,'NB=',I5,9X,'LAMBDA-MIN'
     $         ,F10.2,3X,'LAMBDA-MAX',F10.2,3x,'BIN-WIDTH',I5)

C         IF (LBLANK.EQ.-1) THEN
C         IERR=0
C         CALL WRITMS(3,IPMAX,1,5HIPMAX,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,NBMAX,1,5HNBMAX,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,NBINW,1,5HNBINW,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,ALMIN,1,5HALMIN,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,ALMAX,1,5HALMAX,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,SCAGRI,IPDIM*NBDIM,6HSCAGRI,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,SCAEVT,IPDIM*NBDIM,6HSCAEVT,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         CALL WRITMS(3,ABSEVT,IPDIM*NBDIM,6HABSEVT,-1,IDUMMY,IERR)
C         IF (IERR.LT.0) GOTO 106
C         PRINT *,' LINE BLANKETING TABLE WRITTEN TO MODEL FILE'
C         ENDIF

         LBLAON=1
c         LBLANK=ABS(LBLANK)
         CLOSE (7)
      ELSE IF (LBLANK.GE.1) THEN
         LBLAON=1
         PRINT 999,IPMAX,NBMAX,ALMIN,ALMAX,NBINW
c***  the table is now read in routine READMOD
c	   CNAME='IPMAX'
c         CALL READMSI(3,IPMAX,1,cname,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='NBMAX'
c         CALL READMSI(3,NBMAX,1,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='NBINW'
c         CALL READMSI(3,NBINW,1,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='ALMIN'
c         CALL READMS (3,ALMIN,1,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='ALMAX'
c         CALL READMS (3,ALMAX,1,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='SCAGRI'
c         CALL READMS (3,SCAGRI,IPDIM*NBDIM,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='SCAEVT'
c         CALL READMS (3,SCAEVT,IPDIM*NBDIM,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c	   CNAME='ABSEVT'
c         CALL READMS (3,ABSEVT,IPDIM*NBDIM,CNAME,IERR)
c         IF (IERR .LT. 0) GOTO 107
c         LBLAON=1
c         PRINT *,' LINE BLANKETING TABLE READ FROM MODEL FILE'
c         PRINT *,LBLAON, ALMIN, ALMAX, IPMAX, NBMAX, NBINW
cC         PRINT *,SCAEVT
cC         PRINT *,ABSEVT
      ELSE
         LBLAON=0
         ALMIN=0.
         ALMAX=0.

         PRINT *,' REBLANK: LINE BLANKETING TABLE NOT ACTIVE'
      ENDIF



c***  write the factors into the arreys
      DO 4 K=1,NF


      IF (XLAMBDA(K).GT.ALMAX) GOTO 4

      DO L=1,ND

      DENS = ENTOT(L)*RNE(L)
      CALL LINSCA (XLAMBDA(K),DENS,SCAFAC(L,K),ABSFAC(L,K))

      ENDDO


    4 CONTINUE


      RETURN
  100 PRINT *,' ERROR READING IPMAX FROM FILE LIBLANK'
      STOP ' REBLANK'
  101 PRINT *,' ERROR READING GRID FROM FILE LIBLANK'
      STOP ' REBLANK'
  102 PRINT *,' ERROR READING NBMAX FROM FILE LIBLANK'
      STOP ' REBLANK'
  103 PRINT *,' ERROR READING SCAT EVENTS FROM FILE LIBLANK'
      PRINT *,' NB= ',NB
      STOP ' REBLANK'
  104 PRINT *,' ERROR READING ABS EVENTS FROM FILE LIBLANK'
      PRINT *,' NB= ',NB
      STOP ' REBLANK'
  105 PRINT *,' ERROR DUMMY READ FROM FILE LIBLANK'
      PRINT *,' NB= ',NB
      STOP ' REBLANK'
  106 PRINT *,' ERROR WRITMS TO MODEL FILE'
      PRINT *,' IERR = ',IERR
      STOP ' REBLANK'
  107 PRINT *,' ERROR READMS FROM MODEL FILE'
      PRINT *,' IERR = ',IERR
      STOP ' REBLANK'
      END subroutine
      end module
