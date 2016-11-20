      module MOD_PLOTPOP
      contains
C**********  MODULNAME: PLOTPOP   ******* 24/03/87  21.22.07.******    54 KARTEN
      SUBROUTINE PLOTPOP (LEVELPL,N,ND,LEVEL,ENTOT,POPNUM,MODHEAD,
     $          JOBNUM )
C***  DIRECT TRANSFER OF POPNUMBER PLOT
      USE MOD_PLOTCON
      use MOD_PLOTANF
      implicit real*8(a-h,o-z)

      DIMENSION X(100),Y(100)
      DIMENSION LEVELPL(N),ENTOT(ND),POPNUM(ND,N)
      CHARACTER HEADER*60,MODHEAD*104,LEVEL(N)*10
           
      KANAL=2
      OPEN (KANAL,FILE='PLOT')
C      CALL JSYMSET (2LG2,'TRANSFER')
      write (6,*) 'POP PLOT DATA TO BE ROUTED'
      NPLOT=0
      XMIN=7.
      XMAX=16.
      YMIN=-15.
      YMAX=.0
      XSCALE=2.5
      YSCALE=1.
      XTICK=1.
      YTICK=1.
      XABST=3.
      YABST=5.
     
   10 NPLOT=NPLOT+1
      IF (LEVELPL(NPLOT) .EQ. 0) RETURN
      DECODE (8,5,LEVELPL(NPLOT)) L1,L2,L3
    5 FORMAT (3I2)
      IF (L1.LE.0 .OR. L1 .GT. N) GOTO 10
      ENCODE (60,4,HEADER) JOBNUM
    4 FORMAT (19X,1HJ,I3)
      HEADER (1:18)=MODHEAD(15:32)
      IF (L1 .GT. 0) HEADER(25:34)=LEVEL(L1)
      IF (L2 .GT. 0) HEADER(36:45)=LEVEL(L2)
      IF (L3 .GT. 0) HEADER(47:56)=LEVEL(L3)
      DO 1 L=1,ND
      X(L)=LOG10(ENTOT(L))
    1 Y(L)=LOG10(POPNUM(L,L1))
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,60H      LOG OF NUMBER DENSITY / (CM**-3)                                   
     $ ,60H      LOG OF REL. POP. NUMBERS                                           
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,X,Y,ND,1)
      IF (L2 .LE. 0 .OR. L2 .GT. N) GOTO 10
      DO 2 L=1,ND
    2 Y(L)=LOG10(POPNUM(L,L2))
      CALL PLOTCON (KANAL,X,Y,ND,2)
      IF (L3 .LE. 0 .OR. L3 .GT. N) GOTO 10
      DO 3 L=1,ND
    3 Y(L)=LOG10(POPNUM(L,L3))
      CALL PLOTCON (KANAL,X,Y,ND,3)
      GOTO 10
      END subroutine
      end module