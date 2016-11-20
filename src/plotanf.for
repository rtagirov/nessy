      module MOD_PLOTANF
      contains
C**********  MODULNAME: PLOTANF   ******* 24/03/87  19.38.27.******    21 KARTEN
      SUBROUTINE PLOTANF (KANAL,NPLOT,NHEAD,NX,NY,
     $ B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6,
     $ X,Y,N,ISYMBOL )
      use MOD_PLOTTAB
      implicit real*8(a-h,o-z)
      real*8 X(*),Y(*)
      CHARACTER*60  NHEAD,NPLOT

C*      DIMENSION NPLOT (6),NHEAD(6),NX(6),NY(6)
C*      WRITE (KANAL,1) NPLOT
      WRITE (KANAL,2) NHEAD
C*      WRITE(KANAL,3) NX
C*      WRITE (KANAL,4) NY
C*    1 FORMAT (10H2$$$$$$$$=,/,* PLOT   :*,6A10)
    2 FORMAT (A60)
C*    2 FORMAT (* HEADER :*,6A10)
C*    3 FORMAT (* X-ACHSE:*,6A10)
C*    4 FORMAT (* Y-ACHSE:*,6A10)
C*      WRITE (KANAL,5)
C*     $ B1,B2,B3,B4,B5,B6,C1,C2,C3,C4,C5,C6
C*    5 FORMAT (5X,
C*     $ *MASSTAB    MINIMUM    MAXIMUM    TEILUNGEN  BESCHRIFT. DARUNTER*
C*     $ / ,* X: *,6(G10.4,1X),/,* Y: *,6(G10.4,1X))
      CALL PLOTTAB (KANAL,X,Y,N,ISYMBOL)
      RETURN
      END subroutine
      end module