      module MOD_PLOTT

      contains

      SUBROUTINE PLOTT (ND,R,T,MODHEAD,JOBNUM)

C***  DIRECT TRANSFER OF THE TEMPERATURE STRATIFICATION T(R) VERSUS LOG(R/R*-1)

      use MOD_ISRCHFLT
      use MOD_ISRCHFLE
      use MOD_PLOTANF
      implicit real*8(a-h,o-z)

C***  TRANSFER OF "YMAX" FROM SUBR. DECSTAR
      COMMON / COMPLOT / YMAX
      DIMENSION X(100),Y(100)
      DIMENSION R(ND),T(ND)
      CHARACTER HEAD1*60,HEAD2*60,MODHEAD*104
 
C***  INITIALIZATION
      KANAL=2
      OPEN (KANAL,FILE='PLOT')
C      CALL JSYMSET (2LG2,'TRANSFER')
      write (6,*) 'T-STRATIFICATION TO BE ROUTED'
 
C***  RADIUS RANGE:   2. > LOG(R/R*-1) > -3.5
      LMIN= ISRCHFLT(ND,R(1),1,101.)
      LMAX= ISRCHFLE(ND,R(1),1,1.000316228) - 1
      NDL=LMAX-LMIN+1
      IF (NDL .GT. 100) THEN
         write (6,*) 'T-PLOT IMPOSSIBLE'
         RETURN
         ENDIF
      DO 5 L=1,NDL
      X(L)=LOG10(R(LMAX+1-L)-1.d0)
      Y(L)=T(LMAX+1-L)/1.d3
    5 CONTINUE
C***  PLOT PARAMETER
      XMIN=-3.5
      XMAX=2.
      YMIN=0.
      IF (YMAX .EQ. 0.) YMAX=10.*(INT(Y(1)/10.)+1.d0)
      XSCALE=4.
      YSCALE=15./YMAX
      XTICK=0.5
      YTICK=10.
      XABST=1.
      YABST=10.
      HEAD1=' WR TEMPERATURE STRATIFICATION T(R) VERSUS LOG(R/R*-1) '
c      ENCODE(60,9,HEAD2) JOBNUM
      write (head2,9) jobnum
    9 FORMAT (19X,1HJ,I3)
      HEAD2(1:18)=MODHEAD(15:32)
      HEAD2(24:56)=' TEMPERATURE STRATIFICATION T(R)'
      CALL PLOTANF (KANAL,HEAD1,HEAD2
     $ ,60H      LOG (R/R* - 1)                                         
     $ ,60H      T(R) / KK                                              
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,X,Y,NDL,5)
 
      RETURN
      END subroutine
      end module
