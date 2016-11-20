      module MOD_PLOTFLU
      contains
C**********  MODULNAME: PLOTFLU   ******* 28/05/87  17.26.14.******   129 KARTEN
      SUBROUTINE PLOTFLU (NF,XLAMBDA,EMFLUX,MODHEAD,JOBNUM,
     $           Y0,TEFFE,GRAD,ALDMDT,VINF,BET,PROLIB)
C***  DIRECT TRANSFER OF EMERGENT FLUX PLOT
      use MOD_TRADFUN
      use MOD_PLOTANF
      implicit real*8(a-h,o-z)

      DIMENSION XLAMBDA(NF),EMFLUX(NF)
      DIMENSION X(300),Y(300),Z(300)
      CHARACTER MODHEAD*104,HEADER*60
      LOGICAL PROLIB
     
      IF (NF.gt.300) STOP ' PLOTFLU XYZ-dim'

      KANAL=2
      OPEN (KANAL,FILE='PLOT')
C      CALL JSYMSET (2LG2,'TRANSFER')
      write (6,*) 'FLUX PLOT DATA TO BE ROUTED'

C***  MAX NUMBER OF PLOT-POINTS
      NFPMA=198
     
      XMIN=2.
      XMAX=4.5
      XSCALE=8.
      XTICK=.5
      XABST=1.
     
      IF (PROLIB) THEN
       IF (NF.LT.NFPMA) THEN
         NF3 = (NF/3)*3
        ELSE
         NF3 = NFPMA
       ENDIF
       RWLAE=NF3
       ENCODE (60,5,HEADER) RWLAE,TEFFE,GRAD,ALDMDT,VINF,BET,
     $      MODHEAD(15:22)
    5  FORMAT (F10.2,F11.2,F8.3,F8.3,F8.0,F4.1,A11)
C*    5    FORMAT (F7.1,F7.0,F6.2,F6.2,F6.0,F4.1,A9)
       WRITE(6,6) HEADER(1:48)
    6  FORMAT (' PLOTFLU-PROLIB: ',A48)
       NFX=0
       IF (NF.GT.150) THEN
        DO 11 K=1,NF
         IF (XLAMBDA(K).LT.7080..OR.XLAMBDA(K).GT.8220.) THEN
         IF (XLAMBDA(K).LT.8260..OR.XLAMBDA(K).GT.9330.) THEN
         IF (XLAMBDA(K).LT.9370..OR.XLAMBDA(K).GT.10820.) THEN
         IF (XLAMBDA(K).LT.10850..OR.XLAMBDA(K).GT.12800.) THEN
         IF (XLAMBDA(K).LT.12815..OR.XLAMBDA(K).GT.15200.) THEN
         IF (XLAMBDA(K).LT.15300..OR.XLAMBDA(K).GT.20550.) THEN
         IF (XLAMBDA(K).LT.20600..OR.XLAMBDA(K).GT.30900.) THEN
         IF (XLAMBDA(K).LT.31000..OR.XLAMBDA(K).GT.74300.) THEN
         IF (NFX.LT.NFPMA) THEN
         IF (EMFLUX(K).GT.0.) THEN
          NFX=NFX+1
          YNFX  =LOG10(XLAMBDA(K))
C*        Y(NFX)=YNFX-3.
          Y(NFX)=YNFX
          Z(NFX)=LOG10(EMFLUX(K))-2.d0*YNFX+Y0
          ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF
   11   CONTINUE
       ELSE
C***   case: PROLIB but less points than 150
       NFMIN=NF3-3
       DO 10 K=1,NFMIN
        YNFX  =LOG10(XLAMBDA(K))
C*      Y(K)=YNFX-3.d0
        Y(K)=YNFX
        IF (EMFLUX(K).GT.0.d0) THEN
        Z(K)  =LOG10(EMFLUX(K))-2.*YNFX+Y0
        ELSE
        Z(K)=-30.d0
        ENDIF
   10  CONTINUE
C*** 4th last through second last 
      DO 14 NFX=NFMIN+1,NF3
      NFKK=NF-NF3+NFX-1
      YNFX  =LOG10(XLAMBDA(NFKK))
C*      Y(NFX)=YNFX-3.
      Y(NFX)=YNFX
      Z(NFX)=LOG10(EMFLUX(NFKK))-2.d0*YNFX+Y0
   14 CONTINUE
      ENDIF
C*** reverse order
      DO 12 K=1,NF3
      X(K)=Y(NF3+1-K)
   12 CONTINUE
      DO 13 K=1,NF3
      Y(K)=Z(NF3+1-K)
   13 CONTINUE
      ELSE
C*** case: no PROLIB option
c      ENCODE (60,3,HEADER) JOBNUM
	write (header,3) jobnum
    3 FORMAT ('EMERGENT FLUX OF MODEL',23X,'JOB ',I3)
      HEADER(24:41)=MODHEAD(15:32)
      NF3=NF
      DO 1 K=1,NF
      X(K)=LOG10(XLAMBDA(K))
      Y(K)=TRADFUN(XLAMBDA(K),EMFLUX(K))/1000.
    1 CONTINUE
      ENDIF
     
      YMIN=Y(1)
      DO 2 K=2,NF3
      IF (YMIN .GT. Y(K)) YMIN=Y(K)
    2 CONTINUE
     
      YMIN=INT(YMIN/5.)*5.
      YMAX=YMIN+75.
      YSCALE=.2
      YTICK=5.
      YABST=10.
     
      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,60H      LOG ( LAMBDA/ANGSTROEM )                                           
     $ ,60H       T-RAD / KILO-KELVIN                                               
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,X,Y,NF3,5)
     
      RETURN
      END subroutine
      end module
