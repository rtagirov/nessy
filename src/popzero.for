      MODULE MOD_POPZERO

      CONTAINS

      SUBROUTINE POPZERO(T,RNE,POPNUM,DEPART,ENTOT,ITNE,N,ENLTE,
     $                   WEIGHT,NCHARG,EION,ELEVEL,EN,EINST,LEVEL,
     $                   XLAMBDA,FWEIGHT,XJC,NF,XJL,IFRRA,ITORA,ALPHA,
     $                   SEXPO,AGAUNT,MODHEAD,MODHOLD,JOBNUM,
     $                   LASTIND,ND,LSRAT,SIGMAKI,ALTESUM,COCO,KEYCOL,NOM,NATOM,
     $                   KODAT,NFIRST,NLAST,WAVARR,SIGARR)

      !*******************************************************************************
      !***  IN CASE OF "STARTJOB" (JOB.EQ.'WRSTART'):
      !***  IF (OLDSTART) THEN
      !***  POPNUMBERS ARE READ FROM APPLICABLE OLD MODEL FILE
      !***  ELSE
      !***  CALCULATION OF THE NLTE POP.NUMBERS BY SOLVING THE RATE EQUATIONS.
      !***   --- NO APPROXIMATE LAMBDA OPERATORS ARE INCLUDED HERE -----
      !***  THIS ROUTINE IS SIMILAR TO SUBR. NEWPOP
      !***  IT IS USED BY MAIN PROGRAM STEAL IF ALL GAMMA'S ARE ZERO
      !***  IT HAS THE ADVANTAGE THAT THE RATE COEFFICIENTS ARE PROVIDED
      !***  WHICH CAN BE PRINTED WITH SUBR. PRIRAT
      !***  ENDIF

      !MH   ENLTE:  LTE popnulation number
      !MH   POPNUM: NLTE population numbers
      !MH   DEPART: POPNUM / ENLTE
      !MH   ENE:    absolute electron density
      !MH   RNE:    relative electron density
      !MH   ENTOT:  particle density
      !MH   EN:     POPNUM of level at depth point L
      !*******************************************************************************

      USE MOD_BFCROSS
      use MOD_NLTEPOP
      use MOD_LTEPOP
      use MOD_PRIRAT
      use ABUNDANCES

      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)

      integer,intent(in) :: JOBNUM

      DIMENSION T(ND),ENTOT(ND),RNE(ND),POPNUM(ND, N),ITNE(ND)
      DIMENSION DEPART(ND, N)
      DIMENSION NCHARG(N), EN(N), ENLTE(N)
      DIMENSION NOM(N)
      DIMENSION KODAT(NATOM), NFIRST(NATOM), NLAST(NATOM)
  
      LOGICAL KONVER
      CHARACTER MODHOLD*104,MODHEAD*104,CNAME*10,CHAR*1,TEST*10,JOB*7
      CHARACTER LEVEL(N)*10
      CHARACTER*4 KEYCOL(N, N)
      character*10 cread,help
      DIMENSION WAVARR(N, NF),SIGARR(N, NF)

      real*8	POPHIIL, POPHML, POPHIL

      real*8,dimension(*) :: ALTESUM, COCO, FWEIGHT

      real*8, dimension(N, N) :: CRATE, RRATE, ratco

      real*8,dimension(*) :: WEIGHT,XJC,XJL
      real*8 SIGMAKI(NF,N), ELEVEL(N), EION(N), EINST(N, N)
      real*8 XLAMBDA(NF), ALPHA(*),SEXPO(*)
      character*8 :: agaunt(N)
      COMMON /GIIIERR/  NTUP,NTLOW,NFUP,NFLOW,NHELP
      COMMON / COMNEGI / NEGINTL
      real*8, allocatable::ABXYZ_new(:)

      LOGICAL :: JOB_COND

      NTUP=0
      NTLOW=0
      NFUP=0
      NFLOW=0
      NEGINTL=0

      rewind 99; read (99,'(A7)',err=666) JOB
      print '(A,A)',' file 99 job=',job

      JOB_COND = JOB .EQ. 'wrstart'

      IF (JOB_COND .AND. OLDSTART) THEN
!***  START WITH POPNUMBERS WHICH ARE READ FROM AN APPLICABLE OLD MODEL FILE

        print *,' population numbers copied from an existing model'
        ifl=9
        open (ifl,file='OLDMODEL',status='old',err=777)
        CNAME='MODHEAD'

        READ (ifl,'(A10)') cread
        do while (cread(1:1).eq.' ')
          help=cread(2:10)
          cread=help
        enddo
        ERR_CREAD: if (cread.ne.cname) then
          write (6,*) 'READMSC: KEYWORD MISMATCH'
          write (6,*) cread,char

          stop
        endif ERR_CREAD
        READ (ifl,'(A104)') modhold
        print *, ' POPS FROM OLD MODEL: '
        print ('(A104)'), MODHOLD

        test=' '
        DO WHILE (test.NE.'N')
          READ (ifl,'(A10)',err=777) CNAME
          char=' '
          i=0
          do while (char.eq.' ')
            i=i+1
            ! if (i.ge.11) stop 'error'
            char=cname(i:i)
            if (i.ge.10) char='X'
          enddo
          test=char
        ENDDO
        read (ifl,*) NOLD
        IF (NOLD .NE. N) THEN
          write (6,*) 'POPZERO-N: INAPPLICABLE OLD MODEL FILE'
          STOP 'ERROR'
        ENDIF

        test='  '
        DO WHILE (test(1:2).NE.'ND')
          READ (ifl,'(A10)',err=777) CNAME
          char=' '
          i=0
          do while (char.eq.' ')
            i=i+1
            ! if (i.ge.10) stop 'error'
            char=cname(i:i)
            if (i.ge.9) char='X'
          enddo
          test=cname(i:i+1)
        ENDDO
        read (ifl,*) NDOLD
        IF (NDOLD .NE. ND) THEN
          write (6,*) 'POPZERO-ND: INAPPLICABLE OLD MODEL FILE'
          STOP 'ERROR'
        ENDIF
        test='      '
        DO WHILE (test(1:3).NE.'RNE')
          READ (ifl,'(A10)',err=777) CNAME
          char=' '
          i=0
          do while (char.eq.' ')
            i=i+1
            ! if (i.ge.8) stop 'error'
            char=cname(i:i)
            if (i.ge.8) char='X'
          enddo
          test=cname(i:i+2)
        ENDDO
        READ (ifl,*) (rne(i),i=1,nd)
        ! CALL READMS (ifl,RNE,ND,CNAME,IERR)

        test='      '
        DO WHILE (test(1:6).NE.'POPNUM')
            READ (ifl,'(A10)',err=777) CNAME
            char=' '
            i=0
            do while (char.eq.' ')

              i=i+1
              ! if (i.ge.5) stop 'error'
              char=cname(i:i)
              if (i.ge.5) char='X'
            enddo
            test=cname(i:i+5)
        ENDDO
        READ (ifl,*) ((popnum(i,j),i=1,ND),j=1,N)
        ! CALL READMS (9,POPNUM,ND*N,6HPOPNUM,IERR)
        !***  ARRAYS "ITNE", "DEPART" ARE NEEDED FOR PRINTOUT BY SUBR. PRIPOP
        ITNEL=0
        DO L=1,ND
          ITNE(L)=ITNEL
          TL=T(L)
          ENE=RNE(L)*ENTOT(L)
          if (.NOT. allocated(ABXYZ_new)) allocate(ABXYZ_new(NATOM))
          ABXYZ_new(1:NATOM)=ABXYZn(1:NATOM,L)
          CALL LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $                 ABXYZ_new,NFIRST,NLAST,NATOM)
          DEPART(L,:N)=POPNUM(L,:N)/ENLTE(:N)
        if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)    

        ENDDO

         close (ifl)

      ELSEIF (JOB_COND) THEN

        !***  START WITH CALCULATION OF THE NLTE POPNUMBERS
        !***  EPSNE = ACCURACY LIMIT FOR THE ELECTRON DENSITY ITERATION
        EPSNE=0.001
        !***  MAX. NUMBER OF ITERATIONS
        ITMAX=120
        !***  GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
        !***  SIGMAKI(K,LOW) IN CM**2
        CALL  BFCROSS(SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,
     $                XLAMBDA(1 : NF),ALPHA,SEXPO,AGAUNT,NOM,WAVARR(:, 1 : NF),SIGARR(:, 1 : NF))

        !*** CHANGES BY MARGIT HABERREITER
        PRINT*, 'POPZERO: non-LTE pop numbers replaced by LTE populations.'

        !*** LOOP OVER ALL DEPTH POINTS --------------------------------------
        DO 1 L = 1, ND

          TL = T(L)

          ENTOTL = ENTOT(L)

          RNEL = RNE(L)

          IF (.NOT. ALLOCATED(ABXYZ_new)) ALLOCATE(ABXYZ_new(NATOM))

          ABXYZ_new(1 : NATOM) = ABXYZn(1 : NATOM, L)

          IF (LTE_RUN .OR. CONST_ELEC) THEN ! LTE RUN OR PRE-SET ELECTRON CONCENTRATION

              ENE = RNEL * ENTOTL

              CALL LTEPOP(N, ENLTE, TL, ENE, WEIGHT, NCHARG, EION, ELEVEL, NOM,
     $                    ABXYZ_new, NFIRST, NLAST, NATOM)

              DEPART(L, 1 : N) = 1.0D0
              POPNUM(L, 1 : N) = ENLTE(1 : N)

          ELSEIF (.NOT. LTE_RUN .AND. .NOT. CONST_ELEC) THEN
!          ELSEIF (.NOT. LTE_RUN) THEN

          !***  ITERATION FOR THE ELECTRON DENSITY
              ITNEL = 0
 13           ITNEL = ITNEL + 1

              ENE = RNEL * ENTOTL

              IF (ENE .LT. 0) PRINT*, 'POPZERO: ENE = ', ENE

          !MH proton density
              POPHIIL = POPNUM(L, NLAST(1)) * ENTOTL
              POPHML =  POPNUM(L, 1) *        ENTOTL
              POPHIL =  POPNUM(L, 2) *        ENTOTL

              CALL LTEPOP(N, ENLTE, TL, ENE, WEIGHT, NCHARG, EION, ELEVEL, NOM,
     $                    ABXYZ_new, NFIRST, NLAST, NATOM)

              CALL NLTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,
     $                     EN,EINST,XLAMBDA,FWEIGHT,XJC,NF,ITNEL,L,LEVEL,XJL,ND,LASTIND,
     $                     CRATE,RRATE,RATCO,SIGMAKI,ALTESUM,COCO,
     $                     KEYCOL,NOM,NATOM,ABXYZ_new,KODAT,NFIRST,NLAST,
     $                     POPHIIL, POPHML, POPHIL)

              EN(1 : N) = ENLTE(1 : N)

              RNEOLD = RNEL

              RNEL = SUM(NCHARG(1 : N) * EN(1 : N))

              KONVER = ABS(RNEOLD - RNEL) .LT. EPSNE

              IF (ITNEL .GT. 10) STOP 'electron iteration problem'

              IF (.NOT. KONVER .AND. ITNEL .LT. ITMAX) GOTO 13
              IF (.NOT. KONVER) ITNEL = -ITMAX
     
              ITNE(L) = ITNEL
              RNE(L) = RNEL

              DEPART(L, 1 : N) = EN(1 : N) / ENLTE(1 : N)
              POPNUM(L, 1 : N) = EN(1 : N)

          ENDIF

          IF (ALLOCATED(ABXYZ_new)) DEALLOCATE(ABXYZ_new)

    1 ENDDO

      ELSE

        !***  GAMMA = 0 
        !***  EPSNE = ACCURACY LIMIT FOR THE ELECTRON DENSITY ITERATION
        EPSNE=0.001
        !***  MAX. NUMBER OF ITERATIONS
        ITMAX=10
        !***  GENERATE ONCE FOR ALL PHOTOCROSSSECTIONS AT ALL FREQUENCIES
        !***  SIGMAKI(K,LOW) IN CM**2
         CALL  BFCROSS(SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,
     $                 XLAMBDA(1 : NF),ALPHA,SEXPO,AGAUNT,NOM,WAVARR(:, 1 : NF),SIGARR(:, 1 : NF))

        !***  CHANGES BY MARGIT HABERREITER
        print *,' Gamma=0 branch !!!!!'
        !***  LOOP OVER ALL DEPTH POINTS  --------------------------------------

          MAIN_LOOP: DO L=1,ND

          TL=T(L)
          ENTOTL=ENTOT(L)
          RNEL=RNE(L)
          if (.NOT. allocated(ABXYZ_new)) allocate(ABXYZ_new(NATOM))
          ABXYZ_new(1:NATOM)=ABXYZn(1:NATOM,L)
          !***  ITERATION FOR THE ELECTRON DENSITY
          ITNEL=0
 23       ITNEL=ITNEL+1
         

          ENE = RNEL * ENTOTL

          POPHIIL = POPNUM(L, NLAST(1)) * ENTOTL
          POPHML =  POPNUM(L, 1) *        ENTOTL
          POPHIL =  POPNUM(L, 2) *        ENTOTL

          CALL LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $           ABXYZ_new,NFIRST,NLAST,NATOM)

          CALL NLTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,
     $                 EN,EINST,XLAMBDA,FWEIGHT,XJC,NF,ITNEL,L,LEVEL,
     $                 XJL,ND,LASTIND,
     $                 CRATE,RRATE,RATCO,SIGMAKI,ALTESUM,COCO,
     $                 KEYCOL,NOM,NATOM,ABXYZ_new,KODAT,NFIRST,NLAST,
     $                 POPHIIL, POPHML, POPHIL)

!         for helium always LTE        

!          EN(13:23)=ENLTE(13:23)

!         for helium always LTE 

          RNEOLD=RNEL
          RNEL = sum(NCHARG(:N)*EN(:N))
          KONVER= ABS(RNEOLD-RNEL) .LT. EPSNE
          IF (.NOT. KONVER .AND. ITNEL .LT. ITMAX) GOTO 23
          IF (.NOT. KONVER) ITNEL=-ITMAX
          ITNE(L)=ITNEL
          RNE(L)=RNEL
          DEPART(L,:N)=EN(:N)/ENLTE(:N)
          POPNUM(L,:N)=EN(:N)
     
! test hminus                                                                   
    !         POPNUM(L,1)=ENLTE(1)
                ! test hminus  



!***  ===  PRINTOUT OF RATE COEFFICIENTS ETC.  ===
          IF (LSRAT.NE.-1) THEN
            IF ((L.GE.IFRRA.AND.L.LE.ITORA).OR.ITORA.EQ.0) THEN
              NETTO=0
              LM1=L-1
              IF (IFRRA.GT.0) LM1=L-IFRRA
              IF (((LM1)/LSRAT)*LSRAT.EQ.(LM1).OR.L.EQ.ND)
     $          CALL PRIRAT (IDUMMY,N,LEVEL,L,CRATE,RRATE,RATCO,EN,
     $            IFRRA,MODHEAD,JOBNUM,NETTO )
            ENDIF
          ENDIF
          if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)

        ENDDO MAIN_LOOP
      !***  ENDLOOP  ---------------------------------------------------------
      ENDIF
      if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)
      RETURN
  666 continue
      if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)
      print *,' error reading job control file 99'
      stop ' error popzero'
777   write (6,*) 'POPZERO: ERROR WHEN READING OLD MODEL FILE'
      STOP 'OPMS9'
      END SUBROUTINE
      end module
