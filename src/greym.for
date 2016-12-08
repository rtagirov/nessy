      module MOD_GREYM
      contains
C**********  MODULNAME: GREY      ******* 06/08/87  20.46.41.******   233 KARTEN
      SUBROUTINE GREYM(ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,ENTOT,RNE,RSTAR, 
     $                 ALPHA,SEXPO,AGAUNT,POPNUM,TAUROSS,R23,TTABLE,
     $                 LBKG,XLBKG1,XLBKG2,N, 
     $                 LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,ENLTE,KODAT, 
     $                 NOM,NFIRST,NLAST,NATOM,WAVARR,SIGARR)
C***  CALLED BY WRSTART 
C*********************************************************************** 
C***  COMPUTATION OF THE TEMPERATURE STRUCTURE 
C*********************************************************************** 
C***  IF (TTABLE=.TRUE.), 
C***     THE TEMPERATURE STRUCTURE IS NOT CALCULATED BUT 
C***     ASSUMED TO BE READ FROM FILE. 
C***  ELSE IF (SPHERIC=.TRUE.) 
C***     THE TEMPERATURE STRUCTURE IS CALCULATED 
C***     AS FOR A SPHERICAL, GREY ATMOSPHERE. THE EDDINGTON FACTOR F IS 
C***     APPROXIMATED BY ASSUMING A GEOMETRICALLY DILUTED RADIATION FIELD 
C***     EMERGING FROM A SPHERE OF RADIUS R23 (I.E. THE RADIUS WHERE 
C***     TAU-ROSSELAND= 2/3). 
C***     THE SOLUTION OF THE 1. MOMENT EQUATION IS OBTAINED BY USING THE 
C***     IMSL-ROUTINE IVPRK (DVERK). DTDR IS A USER-PROVIDED COEFFICIENT FUNCTION. 
C***  ELSE 
C***     THE TEMPERATURE STRUCTURE IS CALCULATED AS FOR A PLANE-PARALLEL, 
C***     GREY ATMOSPHERE. 
C***  ENDIF 
C*** 
C***  IF (TMODIFY .NE. .0), THE CALCULATED TEMPERATURE STRUCTURE IS FINALLY 
C***  MODIFIED BY A FACTOR R**TMODIFY . 
C***  IF (TMIN .GT. 0), THE TEMPERATURE STRUCTURE IS MODIFIED NOT TO FALL 
C***  BELOW TMIN. 
C***  FOR THE CALCULATION OF THE LTE IONIZATION STRUCTURE A MINIMUM 
C***  TEMPERATURE OF Teff*0.8112 IS USED
C*********************************************************************** 
      use MOD_LIPO
      use MOD_OPAROSS_M
      use MOD_LTEPOP
      use MOD_RKQS
      use MOD_DTDR
      use MOD_ODEINT
      use MOD_ERROR
      use ABUNDANCES
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /COMDTDR/ OPAMEAN,R23COM,R1COM,R13COM
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC 
      DIMENSION XLAMBDA(NF),FWEIGHT(NF) 
      DIMENSION T(ND),RADIUS(ND),RNE(ND),ENTOT(ND),TAUROSS(ND) 
      DIMENSION NCHARG(N),ENLTE(N)
      DIMENSION POPNUM(ND,N) 
      DIMENSION NOM(N) 
      DIMENSION KODAT(NATOM),NFIRST(NATOM),NLAST(NATOM) 
      
      real*8, allocatable::ABXYZ_new(:)

      !DIMENSION PARAM(50)
      real*8,DIMENSION(*) ::SEXPO,EINST,EION,ELEVEL,WEIGHT,ALPHA
C***      DIMENSION COMVEC(24),WORKSP(1,9) 
CMH  LBKG - KEYWORD FOR NON-LTE OPACITY DISTRIBUTION FUNCTIONS
CMH  XLBKB1, XLBKG2: WAVELENTH RANGE FOR THE ODF
      INTEGER XLBKG1,XLBKG2
      LOGICAL TTABLE, SPHERIC,LBKG 
      CHARACTER*10 LEVEL(N) 
C***  PARAM(4) = NUMBER OF ALLOWED INTEGRATION STEPS
      !DATA PARAM / 3 * 0.,20000, 46 * 0.0 / 
      real*8,parameter :: ak=1.38062259d-16 
      REAL*8 :: RNEL
      character*8 :: agaunt(N)
      real*8 TL_(1)
C******************************************************************************
C***  CHANGES BY MARGIT HABERREITER, 20 MAY, 2002
CMH  LEVLOW NEEDS TO BE DEFINED, AS IT IS USED AS A KEYWORD TO SELECT THE 
CMH  ELEMENT AND LEVEL TO READ THE CONTINUUM OPACITIES FROM AN INPUT TABLE
      !CHARACTER*10 LEVLOW
      DIMENSION WAVARR(N,NF),SIGARR(N,NF)
C****************************************************************************** 
C***************************************************************************** 
C***  ABUNDANCES OF HELIUM AND HYDROGEN FOR THE MODIFIED START APPROXIMATION 
C***************************************************************************** 

            ABHE=0. 
            ABH=0. 
            IF (KODAT(1) .GT. 0) ABHE=ABXYZ_small(KODAT(1)) 
            IF (KODAT(2) .GT. 0) ABH=ABXYZ_small(KODAT(2))

!**********************************************************
!RT
!            PRINT*, 'ACHTUNG: ABHE = ', ABHE, KODAT(1), ABXYZ_small(KODAT(1))
             ABHE = 1.0D0 - 9.11E-01 - 2.57E-05

!*********************************************************

C***************************************************************************** 
C***  FIRST RGRID-POINT:  R=RMAX (L=1) 
C***************************************************************************** 
      TAUROSS(1)=.0 
      IF (.NOT. TTABLE) THEN 
         IF (SPHERIC) THEN 
            T(1)=TEFF*0.7071068/SQRT(RADIUS(1)) 
            ELSE 
            T(1)=TEFF*0.8112 
            ENDIF 
         ENDIF 
C***************************************************************************** 
C***  ITERATION LOOP  ************************************************** 
C***  MAIN PURPOSE IS THE ITERATIONVOF R23 IN THE SPHERICAL CASE 
C***  START VALUE OF THE RADIUS R23  (ONLY USED IN THE SPHERICAL CASE) 
      R23COM=1. 
      R1COM=1. 
      R13COM=1. 
      ITER=0 
  100 ITER=ITER+1 
c     print *, iter,' loop 100; rnel=',rnel 
      DTMAX=.0 
C***  INITIALIZATION OF VARIABLES USED IN THE IMSL-SUBROUTINE IVPRK (DVERK) 
      IDO=1 
      NEQ=1 
      TOL=0.0001
      IND=1 
      NW=1 
      rne(1)=1.e-9
      LOOPCOUNTER=0
C***************************************************************************** 
C***  LOOP OVER ALL DEPTH POINTS  ******************************************* 
C***************************************************************************** 
!       do  j=1, ND
!         do i=1, NATOM
!         print*, 'test', j, i, ABXYZn(i,j)

!         enddo
!       enddo

 !     print*, test3
      DO 10 L=1,ND 
 
        if (.NOT. allocated(ABXYZ_new))     allocate(ABXYZ_new(NATOM))
 
      ABXYZ_new(1:NATOM)=ABXYZn(1:NATOM,L)

!      print*, ABXYZ_new(1), ABXYZ_new(2), ABXYZ_new(3), ABXYZ_new(4)
!      stop
      
 !     do i=1, NATOM
 
 !     print *, 'test', i, l, ABXYZ_new(i), ABXYZn(i,L)
 !     enddo 

 
c     PRINT *, 'GREYM: WITHIN ND LOOP N=', L, ND 
          TL=T(L) 
        IF (ITER .GT. 1) then
          if (L.lt.nd) then
            TOLD=T(L+1)
          else
            TOLD=T(L)
          endif
        endif 
C***************************************************************************** 
C***  USE MODIFIED TEMPERATURE TMIN FOR THE CALCULATION OF THE EL.DENSITY 
C***************************************************************************** 
        TP=TL 
CMH  MARGIT HABERREITER: NO T-MODIFICATION IN CASE OF TTABLE
        IF (.NOT.TTABLE) THEN 
          IF (TP.LT.TMIN) TP=TMIN 
          IF (TP.LT.TEFF*0.8112) TP=TEFF*0.8112 
        ENDIF
C***************************************************************************** 
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY AT POINT L 
C***  FIRST: LTE POPNUMBERS, ITERATION OF ELECTRON DENSITY 
C***************************************************************************** 
        RNEL=RNE(L)
C***************************************************************************** 
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN 
C***************************************************************************** 
c*** first entry for rnel =  0.99997
        IF (ITER .EQ. 1 .AND. TP .LT. 10000.) THEN

            T32=TP*SQRT(TP) 

!          RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TP)*T32/2.07E-16
!     $        *ABHE/ENTOT(L))

          RNEL = ABXYZ_new(3) + ABXYZ_new(4) + ABXYZ_new(5)
!          RNEL = 0.0001

!          PRINT*, 'ACHTUNG 2: RNEL', ITER, L, RNEL

        ENDIF

        !MH   print *,ITER,'RNEL=',rnel
        !!**************************************************************
        n3 = 1
CMH    3 if (RNEL.lt. 0.d0) rnel=abs(rnel)/3.
CMH    3 if (RNEL.lt. 0.d0) rnel=abs(rnel)/1.9
        do  ! while (ABS(RNEDIF/RNEL).GT. 1.e-13 .OR. ABS(RNEDIF).GT. 1.e-13)
          if (RNEL.lt. 0.d0) rnel=abs(rnel)/3d0
          n3=n3+1
          !MH   if (l.gt.30)
          !MH   print *,n3,'loop 3','depth point=',L, ND ,'rnel' ,rnel
          ENE=RNEL*ENTOT(L)
          !BEGIN DEBUG micha
          if(n3 == 1 000 000) then
            print '("greym: RNE has difficulty to converge")'
            print '("greym: L,N,RNEL,RNEDIF=",I4,I4,e10.4," ",e10.4)',
     $             L,N,RNEDIF,RNEDIF/RNEL
          endif  !END DEBUG

          if(n3 > 20 000 000) then
            print '("greym:RNE does not converge L,N,RNEL,RNEDIF=" '//
     $      ',I4,I4,e6.4,e6.4)',
     $             L,N,RNEL,RNEDIF
          call error('greym: RNEL does not to converge')
          endif  !END DEBUG

   
          CALL   LTEPOP (N,ENLTE,TP,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM, 
     $           ABXYZ_new,NFIRST,NLAST,NATOM)
          RNEOLD=RNEL
          RNEL=sum(NCHARG*ENLTE)

!        IF (n3 == 1000000) THEN

!        WRITE(*, '(A,1x,I4,1x,I4,1x,E15.7,12(1x,E15.7))'), 'ACHTUNG 3:', ITER, L, RNEL, ENLTE(1),
!     $                                                                                  ENLTE(2),
!     $                                                                           ENLTE(3),
!     $                                                                           ENLTE(4),
!     $                                                                           ENLTE(5),
!     $                                                                           ENLTE(6), 
!     $                                                                           ENLTE(7), 
!     $                                                                           ENLTE(8), 
!     $                                                                           ENLTE(9),
!     $                                                                           ENLTE(10),
!     $                                                                           ENLTE(11),
!     $                                                                           ENLTE(12)

!         STOP

!         ENDIF

  !        print*, enlte
  !        print*, 'rnel', rnel
 !         stop
          !!************************************************************
          !DO 2 J=1,N
          !2   RNEL=RNEL+NCHARG(J)*ENLTE(J)
          !!************************************************************
          RNEDIF=RNEL-RNEOLD
          !!************************************************************
          IF (ABS(RNEDIF/RNEL) < 1e-13 .and. ABS(RNEDIF) < 1e-13)  exit
        end do

!          print*, rnel
!          stop
C*****************************************************************************
C***  STORE LTE POPNUMBERS TO BE WRITTEN AT THE MODEL FILE (START APPROXIMAT9ON@
C*****************************************************************************
!          DO 5 J=1,N 
!    5     POPNUM(L,J)=ENLTE(J)
          POPNUM(L,1:N)=ENLTE(1:N)
C***************************************************************************** 
C*** ELECTRON DENSITY ALSO STORED 
C*****************************************************************************
        RNE(L)=RNEL 
c*** pressure at depth point L
        press = entot(l)*(1.d0+rnel)*ak*T(L)
C*****************************************************************************
        IF (L .EQ. ND) GOTO 10
C***************************************************************************** 

        CALL OPAROSS_M(OPARL,ENLTE,TP,RNEL,ENTOT(L),RSTAR,N,
     $                 LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST, 
     $                 ALPHA,SEXPO,AGAUNT,NF,XLAMBDA(1 : NF),FWEIGHT(1 : NF),NOM,
     $                 WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2)
C***************************************************************************** 
C***  COMPUTATION OF THE ROSSELAND MEAN OPACITY  AT POINT L+1 
C***  IN THE FIRST ITERATION USING T(L) 
C***  IN THE FOLLOWING ITERATIONS USING TOLD 
C*****************************************************************************
        IF (ITER .EQ. 1) THEN 
          TL1=TL 
        ELSE 
          TL1=TOLD 
        ENDIF
C***************************************************************************** 
C***  USE MODIFIED TEMPERATURE TMIN FOR THE CALCULATION OF THE EL.DENSITY 
C*****************************************************************************
        TP1=TL1 
CMH  MARGIT HABERREITER: NO T-MODIFICATION IN CASE OF TTABLE
        IF (.NOT.TTABLE) THEN 
          IF (TP1.LT.TMIN) TP1=TMIN 
          IF (TP1.LT.TEFF*0.8112) TP1=TEFF*0.8112 
        ENDIF
C***************************************************************************** 
C***  FIRST: LTE POPNUMBERS, ITERATION OF ELECTRON DENSITY 
C***************************************************************************** 
        RNEL=RNE(L+1) 
C***************************************************************************** 
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN 
C***************************************************************************** 
        n13 = 0
c   13  if (rnel.lt.0.d0) rnel=abs(rnel)/1.9
   13   if (rnel.lt.0.d0) rnel=abs(rnel)/1.3
C  13   if (rnel.lt.0.d0) rnel=abs(rnel)/3
        n13 = n13+1
c       if (l.gt.30) 
c       print *,n13,'loop 13: ','depth point',l
c       print *,n13,'rnel:', rnel
c       print *,n13,'RNEDIF/RNEL',RNEDIF/RNEL
c       IF (ITER .EQ. 1 .AND. TP1 .LT. 1.E4) THEN 
c         T32=TP1*SQRT(TP1) 
c         RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TP1)*T32/2.07E-16 
c     $   *ABHE/ENTOT(L+1)) 
c       ENDIF 
C*****************************************************************************
        ENE=RNEL*ENTOT(L+1) 
        CALL  LTEPOP (N,ENLTE,TP1,ENE,WEIGHT,NCHARG,EION,
     $                ELEVEL,NOM,ABXYZ_new,NFIRST,NLAST,NATOM)

        RNEOLD=RNEL 
        RNEL=.0 
C*****************************************************************************
        DO 12 J=1,N 
   12   RNEL=RNEL+NCHARG(J)*ENLTE(J) 
C*****************************************************************************
        RNEDIF=RNEL-RNEOLD 
C***************************************************************************** 
c     IF (ABS(RNEDIF/RNEL).GT.0.01 .OR. ABS(RNEDIF).GT.0.001) 
        IF (ABS(RNEDIF/RNEL).GT. 1.e-13 .OR.
     $    ABS(RNEDIF).GT. 1.e-13) 
     $    GOTO 13 

C***************************************************************************** 
        CALL OPAROSS_M(OPARL1,ENLTE,TP1,RNEL,ENTOT(L+1),RSTAR,N,
     $                 LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST, 
     $                 ALPHA,SEXPO,AGAUNT,NF,XLAMBDA,FWEIGHT,NOM,
     $                 WAVARR,SIGARR,LBKG,XLBKG1,XLBKG2)
C***************************************************************************** 
C***  ARITHMETIC MEAN OF OPARL AND OPARL1 
C***************************************************************************** 
        OPAMEAN=0.5*(OPARL+OPARL1) 
        TAUROSS(L+1)=OPAMEAN*(RADIUS(L)-RADIUS(L+1))+TAUROSS(L) 
C***************************************************************************** 
C***  NOT TABLE
C*****************************************************************************
        IF (.NOT. TTABLE) THEN 
          IF (SPHERIC) THEN 
            RL=RADIUS(L) 
            RL1=RADIUS(L+1) 
C***************************************************************************** 
C***  DVERK (IMSL VERSION 9) IS REPLACED BY IVPRK (VERSION 10) 
C***        CALL DVERK(1,DTDR,RL,TL,RL1,TOL,IND,COMVEC,NW,WORKSP,IER)
C            print *,l,rl,rl1,tl
c***        PC version: IVPRK replaced with ...
c        CALL IVPRK (IDO,NEQ,DTDR,RL,RL1,TOL,PARAM,TL)
c***        ... numerical recipe routine ODEINT (which calls rkqs)
C*****************************************************************************
            h1=abs((rl-rl1)/10.)
            hmin=h1/1000.
            if (NEQ .ne. 1)  ! NEQ is set to 1 at top - *should* be safe
     $         call error('GREYM.for: NEQ not 1 but TL is a scalar')
            TL_=TL
            call odeint(TL_,NEQ,RL,RL1,TOL,h1,hmin,nok,
     $         nbad,DTDR,rkqs)
            TL=TL_(1)
            T(L+1)=TL
          ELSE 
C***        HOPF FUNCTION, C.F. UNSOELD P. 138 
C            PRINT *, L,TAUROSS(l+1) 
            Q=0.6940-0.1167*EXP(-1.9720*TAUROSS(L+1)) 
            T(L+1)=TEFF*(0.75*(TAUROSS(L+1)+Q))**0.25 
          ENDIF 
        ENDIF 
C***************************************************************************** 
C***  MAXIMUM TEMPERATURE CORRECTION 
C***************************************************************************** 
c***  MARGIT HABERREITER: NO T-MODIFICATION IN CASE OF TTABLE
        IF (.NOT.TTABLE) THEN
            IF (ITER .GT. 1) 
     $        DTMAX=DMAX1(DTMAX,ABS(TOLD-T(L+1))) 
        ENDIF 

        if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)

   10 CONTINUE 
C***************************************************************************** 
CC***  GIVE WORKSPACE FREE
c      IF (.NOT. TTABLE .AND. SPHERIC) THEN 
c      IDO=3 
c      CALL IVPRK (IDO,NEQ,DTDR,RL,RL1,TOL,PARAM,TL) 
c      ENDIF 
C***************************************************************************** 
C***  CALCULATE RADIUS R23 WHERE TAUROSS=2/3 
C***************************************************************************** 
      TAU23=0.666666666666 
      IF (TAUROSS(ND) .LT. TAU23) THEN 
         R23=1. 
         ELSE 
         CALL LIPO (R23,TAU23,RADIUS,TAUROSS,ND) 
      ENDIF 
      R23COM=R23 
      TAU1=1. 
      IF (TAUROSS(ND) .LT. TAU1  ) THEN 
        R1 COM=1. 
      ELSE 
         CALL LIPO (R1 COM,TAU1 ,RADIUS,TAUROSS,ND) 
      ENDIF 
      TAU13=0.333333333333 
      IF (TAUROSS(ND) .LT. TAU13 ) THEN 
        R13COM=1. 
      ELSE 
        CALL LIPO (R13COM,TAU13,RADIUS,TAUROSS,ND) 
      ENDIF 
C***************************************************************************** 
      IF (ITER .LE. 1) GOTO 100 
c***  MARGIT HABERREITER: NO T-MODIFICATION IN CASE OF TTABLE
      IF (.NOT.TTABLE) THEN    
        IF (DTMAX .GT. 10.) GOTO 100 
      ENDIF
C***************************************************************************** 
      IF (.NOT. TTABLE) THEN 
         IF (T(ND) .LE. TEFF ) THEN 
            T(ND)=TEFF 
            T(ND-1)=TEFF 
         ENDIF 
      ENDIF 
C***************************************************************************** 
C***  MODIFY THE TEMPERATURE STRATIFICATION BY A FACTOR RADIUS**TMODIFY 
C***  AND/OR BY REQUIRING A MINIMUM TEMPERATURE TMIN 
C***************************************************************************** 
CMH   MARGIT HABERREITER: NO T-MODIFICATION IN CASE OF TTABLE
C***************************************************************************** 
      if (.not.ttable) then
        DO 4 L=1,ND 


          if (.NOT. allocated(ABXYZ_new)) allocate(ABXYZ_new(NATOM))




          ABXYZ_new(1:NATOM)=ABXYZn(1:NATOM,L)

          IF (TMODIFY .EQ. .0 .AND. T(L) .GT. TMIN ) GOTO 4 
          IF (TMODIFY .NE. .0) T(L)=T(L)*RADIUS(L)**TMODIFY 
          IF (T(L) .LT. TMIN) T(L)=TMIN 
C***************************************************************************** 
C***     CALCULATE LTE POP.NUMBERS WITH MODIFIED TEMPERATURE! 
C***************************************************************************** 
          TL=T(L)
C***************************************************************************** 
c*** minimum temperature for the ionization
C***************************************************************************** 
          IF (TL.LT.TEFF*0.8112) TL=TEFF*0.8112 
          T32=TL*SQRT(TL) 
          ENTOTL=ENTOT(L) 
          RNEL=RNE(L) 
C***************************************************************************** 
C***  MODIFIED START APPROXIMATION FOR TEMPERATURE .LT. 10000K: 
C***  SAHA EQUATION FOR HEI/HEII - IONIZED HYDROGEN 
C***************************************************************************** 
          IF (TL.LT.1.E4 .AND. .FALSE.) THEN
            RNEL=ABH/2.+SQRT(ABH*ABH/4.+2.*EXP(-285645./TL)* 
     $      T32/2.07E-16*ABHE/ENTOTL) 
          ENDIF 
   23     ENE=RNEL*ENTOTL
          CALL  LTEPOP (N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,NOM,
     $              ABXYZ_new,NFIRST,NLAST,NATOM)
          RNEOLD=RNEL
          RNEL=.0 
          DO 24 J=1,N 
   24     RNEL=RNEL+NCHARG(J)*ENLTE(J) 
          RNEDIF=RNEL-RNEOLD 
          IF (ABS(RNEDIF/RNEL).GT.0.0001 .OR. 
     $      ABS(RNEDIF).GT.0.00001) 
     $      GOTO 23 
          DO 25 J=1,N 
   25     POPNUM(L,J)=ENLTE(J) 
          RNE(L)=RNEL 
          if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)

    4   CONTINUE 
c***  margit haberreiter, endif .not. ttable
      endif 
C***************************************************************************** 
      if (allocated(ABXYZ_new)) deallocate(ABXYZ_new)
      RETURN 
      END subroutine
      end module
