      MODULE MOD_RGRIDM

      real*8, allocatable :: XNETAB(:)

      CONTAINS

      SUBROUTINE RGRIDM(RADIUS,ENTOT,RHO,T,R,FAL,AINCRIT,RMAX,RSTAR,AMU,ATMEAN,NDDIM,ND)

      USE MOD_ERROR
      USE COMMON_BLOCK

      IMPLICIT REAL*8(A - H, O - Z)

      integer,intent(out) :: ND
      real*8,intent(inout):: T,ENTOT,R,RADIUS,RHO
      real*8,intent(in)   :: RMAX,RSTAR,AMU,ATMEAN
      logical,intent(in)  :: FAL
      CHARACTER*8,intent(in):: AINCRIT(*)
      integer,intent(in) ::NDDIM

      DIMENSION RHO(NDDIM),T(NDDIM),PTAB(NDDIM),
     $ ENTOTN(NDDIM),DELR(NDDIM),
     $ ENTOT(NDDIM),RADIUS(NDDIM),XX(NDDIM),YY(NDDIM),ZZ(NDDIM)
     $ ,height_grid(nddim),TVAR(NDDIM),PVAR(NDDIM),RHOVAR(NDDIM)
     $ ,VDTAB(NDDIM)

      CHARACTER(5) STRL
      COMMON /COMRPAR/ RPAR

      DIMENSION R(NDDIM)
      CHARACTER RPAR*80

      LOGICAL VAR

C***  ak: BOLTZMANN KONSTANTE
      real*8,parameter :: ak=1.38062259d-16
      real*8,parameter :: MUN=1.66054d-24

c**************************************************************************
CMH   ENTOTN(L): (PARTICLE+ELECTRONS)/VOLUME
CMH   ENTOT: HEAVY PARTICLE DENSITY
CMH   RHO(L): MASS COLUMN IN G/CM^2
CMH   XNETAB(L): ELECTRON DENSITY IN ELECTRON/CM^3
CMH   PTAB(L): pressure
CMH   RSUN = 6.960E^10 CM
CMH   RADIUS = HEIGHT IN UNITS OF SOLAR RADII
CMH   R = RADIUS
C**************************************************************************

      IF (.NOT. ALLOCATED(XNETAB)) allocate(XNETAB(NDDIM))
      VAR = .FALSE.
      PRINT *, 'RGRIDM, VAR = ',VAR
      ND = 1
      VDTAB(1:NDDIM) = 0
      IF_FAL:IF (FAL) THEN

              OPEN(9, FILE = 'FAL_VD', STATUS = 'OLD')

  400         READ(9, *, end = 55) height_grid(ND), T(ND), XNETAB(ND), ENTOT(ND), VDTAB(ND)

              RADIUS(ND) = 1.D0 + height_grid(ND) * 1.D5 / RSTAR  ! HEIGTH in km, RSTAR in cm

              R(ND) = RADIUS(ND)

          IF(ND.GT.NDDIM) THEN
              WRITE(6,*) 'RADIUS DIMENSION INSUFFICIENT'
              STOP 'ERROR'
          ENDIF
          ND = ND+1
          GOTO 400

   55     CLOSE (9)
          ND=ND-1
          CONTINUE
      ELSE IF_FAL
CMH   READ ATMOSPHERE FROM 'TABLE'
          OPEN (8,FILE='TABLE', STATUS='OLD')
  500         READ (8,*,END=66) RHO(ND),T(ND),PTAB(ND),XNETAB(ND),
     $ XX(ND),YY(ND),ZZ(ND)

         PTAB(L)=RHO(L)*10.**4.44
  

          IF(ND.GT.NDDIM) THEN
              WRITE(6,*) 'RADIUS DIMENSION INSUFFICIENT'
              STOP 'ERROR'
          ENDIF
          ND = ND + 1
          GOTO 500

   66     CLOSE (8)
          ! Sanity checks, table must be non empty and t and ptab != 0
          ND = ND - 1

          IF (ND .EQ. 0) CALL ERROR('RGRIDM: EMPTY TABLE')
          DO L=1,ND
              WRITE( STRL,'(I5)') L
              IF (T(L).EQ.0.) CALL ERROR('RGRIDM: T(' // STRL // ') = 0')
              IF (PTAB(L).EQ.0.) CALL ERROR('RGRIDM: PTAB(' // STRL  //') = 0')
          END DO
      CONTINUE
C*******************************************************************************
      IFVAR: IF (VAR) THEN
          PRINT *,'RGRIDM: VARIABILITY OF INPUT PARAMETERS ACTIVE!!!!'
          OPEN (9,FILE='VAR',STATUS='OLD')
          DO L=1,ND
              READ (9,*) TVAR(L),PVAR(L),RHOVAR(L)

              T(L)=T(L)*(1.+TVAR(L))

C             NEW - TAKING INTO ACCOUNT TURBULEN PRESSURE
              
              ENTOTN(L) = PTAB(L)/(AK*T(L)+0.5*ATMEAN*MUN*ZZ(L)**2.)
              
              ENTOT(L)  = ENTOTN(L)-XNETAB(L)

              ENTOT(L) = ENTOT(L)*(1.+RHOVAR(L))
          ENDDO
          CLOSE (9)
          CONTINUE
      ELSE IFVAR
          DO L=1,ND

C             NEW - TAKING INTO ACCOUNT TURBULEN PRESSURE
              
              ENTOTN(L) = PTAB(L)/(AK*T(L)+0.5*ATMEAN*MUN*ZZ(L)**2.)
              ENTOT(L)  = ENTOTN(L)-XNETAB(L)

          ENDDO

          PRINT *,'RGRIDM: VARIABILITY NOT ACTIVE!!!!'
      ENDIF IFVAR

      DO L=1,ND-1
          DELR(L) =(2/(AMU*ATMEAN*RSTAR))*(RHO(L+1)-RHO(L))/
     $    (ENTOT(L+1)+ENTOT(L))

      ENDDO

      ND=ND-1
      RADIUS(ND)= 1.
      R(ND)=RADIUS(ND)
      DO K=1,ND-1
          RADIUS(ND-K) = RADIUS(ND-K+1)+ DELR(ND-K)
          R(ND-K) = RADIUS(ND-K)
      ENDDO

      ENDIF IF_FAL

      DPN = ND

      END subroutine

      end module
