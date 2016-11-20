      MODULE MOD_DLIOP

      CONTAINS

      SUBROUTINE DLIOP (I,ENTOTL,DOPAL,DETAL,VDOP,RSTAR,N,
     $                  NDIM,EINST,WEIGHT,ELEVEL,LASTIND,INDLOW,INDNUP)
C*******************************************************************************
C***  DERIVATIVES OF ALL LINE OPACITIES AND EMISSIVITIES WITH RESPECT TO EN(I)
C***  AT CURRENT DEPTH POINT
C*******************************************************************************

      implicit real*8(a-h,o-z)
     
      DIMENSION EINST(NDIM,NDIM),WEIGHT(NDIM),ELEVEL(NDIM)
      DIMENSION INDLOW(LASTIND),INDNUP(LASTIND)
      DIMENSION DETAL(LASTIND),DOPAL(LASTIND)
     
C***  C3 = 4 * PI / H / C   ( CGS UNITS )
      DATA C3 /6.3268d16 /
C***  PI8 = 8*PI
      DATA PI8 /25.1327d0 /
     
C***  LOOP OVER ALL LINE TRANSITIONS  ----------------------------------
      DO 1 IND = 1, LASTIND

      LOW = INDLOW(IND)
      NUP = INDNUP(IND)

      IF (EINST(LOW, NUP) .EQ. -2.d0) GOTO 1

      DETAL(IND) = 0.0D0
      DOPAL(IND) = 0.0D0

      IF (I .NE. LOW .AND. I .NE. NUP) GOTO 1

      XLAMCM = 1.0D0 / (ELEVEL(NUP) - ELEVEL(LOW))

!     DND = DELTA-NUE-DOPPLER (HERTZ)
!     ********************************************
!     RINAT TAGIROV:
!     The description of VDOP is given in LIOP.for
!     ********************************************
      DND = VDOP * 1.0D5 / XLAMCM

      EMINDU=EINST(NUP,LOW)*XLAMCM*XLAMCM/PI8*RSTAR
      IF (I .EQ. LOW) THEN
C***  DERIVATIVE WITH RESPECT TO LOWER LEVEL
      DETAL(IND)=.0
      ABSORP=EMINDU*WEIGHT(NUP)/WEIGHT(LOW)
      DOPAL(IND)=ENTOTL*ABSORP/DND
      ELSE
C***  DERIVATIVE WITH RESPECT TO UPPER LEVEL
      EMSPON=EINST(NUP,LOW)/XLAMCM/C3*RSTAR
      DETAL(IND)=ENTOTL*EMSPON/DND
      DOPAL(IND)=-ENTOTL*EMINDU/DND
      ENDIF
    1 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
     
      RETURN

      END SUBROUTINE

      END MODULE
