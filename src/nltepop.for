      MODULE MOD_NLTEPOP

      CONTAINS

      SUBROUTINE NLTEPOP(N,ENLTE,TL,ENE,WEIGHT,NCHARG,EION,ELEVEL,
     $                   EN,EINST,XLAMBDA,FWEIGHT,XJC,NF,ITNE,L,LEVEL,
     $                   XJL,ND,LASTIND,
     $                   CRATE,RRATE,RATCO,SIGMAKI,ALTESUM,COCO,
     $                   KEYCOL,NOM,NATOM,ABXYZ,KODAT,NFIRST,NLAST,
     $                   POPHIIL, POPHML, POPHIL)

C******************************************************************************
C***  CALCULATION OF NEW NLTE POPULATION NUMBERS
C***  SOLUTION OF THE LINEAR RATE EQUATION SYSTEM
C***  CALLED TWICE BY POPZERO
C***  HMINUS IS LEVEL 1
C******************************************************************************

      use MOD_COLLI
      use MOD_RADIO
      use MOD_INV

      implicit real*8(a-h,o-z)

      DIMENSION EINST(N,N),CRATE(N,N),RRATE(N,N)
      DIMENSION RATCO(N,N)
      DIMENSION ENLTE(N),EN(N),NCHARG(N),WEIGHT(N)
      DIMENSION EION(N),ELEVEL(N)
      DIMENSION NOM(N)
      DIMENSION KODAT(NATOM),NFIRST(NATOM),NLAST(NATOM)
      real*8, allocatable::ABXYZ(:)
      DIMENSION XLAMBDA(NF),FWEIGHT(NF)
      CHARACTER*4 KEYCOL(N,N)
      CHARACTER LEVEL(N)*10
      REAL*8	POPHIIL, POPHML, POPHIL
      real*8    COCO(*),ALTESUM(*),XJC(*),XJL(ND,LASTIND),SIGMAKI(*)

      integer :: fi, si

C***  SET UP THE COEFFICIENT MATRICES CRATE AND RRATE FOR ALL ELEMENTS
      CALL COLLI(N,ENLTE,ENLTE,TL,ENE,NCHARG,ELEVEL,EINST,CRATE,
     $           EION,COCO,KEYCOL,WEIGHT,ALTESUM,NATOM,NOM,KODAT,
     $           POPHIIL, POPHML, POPHIL, LEVEL, 0, L)

CMH - new: POPNUM,ND) - needed to calculate new collision cross sections for Hminus     
      CALL RADIO(N,ENLTE,TL,WEIGHT,NCHARG,EION,ELEVEL,EINST,
     $           RRATE,XLAMBDA,FWEIGHT,XJC,NF,L,XJL,ND,LASTIND,SIGMAKI,NOM)

C***  LOOP FOR EACH ELEMENT  -------------------------------------------
 
      DO 11 NA = 1, NATOM

      NFIRNA = NFIRST(NA)
      NLANA =  NLAST(NA)

      NDELTA = NLANA - NFIRNA
      NDELP1 = NDELTA + 1

!      print*, 'nltepop: ', N, NDELP1; stop

C***  ADDING RADIATIVE AND COLLISIONAL RATES INTO RATCO (=RATE COEFFICIENTS)
      DO 4 I=NFIRNA,NLANA
      ISHIFT=I-NFIRNA+1
      DO 4 J=NFIRNA,NLANA

      JSHIFT=J-NFIRNA+1
      RATCO(ISHIFT,JSHIFT)=CRATE(I,J)+RRATE(I,J)

    4 CONTINUE

C***  DIAGONAL ELEMENTS : - SUM OF THE ROW (I.E. OVER COLUMN INDEX)
      DO 1 I = 1, NDELTA
      SUM=.0
      DO 2 J = 1, NDELP1
    2 SUM = SUM + RATCO(I, J)
    1 RATCO(I, I) = -SUM
     
C***  LAST COLUMN : NUMBER CONSERVATION
      do I = 1, NDELP1; RATCO(I, NDELP1) = 1.0d0; enddo

C***  INVERSION OF RATE COEFFICIENT MATRIX RATCO

!      do fi = nfirna, nlana

!         do si = nfirna, nlana

!            write(*, '(A,3(2x,i4),3(2x,e15.7))'), 'ratco: ', NA, fi, si, ratco(fi, si),
!     $                                             crate(fi, si), rrate(fi, si)

!         enddo

!      enddo

      write(*, *) 'nltepop flag!!!'

      CALL INV(NDELP1, RATCO(1 : ndelp1, 1 : ndelp1))
	
C***  POPULATION NUMBERS EN(J) = LAST ROW OF INVERSE MATRIX
      AB = ABXYZ(NA)

      do J = NFIRNA, NLANA

         EN(J) = AB * RATCO(NDELP1, J - NFIRNA + 1)

      enddo

   11 CONTINUE

      RETURN

      END SUBROUTINE

      END MODULE
