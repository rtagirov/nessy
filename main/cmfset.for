      MODULE MOD_CMFSET

      CONTAINS

      SUBROUTINE CMFSET(PHIK,Z,ND,LMAX,TA,TB,TC,UB,VA,VB,GA,H,S,OPA,OPAL,
     $                  ETA,ETAL,PP,BCORE,DBDR,XIMINUS,DXI)

!     THIS SUBROUTINE IS TO SET UP THE ARRAY ELEMENTS FOR THE CMF FORMALISM

!     RINAT TAGIROV:
!     It follows formalism described in Mihalas, Kunasz & Hummer, 1975, 202: 465 - 489
!     "Solution of the comoving-frame equation of transfer in spherically symmetric flows. 
!     I.Computational method for equivalent-two-level-atom source functions" (MKH),
!     but with the source function of general form, i.e the following equation is solved:
!     T_{k,j}u_{k,j} + U_{k,j}u_{k-1,j} + V_{k,j}v_{k-1,j} = S_{k,j}
!     instead of Eq. (2.39) in MKH.
!     Specifically see their Eq. (2.31) and Appendix C.
!     TA, TB, TC are the diagonals of matrix T_{k,j} (tridiagonal).
!     VA, VB are the diagonals of matrix V_{k,j} (bidiagonal).
!     UB is the matrix U_{k,j} (diagonal).
!     H is the matrix H_{k,j} (diagonal).
!     GA is matrix G_{k,j} which is bidiagonal but its diagonal element 
!     is minus one multiplied by the non-diagonal one, so only GA array is needed here.
!     The signs do not match their formulation because apparently the routine is written
!     for the aforementioned matrix equation having slightly different 
!     form with some signs taken out of the matrix elements
!     and put directly into the equation.

      implicit real*8(a-h,o-z)

      parameter (half=0.5d0, one=1.0d0, two=2.0d0)

      DIMENSION S(ND),OPA(ND),OPAL(ND),ETA(ND),ETAL(ND),VA(ND),VB(ND)
      DIMENSION UB(ND),GA(ND),H(ND)
      DIMENSION PP(ND),Z(ND)

      DIMENSION TA(LMAX), TB(LMAX), TC(LMAX)

      LZ = LMAX - 1
     
!     OUTER BOUNDARY CONDITION - FIRST ORDER
      AK=OPA(1)+PHIK*OPAL(1)
      AZ=half*(AK+OPA(2)+PHIK*OPAL(2))
      TAUZ=Z(1)-Z(2)
      DX=PP(1)/AK
      S(1)=XIMINUS+DX*DXI
      TC(1)=one/AK/TAUZ
      TB(1)=TC(1)+DX+one
      UB(1)=DX
      VB(1)=0.0

!     FOR G AND H THE MATRIX ELEMENT S ARE NOT DIFFERENT FROM INNER POINTS
      DTZM=one/AZ/TAUZ
      DXZM=(PP(1)+PP(2))/AZ/two
      DAZM=DTZM/(one+DXZM)
      DBZM=DXZM/(one+DXZM)
      GA(1)=-DAZM
      H(1)=DBZM
      IF(LZ.LT.2) GOTO 2
     
!     NON - BOUNDARY POINTS
      DO 1 L=2,LZ
      AK=OPA(L)+PHIK*OPAL(L)
      EK=ETA(L)+PHIK*ETAL(L)
      S(L)=EK/AK
      AZ=half*(AK+OPA(L+1)+PHIK*OPAL(L+1))
      AZM=half*(AK+OPA(L-1)+PHIK*OPAL(L-1))
      TAU=half*(Z(L-1)-Z(L+1))
      TAUZ=Z(L)-Z(L+1)
      TAUZM=Z(L-1)-Z(L)
      DT=one/AK/TAU
      DTZ=one/(AZ*TAUZ)
      DTZM=one/(AZM*TAUZM)
      DX=PP(L)/AK
      DXZ=(PP(L)+PP(L+1))*half/AZ
      DXZM=(PP(L)+PP(L-1))*half/AZM
      DAZ=DTZ/(one+DXZ)
      DAZM=DTZM/(one+DXZM)
      DBZ=DXZ/(one+DXZ)
      DBZM=DXZM/(one+DXZM)

      TA(L)=DT*DAZM
      TC(L)=DT*DAZ
      TB(L)=TA(L)+TC(L)+DX+one

      UB(L)=DX

      VA(L)=-DT*DBZM
      VB(L)=DT*DBZ

      GA(L)=-DAZ

      H(L)=DBZ

    1 CONTINUE
     
    2 L = LMAX
      IF (LMAX .LT. ND) GOTO 4
     
!     INNER BOUNDARY CONDITION (CORE RAYS) - ONLY TO FIRST ORDER
      AK=OPA(ND)+PHIK*OPAL(ND)
      S(ND)=BCORE+DBDR*Z(ND)/AK
      TAUZ=Z(ND-1)-Z(ND)
      DT=one/TAUZ/AK
      DX=PP(L)/AK
      TA(L)=DT
      TB(L)=DT+DX+one
      UB(L)=DX
      VA(L)=0.0d0

      RETURN
     
!     INNER BOUNDARY CONDITION (NON-CORE RAYS) - SECOND ORDER
    4 AK=OPA(LMAX)+PHIK*OPAL(LMAX)
      EK=ETA(L)+PHIK*ETAL(L)
      S(L)=EK/AK
      TAUZ=Z(LZ)
      DT=one/AK/TAUZ
      DX=PP(L)/AK
      DA=DT/(one+DX)
      DB=DX/(one+DX)
      TA(L)=two*DT*DA
      TB(L)=TA(L)+DX+one
      UB(L)=DX
      VA(L)=-two*DT*DB

      RETURN

      END SUBROUTINE

      END MODULE
