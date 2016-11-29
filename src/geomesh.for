      MODULE MOD_GEOMESH

      CONTAINS

      SUBROUTINE GEOMESH(RADIUS,ENTOT,T,P,Z,RMAX,RSTAR,AMU,ATMEAN,ND,NDDIM,NP,NPDIM)

      USE MOD_PGRID
      USE MOD_RGRIDM

!     THIS SUBROUTINE GENERATES THE GEOMETRICAL POINT MESH IN RADIUS, P AND Z
!     P and Z mesh is needed for the ray-by-ray solution of the radiative transfer equation in spherical symmetry

      IMPLICIT REAL*8(A - H, O - Z)

      INTEGER, INTENT(OUT) ::    ND, NP
      REAL*8,  INTENT(INOUT) ::  RADIUS, entot(:), T(:), P, Z(NDDIM*NPDIM)
      INTEGER, INTENT(IN) ::     NDDIM, NPDIM
      REAL*8,  INTENT(IN) ::     RMAX, RSTAR, AMU, ATMEAN

      DIMENSION                  RADIUS(NDDIM), P(NPDIM)
      
      CALL RGRIDM(RADIUS, ENTOT, T, RMAX, RSTAR, AMU, ATMEAN, NDDIM, ND)

      CALL PGRID(NPDIM, NP, ND, RADIUS, P)

      do i=1,NDDIM*NPDIM
          z(i)=0.0D0
      enddo

!      print*, 'geomesh 0:', nddim, npdim, size(z)

!      stop

      DO L=1,ND
         RR=RADIUS(L)*RADIUS(L)
        JMAX=NP+1-L
        DO 2 J=1,JMAX
            PJ=P(J)
            PJPJ=PJ*PJ
            I=(J-1)*ND+L   
            Z(I)=SQRT(RR-PJPJ)

            write(*,'(A,2x,3(i4,2x),3(e15.7,2x))') 'geomesh:', l, j, i, rr, pjpj, z(i)
!            write(*,'(A8,2x,3(i4,2x),A,1(2x,d15.7))') 'geomesh:', l, j, i, '                                ', z(i)

    2   enddo
    1 enddo

!      stop

      RETURN

      END SUBROUTINE

      END MODULE
