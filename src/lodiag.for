      MODULE MOD_CALCLAMBDAS

      CONTAINS

      !*** micha: Caclulate the accelerated Lambda Iterator.
      !PURE FUNCTION CALCLAMBDAS(OPA,ETA,THOMSON,R,EDDI,ND)

      FUNCTION CALCLAMBDAS(OPA, R, EDDI, ND)
      ! called by como.for
      ! calculate the LambdaStar operator
      IMPLICIT NONE
      integer, intent(in):: ND                 !Number of Depth/Frequency points, index
      real*8, dimension(ND) :: CALCLAMBDAS     !Approximated Lambda Operator
      real*8, intent(in), dimension(ND)  :: OPA, R  !OPA: Opacity, R: Distance from surface
      !real*8,intent(in),dimension(ND)  :: THOMSON !relative Electron Scattering factor
      real*8, intent(in), dimension(3,ND):: EDDI    !Eddi,sphericity and flux factors, see comment below
      real*8, dimension(ND)  :: A,B, C             !calculate the Lambda Operator
      ! real*8 :: DA, DC, DR,DRM, FL, FM, FP, QL, QM, QP, H, HPLUS     !helper variables
      integer:: L             !index
      real*8 :: dXm, dXl,dX
      !**** EDDI(1,L) is the eddington factor  f_nu = K_nu / J_nu = 2nd order
      !***            acceleration / mean Intensity
      !***  EDDI(2,L) is the sphericity factor q see Mihalas 2nd Edition, 
      !***            7-6 Extended Atmospheres, eq. 7-193 p. 251
      !***  EDDI(3,L) is the eddington factor h / j  (only at the boundaries)
      !***  EDDI(3,ND-1) is the outward flux hplus at the inner boundary
      !***            See Mihalas "Stellar Atmospheres", 2nd Edition, p. 251:
      !***            7-196, 7-198, 7-200) Source function = ETA/OPA
      !*** A,B,C give a triangular matrix of the following form:
      !***                B  -C
      !***               -A   B  -C
      !***                   -A   B  -C
      !***                       -A   B
      !### For the Boundary Conditions see Auer "Improved Boundary Conditions 
      !### for the Feautrier Method" 1967.
      !***=== OUTER BOUNDARY ================================================
      dXl = (OPA(2)+OPA(1))/2d0*(R(2)-R(1))!*RSTAR
      dX = dXl*EDDI(2,1)  ! EDDI(2,L) = q(l)
      A(1) =  0d0
      C(1) =  2./(dX*dXl)
      !***===  NON-BOUNDARY POINTS ==========================================
      !*** from G.B. Rybicki and D.G. Hummer, "An accelerated lambda
      !***   iteration method for multilevel radiative transfer", A&A 254,
      !***   171-181 (1991), [1991A&A...245..171R]
      !*** and R.J. Rutten, "Radiative Transfer in Stellar Atmospheres",
      !*** Intitute of Theoretical Astrophysics Oslo, May 8, 2003 (sphericity
      !*** factor) the complete model is described in Auer, 1967.
      NonBoundary: DO L=2,ND-1
        !*** See Mihalas, Equation 7-196
        !*** Define the Operator "2nd Derivative" (e.g. Rutten, 5.23-25 )
        dXm = dXl
        dXl = (OPA(L+1)+OPA(L))/2d0*(R(l+1)-R(l))!*RSTAR
        dX = (dXm+dXl)/2d0*EDDI(2,L)**2  ! EDDI(2,L) = q(l)
        A(L) =  2d0/(dX*dXm)
        C(L) =  2d0/(dX*dXl)
      ENDDO NonBoundary
      !***===  INNER BOUNDARY (Auer, 1967) ==================================
      !dXm=dXl
      dXl = (OPA(ND)+OPA(ND-1))/2d0*(R(ND)-R(ND-1))!*RSTAR
      !dX = (dXm+dXl)/2d0*EDDI(2,ND)  ! EDDI(2,L) = q(l)
      dX = (dXl)*EDDI(2,ND)  ! EDDI(2,L) = q(l)
      A(ND) =  2d0/(dX*dX)
      C(ND) = 0d0
      B =  A+C+1d0

      CALCLAMBDAS = INVTRIDIAG(A, B, C)

      END FUNCTION


      FUNCTION INVTRIDIAG(A, B, C)

      !***  ==== DIAGONAL ELEMENTS OF THE INVERSE OF A TRIDIAGONAL MATRIX ====
      !*  W(LMAX)=RL*RL*(ETA(LMAX)+two*QL*X*HPLUS/DX)
      !*  A,B,C for T=tri(A,B,C) are calculated, now calculate
      !*  \Lambda* = diag(T^(-1)) according to G.B. Rybicki and D.G. Hummer
      !*  in Astron Astrophysics, 245, 171--181 (1991) (Appendix B, p.180)
      !*  A(1) is one below the top of the matrix
      !*  C(LMAX) is one above the bottom of the matrix
      !*  D(0), E(LMAX+1) are the initial help-vectors for xL
      !*  D(0) = 0d0;z
      !*   Matrix:
      !* /  B1  C1                   \
      !* |  A2  B2  C2               |
      !* |      A3  B3  C3           |
      !* |                           |
      !* |                 AM BM CM  |
      !* \                    AN BN  /, M=N-1, A1=CN=0

      IMPLICIT NONE

      REAL*8, INTENT(IN), DIMENSION(:) :: A, B, C

      INTEGER :: L, LMAX

      REAL*8 :: INVTRIDIAG(size(A))
      REAL*8 :: E(2 : size(A) + 1)  ! for the inverse
      REAL*8 :: DL, DLM             ! DL -> DLM at the end of each iteration

      LMAX = size(A)

      DLM = 0.0d0

      E(LMAX + 1) = 0.0d0
      E(LMAX) = A(LMAX) / B(LMAX)

      DO L = LMAX - 1, 2, -1; E(L) = A(L) / (B(L) - C(L) * E(L + 1)); ENDDO

      DO L = 1, LMAX

         DL = C(L) / (B(L) - A(L) * DLM);

         INVTRIDIAG(L) = 1.0d0 / ((1.0d0 - DL * E(L + 1)) * (B(L) - A(L) * DLM))

         DLM = DL

!        RINAT TAGIROV:
!        Check if INVTRIDIAG(L) is a NaN, in which case the coefficients of the
!        Feautrier matrix are infinite which in turn means that we're in the optically
!        thin case => no need for acceleration.

         if (isnan(INVTRIDIAG(L))) INVTRIDIAG(L) = 0.0D0

      ENDDO

      END FUNCTION

      END MODULE
