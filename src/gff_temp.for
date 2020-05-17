      MODULE MOD_GFF_TEMP

      CONTAINS

!     RINAT TAGIROV:
!     This subroutine was copied from the paper by Christian Janicki, 1990:
!     "A computer program for the free-free and bound-free Gaunt factors of Rydberg systems",
!     Computer Physics Communications, Vol. 60, Issue 3, p. 281-296


      SUBROUTINE GFF_TEMP(IQ, Z, hv, Te, gmoy)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION QX(10), QW(10), QX4(4), QW4(4), QX10(10), QW10(10)

      DATA QX4, QW4 /0.3225478D0, 1.745761D0, 4.536620D0, 9.395072D0,
     $ 0.6031541D0, 0.3574187D0, 3.8887952D-02, 5.3929421D-04/

      DATA QX10, QW10 /0.1377939D0, 0.7294555D0, 1.808344D0, 3.401433D0,
     $ 5.552495D0, 8.330153D0, 11.84379D0, 16.27926D0, 21.99659D0, 
     $ 29.92069D0, 0.3084415D0, 0.4011198D0, 0.2180681D0, 6.2087435D-02,
     $ 9.5015280D-03, 7.5300789D-04, 2.8259257D-05, 4.2493159D-07,
     $ 1.8395685D-09, 9.9118652D-13/

      IF (IQ .EQ. 0) THEN
         
         NQ = 4
         
         DO I = 1, 4

            QX(I) = QX4(I)
            QW(I) = QW4(I)

         ENDDO

      ELSE

         NQ = 10

         DO I = 1, 10

            QX(I) = QX10(I)
            QW(I) = QW10(I)

         ENDDO

      ENDIF

      gmoy = 0

      DO I = 1, NQ

         Ef = Te * QX(I) + hv

         CALL GAUNT_FF(Z, Ef, hv, gff)

         gmoy = gmoy + gff * QW(I)

      ENDDO

      RETURN

      END SUBROUTINE


      SUBROUTINE GAUNT_FF(Z, Ef, hv, gff)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DATA RY /13.60580436D0/ !Rydberg energy (eV)

      Ei = Ef - hv

      ETA_i = DSQRT(RY * Z**2 / Ei)
      ETA_f = DSQRT(RY * Z**2 / Ef)

      CALL F_COUL(0, ETA_i, ETA_f, F0)
      CALL F_COUL(1, ETA_i, ETA_f, F1)

      gff = (1.1027D0 * F0) * (F0 * (ETA_i / ETA_f + ETA_f
     $ / ETA_i + 2.D0 * ETA_i * ETA_f) - 2.D0 * F1 *
     $ DSQRT((1.D0 + ETA_i**2) * (1.D0 + ETA_f**2)))

      RETURN

      END SUBROUTINE


      SUBROUTINE F_COUL(L, ETA_i, ETA_f, FL)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DATA PI /3.1415926536D0/

      COMPLEX*16 C1, C2, C3

      CALL GL_FF(L, ETA_i, ETA_f, G1)

      C1 = DCMPLX(1.D0 + L, ETA_i)
      C2 = DCMPLX(1.D0 + L, ETA_f)
      C3 = DCMPLX(2.D0 * L + 2.D0, 0.D0)

      ARG = 0.5D0 * PI * ABS(ETA_i - ETA_f) + GALNDC(C1) + GALNDC(C2)
     $      - GALNDC(C3)

      FL = 0.25D0 * G1 * DEXP(ARG) * (4.D0 * ETA_i * ETA_f /
     $     (ETA_f - ETA_i)**2)**(L + 1)

      RETURN

      END SUBROUTINE


      SUBROUTINE GL_FF(L, ETA_i, ETA_f, G1)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMPLEX*16 CX0, CX1, CX2, CX3, CX4

      DIMENSION Y(200), Z(200)

      X = -4.D0 * ETA_i * ETA_f / (ETA_i - ETA_f)**2

      IF (DABS(X) .LE. 1.618) THEN

         TX = X / (X - 1.D0)

         A0 = 1.D0
         A1 = ALPHA0(1, L, ETA_i, ETA_f)

         LIMI = -37 / DLOG10(DABS(TX)) + 1

         Y(LIMI + 2) = 0
         Y(LIMI + 1) = 0

         DO I = 1, LIMI

            N = LIMI - I + 1
            Y(N) = ALPHA0(N + 1, L, ETA_i, ETA_f) * Y(N + 1) +
     $             BETA0(N + 2, L, ETA_i, ETA_f) * Y(N + 2) + TX**N

         ENDDO

         G1 = A0 + A0 * BETA0(2, L, ETA_i, ETA_f) * Y(2) + A1 * Y(1)

      ELSE

         YM = -1.D0 / X

         CX1 = DCMPLX(0.D0, ETA_f - ETA_i)
         CX2 = DCMPLX(L + 1.D0, -ETA_i)
         CX3 = DCMPLX(L + 1.D0,  ETA_f)
         CX4 = DCMPLX(2.D0 * L + 2.D0, 0.D0)

         CX0 = GALNDC(CX4) + GALNDC(CX1) - GALNDC(CX2) - GALNDC(CX3)
         CX0 = CDEXP(CX0)

         a0 =  2.d0 * DREAL(CX0)
         b0 = -2.d0 * DIMAG(CX0)

         LIMI = -37 / DLOG10(DABS(YM)) + 1

         Y(LIMI + 2) = 0
         Y(LIMI + 1) = 0
         Z(LIMI + 2) = 0
         Z(LIMI + 1) = 0

         DO I = 1, LIMI + 1

            N = LIMI - I + 1

            AN1 = ALPHA1(N + 1, L, ETA_i, ETA_f)
            BN2 =  BETA1(N + 2, L, ETA_i, ETA_f)

            GN1 =  GAMM1(N + 1, L, ETA_i, ETA_f)
            DN2 = DELTA1(N + 2, L, ETA_i, ETA_f)

            IF (N .NE. 0) THEN

               Y(N) = AN1 * Y(N + 1) + BN2 * Y(N + 2) - GN1 * Z(N + 1)
     $                - DN2 * Z(N + 2) + YM**N

               Z(N) = AN1 * Z(N + 1) + BN2 * Z(N + 2) + GN1 * Y(N + 1)
     $                + DN2 * Y(N + 2) + YM**N

            ELSE

               Y0 = AN1 * Y(1) + BN2 * Y(2) - GN1 * Z(1) - DN2 * Z(2)
     $              + 1.D0

               Z0 = AN1 * Z(1) + BN2 * Z(2) + GN1 * Y(1) + DN2 * Y(2)
     $              + 1.D0

            ENDIF

         ENDDO

         AYM = 0.5D0 * ((a0 - b0) * Y0 + (a0 + b0) * Z0)
         BYM = 0.5D0 * ((a0 + b0) * Y0 - (a0 - b0) * Z0)

         ARG = 0.5D0 * (ETA_i - ETA_f) * DLOG(YM)

         G1 = YM**(L + 1) * (AYM * DCOS(ARG) + BYM * DSIN(ARG))

      ENDIF

      RETURN

      END SUBROUTINE


      DOUBLE PRECISION FUNCTION ALPHA0(N, L, ETA_i, ETA_f)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      E = ETA_i * ETA_f

      AN = N

      ALPHA0 = ((AN - 1) * (2 * AN + 2 * L - 1) - (L + 1)**2 + E) /
     $         (AN * (2 * L + 1 + AN))

      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION BETA0(N, L, ETA_i, ETA_f)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      E2 = (0.5D0 * (ETA_i - ETA_f))**2

      AN = N

      E = ETA_i * ETA_f

      BETA0 = -((AN - 2)**2 + E2 + E) / (AN * (2 * L + 1 + AN))

      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION ALPHA1(N, L, ETA_i, ETA_f)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ALAMBDA = 0.5D0 * (ETA_i - ETA_f)

      AALPHA = L * (L + 1) + 2 * ALAMBDA**2 + ETA_i * ETA_f

      ABETA = L * (L + 1) + ALAMBDA**2

      AN = N

      ALPHA1 = (AN * AALPHA - AN * (AN - 1) * (2 * AN - 1) - 2.D0 *
     $         ALAMBDA**2 * (4 * AN - 3)) / 
     $         (AN * (AN**2 + 4.D0 * ALAMBDA**2))

      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION BETA1(N, L, ETA_i, ETA_f)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!      ALAMBDA = 0.5D0 * (ETA_i _ ETA_f)
      ALAMBDA = 0.5D0 * (ETA_i - ETA_f)

      AALPHA = L * (L + 1) + 2.D0 * ALAMBDA**2 + ETA_i * ETA_f

      ABETA = L * (L + 1) + ALAMBDA**2

      AN = N

      BETA1 = (AN * ABETA - AN * (AN - 1) * (AN - 2) - 2 * ALAMBDA**2 *
     $        (2 * AN - 3)) / (AN * (AN**2 + 4.D0 * ALAMBDA**2))

      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION GAMM1(N, L, ETA_i, ETA_f)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ALAMBDA = 0.5D0 * (ETA_i - ETA_f)

      AALPHA = L * (L + 1) + 2 * ALAMBDA**2 + ETA_i * ETA_f

      AN = N

      GAMM1 = -ALAMBDA * (3 * AN - 2 + 2.D0 * AALPHA) /
     $ (AN * (AN**2 + 4.D0 * ALAMBDA**2))

      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION DELTA1(N, L, ETA_i, ETA_f)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ALAMBDA = 0.5D0 * (ETA_i - ETA_f)

      AALPHA = L * (L + 1) + 2 * ALAMBDA**2 + ETA_i * ETA_f

      ABETA = L * (L + 1) + ALAMBDA**2

      AN = N

      DELTA1 = -ALAMBDA * (3 * AN - 4 + 2.D0 * ABETA)
     $ / (AN * (AN**2 + 4.D0 * ALAMBDA**2))

      RETURN

      END FUNCTION


      DOUBLE COMPLEX FUNCTION GALNDC(XX)

      COMPLEX*16 XX, X, TMP, SER

      REAL*8 COF(6), STP, HALF, ONE, FPF, ZERO

      DATA COF, STP /76.18009173D0, -86.50532033D0, 24.01409822D0,
     $ -1.231739516D0, .120858003D-2, -.536382D-5, 2.50662827465D0/

      DATA HALF, ONE, FPF, ZERO /0.5D0, 1.0D0, 5.5D0, 0.0D0/

      X = XX - DCMPLX(ONE, ZERO)

      TMP = X + DCMPLX(FPF, ZERO)

      TMP = (X + DCMPLX(HALF, ZERO)) * CDLOG(TMP) - TMP

      SER = DCMPLX(ONE, ZERO)

      DO J = 1, 6

         X =   X + DCMPLX(ONE, ZERO)
         SER = SER + DCMPLX(COF(J), ZERO) / X

      ENDDO

      GALNDC = TMP + CDLOG(STP * SER)

      RETURN

      END FUNCTION


      END MODULE
