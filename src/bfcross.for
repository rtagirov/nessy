      MODULE MOD_BFCROSS

      contains

      SUBROUTINE BFCROSS(SIGMAKI,NF,N,NCHARG,ELEVEL,EION,EINST,XLAMBDA,ALPHA,SEXPO,AGAUNT,NOM,WAVARR,SIGARR)

!     THIS ROUTINE PREPARES AN ARRAY SIGMAKI WITH THE BOUND-FREE CROSS SECTIONS
!     (IN CM**2) TO AVOID UNNECCESSARY MULTIPLE CALCULATIONS

      USE MOD_PHOTOCS

      implicit none

      integer, intent(in) ::                   NF, N
      integer, dimension(N), intent(in) ::     NOM, NCHARG

      real*8, intent(in) ::                    EION(N), ELEVEL(N), EINST(N, N)
      real*8, intent(in), dimension(N) ::      ALPHA, SEXPO
      character*8, intent(in) ::               agaunt(N)
      real*8, intent(in), dimension(N, NF) ::  WAVARR, SIGARR
      real*8, intent(in), dimension(NF)   ::   XLAMBDA

      real*8, intent(out), dimension(NF, N) :: SIGMAKI

      integer ::                               NUP, K, LOW, I
      real*8 ::                                SIGMATH, EDGE, WAVENUM, sigma

      do 8 NUP = 2, N

         do 8 LOW = 1, NUP - 1

      !***  LEVELS MUST BELONG TO THE SAME ELEMENT

            if (NOM(LOW) .ne. NOM(NUP)) goto 8

      !***  CHARGE DIFFERENCE MUST BE 1

            if (NCHARG(NUP) .ne. NCHARG(LOW)+1 ) goto 8

      !***  UPPER LEVEL MUST BE GROUND STATE

            if (NCHARG(NUP) .eq. NCHARG(NUP-1)) goto 8

      !***  EINST = THRESHOLD CROSS SECTION IN 10**-18 CM**2

            SIGMATH = EINST(LOW, NUP) * 1.d-18

      !***  EDGE = THRESHOLD ENERGY IN KAYSER *****

            EDGE = EION(LOW) - ELEVEL(LOW)

            do K = 1, NF

               WAVENUM = 1.E8 / XLAMBDA(K)

               if (WAVENUM .lt. EDGE) then

                  SIGMAKI(K, LOW) = 0.0d0

               else

                  call PHOTOCS(SIGMA,SIGMATH,EDGE,WAVENUM,ALPHA,SEXPO,AGAUNT,LOW,
     $                         WAVARR(1 : N, 1 : NF), SIGARR(1 : N, 1 : NF), N, NF)

                  sigmaki(k, low) = sigma

               endif

            enddo

    8 continue

      return

      END SUBROUTINE

      END MODULE
