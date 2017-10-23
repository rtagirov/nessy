      MODULE MOD_CSTABREAD

      contains

      SUBROUTINE cstabread(SIGMA, WAVENUM, LOW, WAVARR, SIGARR, N, NF)

!     CALLED BY PHOTOCS
!     written by Margit Haberreiter to read cross section tables

      USE MOD_INTPL
      USE UTILS

      IMPLICIT NONE

      INTEGER, intent(in) :: LOW, N, NF

      real*8, intent(in) :: WAVENUM, SIGARR(N, NF), WAVARR(N, NF)

      real*8, intent(out) :: SIGMA

      integer :: nup, nlow, k, ii

      real*8 :: WAV1, SIG1, WAV2, SIG2

      SIGMA = 0.

      ! find the index k in the reverse-sorted array WAVARR(low,:)
      ! such that WAVARR(low,k)>= wavenum >= WAVARR(low,k+1)
      ! using a binary search tree. Runtime is O(ln2(NF))
      nlow = 1

      nup = NF

      BIN_SEARCH: DO ii = 1, 999 ! safety guard if list is unorderd

        k = (nup + nlow) / 2

        IF (wavenum .gt. WAVARR(low, k)) THEN

           nup = k

        ELSE

           nlow = k

        ENDIF

        IF (nlow + 1 .ge. nup) exit

      ENDDO BIN_SEARCH

      k = nlow

      !* check for correct result and warn(1) if an overflow occured
      call assert(ii .lt. 900,'cstabread: BIN_SEARCH: WAVARR unsorted?',1)
      call assert((WAVARR(low, k) .ge. wavenum) .and. (WAVARR(low, k + 1) .le. wavenum),
     &            'cstabread: wavenum does not fit, Index not correct')
      !* End Binary search

!      if (.not. ((WAVARR(low, k) .ge. wavenum) .and. (WAVARR(low, k + 1) .le. wavenum))) then

!         write(*, '(A,2x,I4,2x,I4,3(2x,F10.4))'), 'cstabread: ', low, k, 1.0d8 / wavenum, 
!     $                                           1.0d8 / wavarr(low, k), 1.0d8 / wavarr(low, k + 1)

!      endif
      
      WAV1 = WAVARR(low, k)
      SIG1 = SIGARR(low, K)
      WAV2 = WAVARR(low, k + 1)
      SIG2 = SIGARR(low, k + 1)

      call intpl(SIGMA,WAVENUM,WAV1,WAV2,SIG1,SIG2)

      SIGMA = SIGMA * 1.d-18

      END SUBROUTINE

      END MODULE
