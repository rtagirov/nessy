      MODULE MOD_CSTABREAD
      contains
      SUBROUTINE cstabread(SIGMA,WAVENUM,LOW,N,WAVARR,SIGARR,NDIM,NFDIM) 
C*********************************************************
C***  CALLED BY PHOTOCS_M
CMH   written by Margit Haberreiter to read cross section tables
C**************************
C********************************************************************* 
      USE MOD_INTPL
      USE UTILS
      IMPLICIT NONE
      INTEGER, intent(in):: LOW,NDIM,NFDIM,N
      real*8, intent(in)::WAVENUM,SIGARR(NDIM,NFDIM),WAVARR(NDIM,NFDIM)
      real*8, intent(out)::SIGMA
      integer :: nup,nlow,k,ii
      real*8 :: WAV1,SIG1,WAV2,SIG2
      SIGMA=0.

      !* find the index k in the reverse-sorted array WAVARR(low,:)
      !* such that WAVARR(low,k)>= wavenum >= WAVARR(low,k+1)
      !* using a binary search tree. Runtime is O(ln2(NFDIM))
      nlow=1
      nup=NFDIM
      BIN_SEARCH:DO ii=1,999 !* safety guard if list is unorderd
        k=(nup+nlow)/2   
        IF (wavenum > WAVARR(low,k)) THEN 
          nup=k  
        ELSE
          nlow=k
        END IF
        IF (nlow+1 >= nup) exit
      ENDDO BIN_SEARCH
      k=nlow
      !* check for correct result and warn(1) if an overflow occured
      call assert(ii<900,'cstabread: BIN_SEARCH: WAVARR unsorted?',1)
      call assert((WAVARR(low,k)  .ge. wavenum) .and. 
     &            (WAVARR(low,k+1) .le. wavenum),
     &            'cstabread: wavenum does not fit, Index not correct')
      !* End Binary search
      
      WAV1=WAVARR(low,k)
      SIG1=SIGARR(low,K)
      WAV2=WAVARR(low,k+1)
      SIG2=SIGARR(low,k+1)
      !if (sig1 .eq. 0.) then ! IF FIRST CROSS SECTION IS ZERO
      !  sigma=sig2*1.d-18
      !  goto 1
      !end if
      !if (sig2 .eq. 0.) then ! IF SECOND CROSS SECTION IS ZERO
      !  sigma= SIG1*1.d-18
      !  goto 1
      !end if
      call intpl(SIGMA,WAVENUM,WAV1,WAV2,SIG1,SIG2)
      SIGMA=SIGMA*1.d-18     !sigma=sigma/1d18
c     SIGMA=0.
      END SUBROUTINE
      END MODULE

