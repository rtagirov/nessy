      module auxfioss

      contains

      function read_param(par_name) result(par)

      character (len = 2), intent(in) :: par_name

      integer :: par

      character (len = 2) :: ftc

      integer :: fu

      fu = 13745
     
      open(unit = fu, file = 'MODFILE', action = 'read')

      ftc = '  '

      do while (ftc .ne. par_name)

         read(fu, '(A2)') ftc

      enddo

      read(fu, *) par

      close(fu)

      return

      end function read_param

      SUBROUTINE READ_NLTETRAPOP(TI, TPL, TPU, TDL, TDU)

      USE COMMON_BLOCK
      USE FILE_OPERATIONS

      IMPLICIT NONE

      INTEGER, DIMENSION(NTC), INTENT(IN) ::               TI

      REAL*8, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: TPL, TPU, TDL, TDU

      CHARACTER (LEN = 32) :: FUDGE

      INTEGER :: I, J, L, DEPTH_POINTS_NUM

      REAL*8 :: PL, PU, DL, DU

      depth_points_num = num_of_lines(atm_mod_file)

      ALLOCATE(TPL(DEPTH_POINTS_NUM, NTC))
      ALLOCATE(TPU(DEPTH_POINTS_NUM, NTC))
      ALLOCATE(TDL(DEPTH_POINTS_NUM, NTC))
      ALLOCATE(TDU(DEPTH_POINTS_NUM, NTC))

      TPL(:, :) = 0.0D0; TPU(:, :) = 0.0D0; TDL(:, :) = 0.0D0; TDU(:, :) = 0.0D0
 
      DO I = 1, NTC

         OPEN(773, FILE = NTP_FILE)

         IF (TI(I) .NE. 1) THEN

             DO J = 1, (TI(I) - 1) * DEPTH_POINTS_NUM; READ(773, '(A66)') FUDGE; ENDDO

         ENDIF

         DO L = 1, DEPTH_POINTS_NUM

            READ(773, '(E15.7,3(2x,E15.7))') PL, PU, DL, DU

            TPL(L, I) = PL; TPU(L, I) = PU; TDL(L, I) = DL; TDU(L, I) = DU

            WRITE(*, '(A,4(1x,E15.7))') 'TEST READ NLTE:', TPL(L, I), TPU(L, I), TDL(L, I), TDU(L, I)

         ENDDO

         CLOSE(773)

      ENDDO

      END SUBROUTINE


      SUBROUTINE READ_NLTE_TRA(NFE, NUMTR, EID, CH, WL, WU, AUPLOW, WAVTR)

      USE FILE_OPERATIONS

      IMPLICIT NONE

      LOGICAL, INTENT(OUT) ::                            NFE

      INTEGER, INTENT(OUT) ::                            NUMTR

      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: EID, WL, WU, CH

      REAL*8, ALLOCATABLE, DIMENSION(:), INTENT(OUT) ::  WAVTR, AUPLOW

      INTEGER ::                                         I, EI, GU, GL, CHARG

      CHARACTER (LEN = 10) ::                            LL, LU

      REAL*8 ::                                          WAVT, AUL

      INQUIRE(FILE = NTW_FILE, EXIST = NFE)

      IF (NFE .EQ. .FALSE.) RETURN

      NUMTR = NUM_OF_LINES(NTW_FILE)

      ALLOCATE(AUPLOW(NUMTR))
      ALLOCATE(WAVTR(NUMTR))
      ALLOCATE(EID(NUMTR))    ! ELEMENT IDENTIFICATION, I.E. ITS NUMBER IN THE PERIODICAL TABLE
      ALLOCATE(CH(NUMTR))
      ALLOCATE(WL(NUMTR))
      ALLOCATE(WU(NUMTR))

      OPEN(132, FILE = NTW_FILE)

      DO I = 1, NUMTR

         READ(132, '(2(A10,2x),I2,2x,I2,2(2x,I3),2(2x,E15.7))') LL, LU, EI, CHARG, GL, GU, AUL, WAVT

         WAVTR(I) = WAVT; EID(I) = EI; CH(I) = CHARG; WL(I) = GL; WU(I) = GU; AUPLOW(I) = AUL

         WRITE(*, '(A,I2,3(2x,I3),2(2x,E15.7))') 'NLTE TRA:', EID(I), CH(I), WL(I), WU(I), AUPLOW(I), WAVTR(I)

      ENDDO

      END SUBROUTINE

      SUBROUTINE COMPXJ9(ND, NP, JP, NVOPA, DINT, XJK, WLK, IBACK)

C***  Compute mean Intensity (Moment 0)
C***  add for each impact parameter the intensities DINT
C***  WLK are the weigths given by the geometry of the grid
C***  This routine is called within the impact parameter loop
c     modified: 20.1.01  do loop only up tp LMAX

      implicit real*8(a - h, o - z)

      DIMENSION DINT(:,:),XJK(:,:),WLK(ND,NP)
      DIMENSION IBACK(ND,NP)

      LMAX=MIN(NP+1-JP,ND)

      DO 1 L=1,LMAX

        DO 2 K=1,NVOPA
           XJK(L,K)=XJK(L,K)+WLK(L,JP)*((DINT(K,L)
     $              +DINT(K,IBACK(L,JP)))/2.0d0)
    2    ENDDO
    1 ENDDO

      RETURN

      END subroutine

      SUBROUTINE COMPGAU9(XJ,XJK,CWK,XNU,T,DSTEP,ND,NVOPA)

C***  Computation of the Gauss-Profile to fold on to the Thomson Emissivity
C***  The formula for the thermal Gauss Profile is written in
C***  Mihalas : Stellar Atmospheres (Second Ed.) page 420
C***  and has been translated so that the frequencies (XNU) are in
C***  Doppler Units

      implicit real*8(a-h,o-z)

      parameter (PI=3.1415926535898d0)
C***  Constant for thermal Gauss-Profile (= m(e)/(4k)) (cgs?)
      PARAMETER (GAUKONST=1.649538d-12)

      real*4 pot
      DIMENSION XJ(:,:),XJK(:,:)
      DIMENSION CWK(:,:),XNU(:),T(ND)

C***  First Reset the Intensity Array
      XJ(:,:) = 0.0
      
C***  now convolve the mean intensity
C***  with the Gauss profile of thermal electrons
C***  output: array XJ

      DO L =1,ND





        weight = dstep*sqrt(gaukonst/t(L)/pi)
        DO K=1,NVOPA
          DO KK=1,NVOPA
            XNU2=(XNU(K)-XNU(KK))*(XNU(K)-XNU(KK))
            pot=GAUKONST*XNU2/T(L)
            XJ(L,K)=XJ(L,K)+XJK(L,KK)
     $           *exp(-pot)*WEIGHT
          ENDDO


c version with fraction of integral
          XJ(L,K)=XJ(L,K)/CWK(L,K)




        ENDDO
      ENDDO

      RETURN

      END subroutine

      FUNCTION airlambda(vaclambda)
!    translate vacuum to the airlambda lambda in A                             
                                                                                
      IMPLICIT NONE
      real*8 vaclambda, airlambda, sig, n

      sig=1.d4/(vaclambda)
      n=1.d0+6.4328d-5+(2.94981d-2)/(146.-sig**2.)+(2.554d-4)/
     $ (41-sig**2.)
      airlambda=vaclambda/n
      return

      END FUNCTION

      function getLocalTMin(T,ND)
      implicit none
      integer :: getLocalTMin
      integer,intent(in) ::ND
      real*8,intent(in) :: T(ND)
      integer :: NDPMIN,L
      IF(ND==1) THEN
        getLocalTMin=1
        RETURN
      ELSE
      NDPMIN = 0                    ! Sanity check 1
      DO L=2,ND-1                   ! Find local minimum in T
        if ((T(L) .LT. T(L-1)) .AND. (T(L) .LT. T(L+1))) THEN
          NDPMIN = L
        endif
      ENDDO

!      PRINT*, 'TEMP MIN', NDPMIN, T(NDPMIN)

      if (NDPMIN .eq. 0) then       ! Sanity check 1 - finish
        print *,'linop: something wrong! NDPMIN = ',NDPMIN,
     $   ' setting NDPMIN to 1'
        NDPMIN = 1              ! continue, dont abort
      endif
      if ((T(NDPMIN) .gt. 5000d0) .and. (ND .gt. 1)) then !Sanity check 2
        print *,'something wrong with Temperature minimum!',NDPMIN, T(NDPMIN)
        write(6,*)'something wrong with Temperature minimum!',
     $                        NDPMIN,T(NDPMIN)
        pause
      endif
      ENDIF

      getLocalTMin = NDPMIN

      end function getLocalTMin

      function getContIdx(idxFreq)
      use MOD_ERROR
      use UTILS, only: ASSERT
      use SYNTHP_CONT, only: NFCONT
      implicit none
      include '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      integer i,getContIdx,idxFreq
      call assert(LOC(FRXIDX)/=0,'FRXIDX is not yet available')
      call assert(idxFreq>=1,'Frequency Index can not be negative')
      do i=1,NFCONT()-1
        if(idxFreq <= FRXIDX(i)) then
          getContIdx=i
          return
        endif
      enddo
      if(idxFreq<=FRXIDX(NFCONT())) then
        getContIdx=i
        return
      endif
      print *,'inibl:getContIdx:ERROR:',
     &                    idxFreq,NFCONT(),FRXIDX(NFCONT())
      call error('Could not find FrxIdx')
      end function getContIdx

      end module
