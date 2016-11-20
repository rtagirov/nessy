      module MOD_COMPXJ9
      contains
      SUBROUTINE COMPXJ9 (ND,NP,JP,NVOPA,DINT,XJK,WLK,IBACK,NDDIM)
C***********************************************************************
C***  Compute mean Intensity (Moment 0)                                *
C***  add for each impact parameter the intensities DINT               *
C***  WLK are the weigths given by the geometry of the grid            *
C***  This routine is called within the impact parameter loop          *
c     modified: 20.1.01  do loop only up tp LMAX
C***********************************************************************
      implicit real*8(a-h,o-z)

      DIMENSION DINT(:,:),XJK(:,:),WLK(ND,NP)
      DIMENSION IBACK(ND,NP)

      LMAX=MIN(NP+1-JP,ND)
      DO 1 L=1,LMAX
c      DO L=1,ND
        DO 2 K=1,NVOPA
           XJK(L,K)=XJK(L,K)+WLK(L,JP)*((DINT(K,L)
     $              +DINT(K,IBACK(L,JP)))/2.0d0)
    2    ENDDO
    1 ENDDO
      RETURN
      END subroutine
      end module