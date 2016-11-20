      module MOD_COMPRO
      contains
      SUBROUTINE COMPRO (PROFILE,EMINT,NFOBR,W,JFIRST,JLAST)
C*******************************************************************
C***  Update the Profile for each Impact Parameter                 *
C*******************************************************************
      implicit real*8(a-h,o-z)
      DIMENSION PROFILE(NFOBR)
      DIMENSION EMINT(NFOBR)
 
    !  print*, 'w test', W
             
      DO K=1,NFOBR
          IF (JFIRST .EQ. JLAST) THEN
              PROFILE(K)=EMINT(K)
          ELSE
              PROFILE(K)=PROFILE(K)+W*EMINT(K)
          ENDIF
      ENDDO
     
  !    print*, 'test sasha', profile(1), profile(2), profile(3)

      RETURN
      END subroutine
      end module
