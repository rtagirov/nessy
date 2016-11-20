      module MOD_COMPGAU9
      contains
      SUBROUTINE COMPGAU9(XJ,XJK,CWK,XNU,T,DSTEP,NDDIM,ND,NVOPA)
C****************************************************************************
C***  Computation of the Gauss-Profile to fold on to the Thomson Emissivity.*
C***  The formula for the thermal Gauss Profile is written in               *
C***  Mihalas : Stellar Atmospheres (Second Ed.) page 420                   *
C***  and has been translated so that the frequencies (XNU) are in          *
C***  Doppler Units                                                         *
C****************************************************************************
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
c$$$c-pr
c$$$        if (L.eq.40) then
c$$$           print *,' L= ',L
c$$$           print *,' XJK(40,K)',(XJK(40,KK),KK=1,NVOPA)
c$$$        endif
        weight = dstep*sqrt(gaukonst/t(L)/pi)
        DO K=1,NVOPA
          DO KK=1,NVOPA
            XNU2=(XNU(K)-XNU(KK))*(XNU(K)-XNU(KK))
            pot=GAUKONST*XNU2/T(L)
            XJ(L,K)=XJ(L,K)+XJK(L,KK)
     $           *exp(-pot)*WEIGHT
          ENDDO
c version with continuum contribution estimated from elimin
c          XJ(L,K)=XJ(L,K)+CWK(L,K)
c version with fraction of integral
          XJ(L,K)=XJ(L,K)/CWK(L,K)
c$$$c-pr
c$$$          if (K.eq.1 .or. K.eq.NVOPA/2 .or. K.eq.NVOPA) then
c$$$             print *,L,K,XJ(L,K),CWK(L,K)
c$$$          endif
        ENDDO
      ENDDO
      RETURN
      END subroutine
      end module