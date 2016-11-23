      module MOD_CWGAUSS9
      contains
      SUBROUTINE CWGAUSS9(CWK,XNU,ND,NVOPA,T,xjc,dstep)
C***************************************************************
C***  Compute the Integration Weights for the Gauss Profile 
c***  then compare with the theoretical weight and assume the
c***  rest has a intensity equal to the reference continuum
c***  The numbers CWK are used to correct the integration in
c***  routine COMPGAU
C***************************************************************
      implicit real*8(a-h,o-z)
      parameter (PI=3.141592654)
C***  Constant for thermal Gauss-Profile (= m(e)/(4k)) (cgs?)
      PARAMETER (GAUKONST=1.649538E-12)

      DIMENSION CWK(:,:),XNU(:),T(ND),xjc(nd)
      real*4 pot
      
      DO L=1,ND
       weight = dstep*sqrt(gaukonst/t(L)/pi)
       DO K=1,NVOPA
         CWK(L,K)=0.0
         DO KK=1,NVOPA
           XNU2=(XNU(K)-XNU(KK))*(XNU(K)-XNU(KK))
           pot=GAUKONST*XNU2/T(L)
           CWK(L,K)=CWK(L,K)+exp(-pot)
         ENDDO
         frac = weight*CWK(L,K)
         if (frac.gt.1.0000001) then
            print *,' sum .gt. 1:',L,K,frac
         endif
         if (frac.gt.1.) frac=1.

c         CWK(L,K)=(1.-frac)*xjc(L)
c version with fraction of integral
         CWK(L,K)=frac
       ENDDO
      ENDDO
c$$$c-pr
c$$$      L=40
c$$$      print *,' L, T(L) = ',L,T(L)
c$$$      do k=1,nvopa,20
c$$$       diff=xnu(k)*xnu(k)*GAUKONST/T(L)
c$$$       print *,k,xnu(k)/1.e5,diff,cwk(l,k)
c$$$      enddo
      RETURN
      END subroutine
      end module
