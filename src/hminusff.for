      module MOD_HMINUSFF
      contains
      subroutine hminusff (sighmff,xlam,T)
C***  source John 1988 A&A 193,189
c     control table: Stilley & Callaway 1970, ApJ 160, 245
!      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT NONE
      real*8,intent(out) :: sighmff
      real*8, intent(in) :: xlam,T
      real*8 an,bn,cn,dn,en,fn,sum,sumn,theta,pow,sq,x
      integer i,j
      real*8 :: sqrt
      dimension an(6,2),bn(6,2),cn(6,2),dn(6,2),en(6,2),fn(6,2)
C23456789012345678901234567890123456789012345678901234567890123456789012
      data an/0.,     2483.346, -3449.889,  2200.04,  -696.271,  88.283,
     &      518.1021,  473.2636, -482.2089,  115.5291,   0.,      0./,
     &     bn/0.,      285.827, -1158.382,  2427.719,-1841.4,   444.517,
     &     -734.8666, 1443.4137, -737.1616,  169.6374,   0.,      0./,
     &     cn/0.,    -2054.291,  8746.523,-13651.105, 8624.97,-1863.864,
     &     1021.1775,-1977.3395, 1096.8827, -245.649,    0.,      0./,
     &     dn/0.,     2827.776,-11485.632,16755.524,-10051.53, 2095.288,
     &     -479.0721,  922.3575, -521.1341,  114.243,    0.,      0./,
     &     en/0.,    -1341.537,  5303.609, -7510.494, 4400.067,-901.788,
     &       93.1373, -178.9275,  101.7963,  -21.9972,   0.,      0./,
     &     fn/0.,      208.952,  -812.939,  1132.738, -655.02,  132.985,
     &       -6.4285,   12.36,     -7.057,     1.5097,   0.,      0./
C23456789012345678901234567890123456789012345678901234567890123456789012
	x=xlam/1.e+4
      if (x.gt.0.3645d0) then
         i=1
      else if (x.ge.0.1823d0) then
         i=2
      else
         sighmff=0.d0
         return
      endif
      theta=5040.d0/T
      sq = sqrt(theta)
      pow=sq
      sum=0.d0
      do j=1,6
         pow=pow*sq
         sumn=(((((x*an(j,i)+bn(j,i))*x+cn(j,i))*x+dn(j,i))*x+en(j,i))
     &        *x+fn(j,i))/x/x/x
         sum=sum+sumn*pow
      enddo
      sighmff=1.d-29*sum
      return
      end subroutine
      end module
