      module MOD_SDOT
      contains
      function sdot (N,A,N1,B,N2)
C***  INTEGRATION SUM, USING CRAY VECTOR FUNCTION SDOT (SCALAR PRODUCT)
C***  sum a*b
      implicit real*8(a-h,o-z)

      parameter (zero=0.0d0)
      
      DIMENSION A(N),B(N)

	if (n1.ne.n2) then
	   write (6,*) ' SDOT: not coded option N1.ne.N2'
	   stop 'SDOT'
	endif
	if (n1.ne.1) then
	   write (6,*) ' SDOT: not coded option N1.ne.1'
	   stop 'SDOT'
	endif

      sum=zero
      DO 1 I=N1,N
	   sum=sum+a(i)*b(i)
    1 enddo

      sdot=sum

      RETURN
      END function
      end module