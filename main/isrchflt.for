      module MOD_ISRCHFLT
      contains
      function ISRCHFLT (N,array,istep,value)

	real*8 array(n), value, arr

	i=0
	do k=n,1,-istep
	   arr=array(k)
	   if (arr.lt.value) i=k
	enddo

	if (i.gt.0) then
	   isrchflt=i
	   print *,' isrchflt= ',i,array(i),value
	else
	   write (6,*) 'isrchflt failed'
	   stop 'isrchflt: error'
	endif

	return
	end function
      end module