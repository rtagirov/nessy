      module MOD_ISRCHFLE
      contains
      function ISRCHFLE (N,array,istep,value)

	real*8 array(n), value, arr

	i=0
	do k=n,1,-istep
	   arr=array(k)
	   if (arr.le.value) i=k
	enddo

	if (i.gt.0) then
	   isrchfle=i
	   print *,' isrchfle= ',i,array(i),value
	else
	   write (6,*) 'isrchflt failed'
	   stop 'error'
	endif

	return
	end function
      end module