      module MOD_ISAMAX
      contains
      function isamax (n,array,istep)

	real*8 array(n), aval, amax

	i=0
	amax=0.
	do l=1,n,istep
	   aval=abs(array(l))
	   if (aval.gt.amax) then
	      amax=aval
            i=l
	   endif
	enddo
	
	if (i.gt.0) then
         isamax=i
	else
	   write (6,*) ' isamax failed'
	   stop 'error'
	endif

	return
	end function
      end module
