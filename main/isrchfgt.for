      module MOD_ISRCHFGT
      contains
      function ISRCHFGT (N,array,istep,value)
      use MOD_ERROR

      real*8 array(n), value, arr

      if (value.ge.array(n)) then
c	   print *,' ISRCHFGT: arr(n)=',array(n),' val=',value
         isrchfgt=n+1
         return
      endif

      i=0
      do k=n,1,-istep
         arr=array(k)
         if (arr.gt.value) i=k
      enddo

      if (i.gt.0) then
         isrchfgt=i
c        print *,' isrchfgt= ',i
      else
         write (6,*) 'isrchfgt failed'
         pause
         call ERROR( 'isrchfgt failed' )
      endif

      return
      end function
      end module