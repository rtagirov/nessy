      module MOD_WRITMS
      contains
      SUBROUTINE WRITMS (ifl,array,n,cname,iflag,ierr)

      implicit real*8(a-h,o-z)
      integer,intent(in)::ifl,n,iflag
      integer,intent(inout)::ierr
      real*8,intent(in):: array(n)
      character*(*),intent(in):: cname
      character*10 char
      char=cname

      ierr=1
      if (iflag.eq.-1) then
        write (ifl,'(A10)') char
        write (ifl,*) (array(i),i=1,N)
        ierr=0
      else
        write(6,*) ' MODFILE write mode unknown'
        pause
        stop
      endif

	end subroutine

      SUBROUTINE WRITMS1 (ifl,scalar,cname,iflag,ierr)
      implicit none
      integer,intent(in)::ifl,iflag
      integer,intent(inout)::ierr
      real*8,intent(in):: scalar
      character*(*),intent(in):: cname
      integer idummy
      character*10 char
      char=cname

      ierr=1
      if (iflag.eq.-1) then
         write (ifl,'(A10)') char
         write (ifl,*) scalar
         ierr=0
      else
         write(6,*) ' MODFILE write mode unknown'
         pause
         stop
      endif

      end subroutine
      end module