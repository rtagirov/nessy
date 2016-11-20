      module MOD_WRITMSI
      contains
      SUBROUTINE WRITMSI (ifl,irray,n,cname,iflag,ierr)
      !*** Write the array irray(n) to the file ifl under the entry cname 
      !*** INPUT:
      !*** ifl    : File Descriptor to write to
      !*** irray  : Integers to write
      !*** n      : Number of Entries in irray
      !*** cname  : Name of entry
      !*** iflag  : must be set to -1
      !*** IN/OUTPUT:
      !*** ierr   : stores the error if one occured, otherwise 0
      !*** 
      implicit none
      integer,intent(in)::ifl,n,iflag
      integer,intent(inout)::ierr
      integer,intent(in):: irray(n)
      character*(*),intent(in):: cname
      character*10 char
      integer i
      char=cname
      ierr=1
      if (iflag.eq.-1) then
        write (ifl,'(A10)') char
        write (ifl,*) (irray(i),i=1,N)
        ierr=0
      else
        write(6,*) ' MODFILE write mode unknown'
        pause
        stop
      endif

      return
      end subroutine

      SUBROUTINE WRITMSI1 (ifl,iscalar,cname,iflag,ierr)
      !*** Write the scalar iscalar to the file ifl under the entry cname 
      !*** INPUT:
      !*** ifl:     File Descriptor to write to
      !*** iscalar: Integer to write
      !*** cname  : Name of entry
      !*** iflag must be set to -1
      !*** IN/OUTPUT:
      !*** ierr returns an error if one occured, otherwise 0
      !*** 
      implicit none
      integer,intent(in)   :: ifl,iflag
      integer,intent(inout):: ierr
      integer,intent(in)   :: iscalar
      character*(*),intent(in):: cname
      character*10 char
      char=cname
      ierr=1
      if (iflag.eq.-1) then
         write (ifl,'(A10)') char
         write (ifl,*) iscalar
         ierr=0
      else
         write(6,*) ' MODFILE write mode unknown'
         pause
         stop
      endif
      return
      end subroutine
      end module