      module MOD_FORMATS
      implicit none
      character(*),parameter :: FMT_KEY='(A4,I4)'  ! Format for the RADIO* files, etc. E.g. WCHA 111
      contains
      !***************************************
      !* Compare two strings if they are equal, ignoring all blanks
      !* Needed for the keys
      !*       (08/28/2007 - micha  initial version)
      function equalWithoutBlanks(A,B)
      implicit none
      character*(*),intent(in) :: A,B
      logical equalWithoutBlanks
      integer n,m,length 
      n=1;m=1;
      length=min(len(A),len(B))
      equalWithoutBlanks=.false.
      do while(max(n,m)<=length)
        do while(n<length.and.A(n:n)==' '); n=n+1; enddo  ! jump over blanks
        do while(m<length.and.B(m:m)==' '); m=m+1; enddo
        if(A(n:n) .ne. B(m:m)) return                     ! compare the non-blanks
        n=n+1;   m=m+1;
      enddo
      if(A(n:len(A)) .ne. '') return  ! Remainder must be all blanks
      if(B(m:len(B)) .ne. '') return
      equalWithoutBlanks=.true.  ! if we get here all is fine & equal
      end function
      end module