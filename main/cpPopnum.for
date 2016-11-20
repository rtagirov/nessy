      module MOD_CPPOPNUM
      contains
      !*** copy the POPNUM File to the file POPNUM.(JOBNUM)
      !*** only needed for analysis
      !*** option set in CARDS: POPNUM CP <jobnum-steps>
      subroutine cpPopnum(JOBNUM)
      use UTILS
      implicit none
      integer,intent(in) :: JOBNUM
      character*(40):: flname
      character*(220) :: X
      integer :: ifl,ofl
      ifl=getFileUnit(start=10)
      open(unit=ifl,status='old',action='read',name='POPNUM')
      flname = 'POPNUM.'//adjustl(int2str(JOBNUM))
      ofl=getFileUnit(start=10)
      open(unit=ofl,action='write',name=flname)
      do
        read(ifl,'(A)',end=900) X
        write(ofl,'(A)') trim(X)
      enddo
 900  continue 
      close(ifl)
      close(ofl)
      end subroutine
      end module