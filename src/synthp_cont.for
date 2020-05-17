      module SYNTHP_CONT

      implicit none

      integer, private :: NFCONT__ = -1

      real*8, allocatable :: FREQC(:), ABSOC(:), EMISC(:)!, SCATC(:)

      real*8, allocatable :: absoc_rayleigh(:)

      contains

      subroutine setNFCONT(NFCONT_VAL)
        use UTILS
        integer,intent(in) :: NFCONT_VAL
        call assert(NFCONT__==-1,'NFCONT already set')
        call assert(NFCONT_VAL>0,'NFCONT must be > 0!')
        NFCONT__ = NFCONT_VAL
      end subroutine
      
      function NFCONT()

!        use utils

        integer :: NFCONT

!        if (NFCONT == -1) stop 'NFCONT not set'

        NFCONT = NFCONT__

      end function

      subroutine SYNTHP_CONT_INIT(NFCONT_VAL)
        integer,intent(in) :: NFCONT_VAL
        call setNFCONT(NFCONT_VAL)
        print *,'SYNTHP_CONT_INIT: Allocate with NFCONT = ',NFCONT()
        allocate(FREQC(NFCONT()))
        allocate(ABSOC(NFCONT()))
        allocate(EMISC(NFCONT()))

        allocate(absoc_rayleigh(nfcont()))

!        allocate(SCATC(NFCONT()))
      endsubroutine
      endmodule
