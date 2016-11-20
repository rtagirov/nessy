!***  changed by Margit Haberreiter
!*** changed by Micha Schoell
!*** used to be commons declaration, changed to ALLOCATABLES
      module OPINT
      implicit none
      !integer,PARAMETER ::  NFMAX=1003
      !NVDIM  = 100 000,
      !integer,parameter :: NVDIM = 100 000
      integer,private :: NVDIM__ = -1
      integer,private :: NFOBS__ = -1
      integer,PARAMETER :: NMAX =1200,NDDD = 200
      real*8,allocatable :: OPATOT(:,:),ETATOT(:,:),VOPA(:)
     &            ,TAU(:,:),DTAU(:,:),WTAU(:,:)
     &            ,SFINE(:,:),OPAFINE(:,:)
     &            ,XFINE(:),ZFINE(:),SCRAOP(:)
      contains
      function NVDIM()
        use UTILS,only:assert
        implicit none
        integer NVDIM
        call assert(NVDIM__ /= -1, 'NVDIM not yet set!')
        NVDIM = NVDIM__
      end function

!       subroutine setNVDIM(VALUE)
!         use MOD_ERROR
!         use UTILS,only:assert
!         implicit none
!         integer,intent(in) :: VALUE
!         call assert(NVDIM__ .eq. -1, 'NVDIM already set')
!         NVDIM__ = VALUE
!       end subroutine
      subroutine OPINT_INIT(NVDIM_,NFOBS_)
        use UTILS, only:assert
        implicit none
        integer :: NVDIM_,NFOBS_
        PRINT*, 'AHTUNG: NVDIM_', NVDIM_
        call assert(NVDIM_==-1,'NVDIM already set')
        call assert(NFOBS_==-1,'NFMAX already set')
        nvdim__=NVDIM_
        nfobs__=NFOBS_
        call assert(nvdim_>0, 'nvdim_ to small')
!        PRINT*, 'ARRAY 1'
!        PRINT*, 'RIGHT BEFORE THE ALLOCATION:', NVDIM_, NDDD
!        PRINT*, 'OPATOT RIGHT BEFORE THE ALLOCATION:', NVDIM_, NDDD
        allocate(OPATOT(nvdim_,NDDD))
!        PRINT*, 'ARRAY 2'
        allocate(ETATOT(nvdim_,NDDD))
!        PRINT*, 'ARRAY 3'
        allocate(VOPA(nvdim_))
!        PRINT*, 'ARRAY 4'
        allocate(TAU(NFOBS_,NMAX))
!        PRINT*, 'ARRAY 5'
        allocate(DTAU(NFOBS_,NMAX))
!        PRINT*, 'ARRAY 6'
        allocate(WTAU(NFOBS_,NMAX))
!        PRINT*, 'ARRAY 7'
        allocate(SFINE(NFOBS_,NMAX))
!        PRINT*, 'ARRAY 8'
        allocate(OPAFINE(NFOBS_,NMAX))
!        PRINT*, 'ARRAY 9'
        allocate(XFINE(NMAX))
!        PRINT*, 'ARRAY 10'
        allocate(ZFINE(NMAX))
!        PRINT*, 'ARRAY 11'
        allocate(SCRAOP(nvdim__))
!        PRINT*, 'ARRAY 12'
      end subroutine

      subroutine OPINT_EXIT
        deallocate(OPATOT)
        deallocate(ETATOT)
        deallocate(VOPA)
        deallocate(TAU)
        deallocate(DTAU)
        deallocate(WTAU)
        deallocate(SFINE)
        deallocate(OPAFINE)
        deallocate(XFINE)
        deallocate(ZFINE)
        deallocate(SCRAOP)
      end subroutine
      end module
