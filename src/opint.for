!***  changed by Margit Haberreiter
!*** changed by Micha Schoell
!*** used to be commons declaration, changed to ALLOCATABLES
      module OPINT
      implicit none



      integer,private :: NVDIM__ = -1
      integer,private :: NFOBS__ = -1
      integer,PARAMETER :: NMAX =1200,NDDD = 1000
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

        allocate(OPATOT(nvdim_,NDDD))

        allocate(ETATOT(nvdim_,NDDD))

        allocate(VOPA(nvdim_))

        allocate(TAU(NFOBS_,NMAX))

        allocate(DTAU(NFOBS_,NMAX))

        allocate(WTAU(NFOBS_,NMAX))

        allocate(SFINE(NFOBS_,NMAX))

        allocate(OPAFINE(NFOBS_,NMAX))

        allocate(XFINE(NMAX))

        allocate(ZFINE(NMAX))

        allocate(SCRAOP(nvdim__))

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
