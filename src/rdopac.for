      module MOD_RDOPAC
      type ROW
        integer :: LBKG_IDX
        real*8,allocatable :: V(:)
      endtype
      type(ROW),private,allocatable:: ROWS(:)
      contains
      subroutine RDOPAC(XLAM,LINOP,NDPMIN,XLBKG1,XLBKG2)
      !*********************************************
      !***  PROGRAM READS THE BINNED LINE OPACITIES
      !***  AND INCLUDES THEM IN THE HMINUS CODE TO
      !***  CALCULATE THE CONTINUUMS OPACITIES
      !*** L = 1 : MOST OUTWARD DEPTH POINT
      !***  L = ND: MOST INWARD DEPTH POINT
      !*********************************************
      use MOD_ERROR
      implicit none
      real*8, intent(out) :: LINOP(:)
      real*8, intent(in)  :: XLAM
      integer,intent(in)  :: XLBKG1,XLBKG2,NDPMIN
      integer             :: LAM
      character*50        :: FLNAM !,DUMMY
      integer :: i,IDX, nLINOP
      if(.not.allocated(ROWS)) then
        allocate(ROWS(0))
      endif
      !**** ABBIN: BINNED LINE OPACITY -  OPACITY DISTRIBUTION FUNCTION
      !**** XLAM: WAVELENGTH GIVEN IN FGRID INCLUDING EDGE WAVELENGHTS
      LINOP = 0.
      IDX=-1
      IF (XLAM < XLBKG1) RETURN
      IF (XLAM > XLBKG2) RETURN
      nLINOP = size(LINOP)
      LAM = XLAM
      LAM = ((XLAM)/10);   LAM = LAM*10+5

      !*** Search for frequency index
      SEARCH_IDX: do i=1,size(ROWS)
        if(ROWS(i)%LBKG_IDX == LAM) then
          IDX=i
          if(size(ROWS(i)%V) <= nLINOP) exit SEARCH_IDX !** reread the LINOP
          LINOP = ROWS(i)%V(1:nLINOP)
          return !** if found and correct size, return
        endif
      enddo SEARCH_IDX
      write(flnam,'(i6)') LAM
      flnam=adjustl(adjustr(flnam)//'.lbkg')
      open (300,file=flnam,status='old',action='read', err=999)
      read (300,'(A)')           ! read the header
      !***********************************************************
      !*  L = 1 : MOST OUTWARD DEPTH POINT
      !*  L = ND: MOST INWARD DEPTH POINT
      !*  If depth point is outward of temperature minimum,
      !*  use LINOP of temperature minimum
      read(300,'(1pe12.5)') LINOP
      LINOP(1:NDPMIN) = LINOP(NDPMIN)
      close (300)
      !************************************************************
      if(IDX>0) then 
        !*** Change the size of ROWS(IDX)%V and replace with new value
        deallocate(ROWS(IDX)%V)
        allocate(ROWS(IDX)%V(nLINOP))
        ROWS(IDX)%V = LINOP
      else
        call append_ROW(LAM,LINOP)  ! ROW(end+1)%(LBKG_IDX, V) = (LAM, LINOP)
      endif

      if  (any(linop  <  0.)) then
        print *,xlam,'linop= ',linop
        call ERROR('rdopac: negative linop in ' // flnam)
      endif
      RETURN

      !***************************************************************
      !* Error handling
  999 print *, 'RDOPAC: IO Error opening file ' // flnam
      print '("XLBKG1,2 = ", i7," ",i7," xlam =",e10.4," lam=",i6)',
     $    XLBKG1,XLBKG2,xlam,lam
      call ERROR('RDOPAC: IO Error opening file ' // flnam)

      END subroutine


      !****************************************
      !* Append a row to the array of rows.
      subroutine append_ROW(LAM,LINOP)
      integer,intent(in) :: LAM
      real*8,intent(in) :: LINOP(:)
      type(ROW),allocatable :: old_ROWS(:)
      integer :: nrow
      nrow = size(ROWS)
      !* make a copy of the original, reallocate by n+1 and copy back the old one to the first n rows.
      allocate(old_ROWS(nrow))
      old_ROWS = ROWS
      deallocate(ROWS)
      allocate(ROWS(nrow+1))
      ROWS(1:nrow) = old_ROWS

      !* Set Data for the last row (LBKG_IDX, V)
      ROWS(nrow+1)%LBKG_IDX = LAM
      allocate(ROWS(nrow+1)%V(size(LINOP)))
      ROWS(nrow+1)%V = LINOP
      end subroutine append_ROW
      end module
