!***********************************************************************
!***  DECODING INPUT CARDS FOR PROGRAM "FORMAL" FROM FTO5 = $IN (DEFAULT)
!***  THIS ROUTINE IS LEFT WHEN A "LINE"-OPTION IS FOUND,
!***  AND ENTERED AGAIN WHEN THE LINE HAS BEEN COMPUTED.
!***********************************************************************
      module MOD_DECF_SYN
      type :: CARD_PARAMS_
        integer:: ABEMLIN_READ
        integer:: ABEMLIN_WRITE
        integer:: ABEMLIN_AUTO
      end type
      type :: CARD_OPTS
        integer :: ABEMLIN = 0  ! 0 for nothing, 1 for read, 2 for write
        character*80 :: ABEMLIN_PATH = ''
        logical :: PRINT_TAU = .false.
        character*4 :: COMPRESS = ''
        integer     :: CMPRS_TMIN_OFFSET
        integer     :: CMPRS_LENGTH
      end type
      type(CARD_OPTS) :: CARDS
      type(CARD_PARAMS_),parameter :: CARD_PARAMS=CARD_PARAMS_(1,2,3)
      logical :: odf_cards = .FALSE.
      character*50 :: odf_name
C     ************

      real*8,allocatable :: wav_f(:), ffactor(:)
      integer :: Nfudge
C     ***********


      contains
      SUBROUTINE DECF_SYN (KARTE,PLOT,NFOBS,LSOPA,FMAX,FMIN,
     $        MAXITER,felsca,RWLAE,PHEAD,PROLIB,VSINI,SHIFT,
     $        LSPRO,LSDWL,NORM,TRANS,FIN,VDOP,
     $        NPHIP,LPSTI,LPENI,JFIRSI,JLASI,COROT,iTionsel,XLMIN,XLMAX)
      use MOD_ERROR
      use CONSTANTS,only:CLIGHT_SI
      implicit real*8(a-h,o-z)
      !integer,intent(in   ) :: NFODIM
      integer,intent(inout) :: NFOBS,LSOPA,MAXITER,LSDWL,JFIRSI,JLASI
      integer,intent(inout) :: iTionsel,LSPRO,NPHIP,LPSTI,LPENI
      real*8, intent(inout) :: RWLAE,FMAX,FMIN,VDOP,VSINI,SHIFT,felsca
      logical,intent(inout) :: NORM,PLOT,TRANS,FIN,PROLIB,COROT
      character,intent(inout):: KARTE*80,PHEAD*28
      real*8,parameter :: clight = CLIGHT_SI/1e3 !2.997924562d+5
      integer :: l

      REAL*8, INTENT(OUT) :: XLMIN, XLMAX

C    ****************

       Nfudge=161

       if (.not. allocated(wav_f) ) then

       allocate(wav_f(Nfudge))
       allocate(ffactor(Nfudge))

       endif

        OPEN (UNIT=240,FILE='../fudge_As.txt',STATUS='OLD',READONLY)

        do l=1, Nfudge
        READ (240,*) wav_f(l), ffactor(l)
!        ffactor(l)=max(1.,ffactor(l))
!        ffactor(l)=1.
        enddo


        close (240)



C    ****************


      CARDS.ABEMLIN_PATH=''; CARDS.ABEMLIN=0

    1 READ (1,'(A)',END=66) KARTE
      PRINT '(A80)', KARTE
      IF ( KARTE(:4) .EQ. 'LINE' ) THEN
C                          ====
       IF ( KARTE(10:15) .EQ. 'PROLIB' ) THEN
         DECODE(80,'(15X,F10.3)',KARTE) RWLAE
         PROLIB =  (RWLAE.GT.0.)
       ENDIF
       RETURN  ! LINE is last line in cards
      ENDIF

      IF ( KARTE(:7).EQ. 'ABEMLIN') THEN  ! ABEMLIN READ .... or WRITE ....
                        !========
        IF (KARTE(9:12) .EQ. 'READ' ) THEN
          CARDS.ABEMLIN_PATH  = adjustl(karte(13:))
          CARDS.ABEMLIN       = CARD_PARAMS.ABEMLIN_READ
        ELSEIF (KARTE(9:13) .EQ. 'WRITE') THEN
          CARDS.ABEMLIN_PATH  = adjustl(karte(14:))
          CARDS.ABEMLIN       = CARD_PARAMS.ABEMLIN_WRITE
        ELSEIF (KARTE(9:12) .EQ. 'AUTO') THEN
          CARDS.ABEMLIN_PATH  = adjustl(karte(14:))
          CARDS.ABEMLIN       = CARD_PARAMS.ABEMLIN_AUTO
        ELSEIF (KARTE(9:10) .EQ. 'NO' ) THEN
          CARDS.ABEMLIN_PATH  = ''
          CARDS.ABEMLIN       = 0
        ELSE
          call ERROR('decf_syn: ABEMLIN: Must be followed by '//
     $               'either READ, WRITE, AUTO or NO (followed by the path)')
        ENDIF
        print '(A,i1)','decf_syn: ABEMLIN = ', CARDS.ABEMLIN
        GOTO 1
      ENDIF


      IF ( KARTE(:8) .EQ. 'INTERVAL' ) THEN
C                          ========
         read (karte(10:80),*) xlmin,xlmax
         if (xlmax.le.xlmin) then
            print *,' interval not defined'
            print *,' input was: ',xlmin,xlmax
            stop ' error DECF_SYN'
            endif
         print *,' Interval decoded: L-min=',xlmin,' L-max=',xlmax
         rwlae=(xlmin+xlmax)/2.
         FMAX=(xlmax-rwlae)/rwlae*clight
         PROLIB=.TRUE.
         ! PRINT *,'/PFIRST:PLAST=',PFIRST,PLAST
         return
         endif
      IF(KARTE(:8) == 'COMPRESS') THEN
        read(KARTE,'(9x,A4,X,i3,X,i3)') CARDS%COMPRESS
     &                                 ,CARDS%CMPRS_TMIN_OFFSET
     &                                 ,CARDS%CMPRS_LENGTH
        print '("DECFSYN: COMPRESS SET to :",A4,": offset = ",i0,"'//
     &  ', length= ",i6)'  ,CARDS%COMPRESS ,   CARDS%CMPRS_TMIN_OFFSET
     &                     ,CARDS%CMPRS_LENGTH
        GOTO 1
      ENDIF

      IF ( KARTE(:11) .EQ. 'SELECT TION' ) THEN
C                           ===========
            read (karte(12:22),*) asel
		  isel = nint(asel)
            if (isel.gt.0. .and. isel.lt.4) then
               if (isel.eq.1) then
                  print *,' H0/H+ ionization Temp selected'
                  iTionsel=3

               else if (isel.eq.2) then
                  print *,' He0/He+ ionization Temp selected'
                  iTionsel=1
               else if (isel.eq.3) then
                  print *,' He+/He++ ionization Temp selected'
                  iTionsel=2
               endif
            else
               if (isel.eq.4) then
                  print *,' Te ionization Temp selected'
                  iTionsel=-1
               endif
            endif
            GOTO 1
            ENDIF
      IF ( KARTE(:9) .EQ. 'PRINT TAU' ) THEN
!                          ====
            CARDS%PRINT_TAU=.TRUE.
            GOTO 1
      ENDIF
      IF ( KARTE(:7) .EQ. 'MAXITER' ) THEN
C                          =======
            read (karte(9:18),'(F10.8)') FITER
            if (fiter.gt.0.) then
               maxiter=nint(fiter)
            else
               maxiter=1
            endif
            GOTO 1
            ENDIF
      IF ( KARTE(:10) .EQ. 'ELSC_WIDTH' ) THEN
C                          ==========
            read (karte(12:21),'(F10.8)') felsca
            if (felsca.lt.0) felsca=0.
            GOTO 1
            ENDIF
      IF ( KARTE(:4) .EQ. 'PLOT' ) THEN
C                          ====
            PLOT=.TRUE.
            GOTO 1
            ENDIF
      IF ( KARTE(:4) .EQ. 'NORM' ) THEN
C                          ====
            NORM=.TRUE.
            print *,' Normalize option selected'
            GOTO 1
            ENDIF
      IF ( KARTE(:8) .EQ. 'TRANSFER' ) THEN
C                          ========
            TRANS=.TRUE.
            GOTO 1
            ENDIF
      IF (KARTE(:9) .EQ. 'TRANS DWL') THEN
C                         =========
            IF (LSDWL.LT.0) LSDWL=2
            IF (LSDWL.EQ.1) LSDWL=4
            GOTO 1
            ENDIF
      IF ( KARTE(:9) .EQ. 'PRINT DWL') THEN
C                          ========
            IF (LSDWL.LT.0) LSDWL=1
            IF (LSDWL.EQ.2) LSDWL=4
            GOTO 1
            ENDIF
      IF ( KARTE(:6) .EQ. 'PROLIB' ) THEN
C                          ======
         DECODE(80,3,KARTE) PHEAD
    3    FORMAT(10X,A28)
         GOTO 1
         ENDIF
      IF ( KARTE(:5) .EQ. 'SHIFT' ) THEN
C                          =====
            read (karte(7:16),'(F10.8)') shift
            if (shift.lt.0) shift=0.
c            print '(A,F10.0,A)',' Wavelength shift by ',shift,' [km/s]'
            GOTO 1
            ENDIF
      IF ( KARTE(:5) .EQ. 'VSINI' ) THEN
C                          =====
         DECODE(80,6,KARTE) VSINI
    6    FORMAT(6X,F10.0)
         GOTO 1
         ENDIF
      IF ( KARTE(:5) .EQ. 'COROT') THEN
C                          =====
         COROT=.TRUE.
         GOTO 1
         ENDIF
      IF ( KARTE(:4) .EQ. 'NPHI' ) THEN
C                          ====
         DECODE(80,26,KARTE) XNPHI
 26      FORMAT (5X,F10.0)
         NPHIP=NINT(XNPHI)
         GOTO 1
         ENDIF
      IF ( KARTE(:6) .EQ. 'PHISTA' ) THEN
C                          ======
         DECODE(80,5,KARTE) PHISTA
    5    FORMAT(7X,F10.0)
         LPSTI=NINT(PHISTA)
         GOTO 1
         ENDIF
      IF ( KARTE(:6) .EQ. 'PHIEND' ) THEN
C                          ======
         DECODE(80,5,KARTE) PHIEND
         LPENI=NINT(PHIEND)
         GOTO 1
         ENDIF
      IF ( KARTE(:6) .EQ. 'PFIRST' ) THEN
C                          ======
         DECODE(80,5,KARTE) PFIRST
         JFIRSI=NINT(PFIRST)
         GOTO 1
         ENDIF
      IF ( KARTE(:5) .EQ. 'PLAST' ) THEN
C                          =====
         DECODE(80,9,KARTE) PLAST
         JLASI=NINT(PLAST)
         GOTO 1
         ENDIF
      IF ( KARTE(:9) .EQ. 'PRINT PRO') THEN
C                          ========
            LSPRO=1
            GOTO 1
            ENDIF
      IF ( KARTE(:5) .EQ. 'NFOBS' ) THEN
C                          =====
            DECODE (80,9,KARTE) XNFOBS
    9       FORMAT (6X,F10.0)
            NFOBS=NINT(XNFOBS)
c???            IF (NFOBS .GT. NFODIM/3) NFOBS=NFODIM/3
            !IF (NFOBS.GT.NFODIM) NFOBS=NFODIM
            GOTO 1
            ENDIF
      IF ( KARTE(:10) .EQ. 'PRINT OPAL' ) THEN
C                           ==========
            DECODE (80,12,KARTE) XL
   12       FORMAT (10X,F10.0)
            LSOPA=1
            IF (XL.GT.1) LSOPA=NINT(XL)
            GOTO 1
            ENDIF
      IF ( KARTE(:4) .EQ. 'FMAX' ) THEN
C                          ====
            DECODE (80,19,KARTE) FMAX
   19       FORMAT (5X,F10.0)
            GOTO 1
            ENDIF
      IF ( KARTE(:4) .EQ. 'FMIN' ) THEN
C                          ====
            DECODE (80,19,KARTE )FMIN
            FMIN=ABS(FMIN)
            GOTO 1
            ENDIF
      IF ( KARTE(:4) .EQ. 'VDOP' ) THEN
C                          ====
            DECODE (80,22,KARTE) VDOP
   22       FORMAT (5X,F10.0)
            GOTO 1
            ENDIF
      IF ( KARTE(:3) .EQ. 'ODF') THEN
C                          ===
         odf_cards = .TRUE.
         GOTO 1
         ENDIF
      IF (KARTE(:8) .EQ. 'odf_name') THEN
C                         ========
!          READ (1,'(A)',END=66) odf_name
         odf_name = KARTE(10:)
         print*, "odf_namedecfy: ", odf_name
         GOTO 1
         ENDIF


      GOTO 1

66    CONTINUE
C***  EOF REACHED
99    FIN=.TRUE.
      RETURN
c	 pause
      END subroutine
      end module
