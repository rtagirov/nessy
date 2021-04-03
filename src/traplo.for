      module MOD_TRAPLO
      contains
      SUBROUTINE traplo(PROFILE,DLAM,NFOBS,KARTE,MODHEAD,JOBNUM,RWLAE,PHEAD,PROLIB)

      implicit real*8(a-h,o-z)
      real*8, intent(in) :: PROFILE(NFOBS),DLAM(NFOBS),RWLAE
      integer,intent(in) :: NFOBS
      LOGICAL,intent(in) :: PROLIB
      CHARACTER,intent(in):: MODHEAD*104,PHEAD*28,KARTE*80
C***  DIRECT TRANSFER OF PLOT DATA
      real*8,allocatable:: PROF(:),DLAMBD(:),DPROF(:)
      !DIMENSION DLAM(NFOBS),PROFILE(NFOBS)
      CHARACTER HEADER*64,idat*8,itime*8,flnam*18,henam*18

      ! RT: Introduced for calculation of brightness temperatures in radio range
      !*************************************************************************
      REAL*8 :: Wavelength, I_mean
      !*************************************************************************

      IDAT=MODHEAD(15:22)
      ITIME=MODHEAD(25:32)
      NF=NFOBS
c*** change here if more than 63 points should be printed
c     IF (NOUT.gt.63) NOUT=63
c-nout      IF (NOUT.gt.300) NOUT=300
      NOUT=NFOBS

      !if (nfobs.gt.999) stop ' TRAPLO - NFOBS too large'
      allocate(PROF(NFOBS))
      allocate(DLAMBD(NFOBS))
      allocate(DPROF(NFOBS))
      DO 7 I=1,NFOBS
         PROF(I)=PROFILE(NFOBS+1-I)
         DLAMBD(I)=DLAM(NFOBS-I+1)+RWLAE
 7    CONTINUE

c*** reduce the number of profile points
      IF (NOUT.lt.nf) then
 11      CONTINUE
         DO 12 K=2,NF-1
            DLAM1=(DLAMBD(K+1)-DLAMBD(K))/(DLAMBD(K+1)-DLAMBD(K-1))
            DLAM2=1.-DLAM1
            DPROF(K)=ABS(PROF(K-1)*DLAM1+PROF(K+1)*DLAM2-PROF(K))
 12      CONTINUE
         DMAX=1.
         KMIN=NF+1
         DO 13 K=2,NF-1
            IF (DPROF(K).LT.DMAX) then
               KMIN=K
               DMAX=DPROF(K)
            endif
 13      CONTINUE
         IF (KMIN.GE.NF) THEN
            PRINT *,' SOMETHING WRONG,  PROF=',PROF
            STOP
         ENDIF
         DO 14 K=KMIN+1,NF
            PROF(K-1)=PROF(K)
            DLAMBD(K-1)=DLAMBD(K)
 14      CONTINUE
         NF=NF-1
         IF (NF.GT.NOUT) GOTO 11
      endif
  
 
c***  plot file
   
      IF (RWLAE.GE.10000.d0) THEN
         write (header,5)  NOUT,RWLAE,PHEAD,IDAT
c         ENCODE (60,5,HEADER) NOUT,RWLAE,PHEAD,IDAT
    5    FORMAT (I5,1x,F7.0,1X,A28,2X,A8)
      ELSE
         write (header,55) NOUT,RWLAE,PHEAD,IDAT
c     ENCODE (60,55,HEADER) NOUT,RWLAE,PHEAD,IDAT
 55      FORMAT (I5,1x,F6.1,1X,A28,2X,A8)
      ENDIF

   
   

 
      if (rwlae.lt.100.) then
         write (flnam,'(F4.1,6H.mdisp)') rwlae
         write (henam,'(F4.1,6H.title)') rwlae
      else if (rwlae.lt.1000.) then
         write (flnam,'(F4.0,5Hmdisp)') rwlae
         write (henam,'(F4.0,5Htitle)') rwlae
      else if (rwlae.lt.10000.) then
         write (flnam,'(F5.0,5Hmdisp)') rwlae
         write (henam,'(F5.0,5Htitle)') rwlae
      else if (rwlae.lt.100000.) then
         write (flnam,'(F6.0,5Hmdisp)') rwlae
         write (henam,'(F6.0,5Htitle)') rwlae
      else if (rwlae.lt.1000000.) then
         write (flnam,'(F7.0,5Hmdisp)') rwlae
         write (henam,'(F7.0,5Htitle)') rwlae
      else if (rwlae.lt.1.d7) then
         write (flnam,'(F8.0,5Hmdisp)') rwlae
         write (henam,'(F8.0,5Htitle)') rwlae
      else if (rwlae.lt.1.d8) then
         write (flnam,'(F9.0,5Hmdisp)') rwlae
         write (henam,'(F9.0,5Htitle)') rwlae
      else if (rwlae.lt.1.d9) then
         write (flnam,'(F10.0,5Hmdisp)') rwlae
         write (henam,'(F10.0,5Htitle)') rwlae
      else if (rwlae.lt.1.d10) then
         write (flnam,'(F11.0,5Hmdisp)') rwlae
         write (henam,'(F11.0,5Htitle)') rwlae
      else
         print *,' rwlae=',rwlae
         stop 'error'
      endif

      open (10,file=henam,status='unknown')
      WRITE (10,2) header
 2    format ('TITLE= ',A)
      close (10)

      open (10,file=flnam,status='unknown')
       
      if (DLAMBD(1) .lt. 67660.0D0) then

         I_mean = sum(prof) / NF

!         write(10, 3) rwlae, I_mean
! 3       format (F10.4," ",e12.5)

!         write(10, 3) I_mean
! 3       format (e12.5)

      do  I=NF,1,-1

         write (10, 3) DLAMBD(i), prof(i)
 3       format (F10.4," ",e12.5)

      enddo

      else
      DO  I=NF,1,-1
!         WRITE (10,4) DLAMBD(i),prof(i)
! 4       format (F16.4," ",e12.5)

!*******************************************************************************************************************************************
!Changes by RINAT TAGIROV for radio calculations. One can easily revert to the former version of the output using the two lines right above.

         Wavelength = DLAMBD(i) / 1D+8

         WRITE (10,4) Wavelength, prof(i)
 4       format (F8.5,2x,e12.5)

!*******************************************************************************************************************************************

      enddo
      endif

      close (10)
      RETURN
      END subroutine
      end module
