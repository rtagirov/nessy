      module mod_quant

      contains

      subroutine quant(dpn,n,t,rne,entot,popnum,vdop,vdopp,
     $                 nvopa,xlam,Npot,Tion_pot,dil,teff,xjc)
      use UTILS
      use OPINT,only: VOPA
      use SYNTHP_CONT,only:FREQC,ABSOC,EMISC,NFCONT,setNFCONT
      use SYNTHP_CONT,only:SYNTHP_CONT_INIT

!     setup model quantities for Ivan's routines

      implicit none
      integer,intent(in) :: dpn,n,nvopa,Npot
      real*8, intent(in) :: t,rne,entot,popnum,vdop,vdopp,xlam,
     $                      Tion_pot,dil,teff,xjc(:)

      real*8  :: WW,XX,YY,ZZ, TURB,VELORMS, cl8,clkm,VEL
      integer :: I,L,LEV,KOPA
      integer :: nfcont_
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/MODELP.FOR'
      DIMENSION POPNUM(dpn, N), ENTOT(dpn),rne(dpn),t(dpn)
     $         ,Tion_pot(dpn, Npot), dil(dpn)
      DIMENSION VELORMS(dpn), TURB(dpn)
      LOGICAL ::  ADDASPL,ADDTURB
CMH   ADDASPL = .TRUE. THEN THE TURBULENT VELOCITY FROM
CMH   ASPLUND, ET AL. 2000, A&A, 359, P. 729 IS ICLUDED
CMH   ADDTURB = TRUE THEN THE DEPTH DEPENDEND TURBULENCE
CMH   BROADENING IN THE FONTENLA ATMOSPHERES ARE INLCUDED
      ADDASPL = .FALSE.

      ADDTURB = .FALSE.


      ND = DPN

      DO L = 1, ND
        TURB(L)=0.
      ENDDO
      if (nd .gt. MDEPTH) then
         print *,' MDEPTH .lt. nd'
         stop ' adjust grid dimensions - nd'
      endif
      if (nvopa.gt.MFREQ) then
         print *,'nvopa,MFREQ',nvopa,MFREQ
         print *,' MFREQ  .lt. nvopa' 
         stop ' adjust grid dimensions - nvop'
      endif




C***  READING VELOCITY ASPLUND 2000, A&A 359, 729
      IF (ADDASPL) THEN
      open (unit=9999,file='VELO',STATUS='OLD')
      DO 9999 L = 1, ND
        velorms(l)=0.
        READ (UNIT=9999, fmt=*), XX,yy,VELORMS(L)
        write(6,*) L,'quant: VELORMS= ',velorms(L)
 9999 continue
      CLOSE (9999)
      ENDIF

CMH   IF ADDTURB TRUE THEN READ THE TURBULENCE BROADENING FROM atm.inp
      IF (ADDTURB) THEN 
! careful here
! for muram/kurucz format atmosphere or
! if the atmosphere is read from a muram cube slice
! this bit is not going to work
! it works only for FAL format atmosphere

      open (UNIT = 9999, file = 'atm.inp', STATUS = 'OLD')

      DO L=1, ND
        READ (UNIT=9999, fmt=*),WW, XX,YY,ZZ,TURB(L)
        PRINT *, L,'FIOSS: TURB= ',TURB(L)
        write (6,*) L,'FIOSS: TURB= ',TURB(L)

      ENDDO
      CLOSE (9999)
      ENDIF

      do l=1, nd
         temp(l)   = t(l)
         elec(l)   = rne(l)*entot(l)

C***  CHANGES by Margit Hhaberreiter
C***  22.8.03: replace expansion velocity by tabulated 
C*** (horizontally averaged) rms vertical velocity
C***  Reference: Asplund et al. 2000, 358, 729-742
C***  dann noch ganz frech die Turbulenz mit Faktor multipliziert ...

      IF (ADDASPL) THEN
         vturb(l)  = velorms(L)*1.e5

        PRINT*,'INTRFC: ORIGINAL ASPLUND VELOCITY ACTIVE, no tuning!!'
      ElSE IF (ADDTURB) THEN
CMH   ATTENTION:
CMH   IN THE FONTENAL ATMOSPHERES V_T(TURBULENCE BROADENING) 
CMH   AND V_P (THERMAL BROADENING)

c     PRINT *,'FONTENLAS DETPH DEPENDENT BROADENING TAKEN *sqrt(2.)!!!'
      vturb(l) = vturb(l)+(TURB(l)*1.e5)
      ELSE
C***    original:
         vturb(l)  = vdopp*1.e5
c         vturb(l)  = vdopp
        PRINT*,'ASPLUND VELOCITY NOT ACTIVE',vturb(l)
      ENDIF

C***  END OF CHANGES 
C****************************************************************************
         wdil(l)   = dil(l)
c     print *,'2. intrfc: wdil(l) = ', wdil(l)
c****  hier Kappa ausrechnenen!!!
         do lev=1,n
            popul(lev,l) = popnum(l,lev)*entot(l)
         enddo
         trad(1,l)=tion_pot(l,3)
         trad(2,l)=tion_pot(l,1)
         trad(3,l)=tion_pot(l,2)
         xjcon(l)=xjc(l)
      enddo
      cl8 = cl*1.e8
      clkm = vdop*1.e5/cl
      do kopa=1,nvopa
         vel = vopa(kopa)
cws      May-24-96: change from + to - in calculating wlam from vopa
c        loop in FIOSS is from blue=1 to red=nvopa
c        Ivan's routine wants low to high freqency
         wlam(nvopa-kopa+1) = xlam * (1. - vel*clkm)
         freq(nvopa-kopa+1) = cl8/wlam(nvopa-kopa+1)
c         write (90,*) nvopa-kopa+1,wlam(nvopa-kopa+1)
      enddo

      !*** Prepare Continous frequency array
      nfreq=nvopa

      nfcont_=max(2,(nfreq/50)+1)
      if( (nfreq/(nfcont_-1))*(nfcont_-1)==nfreq) nfcont_=nfcont_-1
      call SYNTHP_CONT_INIT(NFCONT_)
      freqc(1)=freq(1)
      forall(i=2:nfcont_-1)
        freqc(i)=freq((nfreq/(nfcont_-1))*(i-1))
      endforall
      freqc(nfcont_)=freq(nfreq)
      call assert(all(freqc(:nfcont_)>0),'FREQC<=0')
      call assert(NFCONT_ <= MFCONT, 'NFCONT > MFCONT!')
      teff0=teff

      teff0=teff 

      return

      end subroutine

      end module
