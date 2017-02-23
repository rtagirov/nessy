      module MOD_INTRFC_M

      contains

      subroutine intrfc_m(ncharg,weight,elevel,eion,
     *                    einst,alpha,sexpo,agaunt,natom,
     *                    symbol,nfirs0,nlast0,
     *                    WAVARR,SIGARR,N,NF)

!     now includes fix of quantum number of HeI lines

      use MOD_INIBL0
      use MOD_SYNSUBM
      use MOD_STATEQ

      implicit none

      integer,intent(in) :: n,ncharg,natom,
     *                      nfirs0,nlast0,NF

      real*8, intent(in) :: weight,elevel,eion,
     *                      einst,alpha,sexpo,
     *                      WAVARR,SIGARR

      character*8, intent(in) :: agaunt(N)
      character*2, intent(in) :: SYMBOL(NATOM)

      integer :: INCODE,ICHEMC,IOPADD,IRSCT,ION,NA,IAT,I,II,IA
      integer :: IOPHMI,IOPH2P,IOPHE1,IOPHE2,IOPFE1,IOPHLI
      integer :: JJ,ID,INMOD,INTRPL,ICHANG,NCH
      real*8  :: SIG

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'

      DIMENSION NCHARG (N), WEIGHT(N),ELEVEL(N)
      DIMENSION EION(N),EINST(N, N)

      DIMENSION ALPHA(N), SEXPO(N)
      DIMENSION NFIRS0(NATOM), NLAST0(NATOM)
      DIMENSION WAVARR(N, NF), SIGARR(N, NF)

      COMMON/INTKEY/INMOD,INTRPL,ICHANG,ICHEMC
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,IOPFE1

      open(555, file = 'input_sun', STATUS='OLD')
      read(555, *) incode

C     original input - routine START and others

      if(incode.eq.0) then
         print *,' ***** original form of input!'

         call inimod
         call tint

C**** CHANGED BY MARGIT HABERREITER
         call inibl0(WAVARR, SIGARR, N, NF)
C***********************************
         call hylset
         call he2set
         call inibla
         return
      end if
C
C     New input - atomic data obtained form Werner's input
C
      grav=4.
      READ(55,*) IMODE,IDSTD,ICHEMC
      call state0
      READ(555,*) IOPADD
      IF(IOPADD.GT.0) THEN
C***  CHANGES BY MARGIT HABERREITER
C***  ADDITIONAL OPACITY SOURCES READ FROM INPUT_SUN (formerly INPUT5)
c         READ(555,*) IRSCT,IOPHMI,IOPH2P,IOPHE1,IOPHE2
        READ(555,*) IRSCT,IOPHMI,IOPH2P,IOPHE1,IOPHE2,IOPFE1
      END IF
      IOPHLI=0
      WRITE(6,605) IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPFE1
  605 FORMAT(1H0//' ADDITIONAL OPACITY SOURCES'/
     *            ' --------------------------'/
     * ' IOPADD                             =',I3/
     * ' IOPHMI  (H-  OPACITY IN LTE)       =',I3/
     * ' IOPH2P  (H2+  OPACITY)             =',I3/
     * ' IRSCT   (RAYLEIGH SCATTERING)      =',I3/
     * ' IOPHLI  (LYMAN LINES WINGS)        =',I3/
     * ' IOPFE1  (EXPLICIT OPACITY SOURCES) =',I3)
C
C     data for explicit atoms and explicit ions
C
      ion=0
      iath=0
      do 20 na=1,natom
        iat=na
        i=0
        if(symbol(na).eq.'H ') numat(iat)=1
        if(symbol(na).eq.'HE') numat(iat)=2
        if(symbol(na).eq.'Li') numat(iat)=3
        if(symbol(na).eq.'Be') numat(iat)=4
        if(symbol(na).eq.'B ') numat(iat)=5
        if(symbol(na).eq.'C ') numat(iat)=6
        if(symbol(na).eq.'N ') numat(iat)=7
        if(symbol(na).eq.'O ') numat(iat)=8
        if(symbol(na).eq.'F ') numat(iat)=9
        if(symbol(na).eq.'Ne') numat(iat)=10
        if(symbol(na).eq.'Na') numat(iat)=11
        if(symbol(na).eq.'NA') numat(iat)=11
        if(symbol(na).eq.'Mg') numat(iat)=12
        if(symbol(na).eq.'MG') numat(iat)=12
        if(symbol(na).eq.'Al') numat(iat)=13
        if(symbol(na).eq.'AL') numat(iat)=13
        if(symbol(na).eq.'Si') numat(iat)=14
        if(symbol(na).eq.'SI') numat(iat)=14
        if(symbol(na).eq.'P ') numat(iat)=15
        if(symbol(na).eq.'S ') numat(iat)=16
        if(symbol(na).eq.'Cl') numat(iat)=17
        if(symbol(na).eq.'Ar') numat(iat)=18
        if(symbol(na).eq.'K ') numat(iat)=19
        if(symbol(na).eq.'Ca') numat(iat)=20
        if(symbol(na).eq.'CA') numat(iat)=20
        if(symbol(na).eq.'Sc') numat(iat)=21
        if(symbol(na).eq.'Ti') numat(iat)=22
        if(symbol(na).eq.'V ') numat(iat)=23
        if(symbol(na).eq.'Cr') numat(iat)=24
        if(symbol(na).eq.'Mn') numat(iat)=25
        if(symbol(na).eq.'Fe') numat(iat)=26
        if(symbol(na).eq.'FE') numat(iat)=26
        if(symbol(na).eq.'Co') numat(iat)=27
        if(symbol(na).eq.'Ni') numat(iat)=28
        if(symbol(na).eq.'Cu') numat(iat)=29
        if(symbol(na).eq.'Zn') numat(iat)=30

c*************************************************
        ia=numat(iat)
c         if(ia.eq.0) call stop(' conflict in interface - numat')
        if(ia.eq.0) then
            print *, 'conflict in interface - numat'
            stop
        end if
        abun(iat)=abnd(ia)
        amass(iat)=amas(ia)*1.67333e-24
        n0a(iat)=nfirs0(na)
        nka(iat)=nlast0(na)
        print *,' iat=',iat,abun(iat),n0a(iat),nka(iat)
        nch=-2
        do 10 ii=n0a(iat),nka(iat)
            iatm(ii)=iat
            ilk(ii)=0
            if(ncharg(ii).gt.nch) then
              nch=ncharg(ii)
              if(ii.lt.nka(iat)) then
                  ion=ion+1
                  iz(ion)=nch+1
                  nfirst(ion)=ii
                  ff(ion)=0.
              end if
              if(ii.gt.n0a(iat)) then
                  if(ii.lt.nka(iat)) then
                    nlast(ion-1)=ii-1
                    nnext(ion-1)=ii
                    ilk(ii)=ion-1
                  else
                    nlast(ion)=ii-1
                    nnext(ion)=ii
                    ilk(ii)=ion
                    iel(ii)=ion
                  end if
              end if
              ifree(ion)=1
              if(ia.eq.1) iath=iat
              if(ia.eq.2) iathe=iat
            end if
  10    continue
  20  continue
      iatref=2
      if(iath.gt.0) iatref=1
C
      nion=ion
      WRITE(6,613)
      do 25 ion=1,nion
         if(iz(ion).eq.1.and.iatm(nfirst(ion)).eq.iath) then
           ielh=ion
        PRINT *, 'IELH=',IELH
           ifree(ielh)=2
         end if
         if(iatm(nfirst(ion)).eq.iathe) then
            if(iz(ion).eq.1) ielhe1=ion
            if(iz(ion).eq.2) then
              ielhe2=ion
              ifree(ielhe2)=2
            end if
         end if
         do 22 ii=nfirst(ion),nlast(ion)
            iel(ii)=ion
   22    continue
         print *, nfirst(ion),nlast(ion),nnext(ion),iz(ion)
         WRITE(6,614) nfirst(ion),nlast(ion),nnext(ion),iz(ion)
   25 continue
C
C     main quantum numbers for hydrogen levels
C
      if(ielh.gt.0) then
      do 30 ii=nfirst(ielh),nlast(ielh)
         nquant(ii)=ii-nfirst(ielh)+1
   30 continue
      n0hn=nfirst(ielh)
      n0h=n0a(iath) 
      n1h=nlast(ielh)
      nkh=nnext(ielh)
      n0m=0
      end if
C
C     main quantum numbers for He I levels
C
      if(ielhe1.gt.0) then
      do 35 ii=nfirst(ielhe1),nlast(ielhe1)
         if (elevel(ii).lt. 500.) then
            nquant(ii)=1
         else if (elevel(ii).lt.175000.) then
            nquant(ii)=2
         else if (elevel(ii).lt.189000.) then
            nquant(ii)=3
         else if (elevel(ii).lt.192000.) then
            nquant(ii)=4
         else if (elevel(ii).lt.194500.) then
            nquant(ii)=5
         else if (elevel(ii).lt.195600.) then
            nquant(ii)=6
         else if (elevel(ii).lt.196300.) then
            nquant(ii)=7
         else if (elevel(ii).lt.196700.) then
            nquant(ii)=8
         else if (elevel(ii).lt.197000.) then
            nquant(ii)=9
         else if (elevel(ii).lt.197300.) then
            nquant(ii)=10
         else if (elevel(ii).lt.197470.) then
            nquant(ii)=11
         else if (elevel(ii).lt.197600.) then
            nquant(ii)=12
         else if (elevel(ii).lt.197700.) then
            nquant(ii)=13
         else if (elevel(ii).lt.197800.) then
            nquant(ii)=14
         else if (elevel(ii).lt.197850.) then
            nquant(ii)=15
         else if (elevel(ii).lt.197900.) then
            nquant(ii)=16
         else if (elevel(ii).lt.197950.) then
            nquant(ii)=17
         else if (elevel(ii).lt.197990.) then
            nquant(ii)=18
         else if (elevel(ii).lt.198020.) then
            nquant(ii)=19
         else if (elevel(ii).lt.198040.) then
            nquant(ii)=20
         else
            print *,' nfirst(ielhe1),nlast(ielhe1) ',
     $           nfirst(ielhe1),nlast(ielhe1) 
            print *,' ii,elevel(ii) ',ii,elevel(ii)
            print *,' main quantum number of Hei not yet coded'
            stop 'update INTRFC'
         endif
 35   continue
      endif
C
C     main quantum numbers for He II levels
C
      if(ielhe2.gt.0) then
      do 40 ii=nfirst(ielhe2),nlast(ielhe2)
         nquant(ii)=ii-nfirst(ielhe2)+1
   40 continue
      end if
C
C     data for explicit levels and continua
C     1.9857e-16: h*c in erg*s*cm/s
      nleve0=n
      do 50 ii=1,nleve0
         enion(ii)=(eion(ii)-elevel(ii))*1.9857e-16
         g(ii)=weight(ii)
         if(nquant(ii).eq.0) nquant(ii)=1
c         if(iel(ii).eq.ilk(ii)) go to 45
         jj=nnext(iel(ii))
         sig=einst(ii,jj)
         if(ii.ne.jj.and.sig.le.0.) then
            write(6,601) ii,jj
  601       format(' *** continuum transition',2i4,' has zero c.-s.')
            go to 50
         end if
         if(iel(ii).eq.ielh) then
            ibf(ii)=1
          else if(iel(ii).eq.ielhe2) then
            ibf(ii)=1 
C*********************************************************
C***  MARGIT HABERREITER
CMH   iel*** extended for additional elements  
CMH     else if(iel(ii).eq.ielc) then
CMH            ibf(ii)=1 
c*********************************************************
          else
            if(agaunt(ii).eq.'KOESTER ') then
               ibf(ii)=5
               s0bf(ii)=sexpo(ii)
               alfbf(ii)=alpha(ii)
               betbf(ii)=sig
CMH Keyword 'SEATON' now used as in hminus-part
CMH used as hydrogenic with explicit Gaunt factor
                 else if(agaunt(ii).eq.'SEATON  ') then
                      ibf(ii)=2
                  s0bf(ii)=3.
                 alfbf(ii)=sig
                 betbf(ii)=1.
c                ibf(ii)=6
c                    s0bf(ii)=sexpo(ii)
c                    alfbf(ii)=sig
c           print *,'SEATON, INTRFC' 
c           print *,ii,ibf(ii),sexpo(ii),sig,nquant(ii)
c           pause
CMH  IFB IS THE SWITCH CONTROLLING THE MODE OF EVALUATION OF THE 
CMH  CROSS SECTION
CMH  IFB = 9: OPACITY PROJECT XS-FITS (INTERPOLATIONS)
          else if(agaunt(ii).eq.'TABLE  ') then
          ibf(ii) = 9
          else
               ibf(ii)=2
               s0bf(ii)=3.
               alfbf(ii)=sig
               betbf(ii)=1.
            end if
         end if
   45    WRITE(6,616) ii,ENION(Ii),G(Ii),NQUANT(Ii),
     *                IEL(Ii),ILK(Ii),IATM(Ii)
     *               ,ibf(ii),s0bf(ii),alfbf(ii),betbf(ii),gambf(ii)
   50 continue
      WRITE(6,603) INMOD,ND,IDSTD,INTRPL,ICHANG,
     *             NATOM,NION,NLEVE0,
     *             IELH,IELHM,IATH

      do 60 id=1,nd
       !   print*, 'vt test', id, vturb(id)
         vturb(id)=vturb(id)*vturb(id)
   60 continue

      print *,' call inimod'
      call inimod
      print *,' call tint'
      call tint
      print *,' call inibl0(WAVARR,SIGARR)'

C***  CHANGED BY MARGIT HABERREITER
      call inibl0(WAVARR, SIGARR, N, NF)
C***********************************
      print *,' call hylset'
      call hylset
      print *,' call he2set'
      call he2set
      print *,' call inibla'
      call inibla
  603 FORMAT(1H0//' BASIC INPUT PARAMETERS'/
     *            ' ----------------------'/
     *            ' INMOD  =',I5/
     *            ' ND     =',I5/
     *            ' IDSTD  =',I5/
     *            ' INTRPL =',I5/
     *            ' ICHANG =',I5/
     *            ' NATOM  =',I5/
     *            ' NION   =',I5/
     *            ' NLEVEL =',I5/
     *            ' IELH   =',I5/
     *            ' IELHM  =',I5/
     *            ' IATH   =',I5)
  613 FORMAT(1H0//' EXPLICIT IONS INCLUDED'/
     *            ' ----------------------'//
     *  ' ION     N0    N1    NK    IZ  IUPSUM ICUP       FF'/)
c  614 FORMAT(1H ,A4,6I6,1PD15.3)
  614 format(5i4)
  616 FORMAT(1H ,i3,1PD15.7,0PF5.2,4I3,i5,1x,1p4e9.2)
      return
      end subroutine
c
c
C
      subroutine quant (nddim,n,t,rne,entot,popnum,vdop,vdopp,
     $                  nvopa,xlam,Npot,Tion_pot,dil,teff,xjc)
      use UTILS
      use OPINT,only: VOPA
      use SYNTHP_CONT,only:FREQC,ABSOC,EMISC,SCATC,NFCONT,setNFCONT
      use SYNTHP_CONT,only:SYNTHP_CONT_INIT

c called by 
c setup model quantities for Ivan's routines
c 
      implicit none
      integer,intent(in   ) :: nddim,n,
     $                  nvopa,Npot
      real*8, intent(in   ) :: t,rne,entot,popnum,vdop,vdopp,
     $                  xlam,
     $                  Tion_pot,dil,teff,xjc(:) !XJC(NDDIM)
      !integer,intent(inout) ::
      !real*8, intent(inout) ::
      real*8  :: WW,XX,YY,ZZ, TURB,VELORMS, cl8,clkm,VEL
      integer :: I,L,LEV,KOPA
      integer :: nfcont_
      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/MODELP.FOR'
      DIMENSION POPNUM(nddim,N),ENTOT(nddim),rne(nddim),t(nddim)
     $         ,Tion_pot(nddim,Npot),dil(nddim)
      DIMENSION VELORMS(NDDIM),TURB(NDDIM)
      LOGICAL ::  ADDASPL,ADDTURB
CMH   ADDASPL = .TRUE. THEN THE TURBULENT VELOCITY FROM
CMH   ASPLUND, ET AL. 2000, A&A, 359, P. 729 IS ICLUDED
CMH   ADDTURB = TRUE THEN THE DEPTH DEPENDEND TURBULENCE
CMH   BROADENING IN THE FONTENLA ATMOSPHERES ARE INLCUDED
      ADDASPL = .FALSE.
c     ADDASPL = .TRUE.
      ADDTURB = .FALSE.
c     ADDTURB = .TRUE.

      DO L=1,NDDIM
        TURB(L)=0.
      ENDDO
      if (nddim.gt.MDEPTH) then
         print *,' MDEPTH .lt. nd'
         stop ' adjust grid dimensions - nd'
      endif
      if (nvopa.gt.MFREQ) then
         print *,'nvopa,MFREQ',nvopa,MFREQ
         print *,' MFREQ  .lt. nvopa' 
         stop ' adjust grid dimensions - nvop'
      endif
cif (n.ne.nlevel) then
c     print *,' n .ne. nlevel' 
cstop ' adjust grid dimensions - nlevel'
cendif
C***  READING VELOCITY ASPLUND 2000, A&A 359, 729
      IF (ADDASPL) THEN
      open (unit=9999,file='VELO',STATUS='OLD')
      DO 9999 L=1,NDDIM
        velorms(l)=0.
        READ (UNIT=9999, fmt=*), XX,yy,VELORMS(L)
        write(6,*) L,'quant: VELORMS= ',velorms(L)
 9999 continue
      CLOSE (9999)
      ENDIF

CMH   IF ADDTURB TRUE THEN READ THE TURBULENCE BROADENING FROM FAL_VD
      IF (ADDTURB) THEN
      open (UNIT=9999,file='FAL_VD',STATUS='OLD')
      DO L=1,NDDIM
        READ (UNIT=9999, fmt=*),WW, XX,YY,ZZ,TURB(L)
        PRINT *, L,'FIOSS: TURB= ',TURB(L)
        write (6,*) L,'FIOSS: TURB= ',TURB(L)

      ENDDO
      CLOSE (9999)
      ENDIF
C**********************************************************************
      do l=1,nddim
         temp(l)   = t(l)
         elec(l)   = rne(l)*entot(l)
C***************************************************************************
C***  CHANGES by Margit Hhaberreiter
C***  22.8.03: replace expansion velocity by tabulated 
C*** (horizontally averaged) rms vertical velocity
C***  Reference: Asplund et al. 2000, 358, 729-742
C***  dann noch ganz frech die Turbulenz mit Faktor multipliziert ...
c        vturb(l)  = velorms(L)*1.e5*1.4
      IF (ADDASPL) THEN
         vturb(l)  = velorms(L)*1.e5
c        vturb(l)  = velorms(L)*1.e5
        PRINT*,'INTRFC: ORIGINAL ASPLUND VELOCITY ACTIVE, no tuning!!'
      ElSE IF (ADDTURB) THEN
CMH   ATTENTION:
CMH   IN THE FONTENAL ATMOSPHERES V_T(TURBULENCE BROADENING) 
CMH   AND V_P (THERMAL BROADENING)
c     vturb(l) = vturb(l)+(TURB(l)*1.e5*sqrt(2.))
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
      nd=nddim
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


      subroutine synopa(WAVARR, SIGARR, N, NF)
      use MOD_SYNSUBM
      use MOD_INIBL0  !contains OPAC
      use SYNTHP_CONT
      use OPINT,  only: OPATOT,ETATOT
      
c     Interface opacity routine

      implicit none
      integer, intent(in) :: N, NF
      real*8,  intent(in) :: WAVARR,SIGARR
      real*8 :: ABSO,EMIS
    
      integer :: ID,IJ, i

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/SYNTHP.FOR'
      INCLUDE '../inc/MODELP.FOR'

      DIMENSION WAVARR(N, NF), SIGARR(N, NF)

CMH   MARGIT HABERREITER
CMH   IJ:    INDEX OF FREQUENCT
CMH   NFREQ: NUMBER OF FREQUENCY POINTS
CMH   INPUT PARAMETERS:

CMH   OUTPUT PARAMETERS:
CMH   OPATOT: NEW ARRAY FOR OPACITY for each frequency nfreq and depth point nd
CMH   ETATOT: NEW ARRAY FOR EMISSIVITY for each frequency nfreq and depth point nd
      dimension abso(mfreq),emis(mfreq) !,scat(2)

       if (.not. allocated(wav_oi)) then
       allocate(wav_oi(20))
       allocate(opac_oi(20))
       allocate(opac_f(20, ND))
       endif


       do 20 id=1,nd

         call opac(id, 1, abso, emis, WAVARR, SIGARR, N, NF)

         opac_f(1:20,ND+1-id)=opac_oi(1:20)

C     loop over frequency points
         do ij=1,nfreq

cws May-24-96: Ivan's routine calculated low to high freqency
c         it is (in QUANT): 
c                wlam(nvopa-kopa+1) = xlam * (1. - vel*clkm)
c                freq(nvopa-kopa+1) = cl8/wlam(nvopa-kopa+1)
            opatot(nfreq-ij+1,id)=abso(ij)
            etatot(nfreq-ij+1,id)=emis(ij)

            if (id.le.1) then
               etatot(nfreq-ij+1,id)=emisc(getContIdx(ij))
            endif


CMH   print line opacity to file
cpr_opa          write(200,310) opatot(nfreq-ij+1,id),etatot(nfreq-ij+1,id)
         enddo

c        write(91,691) id
c        write(91,692) (abso(ij),ij=1,nfreq)

c        write(91,692) (abso(ij),ij=213,216)
c        write(92,691) id
c        write(92,692) (emis(ij),ij=1,nfreq)
c 691    format(i5)
c 692    format(1p8e10.2)
   20 continue
  310 format (1pe12.5,1pe12.5)


!*       open (unit=100, file="../data_opac.txt",access="APPEND") 

!*      do i=20, 1, -1
!*       write(100, *), '**********************' 
!*       write(100, *), airlambda(wav_oi(i))
!*       write(100, '(81e12.5)'), opac_f(i,1:81)

!*       enddo
!*       close(100)
      
      return
      end subroutine
      end module

