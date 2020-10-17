      module MOD_INTRFC

      contains

      subroutine intrfc(ncharg,weight,elevel,eion,
     *                  einst,alpha,sexpo,agaunt,natom,
     *                  symbol,nfirs0,nlast0,
     *                  WAVARR,SIGARR,N,NFDIM)
!     fort.55 is read here
!     now includes fix of quantum number of HeI lines

      use MOD_INIBL
      use MOD_SYNSUBM
      use MOD_STATEQ

      implicit none

      integer,intent(in) :: n,ncharg,natom,
     *                      nfirs0,nlast0,NFDIM

      real*8, intent(in) :: weight,elevel,eion,
     *                      einst,alpha,sexpo,
     *                      WAVARR,SIGARR

      character*8, intent(in) :: agaunt(N)
      character*2, intent(in) :: SYMBOL(NATOM)

      integer :: INCODE,ICHEMC,IOPADD,IRSCT,ION,NA,IAT,I,II,IA
      integer :: IOPHMI,IOPH2P,IOPHE1,IOPHE2,IOPFE1,IOPHLI
      integer :: JJ,ID,INMOD,INTRPL,ICHANG,NCH
      real*8  :: SIG

      real*8 :: inibl_end_1, inibl_start_1, inibl_end_2, inibl_start_2

      INCLUDE '../inc/PARAMS.FOR'
      INCLUDE '../inc/MODELP.FOR'

      DIMENSION NCHARG (N), WEIGHT(N),ELEVEL(N)
      DIMENSION EION(N),EINST(N, N)

      DIMENSION ALPHA(N), SEXPO(N)
      DIMENSION NFIRS0(NATOM), NLAST0(NATOM)
      DIMENSION WAVARR(N, NFDIM), SIGARR(N, NFDIM)

      COMMON/INTKEY/INMOD,INTRPL,ICHANG,ICHEMC
      COMMON/OPCPAR/IOPADD,IOPHMI,IOPH2P,IRSCT,IOPHLI,IOPHE1,IOPHE2,IOPFE1

      open(555, file = 'sun.inp', STATUS='OLD')
      read(555, *) incode

!     original input - routine START and others

      if (incode .eq. 0) then

         print*, 'original form of input'

         call inimod
         call tint

!     changed by Margit Haberreiter
!     *****************************
         call cpu_time(inibl_start_1)

         call inibl(wavarr, sigarr, n, nfdim)

         call cpu_time(inibl_end_1)

         print*, 'intrfc: inibl time 1 = ', inibl_end_1 - inibl_start_1
!     *****************************

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

         vturb(id)=vturb(id)*vturb(id)
   60 continue

      call inimod

      call tint

      call cpu_time(inibl_start_2)

      call inibl(WAVARR, SIGARR, N, NFDIM)

      call cpu_time(inibl_end_2)

      print*, 'intrfc: inibl time 2 = ', inibl_end_2 - inibl_start_2

      call hylset

      call he2set

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

  614 format(5i4)
  616 FORMAT(1H ,i3,1PD15.7,0PF5.2,4I3,i5,1x,1p4e9.2)

      return

      end subroutine

      end module
