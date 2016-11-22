      !IMPLICIT REAL*8 (A-H, O-Z)
C
C     Parameters that specify dimensions of arrays
C

      integer,PARAMETER :: MATEX = 90,
     *           MIOEX = 90,
c     *           MIOEX = 30,
     *           MLEV  = 200,
     *           MDEPTH= 200,
     *           MFREQ = 50000,
     *           MFCONT= 2+MFREQ/50,       !max number of continous frequency points FREQC(MFCONT)
     *           MCROSS= MLEV,
     *           MATOM = 90,
c     *           MATOMN = 30,
     *           MION  = 90,
c     *           MATOM = 30,
c     *           MION  = 30,
     *           MION0 =  9,
     *           MI1   = MION-1,
     *           ibuff = 15
CMH   CHANGED BY MARGIT HABERREITER
      INTEGER NDPMIN
CMH  
C     Basic physical constants
C      
      real*8,PARAMETER :: CL    = 2.99792458D10,
     *           H     = 6.6256D-27,        
     *           BOLK  = 1.38054D-16,   
     *           HK    = 4.79928144D-11,
     *           EH    = 2.17853041D-11,
     *           BN    = 1.4743D-2,     
     *           SIGEL = 6.6516D-25,   
     *           PI4H  = 1.8966D27,   
     *           HMASS = 1.67333D-24
C
C     Basic parameters
C
      integer NATOM0,NION,NLEVE0,ND,NFREQ
      COMMON/BASNUM/NATOM0,
     *              NION, 
     *              NLEVE0,
     *              ND,    
     *              NFREQ
      real*8 TEFF0,GRAV,YTOT,WMM,WMY
      COMMON/INPPAR/TEFF0, 
     *              GRAV,
     *              YTOT,
     *              WMM,
     *              WMY
      integer IMODE,IFREQ,INLTE,IDSTD
      COMMON/BASICM/IMODE,
     *              IFREQ,
     *              INLTE,
     *              IDSTD 
C
C     Parameters for explicit atoms
C
      real*8 AMASS,ABUN
      integer NUMAT,N0A,NKA,IATREF
      COMMON/ATOPAR/AMASS(MATEX),
     *              ABUN (MATEX),
     *              NUMAT(MATEX),
     *              N0A(MATEX), 
     *              NKA(MATEX),
     *              IATREF
C
C     Parameters for explicit ions
C
      real*8 FF

      integer NFIRST,NLAST,NNEXT,IZ,IFREE

      COMMON/IONPAR/FF(MIOEX), 
     *              NFIRST(MIOEX),
     *              NLAST(MIOEX),
     *              NNEXT(MIOEX),
     *              IZ(MIOEX),
     *              IFREE(MIOEX)
C
C     Parameters for explicit levels
C
      real*8 ENION,G
      integer NQUANT,IATM,IEL,ILK,ifwop
      COMMON/LEVPAR/ENION(MLEV),
     *              G(MLEV), 
     *              NQUANT(MLEV),
     *              IATM(MLEV), 
     *              IEL(MLEV), 
     *              ILK(MLEV),
     *              ifwop(mlev) 
C
      CHARACTER*10 TYPLEV(MLEV)
      COMMON/LEVLAB/TYPLEV
C
C     Parameters for all considered transitions
C
      real*8 S0BF,ALFBF,BETBF,GAMBF
      integer IBF
      COMMON/TRAPAR/S0BF(MLEV),
     *              ALFBF(MLEV),
     *              BETBF(MLEV),
     *              GAMBF(MLEV),
     *              IBF(MLEV)
C
C     Auxiliary parameters
C
      !*** NKH: Index of Hydrogen in POPNUM?
      integer IATH,IELH,IELHM,N0H,N1H,NKH,N0HN,N0M,IATHE,IELHE1,IELHE2
      integer IATC,IATMG,IATAL,IATSI,IATFE,IELC1,IELMG1,IELAL1,IELSI1
      integer IATFE1
      COMMON/AUXIND/IATH,IELH,IELHM,N0H,N1H,NKH,N0HN,N0M,
     *              IATHE,IELHE1,IELHE2,
     *              IATC,IATMG,IATAL,IATSI,IATFE,
     *              IELC1,IELMG1,IELAL1,IELSI1,IATFE1
C
C     Parameters for atoms considered in line blanketing opacity
C
      LOGICAL LGR(MATOM),LRM(MATOM)
      real*8 PFSTD,RR,ENEV,AMAS,ABND,TYPAT,RELAB
      integer MODPF,INPOT,NATOMS,IONIZ
      COMMON/PFSTDS/PFSTD(MION,MATOM),MODPF(MATOM)
      COMMON/ADDPOP/RR(MATOM,MION),INPOT(MATOM,MION)
      COMMON/ATOBLN/ENEV(MATOM,MI1),AMAS(MATOM),ABND(MATOM),TYPAT(MATOM)
      COMMON/ATOINI/NATOMS,IONIZ(MATOM),LGR,LRM
      COMMON/POPCHN/RELAB(MLEV)
      COMMON/NDTMIN/NDPMIN
