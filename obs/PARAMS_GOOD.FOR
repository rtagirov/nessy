      IMPLICIT REAL*8 (A-H, O-Z)
C
C     Parameters that specify dimensions of arrays
C
      PARAMETER (MATEX = 20,   
C CHANGED BY MARGIT HABERREITER
     *           FGARR = 2000,
     *           NDIM  = 55, 
     *           MIOEX = 30,
     *           MLEV  = 200,
     *           MDEPTH= 200,
C     *           MDEPTH= 75,
c     *           MFREQ = 5300,
c     *		 MFREQ = 8019,
     *	         MFREQ = 50000,
     *           MFCONT= 2,
     *           MCROSS= MLEV,
     *           MATOM = 30,
     *           MION  = 30,
     *           MION0 =  9,
     *           MI1   = MION-1)
      parameter (ibuff = 15)
C CHANGED BY MARGIT HABERREITER
C      DIMENSION WAVARR(NFDIM),SIGARR(NDIM,NFDIM)	
      DIMENSION WAVARR(37,2000),SIGARR(37,2000)	
c      DIMENSION WAVARR(100,2000),SIGARR(100,2000)
C     Basic physical constants
C      
      PARAMETER (CL    = 2.997925D10,
     *           H     = 6.6256D-27,        
     *           BOLK  = 1.38054D-16,   
     *           HK    = 4.79928144D-11,
     *           EH    = 2.17853041D-11,
     *           BN    = 1.4743D-2,  
     *           SIGEL = 6.6516D-25,   
     *           PI4H  = 1.8966D27,   
     *           HMASS = 1.67333D-24)
C
C     Basic parameters
C
      COMMON/BASNUM/NATOM0,
     *              NION, 
     *              NLEVE0,
     *              ND,    
     *              NFREQ,
     *              NFCONT
      COMMON/INPPAR/TEFF0, 
     *              GRAV,
     *              YTOT,
     *              WMM,
     *              WMY 
      COMMON/BASICM/IMODE,
     *              IFREQ,
     *              INLTE,
     *              IDSTD 
C
C     Parameters for explicit atoms
C
      COMMON/ATOPAR/AMASS(MATEX),
     *              ABUN (MATEX),
     *              NUMAT(MATEX),
     *              N0A(MATEX), 
     *              NKA(MATEX),
     *              IATREF
C
C     Parameters for explicit ions
C
      COMMON/IONPAR/FF(MIOEX), 
     *              NFIRST(MIOEX),
     *              NLAST(MIOEX),
     *              NNEXT(MIOEX),
     *              IZ(MIOEX),
     *              IFREE(MIOEX)
C
C     Parameters for explicit levels
C
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
      COMMON/TRAPAR/S0BF(MLEV),
     *              ALFBF(MLEV),
     *              BETBF(MLEV),
     *              GAMBF(MLEV),
     *              IBF(MLEV)
C
C     Auxiliary parameters
C
      COMMON/AUXIND/IATH,IELH,IELHM,N0H,N1H,NKH,N0HN,N0M,
     *              IATHE,IELHE1,IELHE2,
     *              IATC,IATMG,IATAL,IATSI,IATFE,
     *              IELC1,IELMG1,IELAL1,IELSI1,IATFE1    
C
C     Parameters for atoms considered in line blanketing opacity
C
      LOGICAL LGR(MATOM),LRM(MATOM)
      COMMON/PFSTDS/PFSTD(MION,MATOM),MODPF(MATOM)
      COMMON/ADDPOP/RR(MATOM,MION),INPOT(MATOM,MION)
      COMMON/ATOBLN/ENEV(MATOM,MI1),AMAS(MATOM),ABND(MATOM),TYPAT(MATOM)
      COMMON/ATOINI/NATOMS,IONIZ(MATOM),LGR,LRM
      COMMON/POPCHN/RELAB(MLEV)
