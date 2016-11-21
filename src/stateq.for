      module MOD_STATEQ 
      contains
c ********************************************************************
C
C
      SUBROUTINE STATE0
C     =================
C
C     Initialization of basic parameters for the Saha equations
C
c      INCLUDE 'PARAMS.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
      COMMON/INTKEY/INMOD,INTRPL,ICHANG,ICHEMC
      DIMENSION D(4,MATOM),XI(9,MATOM),iatex(matom)
      !DIMENSION XIFE(10)
      LOGICAL LY,LW
C
C     Standard atomic constants for the first 30 species:
C
C            Element Atomic  Solar    Std.
C                    weight abundance highest
C

       DATA D/'  H ', 1.008, 1.00E0, 2.,
     *       '  HE', 4.003, 1.00E-1, 3.,
     *       '  LI', 6.941, 3.9E-12, 3.,
     *       '  BE', 9.012, 1.1E-11, 3.,
     *       '  B ',10.810, 6.2E-10, 4.,
     *       '  C ',12.011, 3.68E-4, 5.,
     *       '  N ',14.007, 1.14E-4, 5.,
     *       '  O ',16.000, 6.70E-4, 5.,
     *       '  F ',18.918, 3.63E-8, 4.,
     *       '  NE',20.179, 2.79E-5, 4.,
     *       '  NA',22.990, 1.72E-6, 4.,
     *       '  MG',24.305, 3.43E-5, 4.,
     *       ' AL ',26.982, 2.49E-6, 4.,
     *       ' SI ',28.086, 3.51E-5, 5.,
     *       ' P  ',30.974, 2.67E-7, 5.,
     *       ' S  ',32.060, 1.61E-5, 5.,
     *       ' CL ',35.453, 4.42E-7, 5.,
     *       ' AR ',39.948, 4.42E-6, 5.,
     *       ' K  ',39.098, 1.11E-7, 5.,
     *       ' CA ',40.080, 2.12E-6, 5.,
     *       ' SC ',44.956, 1.16E-9, 5.,
     *       ' TI ',47.900, 5.44E-8, 5.,
     *       ' V  ',50.941, 1.25E-8, 5.,
     *       ' CR ',51.996, 4.96E-7, 5.,
     *       ' MN ',54.938, 1.57E-7, 5.,
     *       ' FE ',55.847, 2.49E-5, 6.,
     *       ' CO ',58.933, 3.13E-8, 5.,
     *       ' NI ',58.700, 1.89E-6, 5.,
     *       ' CU ',63.546, 2.79E-8, 5.,
     *       ' ZN ',65.380, 2.60E-8, 5.,
CMH Abundances for GA to TH from Anders & Grevesse
     *       ' GA ',69.723, 7.58E-10,5.,
     *       ' GE ',72.61 , 2.57E-09,5.,
     *       ' AS ',74.922, 2.34E-10,5.,
     *       ' SE ',78.96 , 2.23E-09,5., 
     *       ' BR ',79.904, 4.26E-10,5.,
     *       ' KR ',83.80 , 1.69E-09,5.,
     *       ' RB ',85.468, 3.98E-10,5.,
     *       ' SR ',87.62 , 7.94E-10,5.,
     *       ' Y  ',88.906, 1.73E-10,5.,
     *       ' ZR ',91.224, 3.98E-10,5.,
     *       ' NB ',92.906, 2.63E-11,5.,
     *       ' MO ',95.94 , 8.31E-11,5.,
     *       ' TC ',98.906, 1.00E-12,5.,
     *       ' RU ',101.07, 6.91E-11,5.,
     *       ' RH ',102.91, 1.31E-11,5.,
     *       ' PD ',106.42, 4.89E-11,5.,
     *       ' AG ',107.87, 8.70E-12,5.,
     *       ' CD ',112.41, 7.24E-11,5.,
     *       ' IN ',114.82, 4.57E-11,5.,
     *       ' SN ',118.71, 1.00E-10,5.,
     *       ' SB ',121.76, 1.00E-11,5.,
     *       ' TE ',127.60, 1.73E-10,5.,
     *       ' I  ',126.90, 3.23E-11,5.,
     *       ' XE ',131.29, 1.69E-10,5.,
     *       ' CS ',132.91, 1.31E-11,5.,
     *       ' BA ',137.33, 1.34E-10,5.,
     *       ' LA ',138.91, 1.65E-11,5.,
     *       ' CE ',140.12, 3.54E-11,5.,
     *       ' PR ',140.91, 5.12E-12,5.,
     *       ' ND ',144.24, 3.16E-11,5.,
     *       ' PM ',146.92, 1.00E-12,5.,
     *       ' SM ',150.36, 1.00E-11,5.,
     *       ' EU ',151.96, 3.23E-12,5.,
     *       ' GD ',157.25, 1.31E-11,5.,
     *       ' TB ',158.93, 7.94E-13,5.,
     *       ' DY ',162.50, 1.25E-11,5.,
     *       ' HO ',164.93, 1.81E-12,5.,
     *       ' ER ',167.26, 8.51E-12,5.,
     *       ' TM ',168.93, 1.00E-12,5.,
     *       ' YB ',170.04, 1.20E-11,5.,
     *       ' LU ',174.97, 5.75E-12,5.,
     *       ' HF ',178.49, 7.58E-12,5.,
     *       ' TA ',180.95, 1.34E-12,5.,
     *       ' W  ',183.85, 1.28E-11,5.,
     *       ' RE ',186.21, 1.86E-12,5.,
     *       ' OS ',190.2 , 2.81E-11,5.,
     *       ' IR ',192.22, 2.23E-11,5.,
     *       ' PT ',195.08, 6.30E-11,5.,
     *       ' AU ',196.97, 1.02E-11,5.,
     *       ' HG ',200.59, 1.23E-11,5.,
     *       ' TL ',204.38, 7.94E-12,5.,
     *       ' PB ',207.20, 7.07E-11,5.,
     *       ' BI ',208.98, 5.12E-12,5.,
     *       ' PO ',209.98, 1.00E-12,5.,
     *       ' AT ',209.99, 1.00E-12,5.,
     *       ' RN ',222.02, 1.00E-12,5.,
     *       ' FR ',223.02, 1.00E-12,5.,
     *       ' RA ',226.03, 1.00E-12,5.,
     *       ' AC ',227.03, 1.00E-12,5.,
     *       ' TH ',232.04, 1.31E-12,5./
C
C     Ionization potentials for first 90 species:
CMH	after Allen's Astrophysical Quantities, 4th edition, p.36ff
      DATA XI/
C
C     Element Ionization potentials (eV) 
C              I     II      III     IV       V     VI     VII    VIII
C
     *'  H ',13.595,  0.   ,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *'  HE',24.580, 54.400,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *'  LI', 5.392, 75.619,122.451,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *'  BE', 9.322, 18.206,153.850,217.713,  0.  ,  0.  ,  0.  ,  0.  ,
     *'  B ', 8.296, 25.149, 37.920,259.298,340.22,  0.  ,  0.  ,  0.  ,
     *'  C ',11.264, 24.376, 47.864, 64.476,391.99,489.98,  0.  ,  0.  ,
     *'  N ',14.530, 29.593, 47.426, 77.450, 97.86,551.93,667.03,  0.  ,
     *'  O ',13.614, 35.108, 54.886, 77.394,113.87,138.08,739.11,871.39,
     *'  F ',17.418, 34.980, 62.646, 87.140,114.21,157.12,185.14,953.6 ,
     *'  NE',21.559, 41.070, 63.500, 97.020,126.30,157.91,207.21,239.0 ,
     *'  NA', 5.138, 47.290, 71.650, 98.880,138.37,172.09,208.44,264.16,
     *'  MG', 7.664, 15.030, 80.120,102.290,141.23,186.49,224.9 ,265.96, 
     *' AL ', 5.984, 18.823, 28.440,119.960,153.77,190.42,241.38,284.53, 
     *' SI ', 8.151, 16.350, 33.460, 45.140,166.73,205.11,246.41,303.07, 
     *' P  ',10.484, 19.720, 30.156, 51.354, 65.01,220.41,263.31,309.26,
     *' S  ',10.357, 23.400, 35.000, 47.290, 72.50, 88.03,280.99,328.8 ,
     *' CL ',12.970, 23.800, 39.900, 53.500, 67.80, 96.7 ,114.27,348.3 ,
     *' AR ',15.755, 27.620, 40.900, 59.790, 75.00, 91.3 ,124.0 ,143.46,
     *' K  ', 4.339, 31.810, 46.000, 60.900, 82.6 , 99.7 ,118.0 ,155.0 ,
     *' CA ', 6.111, 11.870, 51.210, 67.700, 84.39,109.0 ,128.0 ,147.0 ,
CMH *' SC ', 6.560, 12.890, 24.750, 73.900, 92.0 ,111.1 ,138.0 ,158.7 ,
     *' SC ', 6.560, 12.800, 24.760, 73.490, 91.7 ,110.7 ,138.0 ,158.7 ,
     *' TI ', 6.830, 13.630, 28.140, 43.240, 99.8 ,120.0 ,140.8 ,168.5 ,
     *' V  ', 6.740, 14.200, 29.700, 48.000, 65.2 ,128.9 ,151.0 ,173.7 ,
     *' CR ', 6.763, 16.490, 30.950, 49.600, 73.0 , 90.6 ,161.1 ,184.7 ,
     *' MN ', 7.432, 15.640, 33.690, 53.000, 76.0 , 97.0 ,119.24,196.46,
     *' FE ', 7.870, 16.183, 30.652, 54.800, 75.0 , 99.1 ,125.0 ,151.06,
     *' CO ', 7.860, 17.060, 33.490, 51.300, 79.5 ,102.0 ,129.0 ,157.0 ,
     *' NI ', 7.635, 18.168, 35.170, 54.900, 75.5 ,108.0 ,133.0 ,162.0 ,
     *' CU ', 7.726, 20.292, 36.830, 55.200, 79.9 ,103.0 ,139.0 ,166.0 ,
     *' ZN ', 9.394, 17.964, 39.722, 59.400, 82.6 ,108.0 ,134.0 ,174.0 ,
     *' GA ', 5.999,20.514, 30.710, 64.000, 87.0  ,116.0 ,140.0 ,170.0, 
     *' GE ', 7.900,15.935, 34.224, 45.71 , 93.5  ,112.0 ,144.0 ,174.0,
     *' AS ', 9.815,18.633, 28.351, 50.13 , 62.63 ,127.6 ,147.0 ,179.0, 
     *' SE ', 9.752,21.19 , 30.820, 42.944, 68.3  , 81.7 ,155.4 ,184.0, 
     *' BR ',11.814,21.8  , 36.    , 47.3  , 59.7  , 88.6 ,103.0 ,192.8, 
     *' KR ',14.000,24.36 , 36.95 , 52.5  , 64.7  , 78.5 ,111.0 ,126.0, 
     *' RB ', 4.177,27.285, 40.0  , 52.6  , 71.0  , 84.4 ,  97.4 ,136.0, 
     *' SR ', 5.685,11.030, 42.89 , 57.0  , 71.6  , 90.8 , 106.0,122.3, 
     *' Y  ', 6.217,12.24 , 20.52 , 61.8  , 77.0  , 93.0 , 116.0,129.0, 
     *' ZR ', 6.634,13.13 , 22.99 , 34.34 , 81.5  , 99.0 , 117.0,140.0, 
     *' NB ', 6.758,14.32 , 25.04 , 38.3  , 50.55 ,102.6 , 125.0,142.0, 
     *' MO ', 7.092,16.16 , 27.13 , 46.4  , 61.2  , 68.0 , 126.8,153.0, 
     *' TC ', 7.28 ,15.26 , 29.54 , 46.0  , 55.0  , 80.0 ,120. ,145.,!!last2values 
     *' RU ', 7.361,16.76 , 28.47 , 50.0  , 60.0  , 92.0 ,116.0,129.0,!!last2values 
     *' RH ', 7.459,18.08 , 31.06 , 48.0  , 65.0  , 97.0 ,117.0,140.0,!!last2values 
     *' PD ', 8.337,18.43 , 32.93 , 53.0  , 62.0  , 90.0 ,110.0,130.0, 
     *' AG ', 7.576,21.49 , 34.83 , 56.0  , 68.0  , 89.0 ,115.0,140.0, 
     *' CD ', 8.994,16.908, 37.48 , 59.0  , 72.0  , 94.0 ,115.0,145.0, 
     *' IN ', 5.786,18.870, 28.03 , 54.4  , 77.0  , 98.0 ,120.0,145.0, 
     *' SN ', 7.343,14.632, 30.503, 40.734, 72.28, 103.0,125.0,150.0, 
     *' SB ', 8.64 ,16.531, 25.3  , 44.2   , 56.0 , 108.0,130.0,155.0, 
     *' TE ', 9.010,18.6  , 27.96 , 37.41  , 58.75, 70.7, 137.0,165.0, 
     *' I  ',10.451,19.131, 33.0  , 42.0   , 66.0 , 81.0, 100.0,170.0, 
     *' XE ',12.130,21.21 , 32.123, 46.0   , 57.0 , 82.0, 100.0,120.0, 
     *' CS ', 3.894,23.157, 35.0  , 46.0   , 62.0 , 74.0, 100.0,120.0, 
     *' BA ', 5.212,10.004, 20.   , 49.0   , 62.0 , 80.0,  95.0,120.0, !!3rd value
     *' LA ', 5.577,11.06 , 19.177, 52.0   , 66.0 , 80.0, 100.0,115.0, 
     *' CE ', 5.538,10.85 , 20.198, 36.72  , 70.0 , 85.0, 100.0,120.0, 
     *' PR ', 5.464,10.55 , 21.624, 38.95  , 57.45, 89.0, 105.0,120.0, 
     *' ND ', 5.525,10.73 , 20.    , 36.    , 70. , 85., 110. , 130.,!3.,4.,5.,6.
     *' PM ', 5.55 ,10.90 , 20.    , 36.    , 70. , 85., 110. ,135., !3.,4.,5.,6.,7
     *' SM ', 5.644,11.07 , 20.    , 36.    , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' EU ', 5.670,11.241, 20.    , 36.    , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' GD ', 6.150,12.09 , 20.    , 36.    , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' TB ', 5.864,11.52, 20.    , 36.  , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' DY ', 5.939,11.67, 20.    , 36.  , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' HO ', 6.022,11.80, 20.    , 36.  , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' ER ', 6.108,11.93, 20.    , 36.  , 70. , 85., 110. ,135., !3.,4.,5.,6,7,8
     *' TM ', 6.184,12.05, 23.68  , 36.  , 70. , 85., 110. ,135., !4.,5.,6,7,8
     *' YB ', 6.254,12.176,25.05  , 36.  , 70. , 85., 110. ,135., !4.,5.,6,7,8 
     *' LU ', 5.425,13.9  , 20.96 , 36.  , 70. , 85., 110. ,135., !4.,5.,6,7,8 
     *' HF ', 6.825,14.9  , 23.3  , 33.3 , 70. , 85., 110. ,135., !5.,6,7,8
     *' TA ', 7.89 ,16.0  , 22.0  , 33.0 , 45.0,85., 110. ,135., !6,7,8
     *' W  ', 7.98 ,18.0  , 24.0  , 35.0 , 48.0,61.0, 110. ,135., !7,8
     *' RE ', 7.88 ,17.0  , 26.0  , 38.0 , 51.0,64.0, 79.0 ,135., !8
     *' OS ', 8.7  ,17.0  , 25.0  , 40.0 , 54.0,68.0, 83.0 ,100., 
     *' IR ', 9.1  ,17.0  , 27.0  , 39.0 , 57.0,72.0, 88.0,105.0, 
     *' PT ', 9.0  ,18.563, 28.0  , 41.0 , 55.0,75.0 , 92.0,110.0, 
     *' AU ', 9.225,20.5  , 30.0  , 44.0 , 58.0,75.0 , 96.0,115.0, 
     *' HG ',10.437,18.756, 34.2 , 46.0 , 61.0, 77.0, 94.0,120.0, 
     *' TL ', 6.108,20.428, 29.83, 50.7 , 64.0, 81.0, 98.0,115.0, 
     *' PB ', 7.417,15.032, 31.937,42.32, 68.8, 84.0,103.0,120.0, 
     *' BI ', 7.289,16.69 , 25.56 ,45.3 , 56.0, 88.3 ,107.0,125.0, 
     *' PO ', 8.417,19.0  , 27.0   ,38.0 , 61.0,73.0  ,112.0,130.0, 
     *' AT ', 9.3  ,20.0   , 29.0   ,41.0 , 51.0,78.0 ,91.0  ,140.0, 
     *' RN ',10.748,21.0    , 29.0  ,44.0 , 55.0,67.0 ,97.0  ,110.0, 
     *' FR ', 4.0  ,22.0   , 33.0   ,43.0 , 59.0,71.0 ,84.0 ,115.0, 
     *' RA ', 5.273,10.147, 34.0   ,46.0 , 58.0,76.0 ,89.0 ,105.0, 
     *' AC ', 5.17 ,12.1   , 20.0   ,49.0 , 62.0,76.0 ,95.0 ,110.0, 
     *' TH ', 6.08 ,11.5   , 20.0   ,28.8 , 65.0,80.0 ,94.0,115.0 /
C
      !DATA XIFE /8*0.,233.6,262.1/
      !DATA NTOTA /30/
      PARAMETER (enh1=13.595, enhe1=24.580, enhe2=54.400)
      LY=YTOT.EQ.0.
      LW=WMY.EQ.0.
C
C     An element (hydrogen through zinc) can be considered in one of
C     the three following options:
C     1. explicitly - some of energy levels of some of its ionization
C                     states are considered explicitly, ie. their
C                     populations follow from the input model
C                     atmosphere
C     2. normally   - ie. the ionization equilibium and level populations
C                     are determined from the Saha-Boltzmann equations;
C                     these elements are allowed to contribute to the
C                     opacity only through the line opacity
C     3. not considered at all
C
C     Input:
C
C     For each element from 1 (hydrogen) to 30 (zinc) the following
C     parameters:
C
C     MA     =  0  - if the element is not considered (option 3)
C            =  1  - if the element is non-explicit (option 2)
C            =  2  - if the element is explicit (option 1)
C     NA0,NAK - have the meaning only for MA=2; indicate that the
C               explicit energy levels of the present species have
C               the indices between NA0 and NAK (NAK is thus the index
C               of the highest ionization state, which is represented
C               as one-level ion).
C     ION    -  has the meaning for MA=1 only;
C               if ION=0, standard number of ionization degrees is
C                         considered
C               if ION>0, then ION ionization degrees is considered
C     MODPF  -  mode of evaluation of partition functions
C            =  0  -  standard evaluation (see procedure PARTF)
C            >  0  -  partition functions assumed constant (independent
C                     of temperature and electron density), and given
C                     as input parameters PFS(J), for J=1 - neutrals,
C                     J=2 - once ionized, etc., up to (possibly) J=5;
C            <  0  -  non-standard evaluation, by user supplied
C                     procedure PFSPEC
C     ABN    -  if ABN=0, solar abundance is assumed (given above;
C                         abundance here is always assumed as relative
C                         to hydrogen by number
C               if ABN>0, non-solar abundance ABN is assumed
C               if ABN<0, non-solar abundance is assumed, abs(ABN) means
C                         now the abundance expressed as a fraction of the
C                         solar abundance
C     PFS    -  see above
C
c      iatref=1
      READ(555,*) NATOMS
cc    READ(81,*) NATOMS
      WRITE(6,600)
      IAT=0
      IREF=0
      DO 10 I=1,MATOM
         DO 10 J=1,MION
   10       RR(I,J)=0.
      DO 20 I=1,MATOM
         TYPAT(I)=D(1,I)
         LGR(I)=.TRUE.
         LRM(I)=.TRUE.
         IATEX(I)=-1
 	 IF(I.LE.NATOMS) THEN
c***   READING FILE 555 input_sun (input5)
            READ(555,*) MA,NA0,NAK,ION,MODPF(I),ABN,(PFSTD(J,I),J=1,5)
cc          READ(81,*) MA,NA0,NAK,ION,MODPF(I),ABN,(PFSTD(J,I),J=1,5)
 	  ELSE IF(IMODE.LE.1) THEN
 	    MA=1
 	    ABN=0.
 	    ION=0
 	    MODPF(I)=0
           ELSE
 	    MA=0
 	 END IF
         IF(MA.EQ.0) GO TO 20
         AMAS(I)=D(2,I)
         ABND(I)=D(3,I)
         if(iref.gt.0) abnd(i)=d(3,i)*abnd(iref)/d(3,iref)
         IONIZ(I)=D(4,I)
C
C        increase the standard highest ionization for Teff larger
C        than 50000 K for N, O, Ne, and Fe
C
         IF(TEFF0.GT.5.D4) THEN
            IF(I.EQ.7) IONIZ(I)=6
            IF(I.EQ.8) IONIZ(I)=7
c           IF(I.EQ.10) IONIZ(I)=9
            IF(I.EQ.26) IONIZ(I)=9
         END IF
C
         DO J=1,9
            IF(J.LE.8) ENEV(I,J)=XI(J+1,I)
            if(enev(i,j).ge.enhe2) then
               inpot(i,j)=3
             else if(enev(i,j).ge.enhe1) then
               inpot(i,j)=2
             else
               inpot(i,j)=1
            end if
         END DO
         LGR(I)=.FALSE.
         IF(ABN.GT.0) ABND(I)=ABN
         IF(ABN.LT.0) ABND(I)=-ABN*D(3,I)
         IF(ION.NE.0) IONIZ(I)=ION
         IF(MA.EQ.1) THEN
            LRM(I)=.FALSE.
            IATEX(I)=0
          ELSE
            IAT=IAT+1
            IATEX(I)=IAT
            IF(IAT.EQ.IATREF) THEN
               IREF=I
               ABNREF=ABND(I)
            END IF
C
C           store parameters for explicit atoms
C
            ABUN (IAT)=ABND(I)
            AMASS(IAT)=AMAS(I)*HMASS
            NUMAT(IAT)=I
            N0A(IAT)=NA0
            NKA(IAT)=NAK
         END IF
         IF(LY) YTOT=YTOT+ABND(I)
         IF(LW) WMY=WMY+ABND(I)*AMAS(I)
         ABN=ABND(I)/D(3,I)
         IF(MA.EQ.1) WRITE(6,601) I,TYPAT(I),ABND(I),ABN
         IF(MA.EQ.2) WRITE(6,602) I,TYPAT(I),ABND(I),ABN,IAT,NA0,NAK
   20 CONTINUE

      IF(IMODE.LE.1) NATOMS=MATOM
      WMM=WMY*HMASS/YTOT
      DO 30 JJ=1,MLEV
   30    RELAB(JJ)=1.
      IF(ICHEMC.NE.1) go to 100
C
C     abundance change with respect to the model atmosphere input
C     (unit 5);
C     this option is switched on by the parameter ICHEMC (read from
C     unit 55), if it is non-zero, an additional input from
C     unit 56 is required
C
C     unit 56 input:
C
C     NCHANG  -  number of chemical elements for which the abundances
C                are going to be changes;
C
C     then there are NCHANG records, each contains:
C
C     I       - atomic number
C     ABN     - new abundance; coded using teh same conventions as in
C               the standard input
C
      READ(56,*,ERR=566,END=566) NCHANG
      WRITE(6,610)
      DO 60 II=1,NCHANG
         READ(56,*) I,ABN
         ABND(I)=D(3,I)
C***  CHANGES BY MARGIT HABERREITER
CMH         IF(ABN.GT.0) ABND(I)=ABN
		IF(ABN.GT.0) then
			ABND(I)=ABN
c		print *,i,abnd(i),'abundance changed!!!'
c		write(6,*) i,IATEX(I),D(3,I),'changed to',abnd(i)
		end if
         IF(ABN.LT.0) ABND(I)=-ABN*D(3,I)
c		print *,i,abnd(i),'abundance changed!!!'
c		write(6,*) i,abnd(i),'abundance changed!!!'
         IF(IATEX(I).GT.0) THEN
            RELA=ABND(I)/ABUN(IATEX(I))
		print *,i,IATEX(I),'RELA = ', RELA
            ABUN(IATEX(I))=ABND(I)
           DO 40 JJ=N0A(IATEX(I)),NKA(IATEX(I))
   40          RELAB(JJ)=RELA
         END IF
         ABN=ABND(I)/D(3,I)
         WRITE(6,601) I,TYPAT(I),ABND(I),ABN
c	pause
   60 CONTINUE
C
C     renormalize abundances to have the standard element abundance
C     equal to unity
C 
  100 IF(IREF.LE.1) RETURN
      ytot=0.
      wmy=0.
      write(6,620)
      DO 70 I=1,NATOMS
         IAT=IATEX(I)
         IF(IAT.LT.0) GO TO 70
         ABND(I)=ABND(I)/ABNREF
         YTOT=YTOT+ABND(I)
         WMY=WMY+ABND(I)*AMAS(I)
         ABNR=ABND(I)/D(3,I)
         IF(IAT.EQ.0) THEN
            WRITE(6,601) I,TYPAT(I),ABND(I),ABNR
          ELSE
            ABUN(IAT)=ABND(I) 
            WRITE(6,602) I,TYPAT(I),ABND(I),ABNR,IAT,N0A(IAT),NKA(IAT)
         END IF
   70 CONTINUE
      WMM=WMY*HMASS/YTOT
      RETURN
C
  566 WRITE(6,656)
      STOP
C
  600 FORMAT(//'    CHEMICAL ELEMENTS INCLUDED'/
     *         '    --------------------------'//
     * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
     * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)
  601 FORMAT(1H ,I4,3X,A5,1P2E14.2)
  602 FORMAT(1H ,I4,3X,A5,1P2E14.2,3X,
     *       'EXPLICIT: IAT=',I3,'  N0A=',I3,'  NKA=',I3)
  610 FORMAT(//'    CHEMICAL ELEMENTS INCLUDED - CHANGED (unit 56)'
     *           /'    --------------------------'//
     * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
     * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/
CMH
     *	'ATTENTION EXPLIXIT LEVELS:'/ 
     *	'ABUNDANCE(ELEMENT)=RELA * SUM(POPNUM(ELEMENT))!!!!!')
  620 FORMAT(1H0//'    CHEMICAL ELEMENTS INCLUDED - RENORMALIZATION'/
     *            '    --------------------------'//
     * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
     * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)
  656 FORMAT(//' CHEMICAL COMPOSITION COULD NOT BE READ FROM '
     *       'UNIT 56',//,' STOP.')
      END subroutine
C
C
C ********************************************************************
C 
C
      !*** called by INIMOD
      SUBROUTINE STATE(ID)
      use MOD_PARTF2
C
C     modified LTE Saha eautions - using radiation temperatures
C     (Schaerer and Schmutz AA 288, 321, 1994)
C
c      INCLUDE 'PARAMS.FOR'
c      INCLUDE 'MODELP.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: ID
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      dimension FFI(MION0)
	!REAL*8 SUM2,UM3
C
      TE=TEMP(ID)
      ANE=ELEC(ID)
CMH-orig      
CMH  NEW: !Loop only over maximal number of nlte levels defined in PARAMS.FOR
	DO 50 I=2,NATOMS
         IF(LGR(I)) GO TO 50
         ION=IONIZ(I)
         RS=1.
         T=TRAD(INPOT(I,1),ID)
!         IF (TE .NE. T) THEN
!         print *, 'XM...', i, TE, T
!         stop
!         ENDIF
!         print *, 'mytest1', i, TE, T
         X=SQRT(T/ANE)
         XMX=2.145E4*SQRT(X)
!         CALL  PARTF(I,1,T,ANE,XMX,UM1)
!         CALL  PARTF2(I,1,T,UM2)
!         print*, 'mytest2', i, UM1, UM2, UM1/UM2
       
CMH	I: number of atom
CMH	J=1: ground state
CMH	ANE: electron density
CMH	UM: returned partition function
         IF (I .LE. 30) THEN
          CALL PARTF(I,1,T,ANE,XMX,UM)
         ELSE 
          CALL PARTF2(I,1,T,UM)
         ENDIF
         PFSTD(1,I)=UM
         JMAX=1
         DO J=2,ION
            T=TRAD(INPOT(I,J),ID)
            TLN=LOG(T)*1.5
            TK=BOLK*T
            THL=11605./T
            X=SQRT(T/ANE)
            XMX=2.145E4*SQRT(X)
            DCH=EH/XMX/XMX/TK
            DCHT=DCH*(J-1)
            FI=36.113+TLN-THL*ENEV(I,J-1)+DCHT
            X=J
            XMAX=XMX*SQRT(X)
CMH	FOR HEAVY ELEMENTS ABOVE AN=30 USE PARTITION FUNCTION
CMH 	ACCORDING TO A.IRWIN, 1981
		IF (I .LE. 30) THEN
	            CALL PARTF(I,J,T,ANE,XMAX,U)
		ELSE 
		    CALL PARTF2(I,J,T,U)
		ENDIF
            PFSTD(J,I)=U
            FI=EXP(FI)*U/UM/ANE*wdil(id)*sqrt(te/t)
            FFI(J)=FI
            IF(FFI(J).GT.1.) JMAX=J
            UM=U

         ENDDO
         IF(JMAX.LT.ION) THEN
         RA=1.
         DO 20 J=JMAX+1,ION
            RA=RA*FFI(J)
            RR(I,J)=RA/PFSTD(J,I)
            RS=RS+RA
   20    CONTINUE
         END IF
         IF(JMAX.GT.1) THEN
         RA=1.
         DO 30 JJ=1,JMAX-1
            J=JMAX-JJ
            RA=RA/FFI(J+1)
            RR(I,J)=RA/PFSTD(J,I)
            RS=RS+RA
   30    CONTINUE
         END IF
         RR(I,JMAX)=ABND(I)/RS
         DO 40 J=1,ION
            IF(J.NE.JMAX) RR(I,J)=RR(I,J)*RR(I,JMAX)
            if(rr(i,j).lt.1.e-35) rr(i,j)=0.
   40    CONTINUE
         RR(I,JMAX)=RR(I,JMAX)/PFSTD(JMAX,I)
   50 CONTINUE
      RETURN
      END SUBROUTINE
C
C
C ********************************************************************
C 
C
      SUBROUTINE PARTF(IAT,IZI,T,ANE,XMAXN,U)
C
C   PARTITION FUNCTIONS AFTER  TRAVING, BASCHEK, AND HOLWEGER, 1966,
C   ABHAND. HAMBURG. STERNWARTE, BAND VIII, NR. 1
C
c      INCLUDE 'PARAMS.FOR'
      implicit real*8(a-h,o-z)
      integer,intent(in   ) :: IAT,IZI
      real*8, intent(in   ) :: T,ANE,XMAXN
      real*8, intent(inout) :: U
	INCLUDE '../inc/PARAMS.FOR'
      REAL*4 AHH( 6),  ALB(12),  AB (11),  AC (19),  AN (30),  AO (49),
     *       AF (34),  ANN(23),  ANA(19),  AMG(15),  AAL(17),  ASI(23),
     *       AP (19),  AS (29),  ACL(28),  AAR(25),  AK (30),  ACA(17),
     *       ASC(24),  ATI(33),  AV (33),  ACR(29),  AMN(28),  AFE(35),
     *       ACO(29),  ANI(23),  ACU(20),  AZN(18)
      REAL*4 GHH( 6),  GLB(12),  GB (11),  GC (19),  GN (30),  GO (49),
     *       GF (34),  GNN(23),  GNA(19),  GMG(15),  GAL(17),  GSI(23),
     *       GP (19),  GS (29),  GCL(28),  GAR(25),  GK (30),  GCA(17),
     *       GSC(24),  GTI(33),  GV (33),  GCR(29),  GMN(28),  GFE(35),
     *       GCO(29),  GNI(23),  GCU(20),  GZN(18)
      REAL*4 XL1(99), XL2(123),  XL(222),
     *       CH1(66),  CH2(72),  CH3(55),  CH4(29),  CHION(222)
      REAL*4 ALF(678), GAM(678), MAX,      MAX2   
      INTEGER*2 II1(5,15),II2(5,15),INDEX0(5,30),
     *          IS1(53),IS2(70),IS(123),INDEXS(123),
     *          IM1(99),IM2(123),IM(222),INDEXM(222),
     *          IGP1(99),IGP2(123),IGPR(222),
     *          IG01(53),IG02(70),IG0(123)
C
      DATA II1      /   1,  -1,   0,   0,   0,
     *                  2,   3,  -1,   0,   0,
     *                  4,   5,  -2,  -1,   0,
     *                  6,   7,  -1,  -2,  -1,
     *                  8,   9,  10,  -1,  -2,
     *                 11,  12,  13,  14,  -1,
     *                 15,  16,  17,  18,  19,
     *                 20,  21,  22,  23,  24,
     *                 25,  26,  27,  28,  -6,
     *                 29,  30,  31,  32,  -9,
     *                 33,  34,  35,  36,  -4,
     *                 37,  38,  39,  40,  -9,
     *                 41,  42,  43,  44,  -6,
     *                 45,  46,  47,  48,  -1,
     *                 49,  50,  51,  52,  53                         /
      DATA II2      /  54,  55,  56,  57,  58,
     *                 59,  60,  61,  62,  63,
     *                 64,  65,  66,  67,  68,
     *                 69,  70,  71,  72,  73,
     *                 74,  75,  76,  77,  -9,
     *                 78,  76,  80,  81,  82,
     *                 83,  84,  85,  86,  87,
     *                 88,  89,  90,  91,  92,
     *                 93,  94,  95,  96,  97,
     *                 98,  99, 100, 101, 102,
     *                103, 104, 105, 106, 107,
     *                108, 109, 110, 111, -25,
     *                112, 113, 114, 115,  -1,
     *                116, 117, 118, 119,  -1,
     *                120, 121, 122, 123,  -1                         /
C
      DATA IG01     /   2,
     *                  1,   2,
     *                  2,   1,
     *                  1,   2,
     *                  2,   1,   2,
     *                  1,   2,   1,   2,
     *                  4,   1,   2,   1,   2,
     *                  5,   4,   1,   2,   1,
     *                  4,   5,   4,   1,
     *                  1,   4,   5,   4,
     *                  2,   1,   4,   5,
     *                  1,   2,   1,   4,
     *                  2,   1,   2,   1,
     *                  1,   2,   1,   2,
     *                  4,   1,   2,   1,   2                         /
      DATA  IG02    /   5,   4,   1,   2,   1,
     *                  4,   5,   4,   1,   2,
     *                  1,   4,   5,   4,   1,
     *                  2,   1,   4,   5,   4,
     *                  1,   2,   1,   4,
     *                  4,   3,   4,   1,   4,
     *                  5,   4,   5,   4,   1,
     *                  4,   1,   4,   5,   4,
     *                  7,   6,   1,   4,   5,
     *                  6,   7,   6,   1,   4,
     *                  9,  10,   9,   6,   1,
     *                 10,   9,  10,  20,
     *                  9,   6,   9,  28,
     *                  2,   1,   6,  21,
     *                  1,   2,   1,  10                              /
C
      DATA IS1      /   1,
     *                  1,   1,
     *                  1,   1,
     *                  2,   1,
     *                  1,   2,   1,
     *                  1,   2,   2,   1,
     *                  2,   2,   3,   2,   1,
     *                  3,   4,   3,   5,   2,
     *                  2,   3,   4,   3,
     *                  2,   2,   3,   2,
     *                  1,   2,   2,   3,
     *                  1,   1,   2,   2,
     *                  2,   2,   1,   2,
     *                  1,   2,   2,   1,
     *                  2,   1,   1,   1,   1                         /
      DATA  IS2     /   3,   2,   1,   2,   2,
     *                  2,   3,   2,   1,   1,
     *                  2,   2,   3,   1,   1,
     *                  1,   2,   3,   3,   2,
     *                  2,   1,   2,   2,
     *                  3,   1,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  2,   3,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  3,   2,   2,   1,   1,
     *                  4,   2,   1,   1,
     *                  2,   2,   1,   1,
     *                  3,   2,   1,   1,
     *                  3,   3,   1,   1                              /
C
      DATA IM1      /   2,
     *                  2,   2,
     *                  2,   2,
     *                  3,   2,   3,
     *                  3,   3,   2,   3,
     *                  4,   3,   3,   3,   3,   3,
     *                  3,   3,   4,   3,   3,   4,   2,   3,   2,   3,
     *                  4,   2,   2,   4,   2,   3,   3,   4,   4,   2,
     *                  3,   4,   2,   2,   2,   3,   3,
     *                  3,   3,   4,   2,   2,
     *                  4,   2,   3,   2,   5,   2,   2,
     *                  2,   2,   3,   2,   4,   2,   2,   4,   2,
     *                  2,   2,   2,   3,   2,   4,   2,   2,
     *                  3,   3,   2,   2,   3,   2,
     *                  3,   2,   3,   2,   3,   2,   2,
     *                  5,   4,   4,   4,   3,   3,
     *                  3,   2,   4,   4,   3,   3                    /
      DATA  IM2     /   4,   2,   2,   4,   2,   5,   4,   2,   3,   1,
     *                  3,   2,   5,   2,   2,   4,   2,   4,   4,
     *                  2,   2,   3,   2,   4,   2,   2,   4,   4,
     *                  3,   2,   3,   3,   2,   3,
     *                  4,   2,   2,   4,   2,
     *                  3,   2,   3,   2,   2,   3,   2,
     *                  4,   3,   3,   5,   4,   2,   3,
     *                  6,   4,   3,   6,   3,   5,   4,   2,
     *                  5,   3,   5,   4,   4,   4,   4,   4,
     *                  3,   3,   3,   4,   4,   4,   4,   4,
     *                  3,   2,   3,   4,   4,   4,   4,   4,
     *                  4,   4,   3,   5,   3,   4,   4,   4,   4,
     *                  5,   3,   3,   3,   5,   4,   5,   1,
     *                  6,   3,   5,   3,   5,   1,
     *                  2,   3,   3,   4,   3,   4,   1,
     *                  2,   2,   2,   3,   3,   2,   3,   1          /
C
      DATA IGP1     /   2,
     *                  4,   2,
     *                  2,   4,
     *                  4,  12,   2,
     *                  2,   4,  12,   2,
     *                 12,   2,  18,   4,  12,   2,
     *                 18,  10,  12,  24,   2,  18,   6,   4,  12,   2,
     *                  8,  20,  12,  18,  10,   2,  10,  12,  24,  20,
     *                  2,  18,   6,  18,  10,   4,  12,
     *                 18,  10,   8,  20,  12,
     *                 18,  10,   2,  10,  12,  24,  20,
     *                  8,   4,  18,  10,   8,  20,  12,  18,  10,
     *                  2,   8,   4,  18,  10,   8,  20,  12,
     *                  4,   2,   8,   4,  18,  10,
     *                  2,  18,   4,  12,   2,   8,   4,
     *                 12,   2,  18,   4,  12,   2,
     *                 18,  10,  12,   2,   4,   2                    /
      DATA  IGP2    /   8,  20,  12,  18,  10,  12,   2,  18,   4,  12,
     *                 18,  10,   8,  20,  12,  18,  10,  12,   2,
     *                  8,   4,  18,  10,   8,  20,  12,  18,  12,
     *                  2,   8,   4,  18,  10,   2,
     *                  8,  20,  12,  18,  10,
     *                  4,  20,   2,   8,   4,  18,  10,
     *                 30,  42,  18,  20,   2,  12,  18,
     *                 56,  56,  28,  42,  10,  20,   2,  12,
     *                 50,  70,  56,  72,  64,  42,  20,   2,
     *                 12,  60,  40,  50,  18,  56,  42,  20,
     *                 14,  10,  50,  12,  72,  50,  56,  42,
     *                 60,  56,  40,  50,  18,  12,  72,  50,  56,
     *                 42,  70,  42,  18,  56,  24,  50,  12,
     *                 20,  56,  42,  18,  56,  50,
     *                  2,  30,  10,  20,  56,  42,  56,
     *                  4,   8,  12,   2,  30,  10,  20,  42          /
C
      DATA XL1      /11.0,
     *                8.0,12.0,
     *                6.0, 6.0,
     *                6.0, 4.0, 8.0,
     *                9.0, 6.0, 4.0, 6.0,
     *                6.0, 6.0, 5.0, 6.1, 5.0, 6.0,
     *                6.1, 4.0, 5.0, 3.9, 6.0, 5.0, 4.0, 6.0, 6.3, 6.0,
     *                8.0, 6.0, 3.4, 6.0, 5.0, 3.9, 3.9, 6.0, 4.9, 4.0,
     *                5.9, 5.0, 4.9, 4.0, 4.0, 6.0, 6.0,
     *                4.0, 4.0, 5.0, 4.0, 4.0,
     *                5.0, 4.0, 3.9, 4.0, 5.0, 5.0, 4.0,
     *                6.0, 6.0, 5.0, 4.0, 3.9, 4.0, 4.0, 5.0, 5.0,
     *                7.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
     *                7.0, 7.0, 5.0, 5.0, 5.0, 5.0,
     *                7.0, 4.0, 7.0, 4.0, 7.0, 5.0, 5.0,
     *                6.1, 5.9, 5.0, 5.0, 5.0, 7.0,
     *                5.0, 5.0, 5.0, 7.0, 8.6, 8.0                    /
      DATA  XL2     / 6.0, 5.0, 5.0, 5.0, 5.0, 3.5, 5.0,14.4, 5.0, 4.0,
     *                6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,
     *                6.0, 6.0, 5.1, 5.0, 5.0, 5.0, 5.0, 5.0, 4.0,
     *                7.0, 5.0, 5.0, 6.0, 6.0, 5.0,
     *                6.0, 5.0, 5.0, 3.6, 4.0,
     *                5.9, 6.0, 7.0, 5.0, 4.9, 5.0, 4.3,
     *                4.9, 4.9, 5.0, 5.0, 6.0, 4.6, 3.8,
     *                5.0, 4.7, 5.0, 5.0, 5.0, 5.0, 6.0, 4.8,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,11.2,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,
     *                6.0, 5.0, 6.0, 7.0, 5.0, 5.0, 5.0, 5.0,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 5.0, 3.6, 3.8,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.0,
     *                5.4, 5.0, 9.0, 5.0, 5.0, 3.0,
     *                8.0, 6.0, 5.0, 7.0, 5.0, 5.0, 2.9,
     *                8.0, 5.0, 5.0, 8.0, 5.0, 5.0, 5.0, 2.8          /
C
C
      DATA CH1      /  13.595 ,
     *                 24.580 ,  54.403 ,
     *                  5.390 , 75.619 ,
     *                  9.320 ,  13.278 ,  18.206 ,
     *                  8.296 ,  25.149 ,  31.146 ,  37.920 ,
     *                 11.256 ,  24.376 ,  30.868 ,  47.871 ,  55.873 ,
     *                 64.476 ,
     *                 14.529 ,  16.428 ,  29.593 ,  36.693 ,
     *                 47.426 ,  55.765 ,  63.626 ,  77.450 ,  87.445 ,
     *                 97.863 ,
     *                 13.614 ,  16.938 ,  18.630 ,
     *                 35.108 ,  37.621 ,  40.461 ,  42.584 ,
     *                 54.886 ,  63.733 ,  70.556 ,
     *                 77.394 ,  87.609 ,  97.077 , 103.911 , 106.116 ,
     *                113.873 , 125.863 ,
     *                 17.418 ,  20.009 ,  34.977 ,  39.204 ,  41.368 ,
     *                 62.646 ,  65.774 ,  69.282 ,  71.882 ,
     *                 87.139 ,  97.852 , 106.089 ,
     *                 21.559 ,  21.656 ,  41.071 ,  44.274 ,
     *                 63.729 ,  68.806 ,  71.434 ,  97.162 , 100.917 /
      DATA CH2      /   5.138 ,  47.290 ,  47.459 ,  71.647 ,  75.504 ,
     *                 98.880 , 104.778 , 107.864 ,
     *                  7.644 ,  15.031 ,  80.117 ,  80.393 ,
     *                109.294 , 113.799 ,
     *                  5.984 ,  10.634 ,  18.823 ,  25.496 ,
     *                 28.441 , 119.957 , 120.383 ,
     *                  8.149 ,  16.339 ,  22.894 ,
     *                 33.459 ,  42.333 ,  45.130 ,
     *                 10.474 ,  11.585 ,  19.720 ,
     *                 30.156 ,  51.354 ,  65.007 ,
     *                 10.357 ,  12.200 ,  13.401 ,  23.405 ,  24.807 ,
     *                 35.047 ,  47.292 ,  57.681 ,  72.474 ,  85.701 ,
     *                 13.014 ,  14.458 ,  23.798 ,  26.041 ,  27.501 ,
     *                 39.904 ,  41.610 ,  53.450 ,  67.801 ,
     *                 15.755 ,  15.933 ,  27.619 ,  29.355 ,
     *                 40.899 ,  42.407 ,  45.234 ,  59.793 ,  75.002 ,
     *                  4.339 ,  31.810 ,  32.079 ,
     *                 45.738 ,  47.768 ,  50.515 ,
     *                 60.897 ,  63.890 ,  65.849 ,  82.799 ,  85.150 /
      DATA CH3      /   6.111 ,   7.808 ,  11.868 ,
     *                 51.207 ,  51.596 ,  67.181 ,  69.536 ,
     *                  6.538 ,   7.147 ,   8.042 ,
     *                 12.891 ,  24.752 ,  74.090 ,  91.847 ,
     *                  6.818 ,  6.953 ,  7.411 ,
     *                 13.635 ,  14.685 ,  28.137 ,  43.236 , 100.083 ,
     *                  6.738 ,   7.101 ,  14.205 ,  15.670 ,  16.277 ,
     *                 29.748 ,  48.464 ,  65.198 ,
     *                  6.763 ,   8.285 ,   9.221 ,
     *                 16.493 ,  18.662 ,  30.950 ,  49.580 ,  73.093 ,
     *                  7.432 ,   8.606 ,   9.240 ,  15.636 ,  18.963 ,
     *                 33.690 ,  53.001 ,  76.006 ,
     *                  7.896 ,   8.195 ,   8.927 ,  16.178 ,  18.662 ,
     *                 30.640 ,  34.607 ,  56.001 ,  79.001           /
      DATA CH4      /   7.863 ,   8.378 ,   9.160 ,   9.519 ,
     *                 17.052 ,  18.958 ,  33.491 ,  53.001 ,
     *                  7.633 ,   8.793 ,  18.147 ,  20.233 ,  35.161 ,
     *                 56.025 ,
     *                  7.724 ,  10.532 ,  10.980 ,
     *                 20.286 ,  27.985 ,  36.826 ,  61.975 ,
     *                  9.391 ,  17.503 ,  17.166 ,
     *                 17.959 ,  27.757 ,  28.310 ,  39.701 ,  65.074 /
C
      DATA AHH      /  20.4976, 747.5023,
     *                 28.1703, 527.8296,  22.2809, 987.7189          /
      DATA GHH      /  10.853 ,  13.342 ,
     *                 21.170 ,  24.125 ,  43.708 ,  53.542           /
C
      DATA ALB      /   8.4915,  97.5015,  23.3299, 192.6701,
     *                  9.1849,  32.9263, 183.8887,  19.9563,  88.0437,
     *                  6.0478,  35.9723, 233.9798                    /
      DATA GLB      /   2.022 ,   4.604 ,  62.032 ,  72.624 ,
     *                  2.735 ,   6.774 ,   8.569 ,  10.750 ,  11.672 ,
     *                  3.967 ,  12.758 ,  16.692                     /
C
      DATA AB       /   4.0086,  19.6741, 402.3110,
     *                  9.7257,  30.9262, 186.3466,  44.1629,  60.8371,
     *                  6.0084,  23.5767,  76.4149                    /
      DATA GB       /   0.002 ,   3.971 ,   7.882 ,
     *                  4.720 ,  13.477 ,  22.103 ,  23.056 ,  24.734 ,
     *                  6.000 ,  24.540 ,  32.300                     /
C
      DATA AC       /   8.0158,   5.8833,  33.7521, 595.3432,
     *                  4.0003,  17.0841,  82.9154,
     *                 15.9808,  48.2044, 435.8093,
     *                 10.0281,  15.7574, 186.2109,
     *                 15.4127,  55.9559, 243.6311,
     *                  6.0057,  23.5757,  76.4185                    /
      DATA GC       /   0.004 ,   1.359 ,   6.454 ,  10.376 ,
     *                  0.008 ,  16.546 ,  21.614 ,
     *                  5.688 ,  15.801 ,  26.269 ,
     *                  6.691 ,  25.034 ,  40.975 ,
     *                 17.604 ,  36.180 ,  47.133 ,
     *                  8.005 ,  40.804 ,  54.492                     /
C
      DATA AN       /  14.0499,  30.8008, 883.1443,
     *                 10.0000,  16.0000,  64.0000,
     *                  8.0462,   6.2669,  17.8696, 282.8084,
     *                  7.3751,  33.1390, 215.4829,
     *                  4.0003,  19.3533,  80.6462,
     *                 13.0998,  19.6425,  94.3035, 370.9539,
     *                 16.0000,  38.0000,
     *                 10.3289,  14.5021, 187.1624, 108.1615, 191.8383,
     *                  6.0044,  23.5612,  76.4344                    /
      DATA GN       /   2.554 ,   9.169 ,  13.651 ,
     *                 12.353 ,  13.784 ,  14.874 ,
     *                  0.014 ,   2.131 ,  15.745 ,  24.949 ,
     *                  6.376 ,  14.246 ,  29.465 ,
     *                  0.022 ,  31.259 ,  41.428 ,
     *                  7.212 ,  15.228 ,  34.387 ,  46.708 ,
     *                 46.475 ,  49.468 ,
     *                  8.693 ,  37.650 ,  65.479 ,  61.155 ,  79.196 ,
     *                  9.999 ,  60.991 ,  82.262                     /
C
      DATA AO       /   4.0029,   5.3656,  36.2853,1044.3447,
     *                131.0217, 868.9779,  14.8533,  93.1466,
     *                 12.7843,   5.6828,  98.0919, 829.4396,
     *                 50.9878, 199.0120,   2.0000,   6.0000,  10.0000,
     *                 10.0000,  30.0000,  50.0000,
     *                  8.0703,   5.7144,  84.1156, 529.0927,
     *                  5.6609,  28.9355, 111.3620, 494.0413,
     *                 45.5249, 134.4751,
     *                  4.0003,  21.2937,  78.7058,
     *                 12.8293,  16.2730, 123.6578, 327.2396,
     *                 48.7883, 102.2117,  20.0060, 161.9903,
     *                 28.4184,  61.5816,
     *                 10.5563,  13.2950, 188.1390,
     *                 14.6560, 129.4922, 470.8512                    /
      DATA GO       /   0.022 ,   2.019 ,   9.812 ,  13.087 ,
     *                 13.804 ,  16.061 ,  14.293 ,  16.114 ,
     *                  3.472 ,   7.437 ,  22.579 ,  32.035 ,
     *                 27.774 ,  33.678 ,  28.118 ,  31.019 ,  34.204 ,
     *                 30.892 ,  33.189 ,  36.181 ,
     *                  0.032 ,   2.760 ,  35.328 ,  48.277 ,
     *                  7.662 ,  16.786 ,  42.657 ,  54.522 ,
     *                 50.204 ,  56.044 ,
     *                  0.048 ,  50.089 ,  66.604 ,
     *                  8.954 ,  18.031 ,  57.755 ,  72.594 ,
     *                 68.388 ,  82.397 ,  31.960 ,  76.876 ,
     *                 75.686 ,  80.388 ,
     *                 10.747 ,  52.323 ,  94.976 ,
     *                 27.405 ,  86.350 , 109.917                     /
C
      DATA AF       /   2.0001,  39.9012, 122.0986,
     *                 10.0000,  30.0000,  50.0000,
     *                  4.0199,   5.5741,  22.1839, 190.2179,
     *                 53.0383, 126.9616,  31.6894,  75.3105,
     *                 13.5014,   7.9936,  55.7981, 298.7039,
     *                 26.2496,  63.7503,   2.0000,   6.0000,  10.0000,
     *                 28.7150,  71.2850,
     *                  8.0153,   6.1931,  21.7287,  48.7780, 278.2782,
     *                178.5560, 421.4435,  51.7632,  95.2368          /
      DATA GF       /   0.050 ,  13.317 ,  15.692 ,
     *                 15.361 ,  17.128 ,  18.498 ,
     *                  0.048 ,   2.735 ,  20.079 ,  30.277 ,
     *                 27.548 ,  32.532 ,  30.391 ,  34.707 ,
     *                  4.479 ,  12.072 ,  31.662 ,  51.432 ,
     *                 44.283 ,  50.964 ,  46.193 ,  50.436 ,  54.880 ,
     *                 50.816 ,  57.479 ,
     *                  0.058 ,   3.434 ,  14.892 ,  37.472 ,  69.883 ,
     *                 67.810 ,  83.105 ,  72.435 ,  79.747           /
C
      DATA ANN      /  34.5080, 365.4919,  16.5768, 183.4231,
     *                  2.0007,  89.5607, 380.4381,  26.4473,  63.5527,
     *                  4.0342,   5.6162,  11.5176,  72.8273,
     *                 48.5684, 131.4315,  31.1710,  76.8290,
     *                 14.0482,  13.3077,  52.7897, 467.8487,
     *                 54.2196, 195.7800                              /
      DATA GNN      /  17.796 ,  20.730 ,  17.879 ,  20.855 ,
     *                  0.097 ,  29.878 ,  37.221 ,  31.913 ,  37.551 ,
     *                  0.092 ,   3.424 ,  24.806 ,  46.616 ,
     *                 45.643 ,  54.147 ,  48.359 ,  57.420 ,
     *                  5.453 ,  18.560 ,  46.583 ,  80.101 ,
     *                 70.337 ,  85.789                               /
C
      DATA ANA      /  11.6348, 158.3593,
     *                 21.0453,  50.9546,  10.1389,  25.8611,
     *                  2.0019,  38.0569, 137.9398,  28.3106,  61.6893,
     *                  4.0334,   5.8560,  18.1786, 208.9142,
     *                 93.6895, 406.3095,  60.4276, 239.5719          /
      DATA GNA      /   2.400 ,   4.552 ,
     *                 34.367 ,  40.566 ,  34.676 ,  40.764 ,
     *                  0.170 ,  44.554 ,  57.142 ,  51.689 ,  60.576 ,
     *                  0.152 ,   4.260 ,  36.635 ,  83.254 ,
     *                 72.561 ,  89.475 ,  75.839 ,  92.582           /
C
      DATA AMG      /  10.7445, 291.5057,  53.7488,
     *                  6.2270,  31.1291, 132.6438,
     *                 40.4379, 159.5618,  20.3845,  79.6154,
     *                  2.0007, 106.8977, 343.1010,  10.1326, 237.8581/
      DATA GMG      /   2.805 ,   6.777 ,   9.254 ,
     *                  4.459 ,   9.789 ,  13.137 ,
     *                 57.413 ,  71.252 ,  58.010 ,  71.660 ,
     *                  0.276 ,  74.440 ,  94.447 ,  54.472 ,  95.858 /
C
      DATA AAL      /   4.0009,  11.7804, 142.2179,  13.6585,  96.3371,
     *                 10.0807,  49.5843, 285.3343,  14.6872,  59.3122,
     *                  6.3277,  29.5086, 134.1634,
     *                 46.3164, 153.6833,  22.9896,  77.0103          /
      DATA GAL      /   0.014 ,   3.841 ,   5.420 ,   3.727 ,   8.833 ,
     *                  4.749 ,  11.902 ,  16.719 ,  11.310 ,  18.268 ,
     *                  6.751 ,  16.681 ,  24.151 ,
     *                 83.551 , 104.787 ,  84.293 , 105.171           /
C
      DATA ASI      /   7.9658,   4.6762,   1.3512, 123.2267, 443.7797,
     *                  4.0000,   7.4186,  24.1754, 60.4060,
     *                 14.4695,  11.9721,  26.5062, 269.0521,
     *                  9.1793,   4.8766,  29.1442,  52.7998,
     *                 13.2674,  36.0417, 180.6910,
     *                  6.4839,  27.6851, 135.8301                    /
      DATA GSI      /   0.020 ,   0.752 ,   1.614 ,   5.831 ,   7.431 ,
     *                  0.036 ,   8.795 ,  11.208 ,  13.835 ,
     *                  5.418 ,   7.825 ,  14.440 ,  19.412 ,
     *                  6.572 ,  11.449 ,  18.424 ,  25.457 ,
     *                 15.682 ,  27.010 ,  34.599 ,
     *                  9.042 ,  24.101 ,  37.445                     /
C
      DATA AP       /  13.5211,  22.2130, 353.2583,  10.0000, 150.0000,
     *                  8.0241,   5.8085,  51.7542, 252.4002,
     *                  4.0021,  20.7985,  62.4194, 200.7786,
     *                 11.7414,  63.5124, 179.7420,
     *                  6.8835,  32.7777, 228.3366                    /
      DATA GP       /   1.514 ,   5.575 ,   9.247 ,   8.076 ,  10.735 ,
     *                  0.043 ,   1.212 ,   8.545 ,  15.525 ,
     *                  0.074 ,   7.674 ,  16.639 ,  25.118 ,
     *                  8.992 ,  24.473 ,  40.704 ,
     *                 11.464 ,  33.732 ,  55.455                     /
C
      DATA AS       /   3.9615,   5.0780,  15.0944, 362.8588,
     *                 51.5995, 268.4002,  12.0000, 276.0000,
     *                 11.4377,   5.5126, 141.0009, 254.0478,
     *                 33.0518, 126.9479,
     *                  4.0707,   4.0637,   5.7245, 144.6376, 106.4909,
     *                  4.0011,  19.2813,  27.5990,  35.1179,
     *                 94.7454, 283.2486,
     *                 10.5474,  28.7137,  65.7378,  24.0000          /
      DATA GS       /   0.053 ,   1.121 ,   5.812 ,   9.425 ,
     *                  8.936 ,  11.277 ,   9.600 ,  12.551 ,
     *                  1.892 ,   3.646 ,  13.550 ,  19.376 ,
     *                 16.253 ,  21.062 ,
     *                  0.043 ,   0.123 ,   1.590 ,  13.712 ,  22.050 ,
     *                  0.118 ,   9.545 ,  18.179 ,  31.441 ,
     *                 30.664 ,  56.150 ,
     *                 10.704 ,  27.075 ,  50.599 ,  43.034           /
C
      DATA ACL      /   2.0007,  62.5048, 669.4942,  29.0259, 130.9740,
     *                  3.9064,   0.3993,   5.3570,  60.3424, 119.9913,
     *                138.1567, 278.8418, 102.3681, 158.6314,
     *                 12.6089,   5.9527, 110.5635, 262.8715,
     *                 69.2035, 100.7960,
     *                  7.3458,   5.6638,  44.1256, 202.7846,
     *                  4.0037,  21.8663,  40.5363,  57.5919          /
      DATA GCL      /   0.110 ,   9.919 ,  12.280 ,  11.017 ,  13.532 ,
     *                  0.092 ,   0.581 ,   1.620 ,  13.121 ,  19.787 ,
     *                 16.365 ,  21.988 ,  18.065 ,  23.594 ,
     *                  2.358 ,   5.708 ,  19.084 ,  30.683 ,
     *                 24.880 ,  33.229 ,
     *                  0.102 ,   1.391 ,  14.709 ,  36.968 ,
     *                  0.185 ,  11.783 ,  25.653 ,  44.698           /
C
      DATA AAR      /  43.6623, 324.3375,  20.8298, 163.1701,
     *                  2.0026, 137.4515, 258.5445,  62.8129, 149.1867,
     *                  4.0495,  14.4466,  46.8234, 124.6651,
     *                151.9828, 268.0157, 101.1302, 150.8691,
     *                 13.3718,   8.6528,  60.4614, 285.5072,
     *                  6.7655,   4.7684,  12.8631,  54.5260          /
      DATA GAR      /  12.638 ,  14.958 ,  12.833 ,  15.139 ,
     *                  0.178 ,  17.522 ,  23.584 ,  20.464 ,  25.150 ,
     *                  0.151 ,   1.561 ,  17.399 ,  30.871 ,
     *                 24.684 ,  33.978 ,  27.091 ,  36.481 ,
     *                  2.810 ,   8.877 ,  24.351 ,  44.489 ,
     *                  0.144 ,   1.160 ,  10.210 ,  27.178           /
C
      DATA AK       /  12.9782, 148.6673,   6.3493,
     *                 66.3444, 101.6553,   4.0001,  13.4465,  46.5534,
     *                  2.0171, 116.4767, 713.4965,  63.5907, 396.4079,
     *                  2.0000,  10.0000,  30.0000,
     *                  4.0702,   5.7791,  52.6795, 327.4539,
     *                 62.8604, 357.1331,  55.9337, 196.0646,
     *                 10.9275,   5.5398,  43.2761,  76.2560,
     *                 42.0000,  18.0000                              /
      DATA GK       /   1.871 ,   3.713 ,  18.172 ,
     *                 21.185 ,  27.705 ,   2.059 ,  23.709 ,  28.542 ,
     *                  0.273 ,  26.709 ,  39.640 ,  31.220 ,  41.865 ,
     *                 29.955 ,  37.557 ,  42.862 ,
     *                  0.228 ,   2.274 ,  21.703 ,  50.191 ,
     *                 32.145 ,  49.262 ,  34.155 ,  51.718 ,
     *                  3.043 ,   5.479 ,  20.547 ,  30.680 ,
     *                 36.275 ,  47.345                               /
C
      DATA ACA      /  18.2366,  27.5012, 149.2617,  94.5242, 705.4711,
     *                 11.8706,  14.0710, 106.0547,
     *                 57.2414, 110.7567,  29.8121,  54.1874,
     *                  2.0184,  97.5784, 282.3939, 209.1871, 252.8129/
      DATA GCA      /   2.050 ,   3.349 ,   5.321 ,   4.873 ,   7.017 ,
     *                  1.769 ,   5.109 ,   9.524 ,
     *                 27.271 ,  41.561 ,  29.172 ,  42.140 ,
     *                  0.394 ,  28.930 ,  52.618 ,  38.593 ,  49.646 /
C
      DATA ASC      /   6.0014,  83.1958,  67.3666, 329.4354,
     *                 44.0793, 169.9969, 533.9195,
     *                 34.1642, 124.8475, 228.9879,
     *                 11.9979,  16.9280,  28.4778,  82.0418, 234.5360,
     *                  6.0042,   2.7101,  13.9801,  65.3039,
     *                 12.0000,  12.0000,
     *                  2.0051,   2.9621,  29.0306                    /
      DATA GSC      /   0.021 ,   2.056 ,   3.551 ,   5.465 ,
     *                  1.535 ,   3.797 ,   6.203 ,
     *                  2.389 ,   4.858 ,   7.141 ,
     *                  0.011 ,   0.430 ,   1.156 ,   3.711 ,   8.863 ,
     *                  0.025 ,   3.499 ,  10.463 ,  18.606 ,
     *                 41.779 ,  57.217 ,
     *                  0.539 ,  24.442 ,  51.079                     /
C
      DATA ATI      /   7.0887,   8.9186,  17.5633, 206.6832, 438.5735,
     *                654.1721,
     *                 38.0462,  69.6271, 364.2845, 832.0408,
     *                 98.8562,  57.9934, 442.1498,
     *                 19.7843,  32.0637,  37.0895, 110.6682, 288.4946,
     *                521.8837,
     *                 10.0000,  34.0000, 120.0000,
     *                 16.1691,  22.3550,  24.1646,  83.5128, 222.7963,
     *                  6.0020,   4.6177,  25.2636,  52.1162,
     *                 12.0000,   8.0000                              /
      DATA GTI      /   0.021 ,   0.048 ,   1.029 ,   2.183 ,   4.109 ,
     *                  5.785 ,
     *                  0.846 ,   1.792 ,   3.836 ,   5.787 ,
     *                  2.561 ,   4.869 ,   6.340 ,
     *                  0.023 ,   0.124 ,   0.774 ,   1.810 ,   4.980 ,
     *                  9.585 ,
     *                  1.082 ,   4.928 ,  11.279 ,
     *                  0.041 ,   1.375 ,   4.768 ,  10.985 ,  19.769 ,
     *                  0.048 ,  11.577 ,  24.531 ,  36.489 ,
     *                 54.436 ,  75.373                               /
C
      DATA AV       /  15.2627,  23.9869,  51.3053, 570.3384,1650.9417,
     *                162.2829, 298.8303, 908.8852,
     *                 23.6736,  37.1624,  86.8011, 300.7440, 864.5880,
     *                 57.8961,  79.4605, 214.9007, 864.7425,
     *                 61.8508,  64.0845, 192.8298, 718.2349,
     *                 23.8116,  68.2495, 135.0613, 536.7632,
     *                 15.9543,  22.5542,  71.4921, 248.9544,
     *                  6.0006,   5.8785,  50.5077,  97.6129          /
      DATA GV       /   0.026 ,   0.145 ,   0.718 ,   2.586 ,   5.458 ,
     *                  2.171 ,   4.153 ,   6.097 ,
     *                  0.009 ,   0.366 ,   1.504 ,   5.294 ,  10.126 ,
     *                  1.796 ,   2.353 ,   6.068 ,  12.269 ,
     *                  2.560 ,   3.674 ,   6.593 ,  12.880 ,
     *                  0.045 ,   1.684 ,   8.162 ,  21.262 ,
     *                  0.065 ,   1.746 ,  15.158 ,  33.141 ,
     *                  0.077 ,  21.229 ,  44.134 ,  60.203           /
C
      DATA ACR      /  30.1842,  79.2847, 149.5293,
     *                215.3696, 119.1974, 741.4321,
     *                184.9946,1352.5038, 784.4937,
     *                 46.6191, 160.1361, 488.0449, 657.1928,
     *                 47.1742, 267.0275, 441.1324, 150.6650,
     *                 24.3768, 122.8359, 285.5092, 794.1654,
     *                 24.2296,  75.0258, 172.9452, 543.6511,
     *                 15.9819,  17.6800,  95.2003, 225.0947          /
      DATA GCR      /   0.993 ,   3.070 ,   5.673 ,
     *                  3.339 ,   4.801 ,   7.198 ,
     *                  2.829 ,   4.990 ,   7.643 ,
     *                  1.645 ,   3.727 ,   7.181 ,  12.299 ,
     *                  2.902 ,   4.273 ,   8.569 ,  14.912 ,
     *                  0.047 ,   2.566 ,   9.441 ,  21.198 ,
     *                  0.078 ,   2.242 ,  15.638 ,  32.725 ,
     *                  0.103 ,   2.146 ,  26.153 ,  49.381           /
C
      DATA AMN      /  53.9107,  81.3931, 546.6945 ,
     *                144.1893, 407.8029,  45.6177, 298.4423,2410.9335,
     *                 22.6382,  93.8419, 183.9367, 907.5765,
     *                137.0409, 168.6783, 329.0287, 773.2513,
     *                 70.1925,  72.3372, 213.9512, 539.5165,
     *                 24.2373,  93.5415, 456.6167, 506.5484,
     *                 24.7687,  66.9896, 264.1853, 484.0161          /
      DATA GMN      /   2.527 ,   4.204 ,   6.602 ,
     *                  4.155 ,   7.321 ,   2.285 ,   5.631 ,   8.448 ,
     *                  1.496 ,   3.839 ,   7.751 ,  13.484 ,
     *                  3.681 ,   6.054 ,   9.934 ,  14.936 ,
     *                  3.531 ,   6.967 ,  15.222 ,  25.069 ,
     *                  0.071 ,   2.896 ,  20.725 ,  37.383 ,
     *                  0.126 ,   2.660 ,  28.528 ,  53.413           /
C
      DATA AFE      /  14.4102,   2.7050, 421.6612, 940.1484,
     *                 36.2187,  22.8883, 239.5997, 825.2919,
     *                110.0242, 992.3040, 640.6715,
     *                 17.0494,  32.3783,  34.3184, 420.9626,1067.2064,
     *                154.0059, 462.1117, 329.8618,
     *                 15.7906,  47.1186, 279.9292, 692.1005,
     *                 91.0206, 206.3082, 706.9927, 836.6689,
     *                 40.0790,  27.6965,  28.2243,  18.0001,
     *                 24.0899,  89.6340,  51.5756, 241.6980          /
      DATA GFE      /   0.066 ,   0.339 ,   2.897 ,   6.585 ,
     *                  0.923 ,   1.679 ,   4.620 ,   7.053 ,
     *                  4.249 ,   5.875 ,   7.781 ,
     *                  0.062 ,   0.283 ,   1.504 ,   5.430 ,  11.210 ,
     *                  2.792 ,   7.627 ,  13.623 ,
     *                  0.077 ,   3.723 ,  12.137 ,  23.700 ,
     *                  2.688 ,   7.595 ,  15.444 ,  25.587 ,
     *                  3.982 ,   4.677 ,   6.453 ,  23.561 ,
     *                  0.102 ,   3.354 ,  22.954 ,  33.796           /
C
      DATA ACO      /  11.9120,  20.4424,  28.3863, 132.5038, 600.7461,
     *                 33.3092, 237.4331, 977.2502,
     *                 55.5396, 318.8169, 619.6366,
     *                 32.6900,  83.8694, 107.4378,
     *                 11.2593,  38.2239,  22.9964, 261.3486, 637.1485,
     *                 23.0233,  41.6599, 264.6460, 181.6699,
     *                 16.0356,   7.8633,  70.3158, 423.3512, 742.3553,
     *                  0.                                            /
      DATA GCO      /   0.112 ,   0.341 ,   0.809 ,   3.808 ,   6.723 ,
     *                  2.057 ,   3.484 ,   7.210 ,
     *                  2.405 ,   5.133 ,   8.097 ,
     *                  2.084 ,   5.291 ,   8.426 ,
     *                  0.135 ,   0.517 ,   1.606 ,   6.772 ,  12.622 ,
     *                  2.512 ,   4.348 ,   8.253 ,  15.377 ,
     *                  0.132 ,   0.863 ,   3.086 ,  11.789 ,  23.263 ,
     *                  0.                                            /
C
      DATA ANI      /   7.1268,  12.4486,  11.9953,  10.0546, 114.1658,
     *                391.2064,
     *                 26.3908, 213.8081, 938.7927,
     *                  4.1421,  37.3781,  25.9712, 333.3397, 311.1633,
     *                 33.1031, 184.1854, 136.7072,
     *                 11.1915,   5.4174,  53.6793, 460.6781, 380.0056,
     *                  0.                                            /
      DATA GNI      /   0.026 ,   0.137 ,   0.315 ,   1.778 ,   4.029 ,
     *                  6.621 ,
     *                  2.249 ,   4.042 ,   7.621 ,
     *                  0.191 ,   1.235 ,   3.358 ,   8.429 ,  17.096 ,
     *                  3.472 ,   9.065 ,  16.556 ,
     *                  0.194 ,   1.305 ,   5.813 ,  14.172 ,  26.169 ,
     *                  0.                                            /
C
      DATA ACU      /  11.0549, 238.9423,  10.3077, 126.2990,1073.3876,
     *                 30.0000,  50.0000,  60.0000,
     *                 19.2984,  50.5974, 240.2021,1216.9016,
     *                 48.3048, 583.2011, 320.4931,
     *                  4.0155,  70.3264, 313.1213, 536.5331,
     *                  0.                                            /
      DATA GCU      /   4.212 ,   7.227 ,   1.493 ,   5.859 ,   9.709 ,
     *                  7.081 ,   9.362 ,  10.130 ,
     *                  2.865 ,   8.260 ,  14.431 ,  18.292 ,
     *                  9.650 ,  14.640 ,  24.320 ,
     *                  0.337 ,   8.520 ,  16.925 ,  28.342 ,
     *                  0.                                            /
C
      DATA AZN      /  15.9880, 484.0042,  18.5863, 123.4134,
     *                  3.0000, 189.0000,
     *                  6.1902,  38.9317, 204.8780,
     *                 10.2588,  89.3771, 370.3640,  30.0000, 128.0000,
     *                 24.6904, 106.7491, 439.5586,
     *                  0.                                            /
      DATA GZN      /   4.546 ,   8.840 ,  10.247 ,  16.620 ,
     *                 11.175 ,  16.321 ,
     *                  6.113 ,  12.964 ,  16.444 ,
     *                  7.926 ,  13.633 ,  24.353 ,  16.286 ,  24.910 ,
     *                 10.291 ,  20.689 ,  32.077 ,
     *                  0.                                            /
C
C
      DATA NIONS,NSS,ICOMP /123,222,0/
C
      EQUIVALENCE   ( AHH(1), ALF(  1)),( ALB(1), ALF(  7)),
     *              ( AB (1), ALF( 19)),
     *              ( AC (1), ALF( 30)),( AN (1), ALF( 49)),
     *              ( AO (1), ALF( 79)),( AF (1), ALF(128)),
     *              ( ANN(1), ALF(162)),( ANA(1), ALF(185)),
     *              ( AMG(1), ALF(204)),( AAL(1), ALF(219)),
     *              ( ASI(1), ALF(236)),( AP (1), ALF(259)),
     *              ( AS (1), ALF(278)),( ACL(1), ALF(307)),
     *              ( AAR(1), ALF(335)),( AK (1), ALF(360)),
     *              ( ACA(1), ALF(390)),( ASC(1), ALF(407)),
     *              ( ATI(1), ALF(431)),( AV (1), ALF(464)),
     *              ( ACR(1), ALF(497)),( AMN(1), ALF(526)),
     *              ( AFE(1), ALF(554)),( ACO(1), ALF(589)),
     *              ( ANI(1), ALF(618)),( ACU(1), ALF(641)),
     *              ( AZN(1), ALF(661))
      EQUIVALENCE   ( GHH(1), GAM(  1)),( GLB(1), GAM(  7)),
     *              ( GB (1), GAM( 19)),
     *              ( GC (1), GAM( 30)),( GN (1), GAM( 49)),
     *              ( GO (1), GAM( 79)),( GF (1), GAM(128)),
     *              ( GNN(1), GAM(162)),( GNA(1), GAM(185)),
     *              ( GMG(1), GAM(204)),( GAL(1), GAM(219)),
     *              ( GSI(1), GAM(236)),( GP (1), GAM(259)),
     *              ( GS (1), GAM(278)),( GCL(1), GAM(307)),
     *              ( GAR(1), GAM(335)),( GK (1), GAM(360)),
     *              ( GCA(1), GAM(390)),( GSC(1), GAM(407)),
     *              ( GTI(1), GAM(431)),( GV (1), GAM(464)),
     *              ( GCR(1), GAM(497)),( GMN(1), GAM(526)),
     *              ( GFE(1), GAM(554)),( GCO(1), GAM(589)),
     *              ( GNI(1), GAM(618)),( GCU(1), GAM(641)),
     *              ( GZN(1), GAM(661))
      EQUIVALENCE   ( CH1(1), CHION(  1)),
     *              ( CH2(1), CHION( 67)),
     *              ( CH3(1), CHION(139)),
     *              ( CH4(1), CHION(194)),
     *              ( XL1(1),    XL(  1)),
     *              ( XL2(1),    XL(100))
      EQUIVALENCE   ( IS1(1),  IS(1)),   ( IS2(1),  IS( 54)),
     *              ( IM1(1),  IM(1)),   ( IM2(1),  IM(100)),
     *              (IGP1(1),IGPR(1)),   (IGP2(1),IGPR(100)),
     *              (IG01(1), IG0(1)),   (IG02(1), IG0( 54)),
     *              (II1(1,1),INDEX0(1,1)),(II2(1,1),INDEX0(1,16))
C
C
      EXP10(X)=EXP(2.3025851*X)
C
      IF(ICOMP.NE.0) GO TO 5
      IND=1
      DO 1 K=1,NIONS
         INDEXS(K)=IND
         IND=IND+IS(K)
    1 CONTINUE
      IND=1
      DO 2 K=1,NSS
         INDEXM(K)=IND
         IND=IND+IM(K)
    2 CONTINUE
      ICOMP=1
    5 CONTINUE
      IF(IAT.EQ.26.AND.IZI.GE.4.AND.IZI.LE.9) GO TO 70
CMH      IF(IZI.LE.0.OR.IZI.GT.5.OR.IAT.LE.0.OR.IAT.GT.30) GO TO 50
CMH  new consider IAT up to 90
      IF(IZI.LE.0.OR.IZI.GT.5.OR.IAT.LE.0.OR.IAT.GT.90) GO TO 50
      MODE=MODPF(IAT)
c      prints *,'part function:',mode,izi,iat
      IF(MODE.LT.0) GO TO 50
      IF(MODE.GT.0) GO TO 60
      I0=INDEX0(IZI,IAT)
      IF(I0) 40,50,10
   10 QZ=IZI
      MAX=XMAXN*SQRT(QZ)
      THET=5040.4/T
      A=31.321*QZ*QZ*THET
      MAX2=MAX*MAX
      QAS1=MAX/3.*(MAX2+1.5*MAX+0.5)
      IS0=INDEXS(I0)
      ISS=IS0+IS(I0)-1
      SU1=0.
      SQA=0.
      DO 30 K=IS0,ISS
         XXL=XL(K)
         GPR=IGPR(K)
         X=CHION(K)*THET
         EX=0.
         IF(X.LT.30) EX=EXP10(-X)
         QAS=(QAS1-XXL/3.*(XXL*XXL+1.5*XXL+0.5)+(MAX-XXL)*(1.+A/2./XXL/
     1        MAX)*A)*GPR*EX
         SQA=SQA+QAS
         M0=INDEXM(K)
         M1=M0+IM(K)-1
         AL1=0.
         DO 20 M=M0,M1
         XG=GAM(M)*THET
         IF(XG.GT.20.) GO TO 20
         XM=EXP10(-XG)*ALF(M)
         AL1=AL1+XM
   20    CONTINUE
         SU1=SU1+AL1
   30 CONTINUE
      U=IG0(I0)
      U=U+SU1+SQA
      IF(U.LT.0.) U=IG0(I0)
      RETURN
   40 U=FLOAT(-I0)
      RETURN
   50 CALL PFSPEC(IAT,IZI,T,ANE,U)
      RETURN
   60 U=PFSTD(IZI,IAT)
      RETURN
   70 call pffe(IZI,T,ANE,U)
      RETURN
      END SUBROUTINE
C
C ********************************************************************
C
C
C
      subroutine pffe(ion,t,ane,pf)
c     =============================
c
c     partition functions for Fe IV to Fe IX
c     after Fischel and Sparks, 1971, NASA SP-3066
c
c     Output:  PF   partition function
c
c      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: ion
      real*8, intent(in   ) :: t,ane
      real*8, intent(inout) :: pf
        
      dimension tt(50),pn(10),nca(9)
      dimension p4a(22),p4b(10,28), 
     *          p5a(30),p5b(10,20),
     *          p6a(37),p6b(10,13),
     *          p7a(40),p7b(10,10),
     *          p8a(41),p8b(10,9),
     *          p9a(45),p9b(10,5)
c
c
      data nca /3*0,22,30,37,40,41,45/
     *     nne /10/
c
      data tt /
     * 3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,
     * 20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,
     * 32.,34.,36.,38.,40.,42.,44.,46.,48.,
     * 50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,125.,150./
c
      data pn /-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7./
c
      data p4a /
     * 0.778, 0.778, 0.778, 0.779, 0.783, 0.789, 0.801, 0.818,
     * 0.842, 0.871, 0.906, 0.945, 0.987, 1.030, 1.074, 1.117,
     * 1.160, 1.201, 1.242, 1.280, 1.317, 1.353/
c
      data p4b /
     * 1.406,1.393,1.389,7*1.387,
     * 1.464,1.434,1.424,1.421,1.420,5*1.419,
     * 1.546,1.483,1.461,1.454,1.451,1.451,4*1.450,
     * 1.665,1.547,1.503,1.488,1.482,1.481,4*1.480,
     * 1.826,1.636,1.553,1.524,1.514,1.510,4*1.509,
     * 2.024,1.755,1.618,1.564,1.546,1.540,1.538,3*1.537,
     * 2.480,2.087,1.814,1.674,1.619,1.599,1.593,1.591,1.590,1.590,
     * 2.945,2.489,2.105,1.846,1.717,1.667,1.649,1.643,1.641,1.640,
     * 3.379,2.897,2.452,2.089,1.859,1.751,1.710,1.696,1.691,1.689,
     * 3.774,3.283,2.808,2.381,2.054,1.864,1.782,1.751,1.741,1.738,
     * 4.133,3.637,3.150,2.688,2.292,2.015,1.871,1.814,1.793,1.786,
     * 4.460,3.962,3.468,2.989,2.549,2.199,1.984,1.886,1.848,1.835,
     * 4.757,4.258,3.762,3.274,2.809,2.406,2.121,1.972,1.908,1.886,
     * 5.029,4.530,4.032,3.539,3.061,2.624,2.279,2.073,1.976,1.939,
     * 5.279,4.780,4.281,3.785,3.299,2.840,2.450,2.189,2.051,1.996,
     * 5.510,5.010,4.511,4.013,3.522,3.050,2.628,2.318,2.136,2.057,
     * 6.014,5.514,5.014,4.515,4.018,3.530,3.065,2.666,2.381,2.228,
     * 6.435,5.935,5.435,4.936,4.437,3.943,3.460,3.022,2.658,2.422,
     * 6.794,6.294,5.794,5.294,4.794,4.297,3.807,3.343,2.939,2.631,
     * 7.102,6.602,6.102,5.602,5.102,4.604,4.110,3.638,3.194,2.845,
     * 7.370,6.870,6.370,5.870,5.370,4.871,4.375,3.892,3.439,3.052,
     * 7.606,7.106,6.606,6.106,5.605,5.106,4.608,4.125,3.661,3.249,
     * 7.815,7.315,6.814,6.314,5.814,5.314,4.816,4.333,3.851,3.418,
     * 8.001,7.501,7.001,6.500,6.000,5.500,5.001,4.511,4.032,3.586,
     * 8.168,7.668,7.168,6.668,6.168,5.667,5.168,4.680,4.197,3.741,
     * 8.319,7.819,7.319,6.819,6.319,5.818,5.319,4.832,4.347,3.884,
     * 8.900,8.399,7.899,7.399,6.899,6.398,5.898,5.405,4.917,4.431,
     * 9.294,8.794,8.294,7.793,7.293,6.793,6.292,5.799,5.306,4.824/
c
      data p5a /
     * 1.235, 1.276, 1.301, 1.321, 1.339, 1.359, 1.381, 1.405,
     * 1.432, 1.460, 1.489, 1.518, 1.546, 1.574, 1.601, 1.627,
     * 1.652, 1.675, 1.697, 1.718, 1.738, 1.757, 1.775, 1.792,
     * 1.808, 1.823, 1.838, 1.851, 1.877, 1.900/
c
      data p5b /
     * 1.943,1.928,1.923,7*1.921,
     * 2.011,1.964,1.947,1.942,1.941,5*1.940,
     * 2.144,2.025,1.980,1.965,1.960,1.958,4*1.957,
     * 2.361,2.137,2.032,1.993,1.980,1.976,1.975,3*1.974,
     * 2.646,2.315,2.121,2.035,2.004,1.994,1.991,1.990,1.989,1.989,
     * 2.960,2.553,2.260,2.102,2.037,2.015,2.007,2.005,2.004,2.004,
     * 3.274,2.823,2.450,2.205,2.086,2.040,2.025,2.020,2.018,2.018,
     * 3.575,3.101,2.674,2.348,2.158,2.075,2.045,2.036,2.032,2.031,
     * 4.251,3.757,3.275,2.829,2.466,2.234,2.124,2.083,2.069,2.064,
     * 4.822,4.324,3.829,3.346,2.895,2.522,2.278,2.161,2.116,2.100,
     * 5.308,4.808,4.310,3.816,3.334,2.888,2.525,2.297,2.187,2.145,
     * 5.725,5.225,4.726,4.228,3.736,3.260,2.828,2.496,2.294,2.206,
     * 6.088,5.589,5.089,4.590,4.093,3.604,3.139,2.733,2.447,2.291,
     * 6.407,5.907,5.407,4.908,4.409,3.915,3.433,2.988,2.629,2.399,
     * 6.689,6.189,5.689,5.189,4.690,4.193,3.704,3.236,2.832,2.535,
     * 6.940,6.440,5.940,5.440,4.941,4.443,3.949,3.469,3.038,2.687,
     * 7.166,6.666,6.166,5.666,5.166,4.667,4.171,3.684,3.237,2.847,
     * 7.370,6.870,6.369,5.869,5.369,4.870,4.373,3.882,3.417,3.008,
     * 8.150,7.649,7.149,6.649,6.149,5.649,5.149,4.651,4.167,3.700,
     * 8.677,8.177,7.676,7.176,6.676,6.176,5.676,5.176,4.687,4.203/
c
      data p6a /
     * 1.218, 1.273, 1.309, 1.335, 1.358, 1.379, 1.400, 1.421,
     * 1.442, 1.463, 1.484, 1.504, 1.523, 1.542, 1.560, 1.577,
     * 1.594, 1.609, 1.624, 1.638, 1.652, 1.664, 1.677, 1.688,
     * 1.699, 1.709, 1.719, 1.729, 1.746, 1.762, 1.777, 1.790,
     * 1.803, 1.814, 1.825, 1.834, 1.843/
c
      data p6b /
     * 1.862,1.855,1.853,7*1.852,
     * 1.958,1.900,1.880,1.874,1.872,5*1.871,
     * 2.264,2.045,1.944,1.906,1.894,1.890,4*1.888,
     * 2.776,2.386,2.119,1.984,1.930,1.912,1.906,1.904,2*1.903,
     * 3.321,2.856,2.453,2.165,2.012,1.949,1.927,1.920,1.918,1.917,
     * 3.821,3.333,2.868,2.465,2.178,2.025,1.963,1.941,1.934,1.932,
     * 4.266,3.771,3.285,2.825,2.434,2.164,2.027,1.972,1.953,1.947,
     * 4.662,4.164,3.670,3.187,2.739,2.372,2.135,2.022,1.980,1.965,
     * 5.015,4.516,4.019,3.527,3.052,2.624,2.295,2.102,2.019,1.988,
     * 5.332,4.832,4.344,3.838,3.351,2.889,2.493,2.217,2.075,2.017,
     * 5.618,5.118,4.619,4.121,3.628,3.149,2.711,2.364,2.155,2.058,
     * 6.710,6.210,5.710,5.210,4.711,4.213,3.719,3.241,2.807,2.462,
     * 7.446,6.946,6.446,5.946,5.446,4.946,4.447,3.952,3.474,3.022/
c
      data p7a /
     * 1.074,1.130,1.167,1.194,1.215,1.234,1.250,1.266,1.280,1.293,
     * 1.306,1.318,1.329,1.340,1.350,1.360,1.369,1.378,1.386,1.394,
     * 1.401,1.408,1.415,1.421,1.427,1.433,1.439,1.444,1.454,1.463,
     * 1.471,1.479,1.486,1.492,1.498,1.504,1.509,1.514,1.525,1.534/
c
      data p7b /
     * 1.555,1.546,1.544,1.543,6*1.542,
     * 1.617,1.572,1.557,1.552,1.550,1.550,4*1.549,
     * 1.798,1.648,1.587,1.566,1.559,1.557,4*1.556,
     * 2.134,1.832,1.666,1.597,1.573,1.565,1.563,1.562,2*1.561,
     * 2.550,2.138,1.836,1.671,1.602,1.578,1.570,1.568,2*1.567,
     * 2.968,2.504,2.102,1.816,1.665,1.603,1.582,1.575,2*1.572,
     * 3.359,2.875,2.419,2.037,1.779,1.651,1.601,1.584,1.579,1.577,
     * 3.718,3.224,2.745,2.305,1.953,1.736,1.636,1.599,1.586,1.582,
     * 5.097,4.598,4.098,3.601,3.110,2.638,2.217,1.899,1.719,1.643,
     * 6.026,5.526,5.026,4.527,4.028,3.531,3.042,2.576,2.170,1.885/
c
      data p8a /
     * 0.809,0.849,0.875,0.894,0.908,0.918,0.927,0.934,0.939,0.944,
     * 0.948,0.952,0.955,0.958,0.960,0.962,0.964,0.966,0.967,0.969,
     * 0.970,0.971,0.973,0.974,0.975,0.975,0.976,0.977,0.978,0.980,
     * 0.981,0.982,0.983,0.984,0.984,0.985,0.986,0.986,0.987,0.988,
     * 0.989/
c
      data p8b /
     * 0.992,0.991,8*0.990,
     * 1.000,0.994,0.992,7*0.991,
     * 1.032,1.005,0.996,0.993,0.992,5*0.991,
     * 1.129,1.040,1.008,0.997,0.993,5*0.992,
     * 1.335,1.132,1.042,1.009,0.998,0.994,0.993,0.993,2*0.992,
     * 1.640,1.312,1.121,1.038,1.007,0.998,0.994,3*0.993,
     * 1.987,1.573,1.269,1.101,1.030,1.005,0.997,2*0.994,0.993,
     * 3.514,3.017,2.526,2.053,1.628,1.305,1.119,1.039,1.010,1.000,
     * 4.569,4.069,3.569,3.072,2.580,2.103,1.671,1.336,1.136,1.048/
c
      data p9a /39*0.000,0.001,0.002,0.005,0.008,0.014,0.021/
c
      data p9b /
     * 2*0.032,8*0.031,
     * 0.048,0.045,8*0.044,
     * 0.076,0.065,0.061,0.060,6*0.059,
     * 1.128,0.722,0.429,0.271,0.207,0.184,0.177,0.174,2*0.173,
     * 2.696,2.200,1.712,1.249,0.848,0.564,0.415,0.354,0.333,0.327/
c
      parameter (xen=2.302585093,xmil=0.001,xmilen=xmil*xen)
      parameter (xbtz=1.38054d-16)
c
      na=nca(ion)
      nb=50-na
      pne=log10(ane*xbtz*t)
      t0=xmil*t
      j=1
      if(pne.lt.pn(1)) go to 15
      if(pne.gt.pn(nne)) then
        j1=nne
        j2=nne
        goto 16
      endif
      do 10 j=1,nne-1      
        if(pne.ge.pn(j).and.pne.lt.pn(j+1)) go to 15
   10 continue
   15 j1=j
      j2=j1+1
      if(pne.lt.pn(1)) j2=1
   16 do 20 i=1,49
         if(t0.ge.tt(i).and.t0.lt.tt(i+1)) go to 25
   20 continue
   25 i1=i
      i2=i+1
      if(t0.gt.tt(50)) then
        i1=50
        i2=50
      endif
      if(i2.le.na) then
         if(ion.eq.4) then
           px1=p4a(i1)
           px2=p4a(i1)
           py1=p4a(i2)
           py2=p4a(i2)
         else if(ion.eq.5) then
           px1=p5a(i1)
           px2=p5a(i1)
           py1=p5a(i2)
           py2=p5a(i2)
         else if(ion.eq.6) then
           px1=p6a(i1)
           px2=p6a(i1)
           py1=p6a(i2)
           py2=p6a(i2)
         else if(ion.eq.7) then
           px1=p7a(i1)
           px2=p7a(i1)
           py1=p7a(i2)
           py2=p7a(i2)
         else if(ion.eq.8) then
           px1=p8a(i1)
           px2=p8a(i1)
           py1=p8a(i2)
           py2=p8a(i2)
         else if(ion.eq.9) then
           px1=p9a(i1)
           px2=p9a(i1)
           py1=p9a(i2)
           py2=p9a(i2)
         endif
       else if(i1.eq.na) then
         if(ion.eq.4) then
           px1=p4a(i1)
           px2=p4a(i1)
           py1=p4b(j1,i2-na)
           py2=p4b(j2,i2-na)
         else if(ion.eq.5) then
           px1=p5a(i1)
           px2=p5a(i1)
           py1=p5b(j1,i2-na)
           py2=p5b(j2,i2-na)
         else if(ion.eq.6) then
           px1=p6a(i1)
           px2=p6a(i1)
           py1=p6b(j1,i2-na)
           py2=p6b(j2,i2-na)
         else if(ion.eq.7) then
           px1=p7a(i1)
           px2=p7a(i1)
           py1=p7b(j1,i2-na)
           py2=p7b(j2,i2-na)
         else if(ion.eq.8) then
           px1=p8a(i1)
           px2=p8a(i1)
           py1=p8b(j1,i2-na)
           py2=p8b(j2,i2-na)
         else if(ion.eq.9) then
           px1=p9a(i1)
           px2=p9a(i1)
           py1=p9b(j1,i2-na)
           py2=p9b(j2,i2-na)
         endif
      else
         if(ion.eq.4) then
           px1=p4b(j1,i1-na)
           px2=p4b(j2,i1-na)
           py1=p4b(j1,i2-na)
           py2=p4b(j2,i2-na)
         else if(ion.eq.5) then
           px1=p5b(j1,i1-na)
           px2=p5b(j2,i1-na)
           py1=p5b(j1,i2-na)
           py2=p5b(j2,i2-na)
         else if(ion.eq.6) then
           px1=p6b(j1,i1-na)
           px2=p6b(j2,i1-na)
           py1=p6b(j1,i2-na)
           py2=p6b(j2,i2-na)
         else if(ion.eq.7) then
           px1=p7b(j1,i1-na)
           px2=p7b(j2,i1-na)
           py1=p7b(j1,i2-na)
           py2=p7b(j2,i2-na)
         else if(ion.eq.8) then
           px1=p8b(j1,i1-na)
           px2=p8b(j2,i1-na)
           py1=p8b(j1,i2-na)
           py2=p8b(j2,i2-na)
         else if(ion.eq.9) then
           px1=p9b(j1,i1-na)
           px2=p9b(j2,i1-na)
           py1=p9b(j1,i2-na)
           py2=p9b(j2,i2-na)
         endif
      end if
      dlgunx=px2-px1
      px=px1+(pne-pn(j1))*dlgunx
      dlguny=py2-py1
      py=py1+(pne-pn(j1))*dlguny
      delt=tt(i2)-tt(i1)
      if(delt.ne.0.) then
        dlgut=(py-px)/delt
        pf=px+(t0-tt(i1))*dlgut
      else
        pf=px
      endif
      pf=exp(xen*pf)
      return
      end subroutine
C
C ********************************************************************
C ********************************************************************
C
      !SUBROUTINE PFSPEC(IAT,IZI,T,ANE,XMAX,U)
      SUBROUTINE PFSPEC(IAT,IZI,T,ANE,U)
C     =======================================
                                                                                
C     Non-standard evaluation of the partition function
C     user supplied procedure
C
C     Input:
C      IAT   - atomic number
C      IZI   - ionic charge (=1 for neutrals, =1 for once ionized, etc)
C      T     - temperature
C      ANE   - electron density
C      XMAX  - principal quantum number of the last bound level
C
C     Output:
C      U     - partition function
C
*
* Modified from the ATMOS related programme 5-April-1990
* as an addition to TLUSTY to allow high ionisation states
* of C, N and O
*
* M.A.Barstow - University of Leicester, Dept of Physics & Astronomy
*
C      INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: IAT,IZI
      real*8, intent(in   ) :: T,ANE
!       integer,intent(inout) :: 
      real*8, intent(inout) :: U
!       integer,intent(  out) :: 
!       real*8, intent(  out) :: 
!       integer :: 
!       real*8  :: 
        
      PARAMETER (MH=100,MHEI=100,MHEII=100,MCI=135,
     +     MCII=157,MCIII=156,MCIV=55,MCV=15,MCVI=100,MNI=228,MNII=122,
     +     MNIII=133,MNIV=73,MNV=51,MNVI=8,MNVII=100,MOI=174,MOII=191,
     +     MOIII=168,MOIV=166,MOV=115,MOVI=52,MOVII=16,MOVIII=100)
      DIMENSION GHYD(MH),SHYD(MH),ENHYD(MH),
     +     GHEL(MH),ENHEL(MH),SHEL(MH),
     +     GCI(MCI),ENCI(MCI),SCI(MCI),
     +     GCII(MCII),ENCII(MCII),SCII(MCII),
     +     GCIII(MCIII),ENCIII(MCIII),SCIII(MCIII),
     +     GCIV(MCIV),ENCIV(MCIV),SCIV(MCIV),
     +     GCV(MCV),ENCV(MCV),SCV(MCV),
     +     GNI(MNI),ENNI(MNI),SNI(MNI),
     +     GNII(MNII),ENNII(MNII),SNII(MNII),
     +     GNIII(MNIII),ENNIII(MNIII),SNIII(MNIII),
     +     GNIV(MNIV),ENNIV(MNIV),SNIV(MNIV),
     +     GNV(MNV),ENNV(MNV),SNV(MNV),
     +     GNVI(MNVI),ENNVI(MNVI),SNVI(MNVI),
     +     GOI(MOI),ENOI(MOI),SOI(MOI),
     +     GOII(MOII),ENOII(MOII),SOII(MOII),
     +     GOIII(MOIII),ENOIII(MOIII),SOIII(MOIII),
     +     GOIV(MOIV),ENOIV(MOIV),SOIV(MOIV),
     +     GOV(MOV),ENOV(MOV),SOV(MOV),
     +     GOVI(MOVI),ENOVI(MOVI),SOVI(MOVI),
     +     GOVII(MOVII),ENOVII(MOVII),SOVII(MOVII)
      INTEGER NHYD(MH),NHEL(MHEI),NCI(MCI),NCII(MCII),
     +     NCIII(MCIII),NCIV(MCIV),NCV(MCV),NNI(MNI),NNII(MNII),
     +     NNIII(MNIII),NNIV(MNIV),NNV(MNV),NNVI(MNVI),NOI(MOI),
     +     NOII(MOII),NOIII(MOIII),NOIV(MOIV),NOV(MOV),NOVI(MOVI),
     +     NOVII(MOVII)
      PARAMETER (HI=13.5878,HEI=24.587,HEII=54.416,CVI=489.84,
     +     NVII=666.83,OVIII=871.12)
      PARAMETER (ZH=1.0,ZHE=2.0,ZC=6.0,ZN=7.0,ZO=8.0)
C                           N***=QUANTUM NO. OF LEVEL
C      DATA FOR IONS        G***=STATISTICAL WEIGHT OF LEVEL
C                           EN***=ENERGY OF LEVEL
C                           S*=SCREENING NO. OF LEVEL
	DATA NHYD/ 1, 2, 3, 4, 5, 6,
     +	 	7, 8, 9,10,11,12,
     +		13,14,15,16,17,18,
     +		19,20,21,22,23,24,
     +		25,26,27,28,29,30,
     +		31,32,33,34,35,36,
     +		37,38,39,40,41,42,
     +		43,44,45,46,47,48,
     +		49,50,51,52,53,54,
     +		55,56,57,58,59,60,
     +		61,62,63,64,65,66,
     +		67,68,69,70,71,72,
     +		73,74,75,76,77,78,
     +		79,80,81,82,83,84,
     +		85,86,87,88,89,90,
     +		91,92,93,94,95,96,
     +		97,98,99, 100 /
	DATA GHYD/ 2.000000000000000, 8.000000000000000, 18.00000000000000,
     +		 32.00000000000000, 50.00000000000000, 72.00000000000000,
     +		 98.00000000000000, 128.0000000000000, 162.0000000000000,
     +		 200.0000000000000, 242.0000000000000, 288.0000000000000,
     +		 338.0000000000000, 392.0000000000000, 450.0000000000000,
     +		 512.0000000000000, 578.0000000000000, 648.0000000000000,
     +		 722.0000000000000, 800.0000000000000, 882.0000000000000,
     +		 968.0000000000000, 1058.000000000000, 1152.000000000000,
     +		 1250.000000000000, 1352.000000000000, 1458.000000000000,
     +		 1568.000000000000, 1682.000000000000, 1800.000000000000,
     +		 1922.000000000000, 2048.000000000000, 2178.000000000000,
     +		 2312.000000000000, 2450.000000000000, 2592.000000000000,
     +		 2738.000000000000, 2888.000000000000, 3042.000000000000,
     +		 3200.000000000000, 3362.000000000000, 3528.000000000000,
     +		 3698.000000000000, 3872.000000000000, 4050.000000000000,
     +		 4232.000000000000, 4418.000000000000, 4608.000000000000,
     +		 4802.000000000000, 5000.000000000000, 5202.000000000000,
     +		 5408.000000000000, 5618.000000000000, 5832.000000000000,
     +		 6050.000000000000, 6272.000000000000, 6498.000000000000,
     +		 6728.000000000000, 6962.000000000000, 7200.000000000000,
     +		 7442.000000000000, 7688.000000000000, 7938.000000000000,
     +		 8192.000000000000, 8450.000000000000, 8712.000000000000,
     +		 8978.000000000000, 9248.000000000000, 9522.000000000000,
     +		 9800.000000000000, 10082.00000000000, 10368.00000000000,
     +		 10658.00000000000, 10952.00000000000, 11250.00000000000,
     +		 11552.00000000000, 11858.00000000000, 12168.00000000000,
     +		 12482.00000000000, 12800.00000000000, 13122.00000000000,
     +		 13448.00000000000, 13778.00000000000, 14112.00000000000,
     +		 14450.00000000000, 14792.00000000000, 15138.00000000000,
     +		 15488.00000000000, 15842.00000000000, 16200.00000000000,
     +		 16562.00000000000, 16928.00000000000, 17298.00000000000,
     +		 17672.00000000000, 18050.00000000000, 18432.00000000000,
     +		 18818.00000000000, 19208.00000000000, 19602.00000000000,
     +		 20000.00000000000/
	DATA ENHYD /0.0000000000000000E+00,10.19085000000000,12.07804444444444,
     +		 12.73856250000000, 13.04428800000000, 13.21036111111111,
     +		 13.31049795918367, 13.37549062500000, 13.42004938271605,
     +		 13.45192200000000, 13.47550413223140, 13.49344027777778,
     +		 13.50739881656805, 13.51847448979592, 13.52740977777778,
     +		 13.53472265625000, 13.54078339100346, 13.54586234567901,
     +		 13.55016066481994, 13.55383050000000, 13.55698866213152,
     +		 13.55972603305785, 13.56211417769376, 13.56421006944444,
     +		 13.56605952000000, 13.56769970414201, 13.56916104252401,
     +		 13.57046862244898, 13.57164328180737, 13.57270244444444,
     +		 13.57366077003122, 13.57453066406250, 13.57532268135905,
     +		 13.57604584775087, 13.57670791836735, 13.57731558641975,
     +		 13.57787465303141, 13.57839016620499, 13.57886653517423,
     +		 13.57930762500000, 13.57971683521713, 13.58009716553288,
     +		 13.58045127095727, 13.58078150826446, 13.58108997530864,
     +		 13.58137854442344, 13.58164889090086, 13.58190251736111,
     +		 13.58214077467722, 13.58236488000000, 13.58257593233372,
     +		 13.58277492603550, 13.58296276254895, 13.58314026063100,
     +		 13.58330816528926, 13.58346715561225, 13.58361785164666,
     +		 13.58376082045184, 13.58389658144211, 13.58402561111111,
     +		 13.58414834721849, 13.58426519250780, 13.58437651801461,
     +		 13.58448266601563, 13.58458395266272, 13.58468067033976,
     +		 13.58477308977501, 13.58486146193772, 13.58494601974375,
     +		 13.58502697959184, 13.58510454274945, 13.58517889660494,
     +		 13.58525021580034, 13.58531866325785, 13.58538439111111,
     +		 13.58544754155125, 13.58550824759656, 13.58556663379356,
     +		 13.58562281685627, 13.58567690625000, 13.58572900472489,
     +		 13.58577920880428, 13.58582760923211, 13.58587429138322,
     +		 13.58591933564014, 13.58596281773932, 13.58600480908971,
     +		 13.58604537706612, 13.58608458527964, 13.58612249382716,
     +		 13.58615915952180, 13.58619463610586, 13.58622897444791,
     +		 13.58626222272522, 13.58629442659280, 13.58632562934028,
     +		 13.58635587203741, 13.58638519366930, 13.58641363126212,
     +		 13.58644122000000/
	DATA SHYD/100*0.0D0/
      DATA NHEL/1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,
     +        5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +        23,24,25,26,27,
     +		28,29,30,31,32,33,
     +		34,35,36,37,38,39,
     +		40,41,42,43,44,45,
     +		46,47,48,49,50,51,
     +		52,53,54,55,56,57,
     +		58,59,60,61,62,63,
     +		64,65,66,67,68,69,
     +		70,71,72,73,74,75,
     +		76,77,78,79,80,81/
      DATA GHEL/1.0D0,3.0D0,1.0D0,5.0D0,3.0D0,1.0D0,3.0D0,
     +          3.0D0,1.0D0,5.0D0,3.0D0,
     +         1.0D0,15.0D0,5.0D0,3.0D0,3.0D0,1.0D0,9.0D0,
     +         15.0D0,5.0D0,21.0D0,7.0D0,
     +         3.0D0,100.0D0,144.0D0,196.0D0,256.0D0,324.0D0,
     +         400.0D0,484.0D0,
     +         576.0D0,676.0D0,784.0D0,900.0D0,1024.0D0,1156.0D0,
     +         1296.0D0,1444.0D0,1600.0D0,1764.0D0,1936.0D0,
     +         2116.0D0,2304.0D0,2500.0D0,2704.0D0,3136.0D0,
     +		3136.000000000000,3364.000000000000,3600.000000000000,
     +		3844.000000000000,4096.000000000000,4356.000000000000,
     +		4624.000000000000,4900.000000000000,5184.000000000000,
     +		5476.000000000000,5776.000000000000,6084.000000000000,
     +		6400.000000000000,6724.000000000000,7056.000000000000,
     +		7396.000000000000,7744.000000000000,8100.000000000000,
     +		8464.000000000000,8836.000000000000,9216.000000000000,
     +		9604.000000000000,10000.00000000000,10404.00000000000,
     +		10816.00000000000,11236.00000000000,11664.00000000000,
     +		12100.00000000000,12544.00000000000,12996.00000000000,
     +		13456.00000000000,13924.00000000000,14400.00000000000,
     +		14884.00000000000,15376.00000000000,15876.00000000000,
     +		16384.00000000000,16900.00000000000,17424.00000000000,
     +		17956.00000000000,18496.00000000000,19044.00000000000,
     +		19600.00000000000,20164.00000000000,20736.00000000000,
     +		21316.00000000000,21904.00000000000,22500.00000000000,
     +		23104.00000000000,23716.00000000000,24336.00000000000,
     +		24964.00000000000,25600.00000000000,26244.00000000000/
      DATA ENHEL/0.0D0,19.819D0,20.615D0,20.964D0,
     +           20.964D0,20.964D0,21.218D0,
     +           22.718D0,22.920D0,23.007D0,23.007D0,
     +           23.007D0,23.073D0,23.074D0,
     +           23.087D0,23.593D0,23.673D0,23.707D0,
     +           23.736D0,23.736D0,23.737D0,
     +           23.737D0,23.742D0,24.028D0,24.201D0,
     +           24.304D0,24.371D0,24.417D0,
     +           24.449D0,24.473D0,24.491D0,24.506D0,
     +           24.517D0,24.526D0,24.534D0,
     +           24.540D0,24.545D0,24.549D0,24.553D0,
     +           24.556D0,24.559D0,24.562D0,
     +           24.564D0,24.566D0,24.568D0,24.570D0,
     +		24.57131951530612,24.57238228299643,24.57334055555556,
     +		24.57420759625390,24.57499462890625,24.57571120293848,
     +		24.57636548442907,24.57696448979592,24.57751427469136,
     +		24.57802008765522,24.57848649584488,24.57891748849441,
     +		24.57931656250000,24.57968679357525,24.58003089569161,
     +		24.58035127095727,24.58065005165289,24.58092913580247,
     +		24.58119021739130,24.58143481213219,24.58166427951389,
     +		24.58187984173261,24.58208260000000,24.58227354863514,
     +		24.58245358727811,24.58262353150587,24.58278412208505,
     +		24.58293603305785,24.58307987882653,24.58321622037550,
     +		24.58334557074911,24.58346839988509,24.58358513888889,
     +		24.58369618382155,24.58380189906348,24.58390262030738,
     +		24.58399865722656,24.58409029585799,24.58417780073462,
     +		24.58426141679661,24.58434137110727,24.58441787439614,
     +		24.58449112244898,24.58456129736163,24.58462856867284,
     +		24.58469309438919,24.58475502191381,24.58481448888889,
     +		24.58487162396122,24.58492654747850,24.58497937212360,
     +		24.58503020349303,24.58507914062500,24.58512627648224/
      DATA SHEL/0.375D0,0.622D0,0.622D0,0.842D0,
     +          0.842D0,0.842D0,0.842D0,0.747D0,
     +          0.747D0,0.912D0,0.912D0,0.912D0,
     +          0.993D0,0.993D0,0.912D0,0.810D0,
     +          0.810D0,0.937D0,0.995D0,0.995D0,
     +          1.000D0,1.000D0,0.937D0,0.949D0,
     +          0.958D0,75*1.000D0/
	DATA NCI/2,2,2,2,2,2,3,3,3,3,2,2,2,3,3,3,3,3,
     +		3,3,3,3,3,2,3,4,4,4,3,3,3,3,3,3,4,3,
     +		3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +		4,4,4,5,4,4,4,4,4,5,5,5,5,5,5,5,5,5,
     +		5,5,5,5,6,5,5,5,5,5,6,6,6,6,6,6,6,7,
     +		6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,
     +		8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,10,10,
     +		10,11,11,11,2,3,3,3,2,2/
	DATA GCI/1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,    
     +		5.0D0,1.0D0,3.0D0,5.0D0,3.0D0,    
     +		7.0D0,5.0D0,3.0D0,3.0D0,3.0D0,    
     +		5.0D0,7.0D0,3.0D0,1.0D0,3.0D0,    
     +		5.0D0,5.0D0,1.0D0,9.0D0,5.0D0,    
     +		1.0D0,3.0D0,5.0D0,5.0D0,7.0D0,    
     +		9.0D0,3.0D0,5.0D0,7.0D0,3.0D0,    
     +		3.0D0,3.0D0,5.0D0,3.0D0,1.0D0,    
     +		3.0D0,5.0D0,7.0D0,3.0D0,3.0D0,    
     +		1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,    
     +		5.0D0,5.0D0,7.0D0,9.0D0,3.0D0,    
     +		5.0D0,7.0D0,3.0D0,7.0D0,3.0D0,	 
     +		5.0D0,3.0D0,1.0D0,3.0D0,3.0D0,	 
     +		5.0D0,7.0D0,5.0D0,1.0D0,5.0D0,	 
     +		5.0D0,7.0D0,9.0D0,3.0D0,5.0D0,	 
     +		7.0D0,3.0D0,7.0D0,3.0D0,5.0D0,	 
     +		3.0D0,1.0D0,5.0D0,5.0D0,7.0D0,	 
     +		9.0D0,3.0D0,5.0D0,7.0D0,3.0D0,	 
     +		7.0D0,5.0D0,3.0D0,1.0D0,3.0D0,	 
     +		5.0D0,7.0D0,9.0D0,3.0D0,5.0D0,	 
     +		7.0D0,7.0D0,3.0D0,5.0D0,3.0D0,	 
     +		1.0D0,9.0D0,7.0D0,5.0D0,3.0D0,	 
     +		5.0D0,7.0D0,7.0D0,5.0D0,3.0D0,	 
     +		1.0D0,9.0D0,7.0D0,5.0D0,3.0D0,	 
     +		5.0D0,7.0D0,7.0D0,3.0D0,5.0D0,	 
     +		7.0D0,3.0D0,5.0D0,7.0D0,5.0D0,	 
     +		3.0D0,5.0D0,7.0D0,3.0D0,3.0D0/
	DATA ENCI/0.0D0,2.0333605D-03,5.3933649D-03,1.263870,2.684086,
     +		4.182672,7.480511,7.482891,7.487915,7.684888,
     +		7.946046,7.946620,7.946474,8.537387,8.640516,
     +		8.643146,8.647287,8.771255,8.846707,8.848247,
     +		8.850785,9.002712,9.171972,9.330682,9.631248,
     +		9.683908,9.685375,9.689256,9.695577,9.697620,
     +		9.701885,9.708156,9.708925,9.710041,9.712769,
     +		9.714380,9.761111,9.833419,9.834406,9.834934,
     +		9.940317,9.942698,9.946449,9.988707,10.05592,
     +		10.08144,10.08328,10.08553,10.13833,10.19809,
     +		10.35278,10.38514,10.38514,10.38514,10.39370,
     +		10.39456,10.39580,10.40021,10.40845,10.41874,
     +		10.42750,10.42990,10.42990,10.52043,10.52041,
     +		10.52041,10.53705,10.58840,10.61635,10.67973,
     +		10.70230,10.70328,10.70328,10.70878,10.70878,
     +		10.71184,10.71407,10.71854,10.72362,10.72523,
     +		10.72684,10.72684,10.86509,10.87426,10.87513,
     +		10.87513,10.87997,10.87997,10.88257,10.88533,
     +		10.88679,10.88964,10.89075,10.89075,10.88980,
     +		10.97789,10.97854,10.97854,10.98597,10.98597,
     +		10.98597,10.98808,10.98913,10.98994,10.98994,
     +		10.98994,11.04474,11.04474,11.04487,11.05280,
     +		11.05280,11.05280,11.05392,11.05429,11.05429,
     +		11.05429,11.09049,11.09049,11.09049,11.09843,
     +		11.09843,11.09843,11.09880,11.13129,11.13129,
     +		11.13129,11.15477,11.15477,11.15477,12.13544,
     +		12.83767,12.84024,12.84331,13.11772,14.86312/
	DATA NCII/2,2,2,2,2,2,2,2,2,2,3,3,3,2,3,3,2,2,4,4,4,3,3,3,
     +		4,4,2,2,4,4,5,5,5,3,3,5,5,5,5,6,3,3,3,3,3,3,6,6,
     +		6,6,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +		4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     +		5,5,5,5,5,5,6,6,6,6,6,6,6/
	DATA GCII/2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		6.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +		2.0D0,2.0D0,4.0D0,4.0D0,4.0D0,
     +		6.0D0,6.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +		4.0D0,6.0D0,6.0D0,8.0D0,2.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +		4.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		10.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		4.0D0,6.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0,8.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,4.0D0,6.0D0,4.0D0,
     +		6.0D0,8.0D0,10.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,4.0D0,6.0D0,6.0D0,
     +		4.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +		6.0D0,8.0D0,10.0D0,6.0D0,8.0D0,
     +		6.0D0,8.0D0,10.0D0,12.0D0,8.0D0,
     +		10.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0,4.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,6.0D0,4.0D0,
     +		2.0D0,6.0D0,8.0D0,4.0D0,6.0D0,
     +		8.0D0,10.0D0,6.0D0,8.0D0,10.0D0,
     +		12.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,6.0D0,
     +		4.0D0,2.0D0/
	DATA ENCII/0.0D0,7.9350658D-03,5.331397,5.334075,5.337658,
     +		9.290338,9.290624,11.96386,13.71590,13.72101,
     +		14.44900,16.33194,16.33332,17.60895,18.04607,
     +		18.04625,18.65519,18.65582,19.49478,20.14995,
     +		20.15068,20.70119,20.70413,20.70971,20.84491,
     +		20.84496,20.92025,20.92256,20.95094,20.95094,
     +		21.49265,21.73314,21.73405,22.09347,22.13075,
     +		22.13075,22.13075,22.18799,22.18799,22.47211,
     +		22.52747,22.52929,22.53239,22.53689,22.56844,
     +		22.57086,22.82136,22.82136,22.85996,22.85996,
     +		22.89870,23.11398,23.11600,23.11878,23.38108,
     +		23.38522,24.12408,24.27024,24.27201,24.27444,
     +		24.27787,24.37010,24.37079,24.37187,24.37315,
     +		24.60198,24.60332,24.65351,24.65617,24.65793,
     +		24.78982,24.79512,25.06741,25.07039,25.98117,
     +		25.98415,25.98986,26.58329,26.58615,26.62689,
     +		26.62867,26.63139,26.63554,26.75178,26.82771,
     +		26.82771,26.83016,26.89454,26.89578,27.22147,
     +		27.22329,27.22585,27.22930,27.29263,27.29263,
     +		27.29378,27.29509,27.35131,27.35294,27.37703,
     +		27.37957,27.38104,27.41188,27.41302,27.41395,
     +		27.41395,27.41395,27.41409,27.46301,27.46301,
     +		27.46810,27.46936,27.47200,27.47561,27.47330,
     +		27.47864,27.48713,27.49096,27.49330,27.49330,
     +		27.48854,27.49412,27.55688,27.56022,27.99752,
     +		27.99752,27.99752,28.25640,28.25640,28.61124,
     +		28.61124,28.61124,28.61124,28.64683,28.64683,
     +		28.64683,28.66803,26.43629,28.66875,28.66875,
     +		28.66875,28.66875,28.70253,28.70253,28.70253,
     +		28.70253,28.70515,28.70515,28.70515,28.70515,
     +		29.31561,29.31561,29.31561,29.31561,29.33557,
     +		29.33557,29.33557/
	DATA NCIII/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,
     +		3,4,4,4,4,3,4,4,4,4,4,4,4,4,3,3,3,4,3,3,3,3,3,3,
     +		3,3,3,3,3,3,5,3,3,3,3,5,5,5,5,3,5,5,5,5,5,5,5,5,
     +		5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,
     +		7,8,8,8,8,9,9,9,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +		4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,
     +		6,6,6,6,6,6,7,7,7,7,7,7/
	DATA GCIII/1.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +		1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +		3.0D0,1.0D0,3.0D0,1.0D0,3.0D0,
     +		5.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +		1.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +		3.0D0,5.0D0,7.0D0,5.0D0,7.0D0,
     +		9.0D0,3.0D0,7.0D0,3.0D0,5.0D0,
     +		7.0D0,5.0D0,3.0D0,1.0D0,3.0D0,
     +		5.0D0,5.0D0,5.0D0,5.0D0,7.0D0,
     +		9.0D0,3.0D0,5.0D0,7.0D0,3.0D0,
     +		5.0D0,3.0D0,1.0D0,7.0D0,3.0D0,
     +		1.0D0,3.0D0,5.0D0,1.0D0,3.0D0,
     +		5.0D0,7.0D0,7.0D0,9.0D0,11.0D0,
     +		9.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +		9.0D0,7.0D0,3.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,7.0D0,9.0D0,11.0D0,
     +		9.0D0,5.0D0,5.0D0,7.0D0,9.0D0,
     +		7.0D0,3.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,3.0D0,5.0D0,7.0D0,1.0D0,
     +		3.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,1.0D0,3.0D0,5.0D0,5.0D0,
     +		5.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,1.0D0,3.0D0,5.0D0,
     +		5.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +		5.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +		7.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,5.0D0,3.0D0,1.0D0,
     +		3.0D0,5.0D0,7.0D0,1.0D0,3.0D0,5.0D0/
	DATA ENCIII/0.0D0,6.486296,6.489148,6.496191,12.69008,
     +		17.03237,17.03602,17.04185,18.08638,22.62984,
     +		29.52845,30.64541,32.10371,32.19328,32.19396,
     +		32.19555,33.47080,33.45866,33.47146,34.27982,
     +		38.20770,38.21183,38.22034,38.36164,38.43612,
     +		38.64882,39.39549,39.39549,39.39611,39.64054,
     +		39.84380,39.84582,39.84874,39.91699,39.91782,
     +		39.91892,39.97328,40.01022,40.05026,40.05341,
     +		40.05822,40.19756,40.57121,40.86969,40.87231,
     +		40.87686,41.24874,41.30157,41.32848,41.33158,
     +		41.33611,41.85783,41.80309,41.86202,42.14028,
     +		42.16117,42.16444,42.16623,42.32471,42.55869,
     +		42.67342,42.67342,42.67342,42.78661,42.83001,
     +		42.83001,42.83001,42.96405,42.96405,42.96416,
     +		42.96405,42.98029,42.98736,43.03527,43.03550,
     +		43.03579,43.25349,43.98952,44.27370,44.39248,
     +		44.39248,44.39248,44.46592,44.46592,44.46600,
     +		44.47219,44.47673,44.48596,44.48596,44.48596,
     +		44.52591,45.07626,45.24178,45.32720,45.32720,
     +		45.32720,45.38200,45.86543,45.92891,45.92891,
     +		45.92891,46.33929,46.33929,46.33929,46.69749,
     +		46.69749,46.69749,47.25143,47.35238,47.35238,
     +		47.35722,47.64920,47.64920,47.65379,47.81342,
     +		47.83558,48.06245,48.06245,48.06245,48.16114,
     +		48.16114,48.16114,48.20208,50.51542,50.55803,
     +		50.55803,50.55803,50.69428,50.69428,50.69428,
     +		50.77264,50.79460,50.90022,50.90022,50.90022,
     +		50.93829,50.93829,50.93829,52.24497,52.24497,
     +		52.24497,52.31775,52.31775,52.31775,52.43107,
     +		52.43107,52.43107,52.45302,52.45302,52.45302,
     +		53.23251,53.23251,53.23251,53.27802,53.27802,
     +		53.27802/
	DATA NCIV/2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,
     +		6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,
     +		8,8,8,8,8,8,8/
	DATA GCIV/2.0D0,2.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +		2.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +		6.0D0,8.0D0,8.0D0,10.0D0,2.0D0,
     +		2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +		8.0D0,8.0D0,10.0D0,10.0D0,12.0D0,
     +		2.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +		6.0D0,8.0D0,8.0D0,10.0D0,10.0D0,
     +		12.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		8.0D0,10.0D0,12.0D0,14.0D0,16.0D0/
	DATA ENCIV/0.0D0,7.995100,8.008378,37.54872,39.68134,
     +		39.68525,40.28040,40.28173,49.76113,50.62434,
     +		50.62599,50.87540,50.87595,50.88784,50.88784,
     +		55.21889,55.65134,55.65221,55.77947,55.77947,
     +		55.78577,55.78578,55.78703,55.78703,58.12002,
     +		58.36774,58.36774,58.44275,58.44275,58.44709,
     +		58.44709,58.44764,58.44764,58.44770,58.44770,
     +		59.84267,60.00038,60.00038,60.04725,60.04725,
     +		60.05156,60.05156,60.05191,60.05191,60.05194,
     +		60.05194,61.05946,61.05946,61.09294,61.09294,
     +		61.09319,61.09319,61.09319,61.09319,61.09319/
	DATA NCV/1,2,2,2,2,2,3,3,3,3,4,5,6,7,8/
	DATA GCV/1.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +		3.0D0,3.0D0,5.0D0,7.0D0,3.0D0,
     +		3.0D0,3.0D0,3.0D0,3.0D0,3.0D0/
	DATA ENCV/0.0D0,298.9618,304.4046,304.4030,304.4199,
     +		307.8855,354.2645,354.2645,354.2645,354.5177,
     +		370.9247,378.5349,382.6710,385.1917,386.6807/
	DATA NNI/2,2,2,2,2,3,3,3,3,3,2,2,2,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,3,3,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,4,4,4,4,4,4,4,4,
     +		4,4,4,4,4,4,4,4,4,3,3,3,3,6,6,6,6,6,5,5,5,5,5,5,
     +		5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,6,6,6,6,6,6,6,6,
     +		6,6,6,6,6,6,6,6,6,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,
     +		7,7,7,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,10,10,10,
     +		10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,11,11,11,11,11,10,
     +		10,10,10,10,10,10,10,10,10,10,10,10,12,12,12,12,12,
     +		11,11,11,11,11,11,11,11,11,11,11,11,11,13,13,12,12,
     +		12,12,12,12,12/
	DATA GNI/4.0D0,6.0D0,4.0D0,4.0D0,2.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		6.0D0,4.0D0,2.0D0,2.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +		6.0D0,4.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,10.0D0,6.0D0,
     +		8.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,4.0D0,6.0D0,
     +		2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		2.0D0,4.0D0,6.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,4.0D0,
     +		6.0D0,8.0D0,10.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,4.0D0,2.0D0,6.0D0,
     +		8.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +		6.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		4.0D0,6.0D0,8.0D0,10.0D0,4.0D0,
     +		2.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,2.0D0,4.0D0,6.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,4.0D0,6.0D0,8.0D0,
     +		10.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		4.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,6.0D0,8.0D0,
     +		4.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +		8.0D0,4.0D0,2.0D0,6.0D0,8.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +		2.0D0,6.0D0,8.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0,8.0D0,2.0D0,4.0D0,6.0D0,
     +		8.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,4.0D0,2.0D0,6.0D0,8.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +		6.0D0,4.0D0,6.0D0/
	DATA ENNI/0.0D0,2.383371,2.384363,3.575739,3.575739,
     +		10.32619,10.33038,10.33617,10.67904,10.69042,
     +		10.92429,10.92973,10.93217,11.60284,11.75037,
     +		11.75317,11.75780,11.76412,11.83769,11.83997,
     +		11.84472,11.99580,12.00032,12.00975,12.12207,
     +		12.12649,12.35701,12.35614,12.84713,12.85333,
     +		12.86185,12.91211,12.92268,12.97078,12.97568,
     +		12.97693,12.97929,12.98350,12.98958,12.99502,
     +		13.00392,13.00161,13.00483,13.00074,13.01686,
     +		13.01822,13.01983,13.02095,13.03344,13.03636,
     +		13.20179,13.23674,13.23917,13.24364,13.25041,
     +		13.26429,13.26623,13.27127,13.32189,13.61527,
     +		13.62076,13.62945,13.64202,13.65185,13.66270,
     +		13.66493,13.66914,13.67609,13.66580,13.67249,
     +		13.67410,13.68043,13.66588,13.66872,13.67695,
     +		13.68464,13.67869,13.68191,13.68836,13.69398,
     +		13.69673,13.70310,13.70607,13.92292,13.92614,
     +		13.95653,13.96207,13.97100,13.97749,13.98841,
     +		13.97948,13.98097,13.98543,13.99324,13.98568,
     +		13.98754,13.98803,13.99674,13.98865,13.98865,
     +		13.98865,13.99696,13.99237,13.99473,13.99944,
     +		14.00155,14.00384,14.13620,14.14326,14.15244,
     +		14.15045,14.15455,14.15417,14.15417,14.15417,
     +		14.15417,14.15690,14.15690,14.15690,14.16508,
     +		14.15827,14.16025,14.15864,14.16843,14.16313,
     +		14.17035,14.16645,14.16645,14.16831,14.23464,
     +		14.24468,14.25113,14.25212,14.25212,14.25683,
     +		14.25683,14.25683,14.25683,14.25882,14.25882,
     +		14.26043,14.26043,14.26545,14.27073,14.27109,
     +		14.27109,14.27109,14.36247,14.36247,14.31821,
     +		14.31821,14.31821,14.32329,14.32329,14.32329,
     +		14.32329,14.32403,14.32403,14.32465,14.32465,
     +		14.33234,14.33544,14.33494,14.33494,14.33494,
     +		14.36272,14.36272,14.36433,14.36433,14.36433,
     +		14.36830,14.36830,14.36830,14.36830,14.36854,
     +		14.36854,14.37016,14.37016,14.37896,14.38119,
     +		14.38107,14.38107,14.38107,14.39557,14.39557,
     +		14.39768,14.39768,14.39768,14.40152,14.40152,
     +		14.40202,14.40202,14.40264,14.40264,14.40264,
     +		14.40264,14.41206,14.41206,14.41442,14.41442,
     +		14.41442,14.42012,14.42012,14.42099,14.42099,
     +		14.42099,14.42583,14.42583,14.42682,14.42682,
     +		14.42781,14.42781,14.42781,14.42781,14.43636,
     +		14.43636,14.43698,14.43698,14.43698,14.46253,
     +		14.44021,14.44455,14.44455,14.45434,14.45434,
     +		14.45434,14.45980,14.45980/
	DATA NNII/2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,2,3,3,3,3,2,3,
     +		3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,3,4,4,4,
     +		4,4,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +		4,4,4,4,4,4,4,4,4,4,3,3,3,5,5,5,5,5,5,5,5,5,5,5,
     +		5,5,5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/
	DATA GNII/1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +		5.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +		3.0D0,1.0D0,5.0D0,1.0D0,3.0D0,
     +		5.0D0,3.0D0,3.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,3.0D0,3.0D0,1.0D0,
     +		3.0D0,5.0D0,5.0D0,1.0D0,5.0D0,
     +		7.0D0,9.0D0,5.0D0,3.0D0,5.0D0,
     +		7.0D0,5.0D0,3.0D0,1.0D0,7.0D0,
     +		3.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +		3.0D0,3.0D0,5.0D0,7.0D0,1.0D0,
     +		3.0D0,5.0D0,3.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,1.0D0,5.0D0,7.0D0,
     +		9.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +		5.0D0,3.0D0,1.0D0,7.0D0,5.0D0,
     +		7.0D0,9.0D0,7.0D0,7.0D0,9.0D0,
     +		11.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +		3.0D0,5.0D0,1.0D0,3.0D0,5.0D0,
     +		1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,5.0D0,7.0D0,9.0D0,
     +		7.0D0,7.0D0,9.0D0,11.0D0,9.0D0,
     +		1.0D0,3.0D0,5.0D0,7.0D0,9.0D0,
     +		3.0D0,5.0D0,7.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,9.0D0,11.0D0,7.0D0,
     +		5.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +		7.0D0,9.0D0/
	DATA ENNII/0.0D0,6.0876831D-03,1.6279284D-02,1.898923,4.052723,
     +		5.848106,11.43604,11.43781,11.43801,13.54146,
     +		13.54146,13.54228,17.87734,18.46259,18.46651,
     +		18.48341,18.49722,19.23384,20.40944,20.64636,
     +		20.65389,20.66582,20.67651,20.94027,21.14861,
     +		21.15298,21.16022,21.59986,22.10340,23.12481,
     +		23.13218,23.14229,23.19670,23.23962,23.24260,
     +		23.24636,23.41565,23.42207,23.42555,23.47490,
     +		23.57225,24.36823,24.37465,24.38944,24.53166,
     +		25.06612,25.13369,25.14001,25.15193,25.18946,
     +		25.19245,25.20124,25.23510,25.46049,25.53877,
     +		25.54572,25.55447,25.58160,25.99668,26.00464,
     +		26.01527,26.02787,26.06667,26.06994,26.07548,
     +		26.12440,26.13011,26.13327,26.16475,26.16510,
     +		26.16800,26.16849,26.17391,26.19663,26.19758,
     +		26.20937,26.20252,26.21087,26.21191,26.21252,
     +		26.22134,26.22182,26.25393,26.25770,26.26368,
     +		26.55921,26.56489,26.58065,26.63554,27.36569,
     +		27.36569,27.36569,27.40948,27.40948,27.40999,
     +		27.41783,27.42901,27.42963,27.43824,27.43947,
     +		27.77609,27.77805,27.78169,27.78704,27.79372,
     +		28.01910,28.02209,28.02755,28.54429,30.17253,
     +		30.17448,30.17763,30.18179,30.18682,30.34387,
     +		30.34864,30.35188,30.41607,30.41652,30.41750,
     +		30.41894,30.42068/
	DATA NNIII/2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,3,3,3,3,3,3,3,
     +		3,4,3,3,3,3,3,3,4,4,3,3,3,3,4,4,4,4,3,3,3,3,3,3,
     +		3,3,3,3,3,5,3,3,3,3,3,3,3,5,5,3,3,5,5,5,5,6,6,6,
     +		6,6,6,4,4,4,3,3,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,
     +		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +		4,4,4,4,4,4,4,3,3,5,5,5,5/
	DATA GNIII/2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		6.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +		4.0D0,6.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,2.0D0,4.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,4.0D0,6.0D0,
     +		6.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,10.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +		6.0D0,6.0D0,4.0D0,2.0D0,6.0D0,
     +		8.0D0,4.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0,8.0D0,8.0D0,10.0D0,4.0D0,
     +		6.0D0,6.0D0,8.0D0,8.0D0,10.0D0,
     +		2.0D0,4.0D0,6.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		8.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +		4.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +		6.0D0,8.0D0,10.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,6.0D0,
     +		4.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +		6.0D0,8.0D0,10.0D0,6.0D0,8.0D0,
     +		6.0D0,8.0D0,10.0D0,12.0D0,8.0D0,
     +		10.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0,4.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0/
	DATA ENNIII/0.0D0,2.1635452D-02,7.180255,7.098413,7.108480,
     +		12.52548,12.52643,16.24252,18.08651,18.10019,
     +		23.16076,25.17799,25.18006,27.43827,28.56680,
     +		28.56730,30.45896,30.46342,33.13367,33.13441,
     +		35.65022,35.65797,35.67233,36.84229,36.85629,
     +		38.44641,38.32793,38.33453,38.39367,38.39807,
     +		38.40689,38.41771,38.64517,38.64825,38.95919,
     +		39.34056,39.34595,39.35325,39.39646,39.40031,
     +		39.71098,39.71098,39.79651,39.80747,40.55027,
     +		40.94474,40.94909,40.95552,40.96437,41.26192,
     +		41.26358,41.26631,41.26982,41.37555,41.47835,
     +		41.48166,41.68555,41.69232,41.69667,42.12335,
     +		42.13715,42.39634,42.39655,42.48893,42.49769,
     +		42.49625,42.49625,42.54757,42.54757,43.95493,
     +		43.95493,44.00932,44.00932,44.04135,44.04135,
     +		45.69180,45.69957,45.71402,46.28896,46.29317,
     +		46.46321,46.47039,46.71232,46.71811,46.72555,
     +		46.73671,46.81577,46.81788,46.85206,46.86286,
     +		46.92110,47.02857,47.03412,47.04068,47.61238,
     +		47.61238,47.61845,47.62763,47.75000,47.75000,
     +		47.77108,47.77108,49.01428,47.77802,47.88887,
     +		47.88887,47.88887,47.97657,47.97913,47.98245,
     +		47.98245,47.98363,47.98760,48.07270,48.08297,
     +		48.11119,48.11662,48.12305,48.13089,48.12993,
     +		48.14229,48.14024,48.14488,48.15087,48.15427,
     +		48.15307,48.16119,49.16950,49.17073,50.71214,
     +		50.71214,50.71214,50.71214/
	DATA NNIV/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,3,3,3,3,4,3,3,3,3,3,4,4,4,3,3,3,3,4,4,4,4,3,
     +		3,3,4,4,4,4,3,4,5,5,5,5,5,5,5,6,6,6,4,4,4,4,5,5,4/
	DATA GNIV/1.0D0,1.0D0,3.0D0,7.0D0,3.0D0,
     +		1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +		3.0D0,1.0D0,1.0D0,3.0D0,5.0D0,
     +		3.0D0,5.0D0,7.0D0,5.0D0,1.0D0,
     +		3.0D0,5.0D0,3.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,3.0D0,1.0D0,3.0D0,
     +		5.0D0,5.0D0,5.0D0,5.0D0,7.0D0,
     +		9.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,7.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,5.0D0,3.0D0,1.0D0,
     +		5.0D0,5.0D0,7.0D0,9.0D0,3.0D0,
     +		7.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +		7.0D0,9.0D0,11.0D0,3.0D0,5.0D0,
     +		7.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +		3.0D0,5.0D0,7.0D0/
	DATA ENNIV/0.0D0,8.323934,8.331770,8.349648,16.20427,
     +		21.75491,21.76399,21.77946,23.41898,29.18244,
     +		46.76804,50.15470,50.32483,50.32679,50.33118,
     +		52.06988,52.07031,52.07132,53.20933,57.68086,
     +		57.69048,57.71067,58.64906,59.62210,60.05779,
     +		60.05779,60.07403,60.44809,61.27855,61.27855,
     +		61.29070,61.78379,61.95650,61.97423,61.97423,
     +		61.97423,62.44215,62.44215,62.44215,62.67301,
     +		62.67685,62.68218,62.77282,62.86333,63.40415,
     +		63.40415,63.40415,63.41109,63.41767,63.41767,
     +		63.80760,64.05482,64.05569,64.05706,64.39976,
     +		64.70402,68.21900,68.53058,68.53058,68.53058,
     +		68.73986,68.73986,68.73986,71.28416,71.28416,
     +		71.28416,73.28070,73.60580,73.60580,73.61063,
     +		78.63129,78.63129,78.63129/
	DATA NNV/2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,
     +		6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8/
	DATA GNV/2.0D0,2.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +		8.0D0,10.0D0,12.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +		8.0D0,10.0D0,12.0D0,14.0D0,2.0D0,
     +		2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +		8.0D0,8.0D0,10.0D0,12.0D0,14.0D0,16.0D0/
	DATA ENNV/0.0D0,9.976473,10.00851,56.55396,59.23740,
     +		59.24660,60.05890,60.06188,75.17694,76.26962,
     +		76.26962,76.61120,76.61120,83.55153,84.09893,
     +		84.09893,84.27598,84.27598,88.02306,88.33514,
     +		88.33514,88.43854,88.43742,88.44214,88.44214,
     +		88.44313,88.44313,88.44313,90.68689,90.88043,
     +		90.88043,90.94527,90.94527,90.94912,90.94912,
     +		90.94974,90.94974,90.94974,90.94974,92.40136,
     +		92.53167,92.53167,92.57358,92.57358,92.57618,
     +		92.57618,92.57668,92.57668,92.57668,92.57668,
     +		92.57668/
	DATA NNVI/1,2,2,2,2,2,3,4/
	DATA GNVI/1.0D0,3.0D0,1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,3.0D0/
	DATA ENNVI/0.0D0,419.8009,426.2953,426.2965,426.3325,
     +		425.7398,497.9737,521.5830/
	DATA NOI/2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,3,3,3,3,3,3,3,3,3,
     +		4,4,4,4,4,4,3,3,3,5,5,3,4,4,4,4,4,4,4,4,5,5,5,6,
     +		6,5,5,5,5,5,5,5,5,6,6,6,7,7,6,6,6,6,6,6,6,6,8,8,
     +		7,7,7,7,7,7,7,7,9,9,8,8,8,8,8,8,8,8,10,10,9,9,9,9,
     +		9,9,9,9,11,11,10,10,10,10,10,10,10,10,3,3,3,3,3,3,
     +		3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,4,4,4,2,2,2,3,
     +		3,3,3,3,5,4,4,4,4,4,4,4,4,4,4,4,3,6,5,5,5,5,5,5,5,
     +		5,5,5,7,6,6,6,2/
	DATA GOI/5.0D0,3.0D0,1.0D0,5.0D0,1.0D0,
     +		5.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +		5.0D0,3.0D0,1.0D0,5.0D0,3.0D0,
     +		9.0D0,7.0D0,5.0D0,5.0D0,3.0D0,
     +		1.0D0,7.0D0,5.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,5.0D0,3.0D0,1.0D0,
     +		7.0D0,5.0D0,3.0D0,5.0D0,3.0D0,
     +		5.0D0,9.0D0,7.0D0,5.0D0,3.0D0,
     +		1.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +		3.0D0,1.0D0,5.0D0,3.0D0,9.0D0,
     +		7.0D0,5.0D0,3.0D0,1.0D0,7.0D0,
     +		5.0D0,3.0D0,5.0D0,3.0D0,1.0D0,
     +		5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +		5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +		5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +		5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +		5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +		7.0D0,5.0D0,3.0D0,9.0D0,7.0D0,
     +		5.0D0,5.0D0,3.0D0,1.0D0,7.0D0,
     +		3.0D0,5.0D0,5.0D0,5.0D0,3.0D0,
     +		1.0D0,9.0D0,7.0D0,5.0D0,9.0D0,
     +		11.0D0,9.0D0,7.0D0,7.0D0,7.0D0,
     +		5.0D0,3.0D0,5.0D0,3.0D0,1.0D0,
     +		7.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +		5.0D0,9.0D0,7.0D0,5.0D0,9.0D0,
     +		11.0D0,9.0D0,7.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,1.0D0,5.0D0,9.0D0,
     +		7.0D0,5.0D0,9.0D0,11.0D0,9.0D0,
     +		7.0D0,5.0D0,3.0D0,1.0D0,5.0D0,
     +		5.0D0,3.0D0,1.0D0,3.0D0/
	DATA ENOI/0.0D0,01.9651687D-02,2.8082693D-02,1.967363,0.4206081,
     +		9.146132,9.521420,10.74028,10.74053,10.74098,
     +		10.98893,10.98886,10.98895,11.83768,11.93056,
     +		12.07869,12.07870,12.07870,12.07872,12.07872,
     +		12.07872,12.08711,12.08711,12.08711,12.28604,
     +		12.28612,12.28627,12.35891,12.35891,12.35891,
     +		12.53927,12.54078,12.54176,12.66092,12.69755,
     +		12.72854,12.75377,12.75377,12.75377,12.75377,
     +		12.75377,12.75911,12.75911,12.75911,12.87829,
     +		12.87829,12.87829,13.02082,13.03891,13.06624,
     +		13.06624,13.06624,13.06624,13.06624,13.06913,
     +		13.06913,13.06913,13.13145,13.13145,13.13145,
     +		13.21004,13.22030,13.23559,13.23559,13.23559,
     +		13.23559,13.23559,13.23740,13.23740,13.23740,
     +		13.32166,13.32807,13.33749,13.33749,13.33749,
     +		13.33749,13.33749,13.33869,13.33869,13.33869,
     +		13.39308,13.39756,13.40353,385.3597,13.40353,
     +		13.40353,13.40353,13.40488,13.40488,13.40488,
     +		13.44262,13.44449,13.44872,13.44872,13.44872,
     +		13.44872,13.44872,13.44966,13.44966,13.44966,
     +		13.47577,13.47812,13.48112,13.48112,13.48112,
     +		13.48112,13.48112,13.48148,13.48148,13.48148,
     +		14.04685,14.04687,14.04730,14.09888,14.09975,
     +		14.10046,14.12320,14.12450,14.12526,14.13382,
     +		14.37218,14.46048,15.22525,15.28698,15.29424,
     +		15.29817,15.40062,15.40062,15.40062,15.40372,
     +		15.40390,15.40622,15.40550,15.41465,15.59420,
     +		15.59514,15.59577,15.65520,15.66431,15.66970,
     +		15.78109,15.78181,15.78222,15.82895,15.94391,
     +		16.01073,16.07676,16.07676,16.07676,16.07836,
     +		16.07844,16.08080,16.08005,16.08545,16.11433,
     +		16.11550,16.11614,16.23505,16.35702,16.35702,
     +		16.35702,16.35702,16.39057,16.39063,16.39308,
     +		16.39308,16.40451,16.40451,16.40451,16.54127,
     +		16.56668,16.56668,16.56668,23.53702/
	DATA NOII/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,2,3,3,3,3,3,3,3,3,
     +		3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,3,3,3,3,3,3,3,4,4,4,4,4,3,4,4,4,4,4,4,4,4,3,
     +		3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,
     +		4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,
     +		5,5,5,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     +		5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3,3,3,4,4,4,4,
     +		4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,5,5,3,3,3,3,3,4/
	DATA GOII/4.0D0,6.0D0,4.0D0,4.0D0,2.0D0,
     +		6.0D0,4.0D0,2.0D0,6.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		2.0D0,2.0D0,2.0D0,4.0D0,6.0D0,
     +		8.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,4.0D0,6.0D0,4.0D0,4.0D0,
     +		2.0D0,2.0D0,4.0D0,2.0D0,6.0D0,
     +		8.0D0,6.0D0,4.0D0,4.0D0,6.0D0,
     +		8.0D0,10.0D0,6.0D0,4.0D0,2.0D0,
     +		2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		8.0D0,6.0D0,8.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,2.0D0,4.0D0,8.0D0,6.0D0,
     +		10.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,8.0D0,10.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +		6.0D0,4.0D0,2.0D0,4.0D0,2.0D0,
     +		6.0D0,8.0D0,2.0D0,6.0D0,4.0D0,
     +		8.0D0,5.80D0,4.0D0,2.0D0,6.0D0,
     +		8.0D0,10.0D0,12.0D0,8.0D0,10.0D0,
     +		4.0D0,6.0D0,4.0D0,6.0D0,8.0D0,
     +		10.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,2.0D0,
     +		4.0D0,6.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,8.0D0,6.0D0,4.0D0,
     +		2.0D0,6.0D0,8.0D0,8.0D0,6.0D0,
     +		4.0D0,2.0D0,6.0D0,8.0D0,10.0D0,
     +		12.0D0,8.0D0,10.0D0,4.0D0,6.0D0,
     +		4.0D0,6.0D0,8.0D0,10.0D0,6.0D0,
     +		8.0D0,4.0D0,6.0D0,8.0D0,6.0D0,
     +		8.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		8.0D0,10.0D0,6.0D0,8.0D0,6.0D0,
     +		4.0D0,2.0D0,4.0D0,6.0D0,10.0D0,
     +		12.0D0,2.0D0,4.0D0,4.0D0,6.20D0,
     +		10.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0/
	DATA ENOII/0.0D0,3.323850,3.326454,5.017305,5.017491,
     +		14.85813,14.87838,14.88860,20.58005,20.57736,
     +		22.96648,22.97954,23.001876,23.41940,23.44172,
     +		24.26523,25.28586,25.63160,25.63849,25.64984,
     +		25.66529,25.66142,25.66154,25.83188,25.83760,
     +		25.84900,26.22564,26.24928,26.30498,26.35845,
     +		26.37943,26.55392,26.56133,28.12621,28.35835,
     +		28.36128,28.51330,28.51270,28.67733,28.68403,
     +		28.69369,28.70637,28.82200,28.83108,28.83932,
     +		28.82414,28.82992,28.85285,28.85711,28.85729,
     +		28.85808,28.86334,28.88355,28.94193,28.95606,
     +		29.06249,29.06893,29.58618,29.59923,29.61924,
     +		29.79726,29.82051,30.42546,30.47162,30.47763,
     +		30.48836,30.50400,30.74951,30.77135,30.80112,
     +		30.81214,31.02747,31.02747,31.14773,31.14812,
     +		31.31967,31.31982,31.37404,31.37430,31.46620,
     +		31.46649,31.55199,31.55199,31.55199,31.56553,
     +		31.61407,31.61407,31.61407,31.61407,31.61407,
     +		31.62925,31.63375,31.63644,31.63766,31.65117,
     +		31.65364,31.67396,31.69345,31.70178,31.71699,
     +		31.70200,31.71709,31.72948,31.72935,31.70999,
     +		31.71043,31.71889,31.73747,31.71911,31.73823,
     +		31.72081,31.72752,31.75062,31.75112,31.75553,
     +		31.75715,31.75586,31.75803,31.95026,31.96318,
     +		31.98375,32.03889,32.06284,32.14771,32.14780,
     +		32.35511,32.35511,32.36540,32.38251,32.39264,
     +		32.39264,32.40412,32.44667,32.46798,32.88345,
     +		32.88345,32.88345,32.88345,32.90963,32.91418,
     +		32.91418,32.92780,32.92780,32.93536,32.94354,
     +		32.95061,32.96264,32.93858,32.94181,32.95049,
     +		32.97082,32.95073,32.97146,32.96227,32.96227,
     +		32.97119,32.97528,32.97826,32.97999,32.97863,
     +		32.97999,33.19875,33.19968,33.20123,34.06365,
     +		34.06901,34.08607,34.08607,34.17174,34.17174,
     +		34.20029,34.20029,34.20504,34.20504,34.21390,
     +		34.21390,34.21960,34.22819,34.22819,34.23350,
     +		34.23350,34.25269,34.25269,34.48530,34.48530,
     +		36.19083,36.18759,36.19109,36.19123,36.19131,37.05294/
	DATA NOIII/2,2,2,2,2,2,2,3,3,2,2,2,2,2,2,3,3,3,3,2,2,2,3,3,
     +		3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		2,3,3,3,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,3,3,3,4,4,
     +		4,4,4,3,3,3,3,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,5,5,
     +		5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		5,5,5,5,5,5,5,5,5,3,3,3,6,6,6,6,7,3,3,4,4,4,3,4,
     +		4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,5/
	DATA GOIII/1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +		5.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +		3.0D0,1.0D0,5.0D0,3.0D0,3.0D0,
     +		1.0D0,3.0D0,5.0D0,3.0D0,5.0D0,
     +		3.0D0,1.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,3.0D0,5.0D0,1.0D0,3.0D0,
     +		5.0D0,5.0D0,1.0D0,5.0D0,7.0D0,
     +		9.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +		5.0D0,3.0D0,1.0D0,7.0D0,3.0D0,
     +		3.0D0,5.0D0,7.0D0,1.0D0,1.0D0,
     +		3.0D0,5.0D0,1.0D0,3.0D0,5.0D0,
     +		3.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +		7.0D0,9.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +		1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +		3.0D0,5.0D0,7.0D0,5.0D0,5.0D0,
     +		7.0D0,9.0D0,5.0D0,5.0D0,3.0D0,
     +		1.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		3.0D0,1.0D0,7.0D0,3.0D0,1.0D0,
     +		3.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,3.0D0,5.0D0,7.0D0,9.0D0,
     +		11.0D0,1.0D0,3.0D0,5.0D0,7.0D0,
     +		9.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +		3.0D0,1.0D0,5.0D0,7.0D0,9.0D0,
     +		5.0D0,7.0D0,9.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,7.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +		7.0D0,7.0D0,7.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,3.0D0,3.0D0,1.0D0,
     +		3.0D0,5.0D0,7.0D0,9.0D0,3.0D0,
     +		5.0D0,7.0D0,3.0D0,5.0D0,7.0D0,
     +		7.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +		9.0D0,3.0D0,5.0D0,7.0D0,1.0D0,
     +		3.0D0,5.0D0,3.0D0/
	DATA ENOIII/0.0D0,1.4059945D-02,3.8038719D-02,2.513308,5.354124,
     +		7.477820,14.88140,14.88477,14.88550,17.65325,
     +		17.65339,17.65514,23.19140,24.43587,26.09378,
     +		33.13600,33.15068,33.18253,33.85794,35.18196,
     +		35.20895,35.22094,36.07438,36.43500,36.45190,
     +		36.47919,36.89279,36.98353,37.22392,37.23410,
     +		37.25028,38.01204,38.90675,40.22861,40.25288,
     +		40.27497,40.26230,40.57149,40.57759,40.58673,
     +		40.84922,40.86335,40.87098,41.14086,41.25951,
     +		41.97723,41.99266,42.14902,42.56451,43.39812,
     +		43.41013,43.43237,44.22956,44.24270,44.27655,
     +		44.46952,45.03978,45.31862,45.32294,45.33144,
     +		45.34384,45.35962,45.34443,45.43903,45.45230,
     +		45.47797,45.62070,45.69189,45.69899,45.71153,
     +		45.91510,45.92614,45.93959,45.98626,46.25228,
     +		46.44183,45.21283,46.46955,46.62690,46.78899,
     +		46.78899,46.78899,46.82767,46.91713,46.91867,
     +		46.92080,47.01923,47.02679,47.03461,47.20199,
     +		47.20199,47.20199,47.21141,47.24910,48.62968,
     +		48.62968,48.62968,48.69874,48.86141,48.86587,
     +		48.87442,48.91428,48.91908,48.92621,48.93560,
     +		48.94701,49.36293,49.36248,49.36198,49.36323,
     +		49.37332,49.40500,49.41368,49.41845,49.63815,
     +		49.65178,49.65844,49.76514,49.77709,49.79367,
     +		49.78386,49.78386,49.78386,49.81572,49.78386,
     +		49.78386,49.78386,50.01249,50.03133,50.31391,
     +		50.31750,50.32357,51.41365,51.47638,51.47638,
     +		51.47638,52.44297,52.69355,52.85969,53.12613,
     +		53.14089,53.16110,53.31682,54.18348,54.33549,
     +		54.33549,54.34320,54.35460,54.36977,54.46407,
     +		54.47044,54.48261,54.88958,54.88958,54.88958,
     +		55.81414,55.82281,55.82951,56.14741,56.14741,
     +		56.14741,56.31095,56.31095,56.31095,56.73994,
     +		56.73994,56.73994,58.73808/
	DATA NOIV/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,
     +		3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,4,4,3,3,3,3,3,3,5,3,3,3,3,5,5,5,5,3,4,4,4,3,
     +		3,4,4,6,6,4,4,3,3,3,3,3,3,3,4,4,7,7,4,4,4,4,4,4,
     +		4,4,4,4,4,4,4,4,4,4,3,8,8,4,4,3,3,3,3,3,3,3,3,3,
     +		3,3,3,3,3,5,5,3,3,5,5,3,3,5,5,5,5,5,5,5,5,5,5,5,
     +		3,3,3,3,3,3,3,3,3,6,6,6,6,4,4,3,4,4,7,7,7,7/
	DATA GOIV/2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		6.0D0,4.0D0,2.0D0,2.0D0,6.0D0,
     +		4.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +		2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,4.0D0,6.0D0,
     +		2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		10.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		4.0D0,6.0D0,6.0D0,4.0D0,2.0D0,
     +		4.0D0,6.0D0,6.0D0,8.0D0,4.0D0,
     +		2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +		6.0D0,8.0D0,2.0D0,2.0D0,4.0D0,
     +		6.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,4.0D0,
     +		6.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +		2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		6.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +		6.0D0,8.0D0,4.0D0,6.0D0,6.0D0,
     +		8.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +		2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,2.0D0,4.0D0,6.0D0,
     +		6.0D0,4.0D0,4.0D0,6.0D0,8.0D0,
     +		2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +		6.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +		6.0D0,8.0D0,4.0D0,2.0D0,6.0D0,
     +		4.0D0,2.0D0,4.0D0,6.0D0,6.0D0,
     +		8.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +		6.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +		4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,8.0D0/
	DATA ENOIV/0.0D0,4.7920357D-02,8.824909,8.841201,8.864076,
     +		15.73825,15.73998,20.37910,22.37705,22.40721,
     +		28.67474,31.63571,31.63934,35.83378,35.83476,
     +		44.33902,48.37428,48.38508,54.37857,54.39532,
     +		54.42593,56.14158,56.17444,57.92984,57.94415,
     +		58.03452,58.04428,58.06108,58.08709,58.79609,
     +		59.33789,59.34961,59.36561,59.84372,59.87542,
     +		60.23497,61.10992,61.36131,61.37108,61.38501,
     +		61.40412,61.93150,61.93509,61.94088,61.94888,
     +		62.18008,62.18691,62.46812,62.48219,62.49133,
     +		63.30199,63.30286,63.32506,63.35387,63.75540,
     +		63.77412,64.30924,64.30999,66.87376,67.85857,
     +		67.86167,68.16618,68.17400,68.44416,68.44416,
     +		68.50069,68.50069,68.74507,70.50282,70.51955,
     +		70.55017,70.76975,70.76975,71.12993,71.15609,
     +		71.21387,71.21387,71.31690,71.33785,71.39315,
     +		71.39737,71.48887,71.50672,71.53300,72.12492,
     +		72.12764,72.47591,72.50269,72.88482,72.88482,
     +		73.16019,73.37047,73.37047,73.37047,73.37047,
     +		73.52322,73.52322,73.52322,73.60108,73.61112,
     +		73.64819,73.65725,73.68911,73.71453,73.93237,
     +		73.95444,74.05078,74.06293,74.06293,74.10930,
     +		74.12628,74.40265,74.40438,74.76035,74.76035,
     +		74.76035,74.76035,75.18896,75.18896,75.18896,
     +		76.30446,76.30806,76.44791,77.47625,77.47625,
     +		77.92433,77.92433,78.12258,78.12258,78.19797,
     +		78.21979,78.41159,78.43242,78.59385,78.59385,
     +		78.59385,78.59385,78.63718,78.63718,78.63718,
     +		78.85769,78.88398,78.91572,78.91572,78.96023,
     +		78.97250,78.98019,80.20107,80.20107,80.72665,
     +		80.72900,81.00314,81.01343,81.37509,81.37509,
     +		81.37509,81.37509,81.42716,81.42716,81.83012,
     +		82.88895,82.88895,83.03365,83.03365,83.03365,83.03365/
	DATA NOV/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +		3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,
     +		4,4,4,4,4,4,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,6,
     +		6,6,6,6,4,4,4,6,4,4,4,4,4,7,7,7,7,7,8,8,8,8,5,5,
     +		5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6/
	DATA GOV/1.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +		1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +		3.0D0,1.0D0,3.0D0,1.0D0,3.0D0,
     +		5.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +		3.0D0,5.0D0,7.0D0,3.0D0,1.0D0,
     +		3.0D0,5.0D0,5.0D0,5.0D0,3.0D0,
     +		5.0D0,7.0D0,1.0D0,5.0D0,3.0D0,
     +		1.0D0,7.0D0,3.0D0,3.0D0,1.0D0,
     +		1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +		5.0D0,7.0D0,5.0D0,7.0D0,3.0D0,
     +		3.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		3.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +		3.0D0,1.0D0,3.0D0,5.0D0,5.0D0,
     +		5.0D0,3.0D0,7.0D0,3.0D0,5.0D0,
     +		7.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		5.0D0,3.0D0,1.0D0,7.0D0,3.0D0,
     +		3.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +		3.0D0,3.0D0,5.0D0,7.0D0,3.0D0,
     +		3.0D0,5.0D0,7.0D0,1.0D0,3.0D0,
     +		5.0D0,5.0D0,5.0D0,3.0D0,5.0D0,
     +		7.0D0,7.0D0,3.0D0,3.0D0,5.0D0,
     +		7.0D0,1.0D0,3.0D0,5.0D0,5.0D0/
	DATA ENOV/0.0D0,10.18183,10.19878,10.23674,19.68863,
     +		26.48845,26.50776,26.54108,28.73015,35.69651,
     +		67.83862,69.59028,72.01395,72.28146,72.28596,
     +		72.29554,74.50599,74.50733,74.50979,75.95557,
     +		80.97483,80.99497,81.03748,82.38657,83.40436,
     +		83.97941,84.00407,84.04314,84.82139,85.49855,
     +		85.51269,85.53633,86.12596,86.43890,87.33036,
     +		87.33829,87.35107,87.73579,87.80076,87.81837,
     +		87.82866,88.39750,89.17985,89.60004,90.71603,
     +		91.26665,91.26665,91.26888,91.48672,92.04689,
     +		92.04763,92.04937,92.59937,92.97132,98.72523,
     +		99.49233,106.7049,106.7049,106.7049,100.2237,
     +		102.1987,102.8568,103.0377,103.0583,103.0944,
     +		103.1870,103.5465,103.5465,103.5676,103.8792,
     +		103.8829,104.1001,104.2509,104.2990,104.2990,
     +		104.2990,104.3064,104.3181,104.3333,104.4087,
     +		104.5556,104.5689,104.5754,105.0316,105.0733,
     +		106.7358,106.9963,106.9963,106.9963,106.9274,
     +		108.4187,108.5325,108.5325,108.5325,111.4108,
     +		111.5461,111.5461,111.5461,111.7535,111.7535,
     +		111.7535,111.8898,111.9082,112.1444,112.1444,
     +		112.1444,112.3809,115.9379,116.0435,116.0435,
     +		116.0435,116.1501,116.1501,116.1501,116.2166/
	DATA NOVI/2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,
     +		6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8/
	DATA GOVI/2.0D0,2.0D0,4.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +		4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +		2.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +		2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +		8.0D0,8.0D0,10.0D0,12.0D0,2.0D0,
     +		2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +		8.0D0,8.0D0,10.0D0,12.0D0,14.0D0,
     +		2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +		8.0D0,10.0D0,12.0D0,14.0D0,16.0D0,4.0D0,6.0D0/
	DATA ENOVI/0.0D0,11.94909,12.01505,79.35559,82.58831,
     +		82.60773,83.64374,83.65008,105.7219,107.0408,
     +		107.0487,107.4805,107.4831,107.5050,107.5062,
     +		117.6237,118.2920,118.2920,118.5122,118.5122,
     +		124.3735,124.3735,124.5034,124.5034,124.5142,
     +		124.5142,124.5156,124.5156,124.5156,127.8017,
     +		128.0311,128.0311,128.1171,128.1171,128.1243,
     +		128.1243,128.1252,128.1252,128.1252,128.1252,
     +		130.2520,130.3984,130.3984,130.4674,130.4674,
     +		130.4680,130.4680,130.4680,130.4680,130.4680,
     +		130.4693,130.4693/
	DATA NOVII/1,2,2,2,2,2,3,3,3,3,3,3,3,4,5,6/
	DATA GOVII/1.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +		3.0D0,1.0D0,3.0D0,5.0D0,7.0D0,
     +		5.0D0,3.0D0,3.0D0,3.0D0,3.0D0,3.0D0/
	DATA ENOVII/0.0D0,561.0761,568.6182,568.6255,568.6938,
     +		573.9532,664.1129,664.1129,664.1129,665.1804,
     +		665.1804,665.1804,665.6218,697.8022,712.7239,720.8449/
	DATA SCI/4.179704,4.179868,4.180140,4.284864,4.411317, 
     +		4.556712,4.417538,4.418036,4.419087,4.460873,
     +		5.012059,5.012145,5.012123,4.656621,4.682271,
     +		4.682931,4.683973,4.715525,4.735114,4.735517, 
     +		4.736181,4.776610,4.823287,5.245868,4.960446, 
     +		4.636463,4.637096,4.638772,4.981131,4.981795, 
     +		4.983181,4.985224,4.985476,4.985839,4.648973, 
     +		4.987257,5.002644,5.026932,5.027268,5.027448, 
     +		4.751991,4.753114,4.754886,4.775015,4.807735, 
     +		4.820394,4.821310,4.822433,4.849118,4.880081, 
     +		4.964531,4.983084,4.983084,4.983084,4.988046, 
     +		4.988550,4.989272,4.739796,4.996660,5.002712, 
     +		5.007890,5.009317,5.009317,4.830776,4.830763, 
     +		4.830763,4.843919,4.885500,4.908801,4.963564, 
     +		4.983777,4.984663,4.984663,4.989659,4.989659, 
     +		4.992449,4.793379,4.998577,5.003253,5.004741, 
     +		5.006231,5.006231,4.972326,4.984214,4.985345, 
     +		4.985345,4.991673,4.991673,4.995098,4.831870, 
     +		5.000666,5.004451,5.005935,5.005935,5.004665, 
     +		4.984616,4.985763,4.985763,4.999064,4.999064, 
     +		4.999064,5.002865,5.004758,5.006233,5.006233, 
     +		5.006233,4.984145,4.984145,4.984432,5.002990, 
     +		5.002990,5.002990,5.005628,5.006508,5.006508, 
     +		5.006508,4.983365,4.983365,4.983365,5.006886, 
     +		5.006886,5.006886,5.008002,5.012074,5.012074, 
     +		5.012074,5.014099,5.014099,5.014099,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000/ 
	DATA SCII/3.322208,3.322644,3.633090,3.633257,3.633480,
     +		3.893420,3.893440,4.089184,4.229172,4.229597, 
     +		3.436720,3.692591,3.692789,4.589103,3.953149, 
     +		3.953178,4.702748,4.702820,3.603431,3.770060, 
     +		3.770254,4.440431,4.441056,4.442240,3.961645, 
     +		3.961659,4.991754,4.992091,3.992479,3.992479, 
     +		3.697579,3.795690,3.796068,4.770876,4.780956, 
     +		3.968260,3.968260,3.994325,3.994325,3.754887, 
     +		4.893884,4.894430,4.895358,4.896706,4.906213, 
     +		4.906944,3.971236,3.971236,3.996578,3.996578, 
     +		5.011171,5.086056,5.086787,5.087796,5.188515, 
     +		5.190205,5.591661,5.735427,5.737656,5.740739, 
     +		5.745144,5.937513,5.941305,5.947732,5.956571, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000,6.000000,6.000000,6.000000, 
     +		6.000000,6.000000/
	DATA SCIII/2.247678,2.511178,2.511298,2.511595,2.783334,
     +		2.988424,2.988601,2.988887,3.040348,3.275480,
     +		2.516355,2.624130,2.770250,2.779440,2.779510,
     +		2.779673,2.913505,2.912204,2.913576,3.001503,
     +		3.471912,3.472452,3.473566,2.656192,3.501991,
     +		2.707107,2.843331,2.843331,2.843447,3.667003,
     +		2.928022,2.928408,2.928967,2.942070,2.942230,
     +		2.942442,2.952919,2.960061,3.725865,3.726323,
     +		3.727023,2.996535,3.802981,3.848412,3.848814,
     +		3.849514,3.907525,3.915897,3.920174,3.920667,
     +		3.921389,4.006182,3.997117,4.006877,2.756046,
     +		4.057183,4.057739,4.058045,4.085242,2.876864,
     +		2.910818,2.910818,2.910818,4.166811,2.957773,
     +		2.957773,2.957773,2.998549,2.998549,2.998583,
     +		2.998549,3.003525,3.005696,3.020441,3.020510,
     +		3.020601,3.088544,2.797248,2.916939,2.968366,
     +		2.968366,2.968366,3.000603,3.000603,3.000640,
     +		3.003372,3.005378,3.009464,3.009464,3.009464,
     +		3.027199,2.830505,2.926040,2.976524,2.976524,
     +		2.976524,3.009360,2.932986,2.982088,2.982088,
     +		2.982088,2.986294,2.986294,2.986294,4.828427,
     +		4.828427,4.828427,5.151012,5.224114,5.224114,
     +		5.227788,5.497263,5.497263,5.502663,5.756044,
     +		5.817122,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000,6.000000,6.000000,6.000000,6.000000,
     +		6.000000/
	DATA SCIV/1.644934,1.923884,1.924364,1.778341,1.948965,
     +		1.949284,1.998202,1.998312,1.838941,1.962835,
     +		1.963075,1.999589,1.999669,2.001418,2.001418,
     +		1.874532,1.972046,1.972243,2.001394,2.001394,
     +		2.002842,2.002845,2.003133,2.003133,1.897882,
     +		1.978616,1.978616,2.003384,2.003384,2.004821,
     +		2.004821,2.005002,2.005002,2.005023,2.005023,
     +		1.913888,1.984033,1.984033,2.005113,2.005113,
     +		2.007061,2.007061,2.007217,2.007217,2.007234,
     +		2.007234,1.989961,1.989961,2.009654,2.009654,
     +		2.009800,2.009800,2.009800,2.009800,2.009800/
	DATASCV/0.6309066,0.7688928,0.9242349,0.9241881,0.9246764,
     +		1.026124,1.003322,1.003322,1.003322,1.020118,
     +		1.021842,1.027042,1.033988,1.051912,1.003001/
	DATA SNI/4.931870,5.108953,5.109031,5.204087,5.204087,
     +		5.329969,5.330800,5.331948,5.401419,5.403777,
     +		5.968683,5.969459,5.969806,5.605717,5.641185,
     +		5.641868,5.642995,5.644538,5.662621,5.663186,
     +		5.664362,5.702335,5.703489,5.705897,5.734946,
     +		5.736104,5.797976,5.797737,5.588642,5.591228,
     +		5.594790,5.615995,5.620492,5.980873,5.982464,
     +		5.982872,5.983638,5.985012,5.986994,5.988774,
     +		5.991692,5.990931,5.991989,5.990646,5.995945,
     +		5.996395,5.996926,5.997295,6.001428,6.002394,
     +		5.745163,5.761658,5.762813,5.764937,5.768167,
     +		5.774817,5.775745,5.778173,5.802792,5.696104,
     +		5.699982,5.706142,5.715098,5.722151,5.983984,
     +		5.985277,5.987724,5.991767,5.985780,5.989671,
     +		5.990610,5.994303,5.985830,5.987479,5.992273,
     +		5.996771,5.993288,5.995173,5.998955,6.002261,
     +		6.003886,6.255743,6.257061,6.360914,6.362586,
     +		5.757127,5.763044,5.772635,5.779661,5.791555,
     +		5.984846,5.986194,5.990250,5.997386,5.990476,
     +		5.992170,5.992623,6.000597,5.993189,5.993189,
     +		5.993189,6.000802,5.996590,5.998751,6.003087,
     +		6.005032,6.007153,5.793718,5.804320,5.818227,
     +		5.815205,5.821445,5.989322,5.989322,5.989322,
     +		5.989322,5.992901,5.992901,5.992901,6.003714,
     +		5.994695,5.997311,5.995185,6.008173,6.001116,
     +		6.010741,6.005529,6.005529,6.008008,5.801159,
     +		5.821038,5.833978,5.835981,5.835981,5.989852,
     +		5.989852,5.989852,5.989852,5.993398,5.993398,
     +		5.996287,5.996287,6.005342,6.014956,6.015613,
     +		6.015613,6.015613,5.971639,5.971639,5.850568,
     +		5.850568,5.850568,5.990061,5.990061,5.990061,
     +		5.990061,5.991796,5.991796,5.993244,5.993244,
     +		6.011374,6.018783,6.017594,6.017594,6.017594,
     +		5.858175,5.858175,5.863377,5.863377,5.863377,
     +		5.988659,5.988659,5.988659,5.988659,5.989390,
     +		5.989390,5.994151,5.994151,6.020563,6.027375,
     +		6.026996,6.026996,6.026996,5.866343,5.866343,
     +		5.874646,5.874646,5.874646,5.990859,5.990859,
     +		5.992667,5.992667,5.994933,5.994933,5.994933,
     +		5.994933,6.030020,6.030020,6.038991,6.038991,
     +		6.038991,5.873279,5.873279,5.877365,5.877365,
     +		5.877365,5.992039,5.992039,5.996431,5.996431,
     +		6.000838,6.000838,6.000838,6.000838,6.039686,
     +		6.039686,6.042562,6.042562,6.042562,6.018730,
     +		5.886326,5.994597,5.994597,6.047575,6.047575,
     +		6.047575,6.078406,6.078406/
	DATA SNII/4.048939,4.049242,4.049750,4.145151,4.258360,
     +		4.356432,4.688145,4.688257,4.688270,4.826217,
     +		4.826217,4.826272,5.142618,4.284333,4.284811,
     +		4.286872,4.288557,5.253337,4.532960,4.564950,
     +		4.565974,4.567596,5.379367,4.605227,4.634192,
     +		4.634804,4.635817,4.698180,4.771750,4.928997,
     +		4.930175,4.931792,4.940517,4.947424,4.947905,
     +		4.948512,4.976005,4.977055,4.977624,4.985716,
     +		5.001774,4.517682,4.519204,4.522714,5.167543,
     +		4.688999,4.706267,4.707887,4.710950,4.720615,
     +		4.721386,4.723658,4.732426,4.791679,5.359476,
     +		5.360877,5.362645,4.824183,4.939473,4.941747,
     +		4.944790,4.948400,4.959554,4.960499,4.962098,
     +		4.976268,4.977930,4.978850,4.988034,4.988136,
     +		4.988983,4.989127,4.990716,4.997378,4.997656,
     +		5.001124,4.999108,5.001566,5.001872,5.002052,
     +		5.004650,5.004791,5.510713,5.511550,5.512880,
     +		4.633615,4.635821,4.641956,4.663452,4.970949,
     +		4.970949,4.970949,4.990886,4.990886,4.991119,
     +		4.994713,4.999842,5.000126,5.004091,5.004655,
     +		5.899771,5.900361,5.901458,5.903069,5.905087,
     +		5.975470,5.976437,5.978201,6.162114,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000/
	DATA SNIII/3.264886,3.265738,3.559230,3.555733,3.556163,
     +		3.795859,3.795903,3.971288,4.062202,4.062887,
     +		4.328298,4.441761,4.441879,3.362787,4.644640,
     +		4.644671,3.648880,3.649321,3.924340,3.924418,
     +		4.208216,4.209135,4.210838,4.353292,4.355043,
     +		3.749472,4.546074,4.546964,4.554955,4.555551,
     +		4.556746,4.558211,3.785649,3.786211,4.632736,
     +		4.686664,4.687436,4.688480,3.926233,3.926969,
     +		3.987034,3.987034,4.752837,4.754452,4.866728,
     +		4.928826,4.929522,4.930549,4.931964,4.980142,
     +		4.980413,4.980860,4.981436,3.664742,5.015918,
     +		5.016470,5.050785,5.051935,5.052674,5.126587,
     +		5.129026,3.959079,3.959144,5.192320,5.193925,
     +		3.989434,3.989434,4.005147,4.005147,3.968565,
     +		3.968565,3.992409,3.992409,4.006539,4.006539,
     +		5.571515,5.574721,5.580696,6.132490,6.134099,
     +		5.935632,5.939608,6.083616,6.087341,6.092149,
     +		6.099411,6.364476,6.365573,6.178216,6.185985,
     +		6.229222,6.316157,6.320952,6.326658,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000,7.000000,7.000000,
     +		7.000000,7.000000,7.000000/
	DATA SNIV/2.226836,2.490623,2.490878,2.491461,2.755431,
     +		2.952339,2.952669,2.953231,3.013266,3.231892,
     +		2.493613,2.749589,2.762857,2.763010,2.763353,
     +		2.901417,2.901452,2.901533,2.994477,3.382731,
     +		3.383611,3.385459,3.472422,3.564919,3.607152,
     +		3.607152,3.608737,3.645438,3.728391,3.728391,
     +		2.639493,3.779903,3.797702,3.799535,3.799535,
     +		3.799535,2.797721,2.797721,2.797721,3.872625,
     +		3.873032,3.873596,3.883205,2.857107,2.934635,
     +		2.934635,2.934635,3.951730,3.952443,3.952443,
     +		2.993446,3.029914,3.030044,3.030246,4.061023,
     +		3.127314,2.880353,2.950477,2.950477,2.950477,
     +		2.998267,2.998267,2.998267,2.959708,2.959708,
     +		2.959708,4.785085,4.873190,4.873190,4.874527,
     +		7.000000,7.000000,7.000000/
	DATA SNV/1.634565,1.915400,1.916327,1.771111,1.943796,
     +		1.944398,1.997854,1.998051,1.833396,1.959357,
     +		1.959357,1.999384,1.999384,1.870467,1.969524,
     +		1.969524,2.001982,2.001982,1.895972,1.977561,
     +		1.977561,2.004889,2.004593,2.005843,2.005843,
     +		2.006106,2.006106,2.006106,1.914798,1.983841,
     +		1.983841,2.007186,2.007186,2.008574,2.008574,
     +		2.008797,2.008797,2.008797,2.008797,1.929893,
     +		1.990742,1.990742,2.010469,2.010469,2.011696,
     +		2.011696,2.011930,2.011930,2.011930,2.011930,
     +		2.011930/
	DATA SNVI/0.6290283,0.7657156,0.9208641,0.9208946,0.9217644,
     +		0.9074407,1.024314,1.024867/
	DATA SOI/5.998809,6.000254,6.000874,6.149045,6.029965,
     +		6.280362,6.354168,6.620857,6.620917,6.621027,
     +		6.681873,6.681856,6.681878,6.554275,6.592577,
     +		6.991942,6.991948,6.991948,6.991953,6.991953,
     +		6.991953,6.994710,6.994710,6.994710,6.749977,
     +		6.750016,6.750087,6.784759,6.784759,6.784759,
     +		7.156593,7.157186,7.157570,6.676266,6.701954,
     +		7.234456,6.993919,6.993919,6.993919,6.993919,
     +		6.993919,6.997045,6.997045,6.997045,6.836974,
     +		6.836974,6.836974,6.746833,6.766088,6.996467,
     +		6.996467,6.996467,6.996467,6.996467,6.999115,
     +		6.999115,6.999115,6.869720,6.869720,6.869720,
     +		6.793480,6.808908,6.999084,6.999084,6.999084,
     +		6.999084,6.999084,7.001479,7.001479,7.001479,
     +		6.826997,6.839929,7.001803,7.001803,7.001803,
     +		7.001803,7.001803,7.003955,7.003955,7.003955,
     +		6.852827,6.864539,7.004705,8.000000,7.004705,
     +		7.004705,7.004705,7.007905,7.007905,7.007905,
     +		6.877349,6.883498,7.007761,7.007761,7.007761,
     +		7.007761,7.007761,7.010591,7.010591,7.010591,
     +		6.890946,6.900387,7.011453,7.011453,7.011453,
     +		7.011453,7.011453,7.012788,7.012788,7.012788,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000/
	DATA SOII/4.784610,4.940430,4.940555,5.022952,5.022962,
     +		5.557054,5.558274,5.558889,5.930025,5.929834,
     +		5.160761,5.162283,5.164577,5.214052,5.216705,
     +		6.210938,5.445367,5.490556,5.491464,5.492962,
     +		5.495002,5.494491,5.494507,5.517109,5.517870,
     +		5.519391,5.570158,5.573380,5.580988,6.392210,
     +		6.394130,5.615287,5.616316,5.844496,5.880436,
     +		5.880893,5.904250,5.904674,5.930839,5.931911,
     +		5.933457,5.935489,5.954107,5.955576,5.956912,
     +		5.954453,5.955388,5.959104,5.959794,5.959825,
     +		5.959952,5.960805,5.964087,5.973599,5.975908,
     +		5.993385,5.994448,5.442262,5.445265,5.449878,
     +		5.491284,5.496742,6.232405,5.654758,5.656267,
     +		5.658962,5.662895,5.725537,5.731195,5.738926,
     +		5.741796,6.348958,6.348958,6.373241,6.373322,
     +		6.408604,6.408635,6.419951,6.420006,6.439373,
     +		6.439435,5.943564,5.943564,5.943564,5.947441,
     +		5.961402,5.961402,5.961402,5.961402,6.471051,
     +		5.965786,5.967088,5.967867,5.968222,5.972136,
     +		5.972852,5.978758,6.488330,5.986873,5.991322,
     +		5.986938,5.991353,5.994984,5.994948,5.989273,
     +		5.989404,5.991879,5.997332,5.991945,5.997554,
     +		5.992443,5.994410,6.001196,6.001346,6.002642,
     +		6.003121,6.002740,6.003380,5.576062,5.580966,
     +		5.588796,5.609913,5.619139,6.121709,6.121740,
     +		5.734796,5.734796,5.738977,5.745944,5.750079,
     +		5.750079,5.754774,5.772264,5.781077,5.960446,
     +		5.960446,5.960446,5.960446,5.972282,5.974347,
     +		5.974347,5.980535,5.980535,5.983981,5.987715,
     +		5.990945,5.996456,5.985451,5.986923,5.990890,
     +		6.000214,5.991003,6.000511,5.996286,5.996286,
     +		6.000386,6.002267,6.003636,6.004436,6.003808,
     +		6.004436,6.864734,6.865004,6.865458,6.871479,
     +		6.874276,6.883227,6.883227,6.929313,6.929313,
     +		6.945119,6.945119,6.947771,6.947771,7.214550,
     +		7.214550,6.955942,6.960794,6.960794,6.963802,
     +		6.963802,6.974759,6.974759,6.897856,6.897856,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000/
	DATA SOIII/3.980091,3.980605,3.981483,4.073126,4.181011,
     +		4.263698,4.567496,2.851461,2.851508,4.688399,
     +		4.688406,4.688483,4.944256,5.004756,5.087306,
     +		4.201648,4.202927,4.205704,4.265077,5.589531,
     +		5.591178,5.591910,4.466920,4.500863,4.502461,
     +		4.505044,4.544430,5.702087,4.576288,4.577272,
     +		4.578837,4.653335,4.743010,4.880211,4.882788,
     +		4.885133,4.883788,4.916797,4.917453,4.918434,
     +		4.946753,4.948285,4.949112,4.978528,4.991553,
     +		5.071567,5.073311,5.091046,6.092469,5.236801,
     +		5.238239,5.240906,4.450987,4.453166,4.458786,
     +		4.490992,5.440957,5.477274,5.477840,5.478956,
     +		5.480584,5.482658,4.640882,4.657492,4.659830,
     +		4.664354,4.689624,5.526725,5.527675,5.529354,
     +		4.742366,4.744359,4.746791,4.755242,4.803841,
     +		5.629193,5.463434,5.633066,5.655169,4.904211,
     +		4.904211,4.904211,4.911571,5.696494,5.696715,
     +		5.697021,4.948280,4.949739,4.951245,4.983719,
     +		4.983719,4.983719,4.985557,4.992922,4.595488,
     +		4.595488,4.595488,4.614187,5.995187,5.995924,
     +		5.997336,6.003933,6.004729,6.005913,6.007472,
     +		6.009368,6.079757,6.079680,6.079593,6.079808,
     +		6.081549,6.087021,6.088523,6.089349,6.127790,
     +		6.130199,6.131379,6.150372,6.152512,6.155484,
     +		4.922875,4.922875,4.922875,4.932409,4.922875,
     +		4.922875,4.922875,4.991952,4.997716,6.251313,
     +		6.251994,6.253141,4.947119,4.974443,4.974443,
     +		4.974443,5.003925,6.782259,6.828280,6.541492,
     +		6.547457,6.555665,6.965416,7.060263,7.160808,
     +		7.160808,7.166230,7.174317,7.185196,7.256399,
     +		7.261456,7.271212,7.771374,7.771374,7.771374,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000/
	DATA SOIV/3.228562,3.230039,3.508826,3.509360,3.510109,
     +		3.741247,3.741307,3.904661,3.977057,3.978160,
     +		4.214302,4.331145,4.331291,4.503491,4.503533,
     +		3.322590,3.617383,3.618198,4.097019,4.098440,
     +		4.101037,4.249484,4.252384,4.410741,4.412061,
     +		4.420406,4.421309,4.422863,4.425270,4.491520,
     +		4.543003,4.544125,4.545658,4.591770,4.594849,
     +		3.506632,4.717019,4.742457,4.743450,4.744866,
     +		4.746809,4.800908,4.801279,4.801878,4.802707,
     +		4.826727,4.827440,4.856910,4.858391,4.859354,
     +		3.927959,3.928085,4.948471,4.951597,4.995503,
     +		4.997566,5.057139,5.057223,3.602068,5.487784,
     +		5.488193,5.528638,5.529685,3.943577,3.943577,
     +		3.956409,3.956409,5.607411,5.152442,5.155901,
     +		5.162243,5.906104,5.906104,5.285100,5.290775,
     +		3.955026,3.955026,5.325924,5.330538,6.007065,
     +		6.007765,6.023023,6.026014,6.030425,6.132527,
     +		6.133009,5.594399,5.600957,3.969004,3.969004,
     +		5.768017,5.824149,5.824149,5.824149,5.824149,
     +		5.865850,5.865850,5.865850,5.887424,5.890223,
     +		5.900586,5.903125,5.912084,5.919259,5.981793,
     +		5.988238,6.512458,4.040435,4.040435,6.034046,
     +		6.039135,6.592915,6.593322,6.679720,6.679720,
     +		6.679720,6.679720,6.791924,6.791924,6.791924,
     +		7.150804,7.152208,7.208681,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000,
     +		8.000000/
	DATA SOV/2.212299,2.477108,2.477559,2.478570,2.736373,
     +		2.929941,2.930501,2.931468,2.995396,3.204502,
     +		2.480139,2.586177,2.736415,2.753261,2.753545,
     +		2.754149,2.895501,2.895588,2.895747,2.990361,
     +		3.333698,3.335127,3.338143,3.434917,3.509306,
     +		3.551885,3.553720,3.556629,3.614975,3.666382,
     +		3.667461,3.669268,3.714561,3.738797,3.808601,
     +		3.809227,3.810236,3.840735,3.845908,3.847311,
     +		3.848131,3.893723,3.957266,2.655746,2.780047,
     +		2.842480,2.842480,2.842734,2.867645,2.932266,
     +		2.932352,2.932553,2.996815,3.040747,2.722722,
     +		2.858081,4.369750,4.369750,4.369750,2.990545,
     +		4.293693,4.399677,4.429359,4.432752,4.438707,
     +		4.454041,4.514207,4.514207,4.517767,4.570812,
     +		4.571450,2.913395,2.952782,2.965416,2.965416,
     +		2.965416,4.644915,4.646958,4.649636,2.994349,
     +		4.688903,4.691261,4.692408,4.774586,4.782194,
     +		2.928606,3.022013,3.022013,3.022013,2.997125,
     +		2.933282,2.986425,2.986425,2.986425,5.872365,
     +		5.931636,5.931636,5.931636,6.025977,6.025977,
     +		6.025977,6.090481,6.099398,6.217293,6.217293,
     +		6.217293,6.343698,8.000000,8.000000,8.000000,
     +		8.000000,8.000000,8.000000,8.000000,8.000000/
	DATA SOVI/1.626749,1.908750,1.910343,1.765577,1.939605,
     +		1.940666,1.997515,1.997865,1.829541,1.956604,
     +		1.957376,1.999562,1.999822,2.001964,2.002083,
     +		1.867336,1.968342,1.968342,2.001995,2.001995,
     +		1.976059,1.976059,2.004681,2.004681,2.007063,
     +		2.007063,2.007365,2.007365,2.007365,1.914098,
     +		1.982390,1.982390,2.008208,2.008208,2.010369,
     +		2.010369,2.010631,2.010631,2.010631,2.010631,
     +		1.930105,1.987143,1.987143,2.014184,2.014184,
     +		2.014431,2.014431,2.014431,2.014431,2.014431,
     +		2.014965,2.014965/
	DATA SOVII/0.6273875,0.7631111,0.9180546,0.9182081,0.9196253,
     +		1.029738,0.9543552,0.9543552,0.9543552,1.004676,
     +		1.004676,1.004676,1.025589,1.027917,1.034432,
     +		1.045346/
*
*	Find index for atom and ion, 10*IAT+IZI
*
c       IF(IAT.EQ.26.AND.IZI.GE.6.AND.IZI.LE.9) GO TO 260
	IF(IAT.GT.2.AND.IAT.LT.6)GO TO 9999
	IF(IAT.LT.1.OR.IAT.GT.8)GO TO 9999
	IND=10*IAT+IZI
	IF(IND.EQ.11) GO TO 11
	IF(IND.EQ.21) GO TO 21
	IF(IND.EQ.22) GO TO 22
	IF(IND.EQ.61) GO TO 61
	IF(IND.EQ.61) GO TO 62
	IF(IND.EQ.63) GO TO 63
	IF(IND.EQ.64) GO TO 64
	IF(IND.EQ.65) GO TO 65
	IF(IND.EQ.66) GO TO 66
	IF(IND.EQ.71) GO TO 71
	IF(IND.EQ.72) GO TO 72
	IF(IND.EQ.73) GO TO 73
	IF(IND.EQ.74) GO TO 74
	IF(IND.EQ.75) GO TO 75
	IF(IND.EQ.76) GO TO 76
	IF(IND.EQ.77) GO TO 77
	IF(IND.EQ.81) GO TO 81
	IF(IND.EQ.82) GO TO 82
	IF(IND.EQ.83) GO TO 83
	IF(IND.EQ.84) GO TO 84
	IF(IND.EQ.85) GO TO 85
	IF(IND.EQ.86) GO TO 86
	IF(IND.EQ.87) GO TO 87
	IF(IND.EQ.88) GO TO 88
* 
*	CALCULATING PARTITION FUNCTIONS FOR HYDROGEN
* 
 11	CALL PARTDV(T,ANE,ZH,MH,NHYD,GHYD,ENHYD,SHYD,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR HEI
* 
 21	CALL PARTDV(T,ANE,ZHE,MHEI,NHEL,GHEL,ENHEL,SHEL,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR HEII
* 
 22	CALL PARTDV(T,ANE,ZHE,MHEII,NHYD,GHYD,ENHYD,SHYD,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR CI
* 
 61	CALL PARTDV(T,ANE,ZC,MCI,NCI,GCI,ENCI,SCI,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR CII
* 
 62	CALL PARTDV(T,ANE,ZC,MCII,NCII,GCII,ENCII,SCII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR CIII
* 
 63	CALL PARTDV(T,ANE,ZC,MCIII,NCIII,GCIII,ENCIII,SCIII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR CIV
* 
 64	CALL PARTDV(T,ANE,ZC,MCIV,NCIV,GCIV,ENCIV,SCIV,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR CV
* 
 65	CALL PARTDV(T,ANE,ZC,MCV,NCV,GCV,ENCV,SCV,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR CVI
* 
 66	CALL PARTDV(T,ANE,ZC,MH,NHYD,GHYD,ENHYD,SHYD,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NI
* 
 71	CALL PARTDV(T,ANE,ZN,MNI,NNI,GNI,ENNI,SNI,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NII
* 
 72	CALL PARTDV(T,ANE,ZN,MNII,NNII,GNII,ENNII,SNII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NIII
* 
 73	CALL PARTDV(T,ANE,ZN,MNIII,NNIII,GNIII,ENNIII,SNIII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NIV
* 
 74	CALL PARTDV(T,ANE,ZN,MNIV,NNIV,GNIV,ENNIV,SNIV,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NV
* 
 75	CALL PARTDV(T,ANE,ZN,MNV,NNV,GNV,ENNV,SNV,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NVI
* 
 76	CALL PARTDV(T,ANE,ZN,MNVI,NNVI,GNVI,ENNVI,SNVI,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR NVII
* 
 77	CALL PARTDV(T,ANE,ZN,MH,NHYD,GHYD,ENHYD,SHYD,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OI
* 
 81	CALL PARTDV(T,ANE,ZO,MOI,NOI,GOI,ENOI,SOI,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OII
* 
 82	CALL PARTDV(T,ANE,ZO,MOII,NOII,GOII,ENOII,SOII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OIII
* 
 83	CALL PARTDV(T,ANE,ZO,MOIII,NOIII,GOIII,ENOIII,SOIII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OIV
* 
 84	CALL PARTDV(T,ANE,ZO,MOIV,NOIV,GOIV,ENOIV,SOIV,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OV
* 
 85	CALL PARTDV(T,ANE,ZO,MOV,NOV,GOV,ENOV,SOV,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OVI
* 
 86	CALL PARTDV(T,ANE,ZO,MOVI,NOVI,GOVI,ENOVI,SOVI,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OVII
* 
 87	CALL PARTDV(T,ANE,ZO,MOVII,NOVII,GOVII,ENOVII,SOVII,U)
	GO TO 8888
* 
*	CALCULATING PARTITION FUNCTIONS FOR OVIII
* 
 88	CALL PARTDV(T,ANE,ZO,MH,NHYD,GHYD,ENHYD,SHYD,U)
	GO TO 8888
C
C 
C	CALCULATING PARTITION FUNCTIONS FOR FE VI - FE IX
C 
C260    CALL PFFE(IZI,T,ANE,U)
 8888	CONTINUE
	RETURN
 9999	U=0
	WRITE(*,*)!! INVALID ATOM IN USER SUPPLIED ROUTINE PARTFUN !!
	STOP
	END SUBROUTINE
C
C     **************************************************************
C
* 
	SUBROUTINE PARTDV(TEMP,DNE,Z,NLEV,NE,GEE,ENRGY,S,U)
* 
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(in   ) :: NLEV,NE(*)
      real*8, intent(in   ) :: TEMP,DNE,Z,GEE,ENRGY,S
      real*8, intent(inout) :: U
	DIMENSION GEE(*),ENRGY(*),S(*)
	U=0.0
	ET=TEMP/11604.8
	P=(14.69D0-0.20-0.6667*LOG10(DNE))
* 
	DO 10 I=1,NLEV
	   U1=FLOAT(NE(I))
	   ZSTAR=Z-S(I)
	   IF (ZSTAR.GT.0)THEN
	       W=P+4.*LOG10(ZSTAR)-4.*LOG10(U1)
	    ELSE
	       W=0.0
	   ENDIF
	   IF (W.GT.1.) W=1.
* 
	   IF ((ENRGY(I)/ET).LT.65.0) THEN
		U1=GEE(I)*W*EXP(-ENRGY(I)/ET)
	    ELSE
		U1=0.0
	   ENDIF
	   U=U+U1
 10	CONTINUE
 	RETURN
	END SUBROUTINE
C
C     **************************************************************
C
!       subroutine pfni(ion,t,ane,pf,dut,dun)
! c     =====================================
! c
! c     partition functions for Ni IV to Ni IX
! c
! c     this routine interpolates within a grid
! c     calculated from all levels predicted by
! c     Kurucz (1992), i.e. over 12,000 levels per ion.
! c     the partition functions depend only on T !
! c     (i.e. no level dissolution with increasing density)
! c     TL  27-DEC-1994, 23-JAN-1995
! c
! c     Output:  PF   partition function
! c              DUT  d(PF)/dT
! c              DUN  d(PF)/d(ANE) (=0 in this version)
! c
!       implicit double precision (a-h,o-z)
!       integer,intent(in   ) :: ion
!       real*8, intent(in   ) :: t,ane
!       real*8, intent(inout) :: pf,dut,dun
!         
! c
!       dimension g0(6)
!       dimension p4a(190),p4b(170)
!       dimension p5a(190),p5b(170)
!       dimension p6a(190),p6b(170)
!       dimension p7a(190),p7b(170)
!       dimension p8a(190),p8b(170)
!       dimension p9a(190),p9b(170)
!       parameter (xen=2.302585093,xmil=0.001)
! c
!       data g0/28.,25.,6.,25.,28.,21./
! c
!       data p4a/
!      .    1.447,1.464,1.482,1.501,1.518,1.535,1.551,1.567,1.582,1.596,
!      .    1.610,1.623,1.636,1.648,1.659,1.671,1.681,1.692,1.702,1.711,
!      .    1.721,1.730,1.739,1.748,1.757,1.765,1.774,1.782,1.791,1.799,
!      .    1.808,1.816,1.824,1.833,1.841,1.850,1.859,1.868,1.877,1.886,
!      .    1.895,1.905,1.914,1.924,1.934,1.945,1.955,1.966,1.977,1.989,
!      .    2.000,2.012,2.025,2.037,2.050,2.063,2.077,2.091,2.105,2.119,
!      .    2.134,2.149,2.164,2.179,2.195,2.211,2.227,2.243,2.260,2.276,
!      .    2.293,2.310,2.327,2.344,2.362,2.379,2.397,2.414,2.432,2.449,
!      .    2.467,2.484,2.502,2.519,2.537,2.554,2.571,2.588,2.606,2.623,
!      .    2.640,2.657,2.674,2.690,2.707,2.723,2.740,2.756,2.772,2.788,
!      .    2.804,2.819,2.835,2.850,2.866,2.881,2.896,2.911,2.925,2.940,
!      .    2.954,2.969,2.983,2.997,3.010,3.024,3.038,3.051,3.064,3.077,
!      .    3.090,3.103,3.116,3.128,3.141,3.153,3.165,3.177,3.189,3.201,
!      .    3.213,3.224,3.235,3.247,3.258,3.269,3.280,3.291,3.301,3.312,
!      .    3.322,3.332,3.343,3.353,3.363,3.373,3.382,3.392,3.402,3.411,
!      .    3.421,3.430,3.439,3.448,3.457,3.466,3.475,3.484,3.492,3.501,
!      .    3.509,3.518,3.526,3.534,3.542,3.550,3.558,3.566,3.574,3.582,
!      .    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,
!      .    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722/
!       data p4b/
!      .    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,
!      .    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722,
!      .    3.728,3.734,3.740,3.747,3.753,3.759,3.765,3.771,3.777,3.782,
!      .    3.788,3.794,3.800,3.805,3.811,3.816,3.822,3.827,3.833,3.838,
!      .    3.843,3.849,3.854,3.859,3.864,3.869,3.874,3.879,3.884,3.889,
!      .    3.894,3.899,3.904,3.909,3.913,3.918,3.923,3.927,3.932,3.936,
!      .    3.941,3.945,3.950,3.954,3.959,3.963,3.967,3.972,3.976,3.980,
!      .    3.984,3.988,3.993,3.997,4.001,4.005,4.009,4.013,4.017,4.021,
!      .    4.024,4.028,4.032,4.036,4.040,4.043,4.047,4.051,4.055,4.058,
!      .    4.062,4.065,4.069,4.072,4.076,4.079,4.083,4.086,4.090,4.093,
!      .    4.097,4.100,4.103,4.107,4.110,4.113,4.116,4.120,4.123,4.126,
!      .    4.129,4.132,4.135,4.138,4.141,4.144,4.148,4.151,4.154,4.157,
!      .    4.159,4.162,4.165,4.168,4.171,4.174,4.177,4.180,4.182,4.185,
!      .    4.188,4.191,4.193,4.196,4.199,4.202,4.204,4.207,4.210,4.212,
!      .    4.215,4.217,4.220,4.223,4.225,4.228,4.230,4.233,4.235,4.238,
!      .    4.240,4.243,4.245,4.247,4.250,4.252,4.255,4.257,4.259,4.262,
!      .    4.264,4.266,4.268,4.271,4.273,4.275,4.278,4.280,4.282,4.284/
!       data p5a/
!      .    1.398,1.408,1.427,1.446,1.466,1.486,1.506,1.526,1.545,1.564,
!      .    1.583,1.601,1.619,1.636,1.652,1.668,1.683,1.698,1.712,1.725,
!      .    1.738,1.751,1.763,1.775,1.786,1.797,1.808,1.818,1.828,1.837,
!      .    1.846,1.855,1.864,1.873,1.881,1.889,1.897,1.904,1.912,1.919,
!      .    1.926,1.933,1.940,1.946,1.953,1.960,1.966,1.972,1.979,1.985,
!      .    1.991,1.997,2.003,2.009,2.016,2.022,2.028,2.034,2.040,2.046,
!      .    2.052,2.058,2.065,2.071,2.077,2.084,2.090,2.097,2.103,2.110,
!      .    2.117,2.124,2.131,2.138,2.145,2.152,2.160,2.167,2.175,2.183,
!      .    2.191,2.199,2.207,2.216,2.224,2.233,2.241,2.250,2.259,2.268,
!      .    2.278,2.287,2.297,2.306,2.316,2.326,2.336,2.346,2.356,2.367,
!      .    2.377,2.387,2.398,2.409,2.419,2.430,2.441,2.452,2.463,2.474,
!      .    2.485,2.497,2.508,2.519,2.530,2.542,2.553,2.564,2.576,2.587,
!      .    2.599,2.610,2.621,2.633,2.644,2.655,2.667,2.678,2.689,2.701,
!      .    2.712,2.723,2.734,2.745,2.757,2.768,2.779,2.790,2.801,2.812,
!      .    2.822,2.833,2.844,2.855,2.865,2.876,2.886,2.897,2.907,2.918,
!      .    2.928,2.938,2.948,2.958,2.968,2.978,2.988,2.998,3.008,3.018,
!      .    3.027,3.037,3.046,3.056,3.065,3.075,3.084,3.093,3.102,3.111,
!      .    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,
!      .    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279/
!       data p5b/
!      .    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,
!      .    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279,
!      .    3.287,3.295,3.303,3.310,3.318,3.325,3.333,3.340,3.347,3.355,
!      .    3.362,3.369,3.376,3.383,3.390,3.397,3.404,3.411,3.417,3.424,
!      .    3.431,3.438,3.444,3.451,3.457,3.464,3.470,3.476,3.483,3.489,
!      .    3.495,3.501,3.507,3.514,3.520,3.526,3.531,3.537,3.543,3.549,
!      .    3.555,3.561,3.566,3.572,3.578,3.583,3.589,3.594,3.600,3.605,
!      .    3.610,3.616,3.621,3.626,3.632,3.637,3.642,3.647,3.652,3.657,
!      .    3.662,3.667,3.672,3.677,3.682,3.687,3.692,3.697,3.701,3.706,
!      .    3.711,3.716,3.720,3.725,3.729,3.734,3.738,3.743,3.747,3.752,
!      .    3.756,3.761,3.765,3.769,3.774,3.778,3.782,3.786,3.790,3.795,
!      .    3.799,3.803,3.807,3.811,3.815,3.819,3.823,3.827,3.831,3.835,
!      .    3.839,3.843,3.846,3.850,3.854,3.858,3.862,3.865,3.869,3.873,
!      .    3.876,3.880,3.884,3.887,3.891,3.894,3.898,3.901,3.905,3.908,
!      .    3.912,3.915,3.918,3.922,3.925,3.929,3.932,3.935,3.939,3.942,
!      .    3.945,3.948,3.951,3.955,3.958,3.961,3.964,3.967,3.970,3.974,
!      .    3.977,3.980,3.983,3.986,3.989,3.992,3.995,3.998,4.001,4.004/
!       data p6a/
!      .    0.778,0.804,0.817,0.834,0.854,0.876,0.901,0.928,0.957,0.987,
!      .    1.017,1.048,1.079,1.109,1.139,1.169,1.197,1.225,1.253,1.279,
!      .    1.304,1.329,1.353,1.376,1.398,1.419,1.440,1.459,1.478,1.497,
!      .    1.515,1.532,1.548,1.564,1.580,1.594,1.609,1.623,1.636,1.649,
!      .    1.662,1.674,1.686,1.698,1.709,1.720,1.730,1.740,1.750,1.760,
!      .    1.769,1.779,1.788,1.796,1.805,1.813,1.821,1.829,1.837,1.845,
!      .    1.852,1.860,1.867,1.874,1.881,1.888,1.894,1.901,1.907,1.914,
!      .    1.920,1.926,1.932,1.938,1.944,1.950,1.956,1.962,1.968,1.974,
!      .    1.979,1.985,1.991,1.996,2.002,2.007,2.013,2.018,2.024,2.029,
!      .    2.035,2.041,2.046,2.052,2.057,2.063,2.068,2.074,2.080,2.086,
!      .    2.091,2.097,2.103,2.109,2.115,2.121,2.127,2.133,2.139,2.145,
!      .    2.152,2.158,2.164,2.171,2.177,2.184,2.190,2.197,2.204,2.211,
!      .    2.218,2.225,2.232,2.239,2.246,2.253,2.261,2.268,2.276,2.283,
!      .    2.291,2.298,2.306,2.314,2.322,2.330,2.338,2.346,2.354,2.362,
!      .    2.370,2.379,2.387,2.395,2.404,2.412,2.420,2.429,2.438,2.446,
!      .    2.455,2.463,2.472,2.481,2.489,2.498,2.507,2.516,2.524,2.533,
!      .    2.542,2.551,2.560,2.569,2.577,2.586,2.595,2.604,2.613,2.622,
!      .    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,
!      .    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796/
!       data p6b/
!      .    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,
!      .    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796,
!      .    2.804,2.812,2.821,2.829,2.838,2.846,2.854,2.862,2.871,2.879,
!      .    2.887,2.895,2.903,2.911,2.919,2.927,2.935,2.943,2.951,2.958,
!      .    2.966,2.974,2.982,2.989,2.997,3.005,3.012,3.020,3.027,3.035,
!      .    3.042,3.049,3.057,3.064,3.071,3.078,3.086,3.093,3.100,3.107,
!      .    3.114,3.121,3.128,3.135,3.141,3.148,3.155,3.162,3.169,3.175,
!      .    3.182,3.188,3.195,3.202,3.208,3.214,3.221,3.227,3.234,3.240,
!      .    3.246,3.252,3.259,3.265,3.271,3.277,3.283,3.289,3.295,3.301,
!      .    3.307,3.313,3.319,3.325,3.330,3.336,3.342,3.348,3.353,3.359,
!      .    3.364,3.370,3.376,3.381,3.386,3.392,3.397,3.403,3.408,3.413,
!      .    3.419,3.424,3.429,3.434,3.440,3.445,3.450,3.455,3.460,3.465,
!      .    3.470,3.475,3.480,3.485,3.490,3.495,3.499,3.504,3.509,3.514,
!      .    3.518,3.523,3.528,3.533,3.537,3.542,3.546,3.551,3.555,3.560,
!      .    3.564,3.569,3.573,3.578,3.582,3.586,3.591,3.595,3.599,3.604,
!      .    3.608,3.612,3.616,3.621,3.625,3.629,3.633,3.637,3.641,3.645,
!      .    3.649,3.653,3.657,3.661,3.665,3.669,3.673,3.677,3.681,3.685/
!       data p7a/
!      .    1.398,1.398,1.398,1.398,1.406,1.425,1.443,1.461,1.480,1.498,
!      .    1.516,1.534,1.551,1.568,1.585,1.601,1.616,1.631,1.646,1.660,
!      .    1.674,1.687,1.700,1.712,1.724,1.736,1.747,1.758,1.768,1.778,
!      .    1.788,1.797,1.806,1.815,1.824,1.832,1.840,1.848,1.855,1.863,
!      .    1.870,1.877,1.883,1.890,1.896,1.902,1.908,1.914,1.920,1.925,
!      .    1.931,1.936,1.941,1.946,1.951,1.956,1.960,1.965,1.969,1.974,
!      .    1.978,1.982,1.986,1.990,1.994,1.998,2.001,2.005,2.009,2.012,
!      .    2.016,2.019,2.022,2.026,2.029,2.032,2.035,2.038,2.041,2.044,
!      .    2.047,2.050,2.053,2.056,2.059,2.061,2.064,2.067,2.069,2.072,
!      .    2.075,2.077,2.080,2.082,2.085,2.088,2.090,2.093,2.095,2.098,
!      .    2.100,2.103,2.105,2.107,2.110,2.112,2.115,2.117,2.120,2.122,
!      .    2.125,2.127,2.130,2.132,2.135,2.137,2.140,2.142,2.145,2.148,
!      .    2.150,2.153,2.155,2.158,2.161,2.163,2.166,2.169,2.172,2.175,
!      .    2.178,2.180,2.183,2.186,2.189,2.192,2.195,2.198,2.202,2.205,
!      .    2.208,2.211,2.215,2.218,2.221,2.225,2.228,2.232,2.235,2.239,
!      .    2.243,2.246,2.250,2.254,2.258,2.261,2.265,2.269,2.273,2.277,
!      .    2.282,2.286,2.290,2.294,2.299,2.303,2.307,2.312,2.316,2.321,
!      .    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,
!      .    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420/
!       data p7b/
!      .    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,
!      .    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420,
!      .    2.426,2.431,2.437,2.442,2.448,2.453,2.459,2.464,2.470,2.476,
!      .    2.481,2.487,2.493,2.498,2.504,2.510,2.516,2.521,2.527,2.533,
!      .    2.539,2.545,2.551,2.556,2.562,2.568,2.574,2.580,2.586,2.592,
!      .    2.598,2.604,2.610,2.616,2.622,2.628,2.634,2.640,2.646,2.652,
!      .    2.658,2.664,2.670,2.676,2.682,2.687,2.693,2.699,2.705,2.711,
!      .    2.717,2.723,2.729,2.735,2.741,2.747,2.753,2.759,2.764,2.770,
!      .    2.776,2.782,2.788,2.794,2.799,2.805,2.811,2.817,2.823,2.828,
!      .    2.834,2.840,2.846,2.851,2.857,2.863,2.868,2.874,2.879,2.885,
!      .    2.891,2.896,2.902,2.907,2.913,2.918,2.924,2.929,2.935,2.940,
!      .    2.945,2.951,2.956,2.962,2.967,2.972,2.978,2.983,2.988,2.993,
!      .    2.999,3.004,3.009,3.014,3.019,3.025,3.030,3.035,3.040,3.045,
!      .    3.050,3.055,3.060,3.065,3.070,3.075,3.080,3.085,3.090,3.095,
!      .    3.099,3.104,3.109,3.114,3.119,3.123,3.128,3.133,3.138,3.142,
!      .    3.147,3.152,3.156,3.161,3.165,3.170,3.175,3.179,3.184,3.188,
!      .    3.193,3.197,3.202,3.206,3.210,3.215,3.219,3.224,3.228,3.232/
!       data p8a/
!      .    1.447,1.447,1.447,1.447,1.447,1.447,1.459,1.475,1.489,1.504,
!      .    1.518,1.531,1.544,1.556,1.568,1.580,1.591,1.602,1.612,1.622,
!      .    1.631,1.640,1.649,1.658,1.666,1.674,1.682,1.689,1.696,1.703,
!      .    1.710,1.716,1.722,1.728,1.734,1.740,1.745,1.751,1.756,1.761,
!      .    1.766,1.770,1.775,1.779,1.784,1.788,1.792,1.796,1.800,1.804,
!      .    1.807,1.811,1.814,1.818,1.821,1.824,1.827,1.831,1.834,1.836,
!      .    1.839,1.842,1.845,1.848,1.850,1.853,1.855,1.858,1.860,1.863,
!      .    1.865,1.867,1.870,1.872,1.874,1.876,1.878,1.880,1.882,1.884,
!      .    1.886,1.888,1.890,1.892,1.894,1.896,1.898,1.900,1.902,1.903,
!      .    1.905,1.907,1.909,1.911,1.912,1.914,1.916,1.917,1.919,1.921,
!      .    1.923,1.924,1.926,1.928,1.929,1.931,1.933,1.934,1.936,1.938,
!      .    1.939,1.941,1.943,1.945,1.946,1.948,1.950,1.951,1.953,1.955,
!      .    1.957,1.959,1.960,1.962,1.964,1.966,1.968,1.970,1.971,1.973,
!      .    1.975,1.977,1.979,1.981,1.983,1.985,1.987,1.989,1.991,1.993,
!      .    1.995,1.998,2.000,2.002,2.004,2.006,2.009,2.011,2.013,2.015,
!      .    2.018,2.020,2.023,2.025,2.027,2.030,2.032,2.035,2.037,2.040,
!      .    2.043,2.045,2.048,2.051,2.053,2.056,2.059,2.062,2.064,2.067,
!      .    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,
!      .    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130/
!       data p8b/
!      .    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,
!      .    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130,
!      .    2.133,2.137,2.140,2.143,2.147,2.151,2.154,2.158,2.161,2.165,
!      .    2.168,2.172,2.176,2.180,2.183,2.187,2.191,2.195,2.198,2.202,
!      .    2.206,2.210,2.214,2.218,2.222,2.226,2.230,2.233,2.237,2.241,
!      .    2.245,2.250,2.254,2.258,2.262,2.266,2.270,2.274,2.278,2.282,
!      .    2.286,2.291,2.295,2.299,2.303,2.307,2.312,2.316,2.320,2.324,
!      .    2.329,2.333,2.337,2.341,2.346,2.350,2.354,2.359,2.363,2.367,
!      .    2.371,2.376,2.380,2.384,2.389,2.393,2.397,2.402,2.406,2.410,
!      .    2.415,2.419,2.423,2.428,2.432,2.436,2.441,2.445,2.449,2.454,
!      .    2.458,2.462,2.467,2.471,2.475,2.480,2.484,2.488,2.493,2.497,
!      .    2.501,2.506,2.510,2.514,2.519,2.523,2.527,2.531,2.536,2.540,
!      .    2.544,2.548,2.553,2.557,2.561,2.565,2.570,2.574,2.578,2.582,
!      .    2.586,2.591,2.595,2.599,2.603,2.607,2.611,2.616,2.620,2.624,
!      .    2.628,2.632,2.636,2.640,2.644,2.648,2.652,2.656,2.661,2.665,
!      .    2.669,2.673,2.677,2.681,2.685,2.689,2.693,2.696,2.700,2.704,
!      .    2.708,2.712,2.716,2.720,2.724,2.728,2.732,2.736,2.739,2.743/
!       data p9a/
!      .    1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.325,
!      .    1.334,1.342,1.351,1.358,1.366,1.373,1.380,1.386,1.392,1.398,
!      .    1.404,1.409,1.415,1.420,1.425,1.429,1.434,1.438,1.442,1.446,
!      .    1.450,1.454,1.457,1.461,1.464,1.467,1.470,1.473,1.476,1.479,
!      .    1.482,1.485,1.487,1.490,1.492,1.495,1.497,1.499,1.501,1.503,
!      .    1.505,1.507,1.509,1.511,1.513,1.515,1.517,1.519,1.520,1.522,
!      .    1.524,1.525,1.527,1.528,1.530,1.531,1.533,1.534,1.535,1.537,
!      .    1.538,1.539,1.541,1.542,1.543,1.545,1.546,1.547,1.548,1.549,
!      .    1.551,1.552,1.553,1.554,1.555,1.556,1.558,1.559,1.560,1.561,
!      .    1.562,1.563,1.565,1.566,1.567,1.568,1.569,1.570,1.571,1.573,
!      .    1.574,1.575,1.576,1.577,1.579,1.580,1.581,1.582,1.584,1.585,
!      .    1.586,1.588,1.589,1.590,1.592,1.593,1.594,1.596,1.597,1.599,
!      .    1.600,1.602,1.603,1.605,1.606,1.608,1.609,1.611,1.612,1.614,
!      .    1.616,1.617,1.619,1.621,1.622,1.624,1.626,1.628,1.630,1.631,
!      .    1.633,1.635,1.637,1.639,1.641,1.643,1.645,1.647,1.649,1.651,
!      .    1.653,1.655,1.657,1.659,1.661,1.664,1.666,1.668,1.670,1.673,
!      .    1.675,1.677,1.679,1.682,1.684,1.686,1.689,1.691,1.694,1.696,
!      .    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,
!      .    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749/
!       data p9b/
!      .    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,
!      .    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749,
!      .    1.751,1.754,1.757,1.760,1.763,1.765,1.768,1.771,1.774,1.777,
!      .    1.780,1.783,1.786,1.789,1.792,1.795,1.798,1.801,1.804,1.807,
!      .    1.810,1.813,1.816,1.819,1.822,1.825,1.828,1.831,1.834,1.837,
!      .    1.840,1.843,1.847,1.850,1.853,1.856,1.859,1.862,1.865,1.869,
!      .    1.872,1.875,1.878,1.881,1.884,1.888,1.891,1.894,1.897,1.901,
!      .    1.904,1.907,1.910,1.913,1.917,1.920,1.923,1.926,1.930,1.933,
!      .    1.936,1.939,1.943,1.946,1.949,1.952,1.956,1.959,1.962,1.965,
!      .    1.969,1.972,1.975,1.978,1.982,1.985,1.988,1.992,1.995,1.998,
!      .    2.001,2.005,2.008,2.011,2.014,2.018,2.021,2.024,2.027,2.031,
!      .    2.034,2.037,2.040,2.044,2.047,2.050,2.053,2.057,2.060,2.063,
!      .    2.066,2.070,2.073,2.076,2.079,2.083,2.086,2.089,2.092,2.095,
!      .    2.099,2.102,2.105,2.108,2.111,2.115,2.118,2.121,2.124,2.127,
!      .    2.131,2.134,2.137,2.140,2.143,2.146,2.149,2.153,2.156,2.159,
!      .    2.162,2.165,2.168,2.171,2.175,2.178,2.181,2.184,2.187,2.190,
!      .    2.193,2.196,2.199,2.202,2.205,2.208,2.212,2.215,2.218,2.221/
! c
!       if(t.lt.12000.) then
!         pf=g0(ion-3)
!         dut=0.
!         dun=0.
!         return
!       endif
! c
!       it=t/1000
!       if(it.ge.350) it=349
!       t1=1000.*it
!       t2=t1+1000.
!       if(ion.eq.4) then
!         if(t.le.200000.) then
!           xu1=p4a(it-10)
!           xu2=p4a(it-9)
!         else
!           xu1=p4b(it-180)
!           xu2=p4b(it-179)
!         endif
!       else if(ion.eq.5) then
!         if(t.le.200000.) then
!           xu1=p5a(it-10)
!           xu2=p5a(it-9)
!         else
!           xu1=p5b(it-180)
!           xu2=p5b(it-179)
!         endif
!       else if(ion.eq.6) then
!         if(t.le.200000.) then
!           xu1=p6a(it-10)
!           xu2=p6a(it-9)
!         else
!           xu1=p6b(it-180)
!           xu2=p6b(it-179)
!         endif
!       else if(ion.eq.7) then
!         if(t.le.200000.) then
!           xu1=p7a(it-10)
!           xu2=p7a(it-9)
!         else
!           xu1=p7b(it-180)
!           xu2=p7b(it-179)
!         endif
!       else if(ion.eq.8) then
!         if(t.le.200000.) then
!           xu1=p8a(it-10)
!           xu2=p8a(it-9)
!         else
!           xu1=p8b(it-180)
!           xu2=p8b(it-179)
!         endif
!       else if(ion.eq.9) then
!         if(t.le.200000.) then
!           xu1=p9a(it-10)
!           xu2=p9a(it-9)
!         else
!           xu1=p9b(it-180)
!           xu2=p9b(it-179)
!         endif
!       endif
! c
!       dxt=xmil*(xu2-xu1)
!       xu=xu1+(t-t1)*dxt
!       pf=exp(xen*xu)
!       dut=xen*pf*dxt
!       dun=0.
!       return
!       end subroutine
c 
c ******************************************************************
c
      subroutine frac1
c     ================
c
C      INCLUDE 'PARAMS.FOR'
C      include 'MODELP.FOR'
      implicit real*8(a-h,o-z)
	INCLUDE '../inc/PARAMS.FOR'
	INCLUDE '../inc/MODELP.FOR'
      parameter (mtemp=100,melec=60,mion1=30)
      dimension xt(mdepth),xne(mdepth)
      dimension kt0(mdepth),kn0(mdepth)
      common/fracop/frac(mtemp,melec,mion1),fracm(mtemp,melec),
     *              itemp(mtemp),ntt
c
      do 10 id=1,nd
         xt(id)=dlog10(temp(id))
         kt0(id)=2*int(20.*xt(id))
         xne(id)=dlog10(elec(id))
         kn0(id)=int(2.*xne(id))
   10 continue
c
      DO 20 IAT=1,30
         iatnum=iat
         call fractn(iatnum)
         if(iatnum.le.0) goto 20
         do 30 id=1,nd
           if(kt0(id).lt.itemp(1)) then
             kt1=1
             write(6,611) id,temp(id)
  611        format(' (FRACOP) Extrapol. in T (low)',i4,f7.0)
             goto 41
           endif
           if(kt0(id).ge.itemp(ntt)) then
             kt1=ntt-1
             write(6,612) id,temp(id)
  612        format(' (FRACOP) Extrapol. in T (high)',i4,f12.0)
             goto 41
           endif
           do 40 it=1,ntt
             if(kt0(id).eq.itemp(it)) then
               kt1=it
               goto 41
             endif
   40      continue
   41      continue
           if(kn0(id).lt.1) then
             kn1=1
             goto 49
           endif
           if(kn0(id).ge.60) then
             kn1=59
             write(6,614) id,xne(id)
  614        format(' (FRACOP) Extrapol. in Ne (high)',i4,f9.4)
             goto 49
           endif
           kn1=kn0(id)
   49      continue
           xt1=0.025*itemp(kt1)
           dxt=0.05
           at1=(xt(id)-xt1)/dxt
           xn1=0.5*kn1
           dxn=0.5
           an1=(xne(id)-xn1)/dxn
           do 50 ion=1,mion
              x11=frac(kt1,kn1,ion)
              x21=frac(kt1+1,kn1,ion)
              x12=frac(kt1,kn1+1,ion)
              x22=frac(kt1+1,kn1+1,ion)
              x1221=x11*x21*x12*x22
              if(x1221.eq.0.) then
                  xx1=x11+at1*(x21-x11)
                  xx2=x12+at1*(x22-x12)
                  rrx=xx1+an1*(xx2-xx1)
              else
                  x11=dlog10(x11)
                  x21=dlog10(x21)
                  x12=dlog10(x12)
                  x22=dlog10(x22)
                  xx1=x11+at1*(x21-x11)
                  xx2=x12+at1*(x22-x12)
                  rrx=xx1+an1*(xx2-xx1)
                  rrx=exp(2.3025851*rrx)
              endif
              rrr(id,ion,iat)=rrx*abnd(iat)*dens(id)/wmm/ytot

!               print *,'abnd(iat)=',abnd(iat),iat
!               print *,'rrx =',rrx
!               print *,'dens =',dens
!               print *,'ytot =',ytot
!               pause

   50      continue
   30    continue
   20 CONTINUE
c
      return
      end subroutine
c 
c ******************************************************************
c
      subroutine fractn(iatnum)
c     =========================
c
C     INCLUDE 'IMPLIC.FOR'
	INCLUDE '../inc/IMPLIC.FOR'
      integer,intent(inout) :: iatnum
      parameter (mtemp=100,
     *           melec= 60,
     *           mion1=30,
     *           mdat = 17)
      parameter (inp=71)
      dimension frac0(-1:mion1),ioo(-1:mion1),idat(mion1)
      dimension gg(mion1,mdat),g0(mion1),z0(-1:mion1)
      dimension uu(mion1,mdat),u0(mion1)
      dimension u6(6),u7(7),u8(8),u10(10),u11(11)
      dimension u12(12),u13(13),u14(14),u16(16),u18(18),u20(20)
      dimension u24(24),u25(25),u26(26),u28(28)
      equivalence (u6(1),uu(1,3)),(u7(1),uu(1,4)),(u8(1),uu(1,5))
      equivalence (u10(1),uu(1,6)),(u11(1),uu(1,7)),(u12(1),uu(1,8))
      equivalence (u13(1),uu(1,9)),(u14(1),uu(1,10)),(u16(1),uu(1,11))
      equivalence (u18(1),uu(1,12)),(u20(1),uu(1,13)),(u24(1),uu(1,14))
      equivalence (u25(1),uu(1,15)),(u26(1),uu(1,16)),(u28(1),uu(1,17))
      common/fracop/frac(mtemp,melec,mion1),fracm(mtemp,melec),
     *              itemp(mtemp),ntt
      !character*30 filop
      data idat   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,
     *              7, 8, 9,10, 0,11, 0,12, 0,13,
     *              0, 0, 0,14,15,16, 0,17, 0, 0/ 
      data gg/2.,29*0.,2.,1.,28*0.,
     *        2.,1.,2.,1.,6.,9.,24*0.,2.,1.,2.,1.,6.,9.,4.,23*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,22*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,20*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,19*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,18*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,17*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,16*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,14*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,
     *        12*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *        6.,1.,2.,1.,10*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,
     *        6.,9.,4.,9.,6.,1.,10.,21.,28.,25.,6.,7.,6*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,7.,6.,5*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,30.,25.,4*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,28.,21.,10.,21.,0.,0./
      data uu(1,1)/109.6787/
      data uu(1,2)/198.3108/
      data uu(2,2)/438.9089/
      data u6/90.82,196.665,386.241,520.178,3162.395,3.952061/
      data u7/117.225,238.751,382.704,624.866,789.537,4452.758,5380.089/
      data u8/109.837,283.24,443.086,624.384,918.657,1114.008,5963.135,
     *        7028.393/
      data u10/173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,
     *         1928.462,9645.005,10986.876/
      data u11/41.449,381.395,577.8,797.8,1116.2,1388.5,1681.5,2130.8,
     *         2418.7,11817.061,13297.676/
      data u12/61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,
     *         2645.2,2964.4,14210.261,15829.951/
      data u13/48.278,151.86,229.446,967.8,1239.8,1536.3,1947.3,2295.4,
     *         2663.4,3214.8,3565.6,16825.022,18584.138/
      data u14/65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,
     *         2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63/
      data u16/83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,
     *         3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,
     *         28182.535/
      data u18/127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,
     *         3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,
     *         33237.173,35699.936/
      data u20/49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,
     *         1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,
     *         8770.,9338.,41366.,44177.41/
      data u24/54.576,132.966,249.7,396.5,560.2,731.02,1291.9,1490.,
     *         1688.,1971.,2184.,2404.,2862.,3098.52,8151.,8850.,
     *         9560.,10480.,11260.,12070.,13180.,13882.,60344.,63675.9/
      data u25/59.959,126.145,271.55,413.,584.,771.1,961.44,1569.,
     *         1789.,2003.,2307.,2536.,2771.,3250.,3509.82,9152.,
     *         9872.,10620.,11590.,12410.,13260.,14420.,15162.,
     *         65660.,69137.4/
      data u26/63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,
     *         1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,
     *         10180.,10985.,11850.,12708.,13620.,14510.,15797.,
     *         16500.,71203.,74829.6/
      data u28/61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,
     *         1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,
     *         4896.2,12430.,13290.,14160.,15280.,16220.,17190.,
     *         18510.,19351.,82984.,86909.4/
c
      if(idat(iatnum).eq.0) then
         write(6,600) iatnum
  600    format(' data for element no. ',i3,' do not exist')
         iatnum=-1
         return
      end if
c
      g0(iatnum+1)=1.
      do i=1,iatnum
        ig0=iatnum-i+1
        g0(ig0)=gg(i,idat(iatnum))
        u0(i)=uu(i,idat(iatnum))*1000.
      enddo
c
      if(iatnum.eq.1) open(inp,file='ioniz.dat',status='old')
      do 10 it=1,mtemp
         do 10 ie=1,melec
            fracm(it,ie)=0.
            do 10 ion=1,mion1
               frac(it,ie,ion)=0.
   10 continue
c
      read(inp,*)
      read(inp,*) it0,it1,itstp
      ntt=(it1-it0)/itstp+1
c
      do 100 it=1,ntt
         read(inp,*) itt,ie0,ie1,iestp
         itemp(it)=itt
         net=(ie1-ie0)/iestp+1
         t=exp(2.3025851*0.025*itt)
         safac0=sqrt(t)*t/2.07d-16
         tkcm=0.69496*t
         do 30 ie=1,net
            read(inp,601) iee,ion0,ion1,
     *                    (ioo(i),frac0(i),i=ion0,min(ion1,ion0+3))
            ane=exp(2.3025851*0.25*iee)
            safac=safac0/ane
            nio=ion1-ion0
            if(nio.ge.3) then
               nlin=nio/4
               do ilin=1,nlin
                  read(inp,602) (ioo(i),frac0(i),
     *                 i=ion0+4*ilin,min(ion1,ion0+4*ilin+3))
               end do
            end if
            ieind=iee/2
            do 20 ion=ion0,ion1
              if(ion.lt.iatnum) then
               if(ion.eq.ion0) then
                  z0(ion)=g0(iatnum-ion)
               else
                  z0(ion)=frac0(ion)/frac0(ion-1)*safac*z0(ion-1)
                  z0(ion)=z0(ion)*exp(-u0(iatnum-ion)/tkcm)
               endif
                  frac(it,ieind,iatnum-ion)=frac0(ion)/z0(ion)
              else
                  u0hm=6090.5
                  z0hm=frac0(ion)/frac0(ion-1)*safac
                  z0hm=z0hm*exp(-u0hm/tkcm)
                  fracm(it,ieind)=frac0(ion)/z0hm
              end if
   20       continue
c           write(6,603) it,ieind,t,ane
c           write(6,604) (frac(it,ieind,i),i=1,mion1)
   30    continue
  100 continue
  601 format(3i4,2x,4(i4,1x,e9.3))
  602 format(14x,4(i4,1x,e9.3))
  603 format(2i5,f15.2,1pe10.3)
  604 format(1p8e10.3)
      return      
      end subroutine
      end module
