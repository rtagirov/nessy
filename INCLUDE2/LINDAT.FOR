c      PARAMETER (MLIN0   = 100000,
      integer,PARAMETER :: MLIN0   = 500000,
     *           MGRIEM  =     10,
c     *           MNLT    = 100000,
     *           MNLT    = 500000,
     *           MSPHE2  =     20,
     *           MLIN    =    500,
     *           MPRF    =  MLIN0
C
      real*8 :: FREQ0,EXCL0,GF0,EXTIN,BNUL
      integer:: INDAT,INDNLT,ILOWN,IUPN,IPOTL,IPOTU,INDLIN,INDLIP
      integer:: NLIN0,NNLT,NLIN,NGRIEM
      COMMON/LINTOT/FREQ0(MLIN0),
     *              EXCL0(MLIN0),
     *              GF0(MLIN0),
     *              EXTIN(MLIN0),
     *              BNUL(MLIN0),
     *              INDAT(MLIN0),
     *              INDNLT(MLIN0),
     *              ILOWN(MLIN0),
     *              IUPN(MLIN0),
     *              IPOTL(MLIN0),
     *              IPOTU(MLIN0),
     *              INDLIN(MLIN),
     *              INDLIP(MLIN),
     *              NLIN0,NLIN,
     *              NNLT,NGRIEM
C
      real*8 :: GAMR0,GS0,GW0,WGR0
      integer:: IPRF0,ISPRF,IGRIEM,ISP0,NSP
      COMMON/LINPRF/GAMR0(MPRF),
     *              GS0(MPRF),
     *              GW0(MPRF),
     *              WGR0(4,MGRIEM),
     *              IPRF0(MPRF),
     *              ISPRF(MPRF),
     *              IGRIEM(MPRF),
     *              ISP0(MSPHE2),NSP
C
      real*8 :: ABCENT,SLIN
      COMMON/LINNLT/ABCENT(MNLT,MDEPTH),
     *              SLIN(MNLT,MDEPTH)
C
      real*8 PLAN,STIM,EXHK
      COMMON/LINDEP/PLAN(MDEPTH),
     *              STIM(MDEPTH),
     *              EXHK(MDEPTH)
