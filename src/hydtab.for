      MODULE MOD_HYDTAB

      real*8 ::             WLINE(4, 22)

      integer, parameter :: MLINH = 4 * 22 ! maxLine
      integer, parameter :: MHWL = 80      ! maxWavelength
      integer, parameter :: MHE = 20       ! max#electrons
      integer, parameter :: MHT = 7        ! max#Temp

      real*8 ::             WL(MHWL, MLINH), XT(MHT, MLINH)
      real*8 ::             XNE(MHE, MLINH), PRF(MHWL, MHT, MHE)

      integer ::            NWLH(MLINH), NTH(MLINH), NEH(MLINH)
      integer ::            ILIN0(4, 22)

      END MODULE
