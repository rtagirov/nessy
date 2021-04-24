      MODULE MOD_WRSTART

      CONTAINS

      SUBROUTINE WRSTART

      use MOD_DATOM
      use MOD_FGRID
      use MOD_GRADIFF
      use MOD_GREYM
      use MOD_PRIBLA
      use MOD_PRIDAT
      use MOD_PRIGH
      use MOD_PRIMOD
!      use MOD_WRITPOP
      use MOD_WRVEL
      use MOD_JSTART
      use MOD_PRICOMP
      use MOD_REBLANK
      use MOD_WRITMOD
      use MOD_TICTOC
      use MOD_chemeq

      use mod_decode
      use vardatom_full
      use vardatom_nlte
      use varhminus
      use geo_mesh
      use init_vel
      use common_block
      use file_operations
      use odf_table
      use phys

!     THIS PROGRAM IS TO INITIALIZE THE MODEL FILE FOR SUBSEQUENT
!     CALCULATION OF THE NON-LTE MULTI-LEVEL LINE FORMATION.
!     IT MAKES USE OF THE ATOMIC DATA (FILE DATOM)
!     AND THE FREQUENCY GRID (FILE fgrid.inp)
!     PRESENT VERSION: MODEL ATMOSPHERE OF HELIUM (CODE NR. "1") WITH
!                                          HYDROGEN         "2"
!     FOR IMPLEMENTATION OF ADDITIONAL ELEMENTS:
!     MODIFY SUBROUTINES  "DATOM", "DECSTAR"
!     INSERT CORRESPONDING ATOMIC DATA INTO SUBR. "COLLI", "PHOTOCS"

      IMPLICIT REAL*8(A - H, O - Z)

      COMMON /VELPAR/  VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE
      COMMON /COMTEFF/ TEFF, TMIN, TMODIFY, SPHERIC

      LOGICAL TTABLE, TPLOT, SPHERIC

      CHARACTER MODHEAD*104

      CHARACTER NAME*10, fstring*24
      integer   timer

      real*8 ATMEAN, AMU

      integer NA
      REAL*8 CSARR(5000,4)

      REAL*8, DIMENSION(:), ALLOCATABLE :: ELEC_CONC, HEAVY_ELEM_CONC

      REAL*8, DIMENSION(:), ALLOCATABLE :: VELO_NE, VELO_E

      CHARACTER(:), ALLOCATABLE :: AMF

      REAL*8 :: H

      real*8, dimension(30) :: atomic_mass ! atomic masses of the first 30 elements in atomic mass units

      data amu /1.660531d-24/

      data atomic_mass /1.0080,
     $                  4.0026,
     $                  6.9400,
     $                  9.0122,
     $                  10.810,
     $                  12.011,
     $                  14.007,
     $                  15.999,
     $                  18.998,
     $                  20.180,
     $                  22.990,
     $                  24.305,
     $                  26.982,
     $                  28.085,
     $                  30.974,
     $                  32.060,
     $                  35.450,
     $                  39.948,
     $                  39.098,
     $                  40.078,
     $                  44.956,
     $                  47.867,
     $                  50.942,
     $                  51.996,
     $                  54.938,
     $                  55.845,
     $                  58.933,
     $                  58.693,
     $                  63.546,
     $                  65.380/

      call FDATE(fstring)
      call TIC(timer)

!     READ ATOMIC DATA FROM FILE DATOM
      call datom('full', N, LEVEL,
     $           NCHARG,WEIGHT,ELEVEL,EION,MAINQN,
     $           EINST,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $           INDNUP,INDLOW,LASTIND,NATOM,
     $           ELEMENT,SYMBOL,NOM,KODAT,ATMASS,STAGE,NFIRST,
     $           NLAST,WAVARR,SIGARR,eleatnum,levatnum,NFDIM) ! NFDIM is known from varhminus module

!      call print_eion(N, level, levatnum, ncharg, eion, elevel, eion - elevel)

!      call print_istageinfo(natom, N, symbol, nfirst, nlast, ncharg)

      allocate(ABXYZ(NATOM))

!     DECODING INPUT DATA
      CALL DECSTAR(MODHEAD,FM,RSTAR,t_eff,glog,xmass,VDOP,TTABLE,TPLOT,NATOM,KODAT,IDAT,LBLANK,ATMEAN)

!      do j = 1, 30

!         print*, 'adundances', j, abxyz(j ), sum(abxyz)

!      enddo

!      stop 'stop abundances'

      apm = atomic_mass_unit * sum(abxyz * atomic_mass) ! average particle mass (cgs)

!     if PRINT DATOM option in CARDS is set, printout the atomic data
      IF (IDAT.EQ.1)
     $CALL PRIDAT(N,LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,
     $            KODAT,ALPHA,SEXPO,AGAUNT,COCO,KEYCOL,ALTESUM,
     $            NATOM,ELEMENT,NOM,ABXYZ,ATMASS)

      !***  PRINTOUT OF THE CHEMICAL COMPOSITION
      CALL PRICOMP(N,EINST,NCHARG,NOM,NATOM,ABXYZ,ATMASS,
     $             STAGE,NFIRST,NLAST,ELEMENT,SYMBOL,LASTIND,
     $             INDLOW,INDNUP)

!     GENERATION OF THE CONTINUOUS FREQUENCY GRID
      call fgrid(nfdim,nf,xlambda,fweight,akey,nom,symbol,natom,n,ncharg,elevel,eion,einst)

      call geomesh(radius, entot, T, P, Z, rstar, amu, atmean, ND, NP)

!     read odf.table and odf.table.grid files
!     odf.table.grid --- grid of temperature and pressures
!     at which the opacities in the odf.table file are given
      if (odf_from_table) call read_odf_table()

!     INITIALISATION OF THE VELOCITY-FIELD PARAMETERS
      call initvel(maxval(radius), t_eff, glog, rstar, xmass)

      allocate(xjc2(ND, NF))
      allocate(EDDI(3, ND))
      allocate(EDDARR(3, ND, NF))
      allocate(TAUROSS(ND))
      allocate(RNE(ND))
      allocate(VELO(ND))
      allocate(GRADI(ND))
      allocate(EMFLUX(NF))
      allocate(POPNUM(ND, N))
      allocate(HTOT(ND))
      allocate(GTOT(ND))
      allocate(XTOT(ND))
      allocate(ETOT(ND))
      allocate(POP1(ND, N))
      allocate(POP2(ND, N))
      allocate(POP3(ND, N))

      allocate(scafac(ND, NF))
      allocate(absfac(ND, NF))

      allocate(ABXYZn(NATOM, ND))

      allocate(U(ND, NP))
      allocate(VJL(NP, ND))
      allocate(VL(NP))
      allocate(HNU(ND))

      allocate(mainpro(ND))
      allocate(mainlev(ND))

      call mol_ab(ABXYZn, ABXYZ, SYMBOL, ENTOT, T, NATOM, ND)

      call mark_nlte()

      print*, 'after mark_nlte'

!     RINAT TAGIROV:
!     Calculation or read-out of the velocity field.
!     In case of the calculation
!     the expansion law is regulated with the boundary
!     condition parameters VMIN and VFINAL
!     which are set in the CARDS file.
!     Using these parameters INITVEL subroutine
!     generates VPAR1 and VPAR2 parameters
!     entering the function WRVEL below which calculates the ultimate velocities.

!     In case of the read-out the velocity is read from the VEL_FIELD_FILE set in FILE_OPERATIONS.FOR
!     and if there is a local extremum in the input data...
!     (local maximum in case of negative velocities, local minimum in case of positive velocities;
!     at least it was like that in the 3D simulation data I worked with and => EXTRAP_VEL_FIELD subroutine
!     below works only for this two cases, however the rest of the code can't handle the negative case
!     because apparently the frequency initial condition in it is not suited for contracting flows,
!     see Mihalas, Kunasz & Hummer, ApJ, 202:465-489, 1975, Appendix B for details about various velocity laws and the
!     type of radiative transfer scheme implemented in the code)
!     ...the law gets extrapolated from the point of the extremum
!     to the innermost point yelding therefore a monotonically increasing/decreasing function.
!     The height grid in the VEL_FIELD_FILE has to be the same as in the atmosphere model file atm.inp.
!     The TABLE string in CARDS file was used before to
!     control the calculation/read-out option but is obsolete now (it is still in the CARDS file though).
!     The logical variable VEL_FIELD_FROM_FILE is declared in comblock.for and set in hminus.for.
!     The units of velocity in the VEL_FIELD_FILE are km/s, height is in km.
!     The first column is height, the second is velocity.

      if (vel_field_from_file) then

         allocate(velo_ne(nd))
         allocate(velo_e(nd))

         open(unit = 1832, file = vel_field_file, action = 'read')

         do i = 1, ND; read(1832, *) h, velo_ne(i); enddo

         close(1832)

!        Extrapolation of the velocity law. If there is no local minimum/maximum then VELO_E = VELO_NE.
         velo_e = extrap_vel_field(velo_ne(1 : ND), ND)

         velo(1 : ND) = velo_e(1 : ND)

      else

         do L = 1, ND; velo(L) = wrvel(radius(L)); enddo

      endif

      call gradiff(ND, velo, gradi, radius)
 
C***  STAPEL: NUMBER OF FREE ELECTRONS PER ATOM
C***  S T A R T   A P P R O X I M A T I O N

      STAPEL = 0.0d0

      DO NA = 1, NATOM; STAPEL = STAPEL + ABXYZ(NA) * (STAGE(NA) - 1.); ENDDO

      RNE(1 : ND) = STAPEL

C***  Read Line-blanketing table
      IF (LBLANK .NE. 0) LBLANK = -2

      CALL REBLANK (LBLANK,NF,XLAMBDA,ND,ENTOT,RNE,SCAFAC,ABSFAC)

      IF (ABS(LBLANK) .EQ. 2) CALL PRIBLA (LBLANK,ENTOT,ND,XLAMBDA,NF,JOBNUM,MODHEAD,SCAFAC,ABSFAC)
 
C***  TEMPERATURE STRATIFICATION AND INITIAL POPNUMBERS (LTE)

      CALL GREYM(ND,T,RADIUS,XLAMBDA,FWEIGHT,NF,NFDIM,ENTOT,RNE,RSTAR,
     $           ALPHA,SEXPO,AGAUNT,POPNUM,TAUROSS,R23,TTABLE,N,
     $           LEVEL,NCHARG,WEIGHT,ELEVEL,EION,EINST,KODAT,
     $           NOM,NFIRST,NLAST,NATOM,WAVARR,SIGARR)

      CALL PRIMOD(ND,RADIUS,RSTAR,ENTOT,T,VELO,GRADI,NP,MODHEAD,JOBNUM,TTABLE,TAUROSS,R23)

!=================================================================
! CONSTANT ELECTRON CONCENTRATION RUN AND/OR LTE RUN

      IF (LTE_RUN) THEN

         ALLOCATE(ELEC_CONC(ND))
         ALLOCATE(HEAVY_ELEM_CONC(ND))

         ELEC_CONC =       read_atm_file_col(3)
         HEAVY_ELEM_CONC = read_atm_file_col(4)

         RNE(1 : ND) = ELEC_CONC(1 : ND) / HEAVY_ELEM_CONC(1 : ND)

         DEALLOCATE(ELEC_CONC)
         DEALLOCATE(HEAVY_ELEM_CONC)

      ENDIF

!=================================================================
 
!     MODFILE: MASS STORAGE FILE IN NAME-INDEX MODE
      TOTOUT =       0.0d0
      TOTIN =        0.0d0

      POP1(:, :) =   0.0d0
      POP2(:, :) =   0.0d0
      POP3(:, :) =   0.0d0

      HTOT(1 : ND) =   0.0d0
      GTOT(1 : ND) =   0.0d0
      XTOT(1 : ND) =   0.0d0
      ETOT(1 : ND) =   0.0d0

      EMFLUX(1 : NF) = 0.0d0

      if (allocated(wcharm)) deallocate(wcharm)

      allocate(wcharm(ND, NF))

      WCHARM(1 : ND, 1 : NF) = 0.0d0

      IFL = 3; open(IFL, file = 'MODFILE', STATUS = 'UNKNOWN')

      JOBNUM = 0

      CALL WRITMOD(IFL,N,ND,TEFF,RADIUS,NP,P,Z,ENTOT,VELO,GRADI,RSTAR,VDOP,NF,
     $             XLAMBDA(1 : NF),FWEIGHT(1 : NF),AKEY(1 : NF),
     $             ABXYZ,NATOM,MODHEAD,JOBNUM)

      CLOSE(ifl)

!      ifl = 3; open(ifl, file = 'POPNUM', status = 'unknown')

!      call writpop(ifl, T, popnum, pop1, pop2, pop3, rne, n, nd, modhead, jobnum)

!      close(ifl)

!     START APPROXIMATION FOR THE RADIATION FIELD
!     JSTART writes the files RADIOC and RADIOL

      EDDI(1 : 3, 1 : ND) = 0.0d0

      CALL JSTART(NF, lastind_nlte, XLAMBDA(1 : NF), ND, T, XJC, XJL,
     $            HTOT, GTOT, XTOT, ETOT, EMFLUX, TOTIN, TOTOUT,
     $            ncharg_nlte, elevel_nlte, EDDI, WCHARM, nom_nlte, N_nlte, einst_nlte,
     $            MODHEAD, JOBNUM, TEFF)

      write(*,  *) 'WRSTART - ', fstring, ' run time: ', TOC(timer)

      open(78, file = 'MODHIST', status = 'unknown')
   
      write(78, *) 'WRSTART - ', fstring, ' run time: ', TOC(timer)
   
      close(78)

      return
 
      end subroutine

      subroutine print_eion(N, level, levatnum, ncharg, eion, elevel, E_ionization)

      implicit none

      integer, intent(in) :: N

      character (len = 10), dimension(N), intent(in) :: level

      integer, dimension(N), intent(in) :: levatnum, ncharg

      real*8,  dimension(N), intent(in) :: eion, elevel, E_ionization

      integer :: i

      open(unit = 2764, file = 'eion.out', action = 'write')

      do i = 1, N

         write(2764, '(i4,2x,A10,2x,i3,2x,i2,3(2x,F8.1)))') i, level(i), levatnum(i), ncharg(i),
     $                                                      eion(i), elevel(i), E_ionization(i)

      enddo

      close(2764)

      end subroutine

      subroutine print_istageinfo(natom, N, symbol, nfirst, nlast, ncharg)

!      use utils

      implicit none

      integer, intent(in)                               :: N, natom

      character (len = 2), intent(in), dimension(natom) :: symbol

      integer, intent(in), dimension(natom)             :: nfirst, nlast

      integer, intent(in), dimension(N)                 :: ncharg

      integer,             dimension(natom)             :: num_i_stages

      integer, allocatable, dimension(:)                :: num_stage_lev

      integer                                           :: un, i, sc, k, msc

      un = 1047; open(unit = un, file = 'istageinfo.out', action = 'write')

      do k = 1, natom

!     number of ionization stsages within each element
          num_i_stages(k) = 1

          do i = nfirst(k), nlast(k) - 1

              if (ncharg(i + 1) /= ncharg(i)) num_i_stages(k) = num_i_stages(k) + 1
    
          enddo

      enddo

      msc = maxval(num_i_stages)

      allocate(num_stage_lev(msc))

      do k = 1, natom

         num_stage_lev(:) = 0

         sc = 1

         do i = nfirst(k), nlast(k)

            num_stage_lev(sc) = num_stage_lev(sc) + 1

            if (i == nlast(k)) exit

            if (ncharg(i + 1) /= ncharg(i)) sc = sc + 1

         enddo

         write(un, '(A2,2x,i2,$)') symbol(k), num_i_stages(k)

         do i = 1, msc

            if (i /= msc) write(un, '(2x,i2,$)') num_stage_lev(i)
            if (i == msc) write(un, '(2x,i2)')   num_stage_lev(i)

         enddo

      enddo

      deallocate(num_stage_lev)

      close(un)

      return

      end subroutine print_istageinfo

      FUNCTION EXTRAP_VEL_FIELD(VEL, ND) RESULT(VEL_E)

      USE FILE_OPERATIONS
      USE COMMON_BLOCK

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::               ND

      REAL*8, DIMENSION(ND), INTENT(IN) :: VEL
      REAL*8, DIMENSION(ND) ::             VEL_E

      INTEGER ::                           I, ML, NI, EXTRLOC

      REAL*8 ::                            A, B

      VEL_E(1 : ND) = VEL(1 : ND)

      IF (VEL(1) .GT. 0) EXTRLOC = MINLOC(VEL, 1)
      IF (VEL(1) .LT. 0) EXTRLOC = MAXLOC(VEL, 1)

!     If there is no local minimum/maximum then no extrapolation is needed
      IF (EXTRLOC .EQ. 1) THEN

         WRITE(*, '(/,A,1x,A,1x,A,/)') 'NO LOCAL MINIMUM/MAXIMUM WAS FOUND IN', VEL_FIELD_FILE,
     $                                 '=> NO EXTRAPOLATION HAS BEEN PERFORMED.
     $                                  PROCESSING THE VELOCITY FIELD AS IT IS.'

         RETURN

      ENDIF

      ML = EXTRLOC; NI = ML - 1

      A = (VEL(ML) - VEL(NI)) / (HEIGHT(ML) - HEIGHT(NI))

      B = VEL(ML) - A * HEIGHT(ML)

      DO I = ML + 1, ND, 1; VEL_E(I) = A * HEIGHT(I) + B; ENDDO

      RETURN

      END FUNCTION EXTRAP_VEL_FIELD

      subroutine mark_nlte()

      use vardatom_full
      use vardatom_lte
      use vardatom_nlte
      use varhminus
      use common_block
!      use file_operations
      use mod_datom

      implicit none

      integer :: i, j, L, k

      real*8  :: eground

      call datom('nlte',
     $           N_nlte,
     $           level_nlte,
     $           ncharg_nlte,
     $           weight_nlte,
     $           elevel_nlte,
     $           eion_nlte,
     $           mainqn_nlte,
     $           einst_nlte,
     $           alpha_nlte,
     $           sexpo_nlte,
     $           agaunt_nlte,
     $           coco_nlte,
     $           keycol_nlte,
     $           altesum_nlte,
     $           indnup_nlte,
     $           indlow_nlte,
     $           lastind_nlte,
     $           natom_nlte,
     $           element_nlte,
     $           symbol_nlte,
     $           nom_nlte,
     $           kodat_nlte,
     $           atmass_nlte,
     $           stage_nlte,
     $           nfirst_nlte,
     $           nlast_nlte,
     $           wavarr_nlte,
     $           sigarr_nlte,
     $           eleatnum_nlte,
     $           levatnum_nlte,
     $           nfdim) ! NFDIM is known from varhminus module

      natom_lte = natom - natom_nlte

      allocate(nlte_lev(N))

      allocate(idx_orig(N_nlte))

      allocate(nlte_ele(NATOM))

      allocate(abxyz_nlte(natom_nlte))

      allocate(ABXYZn_nlte(natom_nlte, dpn))

      if (natom_lte /= 0) then

          allocate(abxyz_lte(natom_lte))

          allocate(abxyzn_lte(natom_lte, dpn))

          allocate(eleatnum_lte(natom_lte))
          allocate(eleisnum_lte(natom_lte))

      endif

      nlte_lev(1 : N) =     .false.

      nlte_ele(1 : NATOM) = .false.

      do i = 1, N

         do j = 1, N_nlte

            if (level_nlte(j) .eq. level(i)) then

                nlte_lev(i) = .true.

                idx_orig(j) = i

                cycle

            endif

         enddo

      enddo

      do i = 1, NATOM

         do j = 1, natom_nlte

            if (element_nlte(j) .eq. element(i)) then

                nlte_ele(i) = .true.

                abxyz_nlte(j) = abxyz(i)

                do L = 1, dpn; ABXYZn_nlte(j, L) = ABXYZn(i, L); enddo

                cycle

            endif

         enddo

      enddo

      if (natom_lte /= 0) then

          j = 1

          do i = 1, NATOM

             if (.not. nlte_ele(i)) then

                abxyz_lte(j) = abxyz(i)

                do L = 1, dpn; abxyzn_lte(j, L) = ABXYZn(i, L); enddo

                j = j + 1

             endif

          enddo

          j = 1

          do i = 1, NATOM

             if (.not. nlte_ele(i)) then

                eleatnum_lte(j) = eleatnum(i)

                j = j + 1

             endif

          enddo

          eleisnum_lte(:) = 0

          do j = 1, natom_lte

             do i = 1, N

                if (levatnum(i) .eq. eleatnum_lte(j)) then

                   if (i .ne. 1 .and. ncharg(i) .eq. ncharg(i - 1)) cycle

                   eleisnum_lte(j) = eleisnum_lte(j) + 1 ! number of ionization stages in lte element j

                endif

             enddo

          enddo

      endif

      do i = 1, N

         write(*, '(A,2x,i4,2x,A,2x,L)') 'mark_nlte, levels:', i, level(i), nlte_lev(i)

      enddo

      j = 1; k = 1

      do i = 1, NATOM

         if (.not. nlte_ele(i)) then 

            write(*, '(A,2x,i4,2x,A,2x,L,2x,ES9.3,2x,A10,2x,ES9.3)')
     $      'mark_nlte, abundances:', i, element(i), nlte_ele(i), ABXYZ(i), '         ', abxyz_lte(k)

            k = k + 1

         else

            write(*, '(A,2x,i4,2x,A,2x,L,2(2x,ES9.3))')
     $      'mark_nlte, abundances:', i, element(i), nlte_ele(i), ABXYZ(i), abxyz_nlte(j)

            j = j + 1

         endif

      enddo

      if (natom_lte /= 0) then

          lis_num = sum(eleisnum_lte) ! number of LTE Ionization Stages (LIS)

          allocate(lis_name(lis_num)) ! name of each ionization stage
          allocate(lis_anum(lis_num)) ! atomic number of each ionization stage
          allocate(lis_cnum(lis_num)) ! charge number of each ionization stage

          j = 1

          do i = 1, N

             if (nlte_lev(i)) cycle

             if (i .ne. 1 .and. ncharg(i) .ne. ncharg(i - 1)) then

                lis_name(j) = level(i)

                lis_anum(j) = levatnum(i)

                lis_cnum(j) = ncharg(i)

                j = j + 1

             endif

          enddo

          allocate(lis_lnum(lis_num)) ! number of levels in each ionization stage

          lis_lnum(1 : lis_num) = 1

          j = 0

          do i = 1, N

             if (nlte_lev(i)) cycle

             if (i .ne. 1 .and. ncharg(i) .ne. ncharg(i - 1)) j = j + 1

             if (i .ne. 1 .and. ncharg(i) .eq. ncharg(i - 1)) lis_lnum(j) = lis_lnum(j) + 1

          enddo

          allocate(lis_weight(lis_num, maxval(lis_lnum)))

          allocate(lis_levien(lis_num, maxval(lis_lnum)))

          lis_weight(:, :) = 0

          lis_levien(:, :) = 0.0d0

          j = 0

          do i = 1, N

             if (nlte_lev(i)) cycle

             if (i .ne. 1 .and. ncharg(i) .ne. ncharg(i - 1)) then

                j = j + 1 ! lte ionization stage counter

                k = 1     ! level counter within the lte ionization stage

                lis_weight(j, k) = weight(i)

                lis_levien(j, k) = eion(i)

                eground          = elevel(i)

             endif

             if (i .ne. 1 .and. ncharg(i) .eq. ncharg(i - 1)) then

                k = k + 1

!               statistical weight of level k within the lte ionization stage j
                lis_weight(j, k) = weight(i)

 !              LEVel Ionization ENergy (LEVIEN) or the ionization energy from excited state
                lis_levien(j, k) = eion(i) - elevel(i) + eground

             endif

          enddo

          do i = 1, lis_num

             write(*, '(I2,2x,A10,3(2x,I2))'), i, lis_name(i), lis_anum(i), lis_lnum(i), lis_cnum(i)

          enddo

          write(*, '(/)')

          do i = 1, lis_num

             write(*, '($,I2,2x,A10,2x,I2,A))'), i, lis_name(i), lis_lnum(i), '        '

             do j = 1, maxval(lis_lnum)

                if (j /= maxval(lis_lnum)) write(*, '($,2x,I2)'),   lis_weight(i, j)

                if (j == maxval(lis_lnum)) write(*, '($,2x,I2,/)'), lis_weight(i, j)

             enddo

          enddo

          write(*, '(/)')

          do i = 1, lis_num

             write(*, '($,I2,2x,A10,2x,I2,A))'), i, lis_name(i), lis_lnum(i), '        '

             do j = 1, maxval(lis_lnum)

                if (j /= maxval(lis_lnum)) write(*, '($,2x,F8.1)'),   lis_levien(i, j)

                if (j == maxval(lis_lnum)) write(*, '($,2x,F8.1,/)'), lis_levien(i, j)

             enddo

          enddo

      endif

      end subroutine

      end module
