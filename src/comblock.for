      module common_block

!     hminus variables

      logical ::                                 damp_acc =        .true.
      logical ::                                 rayleigh =        .false.
      logical ::                                 each_ali =        .false.
      logical ::                                 odf_from_table =  .true. ! if .false. then read ODF from *.lbkg files
      logical ::                                 full_conv_print = .false.
      logical ::                                 hm_lte =          .false. ! if true then hminus_lte is called by hminus
      logical ::                                 lte_run
      logical ::                                 vel_field_from_file
      logical ::                                 oldstart
      logical ::                                 lbkg ! keyword for taking into account opacity distribution function (ODF)

      integer ::                                 ndpmin ! depth point of the temperature minimum
      integer ::                                 dpn
      integer ::                                 xlbkg1, xlbkg2 ! wavelength range for the ODF
      integer ::                                 lambda_iter
      integer ::                                 nline
      integer ::                                 negintl

      real*8 ::                                  cormax_elec

      logical, dimension(:),      allocatable :: nofile
      logical, dimension(:, :),   allocatable :: damp_line, damp_cont

      character(:),               allocatable :: nrrm_file, ncrm_file, ntrm_file

      real*8, dimension(:),       allocatable :: height
      real*8, dimension(:),       allocatable :: ABXYZ
      real*8, dimension(:, :),    allocatable :: ABXYZn
      real*8, dimension(:, :),    allocatable :: llo ! the Local approximate lambda-Operator for all lines (Overall)
      real*8, dimension(:, :),    allocatable :: tau_line, tau_cont
      real*8, dimension(:, :),    allocatable :: xjl_lte
      real*8, dimension(:, :),    allocatable :: xjc_lte
      real*8, dimension(:, :),    allocatable :: pop_lte
      real*8, dimension(:, :, :), allocatable :: arr, acr, rbr
      real*8, dimension(:, :, :), allocatable :: arr_lte

!     nlte marking variables

      logical, dimension(:),      allocatable :: nlte_lev, nlte_ele
      integer, dimension(:),      allocatable :: idx_orig
      real*8,  dimension(:),      allocatable :: ABXYZ_nlte
      real*8,  dimension(:, :),   allocatable :: ABXYZn_nlte

!     variables and parameters for the odf_table module (see also odf_from_table above)

      integer ::                                               numt, nump

      integer, parameter ::                                    nbins = 90
      integer, parameter ::                                    nsubbins = 10

      real*8,            allocatable, dimension(:) ::          wvlgrid
      real*8,            allocatable, dimension(:) ::          tabt, tabp
      integer(kind = 2), allocatable, dimension(:, :, :, :) :: odf

!     fioss variables

      logical ::                               lopa = .false.

      integer ::                               numtra
      integer ::                               ntc

      integer, dimension(:),    allocatable :: nti
      integer, dimension(:),    allocatable :: elid, charge, swl, swu
                                              
      real*8,  dimension(:),    allocatable :: wavtra, aul
      real*8,  dimension(:, :), allocatable :: ntpl, ntpu, ntdl, ntdu

      end module

!      real*8,            allocatable, dimension(:)          :: c1, c2, c3, c4
!      integer,           allocatable, dimension(:)          :: idxp, idxt

!      real*8, allocatable, dimension(:) :: rne1
