MODULE icb_oce
  USE par_oce
  USE lib_mpp
  IMPLICIT NONE
  PUBLIC
  PUBLIC :: icb_alloc
  INTEGER, PUBLIC, PARAMETER :: nclasses = 10
  INTEGER, PUBLIC, PARAMETER :: nkounts = 3
  TYPE, PUBLIC :: icebergs_gridded
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: calving
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: calving_hflx
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: floating_melt
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: maxclass
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: tmp
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: stored_ice
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: stored_heat
  END TYPE icebergs_gridded
  TYPE, PUBLIC :: point
    INTEGER :: year
    REAL(KIND = wp) :: xi, yj
    REAL(KIND = wp) :: e1, e2
    REAL(KIND = wp) :: lon, lat, day
    REAL(KIND = wp) :: mass, thickness, width, length, uvel, vvel
    REAL(KIND = wp) :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi
    REAL(KIND = wp) :: mass_of_bits, heat_density
  END TYPE point
  TYPE, PUBLIC :: iceberg
    TYPE(iceberg), POINTER :: prev => NULL(), next => NULL()
    INTEGER, DIMENSION(nkounts) :: number
    REAL(KIND = wp) :: mass_scaling
    TYPE(point), POINTER :: current_point => NULL()
  END TYPE iceberg
  !TYPE(icebergs_gridded), POINTER :: berg_grid
  TYPE(icebergs_gridded) :: berg_grid ! CDe remove POINTER
  TYPE(iceberg), POINTER :: first_berg => NULL()
  REAL(KIND = wp) :: berg_dt
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE :: first_width, first_length
  LOGICAL :: l_restarted_bergs = .FALSE.
  REAL(KIND = wp), DIMENSION(nclasses), PUBLIC :: class_num = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: uo_e, vo_e
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: ff_e, tt_e, fr_e
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: ua_e, va_e
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: ssh_e
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :), ALLOCATABLE :: tmask_e, umask_e, vmask_e
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_rho_ice = 916.7_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_rho_water = 999.8_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_rho_air = 1.1_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_rho_seawater = 1025._wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_Cd_av = 1.3_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_Cd_ah = 0.0055_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_Cd_wv = 0.9_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_Cd_wh = 0.0012_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: pp_Cd_iv = 0.9_wp
  LOGICAL, PUBLIC :: ln_bergdia
  INTEGER, PUBLIC :: nn_verbose_level
  INTEGER, PUBLIC :: nn_test_icebergs
  REAL(KIND = wp), PUBLIC, DIMENSION(4) :: rn_test_box
  LOGICAL, PUBLIC :: ln_use_calving
  INTEGER, PUBLIC :: nn_sample_rate
  INTEGER, PUBLIC :: nn_verbose_write
  REAL(KIND = wp), PUBLIC :: rn_rho_bergs
  REAL(KIND = wp), PUBLIC :: rn_LoW_ratio
  REAL(KIND = wp), PUBLIC :: rn_bits_erosion_fraction
  REAL(KIND = wp), PUBLIC :: rn_sicn_shift
  LOGICAL, PUBLIC :: ln_operator_splitting
  LOGICAL, PUBLIC :: ln_passive_mode
  LOGICAL, PUBLIC :: ln_time_average_weight
  REAL(KIND = wp), PUBLIC :: rn_speed_limit
  REAL(KIND = wp), DIMENSION(nclasses), PUBLIC :: rn_initial_mass
  REAL(KIND = wp), DIMENSION(nclasses), PUBLIC :: rn_distribution
  REAL(KIND = wp), DIMENSION(nclasses), PUBLIC :: rn_mass_scaling
  REAL(KIND = wp), DIMENSION(nclasses), PUBLIC :: rn_initial_thickness
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: src_calving, src_calving_hflx
  INTEGER, PUBLIC, SAVE :: numicb
  INTEGER, PUBLIC, SAVE, DIMENSION(nkounts) :: num_bergs
  INTEGER, PUBLIC, SAVE :: nicbdi, nicbei, nicbdj, nicbej
  REAL(KIND = wp), PUBLIC, SAVE :: ricb_left, ricb_right
  INTEGER, PUBLIC, SAVE :: nicbpack
  INTEGER, PUBLIC, SAVE :: nktberg, nknberg
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nicbfldpts
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nicbflddest
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nicbfldproc
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nicbfldnsend
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nicbfldexpect
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: nicbfldreq
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: griddata
  CONTAINS
  INTEGER FUNCTION icb_alloc()
    INTEGER :: ill
    icb_alloc = 0
    !ALLOCATE(berg_grid, STAT = ill) ! CDe
    !icb_alloc = icb_alloc + ill ! CDe
    ALLOCATE(berg_grid % calving(jpi, jpj), berg_grid % calving_hflx(jpi, jpj), berg_grid % stored_heat(jpi, jpj), berg_grid % &
&floating_melt(jpi, jpj), berg_grid % maxclass(jpi, jpj), berg_grid % stored_ice(jpi, jpj, nclasses), berg_grid % tmp(jpi, jpj), &
&STAT = ill)
    icb_alloc = icb_alloc + ill
    ALLOCATE(uo_e(0 : jpi + 1, 0 : jpj + 1), ua_e(0 : jpi + 1, 0 : jpj + 1), vo_e(0 : jpi + 1, 0 : jpj + 1), va_e(0 : jpi + 1, 0 : &
&jpj + 1), ff_e(0 : jpi + 1, 0 : jpj + 1), fr_e(0 : jpi + 1, 0 : jpj + 1), tt_e(0 : jpi + 1, 0 : jpj + 1), ssh_e(0 : jpi + 1, 0 : &
&jpj + 1), first_width(nclasses), first_length(nclasses), src_calving(jpi, jpj), src_calving_hflx(jpi, jpj), STAT = ill)
    icb_alloc = icb_alloc + ill
    ALLOCATE(tmask_e(0 : jpi + 1, 0 : jpj + 1), umask_e(0 : jpi + 1, 0 : jpj + 1), vmask_e(0 : jpi + 1, 0 : jpj + 1), STAT = ill)
    icb_alloc = icb_alloc + ill
    ALLOCATE(nicbfldpts(jpi), nicbflddest(jpi), nicbfldproc(jpni), nicbfldnsend(jpni), nicbfldexpect(jpni), nicbfldreq(jpni), STAT &
&= ill)
    icb_alloc = icb_alloc + ill
    ALLOCATE(griddata(jpi, jpj, 1), STAT = ill)
    icb_alloc = icb_alloc + ill
    CALL mpp_sum('icb_oce', icb_alloc)
    IF (icb_alloc > 0) CALL ctl_warn('icb_alloc: allocation of arrays failed')
  END FUNCTION icb_alloc
END MODULE icb_oce
