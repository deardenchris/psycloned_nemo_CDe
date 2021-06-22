MODULE sbc_ice
  USE par_oce
  USE sbc_oce
  USE ice
  USE lib_mpp
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_ice_alloc
  LOGICAL, PUBLIC, PARAMETER :: lk_si3 = .TRUE.
  LOGICAL, PUBLIC, PARAMETER :: lk_cice = .FALSE.
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qns_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qsr_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qla_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: dqla_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: dqns_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: tn_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: alb_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qml_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qcn_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qtr_ice_top
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: utau_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: vtau_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: emp_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: topmelt
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: botmelt
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: evap_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: devap_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qns_oce
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qsr_oce
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qemp_oce
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qemp_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qevap_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qprec_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: emp_oce
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wndm_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sstfrz
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: tsfc_ice
  REAL(KIND = wp), PUBLIC, SAVE :: cldf_ice = 0.81
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: snwice_mass
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: snwice_mass_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: snwice_fmass
  CONTAINS
  INTEGER FUNCTION sbc_ice_alloc()
    INTEGER :: ierr(4)
    !$ACC KERNELS
    ierr(:) = 0
    !$ACC END KERNELS
    ALLOCATE(snwice_mass(jpi, jpj), snwice_mass_b(jpi, jpj), snwice_fmass(jpi, jpj), STAT = ierr(1))
    ALLOCATE(qns_ice(jpi, jpj, jpl), qsr_ice(jpi, jpj, jpl), qla_ice(jpi, jpj, jpl), dqla_ice(jpi, jpj, jpl), dqns_ice(jpi, jpj, jpl), tn_ice(jpi, jpj, jpl), alb_ice(jpi, jpj, jpl), qml_ice(jpi, jpj, jpl), qcn_ice(jpi, jpj, jpl), qtr_ice_top(jpi, jpj, jpl), utau_ice(jpi, jpj), vtau_ice(jpi, jpj), wndm_ice(jpi, jpj), evap_ice(jpi, jpj, jpl), devap_ice(jpi, jpj, jpl), qprec_ice(jpi, jpj), qemp_ice(jpi, jpj), qevap_ice(jpi, jpj, jpl), qemp_oce(jpi, jpj), qns_oce(jpi, jpj), qsr_oce(jpi, jpj), emp_oce(jpi, jpj), emp_ice(jpi, jpj), tsfc_ice(jpi, jpj, jpl), sstfrz(jpi, jpj), STAT = ierr(2))
    sbc_ice_alloc = MAXVAL(ierr)
    CALL mpp_sum('sbc_ice', sbc_ice_alloc)
    IF (sbc_ice_alloc > 0) CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
  END FUNCTION sbc_ice_alloc
END MODULE sbc_ice