MODULE trdmxl_oce
  USE par_oce
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trdmxl_oce_alloc
  INTEGER, PUBLIC, PARAMETER :: jpltrd = 12
  INTEGER, PUBLIC :: jpktrd
  INTEGER, PUBLIC, PARAMETER :: jpmxl_xad = 1
  INTEGER, PUBLIC, PARAMETER :: jpmxl_yad = 2
  INTEGER, PUBLIC, PARAMETER :: jpmxl_zad = 3
  INTEGER, PUBLIC, PARAMETER :: jpmxl_ldf = 4
  INTEGER, PUBLIC, PARAMETER :: jpmxl_zdf = 5
  INTEGER, PUBLIC, PARAMETER :: jpmxl_npc = 6
  INTEGER, PUBLIC, PARAMETER :: jpmxl_bbc = 7
  INTEGER, PUBLIC, PARAMETER :: jpmxl_bbl = 8
  INTEGER, PUBLIC, PARAMETER :: jpmxl_for = 9
  INTEGER, PUBLIC, PARAMETER :: jpmxl_dmp = 10
  INTEGER, PUBLIC, PARAMETER :: jpmxl_zdfp = 11
  INTEGER, PUBLIC, PARAMETER :: jpmxl_atf = 12
  INTEGER, PUBLIC :: nn_ctls = 0
  REAL(KIND = wp), PUBLIC :: rn_rho_c = 0.01
  REAL(KIND = wp), PUBLIC :: rn_ucf = 1.
  CHARACTER(LEN = 32), PUBLIC :: cn_trdrst_in = "restart_mxl"
  CHARACTER(LEN = 32), PUBLIC :: cn_trdrst_out = "restart_mxl"
  LOGICAL, PUBLIC :: ln_trdmxl_instant = .FALSE.
  LOGICAL, PUBLIC :: ln_trdmxl_restart = .FALSE.
  CHARACTER(LEN = 80), PUBLIC :: clname, ctrd(jpltrd + 1, 2)
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: nmxl
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: nbol
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: wkx
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hmxl, tml, sml, tmlb, smlb, tmlbb, smlbb, tmlbn, smlbn, tmltrdm, &
&smltrdm, tml_sum, tml_sumb, tmltrd_atf_sumb, sml_sum, sml_sumb, smltrd_atf_sumb, hmxl_sum, hmxlbn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tmlatfb, tmlatfn, smlatfb, smlatfn, tmlatfm, smlatfm
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :, :) :: tmltrd, smltrd, tmltrd_sum, tmltrd_csum_ln, tmltrd_csum_ub, &
&smltrd_sum, smltrd_csum_ln, smltrd_csum_ub
  CONTAINS
  INTEGER FUNCTION trdmxl_oce_alloc()
    USE lib_mpp
    INTEGER :: ierr(5)
    jpktrd = jpk
    ierr(:) = 0
    ALLOCATE(nmxl(jpi, jpj), nbol(jpi, jpj), wkx(jpi, jpj, jpk), hmxl(jpi, jpj), tml(jpi, jpj), sml(jpi, jpj), tmlb(jpi, jpj), &
&smlb(jpi, jpj), tmlbb(jpi, jpj), smlbb(jpi, jpj), STAT = ierr(1))
    ALLOCATE(tmlbn(jpi, jpj), smlbn(jpi, jpj), tmltrdm(jpi, jpj), smltrdm(jpi, jpj), tml_sum(jpi, jpj), tml_sumb(jpi, jpj), &
&tmltrd_atf_sumb(jpi, jpj), STAT = ierr(2))
    ALLOCATE(sml_sum(jpi, jpj), sml_sumb(jpi, jpj), smltrd_atf_sumb(jpi, jpj), hmxl_sum(jpi, jpj), hmxlbn(jpi, jpj), tmlatfb(jpi, &
&jpj), tmlatfn(jpi, jpj), STAT = ierr(3))
    ALLOCATE(smlatfb(jpi, jpj), smlatfn(jpi, jpj), tmlatfm(jpi, jpj), smlatfm(jpi, jpj), tmltrd(jpi, jpj, jpltrd), smltrd(jpi, &
&jpj, jpltrd), STAT = ierr(4))
    ALLOCATE(tmltrd_sum(jpi, jpj, jpltrd), tmltrd_csum_ln(jpi, jpj, jpltrd), tmltrd_csum_ub(jpi, jpj, jpltrd), smltrd_sum(jpi, &
&jpj, jpltrd), smltrd_csum_ln(jpi, jpj, jpltrd), smltrd_csum_ub(jpi, jpj, jpltrd), STAT = ierr(5))
    trdmxl_oce_alloc = MAXVAL(ierr)
    CALL mpp_sum('trdmxl_oce', trdmxl_oce_alloc)
    IF (trdmxl_oce_alloc /= 0) CALL ctl_stop('STOP', 'trdmxl_oce_alloc: failed to allocate arrays')
  END FUNCTION trdmxl_oce_alloc
END MODULE trdmxl_oce