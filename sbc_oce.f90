MODULE sbc_oce
  USE par_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_oce_alloc
  PUBLIC :: sbc_tau2wnd
  LOGICAL, PUBLIC :: ln_usr
  LOGICAL, PUBLIC :: ln_flx
  LOGICAL, PUBLIC :: ln_blk
  LOGICAL, PUBLIC :: lk_oasis = .FALSE.
  LOGICAL, PUBLIC :: ln_cpl
  LOGICAL, PUBLIC :: ln_mixcpl
  LOGICAL, PUBLIC :: ln_dm2dc
  LOGICAL, PUBLIC :: ln_rnf
  LOGICAL, PUBLIC :: ln_isf
  LOGICAL, PUBLIC :: ln_ssr
  LOGICAL, PUBLIC :: ln_apr_dyn
  INTEGER, PUBLIC :: nn_ice
  LOGICAL, PUBLIC :: ln_ice_embd
  INTEGER, PUBLIC :: nn_components
  INTEGER, PUBLIC :: nn_fwb
  LOGICAL, PUBLIC :: ln_wave
  LOGICAL, PUBLIC :: ln_cdgw
  LOGICAL, PUBLIC :: ln_sdw
  LOGICAL, PUBLIC :: ln_tauwoc
  LOGICAL, PUBLIC :: ln_tauw
  LOGICAL, PUBLIC :: ln_stcor
  INTEGER, PUBLIC :: nn_sdrift
  LOGICAL, PUBLIC :: ln_icebergs
  INTEGER, PUBLIC :: nn_lsm
  INTEGER, PUBLIC :: nn_cats_cpl
  INTEGER, PUBLIC, PARAMETER :: jp_usr = 1
  INTEGER, PUBLIC, PARAMETER :: jp_flx = 2
  INTEGER, PUBLIC, PARAMETER :: jp_blk = 3
  INTEGER, PUBLIC, PARAMETER :: jp_purecpl = 4
  INTEGER, PUBLIC, PARAMETER :: jp_none = 5
  INTEGER, PUBLIC, PARAMETER :: jp_breivik_2014 = 0
  INTEGER, PUBLIC, PARAMETER :: jp_li_2017 = 1
  INTEGER, PUBLIC, PARAMETER :: jp_peakfr = 2
  LOGICAL, PUBLIC :: ll_st_bv2014 = .FALSE.
  LOGICAL, PUBLIC :: ll_st_li2017 = .FALSE.
  LOGICAL, PUBLIC :: ll_st_bv_li = .FALSE.
  LOGICAL, PUBLIC :: ll_st_peakfr = .FALSE.
  INTEGER, PUBLIC, PARAMETER :: jp_iam_nemo = 0
  INTEGER, PUBLIC, PARAMETER :: jp_iam_opa = 1
  INTEGER, PUBLIC, PARAMETER :: jp_iam_sas = 2
  INTEGER, PUBLIC :: ncpl_qsr_freq
  LOGICAL, PUBLIC :: lhftau = .FALSE.
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: utau, utau_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: vtau, vtau_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: taum
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wndm
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qsr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qns, qns_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qsr_tot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qns_tot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: emp, emp_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx, sfx_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: emp_tot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: fmmflx
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: rnf, rnf_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: fwfisf, fwfisf_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: fwficb, fwficb_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sbc_tsc, sbc_tsc_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qsr_hc, qsr_hc_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tprecip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sprecip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: fr_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: atm_co2
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: xcplmask
  INTEGER, PUBLIC :: nn_fsbc
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ssu_m
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ssv_m
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sst_m
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sss_m
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ssh_m
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: e3t_m
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: frq_m
  CONTAINS
  INTEGER FUNCTION sbc_oce_alloc()
    INTEGER :: ierr(5)
    ierr(:) = 0
    ALLOCATE(utau(jpi, jpj), utau_b(jpi, jpj), taum(jpi, jpj), vtau(jpi, jpj), vtau_b(jpi, jpj), wndm(jpi, jpj), STAT = ierr(1))
    ALLOCATE(qns_tot(jpi, jpj), qns(jpi, jpj), qns_b(jpi, jpj), qsr_tot(jpi, jpj), qsr(jpi, jpj), emp(jpi, jpj), emp_b(jpi, jpj), &
&sfx(jpi, jpj), sfx_b(jpi, jpj), emp_tot(jpi, jpj), fmmflx(jpi, jpj), STAT = ierr(2))
    ALLOCATE(fwfisf(jpi, jpj), rnf(jpi, jpj), sbc_tsc(jpi, jpj, jpts), qsr_hc(jpi, jpj, jpk), fwfisf_b(jpi, jpj), rnf_b(jpi, jpj), &
&sbc_tsc_b(jpi, jpj, jpts), qsr_hc_b(jpi, jpj, jpk), fwficb(jpi, jpj), fwficb_b(jpi, jpj), STAT = ierr(3))
    ALLOCATE(tprecip(jpi, jpj), sprecip(jpi, jpj), fr_i(jpi, jpj), atm_co2(jpi, jpj), ssu_m(jpi, jpj), sst_m(jpi, jpj), frq_m(jpi, &
&jpj), ssv_m(jpi, jpj), sss_m(jpi, jpj), ssh_m(jpi, jpj), STAT = ierr(4))
    ALLOCATE(e3t_m(jpi, jpj), STAT = ierr(5))
    sbc_oce_alloc = MAXVAL(ierr)
    CALL mpp_sum('sbc_oce', sbc_oce_alloc)
    IF (sbc_oce_alloc > 0) CALL ctl_warn('sbc_oce_alloc: allocation of arrays failed')
  END FUNCTION sbc_oce_alloc
  SUBROUTINE sbc_tau2wnd
    USE dom_oce
    USE lbclnk
    REAL(KIND = wp) :: zrhoa = 1.22
    REAL(KIND = wp) :: zcdrag = 1.5E-3
    REAL(KIND = wp) :: ztx, zty, ztau, zcoef
    INTEGER :: ji, jj
    zcoef = 0.5 / (zrhoa * zcdrag)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ztx = utau(ji - 1, jj) + utau(ji, jj)
        zty = vtau(ji, jj - 1) + vtau(ji, jj)
        ztau = SQRT(ztx * ztx + zty * zty)
        wndm(ji, jj) = SQRT(ztau * zcoef) * tmask(ji, jj, 1)
      END DO
    END DO
    CALL lbc_lnk('sbc_oce', wndm(:, :), 'T', 1.)
  END SUBROUTINE sbc_tau2wnd
END MODULE sbc_oce