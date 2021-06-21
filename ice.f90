MODULE ice
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_alloc
  INTEGER, PUBLIC :: jpl
  INTEGER, PUBLIC :: nlay_i
  INTEGER, PUBLIC :: nlay_s
  LOGICAL, PUBLIC :: ln_virtual_itd
  LOGICAL, PUBLIC :: ln_icedyn
  LOGICAL, PUBLIC :: ln_icethd
  REAL(KIND = wp), PUBLIC :: rn_amax_n
  REAL(KIND = wp), PUBLIC :: rn_amax_s
  CHARACTER(LEN = 256), PUBLIC :: cn_icerst_in
  CHARACTER(LEN = 256), PUBLIC :: cn_icerst_out
  CHARACTER(LEN = 256), PUBLIC :: cn_icerst_indir
  CHARACTER(LEN = 256), PUBLIC :: cn_icerst_outdir
  REAL(KIND = wp), PUBLIC :: rn_himin
  REAL(KIND = wp), PUBLIC :: rn_ishlat
  LOGICAL, PUBLIC :: ln_landfast_L16
  LOGICAL, PUBLIC :: ln_landfast_home
  REAL(KIND = wp), PUBLIC :: rn_depfra
  REAL(KIND = wp), PUBLIC :: rn_icebfr
  REAL(KIND = wp), PUBLIC :: rn_lfrelax
  REAL(KIND = wp), PUBLIC :: rn_tensile
  REAL(KIND = wp), PUBLIC :: rn_crhg
  LOGICAL, PUBLIC :: ln_aEVP
  REAL(KIND = wp), PUBLIC :: rn_creepl
  REAL(KIND = wp), PUBLIC :: rn_ecc
  INTEGER, PUBLIC :: nn_nevp
  REAL(KIND = wp), PUBLIC :: rn_relast
  LOGICAL, PUBLIC :: ln_adv_Pra
  LOGICAL, PUBLIC :: ln_adv_UMx
  REAL(KIND = wp), PUBLIC :: rn_blow_s
  REAL(KIND = wp), PUBLIC :: rn_cio
  INTEGER, PUBLIC :: nn_flxdist
  LOGICAL, PUBLIC :: ln_cndflx
  LOGICAL, PUBLIC :: ln_cndemulate
  INTEGER, PUBLIC, PARAMETER :: np_cnd_OFF = 0
  INTEGER, PUBLIC, PARAMETER :: np_cnd_ON = 1
  INTEGER, PUBLIC, PARAMETER :: np_cnd_EMU = 2
  LOGICAL, PUBLIC :: ln_cndi_U64
  LOGICAL, PUBLIC :: ln_cndi_P07
  REAL(KIND = wp), PUBLIC :: rn_kappa_i
  REAL(KIND = wp), PUBLIC :: rn_cnd_s
  INTEGER, PUBLIC :: nn_icesal
  REAL(KIND = wp), PUBLIC :: rn_icesal
  REAL(KIND = wp), PUBLIC :: rn_simax
  REAL(KIND = wp), PUBLIC :: rn_simin
  LOGICAL, PUBLIC :: ln_pnd_H12
  LOGICAL, PUBLIC :: ln_pnd_CST
  REAL(KIND = wp), PUBLIC :: rn_apnd
  REAL(KIND = wp), PUBLIC :: rn_hpnd
  LOGICAL, PUBLIC :: ln_pnd_alb
  LOGICAL, PUBLIC :: ln_icediachk
  LOGICAL, PUBLIC :: ln_icediahsb
  LOGICAL, PUBLIC :: ln_icectl
  INTEGER, PUBLIC :: iiceprt
  INTEGER, PUBLIC :: jiceprt
  INTEGER, PUBLIC :: kt_ice
  REAL(KIND = wp), PUBLIC :: rdt_ice
  REAL(KIND = wp), PUBLIC :: r1_rdtice
  REAL(KIND = wp), PUBLIC :: r1_nlay_i
  REAL(KIND = wp), PUBLIC :: r1_nlay_s
  REAL(KIND = wp), PUBLIC :: rswitch
  REAL(KIND = wp), PUBLIC :: rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft
  REAL(KIND = wp), PUBLIC, PARAMETER :: epsi06 = 1.E-06_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: epsi10 = 1.E-10_wp
  REAL(KIND = wp), PUBLIC, PARAMETER :: epsi20 = 1.E-20_wp
  LOGICAL, PUBLIC, DIMENSION(2) :: l_split_advumx = .FALSE.
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: u_oce, v_oce
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ht_i_new
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: strength
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: stress1_i, stress2_i, stress12_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: delta_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: divu_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: shear_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: t_bo
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qlead
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qsb_ice_bot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: fhld
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_snw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_snw_sni
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_snw_sum
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_pnd
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_spr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_sub
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_snw_sub
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_ice_sub
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_snw_dyn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_sni
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_opw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_bog
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_dyn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_bom
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_sum
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_lam
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_res
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: wfx_err_sub
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: afx_tot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_bog
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_bom
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_lam
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_sum
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_sni
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_opw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_bri
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_dyn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_res
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sfx_sub
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_bog
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_bom
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_sum
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_opw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_dif
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_snw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_err_dif
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_err_rem
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qt_atm_oi
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qt_oce_ai
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_sub
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_spr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_thd
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_dyn
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hfx_res
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: rn_amax_2d
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qtr_ice_bot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: t1_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: cnd_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: h_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: a_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: v_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: v_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: h_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: t_su
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: s_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sv_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: o_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: oa_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: bv_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: u_ice, v_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: vt_i, vt_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: at_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ato_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: et_i, et_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tm_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tm_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: bvm_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: sm_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tm_su
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hm_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hm_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: om_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tau_icebfr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: t_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: e_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: t_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: e_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: sz_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: a_ip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: v_ip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: a_ip_frac
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: h_ip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: at_ip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: vt_ip
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: v_s_b, v_i_b, h_s_b, h_i_b, h_ip_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: a_i_b, sv_i_b, oa_i_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: e_s_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: e_i_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: u_ice_b, v_ice_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: at_i_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: hi_max
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: hi_mean
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_trp_vi
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_trp_vs
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_trp_ei
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_trp_es
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_trp_sv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_heat
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_sice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_vice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: diag_vsnw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: t_si
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: tm_si
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qcn_ice_bot
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qcn_ice_top
  CONTAINS
  FUNCTION ice_alloc()
    INTEGER :: ice_alloc
    INTEGER :: ierr(15), ii
    ierr(:) = 0
    ii = 1
    ALLOCATE(u_oce(jpi, jpj), v_oce(jpi, jpj), ht_i_new(jpi, jpj), strength(jpi, jpj), stress1_i(jpi, jpj), stress2_i(jpi, jpj), stress12_i(jpi, jpj), delta_i(jpi, jpj), divu_i(jpi, jpj), shear_i(jpi, jpj), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(t_bo(jpi, jpj), wfx_snw_sni(jpi, jpj), wfx_snw(jpi, jpj), wfx_snw_dyn(jpi, jpj), wfx_snw_sum(jpi, jpj), wfx_snw_sub(jpi, jpj), wfx_ice(jpi, jpj), wfx_sub(jpi, jpj), wfx_ice_sub(jpi, jpj), wfx_lam(jpi, jpj), wfx_pnd(jpi, jpj), wfx_bog(jpi, jpj), wfx_dyn(jpi, jpj), wfx_bom(jpi, jpj), wfx_sum(jpi, jpj), wfx_res(jpi, jpj), wfx_sni(jpi, jpj), wfx_opw(jpi, jpj), wfx_spr(jpi, jpj), afx_tot(jpi, jpj), rn_amax_2d(jpi, jpj), qsb_ice_bot(jpi, jpj), qlead(jpi, jpj), sfx_res(jpi, jpj), sfx_bri(jpi, jpj), sfx_dyn(jpi, jpj), sfx_sub(jpi, jpj), sfx_lam(jpi, jpj), sfx_bog(jpi, jpj), sfx_bom(jpi, jpj), sfx_sum(jpi, jpj), sfx_sni(jpi, jpj), sfx_opw(jpi, jpj), hfx_res(jpi, jpj), hfx_snw(jpi, jpj), hfx_sub(jpi, jpj), qt_atm_oi(jpi, jpj), qt_oce_ai(jpi, jpj), fhld(jpi, jpj), hfx_sum(jpi, jpj), hfx_bom(jpi, jpj), hfx_bog(jpi, jpj), hfx_dif(jpi, jpj), hfx_opw(jpi, jpj), hfx_thd(jpi, jpj), hfx_dyn(jpi, jpj), hfx_spr(jpi, jpj), hfx_err_dif(jpi, jpj), hfx_err_rem(jpi, jpj), wfx_err_sub(jpi, jpj), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(qtr_ice_bot(jpi, jpj, jpl), cnd_ice(jpi, jpj, jpl), t1_ice(jpi, jpj, jpl), h_i(jpi, jpj, jpl), a_i(jpi, jpj, jpl), v_i(jpi, jpj, jpl), v_s(jpi, jpj, jpl), h_s(jpi, jpj, jpl), t_su(jpi, jpj, jpl), s_i(jpi, jpj, jpl), sv_i(jpi, jpj, jpl), o_i(jpi, jpj, jpl), oa_i(jpi, jpj, jpl), bv_i(jpi, jpj, jpl), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(u_ice(jpi, jpj), v_ice(jpi, jpj), vt_i(jpi, jpj), vt_s(jpi, jpj), at_i(jpi, jpj), ato_i(jpi, jpj), et_i(jpi, jpj), et_s(jpi, jpj), tm_i(jpi, jpj), tm_s(jpi, jpj), sm_i(jpi, jpj), tm_su(jpi, jpj), hm_i(jpi, jpj), hm_s(jpi, jpj), om_i(jpi, jpj), bvm_i(jpi, jpj), tau_icebfr(jpi, jpj), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(t_s(jpi, jpj, nlay_s, jpl), e_s(jpi, jpj, nlay_s, jpl), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(t_i(jpi, jpj, nlay_i, jpl), e_i(jpi, jpj, nlay_i, jpl), sz_i(jpi, jpj, nlay_i, jpl), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(a_ip(jpi, jpj, jpl), v_ip(jpi, jpj, jpl), a_ip_frac(jpi, jpj, jpl), h_ip(jpi, jpj, jpl), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(at_ip(jpi, jpj), vt_ip(jpi, jpj), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(v_s_b(jpi, jpj, jpl), v_i_b(jpi, jpj, jpl), h_s_b(jpi, jpj, jpl), h_i_b(jpi, jpj, jpl), h_ip_b(jpi, jpj, jpl), a_i_b(jpi, jpj, jpl), sv_i_b(jpi, jpj, jpl), e_i_b(jpi, jpj, nlay_i, jpl), e_s_b(jpi, jpj, nlay_s, jpl), oa_i_b(jpi, jpj, jpl), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(u_ice_b(jpi, jpj), v_ice_b(jpi, jpj), at_i_b(jpi, jpj), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(hi_max(0 : jpl), hi_mean(jpl), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(diag_trp_vi(jpi, jpj), diag_trp_vs(jpi, jpj), diag_trp_ei(jpi, jpj), diag_trp_es(jpi, jpj), diag_trp_sv(jpi, jpj), diag_heat(jpi, jpj), diag_sice(jpi, jpj), diag_vice(jpi, jpj), diag_vsnw(jpi, jpj), STAT = ierr(ii))
    ii = ii + 1
    ALLOCATE(t_si(jpi, jpj, jpl), tm_si(jpi, jpj), qcn_ice_bot(jpi, jpj, jpl), qcn_ice_top(jpi, jpj, jpl), STAT = ierr(ii))
    ice_alloc = MAXVAL(ierr(:))
    IF (ice_alloc /= 0) CALL ctl_stop('STOP', 'ice_alloc: failed to allocate arrays.')
  END FUNCTION ice_alloc
END MODULE ice