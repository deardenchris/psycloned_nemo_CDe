MODULE icestp
  USE oce
  USE dom_oce
  USE c1d
  USE ice
  USE ice1D
  USE phycst
  USE eosbn2
  USE sbc_oce
  USE sbc_ice
  USE icesbc
  USE icedyn
  USE icethd
  USE icecor
  USE iceupdate
  USE icedia
  USE icewri
  USE icerst
  USE icevar
  USE icectl
  USE iceistate
  USE iceitd
  USE icealb
  USE bdy_oce, ONLY: ln_bdy
  USE bdyice
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE timing
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_stp
  PUBLIC :: ice_init
  CONTAINS
  SUBROUTINE ice_stp(kt, ksbc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: ksbc
    INTEGER :: jl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('ice_stp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ice_stp')
    CALL profile_psy_data0 % PostEnd
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      !$ACC KERNELS
      kt_ice = kt
      u_oce(:, :) = ssu_m(:, :)
      v_oce(:, :) = ssv_m(:, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('ice_stp', 'r1', 0, 0)
      CALL eos_fzp(sss_m(:, :), t_bo(:, :))
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      t_bo(:, :) = (t_bo(:, :) + rt0) * tmask(:, :, 1) + rt0 * (1._wp - tmask(:, :, 1))
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('ice_stp', 'r2', 0, 0)
      CALL store_fields
      CALL ice_sbc_tau(kt, ksbc, utau_ice, vtau_ice)
      CALL diag_set0
      CALL ice_rst_opn(kt)
      IF (ln_icedyn .AND. .NOT. lk_c1d) CALL ice_dyn(kt)
      IF (ln_icethd .AND. ln_bdy) CALL bdy_ice(kt)
      CALL ice_var_glo2eqv
      CALL ice_var_agg(1)
      CALL store_fields
      CALL ice_sbc_flx(kt, ksbc)
      IF (ln_icethd) CALL ice_thd(kt)
      CALL ice_cor(kt, 2)
      CALL ice_var_glo2eqv
      CALL ice_var_agg(2)
      CALL ice_update_flx(kt)
      IF (ln_icediahsb) CALL ice_dia(kt)
      CALL ice_wri(kt)
      IF (lrst_ice) CALL ice_rst_write(kt)
      IF (ln_icectl) CALL ice_ctl(kt)
      CALL profile_psy_data2 % PostEnd
    END IF
    CALL profile_psy_data3 % PreStart('ice_stp', 'r3', 0, 0)
    IF (ln_icedyn) CALL ice_update_tau(kt, ub(:, :, 1), vb(:, :, 1))
    IF (ln_timing) CALL timing_stop('ice_stp')
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE ice_stp
  SUBROUTINE ice_init
    INTEGER :: ji, jj, ierr
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'Sea Ice Model: SI3 (Sea Ice modelling Integrated Initiative)'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~'
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'ice_init: Arrays allocation & Initialization of all routines & init state'
    IF (lwp) WRITE(numout, *) '~~~~~~~~'
    CALL ctl_opn(numnam_ice_ref, 'namelist_ice_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp)
    CALL ctl_opn(numnam_ice_cfg, 'namelist_ice_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp)
    IF (lwm) CALL ctl_opn(numoni, 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, 1)
    CALL par_init
    ierr = ice_alloc()
    ierr = ierr + sbc_ice_alloc()
    ierr = ierr + ice1D_alloc()
    CALL mpp_sum('icestp', ierr)
    IF (ierr /= 0) CALL ctl_stop('STOP', 'ice_init : unable to allocate ice arrays')
    CALL ice_itd_init
    CALL ice_thd_init
    IF (.NOT. ln_rstart) THEN
      CALL ice_istate_init
      CALL ice_istate
    ELSE
      CALL ice_rst_read
    END IF
    CALL ice_var_glo2eqv
    CALL ice_var_agg(1)
    CALL ice_sbc_init
    CALL ice_dyn_init
    CALL ice_update_init
    CALL ice_alb_init
    CALL ice_dia_init
    !$ACC KERNELS
    fr_i(:, :) = at_i(:, :)
    tn_ice(:, :, :) = t_su(:, :, :)
    WHERE (gphit(:, :) > 0._wp)
      rn_amax_2d(:, :) = rn_amax_n
    ELSEWHERE
      rn_amax_2d(:, :) = rn_amax_s
    END WHERE
    !$ACC END KERNELS
    IF (ln_rstart) CALL iom_close(numrir)
  END SUBROUTINE ice_init
  SUBROUTINE par_init
    INTEGER :: ios
    NAMELIST /nampar/ jpl, nlay_i, nlay_s, ln_virtual_itd, ln_icedyn, ln_icethd, rn_amax_n, rn_amax_s, cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, nampar, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nampar in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, nampar, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nampar in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, nampar)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   par_init: ice parameters shared among all the routines'
      WRITE(numout, *) '   ~~~~~~~~'
      WRITE(numout, *) '      Namelist nampar: '
      WRITE(numout, *) '         number of ice  categories                           jpl       = ', jpl
      WRITE(numout, *) '         number of ice  layers                               nlay_i    = ', nlay_i
      WRITE(numout, *) '         number of snow layers                               nlay_s    = ', nlay_s
      WRITE(numout, *) '         virtual ITD param for jpl=1 (T) or not (F)     ln_virtual_itd = ', ln_virtual_itd
      WRITE(numout, *) '         Ice dynamics       (T) or not (F)                   ln_icedyn = ', ln_icedyn
      WRITE(numout, *) '         Ice thermodynamics (T) or not (F)                   ln_icethd = ', ln_icethd
      WRITE(numout, *) '         maximum ice concentration for NH                              = ', rn_amax_n
      WRITE(numout, *) '         maximum ice concentration for SH                              = ', rn_amax_s
    END IF
    rn_amax_n = MIN(rn_amax_n, 1._wp - epsi10)
    rn_amax_s = MIN(rn_amax_s, 1._wp - epsi10)
    IF (jpl > 1 .AND. ln_virtual_itd) THEN
      ln_virtual_itd = .FALSE.
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ln_virtual_itd forced to false as jpl>1, no need with multiple categories to emulate them'
    END IF
    IF (ln_cpl .AND. nn_cats_cpl /= 1 .AND. nn_cats_cpl /= jpl) THEN
      CALL ctl_stop('STOP', 'par_init: in coupled mode, nn_cats_cpl should be either 1 or jpl')
    END IF
    IF (ln_bdy .AND. ln_icediachk) CALL ctl_warn('par_init: online conservation check does not work with BDY')
    rdt_ice = REAL(nn_fsbc) * rdt
    r1_rdtice = 1._wp / rdt_ice
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) '      ice timestep rdt_ice = nn_fsbc*rdt = ', rdt_ice
    r1_nlay_i = 1._wp / REAL(nlay_i, wp)
    r1_nlay_s = 1._wp / REAL(nlay_s, wp)
  END SUBROUTINE par_init
  SUBROUTINE store_fields
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    a_i_b(:, :, :) = a_i(:, :, :)
    v_i_b(:, :, :) = v_i(:, :, :)
    v_s_b(:, :, :) = v_s(:, :, :)
    sv_i_b(:, :, :) = sv_i(:, :, :)
    oa_i_b(:, :, :) = oa_i(:, :, :)
    e_s_b(:, :, :, :) = e_s(:, :, :, :)
    e_i_b(:, :, :, :) = e_i(:, :, :, :)
    WHERE (a_i_b(:, :, :) >= epsi20)
      h_i_b(:, :, :) = v_i_b(:, :, :) / a_i_b(:, :, :)
      h_s_b(:, :, :) = v_s_b(:, :, :) / a_i_b(:, :, :)
    ELSEWHERE
      h_i_b(:, :, :) = 0._wp
      h_s_b(:, :, :) = 0._wp
    END WHERE
    WHERE (a_ip(:, :, :) >= epsi20)
      h_ip_b(:, :, :) = v_ip(:, :, :) / a_ip(:, :, :)
    ELSEWHERE
      h_ip_b(:, :, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('store_fields', 'r0', 0, 0)
    at_i_b(:, :) = SUM(a_i_b(:, :, :), dim = 3)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    u_ice_b(:, :) = u_ice(:, :)
    v_ice_b(:, :) = v_ice(:, :)
    !$ACC END KERNELS
  END SUBROUTINE store_fields
  SUBROUTINE diag_set0
    INTEGER :: ji, jj
    !$ACC KERNELS
    sfx(:, :) = 0._wp
    sfx_bri(:, :) = 0._wp
    sfx_lam(:, :) = 0._wp
    sfx_sni(:, :) = 0._wp
    sfx_opw(:, :) = 0._wp
    sfx_bog(:, :) = 0._wp
    sfx_dyn(:, :) = 0._wp
    sfx_bom(:, :) = 0._wp
    sfx_sum(:, :) = 0._wp
    sfx_res(:, :) = 0._wp
    sfx_sub(:, :) = 0._wp
    wfx_snw(:, :) = 0._wp
    wfx_ice(:, :) = 0._wp
    wfx_sni(:, :) = 0._wp
    wfx_opw(:, :) = 0._wp
    wfx_bog(:, :) = 0._wp
    wfx_dyn(:, :) = 0._wp
    wfx_bom(:, :) = 0._wp
    wfx_sum(:, :) = 0._wp
    wfx_res(:, :) = 0._wp
    wfx_sub(:, :) = 0._wp
    wfx_spr(:, :) = 0._wp
    wfx_lam(:, :) = 0._wp
    wfx_snw_dyn(:, :) = 0._wp
    wfx_snw_sum(:, :) = 0._wp
    wfx_snw_sub(:, :) = 0._wp
    wfx_ice_sub(:, :) = 0._wp
    wfx_snw_sni(:, :) = 0._wp
    wfx_pnd(:, :) = 0._wp
    hfx_thd(:, :) = 0._wp
    hfx_snw(:, :) = 0._wp
    hfx_opw(:, :) = 0._wp
    hfx_bog(:, :) = 0._wp
    hfx_dyn(:, :) = 0._wp
    hfx_bom(:, :) = 0._wp
    hfx_sum(:, :) = 0._wp
    hfx_res(:, :) = 0._wp
    hfx_sub(:, :) = 0._wp
    hfx_spr(:, :) = 0._wp
    hfx_dif(:, :) = 0._wp
    hfx_err_rem(:, :) = 0._wp
    hfx_err_dif(:, :) = 0._wp
    wfx_err_sub(:, :) = 0._wp
    afx_tot(:, :) = 0._wp
    diag_heat(:, :) = 0._wp
    diag_sice(:, :) = 0._wp
    diag_vice(:, :) = 0._wp
    diag_vsnw(:, :) = 0._wp
    qcn_ice_bot(:, :, :) = 0._wp
    qcn_ice_top(:, :, :) = 0._wp
    t_si(:, :, :) = rt0
    tau_icebfr(:, :) = 0._wp
    cnd_ice(:, :, :) = 0._wp
    qcn_ice(:, :, :) = 0._wp
    qtr_ice_bot(:, :, :) = 0._wp
    qsb_ice_bot(:, :) = 0._wp
    diag_trp_vi(:, :) = 0._wp
    diag_trp_vs(:, :) = 0._wp
    diag_trp_ei(:, :) = 0._wp
    diag_trp_es(:, :) = 0._wp
    diag_trp_sv(:, :) = 0._wp
    !$ACC END KERNELS
  END SUBROUTINE diag_set0
END MODULE icestp