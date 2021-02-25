MODULE icedyn_adv
  USE phycst
  USE dom_oce
  USE sbc_oce, ONLY: nn_fsbc
  USE ice
  USE icevar
  USE icedyn_adv_pra
  USE icedyn_adv_umx
  USE icectl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE timing
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_adv
  PUBLIC :: ice_dyn_adv_init
  INTEGER :: nice_adv
  INTEGER, PARAMETER :: np_advPRA = 1
  INTEGER, PARAMETER :: np_advUMx = 2
  INTEGER :: nn_UMx
  CONTAINS
  SUBROUTINE ice_dyn_adv(kt)
    INTEGER, INTENT(IN) :: kt
    IF (ln_timing) CALL timing_start('icedyn_adv')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn_adv: sea-ice advection'
      WRITE(numout, FMT = *) '~~~~~~~~~~~'
    END IF
    SELECT CASE (nice_adv)
    CASE (np_advUMx)
      CALL ice_dyn_adv_umx(nn_UMx, kt, u_ice, v_ice, h_i, h_s, h_ip, ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, e_s, e_i)
    CASE (np_advPRA)
      CALL ice_dyn_adv_pra(kt, u_ice, v_ice, ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, e_s, e_i)
    END SELECT
    diag_trp_ei(:, :) = SUM(SUM(e_i(:, :, 1 : nlay_i, :) - e_i_b(:, :, 1 : nlay_i, :), dim = 4), dim = 3) * r1_rdtice
    diag_trp_es(:, :) = SUM(SUM(e_s(:, :, 1 : nlay_s, :) - e_s_b(:, :, 1 : nlay_s, :), dim = 4), dim = 3) * r1_rdtice
    diag_trp_sv(:, :) = SUM(sv_i(:, :, :) - sv_i_b(:, :, :), dim = 3) * r1_rdtice
    diag_trp_vi(:, :) = SUM(v_i(:, :, :) - v_i_b(:, :, :), dim = 3) * r1_rdtice
    diag_trp_vs(:, :) = SUM(v_s(:, :, :) - v_s_b(:, :, :), dim = 3) * r1_rdtice
    IF (iom_use('icemtrp')) CALL iom_put("icemtrp", diag_trp_vi * rhoi)
    IF (iom_use('snwmtrp')) CALL iom_put("snwmtrp", diag_trp_vs * rhos)
    IF (iom_use('salmtrp')) CALL iom_put("salmtrp", diag_trp_sv * rhoi * 1.E-03)
    IF (iom_use('dihctrp')) CALL iom_put("dihctrp", - diag_trp_ei)
    IF (iom_use('dshctrp')) CALL iom_put("dshctrp", - diag_trp_es)
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (ln_icectl) CALL ice_prt(kt, iiceprt, jiceprt, - 1, ' - ice dyn & trp - ')
    IF (ln_timing) CALL timing_stop('icedyn_adv')
  END SUBROUTINE ice_dyn_adv
  SUBROUTINE ice_dyn_adv_init
    INTEGER :: ios, ioptio
    NAMELIST /namdyn_adv/ ln_adv_Pra, ln_adv_UMx, nn_UMx
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namdyn_adv, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_adv in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namdyn_adv, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_adv in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namdyn_adv)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn_adv_init: ice parameters for ice dynamics '
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namdyn_adv:'
      WRITE(numout, FMT = *) '      type of advection scheme (Prather)             ln_adv_Pra = ', ln_adv_Pra
      WRITE(numout, FMT = *) '      type of advection scheme (Ulimate-Macho)       ln_adv_UMx = ', ln_adv_UMx
      WRITE(numout, FMT = *) '         order of the Ultimate-Macho scheme          nn_UMx     = ', nn_UMx
    END IF
    ioptio = 0
    IF (ln_adv_pra) THEN
      ioptio = ioptio + 1
      nice_adv = np_advpra
    END IF
    IF (ln_adv_umx) THEN
      ioptio = ioptio + 1
      nice_adv = np_advumx
    END IF
    IF (ioptio /= 1) CALL ctl_stop('ice_dyn_adv_init: choose one and only one ice adv. scheme (ln_adv_Pra or ln_adv_UMx)')
    IF (ln_adv_Pra) CALL adv_pra_init
  END SUBROUTINE ice_dyn_adv_init
END MODULE icedyn_adv