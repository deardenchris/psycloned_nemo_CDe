MODULE icethd_dh
  USE dom_oce
  USE phycst
  USE ice
  USE ice1D
  USE icethd_sal
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_dh
  PUBLIC :: ice_thd_snwblow
  INTERFACE ice_thd_snwblow
    MODULE PROCEDURE ice_thd_snwblow_1d, ice_thd_snwblow_2d
  END INTERFACE
  CONTAINS
  SUBROUTINE ice_thd_dh
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jk
    INTEGER :: iter
    REAL(KIND = wp) :: ztmelts
    REAL(KIND = wp) :: zdum
    REAL(KIND = wp) :: zfracs
    REAL(KIND = wp) :: zswi1
    REAL(KIND = wp) :: zswi12
    REAL(KIND = wp) :: zswi2
    REAL(KIND = wp) :: zgrr
    REAL(KIND = wp) :: zt_i_new
    REAL(KIND = wp) :: z1_rho
    REAL(KIND = wp) :: zQm
    REAL(KIND = wp) :: zEi
    REAL(KIND = wp) :: zEw
    REAL(KIND = wp) :: zdE
    REAL(KIND = wp) :: zfmdt
    REAL(KIND = wp), DIMENSION(jpij) :: zqprec
    REAL(KIND = wp), DIMENSION(jpij) :: zq_top
    REAL(KIND = wp), DIMENSION(jpij) :: zq_bot
    REAL(KIND = wp), DIMENSION(jpij) :: zq_rema
    REAL(KIND = wp), DIMENSION(jpij) :: zf_tt
    REAL(KIND = wp), DIMENSION(jpij) :: zevap_rema
    REAL(KIND = wp), DIMENSION(jpij) :: zdh_s_mel
    REAL(KIND = wp), DIMENSION(jpij) :: zdh_s_pre
    REAL(KIND = wp), DIMENSION(jpij) :: zdh_s_sub
    REAL(KIND = wp), DIMENSION(jpij, nlay_s) :: zh_s
    REAL(KIND = wp), DIMENSION(jpij, nlay_i) :: zh_i
    REAL(KIND = wp), DIMENSION(jpij, nlay_i) :: zdeltah
    INTEGER, DIMENSION(jpij, nlay_i) :: icount
    REAL(KIND = wp), DIMENSION(jpij) :: zsnw
    REAL(KIND = wp) :: zswitch_sal
    INTEGER :: num_iter_max
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_thd_dh', 'r0', 0, 0)
    SELECT CASE (nn_icesal)
    CASE (1, 3)
      zswitch_sal = 0._wp
    CASE (2)
      zswitch_sal = 1._wp
    END SELECT
    h_i_old(1 : npti, 0 : nlay_i + 1) = 0._wp
    eh_i_old(1 : npti, 0 : nlay_i + 1) = 0._wp
    DO jk = 1, nlay_i
      DO ji = 1, npti
        h_i_old(ji, jk) = h_i_1d(ji) * r1_nlay_i
        eh_i_old(ji, jk) = e_i_1d(ji, jk) * h_i_old(ji, jk)
      END DO
    END DO
    IF (ln_cndflx .AND. .NOT. ln_cndemulate) THEN
      DO ji = 1, npti
        zq_top(ji) = MAX(0._wp, qml_ice_1d(ji) * rdt_ice)
      END DO
    ELSE
      DO ji = 1, npti
        zdum = qns_ice_1d(ji) + qsr_ice_1d(ji) - qtr_ice_top_1d(ji) - qcn_ice_top_1d(ji)
        qml_ice_1d(ji) = zdum * MAX(0._wp, SIGN(1._wp, t_su_1d(ji) - rt0))
        zq_top(ji) = MAX(0._wp, qml_ice_1d(ji) * rdt_ice)
      END DO
    END IF
    DO ji = 1, npti
      zf_tt(ji) = qcn_ice_bot_1d(ji) + qsb_ice_bot_1d(ji) + fhld_1d(ji)
      zq_bot(ji) = MAX(0._wp, zf_tt(ji) * rdt_ice)
    END DO
    DO jk = 1, nlay_i
      DO ji = 1, npti
        zh_i(ji, jk) = h_i_1d(ji) * r1_nlay_i
      END DO
    END DO
    DO jk = 1, nlay_s
      DO ji = 1, npti
        zh_s(ji, jk) = h_s_1d(ji) * r1_nlay_s
      END DO
    END DO
    DO jk = 1, nlay_s
      DO ji = 1, npti
        IF (t_s_1d(ji, jk) > rt0) THEN
          hfx_res_1d(ji) = hfx_res_1d(ji) + e_s_1d(ji, jk) * zh_s(ji, jk) * a_i_1d(ji) * r1_rdtice
          wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) + rhos * zh_s(ji, jk) * a_i_1d(ji) * r1_rdtice
          dh_s_mlt(ji) = dh_s_mlt(ji) - zh_s(ji, jk)
          h_s_1d(ji) = h_s_1d(ji) - zh_s(ji, jk)
          zh_s(ji, jk) = 0._wp
          e_s_1d(ji, jk) = 0._wp
          t_s_1d(ji, jk) = rt0
        END IF
      END DO
    END DO
    CALL ice_thd_snwblow(1. - at_i_1d(1 : npti), zsnw(1 : npti))
    zdeltah(1 : npti, :) = 0._wp
    DO ji = 1, npti
      IF (sprecip_1d(ji) > 0._wp) THEN
        zdh_s_pre(ji) = zsnw(ji) * sprecip_1d(ji) * rdt_ice * r1_rhos / at_i_1d(ji)
        zqprec(ji) = - qprec_ice_1d(ji)
        hfx_spr_1d(ji) = hfx_spr_1d(ji) + zdh_s_pre(ji) * a_i_1d(ji) * zqprec(ji) * r1_rdtice
        wfx_spr_1d(ji) = wfx_spr_1d(ji) - rhos * a_i_1d(ji) * zdh_s_pre(ji) * r1_rdtice
        rswitch = MAX(0._wp, SIGN(1._wp, zqprec(ji) - epsi20))
        zdeltah(ji, 1) = - rswitch * zq_top(ji) / MAX(zqprec(ji), epsi20)
        zdeltah(ji, 1) = MAX(- zdh_s_pre(ji), zdeltah(ji, 1))
        hfx_snw_1d(ji) = hfx_snw_1d(ji) - zdeltah(ji, 1) * a_i_1d(ji) * zqprec(ji) * r1_rdtice
        wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) - rhos * a_i_1d(ji) * zdeltah(ji, 1) * r1_rdtice
        dh_s_mlt(ji) = dh_s_mlt(ji) + zdeltah(ji, 1)
        zq_top(ji) = MAX(0._wp, zq_top(ji) + zdeltah(ji, 1) * zqprec(ji))
        zdh_s_pre(ji) = zdh_s_pre(ji) + zdeltah(ji, 1)
        h_s_1d(ji) = MAX(0._wp, h_s_1d(ji) + zdh_s_pre(ji))
      ELSE
        zdh_s_pre(ji) = 0._wp
        zqprec(ji) = 0._wp
      END IF
    END DO
    DO jk = 1, nlay_s
      DO ji = 1, npti
        zh_s(ji, jk) = h_s_1d(ji) * r1_nlay_s
      END DO
    END DO
    zdeltah(1 : npti, :) = 0._wp
    zdh_s_mel(1 : npti) = 0._wp
    DO jk = 1, nlay_s
      DO ji = 1, npti
        IF (zh_s(ji, jk) > 0._wp .AND. zq_top(ji) > 0._wp) THEN
          rswitch = MAX(0._wp, SIGN(1._wp, e_s_1d(ji, jk) - epsi20))
          zdeltah(ji, jk) = - rswitch * zq_top(ji) / MAX(e_s_1d(ji, jk), epsi20)
          zdeltah(ji, jk) = MAX(zdeltah(ji, jk), - zh_s(ji, jk))
          zdh_s_mel(ji) = zdh_s_mel(ji) + zdeltah(ji, jk)
          hfx_snw_1d(ji) = hfx_snw_1d(ji) - zdeltah(ji, jk) * a_i_1d(ji) * e_s_1d(ji, jk) * r1_rdtice
          wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) - rhos * a_i_1d(ji) * zdeltah(ji, jk) * r1_rdtice
          dh_s_mlt(ji) = dh_s_mlt(ji) + zdeltah(ji, jk)
          zq_top(ji) = MAX(0._wp, zq_top(ji) + zdeltah(ji, jk) * e_s_1d(ji, jk))
          h_s_1d(ji) = MAX(0._wp, h_s_1d(ji) + zdeltah(ji, jk))
          zh_s(ji, jk) = MAX(0._wp, zh_s(ji, jk) + zdeltah(ji, jk))
        END IF
      END DO
    END DO
    zdeltah(1 : npti, :) = 0._wp
    DO ji = 1, npti
      IF (evap_ice_1d(ji) > 0._wp) THEN
        zdh_s_sub(ji) = MAX(- h_s_1d(ji), - evap_ice_1d(ji) * r1_rhos * rdt_ice)
        zevap_rema(ji) = evap_ice_1d(ji) * rdt_ice + zdh_s_sub(ji) * rhos
        zdeltah(ji, 1) = MAX(zdh_s_sub(ji), - zdh_s_pre(ji))
        hfx_sub_1d(ji) = hfx_sub_1d(ji) + (zdeltah(ji, 1) * zqprec(ji) + (zdh_s_sub(ji) - zdeltah(ji, 1)) * e_s_1d(ji, 1)) * a_i_1d(ji) * r1_rdtice
        wfx_snw_sub_1d(ji) = wfx_snw_sub_1d(ji) - rhos * a_i_1d(ji) * zdh_s_sub(ji) * r1_rdtice
        h_s_1d(ji) = MAX(0._wp, h_s_1d(ji) + zdh_s_sub(ji))
        zdh_s_pre(ji) = zdh_s_pre(ji) + zdeltah(ji, 1)
        zdh_s_sub(ji) = zdh_s_sub(ji) - zdeltah(ji, 1)
      ELSE
        zdh_s_sub(ji) = 0._wp
        zevap_rema(ji) = 0._wp
      END IF
    END DO
    DO ji = 1, npti
      dh_s_tot(ji) = zdh_s_mel(ji) + zdh_s_pre(ji) + zdh_s_sub(ji)
    END DO
    DO jk = 1, nlay_s
      DO ji = 1, npti
        rswitch = MAX(0._wp, SIGN(1._wp, h_s_1d(ji) - epsi20))
        e_s_1d(ji, jk) = rswitch / MAX(h_s_1d(ji), epsi20) * ((zdh_s_pre(ji)) * zqprec(ji) + (h_s_1d(ji) - zdh_s_pre(ji)) * rhos * (rcpi * (rt0 - t_s_1d(ji, jk)) + rLfus))
      END DO
    END DO
    zdeltah(1 : npti, :) = 0._wp
    DO jk = 1, nlay_i
      DO ji = 1, npti
        ztmelts = - rTmlt * sz_i_1d(ji, jk)
        IF (t_i_1d(ji, jk) >= (ztmelts + rt0)) THEN
          zEi = - e_i_1d(ji, jk) * r1_rhoi
          zdE = 0._wp
          zdeltah(ji, jk) = MIN(0._wp, - zh_i(ji, jk))
          zfmdt = - zdeltah(ji, jk) * rhoi
          dh_i_itm(ji) = dh_i_itm(ji) + zdeltah(ji, jk)
          zfmdt = - rhoi * zdeltah(ji, jk)
          hfx_res_1d(ji) = hfx_res_1d(ji) + zfmdt * a_i_1d(ji) * zEi * r1_rdtice
          sfx_res_1d(ji) = sfx_res_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * s_i_1d(ji) * r1_rdtice
          wfx_res_1d(ji) = wfx_res_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * r1_rdtice
        ELSE
          zEi = - e_i_1d(ji, jk) * r1_rhoi
          zEw = rcp * ztmelts
          zdE = zEi - zEw
          zfmdt = - zq_top(ji) / zdE
          zdeltah(ji, jk) = - zfmdt * r1_rhoi
          zdeltah(ji, jk) = MIN(0._wp, MAX(zdeltah(ji, jk), - zh_i(ji, jk)))
          zq_top(ji) = MAX(0._wp, zq_top(ji) - zdeltah(ji, jk) * rhoi * zdE)
          dh_i_sum(ji) = dh_i_sum(ji) + zdeltah(ji, jk)
          zfmdt = - rhoi * zdeltah(ji, jk)
          zQm = zfmdt * zEw
          sfx_sum_1d(ji) = sfx_sum_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * s_i_1d(ji) * r1_rdtice
          hfx_thd_1d(ji) = hfx_thd_1d(ji) + zfmdt * a_i_1d(ji) * zEw * r1_rdtice
          hfx_sum_1d(ji) = hfx_sum_1d(ji) - zfmdt * a_i_1d(ji) * zdE * r1_rdtice
          wfx_sum_1d(ji) = wfx_sum_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * r1_rdtice
        END IF
        zdum = MAX(- (zh_i(ji, jk) + zdeltah(ji, jk)), - zevap_rema(ji) * r1_rhoi)
        zdeltah(ji, jk) = zdeltah(ji, jk) + zdum
        dh_i_sub(ji) = dh_i_sub(ji) + zdum
        sfx_sub_1d(ji) = sfx_sub_1d(ji) - rhoi * a_i_1d(ji) * zdum * s_i_1d(ji) * r1_rdtice
        hfx_sub_1d(ji) = hfx_sub_1d(ji) + zdum * e_i_1d(ji, jk) * a_i_1d(ji) * r1_rdtice
        wfx_ice_sub_1d(ji) = wfx_ice_sub_1d(ji) - rhoi * a_i_1d(ji) * zdum * r1_rdtice
        zevap_rema(ji) = zevap_rema(ji) + zdum * rhoi
        rswitch = MAX(0._wp, SIGN(1._wp, - (zh_i(ji, jk) + zdeltah(ji, jk))))
        icount(ji, jk) = NINT(rswitch)
        zh_i(ji, jk) = MAX(0._wp, zh_i(ji, jk) + zdeltah(ji, jk))
        eh_i_old(ji, jk) = eh_i_old(ji, jk) + zdeltah(ji, jk) * e_i_1d(ji, jk)
        h_i_old(ji, jk) = h_i_old(ji, jk) + zdeltah(ji, jk)
      END DO
    END DO
    DO ji = 1, npti
      h_i_1d(ji) = MAX(0._wp, h_i_1d(ji) + dh_i_sum(ji) + dh_i_itm(ji) + dh_i_sub(ji))
    END DO
    DO ji = 1, npti
      wfx_err_sub_1d(ji) = wfx_err_sub_1d(ji) - zevap_rema(ji) * a_i_1d(ji) * r1_rdtice
    END DO
    num_iter_max = 1
    IF (nn_icesal == 2) num_iter_max = 5
    DO ji = 1, npti
      IF (zf_tt(ji) < 0._wp) THEN
        DO iter = 1, num_iter_max
          zgrr = MIN(1.0E-3, MAX(dh_i_bog(ji) * r1_rdtice, epsi10))
          zswi2 = MAX(0._wp, SIGN(1._wp, zgrr - 3.6E-7))
          zswi12 = MAX(0._wp, SIGN(1._wp, zgrr - 2.0E-8)) * (1.0 - zswi2)
          zswi1 = 1. - zswi2 * zswi12
          zfracs = MIN(zswi1 * 0.12 + zswi12 * (0.8925 + 0.0568 * LOG(100.0 * zgrr)) + zswi2 * 0.26 / (0.26 + 0.74 * EXP(- 724300.0 * zgrr)), 0.5)
          s_i_new(ji) = zswitch_sal * zfracs * sss_1d(ji) + (1. - zswitch_sal) * s_i_1d(ji)
          ztmelts = - rTmlt * s_i_new(ji)
          zt_i_new = zswitch_sal * t_bo_1d(ji) + (1. - zswitch_sal) * t_i_1d(ji, nlay_i)
          zEi = rcpi * (zt_i_new - (ztmelts + rt0)) - rLfus * (1.0 - ztmelts / (zt_i_new - rt0)) + rcp * ztmelts
          zEw = rcp * (t_bo_1d(ji) - rt0)
          zdE = zEi - zEw
          dh_i_bog(ji) = rdt_ice * MAX(0._wp, zf_tt(ji) / (zdE * rhoi))
        END DO
        zfmdt = - rhoi * dh_i_bog(ji)
        hfx_thd_1d(ji) = hfx_thd_1d(ji) + zfmdt * a_i_1d(ji) * zEw * r1_rdtice
        hfx_bog_1d(ji) = hfx_bog_1d(ji) - zfmdt * a_i_1d(ji) * zdE * r1_rdtice
        sfx_bog_1d(ji) = sfx_bog_1d(ji) - rhoi * a_i_1d(ji) * dh_i_bog(ji) * s_i_new(ji) * r1_rdtice
        wfx_bog_1d(ji) = wfx_bog_1d(ji) - rhoi * a_i_1d(ji) * dh_i_bog(ji) * r1_rdtice
        eh_i_old(ji, nlay_i + 1) = eh_i_old(ji, nlay_i + 1) + dh_i_bog(ji) * (- zEi * rhoi)
        h_i_old(ji, nlay_i + 1) = h_i_old(ji, nlay_i + 1) + dh_i_bog(ji)
      END IF
    END DO
    zdeltah(1 : npti, :) = 0._wp
    DO jk = nlay_i, 1, - 1
      DO ji = 1, npti
        IF (zf_tt(ji) > 0._wp .AND. jk > icount(ji, jk)) THEN
          ztmelts = - rTmlt * sz_i_1d(ji, jk)
          IF (t_i_1d(ji, jk) >= (ztmelts + rt0)) THEN
            zEi = - e_i_1d(ji, jk) * r1_rhoi
            zdE = 0._wp
            zdeltah(ji, jk) = MIN(0._wp, - zh_i(ji, jk))
            dh_i_itm(ji) = dh_i_itm(ji) + zdeltah(ji, jk)
            zfmdt = - zdeltah(ji, jk) * rhoi
            hfx_res_1d(ji) = hfx_res_1d(ji) + zfmdt * a_i_1d(ji) * zEi * r1_rdtice
            sfx_res_1d(ji) = sfx_res_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * s_i_1d(ji) * r1_rdtice
            wfx_res_1d(ji) = wfx_res_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * r1_rdtice
            eh_i_old(ji, jk) = eh_i_old(ji, jk) + zdeltah(ji, jk) * e_i_1d(ji, jk)
            h_i_old(ji, jk) = h_i_old(ji, jk) + zdeltah(ji, jk)
          ELSE
            zEi = - e_i_1d(ji, jk) * r1_rhoi
            zEw = rcp * ztmelts
            zdE = zEi - zEw
            zfmdt = - zq_bot(ji) / zdE
            zdeltah(ji, jk) = - zfmdt * r1_rhoi
            zdeltah(ji, jk) = MIN(0._wp, MAX(zdeltah(ji, jk), - zh_i(ji, jk)))
            zq_bot(ji) = MAX(0._wp, zq_bot(ji) - zdeltah(ji, jk) * rhoi * zdE)
            dh_i_bom(ji) = dh_i_bom(ji) + zdeltah(ji, jk)
            zfmdt = - zdeltah(ji, jk) * rhoi
            zQm = zfmdt * zEw
            hfx_thd_1d(ji) = hfx_thd_1d(ji) + zfmdt * a_i_1d(ji) * zEw * r1_rdtice
            hfx_bom_1d(ji) = hfx_bom_1d(ji) - zfmdt * a_i_1d(ji) * zdE * r1_rdtice
            sfx_bom_1d(ji) = sfx_bom_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * s_i_1d(ji) * r1_rdtice
            wfx_bom_1d(ji) = wfx_bom_1d(ji) - rhoi * a_i_1d(ji) * zdeltah(ji, jk) * r1_rdtice
            eh_i_old(ji, jk) = eh_i_old(ji, jk) + zdeltah(ji, jk) * e_i_1d(ji, jk)
            h_i_old(ji, jk) = h_i_old(ji, jk) + zdeltah(ji, jk)
          END IF
        END IF
      END DO
    END DO
    DO ji = 1, npti
      h_i_1d(ji) = MAX(0._wp, h_i_1d(ji) + dh_i_bog(ji) + dh_i_bom(ji))
    END DO
    zdeltah(1 : npti, :) = 0._wp
    DO ji = 1, npti
      zq_rema(ji) = zq_top(ji) + zq_bot(ji)
      rswitch = 1._wp - MAX(0._wp, SIGN(1._wp, - h_s_1d(ji)))
      rswitch = rswitch * MAX(0._wp, SIGN(1._wp, e_s_1d(ji, 1) - epsi20))
      zdeltah(ji, 1) = - rswitch * zq_rema(ji) / MAX(e_s_1d(ji, 1), epsi20)
      zdeltah(ji, 1) = MIN(0._wp, MAX(zdeltah(ji, 1), - h_s_1d(ji)))
      dh_s_tot(ji) = dh_s_tot(ji) + zdeltah(ji, 1)
      h_s_1d(ji) = h_s_1d(ji) + zdeltah(ji, 1)
      zq_rema(ji) = zq_rema(ji) + zdeltah(ji, 1) * e_s_1d(ji, 1)
      hfx_snw_1d(ji) = hfx_snw_1d(ji) - zdeltah(ji, 1) * a_i_1d(ji) * e_s_1d(ji, 1) * r1_rdtice
      wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) - rhos * a_i_1d(ji) * zdeltah(ji, 1) * r1_rdtice
      dh_s_mlt(ji) = dh_s_mlt(ji) + zdeltah(ji, 1)
      qt_oce_ai_1d(ji) = qt_oce_ai_1d(ji) + (zq_rema(ji) * a_i_1d(ji)) * r1_rdtice
      IF (ln_icectl .AND. zq_rema(ji) < 0. .AND. lwp) WRITE(numout, *) 'ALERTE zq_rema <0 = ', zq_rema(ji)
    END DO
    z1_rho = 1._wp / (rhos + rau0 - rhoi)
    DO ji = 1, npti
      dh_snowice(ji) = MAX(0._wp, (rhos * h_s_1d(ji) + (rhoi - rau0) * h_i_1d(ji)) * z1_rho)
      h_i_1d(ji) = h_i_1d(ji) + dh_snowice(ji)
      h_s_1d(ji) = h_s_1d(ji) - dh_snowice(ji)
      zfmdt = (rhos - rhoi) * dh_snowice(ji)
      zEw = rcp * sst_1d(ji)
      zQm = zfmdt * zEw
      hfx_thd_1d(ji) = hfx_thd_1d(ji) + zfmdt * a_i_1d(ji) * zEw * r1_rdtice
      sfx_sni_1d(ji) = sfx_sni_1d(ji) + sss_1d(ji) * a_i_1d(ji) * zfmdt * r1_rdtice
      IF (nn_icesal /= 2) THEN
        sfx_bri_1d(ji) = sfx_bri_1d(ji) - sss_1d(ji) * a_i_1d(ji) * zfmdt * r1_rdtice - s_i_1d(ji) * a_i_1d(ji) * dh_snowice(ji) * rhoi * r1_rdtice
      END IF
      wfx_sni_1d(ji) = wfx_sni_1d(ji) - a_i_1d(ji) * dh_snowice(ji) * rhoi * r1_rdtice
      wfx_snw_sni_1d(ji) = wfx_snw_sni_1d(ji) + a_i_1d(ji) * dh_snowice(ji) * rhos * r1_rdtice
      eh_i_old(ji, 0) = eh_i_old(ji, 0) + dh_snowice(ji) * e_s_1d(ji, 1) + zfmdt * zEw
      h_i_old(ji, 0) = h_i_old(ji, 0) + dh_snowice(ji)
    END DO
    DO ji = 1, npti
      rswitch = 1._wp - MAX(0._wp, SIGN(1._wp, - h_i_1d(ji)))
      t_su_1d(ji) = rswitch * t_su_1d(ji) + (1._wp - rswitch) * rt0
    END DO
    DO jk = 1, nlay_s
      DO ji = 1, npti
        rswitch = (1._wp - MAX(0._wp, SIGN(1._wp, - h_s_1d(ji)))) * (1._wp - MAX(0._wp, SIGN(1._wp, - h_i_1d(ji))))
        hfx_res_1d(ji) = hfx_res_1d(ji) + (1._wp - rswitch) * (e_s_1d(ji, jk) * h_s_1d(ji) * r1_nlay_s * a_i_1d(ji) * r1_rdtice)
        wfx_res_1d(ji) = wfx_res_1d(ji) + (1._wp - rswitch) * (rhos * h_s_1d(ji) * r1_nlay_s * a_i_1d(ji) * r1_rdtice)
        e_s_1d(ji, jk) = rswitch * e_s_1d(ji, jk)
        t_s_1d(ji, jk) = rt0 + rswitch * (- e_s_1d(ji, jk) * r1_rhos * r1_rcpi + rLfus * r1_rcpi)
      END DO
    END DO
    WHERE (h_i_1d(1 : npti) == 0._wp)
      a_i_1d(1 : npti) = 0._wp
      h_s_1d(1 : npti) = 0._wp
    END WHERE
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_thd_dh
  SUBROUTINE ice_thd_snwblow_2d(pin, pout)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pin
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pout
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_thd_snwblow_2d', 'r0', 0, 0)
    pout = (1._wp - (pin) ** rn_blow_s)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_thd_snwblow_2d
  SUBROUTINE ice_thd_snwblow_1d(pin, pout)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pin
    REAL(KIND = wp), DIMENSION(:), INTENT(INOUT) :: pout
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_thd_snwblow_1d', 'r0', 0, 0)
    pout = (1._wp - (pin) ** rn_blow_s)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_thd_snwblow_1d
END MODULE icethd_dh