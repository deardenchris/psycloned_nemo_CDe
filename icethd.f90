MODULE icethd
  USE phycst
  USE dom_oce
  USE ice
  USE sbc_oce, ONLY: sss_m, sst_m, e3t_m, utau, vtau, ssu_m, ssv_m, frq_m, qns_tot, qsr_tot, sprecip, ln_cpl
  USE sbc_ice, ONLY: qsr_oce, qns_oce, qemp_oce, qsr_ice, qns_ice, dqns_ice, evap_ice, qprec_ice, qevap_ice, qml_ice, qcn_ice, &
&qtr_ice_top
  USE ice1D
  USE icethd_zdf
  USE icethd_dh
  USE icethd_da
  USE icethd_sal
  USE icethd_ent
  USE icethd_do
  USE icethd_pnd
  USE iceitd
  USE icetab
  USE icevar
  USE icectl
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd
  PUBLIC :: ice_thd_init
  LOGICAL :: ln_icedH
  LOGICAL :: ln_icedA
  LOGICAL :: ln_icedO
  LOGICAL :: ln_icedS
  CONTAINS
  SUBROUTINE ice_thd(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: zfric_u, zqld, zqfr, zqfr_neg
    REAL(KIND = wp), PARAMETER :: zfric_umin = 0._wp
    REAL(KIND = wp), PARAMETER :: zch = 0.0057_wp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zu_io, zv_io, zfric
    IF (ln_timing) CALL timing_start('icethd')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icethd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_thd: sea-ice thermodynamics'
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    IF (ln_icedyn) THEN
      zu_io(:, :) = u_ice(:, :) - ssu_m(:, :)
      zv_io(:, :) = v_ice(:, :) - ssv_m(:, :)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfric(ji, jj) = rn_cio * (0.5_wp * (zu_io(ji, jj) * zu_io(ji, jj) + zu_io(ji - 1, jj) * zu_io(ji - 1, jj) + zv_io(ji, &
&jj) * zv_io(ji, jj) + zv_io(ji, jj - 1) * zv_io(ji, jj - 1))) * tmask(ji, jj, 1)
        END DO
      END DO
    ELSE
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfric(ji, jj) = r1_rau0 * SQRT(0.5_wp * (utau(ji, jj) * utau(ji, jj) + utau(ji - 1, jj) * utau(ji - 1, jj) + vtau(ji, &
&jj) * vtau(ji, jj) + vtau(ji, jj - 1) * vtau(ji, jj - 1))) * tmask(ji, jj, 1)
        END DO
      END DO
    END IF
    CALL lbc_lnk('icethd', zfric, 'T', 1.)
    DO jj = 1, jpj
      DO ji = 1, jpi
        rswitch = tmask(ji, jj, 1) * MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi10))
        zqld = tmask(ji, jj, 1) * rdt_ice * ((1._wp - at_i_b(ji, jj)) * qsr_oce(ji, jj) * frq_m(ji, jj) + (1._wp - at_i_b(ji, jj)) &
&* qns_oce(ji, jj) + qemp_oce(ji, jj))
        zqfr = rau0 * rcp * e3t_m(ji, jj) * (t_bo(ji, jj) - (sst_m(ji, jj) + rt0)) * tmask(ji, jj, 1)
        zqfr_neg = MIN(zqfr, 0._wp)
        zfric_u = MAX(SQRT(zfric(ji, jj)), zfric_umin)
        qsb_ice_bot(ji, jj) = rswitch * rau0 * rcp * zch * zfric_u * ((sst_m(ji, jj) + rt0) - t_bo(ji, jj))
        qsb_ice_bot(ji, jj) = rswitch * MIN(qsb_ice_bot(ji, jj), - zqfr_neg * r1_rdtice / MAX(at_i(ji, jj), epsi10))
        qlead(ji, jj) = MIN(0._wp, zqld - (qsb_ice_bot(ji, jj) * at_i(ji, jj) * rdt_ice) - zqfr)
        IF ((zqld >= 0._wp .AND. at_i(ji, jj) > 0._wp) .OR. at_i(ji, jj) >= (1._wp - epsi10)) THEN
          fhld(ji, jj) = rswitch * zqld * r1_rdtice / MAX(at_i(ji, jj), epsi10)
          qlead(ji, jj) = 0._wp
        ELSE
          fhld(ji, jj) = 0._wp
        END IF
        qt_atm_oi(ji, jj) = qns_tot(ji, jj) + qsr_tot(ji, jj)
      END DO
    END DO
    IF (.NOT. ln_icedO) qlead(:, :) = 0._wp
    IF (.NOT. ln_icedH) THEN
      qsb_ice_bot(:, :) = 0._wp
      fhld(:, :) = 0._wp
    END IF
    qt_oce_ai(:, :) = (1._wp - at_i_b(:, :)) * qns_oce(:, :) + qemp_oce(:, :) - qlead(:, :) * r1_rdtice - at_i(:, :) * &
&qsb_ice_bot(:, :) - at_i(:, :) * fhld(:, :)
    DO jl = 1, jpl
      npti = 0
      nptidx(:) = 0
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (a_i(ji, jj, jl) > epsi10) THEN
            npti = npti + 1
            nptidx(npti) = (jj - 1) * jpi + ji
          END IF
        END DO
      END DO
      IF (npti > 0) THEN
        CALL ice_thd_1d2d(jl, 1)
        s_i_new(1 : npti) = 0._wp
        dh_s_tot(1 : npti) = 0._wp
        dh_i_sum(1 : npti) = 0._wp
        dh_i_bom(1 : npti) = 0._wp
        dh_i_itm(1 : npti) = 0._wp
        dh_i_sub(1 : npti) = 0._wp
        dh_i_bog(1 : npti) = 0._wp
        dh_snowice(1 : npti) = 0._wp
        dh_s_mlt(1 : npti) = 0._wp
        CALL ice_thd_zdf
        IF (ln_icedH) THEN
          CALL ice_thd_dh
          CALL ice_thd_pnd
          CALL ice_thd_ent(e_i_1d(1 : npti, :))
        END IF
        CALL ice_thd_sal(ln_icedS)
        CALL ice_thd_temp
        IF (ln_icedH .AND. ln_virtual_itd) CALL ice_thd_mono
        IF (ln_icedA) CALL ice_thd_da
        CALL ice_thd_1d2d(jl, 2)
      END IF
    END DO
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'icethd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (jpl > 1) CALL ice_itd_rem(kt)
    IF (ln_icedO) CALL ice_thd_do
    IF (ln_icectl) CALL ice_prt(kt, iiceprt, jiceprt, 1, ' - ice thermodyn. - ')
    IF (ln_ctl) CALL ice_prt3D('icethd')
    IF (ln_timing) CALL timing_stop('icethd')
  END SUBROUTINE ice_thd
  SUBROUTINE ice_thd_temp
    INTEGER :: ji, jk
    REAL(KIND = wp) :: ztmelts, zbbb, zccc
    DO jk = 1, nlay_i
      DO ji = 1, npti
        ztmelts = - rTmlt * sz_i_1d(ji, jk)
        zbbb = (rcp - rcpi) * ztmelts + e_i_1d(ji, jk) * r1_rhoi - rLfus
        zccc = SQRT(MAX(zbbb * zbbb - 4._wp * rcpi * rLfus * ztmelts, 0._wp))
        t_i_1d(ji, jk) = rt0 - (zbbb + zccc) * 0.5_wp * r1_rcpi
        rswitch = 1._wp - MAX(0._wp, SIGN(1._wp, - h_i_1d(ji)))
        t_i_1d(ji, jk) = rswitch * t_i_1d(ji, jk) + (1._wp - rswitch) * rt0
      END DO
    END DO
  END SUBROUTINE ice_thd_temp
  SUBROUTINE ice_thd_mono
    INTEGER :: ji
    REAL(KIND = wp) :: zhi_bef
    REAL(KIND = wp) :: zdh_mel, zda_mel
    REAL(KIND = wp) :: zvi, zvs
    DO ji = 1, npti
      zdh_mel = MIN(0._wp, dh_i_itm(ji) + dh_i_sum(ji) + dh_i_bom(ji) + dh_snowice(ji) + dh_i_sub(ji))
      IF (zdh_mel < 0._wp .AND. a_i_1d(ji) > 0._wp) THEN
        zvi = a_i_1d(ji) * h_i_1d(ji)
        zvs = a_i_1d(ji) * h_s_1d(ji)
        zhi_bef = h_i_1d(ji) - zdh_mel
        rswitch = MAX(0._wp, SIGN(1._wp, zhi_bef - epsi20))
        zda_mel = rswitch * a_i_1d(ji) * zdh_mel / (2._wp * MAX(zhi_bef, epsi20))
        a_i_1d(ji) = MAX(epsi20, a_i_1d(ji) + zda_mel)
        h_i_1d(ji) = zvi / a_i_1d(ji)
        h_s_1d(ji) = zvs / a_i_1d(ji)
        at_i_1d(ji) = a_i_1d(ji)
      END IF
    END DO
  END SUBROUTINE ice_thd_mono
  SUBROUTINE ice_thd_1d2d(kl, kn)
    INTEGER, INTENT(IN) :: kl
    INTEGER, INTENT(IN) :: kn
    INTEGER :: jk
    SELECT CASE (kn)
    CASE (1)
      CALL tab_2d_1d(npti, nptidx(1 : npti), at_i_1d(1 : npti), at_i)
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), h_i_1d(1 : npti), h_i(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), h_s_1d(1 : npti), h_s(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), t_su_1d(1 : npti), t_su(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), s_i_1d(1 : npti), s_i(:, :, kl))
      DO jk = 1, nlay_s
        CALL tab_2d_1d(npti, nptidx(1 : npti), t_s_1d(1 : npti, jk), t_s(:, :, jk, kl))
        CALL tab_2d_1d(npti, nptidx(1 : npti), e_s_1d(1 : npti, jk), e_s(:, :, jk, kl))
      END DO
      DO jk = 1, nlay_i
        CALL tab_2d_1d(npti, nptidx(1 : npti), t_i_1d(1 : npti, jk), t_i(:, :, jk, kl))
        CALL tab_2d_1d(npti, nptidx(1 : npti), e_i_1d(1 : npti, jk), e_i(:, :, jk, kl))
        CALL tab_2d_1d(npti, nptidx(1 : npti), sz_i_1d(1 : npti, jk), sz_i(:, :, jk, kl))
      END DO
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_ip_1d(1 : npti), a_ip(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), h_ip_1d(1 : npti), h_ip(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), a_ip_frac_1d(1 : npti), a_ip_frac(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), qprec_ice_1d(1 : npti), qprec_ice)
      CALL tab_2d_1d(npti, nptidx(1 : npti), qsr_ice_1d(1 : npti), qsr_ice(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), qns_ice_1d(1 : npti), qns_ice(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), evap_ice_1d(1 : npti), evap_ice(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), dqns_ice_1d(1 : npti), dqns_ice(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), t_bo_1d(1 : npti), t_bo)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sprecip_1d(1 : npti), sprecip)
      CALL tab_2d_1d(npti, nptidx(1 : npti), qsb_ice_bot_1d(1 : npti), qsb_ice_bot)
      CALL tab_2d_1d(npti, nptidx(1 : npti), fhld_1d(1 : npti), fhld)
      CALL tab_2d_1d(npti, nptidx(1 : npti), qml_ice_1d(1 : npti), qml_ice(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), qcn_ice_1d(1 : npti), qcn_ice(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), qtr_ice_top_1d(1 : npti), qtr_ice_top(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_snw_sni_1d(1 : npti), wfx_snw_sni)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_snw_sum_1d(1 : npti), wfx_snw_sum)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_sub_1d(1 : npti), wfx_sub)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_snw_sub_1d(1 : npti), wfx_snw_sub)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_ice_sub_1d(1 : npti), wfx_ice_sub)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_err_sub_1d(1 : npti), wfx_err_sub)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_bog_1d(1 : npti), wfx_bog)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_bom_1d(1 : npti), wfx_bom)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_sum_1d(1 : npti), wfx_sum)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_sni_1d(1 : npti), wfx_sni)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_res_1d(1 : npti), wfx_res)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_spr_1d(1 : npti), wfx_spr)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_lam_1d(1 : npti), wfx_lam)
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_pnd_1d(1 : npti), wfx_pnd)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_bog_1d(1 : npti), sfx_bog)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_bom_1d(1 : npti), sfx_bom)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_sum_1d(1 : npti), sfx_sum)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_sni_1d(1 : npti), sfx_sni)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_bri_1d(1 : npti), sfx_bri)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_res_1d(1 : npti), sfx_res)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_sub_1d(1 : npti), sfx_sub)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_lam_1d(1 : npti), sfx_lam)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_thd_1d(1 : npti), hfx_thd)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_spr_1d(1 : npti), hfx_spr)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_sum_1d(1 : npti), hfx_sum)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_bom_1d(1 : npti), hfx_bom)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_bog_1d(1 : npti), hfx_bog)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_dif_1d(1 : npti), hfx_dif)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_opw_1d(1 : npti), hfx_opw)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_snw_1d(1 : npti), hfx_snw)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_sub_1d(1 : npti), hfx_sub)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_res_1d(1 : npti), hfx_res)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_err_dif_1d(1 : npti), hfx_err_dif)
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_err_rem_1d(1 : npti), hfx_err_rem)
      CALL tab_2d_1d(npti, nptidx(1 : npti), qt_oce_ai_1d(1 : npti), qt_oce_ai)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sst_1d(1 : npti), sst_m)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sss_1d(1 : npti), sss_m)
      CALL tab_2d_1d(npti, nptidx(1 : npti), o_i_1d(1 : npti), o_i(:, :, kl))
      CALL tab_2d_1d(npti, nptidx(1 : npti), oa_i_1d(1 : npti), oa_i(:, :, kl))
      DO jk = 1, nlay_i
        WHERE (h_i_1d(1 : npti) > 0._wp) e_i_1d(1 : npti, jk) = e_i_1d(1 : npti, jk) / (h_i_1d(1 : npti) * a_i_1d(1 : npti)) * &
&nlay_i
      END DO
      DO jk = 1, nlay_s
        WHERE (h_s_1d(1 : npti) > 0._wp) e_s_1d(1 : npti, jk) = e_s_1d(1 : npti, jk) / (h_s_1d(1 : npti) * a_i_1d(1 : npti)) * &
&nlay_s
      END DO
    CASE (2)
      DO jk = 1, nlay_i
        e_i_1d(1 : npti, jk) = e_i_1d(1 : npti, jk) * h_i_1d(1 : npti) * a_i_1d(1 : npti) * r1_nlay_i
      END DO
      DO jk = 1, nlay_s
        e_s_1d(1 : npti, jk) = e_s_1d(1 : npti, jk) * h_s_1d(1 : npti) * a_i_1d(1 : npti) * r1_nlay_s
      END DO
      v_i_1d(1 : npti) = h_i_1d(1 : npti) * a_i_1d(1 : npti)
      v_s_1d(1 : npti) = h_s_1d(1 : npti) * a_i_1d(1 : npti)
      sv_i_1d(1 : npti) = s_i_1d(1 : npti) * v_i_1d(1 : npti)
      v_ip_1d(1 : npti) = h_ip_1d(1 : npti) * a_ip_1d(1 : npti)
      oa_i_1d(1 : npti) = o_i_1d(1 : npti) * a_i_1d(1 : npti)
      CALL tab_1d_2d(npti, nptidx(1 : npti), at_i_1d(1 : npti), at_i)
      CALL tab_1d_2d(npti, nptidx(1 : npti), a_i_1d(1 : npti), a_i(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), h_i_1d(1 : npti), h_i(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), h_s_1d(1 : npti), h_s(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), t_su_1d(1 : npti), t_su(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), s_i_1d(1 : npti), s_i(:, :, kl))
      DO jk = 1, nlay_s
        CALL tab_1d_2d(npti, nptidx(1 : npti), t_s_1d(1 : npti, jk), t_s(:, :, jk, kl))
        CALL tab_1d_2d(npti, nptidx(1 : npti), e_s_1d(1 : npti, jk), e_s(:, :, jk, kl))
      END DO
      DO jk = 1, nlay_i
        CALL tab_1d_2d(npti, nptidx(1 : npti), t_i_1d(1 : npti, jk), t_i(:, :, jk, kl))
        CALL tab_1d_2d(npti, nptidx(1 : npti), e_i_1d(1 : npti, jk), e_i(:, :, jk, kl))
        CALL tab_1d_2d(npti, nptidx(1 : npti), sz_i_1d(1 : npti, jk), sz_i(:, :, jk, kl))
      END DO
      CALL tab_1d_2d(npti, nptidx(1 : npti), a_ip_1d(1 : npti), a_ip(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), h_ip_1d(1 : npti), h_ip(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), a_ip_frac_1d(1 : npti), a_ip_frac(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_snw_sni_1d(1 : npti), wfx_snw_sni)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_snw_sum_1d(1 : npti), wfx_snw_sum)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_sub_1d(1 : npti), wfx_sub)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_snw_sub_1d(1 : npti), wfx_snw_sub)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_ice_sub_1d(1 : npti), wfx_ice_sub)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_err_sub_1d(1 : npti), wfx_err_sub)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_bog_1d(1 : npti), wfx_bog)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_bom_1d(1 : npti), wfx_bom)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_sum_1d(1 : npti), wfx_sum)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_sni_1d(1 : npti), wfx_sni)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_res_1d(1 : npti), wfx_res)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_spr_1d(1 : npti), wfx_spr)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_lam_1d(1 : npti), wfx_lam)
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_pnd_1d(1 : npti), wfx_pnd)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_bog_1d(1 : npti), sfx_bog)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_bom_1d(1 : npti), sfx_bom)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_sum_1d(1 : npti), sfx_sum)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_sni_1d(1 : npti), sfx_sni)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_bri_1d(1 : npti), sfx_bri)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_res_1d(1 : npti), sfx_res)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_sub_1d(1 : npti), sfx_sub)
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_lam_1d(1 : npti), sfx_lam)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_thd_1d(1 : npti), hfx_thd)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_spr_1d(1 : npti), hfx_spr)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_sum_1d(1 : npti), hfx_sum)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_bom_1d(1 : npti), hfx_bom)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_bog_1d(1 : npti), hfx_bog)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_dif_1d(1 : npti), hfx_dif)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_opw_1d(1 : npti), hfx_opw)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_snw_1d(1 : npti), hfx_snw)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_sub_1d(1 : npti), hfx_sub)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_res_1d(1 : npti), hfx_res)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_err_dif_1d(1 : npti), hfx_err_dif)
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_err_rem_1d(1 : npti), hfx_err_rem)
      CALL tab_1d_2d(npti, nptidx(1 : npti), qt_oce_ai_1d(1 : npti), qt_oce_ai)
      CALL tab_1d_2d(npti, nptidx(1 : npti), qns_ice_1d(1 : npti), qns_ice(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), qtr_ice_bot_1d(1 : npti), qtr_ice_bot(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), cnd_ice_1d(1 : npti), cnd_ice(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), t1_ice_1d(1 : npti), t1_ice(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), t_si_1d(1 : npti), t_si(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), qcn_ice_bot_1d(1 : npti), qcn_ice_bot(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), qcn_ice_top_1d(1 : npti), qcn_ice_top(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), v_i_1d(1 : npti), v_i(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), v_s_1d(1 : npti), v_s(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), sv_i_1d(1 : npti), sv_i(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), v_ip_1d(1 : npti), v_ip(:, :, kl))
      CALL tab_1d_2d(npti, nptidx(1 : npti), oa_i_1d(1 : npti), oa_i(:, :, kl))
    END SELECT
  END SUBROUTINE ice_thd_1d2d
  SUBROUTINE ice_thd_init
    INTEGER :: ios
    NAMELIST /namthd/ ln_icedH, ln_icedA, ln_icedO, ln_icedS
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namthd, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namthd in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namthd, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namthd in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namthd)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_thd_init: Ice Thermodynamics'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namthd:'
      WRITE(numout, FMT = *) '      activate ice thick change from top/bot (T) or not (F)   ln_icedH  = ', ln_icedH
      WRITE(numout, FMT = *) '      activate lateral melting (T) or not (F)                 ln_icedA  = ', ln_icedA
      WRITE(numout, FMT = *) '      activate ice growth in open-water (T) or not (F)        ln_icedO  = ', ln_icedO
      WRITE(numout, FMT = *) '      activate gravity drainage and flushing (T) or not (F)   ln_icedS  = ', ln_icedS
    END IF
    CALL ice_thd_zdf_init
    IF (ln_icedA) CALL ice_thd_da_init
    IF (ln_icedO) CALL ice_thd_do_init
    CALL ice_thd_sal_init
    CALL ice_thd_pnd_init
  END SUBROUTINE ice_thd_init
END MODULE icethd