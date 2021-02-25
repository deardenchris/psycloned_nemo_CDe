MODULE icewri
  USE dianam
  USE phycst
  USE dom_oce
  USE sbc_oce
  USE sbc_ice
  USE ice
  USE icevar
  USE ioipsl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_wri
  PUBLIC :: ice_wri_state
  CONTAINS
  SUBROUTINE ice_wri(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: z2da, z2db, zrho1, zrho2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d, zfast
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmsk00, zmsk05, zmsk15, zmsksn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zmsk00l, zmsksnl
    REAL(KIND = wp) :: zdiag_area_nh, zdiag_extt_nh, zdiag_volu_nh
    REAL(KIND = wp) :: zdiag_area_sh, zdiag_extt_sh, zdiag_volu_sh
    IF (ln_timing) CALL timing_start('icewri')
    CALL ice_var_bv
    DO jj = 1, jpj
      DO ji = 1, jpi
        zmsk00(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi06))
        zmsk05(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - 0.05_wp))
        zmsk15(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - 0.15_wp))
        zmsksn(ji, jj) = MAX(0._wp, SIGN(1._wp, vt_s(ji, jj) - epsi06))
      END DO
    END DO
    DO jl = 1, jpl
      DO jj = 1, jpj
        DO ji = 1, jpi
          zmsk00l(ji, jj, jl) = MAX(0._wp, SIGN(1._wp, a_i(ji, jj, jl) - epsi06))
          zmsksnl(ji, jj, jl) = MAX(0._wp, SIGN(1._wp, v_s(ji, jj, jl) - epsi06))
        END DO
      END DO
    END DO
    zrho1 = (rau0 - rhoi) * r1_rau0
    zrho2 = rhos * r1_rau0
    IF (iom_use('icemask')) CALL iom_put("icemask", zmsk00)
    IF (iom_use('icemask05')) CALL iom_put("icemask05", zmsk05)
    IF (iom_use('icemask15')) CALL iom_put("icemask15", zmsk15)
    IF (iom_use('icemass')) CALL iom_put("icemass", rhoi * vt_i * zmsk00)
    IF (iom_use('snwmass')) CALL iom_put("snwmass", rhos * vt_s * zmsksn)
    IF (iom_use('icepres')) CALL iom_put("icepres", zmsk00)
    IF (iom_use('iceconc')) CALL iom_put("iceconc", at_i * zmsk00)
    IF (iom_use('icevolu')) CALL iom_put("icevolu", vt_i * zmsk00)
    IF (iom_use('icethic')) CALL iom_put("icethic", hm_i * zmsk00)
    IF (iom_use('snwthic')) CALL iom_put("snwthic", hm_s * zmsk00)
    IF (iom_use('icebrv')) CALL iom_put("icebrv", bvm_i * zmsk00 * 100.)
    IF (iom_use('iceage')) CALL iom_put("iceage", om_i * zmsk15 / rday)
    IF (iom_use('icehnew')) CALL iom_put("icehnew", ht_i_new)
    IF (iom_use('snwvolu')) CALL iom_put("snwvolu", vt_s * zmsksn)
    IF (iom_use('icefrb')) THEN
      z2d(:, :) = (zrho1 * hm_i(:, :) - zrho2 * hm_s(:, :))
      WHERE (z2d < 0._wp) z2d = 0._wp
      CALL iom_put("icefrb", z2d * zmsk00)
    END IF
    IF (iom_use('iceapnd')) CALL iom_put("iceapnd", at_ip * zmsk00)
    IF (iom_use('icevpnd')) CALL iom_put("icevpnd", vt_ip * zmsk00)
    IF (iom_use('icesalt')) CALL iom_put("icesalt", sm_i * zmsk00)
    IF (iom_use('icesalm')) CALL iom_put("icesalm", SUM(sv_i, DIM = 3) * rhoi * 1.0E-3 * zmsk00)
    IF (iom_use('icetemp')) CALL iom_put("icetemp", (tm_i - rt0) * zmsk00)
    IF (iom_use('snwtemp')) CALL iom_put("snwtemp", (tm_s - rt0) * zmsksn)
    IF (iom_use('icettop')) CALL iom_put("icettop", (tm_su - rt0) * zmsk00)
    IF (iom_use('icetbot')) CALL iom_put("icetbot", (t_bo - rt0) * zmsk00)
    IF (iom_use('icetsni')) CALL iom_put("icetsni", (tm_si - rt0) * zmsk00)
    IF (iom_use('icehc')) CALL iom_put("icehc", - et_i * zmsk00)
    IF (iom_use('snwhc')) CALL iom_put("snwhc", - et_s * zmsksn)
    IF (iom_use('uice')) CALL iom_put("uice", u_ice)
    IF (iom_use('vice')) CALL iom_put("vice", v_ice)
    IF (iom_use('utau_ai')) CALL iom_put("utau_ai", utau_ice * zmsk00)
    IF (iom_use('vtau_ai')) CALL iom_put("vtau_ai", vtau_ice * zmsk00)
    IF (iom_use('icevel') .OR. iom_use('fasticepres')) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          z2da = (u_ice(ji, jj) + u_ice(ji - 1, jj))
          z2db = (v_ice(ji, jj) + v_ice(ji, jj - 1))
          z2d(ji, jj) = 0.5_wp * SQRT(z2da * z2da + z2db * z2db)
        END DO
      END DO
      CALL lbc_lnk('icewri', z2d, 'T', 1.)
      IF (iom_use('icevel')) CALL iom_put("icevel", z2d)
      WHERE (z2d(:, :) < 5.E-04_wp .AND. zmsk15(:, :) == 1._wp)
        zfast(:, :) = 1._wp
      ELSEWHERE
        zfast(:, :) = 0._wp
      END WHERE
      IF (iom_use('fasticepres')) CALL iom_put("fasticepres", zfast)
    END IF
    IF (iom_use('icemask_cat')) CALL iom_put("icemask_cat", zmsk00l)
    IF (iom_use('iceconc_cat')) CALL iom_put("iceconc_cat", a_i * zmsk00l)
    IF (iom_use('icethic_cat')) CALL iom_put("icethic_cat", h_i * zmsk00l)
    IF (iom_use('snwthic_cat')) CALL iom_put("snwthic_cat", h_s * zmsksnl)
    IF (iom_use('icesalt_cat')) CALL iom_put("icesalt_cat", s_i * zmsk00l)
    IF (iom_use('iceage_cat')) CALL iom_put("iceage_cat", o_i * zmsk00l / rday)
    IF (iom_use('icetemp_cat')) CALL iom_put("icetemp_cat", (SUM(t_i(:, :, :, :), dim = 3) * r1_nlay_i - rt0) * zmsk00l)
    IF (iom_use('snwtemp_cat')) CALL iom_put("snwtemp_cat", (SUM(t_s(:, :, :, :), dim = 3) * r1_nlay_s - rt0) * zmsksnl)
    IF (iom_use('icettop_cat')) CALL iom_put("icettop_cat", (t_su - rt0) * zmsk00l)
    IF (iom_use('icebrv_cat')) CALL iom_put("icebrv_cat", bv_i * 100. * zmsk00l)
    IF (iom_use('iceapnd_cat')) CALL iom_put("iceapnd_cat", a_ip * zmsk00l)
    IF (iom_use('icehpnd_cat')) CALL iom_put("icehpnd_cat", h_ip * zmsk00l)
    IF (iom_use('iceafpnd_cat')) CALL iom_put("iceafpnd_cat", a_ip_frac * zmsk00l)
    IF (iom_use('dmithd')) CALL iom_put("dmithd", - wfx_bog - wfx_bom - wfx_sum - wfx_sni - wfx_opw - wfx_lam - wfx_res)
    IF (iom_use('dmidyn')) CALL iom_put("dmidyn", - wfx_dyn + rhoi * diag_trp_vi)
    IF (iom_use('dmiopw')) CALL iom_put("dmiopw", - wfx_opw)
    IF (iom_use('dmibog')) CALL iom_put("dmibog", - wfx_bog)
    IF (iom_use('dmisni')) CALL iom_put("dmisni", - wfx_sni)
    IF (iom_use('dmisum')) CALL iom_put("dmisum", - wfx_sum)
    IF (iom_use('dmibom')) CALL iom_put("dmibom", - wfx_bom)
    IF (iom_use('dmtsub')) CALL iom_put("dmtsub", - wfx_sub)
    IF (iom_use('dmssub')) CALL iom_put("dmssub", - wfx_snw_sub)
    IF (iom_use('dmisub')) CALL iom_put("dmisub", - wfx_ice_sub)
    IF (iom_use('dmsspr')) CALL iom_put("dmsspr", - wfx_spr)
    IF (iom_use('dmsssi')) CALL iom_put("dmsssi", wfx_sni * rhos * r1_rhoi)
    IF (iom_use('dmsmel')) CALL iom_put("dmsmel", - wfx_snw_sum)
    IF (iom_use('dmsdyn')) CALL iom_put("dmsdyn", - wfx_snw_dyn + rhos * diag_trp_vs)
    IF (iom_use('NH_icearea') .OR. iom_use('NH_icevolu') .OR. iom_use('NH_iceextt')) THEN
      WHERE (ff_t > 0._wp)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_area_nh = glob_sum('icewri', at_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      zdiag_volu_nh = glob_sum('icewri', vt_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      WHERE (ff_t > 0._wp .AND. at_i > 0.15)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_extt_nh = glob_sum('icewri', zmsk00(:, :) * e1e2t(:, :))
      IF (iom_use('NH_icearea')) CALL iom_put("NH_icearea", zdiag_area_nh)
      IF (iom_use('NH_icevolu')) CALL iom_put("NH_icevolu", zdiag_volu_nh)
      IF (iom_use('NH_iceextt')) CALL iom_put("NH_iceextt", zdiag_extt_nh)
    END IF
    IF (iom_use('SH_icearea') .OR. iom_use('SH_icevolu') .OR. iom_use('SH_iceextt')) THEN
      WHERE (ff_t < 0._wp)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_area_sh = glob_sum('icewri', at_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      zdiag_volu_sh = glob_sum('icewri', vt_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      WHERE (ff_t < 0._wp .AND. at_i > 0.15)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_extt_sh = glob_sum('icewri', zmsk00(:, :) * e1e2t(:, :))
      IF (iom_use('SH_icearea')) CALL iom_put("SH_icearea", zdiag_area_sh)
      IF (iom_use('SH_icevolu')) CALL iom_put("SH_icevolu", zdiag_volu_sh)
      IF (iom_use('SH_iceextt')) CALL iom_put("SH_iceextt", zdiag_extt_sh)
    END IF
    IF (ln_timing) CALL timing_stop('icewri')
  END SUBROUTINE ice_wri
  SUBROUTINE ice_wri_state(kid)
    INTEGER, INTENT(IN) :: kid
    CALL iom_rstput(0, 0, kid, 'sithic', hm_i)
    CALL iom_rstput(0, 0, kid, 'siconc', at_i)
    CALL iom_rstput(0, 0, kid, 'sitemp', tm_i - rt0)
    CALL iom_rstput(0, 0, kid, 'sivelu', u_ice)
    CALL iom_rstput(0, 0, kid, 'sivelv', v_ice)
    CALL iom_rstput(0, 0, kid, 'sistru', utau_ice)
    CALL iom_rstput(0, 0, kid, 'sistrv', vtau_ice)
    CALL iom_rstput(0, 0, kid, 'sisflx', qsr)
    CALL iom_rstput(0, 0, kid, 'sinflx', qns)
    CALL iom_rstput(0, 0, kid, 'snwpre', sprecip)
    CALL iom_rstput(0, 0, kid, 'sisali', sm_i)
    CALL iom_rstput(0, 0, kid, 'sivolu', vt_i)
    CALL iom_rstput(0, 0, kid, 'sidive', divu_i * 1.0E8)
    CALL iom_rstput(0, 0, kid, 'si_amp', at_ip)
    CALL iom_rstput(0, 0, kid, 'si_vmp', vt_ip)
    CALL iom_rstput(0, 0, kid, 'sithicat', h_i)
    CALL iom_rstput(0, 0, kid, 'siconcat', a_i)
    CALL iom_rstput(0, 0, kid, 'sisalcat', s_i)
    CALL iom_rstput(0, 0, kid, 'snthicat', h_s)
  END SUBROUTINE ice_wri_state
END MODULE icewri