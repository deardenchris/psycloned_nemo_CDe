MODULE crsfld
  USE crs
  USE crsdom
  USE crslbclnk
  USE oce
  USE dom_oce
  USE sbc_oce
  USE zdf_oce
  USE ldftra
  USE zdfddm
  USE in_out_manager
  USE iom
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: crs_fld
  CONTAINS
  SUBROUTINE crs_fld(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: z2dcrsu, z2dcrsv
    REAL(KIND = wp) :: zztmp
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ze3t, ze3u, ze3v, ze3w
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zt, zs, z3d
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk) :: zt_crs, zs_crs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    IF (ln_timing) CALL timing_start('crs_fld')
    !$ACC KERNELS
    ze3t(:, :, :) = e3t_n(:, :, :)
    ze3u(:, :, :) = e3u_n(:, :, :)
    ze3v(:, :, :) = e3v_n(:, :, :)
    ze3w(:, :, :) = e3w_n(:, :, :)
    IF (kt == nit000) THEN
      tsn_crs(:, :, :, :) = 0._wp
      un_crs(:, :, :) = 0._wp
      vn_crs(:, :, :) = 0._wp
      wn_crs(:, :, :) = 0._wp
      avs_crs(:, :, :) = 0._wp
      hdivn_crs(:, :, :) = 0._wp
      sshn_crs(:, :) = 0._wp
      utau_crs(:, :) = 0._wp
      vtau_crs(:, :) = 0._wp
      wndm_crs(:, :) = 0._wp
      qsr_crs(:, :) = 0._wp
      emp_crs(:, :) = 0._wp
      emp_b_crs(:, :) = 0._wp
      rnf_crs(:, :) = 0._wp
      fr_i_crs(:, :) = 0._wp
    END IF
    !$ACC END KERNELS
    CALL iom_swap("nemo_crs")
    !$ACC KERNELS
    zt(:, :, :) = tsn(:, :, :, jp_tem)
    zt_crs(:, :, :) = 0._wp
    !$ACC END KERNELS
    CALL crs_dom_ope(zt, 'VOL', 'T', tmask, zt_crs, p_e12 = e1e2t, p_e3 = ze3t, psgn = 1.0)
    !$ACC KERNELS
    tsn_crs(:, :, :, jp_tem) = zt_crs(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('crs_fld', 'r0', 0, 0)
    CALL iom_put("toce", tsn_crs(:, :, :, jp_tem))
    CALL iom_put("sst", tsn_crs(:, :, 1, jp_tem))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zs(:, :, :) = tsn(:, :, :, jp_sal)
    zs_crs(:, :, :) = 0._wp
    !$ACC END KERNELS
    CALL crs_dom_ope(zs, 'VOL', 'T', tmask, zs_crs, p_e12 = e1e2t, p_e3 = ze3t, psgn = 1.0)
    !$ACC KERNELS
    tsn_crs(:, :, :, jp_sal) = zt_crs(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('crs_fld', 'r1', 0, 0)
    CALL iom_put("soce", tsn_crs(:, :, :, jp_sal))
    CALL iom_put("sss", tsn_crs(:, :, 1, jp_sal))
    CALL crs_dom_ope(un, 'SUM', 'U', umask, un_crs, p_e12 = e2u, p_e3 = ze3u, p_surf_crs = e2e3u_msk, psgn = - 1.0)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    zt(:, :, :) = 0._wp
    zs(:, :, :) = 0._wp
    zt_crs(:, :, :) = 0._wp
    zs_crs(:, :, :) = 0._wp
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zt(ji, jj, jk) = un(ji, jj, jk) * 0.5 * (tsn(ji, jj, jk, jp_tem) + tsn(ji + 1, jj, jk, jp_tem))
          zs(ji, jj, jk) = un(ji, jj, jk) * 0.5 * (tsn(ji, jj, jk, jp_sal) + tsn(ji + 1, jj, jk, jp_sal))
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('crs_fld', 'r2', 0, 0)
    CALL crs_dom_ope(zt, 'SUM', 'U', umask, zt_crs, p_e12 = e2u, p_e3 = ze3u, p_surf_crs = e2e3u_msk, psgn = - 1.0)
    CALL crs_dom_ope(zs, 'SUM', 'U', umask, zs_crs, p_e12 = e2u, p_e3 = ze3u, p_surf_crs = e2e3u_msk, psgn = - 1.0)
    CALL iom_put("uoce", un_crs)
    CALL iom_put("uocet", zt_crs)
    CALL iom_put("uoces", zs_crs)
    CALL crs_dom_ope(vn, 'SUM', 'V', vmask, vn_crs, p_e12 = e1v, p_e3 = ze3v, p_surf_crs = e1e3v_msk, psgn = - 1.0)
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    zt(:, :, :) = 0._wp
    zs(:, :, :) = 0._wp
    zt_crs(:, :, :) = 0._wp
    zs_crs(:, :, :) = 0._wp
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zt(ji, jj, jk) = vn(ji, jj, jk) * 0.5 * (tsn(ji, jj, jk, jp_tem) + tsn(ji, jj + 1, jk, jp_tem))
          zs(ji, jj, jk) = vn(ji, jj, jk) * 0.5 * (tsn(ji, jj, jk, jp_sal) + tsn(ji, jj + 1, jk, jp_sal))
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('crs_fld', 'r3', 0, 0)
    CALL crs_dom_ope(zt, 'SUM', 'V', vmask, zt_crs, p_e12 = e1v, p_e3 = ze3v, p_surf_crs = e1e3v_msk, psgn = - 1.0)
    CALL crs_dom_ope(zs, 'SUM', 'V', vmask, zs_crs, p_e12 = e1v, p_e3 = ze3v, p_surf_crs = e1e3v_msk, psgn = - 1.0)
    CALL iom_put("voce", vn_crs)
    CALL iom_put("vocet", zt_crs)
    CALL iom_put("voces", zs_crs)
    CALL profile_psy_data3 % PostEnd
    IF (iom_use("eken")) THEN
      !$ACC KERNELS
      z3d(:, :, jk) = 0._wp
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zztmp = r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
            z3d(ji, jj, jk) = 0.25_wp * zztmp * (un(ji - 1, jj, jk) ** 2 * e2u(ji - 1, jj) * e3u_n(ji - 1, jj, jk) + un(ji, jj, &
&jk) ** 2 * e2u(ji, jj) * e3u_n(ji, jj, jk) + vn(ji, jj - 1, jk) ** 2 * e1v(ji, jj - 1) * e3v_n(ji, jj - 1, jk) + vn(ji, jj, jk) &
&** 2 * e1v(ji, jj) * e3v_n(ji, jj, jk))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data4 % PreStart('crs_fld', 'r4', 0, 0)
      CALL lbc_lnk('crsfld', z3d, 'T', 1.)
      CALL crs_dom_ope(z3d, 'VOL', 'T', tmask, zt_crs, p_e12 = e1e2t, p_e3 = ze3t, psgn = 1.0)
      CALL iom_put("eken", zt_crs)
      CALL profile_psy_data4 % PostEnd
    END IF
    !$ACC KERNELS
    DO jk = 1, jpkm1
      DO ji = 2, jpi_crsm1
        DO jj = 2, jpj_crsm1
          IF (tmask_crs(ji, jj, jk) > 0) THEN
            z2dcrsu = (un_crs(ji, jj, jk) * crs_surfu_wgt(ji, jj, jk)) - (un_crs(ji - 1, jj, jk) * crs_surfu_wgt(ji - 1, jj, jk))
            z2dcrsv = (vn_crs(ji, jj, jk) * crs_surfv_wgt(ji, jj, jk)) - (vn_crs(ji, jj - 1, jk) * crs_surfv_wgt(ji, jj - 1, jk))
            hdivn_crs(ji, jj, jk) = (z2dcrsu + z2dcrsv) / crs_volt_wgt(ji, jj, jk)
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data5 % PreStart('crs_fld', 'r5', 0, 0)
    CALL crs_lbc_lnk(hdivn_crs, 'T', 1.0)
    CALL iom_put("hdiv", hdivn_crs)
    CALL profile_psy_data5 % PostEnd
    IF (ln_crs_wn) THEN
      CALL crs_dom_ope(wn, 'SUM', 'W', tmask, wn_crs, p_e12 = e1e2t, p_surf_crs = e1e2w_msk, psgn = 1.0)
    ELSE
      !$ACC KERNELS
      wn_crs(:, :, jpk) = 0._wp
      DO jk = jpkm1, 1, - 1
        wn_crs(:, :, jk) = wn_crs(:, :, jk + 1) - e3t_crs(:, :, jk) * hdivn_crs(:, :, jk)
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data6 % PreStart('crs_fld', 'r6', 0, 0)
    CALL iom_put("woce", wn_crs)
    SELECT CASE (nn_crs_kz)
    CASE (0)
      CALL crs_dom_ope(avt, 'VOL', 'W', tmask, avt_crs, p_e12 = e1e2t, p_e3 = ze3w, psgn = 1.0)
      CALL crs_dom_ope(avs, 'VOL', 'W', tmask, avs_crs, p_e12 = e1e2t, p_e3 = ze3w, psgn = 1.0)
    CASE (1)
      CALL crs_dom_ope(avt, 'MAX', 'W', tmask, avt_crs, p_e12 = e1e2t, p_e3 = ze3w, psgn = 1.0)
      CALL crs_dom_ope(avs, 'MAX', 'W', tmask, avs_crs, p_e12 = e1e2t, p_e3 = ze3w, psgn = 1.0)
    CASE (2)
      CALL crs_dom_ope(avt, 'MIN', 'W', tmask, avt_crs, p_e12 = e1e2t, p_e3 = ze3w, psgn = 1.0)
      CALL crs_dom_ope(avs, 'MIN', 'W', tmask, avs_crs, p_e12 = e1e2t, p_e3 = ze3w, psgn = 1.0)
    END SELECT
    CALL iom_put("avt", avt_crs)
    CALL iom_put("avs", avs_crs)
    CALL crs_dom_ope(sshn, 'VOL', 'T', tmask, sshn_crs, p_e12 = e1e2t, p_e3 = ze3t, psgn = 1.0)
    CALL crs_dom_ope(utau, 'SUM', 'U', umask, utau_crs, p_e12 = e2u, p_surf_crs = e2u_crs, psgn = 1.0)
    CALL crs_dom_ope(vtau, 'SUM', 'V', vmask, vtau_crs, p_e12 = e1v, p_surf_crs = e1v_crs, psgn = 1.0)
    CALL crs_dom_ope(wndm, 'SUM', 'T', tmask, wndm_crs, p_e12 = e1e2t, p_surf_crs = e1e2t_crs, psgn = 1.0)
    CALL crs_dom_ope(rnf, 'MAX', 'T', tmask, rnf_crs, psgn = 1.0)
    CALL crs_dom_ope(qsr, 'SUM', 'T', tmask, qsr_crs, p_e12 = e1e2t, p_surf_crs = e1e2t_crs, psgn = 1.0)
    CALL crs_dom_ope(emp_b, 'SUM', 'T', tmask, emp_b_crs, p_e12 = e1e2t, p_surf_crs = e1e2t_crs, psgn = 1.0)
    CALL crs_dom_ope(emp, 'SUM', 'T', tmask, emp_crs, p_e12 = e1e2t, p_surf_crs = e1e2t_crs, psgn = 1.0)
    CALL crs_dom_ope(sfx, 'SUM', 'T', tmask, sfx_crs, p_e12 = e1e2t, p_surf_crs = e1e2t_crs, psgn = 1.0)
    CALL crs_dom_ope(fr_i, 'SUM', 'T', tmask, fr_i_crs, p_e12 = e1e2t, p_surf_crs = e1e2t_crs, psgn = 1.0)
    CALL iom_put("ssh", sshn_crs)
    CALL iom_put("utau", utau_crs)
    CALL iom_put("vtau", vtau_crs)
    CALL iom_put("wspd", wndm_crs)
    CALL iom_put("runoffs", rnf_crs)
    CALL iom_put("qsr", qsr_crs)
    CALL iom_put("empmr", emp_crs)
    CALL iom_put("saltflx", sfx_crs)
    CALL iom_put("ice_cover", fr_i_crs)
    CALL iom_swap("nemo")
    IF (ln_timing) CALL timing_stop('crs_fld')
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE crs_fld
END MODULE crsfld