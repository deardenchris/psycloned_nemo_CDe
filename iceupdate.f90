MODULE iceupdate
  USE oce, ONLY: sshn, sshb
  USE phycst
  USE dom_oce
  USE ice
  USE sbc_ice
  USE sbc_oce
  USE sbccpl
  USE icealb
  USE traqsr
  USE icectl
  USE bdy_oce, ONLY: ln_bdy
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_update_init
  PUBLIC :: ice_update_flx
  PUBLIC :: ice_update_tau
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: utau_oce, vtau_oce
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: tmod_io
  CONTAINS
  INTEGER FUNCTION ice_update_alloc()
    ALLOCATE(utau_oce(jpi, jpj), vtau_oce(jpi, jpj), tmod_io(jpi, jpj), STAT = ice_update_alloc)
    CALL mpp_sum('iceupdate', ice_update_alloc)
    IF (ice_update_alloc /= 0) CALL ctl_stop('STOP', 'ice_update_alloc: failed to allocate arrays')
  END FUNCTION ice_update_alloc
  SUBROUTINE ice_update_flx(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jl, jk
    REAL(KIND = wp) :: zqmass
    REAL(KIND = wp) :: zqsr
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zalb_cs, zalb_os
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('ice_update_flx', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ice_update')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_update_flx: update fluxes (mass, salt and heat) at the ice-ocean interface'
      WRITE(numout, *) '~~~~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. ln_icethd) THEN
      !$ACC KERNELS
      qt_atm_oi(:, :) = (1._wp - at_i_b(:, :)) * (qns_oce(:, :) + qsr_oce(:, :)) + qemp_oce(:, :)
      qt_oce_ai(:, :) = (1._wp - at_i_b(:, :)) * qns_oce(:, :) + qemp_oce(:, :)
      emp_ice(:, :) = 0._wp
      qemp_ice(:, :) = 0._wp
      qevap_ice(:, :, :) = 0._wp
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zqsr = qsr_tot(ji, jj) - SUM(a_i_b(ji, jj, :) * (qsr_ice(ji, jj, :) - qtr_ice_bot(ji, jj, :)))
        zqmass = hfx_thd(ji, jj) + hfx_dyn(ji, jj) + hfx_res(ji, jj)
        qt_oce_ai(ji, jj) = qt_oce_ai(ji, jj) + zqmass + zqsr
        qt_oce_ai(ji, jj) = qt_oce_ai(ji, jj) + hfx_err_dif(ji, jj) + (hfx_sub(ji, jj) - SUM(qevap_ice(ji, jj, :) * a_i_b(ji, jj, :)))
        qsr(ji, jj) = zqsr
        qns(ji, jj) = qt_oce_ai(ji, jj) - zqsr
        wfx_sub(ji, jj) = wfx_snw_sub(ji, jj) + wfx_ice_sub(ji, jj)
        wfx_ice(ji, jj) = wfx_bog(ji, jj) + wfx_bom(ji, jj) + wfx_sum(ji, jj) + wfx_sni(ji, jj) + wfx_opw(ji, jj) + wfx_dyn(ji, jj) + wfx_res(ji, jj) + wfx_lam(ji, jj) + wfx_pnd(ji, jj)
        wfx_snw(ji, jj) = wfx_snw_sni(ji, jj) + wfx_snw_dyn(ji, jj) + wfx_snw_sum(ji, jj)
        fmmflx(ji, jj) = - (wfx_ice(ji, jj) + wfx_snw(ji, jj) + wfx_err_sub(ji, jj))
        emp(ji, jj) = emp_oce(ji, jj) - wfx_ice(ji, jj) - wfx_snw(ji, jj) - wfx_err_sub(ji, jj)
        sfx(ji, jj) = sfx_bog(ji, jj) + sfx_bom(ji, jj) + sfx_sum(ji, jj) + sfx_sni(ji, jj) + sfx_opw(ji, jj) + sfx_res(ji, jj) + sfx_dyn(ji, jj) + sfx_bri(ji, jj) + sfx_sub(ji, jj) + sfx_lam(ji, jj)
        snwice_mass_b(ji, jj) = snwice_mass(ji, jj)
        snwice_mass(ji, jj) = tmask(ji, jj, 1) * (rhos * vt_s(ji, jj) + rhoi * vt_i(ji, jj))
        snwice_fmass(ji, jj) = (snwice_mass(ji, jj) - snwice_mass_b(ji, jj)) * r1_rdtice
      END DO
    END DO
    fr_i(:, :) = at_i(:, :)
    tn_ice(:, :, :) = t_su(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_update_flx', 'r1', 0, 0)
    CALL ice_alb(t_su, h_i, h_s, ln_pnd_alb, a_ip_frac, h_ip, zalb_cs, zalb_os)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    alb_ice(:, :, :) = (1._wp - cldf_ice) * zalb_cs(:, :, :) + cldf_ice * zalb_os(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('ice_update_flx', 'r2', 0, 0)
    IF (lrst_ice) THEN
      CALL update_rst('WRITE', kt)
    END IF
    IF (iom_use('sfxice')) CALL iom_put("sfxice", sfx * 1.E-03)
    IF (iom_use('sfxbog')) CALL iom_put("sfxbog", sfx_bog * 1.E-03)
    IF (iom_use('sfxbom')) CALL iom_put("sfxbom", sfx_bom * 1.E-03)
    IF (iom_use('sfxsum')) CALL iom_put("sfxsum", sfx_sum * 1.E-03)
    IF (iom_use('sfxlam')) CALL iom_put("sfxlam", sfx_lam * 1.E-03)
    IF (iom_use('sfxsni')) CALL iom_put("sfxsni", sfx_sni * 1.E-03)
    IF (iom_use('sfxopw')) CALL iom_put("sfxopw", sfx_opw * 1.E-03)
    IF (iom_use('sfxdyn')) CALL iom_put("sfxdyn", sfx_dyn * 1.E-03)
    IF (iom_use('sfxbri')) CALL iom_put("sfxbri", sfx_bri * 1.E-03)
    IF (iom_use('sfxres')) CALL iom_put("sfxres", sfx_res * 1.E-03)
    IF (iom_use('sfxsub')) CALL iom_put("sfxsub", sfx_sub * 1.E-03)
    IF (iom_use('emp_oce')) CALL iom_put("emp_oce", emp_oce)
    IF (iom_use('emp_ice')) CALL iom_put("emp_ice", emp_ice)
    IF (iom_use('vfxice')) CALL iom_put("vfxice", wfx_ice)
    IF (iom_use('vfxbog')) CALL iom_put("vfxbog", wfx_bog)
    IF (iom_use('vfxbom')) CALL iom_put("vfxbom", wfx_bom)
    IF (iom_use('vfxsum')) CALL iom_put("vfxsum", wfx_sum)
    IF (iom_use('vfxlam')) CALL iom_put("vfxlam", wfx_lam)
    IF (iom_use('vfxsni')) CALL iom_put("vfxsni", wfx_sni)
    IF (iom_use('vfxopw')) CALL iom_put("vfxopw", wfx_opw)
    IF (iom_use('vfxdyn')) CALL iom_put("vfxdyn", wfx_dyn)
    IF (iom_use('vfxres')) CALL iom_put("vfxres", wfx_res)
    IF (iom_use('vfxpnd')) CALL iom_put("vfxpnd", wfx_pnd)
    IF (iom_use('vfxsub')) CALL iom_put("vfxsub", wfx_ice_sub)
    IF (iom_use('vfxsub_err')) CALL iom_put("vfxsub_err", wfx_err_sub)
    IF (iom_use("vfxthin")) THEN
      WHERE (hm_i(:, :) < 0.2 .AND. hm_i(:, :) > 0.)
        z2d = wfx_bog
      ELSEWHERE
        z2d = 0._wp
      END WHERE
      CALL iom_put("vfxthin", wfx_opw + z2d)
    END IF
    IF (iom_use('vfxsnw')) CALL iom_put("vfxsnw", wfx_snw)
    IF (iom_use('vfxsnw_sum')) CALL iom_put("vfxsnw_sum", wfx_snw_sum)
    IF (iom_use('vfxsnw_sni')) CALL iom_put("vfxsnw_sni", wfx_snw_sni)
    IF (iom_use('vfxsnw_dyn')) CALL iom_put("vfxsnw_dyn", wfx_snw_dyn)
    IF (iom_use('vfxsnw_sub')) CALL iom_put("vfxsnw_sub", wfx_snw_sub)
    IF (iom_use('vfxsnw_pre')) CALL iom_put("vfxsnw_pre", wfx_spr)
    IF (iom_use('qsr_oce')) CALL iom_put("qsr_oce", qsr_oce * (1._wp - at_i_b))
    IF (iom_use('qns_oce')) CALL iom_put("qns_oce", qns_oce * (1._wp - at_i_b) + qemp_oce)
    IF (iom_use('qsr_ice')) CALL iom_put("qsr_ice", SUM(qsr_ice * a_i_b, dim = 3))
    IF (iom_use('qns_ice')) CALL iom_put("qns_ice", SUM(qns_ice * a_i_b, dim = 3) + qemp_ice)
    IF (iom_use('qtr_ice_bot')) CALL iom_put("qtr_ice_bot", SUM(qtr_ice_bot * a_i_b, dim = 3))
    IF (iom_use('qtr_ice_top')) CALL iom_put("qtr_ice_top", SUM(qtr_ice_top * a_i_b, dim = 3))
    IF (iom_use('qt_oce')) CALL iom_put("qt_oce", (qsr_oce + qns_oce) * (1._wp - at_i_b) + qemp_oce)
    IF (iom_use('qt_ice')) CALL iom_put("qt_ice", SUM((qns_ice + qsr_ice) * a_i_b, dim = 3) + qemp_ice)
    IF (iom_use('qt_oce_ai')) CALL iom_put("qt_oce_ai", qt_oce_ai * tmask(:, :, 1))
    IF (iom_use('qt_atm_oi')) CALL iom_put("qt_atm_oi", qt_atm_oi * tmask(:, :, 1))
    IF (iom_use('qemp_oce')) CALL iom_put("qemp_oce", qemp_oce)
    IF (iom_use('qemp_ice')) CALL iom_put("qemp_ice", qemp_ice)
    IF (iom_use('hfxbog')) CALL iom_put("hfxbog", hfx_bog)
    IF (iom_use('hfxbom')) CALL iom_put("hfxbom", hfx_bom)
    IF (iom_use('hfxsum')) CALL iom_put("hfxsum", hfx_sum)
    IF (iom_use('hfxopw')) CALL iom_put("hfxopw", hfx_opw)
    IF (iom_use('hfxdif')) CALL iom_put("hfxdif", hfx_dif)
    IF (iom_use('hfxsnw')) CALL iom_put("hfxsnw", hfx_snw)
    IF (iom_use('hfxerr')) CALL iom_put("hfxerr", hfx_err_dif)
    IF (iom_use('hfxthd')) CALL iom_put("hfxthd", hfx_thd)
    IF (iom_use('hfxdyn')) CALL iom_put("hfxdyn", hfx_dyn)
    IF (iom_use('hfxres')) CALL iom_put("hfxres", hfx_res)
    IF (iom_use('hfxsub')) CALL iom_put("hfxsub", hfx_sub)
    IF (iom_use('hfxspr')) CALL iom_put("hfxspr", hfx_spr)
    IF (iom_use('hfxsensib')) CALL iom_put("hfxsensib", - qsb_ice_bot * at_i_b)
    IF (iom_use('hfxcndbot')) CALL iom_put("hfxcndbot", SUM(qcn_ice_bot * a_i_b, dim = 3))
    IF (iom_use('hfxcndtop')) CALL iom_put("hfxcndtop", SUM(qcn_ice_top * a_i_b, dim = 3))
    IF (iom_use('hfxdhc')) CALL iom_put("hfxdhc", diag_heat)
    IF (ln_icediachk .AND. .NOT. ln_bdy) CALL ice_cons_final('iceupdate')
    IF (ln_icectl) CALL ice_prt(kt, iiceprt, jiceprt, 3, 'Final state ice_update')
    IF (ln_ctl) CALL ice_prt3D('iceupdate')
    IF (ln_timing) CALL timing_stop('ice_update')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE ice_update_flx
  SUBROUTINE ice_update_tau(kt, pu_oce, pv_oce)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pu_oce, pv_oce
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zat_u, zutau_ice, zu_t, zmodt
    REAL(KIND = wp) :: zat_v, zvtau_ice, zv_t, zrhoco
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('ice_update_tau', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ice_update_tau')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_update_tau: update stress at the ice-ocean interface'
      WRITE(numout, *) '~~~~~~~~~~~~~~'
    END IF
    zrhoco = rau0 * rn_cio
    CALL profile_psy_data0 % PostEnd
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu_t = u_ice(ji, jj) + u_ice(ji - 1, jj) - u_oce(ji, jj) - u_oce(ji - 1, jj)
          zv_t = v_ice(ji, jj) + v_ice(ji, jj - 1) - v_oce(ji, jj) - v_oce(ji, jj - 1)
          zmodt = 0.25_wp * (zu_t * zu_t + zv_t * zv_t)
          taum(ji, jj) = (1._wp - at_i(ji, jj)) * taum(ji, jj) + at_i(ji, jj) * zrhoco * zmodt
          tmod_io(ji, jj) = zrhoco * SQRT(zmodt)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('ice_update_tau', 'r1', 0, 0)
      CALL lbc_lnk_multi('iceupdate', taum, 'T', 1., tmod_io, 'T', 1.)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      utau_oce(:, :) = utau(:, :)
      vtau_oce(:, :) = vtau(:, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zat_u = (at_i(ji, jj) * tmask(ji, jj, 1) + at_i(ji + 1, jj) * tmask(ji + 1, jj, 1)) / MAX(1.0_wp, tmask(ji, jj, 1) + tmask(ji + 1, jj, 1))
        zat_v = (at_i(ji, jj) * tmask(ji, jj, 1) + at_i(ji, jj + 1) * tmask(ji, jj + 1, 1)) / MAX(1.0_wp, tmask(ji, jj, 1) + tmask(ji, jj + 1, 1))
        zutau_ice = 0.5_wp * (tmod_io(ji, jj) + tmod_io(ji + 1, jj)) * (u_ice(ji, jj) - pu_oce(ji, jj))
        zvtau_ice = 0.5_wp * (tmod_io(ji, jj) + tmod_io(ji, jj + 1)) * (v_ice(ji, jj) - pv_oce(ji, jj))
        utau(ji, jj) = (1._wp - zat_u) * utau_oce(ji, jj) + zat_u * zutau_ice
        vtau(ji, jj) = (1._wp - zat_v) * vtau_oce(ji, jj) + zat_v * zvtau_ice
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('ice_update_tau', 'r2', 0, 0)
    CALL lbc_lnk_multi('iceupdate', utau, 'U', - 1., vtau, 'V', - 1.)
    IF (ln_timing) CALL timing_stop('ice_update_tau')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE ice_update_tau
  SUBROUTINE ice_update_init
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoefu, zcoefv, zcoeff
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'ice_update_init: ice-ocean stress init'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~~'
    IF (ice_update_alloc() /= 0) CALL ctl_stop('STOP', 'ice_update_init : unable to allocate standard arrays')
    CALL update_rst('READ')
  END SUBROUTINE ice_update_init
  SUBROUTINE update_rst(cdrw, kt)
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER, OPTIONAL, INTENT(IN) :: kt
    INTEGER :: iter
    INTEGER :: id1
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numrir, 'snwice_mass', ldstop = .FALSE.)
        IF (id1 > 0) THEN
          CALL iom_get(numrir, jpdom_autoglo, 'snwice_mass', snwice_mass)
          CALL iom_get(numrir, jpdom_autoglo, 'snwice_mass_b', snwice_mass_b)
        ELSE
          IF (lwp) WRITE(numout, *) '   ==>>   previous run without snow-ice mass output then set it'
          !$ACC KERNELS
          snwice_mass(:, :) = tmask(:, :, 1) * (rhos * vt_s(:, :) + rhoi * vt_i(:, :))
          snwice_mass_b(:, :) = snwice_mass(:, :)
          !$ACC END KERNELS
        END IF
      ELSE
        IF (lwp) WRITE(numout, *) '   ==>>   start from rest: set the snow-ice mass'
        !$ACC KERNELS
        snwice_mass(:, :) = tmask(:, :, 1) * (rhos * vt_s(:, :) + rhoi * vt_i(:, :))
        snwice_mass_b(:, :) = snwice_mass(:, :)
        !$ACC END KERNELS
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, *) '---- update-rst ----'
      iter = kt + nn_fsbc - 1
      CALL iom_rstput(iter, nitrst, numriw, 'snwice_mass', snwice_mass)
      CALL iom_rstput(iter, nitrst, numriw, 'snwice_mass_b', snwice_mass_b)
    END IF
  END SUBROUTINE update_rst
END MODULE iceupdate