MODULE domvvl
  USE oce
  USE phycst
  USE dom_oce
  USE sbc_oce
  USE wet_dry
  USE usrdef_istate
  USE restart
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_vvl_init
  PUBLIC :: dom_vvl_sf_nxt
  PUBLIC :: dom_vvl_sf_swp
  PUBLIC :: dom_vvl_interpol
  LOGICAL, PUBLIC :: ln_vvl_zstar = .FALSE.
  LOGICAL, PUBLIC :: ln_vvl_ztilde = .FALSE.
  LOGICAL, PUBLIC :: ln_vvl_layer = .FALSE.
  LOGICAL, PUBLIC :: ln_vvl_ztilde_as_zstar = .FALSE.
  LOGICAL, PUBLIC :: ln_vvl_zstar_at_eqtor = .FALSE.
  LOGICAL, PUBLIC :: ln_vvl_kepe = .FALSE.
  REAL(KIND = wp) :: rn_ahe3
  REAL(KIND = wp) :: rn_rst_e3t
  REAL(KIND = wp) :: rn_lf_cutoff
  REAL(KIND = wp) :: rn_zdef_max
  LOGICAL, PUBLIC :: ln_vvl_dbg = .FALSE.
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: un_td, vn_td
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: hdiv_lf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: tilde_e3t_b, tilde_e3t_n
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: tilde_e3t_a, dtilde_e3t_a
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: frq_rst_e3t
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: frq_rst_hdv
  CONTAINS
  INTEGER FUNCTION dom_vvl_alloc()
    IF (ln_vvl_zstar) dom_vvl_alloc = 0
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      ALLOCATE(tilde_e3t_b(jpi, jpj, jpk), tilde_e3t_n(jpi, jpj, jpk), tilde_e3t_a(jpi, jpj, jpk), dtilde_e3t_a(jpi, jpj, jpk), un_td(jpi, jpj, jpk), vn_td(jpi, jpj, jpk), STAT = dom_vvl_alloc)
      CALL mpp_sum('domvvl', dom_vvl_alloc)
      IF (dom_vvl_alloc /= 0) CALL ctl_stop('STOP', 'dom_vvl_alloc: failed to allocate arrays')
      un_td = 0._wp
      vn_td = 0._wp
    END IF
    IF (ln_vvl_ztilde) THEN
      ALLOCATE(frq_rst_e3t(jpi, jpj), frq_rst_hdv(jpi, jpj), hdiv_lf(jpi, jpj, jpk), STAT = dom_vvl_alloc)
      CALL mpp_sum('domvvl', dom_vvl_alloc)
      IF (dom_vvl_alloc /= 0) CALL ctl_stop('STOP', 'dom_vvl_alloc: failed to allocate arrays')
    END IF
  END FUNCTION dom_vvl_alloc
  SUBROUTINE dom_vvl_init
    INTEGER :: ji, jj, jk
    INTEGER :: ii0, ii1, ij0, ij1
    REAL(KIND = wp) :: zcoef
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) 'dom_vvl_init : Variable volume activated'
    IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~'
    CALL dom_vvl_ctl
    IF (dom_vvl_alloc() /= 0) CALL ctl_stop('STOP', 'dom_vvl_init : unable to allocate arrays')
    CALL dom_vvl_rst(nit000, 'READ')
    !$ACC KERNELS
    e3t_a(:, :, jpk) = e3t_0(:, :, jpk)
    !$ACC END KERNELS
    CALL dom_vvl_interpol(e3t_b(:, :, :), e3u_b(:, :, :), 'U')
    CALL dom_vvl_interpol(e3t_n(:, :, :), e3u_n(:, :, :), 'U')
    CALL dom_vvl_interpol(e3t_b(:, :, :), e3v_b(:, :, :), 'V')
    CALL dom_vvl_interpol(e3t_n(:, :, :), e3v_n(:, :, :), 'V')
    CALL dom_vvl_interpol(e3u_n(:, :, :), e3f_n(:, :, :), 'F')
    CALL dom_vvl_interpol(e3t_n(:, :, :), e3w_n(:, :, :), 'W')
    CALL dom_vvl_interpol(e3t_b(:, :, :), e3w_b(:, :, :), 'W')
    CALL dom_vvl_interpol(e3u_n(:, :, :), e3uw_n(:, :, :), 'UW')
    CALL dom_vvl_interpol(e3u_b(:, :, :), e3uw_b(:, :, :), 'UW')
    CALL dom_vvl_interpol(e3v_n(:, :, :), e3vw_n(:, :, :), 'VW')
    CALL dom_vvl_interpol(e3v_b(:, :, :), e3vw_b(:, :, :), 'VW')
    !$ACC KERNELS
    e3t_a(:, :, :) = e3t_n(:, :, :)
    e3u_a(:, :, :) = e3u_n(:, :, :)
    e3v_a(:, :, :) = e3v_n(:, :, :)
    gdept_n(:, :, 1) = 0.5_wp * e3w_n(:, :, 1)
    gdepw_n(:, :, 1) = 0.0_wp
    gde3w_n(:, :, 1) = gdept_n(:, :, 1) - sshn(:, :)
    gdept_b(:, :, 1) = 0.5_wp * e3w_b(:, :, 1)
    gdepw_b(:, :, 1) = 0.0_wp
    DO jk = 2, jpk
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zcoef = (tmask(ji, jj, jk) - wmask(ji, jj, jk))
          gdepw_n(ji, jj, jk) = gdepw_n(ji, jj, jk - 1) + e3t_n(ji, jj, jk - 1)
          gdept_n(ji, jj, jk) = zcoef * (gdepw_n(ji, jj, jk) + 0.5 * e3w_n(ji, jj, jk)) + (1 - zcoef) * (gdept_n(ji, jj, jk - 1) + e3w_n(ji, jj, jk))
          gde3w_n(ji, jj, jk) = gdept_n(ji, jj, jk) - sshn(ji, jj)
          gdepw_b(ji, jj, jk) = gdepw_b(ji, jj, jk - 1) + e3t_b(ji, jj, jk - 1)
          gdept_b(ji, jj, jk) = zcoef * (gdepw_b(ji, jj, jk) + 0.5 * e3w_b(ji, jj, jk)) + (1 - zcoef) * (gdept_b(ji, jj, jk - 1) + e3w_b(ji, jj, jk))
        END DO
      END DO
    END DO
    ht_n(:, :) = e3t_n(:, :, 1) * tmask(:, :, 1)
    hu_b(:, :) = e3u_b(:, :, 1) * umask(:, :, 1)
    hu_n(:, :) = e3u_n(:, :, 1) * umask(:, :, 1)
    hv_b(:, :) = e3v_b(:, :, 1) * vmask(:, :, 1)
    hv_n(:, :) = e3v_n(:, :, 1) * vmask(:, :, 1)
    DO jk = 2, jpkm1
      ht_n(:, :) = ht_n(:, :) + e3t_n(:, :, jk) * tmask(:, :, jk)
      hu_b(:, :) = hu_b(:, :) + e3u_b(:, :, jk) * umask(:, :, jk)
      hu_n(:, :) = hu_n(:, :) + e3u_n(:, :, jk) * umask(:, :, jk)
      hv_b(:, :) = hv_b(:, :) + e3v_b(:, :, jk) * vmask(:, :, jk)
      hv_n(:, :) = hv_n(:, :) + e3v_n(:, :, jk) * vmask(:, :, jk)
    END DO
    r1_hu_b(:, :) = ssumask(:, :) / (hu_b(:, :) + 1._wp - ssumask(:, :))
    r1_hu_n(:, :) = ssumask(:, :) / (hu_n(:, :) + 1._wp - ssumask(:, :))
    r1_hv_b(:, :) = ssvmask(:, :) / (hv_b(:, :) + 1._wp - ssvmask(:, :))
    r1_hv_n(:, :) = ssvmask(:, :) / (hv_n(:, :) + 1._wp - ssvmask(:, :))
    !$ACC END KERNELS
    IF (ln_vvl_ztilde) THEN
      !$ACC KERNELS
      frq_rst_e3t(:, :) = 2._wp * rpi / (MAX(rn_rst_e3t, rsmall) * 86400.0_wp)
      frq_rst_hdv(:, :) = 2._wp * rpi / (MAX(rn_lf_cutoff, rsmall) * 86400.0_wp)
      !$ACC END KERNELS
      IF (ln_vvl_ztilde_as_zstar) THEN
        !$ACC KERNELS
        frq_rst_e3t(:, :) = 0._wp
        frq_rst_hdv(:, :) = 1._wp / rdt
        !$ACC END KERNELS
      END IF
      IF (ln_vvl_zstar_at_eqtor) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (ABS(gphit(ji, jj)) >= 6.) THEN
              frq_rst_e3t(ji, jj) = 2.0_wp * rpi / (MAX(rn_rst_e3t, rsmall) * 86400.E0_wp)
              frq_rst_hdv(ji, jj) = 2.0_wp * rpi / (MAX(rn_lf_cutoff, rsmall) * 86400.E0_wp)
            ELSE IF (ABS(gphit(ji, jj)) <= 2.5) THEN
              frq_rst_e3t(ji, jj) = 0.0_wp
              frq_rst_hdv(ji, jj) = 1.0_wp / rdt
            ELSE
              frq_rst_e3t(ji, jj) = 0.0_wp + (frq_rst_e3t(ji, jj) - 0.0_wp) * 0.5_wp * (1.0_wp - COS(rad * (ABS(gphit(ji, jj)) - 2.5_wp) * 180._wp / 3.5_wp))
              frq_rst_hdv(ji, jj) = (1.0_wp / rdt) + (frq_rst_hdv(ji, jj) - (1.E0_wp / rdt)) * 0.5_wp * (1._wp - COS(rad * (ABS(gphit(ji, jj)) - 2.5_wp) * 180._wp / 3.5_wp))
            END IF
          END DO
        END DO
        IF (cn_cfg == "orca" .OR. cn_cfg == "ORCA") THEN
          IF (nn_cfg == 3) THEN
            ii0 = 103
            ii1 = 111
            ij0 = 128
            ij1 = 135
            frq_rst_e3t(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1)) = 0.0_wp
            frq_rst_hdv(mi0(ii0) : mi1(ii1), mj0(ij0) : mj1(ij1)) = 1.E0_wp / rdt
          END IF
        END IF
        !$ACC END KERNELS
      END IF
    END IF
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('e3t_b')
      CALL iom_set_rstw_var_active('e3t_n')
      IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
        CALL iom_set_rstw_var_active('tilde_e3t_b')
        CALL iom_set_rstw_var_active('tilde_e3t_n')
      END IF
      IF (ln_vvl_ztilde) THEN
        CALL iom_set_rstw_var_active('hdiv_lf')
      END IF
    END IF
  END SUBROUTINE dom_vvl_init
  SUBROUTINE dom_vvl_sf_nxt(kt, kcall)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN), OPTIONAL :: kcall
    INTEGER :: ji, jj, jk
    INTEGER, DIMENSION(3) :: ijk_max, ijk_min
    REAL(KIND = wp) :: z2dt, z_tmin, z_tmax
    LOGICAL :: ll_do_bclinic
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zht, z_scale, zwu, zwv, zhdiv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ze3t
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    IF (ln_linssh) RETURN
    CALL profile_psy_data0 % PreStart('dom_vvl_sf_nxt', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dom_vvl_sf_nxt')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dom_vvl_sf_nxt : compute after scale factors'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~'
    END IF
    ll_do_bclinic = .TRUE.
    IF (PRESENT(kcall)) THEN
      IF (kcall == 2 .AND. ln_vvl_ztilde) ll_do_bclinic = .FALSE.
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    z_scale(:, :) = (ssha(:, :) - sshb(:, :)) * ssmask(:, :) / (ht_0(:, :) + sshn(:, :) + 1. - ssmask(:, :))
    DO jk = 1, jpkm1
      e3t_a(:, :, jk) = e3t_b(:, :, jk) + e3t_n(:, :, jk) * z_scale(:, :) * tmask(:, :, jk)
    END DO
    !$ACC END KERNELS
    IF (ln_vvl_ztilde .OR. ln_vvl_layer .AND. ll_do_bclinic) THEN
      !$ACC KERNELS
      zhdiv(:, :) = 0._wp
      zht(:, :) = 0._wp
      DO jk = 1, jpkm1
        zhdiv(:, :) = zhdiv(:, :) + e3t_n(:, :, jk) * hdivn(:, :, jk)
        zht(:, :) = zht(:, :) + e3t_n(:, :, jk) * tmask(:, :, jk)
      END DO
      zhdiv(:, :) = zhdiv(:, :) / (zht(:, :) + 1. - tmask_i(:, :))
      !$ACC END KERNELS
      IF (ln_vvl_ztilde) THEN
        !$ACC KERNELS
        IF (kt > nit000) THEN
          DO jk = 1, jpkm1
            hdiv_lf(:, :, jk) = hdiv_lf(:, :, jk) - rdt * frq_rst_hdv(:, :) * (hdiv_lf(:, :, jk) - e3t_n(:, :, jk) * (hdivn(:, :, jk) - zhdiv(:, :)))
          END DO
        END IF
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      tilde_e3t_a(:, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_vvl_ztilde) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          tilde_e3t_a(:, :, jk) = tilde_e3t_a(:, :, jk) - (e3t_n(:, :, jk) * (hdivn(:, :, jk) - zhdiv(:, :)) - hdiv_lf(:, :, jk))
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        DO jk = 1, jpkm1
          tilde_e3t_a(:, :, jk) = tilde_e3t_a(:, :, jk) - e3t_n(:, :, jk) * (hdivn(:, :, jk) - zhdiv(:, :)) * tmask(:, :, jk)
        END DO
        !$ACC END KERNELS
      END IF
      IF (ln_vvl_ztilde) THEN
        !$ACC KERNELS
        DO jk = 1, jpk
          tilde_e3t_a(:, :, jk) = tilde_e3t_a(:, :, jk) - frq_rst_e3t(:, :) * tilde_e3t_b(:, :, jk)
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      zwu(:, :) = 0._wp
      zwv(:, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            un_td(ji, jj, jk) = rn_ahe3 * umask(ji, jj, jk) * e2_e1u(ji, jj) * (tilde_e3t_b(ji, jj, jk) - tilde_e3t_b(ji + 1, jj, jk))
            vn_td(ji, jj, jk) = rn_ahe3 * vmask(ji, jj, jk) * e1_e2v(ji, jj) * (tilde_e3t_b(ji, jj, jk) - tilde_e3t_b(ji, jj + 1, jk))
            zwu(ji, jj) = zwu(ji, jj) + un_td(ji, jj, jk)
            zwv(ji, jj) = zwv(ji, jj) + vn_td(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          un_td(ji, jj, mbku(ji, jj)) = un_td(ji, jj, mbku(ji, jj)) - zwu(ji, jj)
          vn_td(ji, jj, mbkv(ji, jj)) = vn_td(ji, jj, mbkv(ji, jj)) - zwv(ji, jj)
        END DO
      END DO
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            tilde_e3t_a(ji, jj, jk) = tilde_e3t_a(ji, jj, jk) + (un_td(ji - 1, jj, jk) - un_td(ji, jj, jk) + vn_td(ji, jj - 1, jk) - vn_td(ji, jj, jk)) * r1_e1e2t(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dom_vvl_sf_nxt', 'r1', 0, 0)
      CALL lbc_lnk_multi('domvvl', un_td, 'U', - 1._wp, vn_td, 'V', - 1._wp)
      IF (neuler == 0 .AND. kt == nit000) THEN
        z2dt = rdt
      ELSE
        z2dt = 2.0_wp * rdt
      END IF
      CALL lbc_lnk('domvvl', tilde_e3t_a(:, :, :), 'T', 1._wp)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      tilde_e3t_a(:, :, :) = tilde_e3t_b(:, :, :) + z2dt * tmask(:, :, :) * tilde_e3t_a(:, :, :)
      ze3t(:, :, jpk) = 0._wp
      DO jk = 1, jpkm1
        ze3t(:, :, jk) = tilde_e3t_a(:, :, jk) / e3t_0(:, :, jk) * tmask(:, :, jk) * tmask_i(:, :)
      END DO
      z_tmax = MAXVAL(ze3t(:, :, :))
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('dom_vvl_sf_nxt', 'r2', 0, 0)
      CALL mpp_max('domvvl', z_tmax)
      z_tmin = MINVAL(ze3t(:, :, :))
      CALL mpp_min('domvvl', z_tmin)
      IF ((z_tmax > rn_zdef_max) .OR. (z_tmin < - rn_zdef_max)) THEN
        IF (lk_mpp) THEN
          CALL mpp_maxloc('domvvl', ze3t, tmask, z_tmax, ijk_max)
          CALL mpp_minloc('domvvl', ze3t, tmask, z_tmin, ijk_min)
        ELSE
          ijk_max = MAXLOC(ze3t(:, :, :))
          ijk_max(1) = ijk_max(1) + nimpp - 1
          ijk_max(2) = ijk_max(2) + njmpp - 1
          ijk_min = MINLOC(ze3t(:, :, :))
          ijk_min(1) = ijk_min(1) + nimpp - 1
          ijk_min(2) = ijk_min(2) + njmpp - 1
        END IF
        IF (lwp) THEN
          WRITE(numout, *) 'MAX( tilde_e3t_a(:,:,:) / e3t_0(:,:,:) ) =', z_tmax
          WRITE(numout, *) 'at i, j, k=', ijk_max
          WRITE(numout, *) 'MIN( tilde_e3t_a(:,:,:) / e3t_0(:,:,:) ) =', z_tmin
          WRITE(numout, *) 'at i, j, k=', ijk_min
          CALL ctl_stop('STOP', 'MAX( ABS( tilde_e3t_a(:,:,: ) ) / e3t_0(:,:,:) ) too high')
        END IF
      END IF
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      tilde_e3t_a(:, :, :) = MIN(tilde_e3t_a(:, :, :), rn_zdef_max * e3t_0(:, :, :))
      tilde_e3t_a(:, :, :) = MAX(tilde_e3t_a(:, :, :), - rn_zdef_max * e3t_0(:, :, :))
      DO jk = 1, jpkm1
        dtilde_e3t_a(:, :, jk) = tilde_e3t_a(:, :, jk) - tilde_e3t_b(:, :, jk)
      END DO
      zht(:, :) = 0.
      DO jk = 1, jpkm1
        zht(:, :) = zht(:, :) + tilde_e3t_a(:, :, jk) * tmask(:, :, jk)
      END DO
      z_scale(:, :) = - zht(:, :) / (ht_0(:, :) + sshn(:, :) + 1. - ssmask(:, :))
      DO jk = 1, jpkm1
        dtilde_e3t_a(:, :, jk) = dtilde_e3t_a(:, :, jk) + e3t_n(:, :, jk) * z_scale(:, :) * tmask(:, :, jk)
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      DO jk = 1, jpkm1
        e3t_a(:, :, jk) = e3t_a(:, :, jk) + dtilde_e3t_a(:, :, jk) * tmask(:, :, jk)
      END DO
    END IF
    !$ACC END KERNELS
    IF (ln_vvl_dbg .AND. .NOT. ll_do_bclinic) THEN
      CALL profile_psy_data3 % PreStart('dom_vvl_sf_nxt', 'r3', 0, 0)
      IF (lwp) WRITE(numout, *) 'kt =', kt
      IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
        z_tmax = MAXVAL(tmask(:, :, 1) * tmask_i(:, :) * ABS(zht(:, :)))
        CALL mpp_max('domvvl', z_tmax)
        IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(SUM(tilde_e3t_a))) =', z_tmax
      END IF
      CALL profile_psy_data3 % PostEnd
      !$ACC KERNELS
      zht(:, :) = 0.0_wp
      DO jk = 1, jpkm1
        zht(:, :) = zht(:, :) + e3t_n(:, :, jk) * tmask(:, :, jk)
      END DO
      z_tmax = MAXVAL(tmask(:, :, 1) * tmask_i(:, :) * ABS(ht_0(:, :) + sshn(:, :) - zht(:, :)))
      !$ACC END KERNELS
      CALL profile_psy_data4 % PreStart('dom_vvl_sf_nxt', 'r4', 0, 0)
      CALL mpp_max('domvvl', z_tmax)
      IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(ht_0+sshn-SUM(e3t_n))) =', z_tmax
      CALL profile_psy_data4 % PostEnd
      !$ACC KERNELS
      zht(:, :) = 0.0_wp
      DO jk = 1, jpkm1
        zht(:, :) = zht(:, :) + e3t_a(:, :, jk) * tmask(:, :, jk)
      END DO
      z_tmax = MAXVAL(tmask(:, :, 1) * tmask_i(:, :) * ABS(ht_0(:, :) + ssha(:, :) - zht(:, :)))
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('dom_vvl_sf_nxt', 'r5', 0, 0)
      CALL mpp_max('domvvl', z_tmax)
      IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(ht_0+ssha-SUM(e3t_a))) =', z_tmax
      CALL profile_psy_data5 % PostEnd
      !$ACC KERNELS
      zht(:, :) = 0.0_wp
      DO jk = 1, jpkm1
        zht(:, :) = zht(:, :) + e3t_b(:, :, jk) * tmask(:, :, jk)
      END DO
      z_tmax = MAXVAL(tmask(:, :, 1) * tmask_i(:, :) * ABS(ht_0(:, :) + sshb(:, :) - zht(:, :)))
      !$ACC END KERNELS
      CALL profile_psy_data6 % PreStart('dom_vvl_sf_nxt', 'r6', 0, 0)
      CALL mpp_max('domvvl', z_tmax)
      IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(ht_0+sshb-SUM(e3t_b))) =', z_tmax
      z_tmax = MAXVAL(tmask(:, :, 1) * ABS(sshb(:, :)))
      CALL mpp_max('domvvl', z_tmax)
      IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(sshb))) =', z_tmax
      z_tmax = MAXVAL(tmask(:, :, 1) * ABS(sshn(:, :)))
      CALL mpp_max('domvvl', z_tmax)
      IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(sshn))) =', z_tmax
      z_tmax = MAXVAL(tmask(:, :, 1) * ABS(ssha(:, :)))
      CALL mpp_max('domvvl', z_tmax)
      IF (lwp) WRITE(numout, *) kt, ' MAXVAL(abs(ssha))) =', z_tmax
      CALL profile_psy_data6 % PostEnd
    END IF
    CALL profile_psy_data7 % PreStart('dom_vvl_sf_nxt', 'r7', 0, 0)
    CALL dom_vvl_interpol(e3t_a(:, :, :), e3u_a(:, :, :), 'U')
    CALL dom_vvl_interpol(e3t_a(:, :, :), e3v_a(:, :, :), 'V')
    CALL profile_psy_data7 % PostEnd
    !$ACC KERNELS
    hu_a(:, :) = e3u_a(:, :, 1) * umask(:, :, 1)
    hv_a(:, :) = e3v_a(:, :, 1) * vmask(:, :, 1)
    DO jk = 2, jpkm1
      hu_a(:, :) = hu_a(:, :) + e3u_a(:, :, jk) * umask(:, :, jk)
      hv_a(:, :) = hv_a(:, :) + e3v_a(:, :, jk) * vmask(:, :, jk)
    END DO
    r1_hu_a(:, :) = ssumask(:, :) / (hu_a(:, :) + 1._wp - ssumask(:, :))
    r1_hv_a(:, :) = ssvmask(:, :) / (hv_a(:, :) + 1._wp - ssvmask(:, :))
    !$ACC END KERNELS
    CALL profile_psy_data8 % PreStart('dom_vvl_sf_nxt', 'r8', 0, 0)
    IF (ln_timing) CALL timing_stop('dom_vvl_sf_nxt')
    CALL profile_psy_data8 % PostEnd
  END SUBROUTINE dom_vvl_sf_nxt
  SUBROUTINE dom_vvl_sf_swp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoef
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_linssh) RETURN
    CALL profile_psy_data0 % PreStart('dom_vvl_sf_swp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dom_vvl_sf_swp')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dom_vvl_sf_swp : - time filter and swap of scale factors'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~   - interpolate scale factors and compute depths for next time step'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      IF (neuler == 0 .AND. kt == nit000) THEN
        tilde_e3t_b(:, :, :) = tilde_e3t_n(:, :, :)
      ELSE
        tilde_e3t_b(:, :, :) = tilde_e3t_n(:, :, :) + atfp * (tilde_e3t_b(:, :, :) - 2.0_wp * tilde_e3t_n(:, :, :) + tilde_e3t_a(:, :, :))
      END IF
      tilde_e3t_n(:, :, :) = tilde_e3t_a(:, :, :)
    END IF
    gdept_b(:, :, :) = gdept_n(:, :, :)
    gdepw_b(:, :, :) = gdepw_n(:, :, :)
    e3t_n(:, :, :) = e3t_a(:, :, :)
    e3u_n(:, :, :) = e3u_a(:, :, :)
    e3v_n(:, :, :) = e3v_a(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dom_vvl_sf_swp', 'r1', 0, 0)
    CALL dom_vvl_interpol(e3u_n(:, :, :), e3f_n(:, :, :), 'F')
    CALL dom_vvl_interpol(e3t_n(:, :, :), e3w_n(:, :, :), 'W')
    CALL dom_vvl_interpol(e3u_n(:, :, :), e3uw_n(:, :, :), 'UW')
    CALL dom_vvl_interpol(e3v_n(:, :, :), e3vw_n(:, :, :), 'VW')
    CALL dom_vvl_interpol(e3t_b(:, :, :), e3w_b(:, :, :), 'W')
    CALL dom_vvl_interpol(e3u_b(:, :, :), e3uw_b(:, :, :), 'UW')
    CALL dom_vvl_interpol(e3v_b(:, :, :), e3vw_b(:, :, :), 'VW')
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    gdept_n(:, :, 1) = 0.5_wp * e3w_n(:, :, 1)
    gdepw_n(:, :, 1) = 0.0_wp
    gde3w_n(:, :, 1) = gdept_n(:, :, 1) - sshn(:, :)
    DO jk = 2, jpk
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zcoef = (tmask(ji, jj, jk) - wmask(ji, jj, jk))
          gdepw_n(ji, jj, jk) = gdepw_n(ji, jj, jk - 1) + e3t_n(ji, jj, jk - 1)
          gdept_n(ji, jj, jk) = zcoef * (gdepw_n(ji, jj, jk) + 0.5 * e3w_n(ji, jj, jk)) + (1 - zcoef) * (gdept_n(ji, jj, jk - 1) + e3w_n(ji, jj, jk))
          gde3w_n(ji, jj, jk) = gdept_n(ji, jj, jk) - sshn(ji, jj)
        END DO
      END DO
    END DO
    hu_n(:, :) = hu_a(:, :)
    r1_hu_n(:, :) = r1_hu_a(:, :)
    hv_n(:, :) = hv_a(:, :)
    r1_hv_n(:, :) = r1_hv_a(:, :)
    ht_n(:, :) = e3t_n(:, :, 1) * tmask(:, :, 1)
    DO jk = 2, jpkm1
      ht_n(:, :) = ht_n(:, :) + e3t_n(:, :, jk) * tmask(:, :, jk)
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('dom_vvl_sf_swp', 'r2', 0, 0)
    IF (lrst_oce) CALL dom_vvl_rst(kt, 'WRITE')
    IF (ln_timing) CALL timing_stop('dom_vvl_sf_swp')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE dom_vvl_sf_swp
  SUBROUTINE dom_vvl_interpol(pe3_in, pe3_out, pout)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pe3_in
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pe3_out
    CHARACTER(LEN = *), INTENT(IN) :: pout
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zlnwd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('dom_vvl_interpol', 'r0', 0, 0)
    IF (ln_wd_il) THEN
      zlnwd = 1.0_wp
    ELSE
      zlnwd = 0.0_wp
    END IF
    CALL profile_psy_data0 % PostEnd
    SELECT CASE (pout)
    CASE ('U')
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pe3_out(ji, jj, jk) = 0.5_wp * (umask(ji, jj, jk) * (1.0_wp - zlnwd) + zlnwd) * r1_e1e2u(ji, jj) * (e1e2t(ji, jj) * (pe3_in(ji, jj, jk) - e3t_0(ji, jj, jk)) + e1e2t(ji + 1, jj) * (pe3_in(ji + 1, jj, jk) - e3t_0(ji + 1, jj, jk)))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dom_vvl_interpol', 'r1', 0, 0)
      CALL lbc_lnk('domvvl', pe3_out(:, :, :), 'U', 1._wp)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      pe3_out(:, :, :) = pe3_out(:, :, :) + e3u_0(:, :, :)
      !$ACC END KERNELS
    CASE ('V')
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pe3_out(ji, jj, jk) = 0.5_wp * (vmask(ji, jj, jk) * (1.0_wp - zlnwd) + zlnwd) * r1_e1e2v(ji, jj) * (e1e2t(ji, jj) * (pe3_in(ji, jj, jk) - e3t_0(ji, jj, jk)) + e1e2t(ji, jj + 1) * (pe3_in(ji, jj + 1, jk) - e3t_0(ji, jj + 1, jk)))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('dom_vvl_interpol', 'r2', 0, 0)
      CALL lbc_lnk('domvvl', pe3_out(:, :, :), 'V', 1._wp)
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      pe3_out(:, :, :) = pe3_out(:, :, :) + e3v_0(:, :, :)
      !$ACC END KERNELS
    CASE ('F')
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            pe3_out(ji, jj, jk) = 0.5_wp * (umask(ji, jj, jk) * umask(ji, jj + 1, jk) * (1.0_wp - zlnwd) + zlnwd) * r1_e1e2f(ji, jj) * (e1e2u(ji, jj) * (pe3_in(ji, jj, jk) - e3u_0(ji, jj, jk)) + e1e2u(ji, jj + 1) * (pe3_in(ji, jj + 1, jk) - e3u_0(ji, jj + 1, jk)))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('dom_vvl_interpol', 'r3', 0, 0)
      CALL lbc_lnk('domvvl', pe3_out(:, :, :), 'F', 1._wp)
      CALL profile_psy_data3 % PostEnd
      !$ACC KERNELS
      pe3_out(:, :, :) = pe3_out(:, :, :) + e3f_0(:, :, :)
      !$ACC END KERNELS
    CASE ('W')
      !$ACC KERNELS
      pe3_out(:, :, 1) = e3w_0(:, :, 1) + pe3_in(:, :, 1) - e3t_0(:, :, 1)
      DO jk = 2, jpk
        pe3_out(:, :, jk) = e3w_0(:, :, jk) + (1.0_wp - 0.5_wp * (tmask(:, :, jk) * (1.0_wp - zlnwd) + zlnwd)) * (pe3_in(:, :, jk - 1) - e3t_0(:, :, jk - 1)) + 0.5_wp * (tmask(:, :, jk) * (1.0_wp - zlnwd) + zlnwd) * (pe3_in(:, :, jk) - e3t_0(:, :, jk))
      END DO
      !$ACC END KERNELS
    CASE ('UW')
      !$ACC KERNELS
      pe3_out(:, :, 1) = e3uw_0(:, :, 1) + pe3_in(:, :, 1) - e3u_0(:, :, 1)
      DO jk = 2, jpk
        pe3_out(:, :, jk) = e3uw_0(:, :, jk) + (1.0_wp - 0.5_wp * (umask(:, :, jk) * (1.0_wp - zlnwd) + zlnwd)) * (pe3_in(:, :, jk - 1) - e3u_0(:, :, jk - 1)) + 0.5_wp * (umask(:, :, jk) * (1.0_wp - zlnwd) + zlnwd) * (pe3_in(:, :, jk) - e3u_0(:, :, jk))
      END DO
      !$ACC END KERNELS
    CASE ('VW')
      !$ACC KERNELS
      pe3_out(:, :, 1) = e3vw_0(:, :, 1) + pe3_in(:, :, 1) - e3v_0(:, :, 1)
      DO jk = 2, jpk
        pe3_out(:, :, jk) = e3vw_0(:, :, jk) + (1.0_wp - 0.5_wp * (vmask(:, :, jk) * (1.0_wp - zlnwd) + zlnwd)) * (pe3_in(:, :, jk - 1) - e3v_0(:, :, jk - 1)) + 0.5_wp * (vmask(:, :, jk) * (1.0_wp - zlnwd) + zlnwd) * (pe3_in(:, :, jk) - e3v_0(:, :, jk))
      END DO
      !$ACC END KERNELS
    END SELECT
  END SUBROUTINE dom_vvl_interpol
  SUBROUTINE dom_vvl_rst(kt, cdrw)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER :: ji, jj, jk
    INTEGER :: id1, id2, id3, id4, id5
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        CALL rst_read_open
        CALL iom_get(numror, jpdom_autoglo, 'sshn', sshn, ldxios = lrxios)
        id1 = iom_varid(numror, 'e3t_b', ldstop = .FALSE.)
        id2 = iom_varid(numror, 'e3t_n', ldstop = .FALSE.)
        id3 = iom_varid(numror, 'tilde_e3t_b', ldstop = .FALSE.)
        id4 = iom_varid(numror, 'tilde_e3t_n', ldstop = .FALSE.)
        id5 = iom_varid(numror, 'hdiv_lf', ldstop = .FALSE.)
        IF (MIN(id1, id2) > 0) THEN
          CALL iom_get(numror, jpdom_autoglo, 'e3t_b', e3t_b(:, :, :), ldxios = lrxios)
          CALL iom_get(numror, jpdom_autoglo, 'e3t_n', e3t_n(:, :, :), ldxios = lrxios)
          IF (lwp) WRITE(numout, *) 'dom_vvl_rst : e3t_b and e3t_n found in restart files'
          !$ACC KERNELS
          WHERE (tmask(:, :, :) == 0.0_wp)
            e3t_n(:, :, :) = e3t_0(:, :, :)
            e3t_b(:, :, :) = e3t_0(:, :, :)
          END WHERE
          IF (neuler == 0) THEN
            e3t_b(:, :, :) = e3t_n(:, :, :)
          END IF
          !$ACC END KERNELS
        ELSE IF (id1 > 0) THEN
          IF (lwp) WRITE(numout, *) 'dom_vvl_rst WARNING : e3t_n not found in restart files'
          IF (lwp) WRITE(numout, *) 'e3t_n set equal to e3t_b.'
          IF (lwp) WRITE(numout, *) 'neuler is forced to 0'
          CALL iom_get(numror, jpdom_autoglo, 'e3t_b', e3t_b(:, :, :), ldxios = lrxios)
          !$ACC KERNELS
          e3t_n(:, :, :) = e3t_b(:, :, :)
          neuler = 0
          !$ACC END KERNELS
        ELSE IF (id2 > 0) THEN
          IF (lwp) WRITE(numout, *) 'dom_vvl_rst WARNING : e3t_b not found in restart files'
          IF (lwp) WRITE(numout, *) 'e3t_b set equal to e3t_n.'
          IF (lwp) WRITE(numout, *) 'neuler is forced to 0'
          CALL iom_get(numror, jpdom_autoglo, 'e3t_n', e3t_n(:, :, :), ldxios = lrxios)
          !$ACC KERNELS
          e3t_b(:, :, :) = e3t_n(:, :, :)
          neuler = 0
          !$ACC END KERNELS
        ELSE
          IF (lwp) WRITE(numout, *) 'dom_vvl_rst WARNING : e3t_n not found in restart file'
          IF (lwp) WRITE(numout, *) 'Compute scale factor from sshn'
          IF (lwp) WRITE(numout, *) 'neuler is forced to 0'
          !$ACC KERNELS
          DO jk = 1, jpk
            e3t_n(:, :, jk) = e3t_0(:, :, jk) * (ht_0(:, :) + sshn(:, :)) / (ht_0(:, :) + 1._wp - ssmask(:, :)) * tmask(:, :, jk) + e3t_0(:, :, jk) * (1._wp - tmask(:, :, jk))
          END DO
          e3t_b(:, :, :) = e3t_n(:, :, :)
          neuler = 0
          !$ACC END KERNELS
        END IF
        IF (ln_vvl_zstar) THEN
          IF (MIN(id3, id4) > 0) THEN
            CALL ctl_stop('dom_vvl_rst: z_star cannot restart from a z_tilde or layer run')
          END IF
        ELSE
          IF (MIN(id3, id4) > 0) THEN
            CALL iom_get(numror, jpdom_autoglo, 'tilde_e3t_b', tilde_e3t_b(:, :, :), ldxios = lrxios)
            CALL iom_get(numror, jpdom_autoglo, 'tilde_e3t_n', tilde_e3t_n(:, :, :), ldxios = lrxios)
          ELSE
            !$ACC KERNELS
            tilde_e3t_b(:, :, :) = 0.0_wp
            tilde_e3t_n(:, :, :) = 0.0_wp
            !$ACC END KERNELS
          END IF
          IF (ln_vvl_ztilde) THEN
            IF (id5 > 0) THEN
              CALL iom_get(numror, jpdom_autoglo, 'hdiv_lf', hdiv_lf(:, :, :), ldxios = lrxios)
            ELSE
              !$ACC KERNELS
              hdiv_lf(:, :, :) = 0.0_wp
              !$ACC END KERNELS
            END IF
          END IF
        END IF
      ELSE
        IF (ll_wd) THEN
          IF (cn_cfg == 'wad') THEN
            CALL usr_def_istate(gdept_b, tmask, tsb, ub, vb, sshb)
            !$ACC KERNELS
            tsn(:, :, :, :) = tsb(:, :, :, :)
            sshn(:, :) = sshb(:, :)
            un(:, :, :) = ub(:, :, :)
            vn(:, :, :) = vb(:, :, :)
            !$ACC END KERNELS
          ELSE
            !$ACC KERNELS
            sshn(:, :) = - ssh_ref
            sshb(:, :) = - ssh_ref
            !$ACC loop independent collapse(2)
            DO jj = 1, jpj
              DO ji = 1, jpi
                IF (ht_0(ji, jj) - ssh_ref < rn_wdmin1) THEN
                  sshb(ji, jj) = rn_wdmin1 - (ht_0(ji, jj))
                  sshn(ji, jj) = rn_wdmin1 - (ht_0(ji, jj))
                  ssha(ji, jj) = rn_wdmin1 - (ht_0(ji, jj))
                END IF
              END DO
            END DO
            !$ACC END KERNELS
          END IF
          !$ACC KERNELS
          DO jk = 1, jpk
            e3t_n(:, :, jk) = e3t_0(:, :, jk) * (ht_0(:, :) + sshn(:, :)) / (ht_0(:, :) + 1._wp - ssmask(:, :)) * tmask(:, :, jk) + e3t_0(:, :, jk) * (1._wp - tmask(:, :, jk))
          END DO
          e3t_b(:, :, :) = e3t_n(:, :, :)
          !$ACC END KERNELS
          DO ji = 1, jpi
            DO jj = 1, jpj
              IF (ht_0(ji, jj) .LE. 0.0 .AND. NINT(ssmask(ji, jj)) .EQ. 1) THEN
                CALL ctl_stop('dom_vvl_rst: ht_0 must be positive at potentially wet points')
              END IF
            END DO
          END DO
        ELSE
          !$ACC KERNELS
          sshn(:, :) = 0._wp
          e3t_n(:, :, :) = e3t_0(:, :, :)
          e3t_b(:, :, :) = e3t_0(:, :, :)
          !$ACC END KERNELS
        END IF
        IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
          !$ACC KERNELS
          tilde_e3t_b(:, :, :) = 0._wp
          tilde_e3t_n(:, :, :) = 0._wp
          !$ACC END KERNELS
          IF (ln_vvl_ztilde) hdiv_lf(:, :, :) = 0._wp
        END IF
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, *) '---- dom_vvl_rst ----'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'e3t_b', e3t_b(:, :, :), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'e3t_n', e3t_n(:, :, :), ldxios = lwxios)
      IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
        CALL iom_rstput(kt, nitrst, numrow, 'tilde_e3t_b', tilde_e3t_b(:, :, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'tilde_e3t_n', tilde_e3t_n(:, :, :), ldxios = lwxios)
      END IF
      IF (ln_vvl_ztilde) THEN
        CALL iom_rstput(kt, nitrst, numrow, 'hdiv_lf', hdiv_lf(:, :, :), ldxios = lwxios)
      END IF
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
  END SUBROUTINE dom_vvl_rst
  SUBROUTINE dom_vvl_ctl
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ioptio, ios
    NAMELIST /nam_vvl/ ln_vvl_zstar, ln_vvl_ztilde, ln_vvl_layer, ln_vvl_ztilde_as_zstar, ln_vvl_zstar_at_eqtor, rn_ahe3, rn_rst_e3t, rn_lf_cutoff, rn_zdef_max, ln_vvl_dbg
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dom_vvl_ctl', 'r0', 0, 0)
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nam_vvl, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nam_vvl in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nam_vvl, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nam_vvl in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nam_vvl)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dom_vvl_ctl : choice/control of the variable vertical coordinate'
      WRITE(numout, *) '~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist nam_vvl : chose a vertical coordinate'
      WRITE(numout, *) '      zstar                      ln_vvl_zstar   = ', ln_vvl_zstar
      WRITE(numout, *) '      ztilde                     ln_vvl_ztilde  = ', ln_vvl_ztilde
      WRITE(numout, *) '      layer                      ln_vvl_layer   = ', ln_vvl_layer
      WRITE(numout, *) '      ztilde as zstar   ln_vvl_ztilde_as_zstar  = ', ln_vvl_ztilde_as_zstar
      WRITE(numout, *) '      ztilde near the equator    ln_vvl_zstar_at_eqtor  = ', ln_vvl_zstar_at_eqtor
      WRITE(numout, *) '      !'
      WRITE(numout, *) '      thickness diffusion coefficient                      rn_ahe3      = ', rn_ahe3
      WRITE(numout, *) '      maximum e3t deformation fractional change            rn_zdef_max  = ', rn_zdef_max
      IF (ln_vvl_ztilde_as_zstar) THEN
        WRITE(numout, *) '      ztilde running in zstar emulation mode (ln_vvl_ztilde_as_zstar=T) '
        WRITE(numout, *) '         ignoring namelist timescale parameters and using:'
        WRITE(numout, *) '            hard-wired : z-tilde to zstar restoration timescale (days)'
        WRITE(numout, *) '                         rn_rst_e3t     = 0.e0'
        WRITE(numout, *) '            hard-wired : z-tilde cutoff frequency of low-pass filter (days)'
        WRITE(numout, *) '                         rn_lf_cutoff   = 1.0/rdt'
      ELSE
        WRITE(numout, *) '      z-tilde to zstar restoration timescale (days)        rn_rst_e3t   = ', rn_rst_e3t
        WRITE(numout, *) '      z-tilde cutoff frequency of low-pass filter (days)   rn_lf_cutoff = ', rn_lf_cutoff
      END IF
      WRITE(numout, *) '         debug prints flag                                 ln_vvl_dbg   = ', ln_vvl_dbg
    END IF
    ioptio = 0
    IF (ln_vvl_ztilde_as_zstar) ln_vvl_ztilde = .TRUE.
    IF (ln_vvl_zstar) ioptio = ioptio + 1
    IF (ln_vvl_ztilde) ioptio = ioptio + 1
    IF (ln_vvl_layer) ioptio = ioptio + 1
    IF (ioptio /= 1) CALL ctl_stop('Choose ONE vertical coordinate in namelist nam_vvl')
    IF (.NOT. ln_vvl_zstar .AND. ln_isf) CALL ctl_stop('Only vvl_zstar has been tested with ice shelf cavity')
    IF (lwp) THEN
      WRITE(numout, *)
      IF (ln_vvl_zstar) WRITE(numout, *) '      ==>>>   zstar vertical coordinate is used'
      IF (ln_vvl_ztilde) WRITE(numout, *) '      ==>>>   ztilde vertical coordinate is used'
      IF (ln_vvl_layer) WRITE(numout, *) '      ==>>>   layer vertical coordinate is used'
      IF (ln_vvl_ztilde_as_zstar) WRITE(numout, *) '      ==>>>   to emulate a zstar coordinate'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dom_vvl_ctl
END MODULE domvvl