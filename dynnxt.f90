MODULE dynnxt
  USE oce
  USE dom_oce
  USE sbc_oce
  USE sbcrnf
  USE sbcisf
  USE phycst
  USE dynadv
  USE dynspg_ts
  USE domvvl
  USE bdy_oce, ONLY: ln_bdy
  USE bdydta
  USE bdydyn
  USE bdyvol
  USE trd_oce
  USE trddyn
  USE trdken
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_nxt
  CONTAINS
  SUBROUTINE dyn_nxt(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikt
    REAL(KIND = wp) :: zue3a, zue3n, zue3b, zuf, zcoef
    REAL(KIND = wp) :: zve3a, zve3n, zve3b, zvf, z1_2dt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zue, zve
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ze3u_f, ze3v_f, zua, zva
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    CALL profile_psy_data0 % PreStart('dyn_nxt', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dyn_nxt')
    IF (ln_dynspg_ts) ALLOCATE(zue(jpi, jpj), zve(jpi, jpj))
    IF (l_trddyn) ALLOCATE(zua(jpi, jpj, jpk), zva(jpi, jpj, jpk))
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_nxt : time stepping'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_dynspg_ts) THEN
      !$ACC KERNELS
      zue(:, :) = e3u_a(:, :, 1) * ua(:, :, 1) * umask(:, :, 1)
      zve(:, :) = e3v_a(:, :, 1) * va(:, :, 1) * vmask(:, :, 1)
      !$ACC END KERNELS
      DO jk = 2, jpkm1
        !$ACC KERNELS
        zue(:, :) = zue(:, :) + e3u_a(:, :, jk) * ua(:, :, jk) * umask(:, :, jk)
        zve(:, :) = zve(:, :) + e3v_a(:, :, jk) * va(:, :, jk) * vmask(:, :, jk)
        !$ACC END KERNELS
      END DO
      DO jk = 1, jpkm1
        !$ACC KERNELS
        ua(:, :, jk) = (ua(:, :, jk) - zue(:, :) * r1_hu_a(:, :) + ua_b(:, :)) * umask(:, :, jk)
        va(:, :, jk) = (va(:, :, jk) - zve(:, :) * r1_hv_a(:, :) + va_b(:, :)) * vmask(:, :, jk)
        !$ACC END KERNELS
      END DO
      IF (.NOT. ln_bt_fw) THEN
        DO jk = 1, jpkm1
          !$ACC KERNELS
          un(:, :, jk) = (un(:, :, jk) - un_adv(:, :) * r1_hu_n(:, :) + un_b(:, :)) * umask(:, :, jk)
          vn(:, :, jk) = (vn(:, :, jk) - vn_adv(:, :) * r1_hv_n(:, :) + vn_b(:, :)) * vmask(:, :, jk)
          !$ACC END KERNELS
        END DO
      END IF
    END IF
    CALL profile_psy_data1 % PreStart('dyn_nxt', 'r1', 0, 0)
    CALL lbc_lnk_multi('dynnxt', ua, 'U', - 1., va, 'V', - 1.)
    IF (ln_bdy .AND. ln_dynspg_exp) CALL bdy_dyn(kt)
    IF (ln_bdy .AND. ln_dynspg_ts) CALL bdy_dyn(kt, dyn3d_only = .TRUE.)
    CALL profile_psy_data1 % PostEnd
    IF (l_trddyn) THEN
      CALL profile_psy_data2 % PreStart('dyn_nxt', 'r2', 0, 0)
      z1_2dt = 1._wp / (2. * rdt)
      IF (neuler == 0 .AND. kt == nit000) z1_2dt = 1._wp / rdt
      IF (ln_KE_trd) CALL trd_dyn(ua, va, jpdyn_ken, kt)
      CALL profile_psy_data2 % PostEnd
      IF (ln_dyn_trd) THEN
        !$ACC KERNELS
        zua(:, :, :) = (ua(:, :, :) - ub(:, :, :)) * z1_2dt
        zva(:, :, :) = (va(:, :, :) - vb(:, :, :)) * z1_2dt
        !$ACC END KERNELS
        CALL profile_psy_data3 % PreStart('dyn_nxt', 'r3', 0, 0)
        CALL iom_put("utrd_tot", zua)
        CALL iom_put("vtrd_tot", zva)
        CALL profile_psy_data3 % PostEnd
      END IF
      !$ACC KERNELS
      zua(:, :, :) = un(:, :, :)
      zva(:, :, :) = vn(:, :, :)
      !$ACC END KERNELS
    END IF
    IF (neuler == 0 .AND. kt == nit000) THEN
      DO jk = 1, jpkm1
        !$ACC KERNELS
        un(:, :, jk) = ua(:, :, jk)
        vn(:, :, jk) = va(:, :, jk)
        !$ACC END KERNELS
      END DO
      IF (.NOT. ln_linssh) THEN
        DO jk = 1, jpkm1
          !$ACC KERNELS
          e3t_n(:, :, jk) = e3t_a(:, :, jk)
          e3u_n(:, :, jk) = e3u_a(:, :, jk)
          e3v_n(:, :, jk) = e3v_a(:, :, jk)
          !$ACC END KERNELS
        END DO
      END IF
    ELSE
      IF (ln_linssh) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              zuf = un(ji, jj, jk) + atfp * (ub(ji, jj, jk) - 2._wp * un(ji, jj, jk) + ua(ji, jj, jk))
              zvf = vn(ji, jj, jk) + atfp * (vb(ji, jj, jk) - 2._wp * vn(ji, jj, jk) + va(ji, jj, jk))
              ub(ji, jj, jk) = zuf
              vb(ji, jj, jk) = zvf
              un(ji, jj, jk) = ua(ji, jj, jk)
              vn(ji, jj, jk) = va(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        DO jk = 1, jpkm1
          e3t_b(:, :, jk) = e3t_n(:, :, jk) + atfp * (e3t_b(:, :, jk) - 2._wp * e3t_n(:, :, jk) + e3t_a(:, :, jk))
        END DO
        zcoef = atfp * rdt * r1_rau0
        e3t_b(:, :, 1) = e3t_b(:, :, 1) - zcoef * (emp_b(:, :) - emp(:, :)) * tmask(:, :, 1)
        !$ACC END KERNELS
        IF (ln_rnf) THEN
          IF (ln_rnf_depth) THEN
            !$ACC KERNELS
            DO jk = 1, jpkm1
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, jpj
                DO ji = 1, jpi
                  IF (jk <= nk_rnf(ji, jj)) THEN
                    e3t_b(ji, jj, jk) = e3t_b(ji, jj, jk) - zcoef * (- rnf_b(ji, jj) + rnf(ji, jj)) * (e3t_n(ji, jj, jk) / &
&h_rnf(ji, jj)) * tmask(ji, jj, jk)
                  END IF
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          ELSE
            !$ACC KERNELS
            e3t_b(:, :, 1) = e3t_b(:, :, 1) - zcoef * (- rnf_b(:, :) + rnf(:, :)) * tmask(:, :, 1)
            !$ACC END KERNELS
          END IF
        END IF
        IF (ln_isf) THEN
          !$ACC KERNELS
          DO jk = 1, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpj
              DO ji = 1, jpi
                IF (misfkt(ji, jj) <= jk .AND. jk < misfkb(ji, jj)) THEN
                  e3t_b(ji, jj, jk) = e3t_b(ji, jj, jk) - zcoef * (fwfisf_b(ji, jj) - fwfisf(ji, jj)) * (e3t_n(ji, jj, jk) * &
&r1_hisf_tbl(ji, jj)) * tmask(ji, jj, jk)
                ELSE IF (jk == misfkb(ji, jj)) THEN
                  e3t_b(ji, jj, jk) = e3t_b(ji, jj, jk) - zcoef * (fwfisf_b(ji, jj) - fwfisf(ji, jj)) * (e3t_n(ji, jj, jk) * &
&r1_hisf_tbl(ji, jj)) * ralpha(ji, jj) * tmask(ji, jj, jk)
                END IF
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        IF (ln_dynadv_vec) THEN
          CALL profile_psy_data4 % PreStart('dyn_nxt', 'r4', 0, 0)
          CALL dom_vvl_interpol(e3t_b(:, :, :), e3u_b(:, :, :), 'U')
          CALL dom_vvl_interpol(e3t_b(:, :, :), e3v_b(:, :, :), 'V')
          CALL profile_psy_data4 % PostEnd
          !$ACC KERNELS
          DO jk = 1, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpj
              DO ji = 1, jpi
                zuf = un(ji, jj, jk) + atfp * (ub(ji, jj, jk) - 2._wp * un(ji, jj, jk) + ua(ji, jj, jk))
                zvf = vn(ji, jj, jk) + atfp * (vb(ji, jj, jk) - 2._wp * vn(ji, jj, jk) + va(ji, jj, jk))
                ub(ji, jj, jk) = zuf
                vb(ji, jj, jk) = zvf
                un(ji, jj, jk) = ua(ji, jj, jk)
                vn(ji, jj, jk) = va(ji, jj, jk)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        ELSE
          CALL profile_psy_data5 % PreStart('dyn_nxt', 'r5', 0, 0)
          ALLOCATE(ze3u_f(jpi, jpj, jpk), ze3v_f(jpi, jpj, jpk))
          CALL dom_vvl_interpol(e3t_b(:, :, :), ze3u_f, 'U')
          CALL dom_vvl_interpol(e3t_b(:, :, :), ze3v_f, 'V')
          CALL profile_psy_data5 % PostEnd
          !$ACC KERNELS
          DO jk = 1, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpj
              DO ji = 1, jpi
                zue3a = e3u_a(ji, jj, jk) * ua(ji, jj, jk)
                zve3a = e3v_a(ji, jj, jk) * va(ji, jj, jk)
                zue3n = e3u_n(ji, jj, jk) * un(ji, jj, jk)
                zve3n = e3v_n(ji, jj, jk) * vn(ji, jj, jk)
                zue3b = e3u_b(ji, jj, jk) * ub(ji, jj, jk)
                zve3b = e3v_b(ji, jj, jk) * vb(ji, jj, jk)
                zuf = (zue3n + atfp * (zue3b - 2._wp * zue3n + zue3a)) / ze3u_f(ji, jj, jk)
                zvf = (zve3n + atfp * (zve3b - 2._wp * zve3n + zve3a)) / ze3v_f(ji, jj, jk)
                ub(ji, jj, jk) = zuf
                vb(ji, jj, jk) = zvf
                un(ji, jj, jk) = ua(ji, jj, jk)
                vn(ji, jj, jk) = va(ji, jj, jk)
              END DO
            END DO
          END DO
          e3u_b(:, :, 1 : jpkm1) = ze3u_f(:, :, 1 : jpkm1)
          e3v_b(:, :, 1 : jpkm1) = ze3v_f(:, :, 1 : jpkm1)
          !$ACC END KERNELS
          DEALLOCATE(ze3u_f, ze3v_f)
        END IF
      END IF
      IF (ln_dynspg_ts .AND. ln_bt_fw) THEN
        !$ACC KERNELS
        zue(:, :) = e3u_b(:, :, 1) * ub(:, :, 1) * umask(:, :, 1)
        zve(:, :) = e3v_b(:, :, 1) * vb(:, :, 1) * vmask(:, :, 1)
        !$ACC END KERNELS
        DO jk = 2, jpkm1
          !$ACC KERNELS
          zue(:, :) = zue(:, :) + e3u_b(:, :, jk) * ub(:, :, jk) * umask(:, :, jk)
          zve(:, :) = zve(:, :) + e3v_b(:, :, jk) * vb(:, :, jk) * vmask(:, :, jk)
          !$ACC END KERNELS
        END DO
        DO jk = 1, jpkm1
          !$ACC KERNELS
          ub(:, :, jk) = ub(:, :, jk) - (zue(:, :) * r1_hu_n(:, :) - un_b(:, :)) * umask(:, :, jk)
          vb(:, :, jk) = vb(:, :, jk) - (zve(:, :) * r1_hv_n(:, :) - vn_b(:, :)) * vmask(:, :, jk)
          !$ACC END KERNELS
        END DO
      END IF
    END IF
    IF (.NOT. ln_linssh) THEN
      !$ACC KERNELS
      hu_b(:, :) = e3u_b(:, :, 1) * umask(:, :, 1)
      hv_b(:, :) = e3v_b(:, :, 1) * vmask(:, :, 1)
      !$ACC END KERNELS
      DO jk = 2, jpkm1
        !$ACC KERNELS
        hu_b(:, :) = hu_b(:, :) + e3u_b(:, :, jk) * umask(:, :, jk)
        hv_b(:, :) = hv_b(:, :) + e3v_b(:, :, jk) * vmask(:, :, jk)
        !$ACC END KERNELS
      END DO
      !$ACC KERNELS
      r1_hu_b(:, :) = ssumask(:, :) / (hu_b(:, :) + 1._wp - ssumask(:, :))
      r1_hv_b(:, :) = ssvmask(:, :) / (hv_b(:, :) + 1._wp - ssvmask(:, :))
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    un_b(:, :) = e3u_a(:, :, 1) * un(:, :, 1) * umask(:, :, 1)
    ub_b(:, :) = e3u_b(:, :, 1) * ub(:, :, 1) * umask(:, :, 1)
    vn_b(:, :) = e3v_a(:, :, 1) * vn(:, :, 1) * vmask(:, :, 1)
    vb_b(:, :) = e3v_b(:, :, 1) * vb(:, :, 1) * vmask(:, :, 1)
    !$ACC END KERNELS
    DO jk = 2, jpkm1
      !$ACC KERNELS
      un_b(:, :) = un_b(:, :) + e3u_a(:, :, jk) * un(:, :, jk) * umask(:, :, jk)
      ub_b(:, :) = ub_b(:, :) + e3u_b(:, :, jk) * ub(:, :, jk) * umask(:, :, jk)
      vn_b(:, :) = vn_b(:, :) + e3v_a(:, :, jk) * vn(:, :, jk) * vmask(:, :, jk)
      vb_b(:, :) = vb_b(:, :) + e3v_b(:, :, jk) * vb(:, :, jk) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    un_b(:, :) = un_b(:, :) * r1_hu_a(:, :)
    vn_b(:, :) = vn_b(:, :) * r1_hv_a(:, :)
    ub_b(:, :) = ub_b(:, :) * r1_hu_b(:, :)
    vb_b(:, :) = vb_b(:, :) * r1_hv_b(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data6 % PreStart('dyn_nxt', 'r6', 0, 0)
    IF (.NOT. ln_dynspg_ts) THEN
      CALL iom_put("ubar", un_b(:, :))
      CALL iom_put("vbar", vn_b(:, :))
    END IF
    CALL profile_psy_data6 % PostEnd
    IF (l_trddyn) THEN
      !$ACC KERNELS
      zua(:, :, :) = (ub(:, :, :) - zua(:, :, :)) * z1_2dt
      zva(:, :, :) = (vb(:, :, :) - zva(:, :, :)) * z1_2dt
      !$ACC END KERNELS
      CALL trd_dyn(zua, zva, jpdyn_atf, kt)
    END IF
    CALL profile_psy_data7 % PreStart('dyn_nxt', 'r7', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = un, clinfo1 = ' nxt  - Un: ', mask1 = umask, tab3d_2 = vn, clinfo2 = ' Vn: ', mask2 = vmask)
    IF (ln_dynspg_ts) DEALLOCATE(zue, zve)
    IF (l_trddyn) DEALLOCATE(zua, zva)
    IF (ln_timing) CALL timing_stop('dyn_nxt')
    CALL profile_psy_data7 % PostEnd
  END SUBROUTINE dyn_nxt
END MODULE dynnxt