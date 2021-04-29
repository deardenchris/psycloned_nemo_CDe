MODULE dynldf_iso
  USE oce
  USE dom_oce
  USE ldfdyn
  USE ldftra
  USE zdf_oce
  USE ldfslp
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_ldf_iso
  PUBLIC :: dyn_ldf_iso_alloc
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: akzu, akzv
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: zfuw, zdiu, zdju, zdj1u
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: zfvw, zdiv, zdjv, zdj1v
  CONTAINS
  INTEGER FUNCTION dyn_ldf_iso_alloc()
    ALLOCATE(akzu(jpi, jpj, jpk), zfuw(jpi, jpk), zdiu(jpi, jpk), zdju(jpi, jpk), zdj1u(jpi, jpk), akzv(jpi, jpj, jpk), zfvw(jpi, jpk), zdiv(jpi, jpk), zdjv(jpi, jpk), zdj1v(jpi, jpk), STAT = dyn_ldf_iso_alloc)
    IF (dyn_ldf_iso_alloc /= 0) CALL ctl_warn('dyn_ldf_iso_alloc: array allocate failed.')
  END FUNCTION dyn_ldf_iso_alloc
  SUBROUTINE dyn_ldf_iso(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zabe1, zmskt, zmkt, zuav, zuwslpi, zuwslpj
    REAL(KIND = wp) :: zabe2, zmskf, zmkf, zvav, zvwslpi, zvwslpj
    REAL(KIND = wp) :: zcof0, zcof1, zcof2, zcof3, zcof4, zaht_0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ziut, zivf, zdku, zdk1u
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zjuf, zjvt, zdkv, zdk1v
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('dyn_ldf_iso', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn_ldf_iso : iso-neutral laplacian diffusive operator or '
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~   s-coordinate horizontal diffusive operator'
      IF (dyn_ldf_iso_alloc() /= 0) CALL ctl_stop('STOP', 'dyn_ldf_iso: failed to allocate arrays')
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_dynldf_hor .AND. ln_traldf_iso) THEN
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            uslp(ji, jj, jk) = - (gdept_b(ji + 1, jj, jk) - gdept_b(ji, jj, jk)) * r1_e1u(ji, jj) * umask(ji, jj, jk)
            vslp(ji, jj, jk) = - (gdept_b(ji, jj + 1, jk) - gdept_b(ji, jj, jk)) * r1_e2v(ji, jj) * vmask(ji, jj, jk)
            wslpi(ji, jj, jk) = - (gdepw_b(ji + 1, jj, jk) - gdepw_b(ji - 1, jj, jk)) * r1_e1t(ji, jj) * tmask(ji, jj, jk) * 0.5
            wslpj(ji, jj, jk) = - (gdepw_b(ji, jj + 1, jk) - gdepw_b(ji, jj - 1, jk)) * r1_e2t(ji, jj) * tmask(ji, jj, jk) * 0.5
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dyn_ldf_iso', 'r1', 0, 0)
      CALL lbc_lnk_multi('dynldf_iso', uslp, 'U', - 1., vslp, 'V', - 1., wslpi, 'W', - 1., wslpj, 'W', - 1.)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('dyn_ldf_iso', 'r2', 0, 0)
    zaht_0 = 0.5_wp * rn_Ud * rn_Ld
    CALL profile_psy_data2 % PostEnd
    DO jk = 1, jpkm1
      !$ACC KERNELS
      zdk1u(:, :) = (ub(:, :, jk) - ub(:, :, jk + 1)) * umask(:, :, jk + 1)
      zdk1v(:, :) = (vb(:, :, jk) - vb(:, :, jk + 1)) * vmask(:, :, jk + 1)
      IF (jk == 1) THEN
        zdku(:, :) = zdk1u(:, :)
        zdkv(:, :) = zdk1v(:, :)
      ELSE
        zdku(:, :) = (ub(:, :, jk - 1) - ub(:, :, jk)) * umask(:, :, jk)
        zdkv(:, :) = (vb(:, :, jk - 1) - vb(:, :, jk)) * vmask(:, :, jk)
      END IF
      !$ACC END KERNELS
      IF (ln_zps) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpi
            zabe1 = (ahmt(ji, jj, jk) + rn_ahm_b) * e2t(ji, jj) * MIN(e3u_n(ji, jj, jk), e3u_n(ji - 1, jj, jk)) * r1_e1t(ji, jj)
            zmskt = 1._wp / MAX(umask(ji - 1, jj, jk) + umask(ji, jj, jk + 1) + umask(ji - 1, jj, jk + 1) + umask(ji, jj, jk), 1._wp)
            zcof1 = - zaht_0 * e2t(ji, jj) * zmskt * 0.5 * (uslp(ji - 1, jj, jk) + uslp(ji, jj, jk))
            ziut(ji, jj) = (zabe1 * (ub(ji, jj, jk) - ub(ji - 1, jj, jk)) + zcof1 * (zdku(ji, jj) + zdk1u(ji - 1, jj) + zdk1u(ji, jj) + zdku(ji - 1, jj))) * tmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpi
            zabe1 = (ahmt(ji, jj, jk) + rn_ahm_b) * e2t(ji, jj) * e3t_n(ji, jj, jk) * r1_e1t(ji, jj)
            zmskt = 1._wp / MAX(umask(ji - 1, jj, jk) + umask(ji, jj, jk + 1) + umask(ji - 1, jj, jk + 1) + umask(ji, jj, jk), 1._wp)
            zcof1 = - zaht_0 * e2t(ji, jj) * zmskt * 0.5 * (uslp(ji - 1, jj, jk) + uslp(ji, jj, jk))
            ziut(ji, jj) = (zabe1 * (ub(ji, jj, jk) - ub(ji - 1, jj, jk)) + zcof1 * (zdku(ji, jj) + zdk1u(ji - 1, jj) + zdk1u(ji, jj) + zdku(ji - 1, jj))) * tmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zabe2 = (ahmf(ji, jj, jk) + rn_ahm_b) * e1f(ji, jj) * e3f_n(ji, jj, jk) * r1_e2f(ji, jj)
          zmskf = 1._wp / MAX(umask(ji, jj + 1, jk) + umask(ji, jj, jk + 1) + umask(ji, jj + 1, jk + 1) + umask(ji, jj, jk), 1._wp)
          zcof2 = - zaht_0 * e1f(ji, jj) * zmskf * 0.5 * (vslp(ji + 1, jj, jk) + vslp(ji, jj, jk))
          zjuf(ji, jj) = (zabe2 * (ub(ji, jj + 1, jk) - ub(ji, jj, jk)) + zcof2 * (zdku(ji, jj + 1) + zdk1u(ji, jj) + zdk1u(ji, jj + 1) + zdku(ji, jj))) * fmask(ji, jj, jk)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          zabe1 = (ahmf(ji, jj, jk) + rn_ahm_b) * e2f(ji, jj) * e3f_n(ji, jj, jk) * r1_e1f(ji, jj)
          zmskf = 1._wp / MAX(vmask(ji + 1, jj, jk) + vmask(ji, jj, jk + 1) + vmask(ji + 1, jj, jk + 1) + vmask(ji, jj, jk), 1._wp)
          zcof1 = - zaht_0 * e2f(ji, jj) * zmskf * 0.5 * (uslp(ji, jj + 1, jk) + uslp(ji, jj, jk))
          zivf(ji, jj) = (zabe1 * (vb(ji + 1, jj, jk) - vb(ji, jj, jk)) + zcof1 * (zdkv(ji, jj) + zdk1v(ji + 1, jj) + zdk1v(ji, jj) + zdkv(ji + 1, jj))) * fmask(ji, jj, jk)
        END DO
      END DO
      !$ACC END KERNELS
      IF (ln_zps) THEN
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 1, jpim1
            zabe2 = (ahmt(ji, jj, jk) + rn_ahm_b) * e1t(ji, jj) * MIN(e3v_n(ji, jj, jk), e3v_n(ji, jj - 1, jk)) * r1_e2t(ji, jj)
            zmskt = 1._wp / MAX(vmask(ji, jj - 1, jk) + vmask(ji, jj, jk + 1) + vmask(ji, jj - 1, jk + 1) + vmask(ji, jj, jk), 1._wp)
            zcof2 = - zaht_0 * e1t(ji, jj) * zmskt * 0.5 * (vslp(ji, jj - 1, jk) + vslp(ji, jj, jk))
            zjvt(ji, jj) = (zabe2 * (vb(ji, jj, jk) - vb(ji, jj - 1, jk)) + zcof2 * (zdkv(ji, jj - 1) + zdk1v(ji, jj) + zdk1v(ji, jj - 1) + zdkv(ji, jj))) * tmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        !$ACC loop independent collapse(2)
        DO jj = 2, jpj
          DO ji = 1, jpim1
            zabe2 = (ahmt(ji, jj, jk) + rn_ahm_b) * e1t(ji, jj) * e3t_n(ji, jj, jk) * r1_e2t(ji, jj)
            zmskt = 1. / MAX(vmask(ji, jj - 1, jk) + vmask(ji, jj, jk + 1) + vmask(ji, jj - 1, jk + 1) + vmask(ji, jj, jk), 1.)
            zcof2 = - zaht_0 * e1t(ji, jj) * zmskt * 0.5 * (vslp(ji, jj - 1, jk) + vslp(ji, jj, jk))
            zjvt(ji, jj) = (zabe2 * (vb(ji, jj, jk) - vb(ji, jj - 1, jk)) + zcof2 * (zdkv(ji, jj - 1) + zdk1v(ji, jj) + zdk1v(ji, jj - 1) + zdkv(ji, jj))) * tmask(ji, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) + (ziut(ji + 1, jj) - ziut(ji, jj) + zjuf(ji, jj) - zjuf(ji, jj - 1)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) + (zivf(ji, jj) - zivf(ji - 1, jj) + zjvt(ji, jj + 1) - zjvt(ji, jj)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    CALL profile_psy_data3 % PreStart('dyn_ldf_iso', 'r3', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' ldfh - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = vmask, clinfo3 = 'dyn')
    CALL profile_psy_data3 % PostEnd
    DO jj = 2, jpjm1
      !$ACC KERNELS
      DO jk = 1, jpk
        DO ji = 2, jpi
          zdiu(ji, jk) = tmask(ji, jj, jk) * (ub(ji, jj, jk) - ub(ji - 1, jj, jk))
          zdju(ji, jk) = fmask(ji, jj, jk) * (ub(ji, jj + 1, jk) - ub(ji, jj, jk))
          zdjv(ji, jk) = tmask(ji, jj, jk) * (vb(ji, jj, jk) - vb(ji, jj - 1, jk))
          zdj1u(ji, jk) = fmask(ji, jj - 1, jk) * (ub(ji, jj, jk) - ub(ji, jj - 1, jk))
          zdj1v(ji, jk) = tmask(ji, jj + 1, jk) * (vb(ji, jj + 1, jk) - vb(ji, jj, jk))
        END DO
      END DO
      DO jk = 1, jpk
        DO ji = 1, jpim1
          zdiv(ji, jk) = fmask(ji, jj, jk) * (vb(ji + 1, jj, jk) - vb(ji, jj, jk))
        END DO
      END DO
      DO ji = 1, jpi
        zfuw(ji, 1) = 0.E0
        zfvw(ji, 1) = 0.E0
        zfuw(ji, jpk) = 0.E0
        zfvw(ji, jpk) = 0.E0
      END DO
      DO jk = 2, jpkm1
        DO ji = 2, jpim1
          zcof0 = 0.5_wp * zaht_0 * umask(ji, jj, jk)
          zuwslpi = zcof0 * (wslpi(ji + 1, jj, jk) + wslpi(ji, jj, jk))
          zuwslpj = zcof0 * (wslpj(ji + 1, jj, jk) + wslpj(ji, jj, jk))
          zmkt = 1. / MAX(tmask(ji, jj, jk - 1) + tmask(ji + 1, jj, jk - 1) + tmask(ji, jj, jk) + tmask(ji + 1, jj, jk), 1.)
          zmkf = 1. / MAX(fmask(ji, jj - 1, jk - 1) + fmask(ji, jj, jk - 1) + fmask(ji, jj - 1, jk) + fmask(ji, jj, jk), 1.)
          zcof3 = - e2u(ji, jj) * zmkt * zuwslpi
          zcof4 = - e1u(ji, jj) * zmkf * zuwslpj
          zfuw(ji, jk) = zcof3 * (zdiu(ji, jk - 1) + zdiu(ji + 1, jk - 1) + zdiu(ji, jk) + zdiu(ji + 1, jk)) + zcof4 * (zdj1u(ji, jk - 1) + zdju(ji, jk - 1) + zdj1u(ji, jk) + zdju(ji, jk))
          akzu(ji, jj, jk) = (zuwslpi * zuwslpi + zuwslpj * zuwslpj) / zaht_0
        END DO
      END DO
      DO jk = 2, jpkm1
        DO ji = 2, jpim1
          zcof0 = 0.5_wp * zaht_0 * vmask(ji, jj, jk)
          zvwslpi = zcof0 * (wslpi(ji, jj + 1, jk) + wslpi(ji, jj, jk))
          zvwslpj = zcof0 * (wslpj(ji, jj + 1, jk) + wslpj(ji, jj, jk))
          zmkf = 1. / MAX(fmask(ji - 1, jj, jk - 1) + fmask(ji, jj, jk - 1) + fmask(ji - 1, jj, jk) + fmask(ji, jj, jk), 1.)
          zmkt = 1. / MAX(tmask(ji, jj, jk - 1) + tmask(ji, jj + 1, jk - 1) + tmask(ji, jj, jk) + tmask(ji, jj + 1, jk), 1.)
          zcof3 = - e2v(ji, jj) * zmkf * zvwslpi
          zcof4 = - e1v(ji, jj) * zmkt * zvwslpj
          zfvw(ji, jk) = zcof3 * (zdiv(ji, jk - 1) + zdiv(ji - 1, jk - 1) + zdiv(ji, jk) + zdiv(ji - 1, jk)) + zcof4 * (zdjv(ji, jk - 1) + zdj1v(ji, jk - 1) + zdjv(ji, jk) + zdj1v(ji, jk))
          akzv(ji, jj, jk) = (zvwslpi * zvwslpi + zvwslpj * zvwslpj) / zaht_0
        END DO
      END DO
      DO jk = 1, jpkm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) + (zfuw(ji, jk) - zfuw(ji, jk + 1)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) + (zfvw(ji, jk) - zfvw(ji, jk + 1)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE dyn_ldf_iso
END MODULE dynldf_iso