MODULE dynadv_ubs
  USE oce
  USE dom_oce
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE prtctl
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  REAL(KIND = wp), PARAMETER :: gamma1 = 1._wp / 3._wp
  REAL(KIND = wp), PARAMETER :: gamma2 = 1._wp / 32._wp
  PUBLIC :: dyn_adv_ubs
  CONTAINS
  SUBROUTINE dyn_adv_ubs(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zui, zvj, zfuj, zfvi, zl_u, zl_v
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zfu_t, zfu_f, zfu_uw, zfu
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zfv_t, zfv_f, zfv_vw, zfv, zfw
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, 2) :: zlu_uu, zlu_uv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, 2) :: zlv_vv, zlv_vu
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('dyn_adv_ubs', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn_adv_ubs : UBS flux form momentum advection'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zfu_t(:, :, :) = 0._wp
    zfv_t(:, :, :) = 0._wp
    zfu_f(:, :, :) = 0._wp
    zfv_f(:, :, :) = 0._wp
    zlu_uu(:, :, :, :) = 0._wp
    zlv_vv(:, :, :, :) = 0._wp
    zlu_uv(:, :, :, :) = 0._wp
    zlv_vu(:, :, :, :) = 0._wp
    !$ACC END KERNELS
    IF (l_trddyn) THEN
      !$ACC KERNELS
      zfu_uw(:, :, :) = ua(:, :, :)
      zfv_vw(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpkm1
      zfu(:, :, jk) = e2u(:, :) * e3u_n(:, :, jk) * un(:, :, jk)
      zfv(:, :, jk) = e1v(:, :) * e3v_n(:, :, jk) * vn(:, :, jk)
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zlu_uu(ji, jj, jk, 1) = (ub(ji + 1, jj, jk) - 2. * ub(ji, jj, jk) + ub(ji - 1, jj, jk)) * umask(ji, jj, jk)
          zlv_vv(ji, jj, jk, 1) = (vb(ji, jj + 1, jk) - 2. * vb(ji, jj, jk) + vb(ji, jj - 1, jk)) * vmask(ji, jj, jk)
          zlu_uv(ji, jj, jk, 1) = (ub(ji, jj + 1, jk) - ub(ji, jj, jk)) * fmask(ji, jj, jk) - (ub(ji, jj, jk) - ub(ji, jj - 1, jk)) * fmask(ji, jj - 1, jk)
          zlv_vu(ji, jj, jk, 1) = (vb(ji + 1, jj, jk) - vb(ji, jj, jk)) * fmask(ji, jj, jk) - (vb(ji, jj, jk) - vb(ji - 1, jj, jk)) * fmask(ji - 1, jj, jk)
          zlu_uu(ji, jj, jk, 2) = (zfu(ji + 1, jj, jk) - 2. * zfu(ji, jj, jk) + zfu(ji - 1, jj, jk)) * umask(ji, jj, jk)
          zlv_vv(ji, jj, jk, 2) = (zfv(ji, jj + 1, jk) - 2. * zfv(ji, jj, jk) + zfv(ji, jj - 1, jk)) * vmask(ji, jj, jk)
          zlu_uv(ji, jj, jk, 2) = (zfu(ji, jj + 1, jk) - zfu(ji, jj, jk)) * fmask(ji, jj, jk) - (zfu(ji, jj, jk) - zfu(ji, jj - 1, jk)) * fmask(ji, jj - 1, jk)
          zlv_vu(ji, jj, jk, 2) = (zfv(ji + 1, jj, jk) - zfv(ji, jj, jk)) * fmask(ji, jj, jk) - (zfv(ji, jj, jk) - zfv(ji - 1, jj, jk)) * fmask(ji - 1, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dyn_adv_ubs', 'r1', 0, 0)
    CALL lbc_lnk_multi('dynadv_ubs', zlu_uu(:, :, :, 1), 'U', 1., zlu_uv(:, :, :, 1), 'U', 1., zlu_uu(:, :, :, 2), 'U', 1., zlu_uv(:, :, :, 2), 'U', 1., zlv_vv(:, :, :, 1), 'V', 1., zlv_vu(:, :, :, 1), 'V', 1., zlv_vv(:, :, :, 2), 'V', 1., zlv_vu(:, :, :, 2), 'V', 1.)
    CALL profile_psy_data1 % PostEnd
    DO jk = 1, jpkm1
      !$ACC KERNELS
      zfu(:, :, jk) = 0.25_wp * e2u(:, :) * e3u_n(:, :, jk) * un(:, :, jk)
      zfv(:, :, jk) = 0.25_wp * e1v(:, :) * e3v_n(:, :, jk) * vn(:, :, jk)
      !$ACC loop independent collapse(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zui = (un(ji, jj, jk) + un(ji + 1, jj, jk))
          zvj = (vn(ji, jj, jk) + vn(ji, jj + 1, jk))
          IF (zui > 0) THEN
            zl_u = zlu_uu(ji, jj, jk, 1)
          ELSE
            zl_u = zlu_uu(ji + 1, jj, jk, 1)
          END IF
          IF (zvj > 0) THEN
            zl_v = zlv_vv(ji, jj, jk, 1)
          ELSE
            zl_v = zlv_vv(ji, jj + 1, jk, 1)
          END IF
          zfu_t(ji + 1, jj, jk) = (zfu(ji, jj, jk) + zfu(ji + 1, jj, jk) - gamma2 * (zlu_uu(ji, jj, jk, 2) + zlu_uu(ji + 1, jj, jk, 2))) * (zui - gamma1 * zl_u)
          zfv_t(ji, jj + 1, jk) = (zfv(ji, jj, jk) + zfv(ji, jj + 1, jk) - gamma2 * (zlv_vv(ji, jj, jk, 2) + zlv_vv(ji, jj + 1, jk, 2))) * (zvj - gamma1 * zl_v)
          zfuj = (zfu(ji, jj, jk) + zfu(ji, jj + 1, jk))
          zfvi = (zfv(ji, jj, jk) + zfv(ji + 1, jj, jk))
          IF (zfuj > 0) THEN
            zl_v = zlv_vu(ji, jj, jk, 1)
          ELSE
            zl_v = zlv_vu(ji + 1, jj, jk, 1)
          END IF
          IF (zfvi > 0) THEN
            zl_u = zlu_uv(ji, jj, jk, 1)
          ELSE
            zl_u = zlu_uv(ji, jj + 1, jk, 1)
          END IF
          zfv_f(ji, jj, jk) = (zfvi - gamma2 * (zlv_vu(ji, jj, jk, 2) + zlv_vu(ji + 1, jj, jk, 2))) * (un(ji, jj, jk) + un(ji, jj + 1, jk) - gamma1 * zl_u)
          zfu_f(ji, jj, jk) = (zfuj - gamma2 * (zlu_uv(ji, jj, jk, 2) + zlu_uv(ji, jj + 1, jk, 2))) * (vn(ji, jj, jk) + vn(ji + 1, jj, jk) - gamma1 * zl_v)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zfu_t(ji + 1, jj, jk) - zfu_t(ji, jj, jk) + zfv_f(ji, jj, jk) - zfv_f(ji, jj - 1, jk)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (zfu_f(ji, jj, jk) - zfu_f(ji - 1, jj, jk) + zfv_t(ji, jj + 1, jk) - zfv_t(ji, jj, jk)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    IF (l_trddyn) THEN
      !$ACC KERNELS
      zfu_uw(:, :, :) = ua(:, :, :) - zfu_uw(:, :, :)
      zfv_vw(:, :, :) = va(:, :, :) - zfv_vw(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('dyn_adv_ubs', 'r2', 0, 0)
      CALL trd_dyn(zfu_uw, zfv_vw, jpdyn_keg, kt)
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      zfu_t(:, :, :) = ua(:, :, :)
      zfv_t(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zfu_uw(ji, jj, jpk) = 0._wp
        zfv_vw(ji, jj, jpk) = 0._wp
        zfu_uw(ji, jj, 1) = 0._wp
        zfv_vw(ji, jj, 1) = 0._wp
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_linssh) THEN
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfu_uw(ji, jj, 1) = 0.5_wp * (e1e2t(ji, jj) * wn(ji, jj, 1) + e1e2t(ji + 1, jj) * wn(ji + 1, jj, 1)) * un(ji, jj, 1)
          zfv_vw(ji, jj, 1) = 0.5_wp * (e1e2t(ji, jj) * wn(ji, jj, 1) + e1e2t(ji, jj + 1) * wn(ji, jj + 1, 1)) * vn(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    DO jk = 2, jpkm1
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          zfw(ji, jj, jk) = 0.25_wp * e1e2t(ji, jj) * wn(ji, jj, jk)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfu_uw(ji, jj, jk) = (zfw(ji, jj, jk) + zfw(ji + 1, jj, jk)) * (un(ji, jj, jk) + un(ji, jj, jk - 1))
          zfv_vw(ji, jj, jk) = (zfw(ji, jj, jk) + zfw(ji, jj + 1, jk)) * (vn(ji, jj, jk) + vn(ji, jj, jk - 1))
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zfu_uw(ji, jj, jk) - zfu_uw(ji, jj, jk + 1)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (zfv_vw(ji, jj, jk) - zfv_vw(ji, jj, jk + 1)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (l_trddyn) THEN
      !$ACC KERNELS
      zfu_t(:, :, :) = ua(:, :, :) - zfu_t(:, :, :)
      zfv_t(:, :, :) = va(:, :, :) - zfv_t(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('dyn_adv_ubs', 'r3', 0, 0)
      CALL trd_dyn(zfu_t, zfv_t, jpdyn_zad, kt)
      CALL profile_psy_data3 % PostEnd
    END IF
    CALL profile_psy_data4 % PreStart('dyn_adv_ubs', 'r4', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' ubs2 adv - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = vmask, clinfo3 = 'dyn')
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE dyn_adv_ubs
END MODULE dynadv_ubs