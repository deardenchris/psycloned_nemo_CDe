MODULE zpshde
  USE oce
  USE dom_oce
  USE phycst
  USE eosbn2
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zps_hde
  PUBLIC :: zps_hde_isf
  CONTAINS
  SUBROUTINE zps_hde(kt, kjpt, pta, pgtu, pgtv, prd, pgru, pgrv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: pta
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(OUT) :: pgtu, pgtv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: prd
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT), OPTIONAL :: pgru, pgrv
    INTEGER :: ji, jj, jn
    INTEGER :: iku, ikv, ikum1, ikvm1
    REAL(KIND = wp) :: ze3wu, ze3wv, zmaxu, zmaxv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zri, zrj, zhi, zhj
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt) :: zti, ztj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (ln_timing) CALL timing_start('zps_hde')
    !$ACC KERNELS
    pgtu(:, :, :) = 0._wp
    zti(:, :, :) = 0._wp
    zhi(:, :) = 0._wp
    pgtv(:, :, :) = 0._wp
    ztj(:, :, :) = 0._wp
    zhj(:, :) = 0._wp
    DO jn = 1, kjpt
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = mbku(ji, jj)
          ikum1 = MAX(iku - 1, 1)
          ikv = mbkv(ji, jj)
          ikvm1 = MAX(ikv - 1, 1)
          ze3wu = e3w_n(ji + 1, jj, iku) - e3w_n(ji, jj, iku)
          ze3wv = e3w_n(ji, jj + 1, ikv) - e3w_n(ji, jj, ikv)
          IF (ze3wu >= 0._wp) THEN
            zmaxu = ze3wu / e3w_n(ji + 1, jj, iku)
            zti(ji, jj, jn) = pta(ji + 1, jj, iku, jn) + zmaxu * (pta(ji + 1, jj, ikum1, jn) - pta(ji + 1, jj, iku, jn))
            pgtu(ji, jj, jn) = umask(ji, jj, 1) * (zti(ji, jj, jn) - pta(ji, jj, iku, jn))
          ELSE
            zmaxu = - ze3wu / e3w_n(ji, jj, iku)
            zti(ji, jj, jn) = pta(ji, jj, iku, jn) + zmaxu * (pta(ji, jj, ikum1, jn) - pta(ji, jj, iku, jn))
            pgtu(ji, jj, jn) = umask(ji, jj, 1) * (pta(ji + 1, jj, iku, jn) - zti(ji, jj, jn))
          END IF
          IF (ze3wv >= 0._wp) THEN
            zmaxv = ze3wv / e3w_n(ji, jj + 1, ikv)
            ztj(ji, jj, jn) = pta(ji, jj + 1, ikv, jn) + zmaxv * (pta(ji, jj + 1, ikvm1, jn) - pta(ji, jj + 1, ikv, jn))
            pgtv(ji, jj, jn) = vmask(ji, jj, 1) * (ztj(ji, jj, jn) - pta(ji, jj, ikv, jn))
          ELSE
            zmaxv = - ze3wv / e3w_n(ji, jj, ikv)
            ztj(ji, jj, jn) = pta(ji, jj, ikv, jn) + zmaxv * (pta(ji, jj, ikvm1, jn) - pta(ji, jj, ikv, jn))
            pgtv(ji, jj, jn) = vmask(ji, jj, 1) * (pta(ji, jj + 1, ikv, jn) - ztj(ji, jj, jn))
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('zpshde', pgtu(:, :, :), 'U', - 1., pgtv(:, :, :), 'V', - 1.)
    IF (PRESENT(prd)) THEN
      !$ACC KERNELS
      pgru(:, :) = 0._wp
      pgrv(:, :) = 0._wp
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = mbku(ji, jj)
          ikv = mbkv(ji, jj)
          ze3wu = e3w_n(ji + 1, jj, iku) - e3w_n(ji, jj, iku)
          ze3wv = e3w_n(ji, jj + 1, ikv) - e3w_n(ji, jj, ikv)
          IF (ze3wu >= 0._wp) THEN
            zhi(ji, jj) = gdept_n(ji, jj, iku)
          ELSE
            zhi(ji, jj) = gdept_n(ji + 1, jj, iku)
          END IF
          IF (ze3wv >= 0._wp) THEN
            zhj(ji, jj) = gdept_n(ji, jj, ikv)
          ELSE
            zhj(ji, jj) = gdept_n(ji, jj + 1, ikv)
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('zps_hde', 'r0', 0, 0)
      CALL eos(zti, zhi, zri)
      CALL eos(ztj, zhj, zrj)
      CALL profile_psy_data0 % PostEnd
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = mbku(ji, jj)
          ikv = mbkv(ji, jj)
          ze3wu = e3w_n(ji + 1, jj, iku) - e3w_n(ji, jj, iku)
          ze3wv = e3w_n(ji, jj + 1, ikv) - e3w_n(ji, jj, ikv)
          IF (ze3wu >= 0._wp) THEN
            pgru(ji, jj) = umask(ji, jj, 1) * (zri(ji, jj) - prd(ji, jj, iku))
          ELSE
            pgru(ji, jj) = umask(ji, jj, 1) * (prd(ji + 1, jj, iku) - zri(ji, jj))
          END IF
          IF (ze3wv >= 0._wp) THEN
            pgrv(ji, jj) = vmask(ji, jj, 1) * (zrj(ji, jj) - prd(ji, jj, ikv))
          ELSE
            pgrv(ji, jj) = vmask(ji, jj, 1) * (prd(ji, jj + 1, ikv) - zrj(ji, jj))
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi('zpshde', pgru, 'U', - 1., pgrv, 'V', - 1.)
    END IF
    IF (ln_timing) CALL timing_stop('zps_hde')
  END SUBROUTINE zps_hde
  SUBROUTINE zps_hde_isf(kt, kjpt, pta, pgtu, pgtv, pgtui, pgtvi, prd, pgru, pgrv, pgrui, pgrvi)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kjpt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: pta
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(OUT) :: pgtu, pgtv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(OUT) :: pgtui, pgtvi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: prd
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT), OPTIONAL :: pgru, pgrv
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT), OPTIONAL :: pgrui, pgrvi
    INTEGER :: ji, jj, jn
    INTEGER :: iku, ikv, ikum1, ikvm1, ikup1, ikvp1
    REAL(KIND = wp) :: ze3wu, ze3wv, zmaxu, zmaxv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zri, zrj, zhi, zhj
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt) :: zti, ztj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    IF (ln_timing) CALL timing_start('zps_hde_isf')
    !$ACC KERNELS
    pgtu(:, :, :) = 0._wp
    pgtv(:, :, :) = 0._wp
    pgtui(:, :, :) = 0._wp
    pgtvi(:, :, :) = 0._wp
    zti(:, :, :) = 0._wp
    ztj(:, :, :) = 0._wp
    zhi(:, :) = 0._wp
    zhj(:, :) = 0._wp
    DO jn = 1, kjpt
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = mbku(ji, jj)
          ikum1 = MAX(iku - 1, 1)
          ikv = mbkv(ji, jj)
          ikvm1 = MAX(ikv - 1, 1)
          ze3wu = gdept_n(ji + 1, jj, iku) - gdept_n(ji, jj, iku)
          ze3wv = gdept_n(ji, jj + 1, ikv) - gdept_n(ji, jj, ikv)
          IF (ze3wu >= 0._wp) THEN
            zmaxu = ze3wu / e3w_n(ji + 1, jj, iku)
            zti(ji, jj, jn) = pta(ji + 1, jj, iku, jn) + zmaxu * (pta(ji + 1, jj, ikum1, jn) - pta(ji + 1, jj, iku, jn))
            pgtu(ji, jj, jn) = ssumask(ji, jj) * (zti(ji, jj, jn) - pta(ji, jj, iku, jn))
          ELSE
            zmaxu = - ze3wu / e3w_n(ji, jj, iku)
            zti(ji, jj, jn) = pta(ji, jj, iku, jn) + zmaxu * (pta(ji, jj, ikum1, jn) - pta(ji, jj, iku, jn))
            pgtu(ji, jj, jn) = ssumask(ji, jj) * (pta(ji + 1, jj, iku, jn) - zti(ji, jj, jn))
          END IF
          IF (ze3wv >= 0._wp) THEN
            zmaxv = ze3wv / e3w_n(ji, jj + 1, ikv)
            ztj(ji, jj, jn) = pta(ji, jj + 1, ikv, jn) + zmaxv * (pta(ji, jj + 1, ikvm1, jn) - pta(ji, jj + 1, ikv, jn))
            pgtv(ji, jj, jn) = ssvmask(ji, jj) * (ztj(ji, jj, jn) - pta(ji, jj, ikv, jn))
          ELSE
            zmaxv = - ze3wv / e3w_n(ji, jj, ikv)
            ztj(ji, jj, jn) = pta(ji, jj, ikv, jn) + zmaxv * (pta(ji, jj, ikvm1, jn) - pta(ji, jj, ikv, jn))
            pgtv(ji, jj, jn) = ssvmask(ji, jj) * (pta(ji, jj + 1, ikv, jn) - ztj(ji, jj, jn))
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('zpshde', pgtu(:, :, :), 'U', - 1., pgtv(:, :, :), 'V', - 1.)
    IF (PRESENT(prd)) THEN
      !$ACC KERNELS
      pgru(:, :) = 0.0_wp
      pgrv(:, :) = 0.0_wp
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = mbku(ji, jj)
          ikv = mbkv(ji, jj)
          ze3wu = gdept_n(ji + 1, jj, iku) - gdept_n(ji, jj, iku)
          ze3wv = gdept_n(ji, jj + 1, ikv) - gdept_n(ji, jj, ikv)
          IF (ze3wu >= 0._wp) THEN
            zhi(ji, jj) = gdept_n(ji, jj, iku)
          ELSE
            zhi(ji, jj) = gdept_n(ji + 1, jj, iku)
          END IF
          IF (ze3wv >= 0._wp) THEN
            zhj(ji, jj) = gdept_n(ji, jj, ikv)
          ELSE
            zhj(ji, jj) = gdept_n(ji, jj + 1, ikv)
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('zps_hde_isf', 'r0', 0, 0)
      CALL eos(zti, zhi, zri)
      CALL eos(ztj, zhj, zrj)
      CALL profile_psy_data0 % PostEnd
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = mbku(ji, jj)
          ikv = mbkv(ji, jj)
          ze3wu = gdept_n(ji + 1, jj, iku) - gdept_n(ji, jj, iku)
          ze3wv = gdept_n(ji, jj + 1, ikv) - gdept_n(ji, jj, ikv)
          IF (ze3wu >= 0._wp) THEN
            pgru(ji, jj) = ssumask(ji, jj) * (zri(ji, jj) - prd(ji, jj, iku))
          ELSE
            pgru(ji, jj) = ssumask(ji, jj) * (prd(ji + 1, jj, iku) - zri(ji, jj))
          END IF
          IF (ze3wv >= 0._wp) THEN
            pgrv(ji, jj) = ssvmask(ji, jj) * (zrj(ji, jj) - prd(ji, jj, ikv))
          ELSE
            pgrv(ji, jj) = ssvmask(ji, jj) * (prd(ji, jj + 1, ikv) - zrj(ji, jj))
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi('zpshde', pgru, 'U', - 1., pgrv, 'V', - 1.)
    END IF
    !$ACC KERNELS
    DO jn = 1, kjpt
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = miku(ji, jj)
          ikup1 = miku(ji, jj) + 1
          ikv = mikv(ji, jj)
          ikvp1 = mikv(ji, jj) + 1
          ze3wu = gdept_n(ji, jj, iku) - gdept_n(ji + 1, jj, iku)
          ze3wv = gdept_n(ji, jj, ikv) - gdept_n(ji, jj + 1, ikv)
          IF (ze3wu >= 0._wp) THEN
            zmaxu = ze3wu / e3w_n(ji + 1, jj, ikup1)
            zti(ji, jj, jn) = pta(ji + 1, jj, iku, jn) + zmaxu * (pta(ji + 1, jj, ikup1, jn) - pta(ji + 1, jj, iku, jn))
            pgtui(ji, jj, jn) = ssumask(ji, jj) * (zti(ji, jj, jn) - pta(ji, jj, iku, jn))
          ELSE
            zmaxu = - ze3wu / e3w_n(ji, jj, ikup1)
            zti(ji, jj, jn) = pta(ji, jj, iku, jn) + zmaxu * (pta(ji, jj, ikup1, jn) - pta(ji, jj, iku, jn))
            pgtui(ji, jj, jn) = ssumask(ji, jj) * (pta(ji + 1, jj, iku, jn) - zti(ji, jj, jn))
          END IF
          IF (ze3wv >= 0._wp) THEN
            zmaxv = ze3wv / e3w_n(ji, jj + 1, ikvp1)
            ztj(ji, jj, jn) = pta(ji, jj + 1, ikv, jn) + zmaxv * (pta(ji, jj + 1, ikvp1, jn) - pta(ji, jj + 1, ikv, jn))
            pgtvi(ji, jj, jn) = ssvmask(ji, jj) * (ztj(ji, jj, jn) - pta(ji, jj, ikv, jn))
          ELSE
            zmaxv = - ze3wv / e3w_n(ji, jj, ikvp1)
            ztj(ji, jj, jn) = pta(ji, jj, ikv, jn) + zmaxv * (pta(ji, jj, ikvp1, jn) - pta(ji, jj, ikv, jn))
            pgtvi(ji, jj, jn) = ssvmask(ji, jj) * (pta(ji, jj + 1, ikv, jn) - ztj(ji, jj, jn))
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('zpshde', pgtui(:, :, :), 'U', - 1., pgtvi(:, :, :), 'V', - 1.)
    IF (PRESENT(prd)) THEN
      !$ACC KERNELS
      pgrui(:, :) = 0.0_wp
      pgrvi(:, :) = 0.0_wp
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = miku(ji, jj)
          ikv = mikv(ji, jj)
          ze3wu = gdept_n(ji, jj, iku) - gdept_n(ji + 1, jj, iku)
          ze3wv = gdept_n(ji, jj, ikv) - gdept_n(ji, jj + 1, ikv)
          IF (ze3wu >= 0._wp) THEN
            zhi(ji, jj) = gdept_n(ji, jj, iku)
          ELSE
            zhi(ji, jj) = gdept_n(ji + 1, jj, iku)
          END IF
          IF (ze3wv >= 0._wp) THEN
            zhj(ji, jj) = gdept_n(ji, jj, ikv)
          ELSE
            zhj(ji, jj) = gdept_n(ji, jj + 1, ikv)
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('zps_hde_isf', 'r1', 0, 0)
      CALL eos(zti, zhi, zri)
      CALL eos(ztj, zhj, zrj)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          iku = miku(ji, jj)
          ikv = mikv(ji, jj)
          ze3wu = gdept_n(ji, jj, iku) - gdept_n(ji + 1, jj, iku)
          ze3wv = gdept_n(ji, jj, ikv) - gdept_n(ji, jj + 1, ikv)
          IF (ze3wu >= 0._wp) THEN
            pgrui(ji, jj) = ssumask(ji, jj) * (zri(ji, jj) - prd(ji, jj, iku))
          ELSE
            pgrui(ji, jj) = ssumask(ji, jj) * (prd(ji + 1, jj, iku) - zri(ji, jj))
          END IF
          IF (ze3wv >= 0._wp) THEN
            pgrvi(ji, jj) = ssvmask(ji, jj) * (zrj(ji, jj) - prd(ji, jj, ikv))
          ELSE
            pgrvi(ji, jj) = ssvmask(ji, jj) * (prd(ji, jj + 1, ikv) - zrj(ji, jj))
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi('zpshde', pgrui, 'U', - 1., pgrvi, 'V', - 1.)
    END IF
    IF (ln_timing) CALL timing_stop('zps_hde_isf')
  END SUBROUTINE zps_hde_isf
END MODULE zpshde