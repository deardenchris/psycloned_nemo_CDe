MODULE zdfddm
  USE oce
  USE dom_oce
  USE zdf_oce
  USE eosbn2
  USE in_out_manager
  USE lbclnk
  USE prtctl
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_ddm
  CONTAINS
  SUBROUTINE zdf_ddm(kt, p_avm, p_avt, p_avs)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: p_avs
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zaw, zbw, zrw
    REAL(KIND = wp) :: zdt, zds
    REAL(KIND = wp) :: zinr, zrr
    REAL(KIND = wp) :: zavft, zavfs
    REAL(KIND = wp) :: zavdt, zavds
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zrau, zmsks, zmskf, zmskd1, zmskd2, zmskd3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    DO jk = 2, jpkm1
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zrw = (gdepw_n(ji, jj, jk) - gdept_n(ji, jj, jk)) / (gdept_n(ji, jj, jk - 1) - gdept_n(ji, jj, jk))
          zaw = (rab_n(ji, jj, jk, jp_tem) * (1. - zrw) + rab_n(ji, jj, jk - 1, jp_tem) * zrw) * tmask(ji, jj, jk) * tmask(ji, jj, jk - 1)
          zbw = (rab_n(ji, jj, jk, jp_sal) * (1. - zrw) + rab_n(ji, jj, jk - 1, jp_sal) * zrw) * tmask(ji, jj, jk) * tmask(ji, jj, jk - 1)
          zdt = zaw * (tsn(ji, jj, jk - 1, jp_tem) - tsn(ji, jj, jk, jp_tem))
          zds = zbw * (tsn(ji, jj, jk - 1, jp_sal) - tsn(ji, jj, jk, jp_sal))
          IF (ABS(zds) <= 1.E-20_wp) zds = 1.E-20_wp
          zrau(ji, jj) = MAX(1.E-20, zdt / zds)
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (rn2(ji, jj, jk) + 1.E-12 <= 0.) THEN
            zmsks(ji, jj) = 0._wp
          ELSE
            zmsks(ji, jj) = 1._wp
          END IF
          IF (zrau(ji, jj) <= 1.) THEN
            zmskf(ji, jj) = 0._wp
          ELSE
            zmskf(ji, jj) = 1._wp
          END IF
          IF (zrau(ji, jj) >= 1.) THEN
            zmskd1(ji, jj) = 0._wp
          ELSE
            zmskd1(ji, jj) = 1._wp
          END IF
          IF (zrau(ji, jj) >= 0.5) THEN
            zmskd2(ji, jj) = 0._wp
          ELSE
            zmskd2(ji, jj) = 1._wp
          END IF
          IF (zrau(ji, jj) <= 0.5 .OR. zrau(ji, jj) >= 1.) THEN
            zmskd3(ji, jj) = 0._wp
          ELSE
            zmskd3(ji, jj) = 1._wp
          END IF
        END DO
      END DO
      zmsks(:, :) = zmsks(:, :) * wmask(:, :, jk)
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zinr = 1._wp / zrau(ji, jj)
          zrr = zrau(ji, jj) / rn_hsbfr
          zrr = zrr * zrr
          zavfs = rn_avts / (1 + zrr * zrr * zrr) * zmsks(ji, jj) * zmskf(ji, jj)
          zavft = 0.7 * zavfs * zinr
          zavdt = 1.3635E-6 * EXP(4.6 * EXP(- 0.54 * (zinr - 1.))) * zmsks(ji, jj) * zmskd1(ji, jj)
          zavds = zavdt * zmsks(ji, jj) * ((1.85 * zrau(ji, jj) - 0.85) * zmskd3(ji, jj) + 0.15 * zrau(ji, jj) * zmskd2(ji, jj))
          p_avs(ji, jj, jk) = p_avt(ji, jj, jk) + zavfs + zavds
          p_avt(ji, jj, jk) = p_avt(ji, jj, jk) + zavft + zavdt
          p_avm(ji, jj, jk) = p_avm(ji, jj, jk) + MAX(zavft + zavdt, zavfs + zavds)
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    CALL profile_psy_data0 % PreStart('zdf_ddm', 'r0', 0, 0)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab3d_1 = avt, clinfo1 = ' ddm  - t: ', tab3d_2 = avs, clinfo2 = ' s: ', kdim = jpk)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE zdf_ddm
END MODULE zdfddm