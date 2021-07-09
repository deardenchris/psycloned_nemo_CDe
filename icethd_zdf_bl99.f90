MODULE icethd_zdf_BL99
  USE dom_oce
  USE phycst
  USE ice
  USE ice1D
  USE icevar
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_zdf_BL99
  CONTAINS
  SUBROUTINE ice_thd_zdf_BL99(k_cnd)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: k_cnd
    INTEGER :: ji, jk
    INTEGER :: jm
    INTEGER :: jm_mint, jm_maxt
    INTEGER :: iconv
    INTEGER :: iconv_max = 50
    INTEGER, DIMENSION(jpij) :: jm_min
    INTEGER, DIMENSION(jpij) :: jm_max
    LOGICAL, DIMENSION(jpij) :: l_T_converged
    REAL(KIND = wp) :: zg1s = 2._wp
    REAL(KIND = wp) :: zg1 = 2._wp
    REAL(KIND = wp) :: zgamma = 18009._wp
    REAL(KIND = wp) :: zbeta = 0.117_wp
    REAL(KIND = wp) :: zraext_s = 10._wp
    REAL(KIND = wp) :: zkimin = 0.10_wp
    REAL(KIND = wp) :: ztsu_err = 1.E-5_wp
    REAL(KIND = wp) :: zdti_bnd = 1.E-4_wp
    REAL(KIND = wp) :: zhs_min = 0.01_wp
    REAL(KIND = wp) :: ztmelts
    REAL(KIND = wp) :: zdti_max
    REAL(KIND = wp) :: zcpi
    REAL(KIND = wp) :: zhfx_err, zdq
    REAL(KIND = wp) :: zfac
    REAL(KIND = wp), DIMENSION(jpij) :: isnow
    REAL(KIND = wp), DIMENSION(jpij) :: ztsub
    REAL(KIND = wp), DIMENSION(jpij) :: zh_i, z1_h_i
    REAL(KIND = wp), DIMENSION(jpij) :: zh_s, z1_h_s
    REAL(KIND = wp), DIMENSION(jpij) :: zqns_ice_b
    REAL(KIND = wp), DIMENSION(jpij) :: zfnet
    REAL(KIND = wp), DIMENSION(jpij) :: zdqns_ice_b
    REAL(KIND = wp), DIMENSION(jpij) :: ztsuold
    REAL(KIND = wp), DIMENSION(jpij, nlay_i) :: ztiold
    REAL(KIND = wp), DIMENSION(jpij, nlay_s) :: ztsold
    REAL(KIND = wp), DIMENSION(jpij, nlay_i) :: ztib
    REAL(KIND = wp), DIMENSION(jpij, nlay_s) :: ztsb
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: ztcond_i
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: ztcond_i_cp
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: zradtr_i
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: zradab_i
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: zkappa_i
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: zeta_i
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_s) :: zradtr_s
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_s) :: zradab_s
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_s) :: zkappa_s
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_s) :: zeta_s
    REAL(KIND = wp), DIMENSION(jpij, nlay_i + 3) :: zindterm
    REAL(KIND = wp), DIMENSION(jpij, nlay_i + 3) :: zindtbis
    REAL(KIND = wp), DIMENSION(jpij, nlay_i + 3) :: zdiagbis
    REAL(KIND = wp), DIMENSION(jpij, nlay_i + 3, 3) :: ztrid
    REAL(KIND = wp), DIMENSION(jpij) :: zq_ini
    REAL(KIND = wp), DIMENSION(jpij) :: zghe
    REAL(KIND = wp) :: zepsilon
    REAL(KIND = wp) :: zhe
    REAL(KIND = wp) :: zcnd_i
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    !$ACC KERNELS
    DO ji = 1, npti
      zq_ini(ji) = (SUM(e_i_1d(ji, 1 : nlay_i)) * h_i_1d(ji) * r1_nlay_i + SUM(e_s_1d(ji, 1 : nlay_s)) * h_s_1d(ji) * r1_nlay_s)
    END DO
    DO ji = 1, npti
      isnow(ji) = 1._wp - MAX(0._wp, SIGN(1._wp, - h_s_1d(ji)))
      zh_i(ji) = h_i_1d(ji) * r1_nlay_i
      zh_s(ji) = h_s_1d(ji) * r1_nlay_s
    END DO
!    !$ACC END KERNELS
    ! CDe continue kernels region here rather than split it up
!    CALL profile_psy_data0 % PreStart('ice_thd_zdf_bl99', 'r0', 0, 0)
    WHERE (zh_i(1 : npti) >= epsi10)
      z1_h_i(1 : npti) = 1._wp / zh_i(1 : npti)
    ELSEWHERE
      z1_h_i(1 : npti) = 0._wp
    END WHERE
    WHERE (zh_s(1 : npti) > 0._wp) zh_s(1 : npti) = MAX(zhs_min * r1_nlay_s, zh_s(1 : npti))
    WHERE (zh_s(1 : npti) > 0._wp)
      z1_h_s(1 : npti) = 1._wp / zh_s(1 : npti)
    ELSEWHERE
      z1_h_s(1 : npti) = 0._wp
    END WHERE
!    CALL profile_psy_data0 % PostEnd
!    !$ACC KERNELS
    IF (k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU) THEN
      ztsub(1 : npti) = t_su_1d(1 : npti)
      ztsuold(1 : npti) = t_su_1d(1 : npti)
      t_su_1d(1 : npti) = MIN(t_su_1d(1 : npti), rt0 - ztsu_err)
      zdqns_ice_b(1 : npti) = dqns_ice_1d(1 : npti)
      zqns_ice_b(1 : npti) = qns_ice_1d(1 : npti)
    END IF
    ztsold(1 : npti, :) = t_s_1d(1 : npti, :)
    ztiold(1 : npti, :) = t_i_1d(1 : npti, :)
    zradtr_s(1 : npti, 0) = qtr_ice_top_1d(1 : npti)
    DO jk = 1, nlay_s
      DO ji = 1, npti
        zradtr_s(ji, jk) = zradtr_s(ji, 0) * EXP(- zraext_s * h_s_1d(ji) * r1_nlay_s * REAL(jk))
        zradab_s(ji, jk) = zradtr_s(ji, jk - 1) - zradtr_s(ji, jk)
      END DO
    END DO
    zradtr_i(1 : npti, 0) = zradtr_s(1 : npti, nlay_s) * isnow(1 : npti) + qtr_ice_top_1d(1 : npti) * (1._wp - isnow(1 : npti))
    DO jk = 1, nlay_i
      DO ji = 1, npti
        zradtr_i(ji, jk) = zradtr_i(ji, 0) * EXP(- rn_kappa_i * zh_i(ji) * REAL(jk))
        zradab_i(ji, jk) = zradtr_i(ji, jk - 1) - zradtr_i(ji, jk)
      END DO
    END DO
    qtr_ice_bot_1d(1 : npti) = zradtr_i(1 : npti, nlay_i)
    iconv = 0
    l_T_converged(:) = .FALSE.
    !$ACC END KERNELS
    ! CALL profile_psy_data1 % PreStart('ice_thd_zdf_bl99', 'r1', 0, 0)
    ! !$ACC KERNELS 
    DO WHILE ((.NOT. ALL(l_T_converged(1 : npti))) .AND. iconv < iconv_max)
      iconv = iconv + 1
      ztib(1 : npti, :) = t_i_1d(1 : npti, :)
      ztsb(1 : npti, :) = t_s_1d(1 : npti, :)
      IF (ln_cndi_U64) THEN
        DO ji = 1, npti
          ztcond_i_cp(ji, 0) = rcnd_i + zbeta * sz_i_1d(ji, 1) / MIN(- epsi10, t_i_1d(ji, 1) - rt0)
          ztcond_i_cp(ji, nlay_i) = rcnd_i + zbeta * sz_i_1d(ji, nlay_i) / MIN(- epsi10, t_bo_1d(ji) - rt0)
        END DO
        DO jk = 1, nlay_i - 1
          DO ji = 1, npti
            ztcond_i_cp(ji, jk) = rcnd_i + zbeta * 0.5_wp * (sz_i_1d(ji, jk) + sz_i_1d(ji, jk + 1)) / MIN(- epsi10, 0.5_wp * (t_i_1d(ji, jk) + t_i_1d(ji, jk + 1)) - rt0)
          END DO
        END DO
      ELSE IF (ln_cndi_P07) THEN
        DO ji = 1, npti
          ztcond_i_cp(ji, 0) = rcnd_i + 0.09_wp * sz_i_1d(ji, 1) / MIN(- epsi10, t_i_1d(ji, 1) - rt0) - 0.011_wp * (t_i_1d(ji, 1) - rt0)
          ztcond_i_cp(ji, nlay_i) = rcnd_i + 0.09_wp * sz_i_1d(ji, nlay_i) / MIN(- epsi10, t_bo_1d(ji) - rt0) - 0.011_wp * (t_bo_1d(ji) - rt0)
        END DO
        DO jk = 1, nlay_i - 1
          DO ji = 1, npti
            ztcond_i_cp(ji, jk) = rcnd_i + 0.09_wp * 0.5_wp * (sz_i_1d(ji, jk) + sz_i_1d(ji, jk + 1)) / MIN(- epsi10, 0.5_wp * (t_i_1d(ji, jk) + t_i_1d(ji, jk + 1)) - rt0) - 0.011_wp * (0.5_wp * (t_i_1d(ji, jk) + t_i_1d(ji, jk + 1)) - rt0)
          END DO
        END DO
      END IF
      DO ji = 1, npti
        IF (.NOT. l_T_converged(ji)) ztcond_i(ji, :) = MAX(zkimin, ztcond_i_cp(ji, :))
      END DO
      zghe(1 : npti) = 1._wp
      IF (ln_virtual_itd) THEN
        zepsilon = 0.1_wp
        DO ji = 1, npti
          zcnd_i = SUM(ztcond_i(ji, :)) / REAL(nlay_i + 1, wp)
          zhe = (rn_cnd_s * h_i_1d(ji) + zcnd_i * h_s_1d(ji)) / (rn_cnd_s + zcnd_i)
          IF (zhe >= zepsilon * 0.5_wp * EXP(1._wp)) zghe(ji) = MIN(2._wp, 0.5_wp * (1._wp + LOG(2._wp * zhe / zepsilon)))
        END DO
      END IF
      DO jk = 0, nlay_s - 1
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) zkappa_s(ji, jk) = zghe(ji) * rn_cnd_s * z1_h_s(ji)
        END DO
      END DO
      DO ji = 1, npti
        IF (.NOT. l_T_converged(ji)) THEN
          zfac = 0.5_wp * (ztcond_i(ji, 0) * zh_s(ji) + rn_cnd_s * zh_i(ji))
          IF (zfac > epsi10) THEN
            zkappa_s(ji, nlay_s) = zghe(ji) * rn_cnd_s * ztcond_i(ji, 0) / zfac
          ELSE
            zkappa_s(ji, nlay_s) = 0._wp
          END IF
        END IF
      END DO
      DO jk = 0, nlay_i
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) zkappa_i(ji, jk) = zghe(ji) * ztcond_i(ji, jk) * z1_h_i(ji)
        END DO
      END DO
      DO ji = 1, npti
        IF (.NOT. l_T_converged(ji)) zkappa_i(ji, 0) = zkappa_s(ji, nlay_s) * isnow(ji) + zkappa_i(ji, 0) * (1._wp - isnow(ji))
      END DO
      DO jk = 1, nlay_i
        DO ji = 1, npti
          zcpi = rcpi + zgamma * sz_i_1d(ji, jk) / MAX((t_i_1d(ji, jk) - rt0) * (ztiold(ji, jk) - rt0), epsi10)
          zeta_i(ji, jk) = rdt_ice * r1_rhoi * z1_h_i(ji) / MAX(epsi10, zcpi)
        END DO
      END DO
      DO jk = 1, nlay_s
        DO ji = 1, npti
          zeta_s(ji, jk) = rdt_ice * r1_rhos * r1_rcpi * z1_h_s(ji)
        END DO
      END DO
      IF (k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU) THEN
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) qns_ice_1d(ji) = qns_ice_1d(ji) + dqns_ice_1d(ji) * (t_su_1d(ji) - ztsub(ji))
        END DO
        DO ji = 1, npti
          zfnet(ji) = qsr_ice_1d(ji) - qtr_ice_top_1d(ji) + qns_ice_1d(ji)
        END DO
        ztrid(1 : npti, :, :) = 0._wp
        zindterm(1 : npti, :) = 0._wp
        zindtbis(1 : npti, :) = 0._wp
        zdiagbis(1 : npti, :) = 0._wp
        DO jm = nlay_s + 2, nlay_s + nlay_i
          DO ji = 1, npti
            jk = jm - nlay_s - 1
            ztrid(ji, jm, 1) = - zeta_i(ji, jk) * zkappa_i(ji, jk - 1)
            ztrid(ji, jm, 2) = 1._wp + zeta_i(ji, jk) * (zkappa_i(ji, jk - 1) + zkappa_i(ji, jk))
            ztrid(ji, jm, 3) = - zeta_i(ji, jk) * zkappa_i(ji, jk)
            zindterm(ji, jm) = ztiold(ji, jk) + zeta_i(ji, jk) * zradab_i(ji, jk)
          END DO
        END DO
        jm = nlay_s + nlay_i + 1
        DO ji = 1, npti
          ztrid(ji, jm, 1) = - zeta_i(ji, nlay_i) * zkappa_i(ji, nlay_i - 1)
          ztrid(ji, jm, 2) = 1._wp + zeta_i(ji, nlay_i) * (zkappa_i(ji, nlay_i - 1) + zkappa_i(ji, nlay_i) * zg1)
          ztrid(ji, jm, 3) = 0._wp
          zindterm(ji, jm) = ztiold(ji, nlay_i) + zeta_i(ji, nlay_i) * (zradab_i(ji, nlay_i) + zkappa_i(ji, nlay_i) * zg1 * t_bo_1d(ji))
        END DO
        DO ji = 1, npti
          IF (h_s_1d(ji) > 0._wp) THEN
            DO jm = 3, nlay_s + 1
              jk = jm - 1
              ztrid(ji, jm, 1) = - zeta_s(ji, jk) * zkappa_s(ji, jk - 1)
              ztrid(ji, jm, 2) = 1._wp + zeta_s(ji, jk) * (zkappa_s(ji, jk - 1) + zkappa_s(ji, jk))
              ztrid(ji, jm, 3) = - zeta_s(ji, jk) * zkappa_s(ji, jk)
              zindterm(ji, jm) = ztsold(ji, jk) + zeta_s(ji, jk) * zradab_s(ji, jk)
            END DO
            IF (nlay_i == 1) THEN
              ztrid(ji, nlay_s + 2, 3) = 0._wp
              zindterm(ji, nlay_s + 2) = zindterm(ji, nlay_s + 2) + zeta_i(ji, 1) * zkappa_i(ji, 1) * t_bo_1d(ji)
            END IF
            IF (t_su_1d(ji) < rt0) THEN
              jm_min(ji) = 1
              jm_max(ji) = nlay_i + nlay_s + 1
              ztrid(ji, 1, 1) = 0._wp
              ztrid(ji, 1, 2) = zdqns_ice_b(ji) - zg1s * zkappa_s(ji, 0)
              ztrid(ji, 1, 3) = zg1s * zkappa_s(ji, 0)
              zindterm(ji, 1) = zdqns_ice_b(ji) * t_su_1d(ji) - zfnet(ji)
              ztrid(ji, 2, 1) = - zeta_s(ji, 1) * zkappa_s(ji, 0) * zg1s
              ztrid(ji, 2, 2) = 1._wp + zeta_s(ji, 1) * (zkappa_s(ji, 1) + zkappa_s(ji, 0) * zg1s)
              ztrid(ji, 2, 3) = - zeta_s(ji, 1) * zkappa_s(ji, 1)
              zindterm(ji, 2) = ztsold(ji, 1) + zeta_s(ji, 1) * zradab_s(ji, 1)
            ELSE
              jm_min(ji) = 2
              jm_max(ji) = nlay_i + nlay_s + 1
              ztrid(ji, 2, 1) = 0._wp
              ztrid(ji, 2, 2) = 1._wp + zeta_s(ji, 1) * (zkappa_s(ji, 1) + zkappa_s(ji, 0) * zg1s)
              ztrid(ji, 2, 3) = - zeta_s(ji, 1) * zkappa_s(ji, 1)
              zindterm(ji, 2) = ztsold(ji, 1) + zeta_s(ji, 1) * (zradab_s(ji, 1) + zkappa_s(ji, 0) * zg1s * t_su_1d(ji))
            END IF
          ELSE
            IF (t_su_1d(ji) < rt0) THEN
              jm_min(ji) = nlay_s + 1
              jm_max(ji) = nlay_i + nlay_s + 1
              ztrid(ji, jm_min(ji), 1) = 0._wp
              ztrid(ji, jm_min(ji), 2) = zdqns_ice_b(ji) - zkappa_i(ji, 0) * zg1
              ztrid(ji, jm_min(ji), 3) = zkappa_i(ji, 0) * zg1
              zindterm(ji, jm_min(ji)) = zdqns_ice_b(ji) * t_su_1d(ji) - zfnet(ji)
              ztrid(ji, jm_min(ji) + 1, 1) = - zeta_i(ji, 1) * zkappa_i(ji, 0) * zg1
              ztrid(ji, jm_min(ji) + 1, 2) = 1._wp + zeta_i(ji, 1) * (zkappa_i(ji, 1) + zkappa_i(ji, 0) * zg1)
              ztrid(ji, jm_min(ji) + 1, 3) = - zeta_i(ji, 1) * zkappa_i(ji, 1)
              zindterm(ji, jm_min(ji) + 1) = ztiold(ji, 1) + zeta_i(ji, 1) * zradab_i(ji, 1)
              IF (nlay_i == 1) THEN
                ztrid(ji, jm_min(ji), 1) = 0._wp
                ztrid(ji, jm_min(ji), 2) = zdqns_ice_b(ji) - zkappa_i(ji, 0) * 2._wp
                ztrid(ji, jm_min(ji), 3) = zkappa_i(ji, 0) * 2._wp
                ztrid(ji, jm_min(ji) + 1, 1) = - zeta_i(ji, 1) * zkappa_i(ji, 0) * 2._wp
                ztrid(ji, jm_min(ji) + 1, 2) = 1._wp + zeta_i(ji, 1) * (zkappa_i(ji, 0) * 2._wp + zkappa_i(ji, 1))
                ztrid(ji, jm_min(ji) + 1, 3) = 0._wp
                zindterm(ji, jm_min(ji) + 1) = ztiold(ji, 1) + zeta_i(ji, 1) * (zradab_i(ji, 1) + zkappa_i(ji, 1) * t_bo_1d(ji))
              END IF
            ELSE
              jm_min(ji) = nlay_s + 2
              jm_max(ji) = nlay_i + nlay_s + 1
              ztrid(ji, jm_min(ji), 1) = 0._wp
              ztrid(ji, jm_min(ji), 2) = 1._wp + zeta_i(ji, 1) * (zkappa_i(ji, 1) + zkappa_i(ji, 0) * zg1)
              ztrid(ji, jm_min(ji), 3) = - zeta_i(ji, 1) * zkappa_i(ji, 1)
              zindterm(ji, jm_min(ji)) = ztiold(ji, 1) + zeta_i(ji, 1) * (zradab_i(ji, 1) + zkappa_i(ji, 0) * zg1 * t_su_1d(ji))
              IF (nlay_i == 1) THEN
                ztrid(ji, jm_min(ji), 1) = 0._wp
                ztrid(ji, jm_min(ji), 2) = 1._wp + zeta_i(ji, 1) * (zkappa_i(ji, 0) * 2._wp + zkappa_i(ji, 1))
                ztrid(ji, jm_min(ji), 3) = 0._wp
                zindterm(ji, jm_min(ji)) = ztiold(ji, 1) + zeta_i(ji, 1) * (zradab_i(ji, 1) + zkappa_i(ji, 1) * t_bo_1d(ji)) + t_su_1d(ji) * zeta_i(ji, 1) * zkappa_i(ji, 0) * 2._wp
              END IF
            END IF
          END IF
          zindtbis(ji, jm_min(ji)) = zindterm(ji, jm_min(ji))
          zdiagbis(ji, jm_min(ji)) = ztrid(ji, jm_min(ji), 2)
        END DO
        jm_maxt = 0
        jm_mint = nlay_i + 5
        DO ji = 1, npti
          jm_mint = MIN(jm_min(ji), jm_mint)
          jm_maxt = MAX(jm_max(ji), jm_maxt)
        END DO
        DO jk = jm_mint + 1, jm_maxt
          DO ji = 1, npti
            jm = MIN(MAX(jm_min(ji) + 1, jk), jm_max(ji))
            zdiagbis(ji, jm) = ztrid(ji, jm, 2) - ztrid(ji, jm, 1) * ztrid(ji, jm - 1, 3) / zdiagbis(ji, jm - 1)
            zindtbis(ji, jm) = zindterm(ji, jm) - ztrid(ji, jm, 1) * zindtbis(ji, jm - 1) / zdiagbis(ji, jm - 1)
          END DO
        END DO
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) t_i_1d(ji, nlay_i) = zindtbis(ji, jm_max(ji)) / zdiagbis(ji, jm_max(ji))
        END DO
        DO jm = nlay_i + nlay_s, nlay_s + 2, - 1
          DO ji = 1, npti
            jk = jm - nlay_s - 1
            IF (.NOT. l_T_converged(ji)) t_i_1d(ji, jk) = (zindtbis(ji, jm) - ztrid(ji, jm, 3) * t_i_1d(ji, jk + 1)) / zdiagbis(ji, jm)
          END DO
        END DO
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) THEN
            IF (h_s_1d(ji) > 0._wp) THEN
              t_s_1d(ji, nlay_s) = (zindtbis(ji, nlay_s + 1) - ztrid(ji, nlay_s + 1, 3) * t_i_1d(ji, 1)) / zdiagbis(ji, nlay_s + 1)
            END IF
            ztsub(ji) = t_su_1d(ji)
            IF (t_su_1d(ji) < rt0) THEN
              t_su_1d(ji) = (zindtbis(ji, jm_min(ji)) - ztrid(ji, jm_min(ji), 3) * (isnow(ji) * t_s_1d(ji, 1) + (1._wp - isnow(ji)) * t_i_1d(ji, 1))) / zdiagbis(ji, jm_min(ji))
            END IF
          END IF
        END DO
        DO ji = 1, npti
          zdti_max = 0._wp
          IF (.NOT. l_T_converged(ji)) THEN
            t_su_1d(ji) = MAX(MIN(t_su_1d(ji), rt0), rt0 - 100._wp)
            zdti_max = MAX(zdti_max, ABS(t_su_1d(ji) - ztsub(ji)))
            t_s_1d(ji, 1 : nlay_s) = MAX(MIN(t_s_1d(ji, 1 : nlay_s), rt0), rt0 - 100._wp)
            zdti_max = MAX(zdti_max, MAXVAL(ABS(t_s_1d(ji, 1 : nlay_s) - ztsb(ji, 1 : nlay_s))))
            DO jk = 1, nlay_i
              ztmelts = - rTmlt * sz_i_1d(ji, jk) + rt0
              t_i_1d(ji, jk) = MAX(MIN(t_i_1d(ji, jk), ztmelts), rt0 - 100._wp)
              zdti_max = MAX(zdti_max, ABS(t_i_1d(ji, jk) - ztib(ji, jk)))
            END DO
            IF (zdti_max < zdti_bnd) l_T_converged(ji) = .TRUE.
          END IF
        END DO
      ELSE IF (k_cnd == np_cnd_ON) THEN
        ztrid(1 : npti, :, :) = 0._wp
        zindterm(1 : npti, :) = 0._wp
        zindtbis(1 : npti, :) = 0._wp
        zdiagbis(1 : npti, :) = 0._wp
        DO jm = nlay_s + 2, nlay_s + nlay_i
          DO ji = 1, npti
            jk = jm - nlay_s - 1
            ztrid(ji, jm, 1) = - zeta_i(ji, jk) * zkappa_i(ji, jk - 1)
            ztrid(ji, jm, 2) = 1._wp + zeta_i(ji, jk) * (zkappa_i(ji, jk - 1) + zkappa_i(ji, jk))
            ztrid(ji, jm, 3) = - zeta_i(ji, jk) * zkappa_i(ji, jk)
            zindterm(ji, jm) = ztiold(ji, jk) + zeta_i(ji, jk) * zradab_i(ji, jk)
          END DO
        END DO
        jm = nlay_s + nlay_i + 1
        DO ji = 1, npti
          ztrid(ji, jm, 1) = - zeta_i(ji, nlay_i) * zkappa_i(ji, nlay_i - 1)
          ztrid(ji, jm, 2) = 1._wp + zeta_i(ji, nlay_i) * (zkappa_i(ji, nlay_i - 1) + zkappa_i(ji, nlay_i) * zg1)
          ztrid(ji, jm, 3) = 0._wp
          zindterm(ji, jm) = ztiold(ji, nlay_i) + zeta_i(ji, nlay_i) * (zradab_i(ji, nlay_i) + zkappa_i(ji, nlay_i) * zg1 * t_bo_1d(ji))
        END DO
        DO ji = 1, npti
          IF (h_s_1d(ji) > 0._wp) THEN
            DO jm = 3, nlay_s + 1
              jk = jm - 1
              ztrid(ji, jm, 1) = - zeta_s(ji, jk) * zkappa_s(ji, jk - 1)
              ztrid(ji, jm, 2) = 1._wp + zeta_s(ji, jk) * (zkappa_s(ji, jk - 1) + zkappa_s(ji, jk))
              ztrid(ji, jm, 3) = - zeta_s(ji, jk) * zkappa_s(ji, jk)
              zindterm(ji, jm) = ztsold(ji, jk) + zeta_s(ji, jk) * zradab_s(ji, jk)
            END DO
            IF (nlay_i == 1) THEN
              ztrid(ji, nlay_s + 2, 3) = 0._wp
              zindterm(ji, nlay_s + 2) = zindterm(ji, nlay_s + 2) + zeta_i(ji, 1) * zkappa_i(ji, 1) * t_bo_1d(ji)
            END IF
            jm_min(ji) = 2
            jm_max(ji) = nlay_i + nlay_s + 1
            ztrid(ji, 2, 1) = 0._wp
            ztrid(ji, 2, 2) = 1._wp + zeta_s(ji, 1) * zkappa_s(ji, 1)
            ztrid(ji, 2, 3) = - zeta_s(ji, 1) * zkappa_s(ji, 1)
            zindterm(ji, 2) = ztsold(ji, 1) + zeta_s(ji, 1) * (zradab_s(ji, 1) + qcn_ice_1d(ji))
          ELSE
            jm_min(ji) = nlay_s + 2
            jm_max(ji) = nlay_i + nlay_s + 1
            ztrid(ji, jm_min(ji), 1) = 0._wp
            ztrid(ji, jm_min(ji), 2) = 1._wp + zeta_i(ji, 1) * zkappa_i(ji, 1)
            ztrid(ji, jm_min(ji), 3) = - zeta_i(ji, 1) * zkappa_i(ji, 1)
            zindterm(ji, jm_min(ji)) = ztiold(ji, 1) + zeta_i(ji, 1) * (zradab_i(ji, 1) + qcn_ice_1d(ji))
            IF (nlay_i == 1) THEN
              ztrid(ji, jm_min(ji), 1) = 0._wp
              ztrid(ji, jm_min(ji), 2) = 1._wp + zeta_i(ji, 1) * zkappa_i(ji, 1)
              ztrid(ji, jm_min(ji), 3) = 0._wp
              zindterm(ji, jm_min(ji)) = ztiold(ji, 1) + zeta_i(ji, 1) * (zradab_i(ji, 1) + zkappa_i(ji, 1) * t_bo_1d(ji) + qcn_ice_1d(ji))
            END IF
          END IF
          zindtbis(ji, jm_min(ji)) = zindterm(ji, jm_min(ji))
          zdiagbis(ji, jm_min(ji)) = ztrid(ji, jm_min(ji), 2)
        END DO
        jm_maxt = 0
        jm_mint = nlay_i + 5
        DO ji = 1, npti
          jm_mint = MIN(jm_min(ji), jm_mint)
          jm_maxt = MAX(jm_max(ji), jm_maxt)
        END DO
        DO jk = jm_mint + 1, jm_maxt
          DO ji = 1, npti
            jm = MIN(MAX(jm_min(ji) + 1, jk), jm_max(ji))
            zdiagbis(ji, jm) = ztrid(ji, jm, 2) - ztrid(ji, jm, 1) * ztrid(ji, jm - 1, 3) / zdiagbis(ji, jm - 1)
            zindtbis(ji, jm) = zindterm(ji, jm) - ztrid(ji, jm, 1) * zindtbis(ji, jm - 1) / zdiagbis(ji, jm - 1)
          END DO
        END DO
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) t_i_1d(ji, nlay_i) = zindtbis(ji, jm_max(ji)) / zdiagbis(ji, jm_max(ji))
        END DO
        DO jm = nlay_i + nlay_s, nlay_s + 2, - 1
          DO ji = 1, npti
            IF (.NOT. l_T_converged(ji)) THEN
              jk = jm - nlay_s - 1
              t_i_1d(ji, jk) = (zindtbis(ji, jm) - ztrid(ji, jm, 3) * t_i_1d(ji, jk + 1)) / zdiagbis(ji, jm)
            END IF
          END DO
        END DO
        DO ji = 1, npti
          IF (.NOT. l_T_converged(ji)) THEN
            IF (h_s_1d(ji) > 0._wp) THEN
              t_s_1d(ji, nlay_s) = (zindtbis(ji, nlay_s + 1) - ztrid(ji, nlay_s + 1, 3) * t_i_1d(ji, 1)) / zdiagbis(ji, nlay_s + 1)
            END IF
          END IF
        END DO
        DO ji = 1, npti
          zdti_max = 0._wp
          IF (.NOT. l_T_converged(ji)) THEN
            t_s_1d(ji, 1 : nlay_s) = MAX(MIN(t_s_1d(ji, 1 : nlay_s), rt0), rt0 - 100._wp)
            zdti_max = MAX(zdti_max, MAXVAL(ABS(t_s_1d(ji, 1 : nlay_s) - ztsb(ji, 1 : nlay_s))))
            DO jk = 1, nlay_i
              ztmelts = - rTmlt * sz_i_1d(ji, jk) + rt0
              t_i_1d(ji, jk) = MAX(MIN(t_i_1d(ji, jk), ztmelts), rt0 - 100._wp)
              zdti_max = MAX(zdti_max, ABS(t_i_1d(ji, jk) - ztib(ji, jk)))
            END DO
            IF (zdti_max < zdti_bnd) l_T_converged(ji) = .TRUE.
          END IF
        END DO
      END IF
    END DO
    ! !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_thd_zdf_bl99', 'r1', 0, 0)
    IF (ln_icectl .AND. lwp) THEN
      WRITE(numout, *) ' zdti_max : ', zdti_max
      WRITE(numout, *) ' iconv    : ', iconv
    END IF
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    DO ji = 1, npti
      qcn_ice_top_1d(ji) = - isnow(ji) * zkappa_s(ji, 0) * zg1s * (t_s_1d(ji, 1) - t_su_1d(ji)) - (1._wp - isnow(ji)) * zkappa_i(ji, 0) * zg1 * (t_i_1d(ji, 1) - t_su_1d(ji))
      qcn_ice_bot_1d(ji) = - zkappa_i(ji, nlay_i) * zg1 * (t_bo_1d(ji) - t_i_1d(ji, nlay_i))
    END DO
    IF (k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU) THEN
      DO ji = 1, npti
        hfx_err_dif_1d(ji) = hfx_err_dif_1d(ji) - (qns_ice_1d(ji) - zqns_ice_b(ji)) * a_i_1d(ji)
      END DO
    ELSE IF (k_cnd == np_cnd_ON) THEN
      DO ji = 1, npti
        hfx_err_dif_1d(ji) = hfx_err_dif_1d(ji) - (qcn_ice_top_1d(ji) - qcn_ice_1d(ji)) * a_i_1d(ji)
      END DO
    END IF
    !$ACC END KERNELS
    !CALL profile_psy_data2 % PreStart('ice_thd_zdf_bl99', 'r2', 0, 0)
    IF (k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_ON) THEN
      CALL ice_var_enthalpy
      !$ACC KERNELS ! CDe
      !$ACC loop private(zhfx_err) ! CDe
      DO ji = 1, npti
        zdq = - zq_ini(ji) + (SUM(e_i_1d(ji, 1 : nlay_i)) * h_i_1d(ji) * r1_nlay_i + SUM(e_s_1d(ji, 1 : nlay_s)) * h_s_1d(ji) * r1_nlay_s)
        IF (k_cnd == np_cnd_OFF) THEN
          IF (t_su_1d(ji) < rt0) THEN
            zhfx_err = (qns_ice_1d(ji) + qsr_ice_1d(ji) - zradtr_i(ji, nlay_i) - qcn_ice_bot_1d(ji) + zdq * r1_rdtice) * a_i_1d(ji)
          ELSE
            zhfx_err = (qcn_ice_top_1d(ji) + qtr_ice_top_1d(ji) - zradtr_i(ji, nlay_i) - qcn_ice_bot_1d(ji) + zdq * r1_rdtice) * a_i_1d(ji)
          END IF
        ELSE IF (k_cnd == np_cnd_ON) THEN
          zhfx_err = (qcn_ice_top_1d(ji) + qtr_ice_top_1d(ji) - zradtr_i(ji, nlay_i) - qcn_ice_bot_1d(ji) + zdq * r1_rdtice) * a_i_1d(ji)
        END IF
        hfx_err_dif_1d(ji) = hfx_err_dif_1d(ji) + zhfx_err
        hfx_dif_1d(ji) = hfx_dif_1d(ji) - zdq * r1_rdtice * a_i_1d(ji)
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS ! CDe
    DO ji = 1, npti
      IF (h_s_1d(ji) > 0.1_wp) THEN
        cnd_ice_1d(ji) = 2._wp * zkappa_s(ji, 0)
      ELSE
        IF (h_i_1d(ji) > 0.1_wp) THEN
          cnd_ice_1d(ji) = 2._wp * zkappa_i(ji, 0)
        ELSE
          cnd_ice_1d(ji) = 2._wp * ztcond_i(ji, 0) * 10._wp
        END IF
      END IF
      t1_ice_1d(ji) = isnow(ji) * t_s_1d(ji, 1) + (1._wp - isnow(ji)) * t_i_1d(ji, 1)
    END DO
    !$ACC END KERNELS 
    !CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    IF (k_cnd == np_cnd_EMU) THEN
      t_s_1d(1 : npti, :) = ztsold(1 : npti, :)
      t_i_1d(1 : npti, :) = ztiold(1 : npti, :)
      qcn_ice_1d(1 : npti) = qcn_ice_top_1d(1 : npti)
    END IF
    !CALL profile_psy_data3 % PreStart('ice_thd_zdf_bl99', 'r3', 0, 0)
    DO ji = 1, npti
      zfac = rn_cnd_s * zh_i(ji) + ztcond_i(ji, 1) * zh_s(ji)
      IF (h_s_1d(ji) >= zhs_min) THEN
        t_si_1d(ji) = (rn_cnd_s * zh_i(ji) * t_s_1d(ji, 1) + ztcond_i(ji, 1) * zh_s(ji) * t_i_1d(ji, 1)) / MAX(epsi10, zfac)
      ELSE
        t_si_1d(ji) = t_su_1d(ji)
      END IF
    END DO
    !$ACC END KERNELS
    !CALL profile_psy_data3 % PostEnd
  END SUBROUTINE ice_thd_zdf_BL99
END MODULE icethd_zdf_BL99
