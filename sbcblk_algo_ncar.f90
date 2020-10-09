MODULE sbcblk_algo_ncar
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE sbcwave, ONLY: cdn_wave
  USE iom
  USE lib_mpp
  USE in_out_manager
  USE prtctl
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TURB_NCAR
  REAL(KIND = wp), PARAMETER :: rctv0 = 0.608
  CONTAINS
  SUBROUTINE turb_ncar(zt, zu, sst, t_zt, ssq, q_zt, U_zu, Cd, Ch, Ce, t_zu, q_zu, U_blk, Cdn, Chn, Cen)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: zt
    REAL(KIND = wp), INTENT(IN) :: zu
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: sst
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: t_zt
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: ssq
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: q_zt
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: U_zu
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: Cd
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: Ch
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: Ce
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: t_zu
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: q_zu
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: U_blk
    REAL(KIND = wp), INTENT(OUT), DIMENSION(jpi, jpj) :: Cdn, Chn, Cen
    INTEGER :: j_itt
    LOGICAL :: l_zt_equal_zu = .FALSE.
    INTEGER, PARAMETER :: nb_itt = 4
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: Cx_n10
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: sqrt_Cd_n10
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zeta_u
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zpsi_h_u
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmp0, ztmp1, ztmp2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: stab
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    CALL profile_psy_data0 % PreStart('turb_ncar', 'r0', 0, 0)
    l_zt_equal_zu = .FALSE.
    IF (ABS(zu - zt) < 0.01) l_zt_equal_zu = .TRUE.
    U_blk = MAX(0.5, U_zu)
    ztmp0 = t_zt * (1. + rctv0 * q_zt) - sst * (1. + rctv0 * ssq)
    stab = 0.5 + SIGN(0.5, ztmp0)
    CALL profile_psy_data0 % PostEnd
    IF (ln_cdgw) THEN
      !$ACC KERNELS
      cdn_wave(:, :) = cdn_wave(:, :) + rsmall * (1._wp - tmask(:, :, 1))
      ztmp0(:, :) = cdn_wave(:, :)
      !$ACC END KERNELS
    ELSE
      CALL profile_psy_data1 % PreStart('turb_ncar', 'r1', 0, 0)
      ztmp0 = cd_neutral_10m(U_blk)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('turb_ncar', 'r2', 0, 0)
    sqrt_Cd_n10 = SQRT(ztmp0)
    Cd = ztmp0
    Ce = 1.E-3 * (34.6 * sqrt_Cd_n10)
    Ch = 1.E-3 * sqrt_Cd_n10 * (18. * stab + 32.7 * (1. - stab))
    stab = sqrt_Cd_n10
    IF (ln_cdgw) cen = ce
    chn = ch
    t_zu = t_zt
    q_zu = q_zt
    CALL profile_psy_data2 % PostEnd
    DO j_itt = 1, nb_itt
      CALL profile_psy_data3 % PreStart('turb_ncar', 'r3', 0, 0)
      ztmp1 = t_zu - sst
      ztmp2 = q_zu - ssq
      ztmp1 = Ch / stab * ztmp1
      ztmp2 = Ce / stab * ztmp2
      ztmp0 = 1. + rctv0 * q_zu
      ztmp0 = (grav * vkarmn / (t_zu * ztmp0) * (ztmp1 * ztmp0 + rctv0 * t_zu * ztmp2)) / (Cd * U_blk * U_blk)
      zeta_u = zu * ztmp0
      zeta_u = SIGN(MIN(ABS(zeta_u), 10.0), zeta_u)
      zpsi_h_u = psi_h(zeta_u)
      IF (.NOT. l_zt_equal_zu) THEN
        stab = zt * ztmp0
        stab = SIGN(MIN(ABS(stab), 10.0), stab)
        stab = LOG(zt / zu) + zpsi_h_u - psi_h(stab)
        t_zu = t_zt - ztmp1 / vkarmn * stab
        q_zu = q_zt - ztmp2 / vkarmn * stab
        q_zu = MAX(0., q_zu)
      END IF
      ztmp2 = psi_m(zeta_u)
      CALL profile_psy_data3 % PostEnd
      IF (ln_cdgw) THEN
        CALL profile_psy_data4 % PreStart('turb_ncar', 'r4', 0, 0)
        stab = vkarmn / (vkarmn / sqrt_Cd_n10 - ztmp2)
        Cd = stab * stab
        ztmp0 = (LOG(zu / 10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
        ztmp2 = stab / sqrt_Cd_n10
        ztmp1 = 1. + Chn * ztmp0
        Ch = Chn * ztmp2 / ztmp1
        ztmp1 = 1. + Cen * ztmp0
        Ce = Cen * ztmp2 / ztmp1
        CALL profile_psy_data4 % PostEnd
      ELSE
        CALL profile_psy_data5 % PreStart('turb_ncar', 'r5', 0, 0)
        ztmp0 = MAX(0.25, U_blk / (1. + sqrt_Cd_n10 / vkarmn * (LOG(zu / 10.) - ztmp2)))
        ztmp0 = cd_neutral_10m(ztmp0)
        CALL profile_psy_data5 % PostEnd
        !$ACC KERNELS
        Cdn(:, :) = ztmp0
        sqrt_Cd_n10 = SQRT(ztmp0)
        stab = 0.5 + SIGN(0.5, zeta_u)
        Cx_n10 = 1.E-3 * sqrt_Cd_n10 * (18. * stab + 32.7 * (1. - stab))
        Chn(:, :) = Cx_n10
        ztmp1 = 1. + sqrt_Cd_n10 / vkarmn * (LOG(zu / 10.) - ztmp2)
        Cd = ztmp0 / (ztmp1 * ztmp1)
        stab = SQRT(Cd)
        ztmp0 = (LOG(zu / 10.) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
        ztmp2 = stab / sqrt_Cd_n10
        ztmp1 = 1. + Cx_n10 * ztmp0
        Ch = Cx_n10 * ztmp2 / ztmp1
        Cx_n10 = 1.E-3 * (34.6 * sqrt_Cd_n10)
        Cen(:, :) = Cx_n10
        ztmp1 = 1. + Cx_n10 * ztmp0
        Ce = Cx_n10 * ztmp2 / ztmp1
        !$ACC END KERNELS
      END IF
    END DO
  END SUBROUTINE turb_ncar
  FUNCTION cd_neutral_10m(pw10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pw10
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: cd_neutral_10m
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zgt33, zw, zw6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('cd_neutral_10m', 'r0', 0, 0)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zw = pw10(ji, jj)
        zw6 = zw * zw * zw
        zw6 = zw6 * zw6
        zgt33 = 0.5 + SIGN(0.5, (zw - 33.))
        cd_neutral_10m(ji, jj) = 1.E-3 * ((1. - zgt33) * (2.7 / zw + 0.142 + zw / 13.09 - 3.14807E-10 * zw6) + zgt33 * 2.34)
        cd_neutral_10m(ji, jj) = MAX(cd_neutral_10m(ji, jj), 1.E-6)
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END FUNCTION cd_neutral_10m
  FUNCTION psi_m(pzeta)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pzeta
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: psi_m
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zx2, zx, zstab
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zx2 = SQRT(ABS(1. - 16. * pzeta(ji, jj)))
        zx2 = MAX(zx2, 1.)
        zx = SQRT(zx2)
        zstab = 0.5 + SIGN(0.5, pzeta(ji, jj))
        psi_m(ji, jj) = zstab * (- 5. * pzeta(ji, jj)) + (1. - zstab) * (2. * LOG((1. + zx) * 0.5) + LOG((1. + zx2) * 0.5) - 2. * &
&ATAN(zx) + rpi * 0.5)
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION psi_m
  FUNCTION psi_h(pzeta)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pzeta
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: psi_h
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zx2, zstab
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zx2 = SQRT(ABS(1. - 16. * pzeta(ji, jj)))
        zx2 = MAX(zx2, 1.)
        zstab = 0.5 + SIGN(0.5, pzeta(ji, jj))
        psi_h(ji, jj) = zstab * (- 5. * pzeta(ji, jj)) + (1. - zstab) * (2. * LOG((1. + zx2) * 0.5))
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION psi_h
END MODULE sbcblk_algo_ncar