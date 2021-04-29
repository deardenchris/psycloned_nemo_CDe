MODULE sbcblk_algo_ecmwf
  USE oce
  USE dom_oce
  USE phycst
  USE iom
  USE lib_mpp
  USE in_out_manager
  USE prtctl
  USE sbcwave, ONLY: cdn_wave
  USE lib_fortran
  USE sbc_oce
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: TURB_ECMWF
  REAL(KIND = wp), PARAMETER :: charn0 = 0.018
  REAL(KIND = wp), PARAMETER :: zi0 = 1000.
  REAL(KIND = wp), PARAMETER :: Beta0 = 1.
  REAL(KIND = wp), PARAMETER :: rctv0 = 0.608
  REAL(KIND = wp), PARAMETER :: Cp_dry = 1005.0
  REAL(KIND = wp), PARAMETER :: Cp_vap = 1860.0
  REAL(KIND = wp), PARAMETER :: alpha_M = 0.11
  REAL(KIND = wp), PARAMETER :: alpha_H = 0.40
  REAL(KIND = wp), PARAMETER :: alpha_Q = 0.62
  CONTAINS
  SUBROUTINE TURB_ECMWF(zt, zu, sst, t_zt, ssq, q_zt, U_zu, Cd, Ch, Ce, t_zu, q_zu, U_blk, Cdn, Chn, Cen)
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
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: u_star, t_star, q_star, dt_zu, dq_zu, znu_a, Linv, z0, z0t, z0q
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: func_m, func_h
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmp0, ztmp1, ztmp2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('turb_ecmwf', 'r0', 0, 0)
    l_zt_equal_zu = .FALSE.
    IF (ABS(zu - zt) < 0.01) l_zt_equal_zu = .TRUE.
    t_zu = MAX(t_zt, 0.0)
    q_zu = MAX(q_zt, 1.E-6)
    dt_zu = t_zu - sst
    dt_zu = SIGN(MAX(ABS(dt_zu), 1.E-6), dt_zu)
    dq_zu = q_zu - ssq
    dq_zu = SIGN(MAX(ABS(dq_zu), 1.E-9), dq_zu)
    znu_a = visc_air(t_zt)
    ztmp2 = 0.5 * 0.5
    U_blk = SQRT(U_zu * U_zu + ztmp2)
    ztmp2 = 10000.
    ztmp0 = LOG(zu * ztmp2)
    ztmp1 = LOG(10. * ztmp2)
    u_star = 0.035 * U_blk * ztmp1 / ztmp0
    z0 = charn0 * u_star * u_star / grav + 0.11 * znu_a / u_star
    z0t = 0.1 * EXP(vkarmn / (0.00115 / (vkarmn / ztmp1)))
    Cd = (vkarmn / ztmp0) ** 2
    ztmp0 = vkarmn * vkarmn / LOG(zt * z0t) / Cd
    ztmp2 = Ri_bulk(zu, t_zu, dt_zu, q_zu, dq_zu, U_blk)
    ztmp1 = 0.5 + SIGN(0.5, ztmp2)
    func_m = ztmp0 * ztmp2
    func_h = (1. - ztmp1) * (func_m / (1. + ztmp2 / (- zu / (zi0 * 0.004 * Beta0 ** 3)))) + ztmp1 * (func_m * (1. + 27. / 9. * ztmp2 / ztmp0))
    ztmp0 = vkarmn / (LOG(zu * z0t) - psi_h_ecmwf(func_h))
    u_star = U_blk * vkarmn / (LOG(zu) - LOG(z0) - psi_m_ecmwf(func_h))
    t_star = dt_zu * ztmp0
    q_star = dq_zu * ztmp0
    IF (.NOT. l_zt_equal_zu) THEN
      ztmp0 = psi_h_ecmwf(func_h) - psi_h_ecmwf(zt * func_h / zu)
      ztmp1 = LOG(zt / zu) + ztmp0
      t_zu = t_zt - t_star / vkarmn * ztmp1
      q_zu = q_zt - q_star / vkarmn * ztmp1
      q_zu = (0.5 + SIGN(0.5, q_zu)) * q_zu
      dt_zu = t_zu - sst
      dt_zu = SIGN(MAX(ABS(dt_zu), 1.E-6), dt_zu)
      dq_zu = q_zu - ssq
      dq_zu = SIGN(MAX(ABS(dq_zu), 1.E-9), dq_zu)
    END IF
    ztmp0 = (1. + rctv0 * q_zu)
    Linv = grav * vkarmn * (t_star * ztmp0 + rctv0 * t_zu * q_star) / (u_star * u_star * t_zu * ztmp0)
    ztmp1 = zu + z0
    ztmp0 = ztmp1 * Linv
    func_m = LOG(ztmp1) - LOG(z0) - psi_m_ecmwf(ztmp0) + psi_m_ecmwf(z0 * Linv)
    func_h = LOG(ztmp1 * z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(1. / z0t * Linv)
    DO j_itt = 1, nb_itt
      ztmp0 = Ri_bulk(zu, t_zu, dt_zu, q_zu, dq_zu, U_blk)
      Linv = ztmp0 * func_m * func_m / func_h / zu
      ztmp1 = zu + z0
      func_m = LOG(ztmp1) - LOG(z0) - psi_m_ecmwf(ztmp1 * Linv) + psi_m_ecmwf(z0 * Linv)
      u_star = U_blk * vkarmn / func_m
      ztmp2 = u_star * u_star
      ztmp1 = znu_a / u_star
      z0 = alpha_M * ztmp1 + charn0 * ztmp2 / grav
      z0t = alpha_H * ztmp1
      z0q = alpha_Q * ztmp1
      ztmp2 = ztmp2 * (MAX(- zi0 * Linv / vkarmn, 0.)) ** (2. / 3.)
      U_blk = MAX(SQRT(U_zu * U_zu + ztmp2), 0.2)
      IF (.NOT. l_zt_equal_zu) THEN
        func_h = psi_h_ecmwf((zu + z0) * Linv)
        func_m = psi_h_ecmwf((zt + z0) * Linv)
        ztmp2 = psi_h_ecmwf(z0t * Linv)
        ztmp0 = func_h - ztmp2
        ztmp1 = vkarmn / (LOG(zu + z0) - LOG(z0t) - ztmp0)
        t_star = dt_zu * ztmp1
        ztmp2 = ztmp0 - func_m + ztmp2
        ztmp1 = LOG(zt / zu) + ztmp2
        t_zu = t_zt - t_star / vkarmn * ztmp1
        ztmp2 = psi_h_ecmwf(z0q * Linv)
        ztmp0 = func_h - ztmp2
        ztmp1 = vkarmn / (LOG(zu + z0) - LOG(z0q) - ztmp0)
        q_star = dq_zu * ztmp1
        ztmp2 = ztmp0 - func_m + ztmp2
        ztmp1 = LOG(zt / zu) + ztmp2
        q_zu = q_zt - q_star / vkarmn * ztmp1
        dt_zu = t_zu - sst
        dt_zu = SIGN(MAX(ABS(dt_zu), 1.E-6), dt_zu)
        dq_zu = q_zu - ssq
        dq_zu = SIGN(MAX(ABS(dq_zu), 1.E-9), dq_zu)
      END IF
      ztmp1 = zu + z0
      ztmp0 = ztmp1 * Linv
      func_m = LOG(ztmp1) - LOG(z0) - psi_m_ecmwf(ztmp0) + psi_m_ecmwf(z0 * Linv)
      func_h = LOG(ztmp1) - LOG(z0t) - psi_h_ecmwf(ztmp0) + psi_h_ecmwf(z0t * Linv)
    END DO
    Cd = vkarmn * vkarmn / (func_m * func_m)
    Ch = vkarmn * vkarmn / (func_m * func_h)
    ztmp1 = LOG((zu + z0) / z0q) - psi_h_ecmwf((zu + z0) * Linv) + psi_h_ecmwf(z0q * Linv)
    Ce = vkarmn * vkarmn / (func_m * ztmp1)
    ztmp1 = zu + z0
    Cdn = vkarmn * vkarmn / (LOG(ztmp1 / z0) * LOG(ztmp1 / z0))
    Chn = vkarmn * vkarmn / (LOG(ztmp1 / z0t) * LOG(ztmp1 / z0t))
    Cen = vkarmn * vkarmn / (LOG(ztmp1 / z0q) * LOG(ztmp1 / z0q))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE TURB_ECMWF
  FUNCTION psi_m_ecmwf(pzeta)
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: psi_m_ecmwf
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pzeta
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zzeta, zx, ztmp, psi_unst, psi_stab, stab
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zzeta = MIN(pzeta(ji, jj), 5.)
        zx = SQRT(ABS(1. - 16. * zzeta))
        ztmp = 1. + SQRT(zx)
        ztmp = ztmp * ztmp
        psi_unst = LOG(0.125 * ztmp * (1. + zx)) - 2. * ATAN(SQRT(zx)) + 0.5 * rpi
        psi_stab = - 2. / 3. * (zzeta - 5. / 0.35) * EXP(- 0.35 * zzeta) - zzeta - 2. / 3. * 5. / 0.35
        stab = 0.5 + SIGN(0.5, zzeta)
        psi_m_ecmwf(ji, jj) = (1. - stab) * psi_unst + stab * psi_stab
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION psi_m_ecmwf
  FUNCTION psi_h_ecmwf(pzeta)
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: psi_h_ecmwf
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pzeta
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zzeta, zx, psi_unst, psi_stab, stab
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zzeta = MIN(pzeta(ji, jj), 5.)
        zx = ABS(1. - 16. * zzeta) ** .25
        psi_unst = 2. * LOG(0.5 * (1. + zx * zx))
        psi_stab = - 2. / 3. * (zzeta - 5. / 0.35) * EXP(- 0.35 * zzeta) - ABS(1. + 2. / 3. * zzeta) ** 1.5 - 2. / 3. * 5. / 0.35 + 1.
        stab = 0.5 + SIGN(0.5, zzeta)
        psi_h_ecmwf(ji, jj) = (1. - stab) * psi_unst + stab * psi_stab
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION psi_h_ecmwf
  FUNCTION Ri_bulk(pz, ptz, pdt, pqz, pdq, pub)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: Ri_bulk
    REAL(KIND = wp), INTENT(IN) :: pz
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptz
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqz
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pdq
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pub
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ri_bulk', 'r0', 0, 0)
    Ri_bulk = grav * pz / (pub * pub) * (pdt / (ptz - 0.5_wp * (pdt + grav * pz / (Cp_dry + Cp_vap * pqz))) + rctv0 * pdq)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION Ri_bulk
  FUNCTION visc_air(ptak)
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: visc_air
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptak
    INTEGER :: ji, jj
    REAL(KIND = wp) :: ztc, ztc2
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        ztc = ptak(ji, jj) - rt0
        ztc2 = ztc * ztc
        visc_air(ji, jj) = 1.326E-5 * (1. + 6.542E-3 * ztc + 8.301E-6 * ztc2 - 4.84E-9 * ztc2 * ztc)
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION visc_air
END MODULE sbcblk_algo_ecmwf