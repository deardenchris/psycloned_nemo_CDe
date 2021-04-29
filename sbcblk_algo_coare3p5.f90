MODULE sbcblk_algo_coare3p5
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
  PUBLIC :: TURB_COARE3P5
  REAL(KIND = wp), PARAMETER :: charn0_max = 0.028
  REAL(KIND = wp), PARAMETER :: zi0 = 600.
  REAL(KIND = wp), PARAMETER :: Beta0 = 1.25
  REAL(KIND = wp), PARAMETER :: rctv0 = 0.608
  CONTAINS
  SUBROUTINE turb_coare3p5(zt, zu, sst, t_zt, ssq, q_zt, U_zu, Cd, Ch, Ce, t_zu, q_zu, U_blk, Cdn, Chn, Cen)
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
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: u_star, t_star, q_star, dt_zu, dq_zu, znu_a, z0, z0t
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zeta_u
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmp0, ztmp1, ztmp2
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zeta_t
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('turb_coare3p5', 'r0', 0, 0)
    l_zt_equal_zu = .FALSE.
    IF (ABS(zu - zt) < 0.01) l_zt_equal_zu = .TRUE.
    IF (.NOT. l_zt_equal_zu) ALLOCATE(zeta_t(jpi, jpj))
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
    ztmp2 = MIN(0.0017 * U_zu - 0.005, charn0_max)
    ztmp2 = MAX(ztmp2, 0.)
    z0 = ztmp2 * u_star * u_star / grav + 0.11 * znu_a / u_star
    z0t = 0.1 * EXP(vkarmn / (0.00115 / (vkarmn / ztmp1)))
    ztmp2 = vkarmn / ztmp0
    Cd = ztmp2 * ztmp2
    ztmp0 = vkarmn * vkarmn / LOG(zt * z0t) / Cd
    ztmp2 = grav * zu * (dt_zu + rctv0 * t_zu * dq_zu) / (t_zu * U_blk * U_blk)
    ztmp1 = 0.5 + SIGN(0.5, ztmp2)
    ztmp0 = ztmp0 * ztmp2
    zeta_u = (1. - ztmp1) * (ztmp0 / (1. + ztmp2 / (- zu / (zi0 * 0.004 * Beta0 ** 3)))) + ztmp1 * (ztmp0 * (1. + 27. / 9. * ztmp2 / ztmp0))
    ztmp0 = vkarmn / (LOG(zu * z0t) - psi_h_coare(zeta_u))
    u_star = U_blk * vkarmn / (LOG(zu) - LOG(z0) - psi_m_coare(zeta_u))
    t_star = dt_zu * ztmp0
    q_star = dq_zu * ztmp0
    IF (.NOT. l_zt_equal_zu) THEN
      zeta_t = zt * zeta_u / zu
      ztmp0 = psi_h_coare(zeta_u) - psi_h_coare(zeta_t)
      ztmp1 = LOG(zt / zu) + ztmp0
      t_zu = t_zt - t_star / vkarmn * ztmp1
      q_zu = q_zt - q_star / vkarmn * ztmp1
      q_zu = (0.5 + SIGN(0.5, q_zu)) * q_zu
      dt_zu = t_zu - sst
      dt_zu = SIGN(MAX(ABS(dt_zu), 1.E-6), dt_zu)
      dq_zu = q_zu - ssq
      dq_zu = SIGN(MAX(ABS(dq_zu), 1.E-9), dq_zu)
    END IF
    DO j_itt = 1, nb_itt
      ztmp0 = One_on_L(t_zu, q_zu, u_star, t_star, q_star)
      ztmp1 = u_star * u_star
      ztmp2 = Beta0 * Beta0 * ztmp1 * (MAX(- zi0 * ztmp0 / vkarmn, 0.)) ** (2. / 3.)
      U_blk = MAX(SQRT(U_zu * U_zu + ztmp2), 0.2)
      ztmp2 = u_star / vkarmn * LOG(10. / z0)
      ztmp2 = MIN(0.0017 * ztmp2 - 0.005, charn0_max)
      ztmp2 = MAX(ztmp2, 0.)
      z0 = ztmp2 * ztmp1 / grav + 0.11 * znu_a / u_star
      ztmp1 = z0 * u_star / znu_a
      z0t = MIN(1.6E-4, 5.8E-5 * ztmp1 ** (- 0.72))
      zeta_u = zu * ztmp0
      zeta_u = SIGN(MIN(ABS(zeta_u), 50.0), zeta_u)
      IF (.NOT. l_zt_equal_zu) THEN
        zeta_t = zt * ztmp0
        zeta_t = SIGN(MIN(ABS(zeta_t), 50.0), zeta_t)
      END IF
      ztmp0 = psi_h_coare(zeta_u)
      ztmp1 = vkarmn / (LOG(zu) - LOG(z0t) - ztmp0)
      t_star = dt_zu * ztmp1
      q_star = dq_zu * ztmp1
      u_star = U_blk * vkarmn / (LOG(zu) - LOG(z0) - psi_m_coare(zeta_u))
      IF (.NOT. l_zt_equal_zu) THEN
        ztmp2 = ztmp0 - psi_h_coare(zeta_t)
        ztmp1 = LOG(zt / zu) + ztmp2
        t_zu = t_zt - t_star / vkarmn * ztmp1
        q_zu = q_zt - q_star / vkarmn * ztmp1
        dt_zu = t_zu - sst
        dt_zu = SIGN(MAX(ABS(dt_zu), 1.E-6), dt_zu)
        dq_zu = q_zu - ssq
        dq_zu = SIGN(MAX(ABS(dq_zu), 1.E-9), dq_zu)
      END IF
    END DO
    ztmp0 = u_star / U_blk
    Cd = ztmp0 * ztmp0
    Ch = ztmp0 * t_star / dt_zu
    Ce = ztmp0 * q_star / dq_zu
    ztmp1 = zu + z0
    Cdn = vkarmn * vkarmn / (LOG(ztmp1 / z0) * LOG(ztmp1 / z0))
    Chn = vkarmn * vkarmn / (LOG(ztmp1 / z0t) * LOG(ztmp1 / z0t))
    Cen = Chn
    IF (.NOT. l_zt_equal_zu) DEALLOCATE(zeta_t)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE turb_coare3p5
  FUNCTION One_on_L(ptha, pqa, pus, pts, pqs)
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: One_on_L
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptha, pqa, pus, pts, pqs
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zqa
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zqa = (1. + rctv0 * pqa(ji, jj))
        One_on_L(ji, jj) = grav * vkarmn * (pts(ji, jj) * zqa + rctv0 * ptha(ji, jj) * pqs(ji, jj)) / (pus(ji, jj) * pus(ji, jj) * ptha(ji, jj) * zqa)
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION One_on_L
  FUNCTION psi_m_coare(pzeta)
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: psi_m_coare
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pzeta
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zta, zphi_m, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zta = pzeta(ji, jj)
        zphi_m = ABS(1. - 15. * zta) ** .25
        zpsi_k = 2. * LOG((1. + zphi_m) / 2.) + LOG((1. + zphi_m * zphi_m) / 2.) - 2. * ATAN(zphi_m) + 0.5 * rpi
        zphi_c = ABS(1. - 10.15 * zta) ** .3333
        zpsi_c = 1.5 * LOG((1. + zphi_c + zphi_c * zphi_c) / 3.) - 1.7320508 * ATAN((1. + 2. * zphi_c) / 1.7320508) + 1.813799447
        zf = zta * zta
        zf = zf / (1. + zf)
        zc = MIN(50., 0.35 * zta)
        zstab = 0.5 + SIGN(0.5, zta)
        psi_m_coare(ji, jj) = (1. - zstab) * ((1. - zf) * zpsi_k + zf * zpsi_c) - zstab * (1. + 1. * zta + 0.6667 * (zta - 14.28) / EXP(zc) + 8.525)
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION psi_m_coare
  FUNCTION psi_h_coare(pzeta)
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: psi_h_coare
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pzeta
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zta, zphi_h, zphi_c, zpsi_k, zpsi_c, zf, zc, zstab
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zta = pzeta(ji, jj)
        zphi_h = (ABS(1. - 15. * zta)) ** .5
        zpsi_k = 2. * LOG((1. + zphi_h) / 2.)
        zphi_c = (ABS(1. - 34.15 * zta)) ** .3333
        zpsi_c = 1.5 * LOG((1. + zphi_c + zphi_c * zphi_c) / 3.) - 1.7320508 * ATAN((1. + 2. * zphi_c) / 1.7320508) + 1.813799447
        zf = zta * zta
        zf = zf / (1. + zf)
        zc = MIN(50., 0.35 * zta)
        zstab = 0.5 + SIGN(0.5, zta)
        psi_h_coare(ji, jj) = (1. - zstab) * ((1. - zf) * zpsi_k + zf * zpsi_c) - zstab * ((ABS(1. + 2. * zta / 3.)) ** 1.5 + .6667 * (zta - 14.28) / EXP(zc) + 8.525)
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION psi_h_coare
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
END MODULE sbcblk_algo_coare3p5