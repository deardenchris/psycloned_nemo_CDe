MODULE icbdyn
  USE par_oce
  USE dom_oce
  USE phycst
  USE icb_oce
  USE icbutl
  USE icbdia
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_dyn
  CONTAINS
  SUBROUTINE icb_dyn(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    LOGICAL :: ll_bounced
    REAL(KIND = wp) :: zuvel1, zvvel1, zu1, zv1, zax1, zay1, zxi1, zyj1
    REAL(KIND = wp) :: zuvel2, zvvel2, zu2, zv2, zax2, zay2, zxi2, zyj2
    REAL(KIND = wp) :: zuvel3, zvvel3, zu3, zv3, zax3, zay3, zxi3, zyj3
    REAL(KIND = wp) :: zuvel4, zvvel4, zu4, zv4, zax4, zay4, zxi4, zyj4
    REAL(KIND = wp) :: zuvel_n, zvvel_n, zxi_n, zyj_n
    REAL(KIND = wp) :: zdt, zdt_2, zdt_6, ze1, ze2
    TYPE(iceberg), POINTER :: berg
    TYPE(point), POINTER :: pt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_dyn', 'r0', 0, 0)
    zdt = berg_dt
    zdt_2 = zdt * 0.5_wp
    zdt_6 = zdt / 6._wp
    berg => first_berg
    DO WHILE (ASSOCIATED(berg))
      pt => berg % current_point
      ll_bounced = .FALSE.
      zxi1 = pt % xi
      zuvel1 = pt % uvel
      zyj1 = pt % yj
      zvvel1 = pt % vvel
      CALL icb_accel(berg, zxi1, ze1, zuvel1, zuvel1, zax1, zyj1, ze2, zvvel1, zvvel1, zay1, zdt_2)
      zu1 = zuvel1 / ze1
      zv1 = zvvel1 / ze2
      zxi2 = zxi1 + zdt_2 * zu1
      zuvel2 = zuvel1 + zdt_2 * zax1
      zyj2 = zyj1 + zdt_2 * zv1
      zvvel2 = zvvel1 + zdt_2 * zay1
      CALL icb_ground(zxi2, zxi1, zu1, zyj2, zyj1, zv1, ll_bounced)
      CALL icb_accel(berg, zxi2, ze1, zuvel2, zuvel1, zax2, zyj2, ze2, zvvel2, zvvel1, zay2, zdt_2)
      zu2 = zuvel2 / ze1
      zv2 = zvvel2 / ze2
      zxi3 = zxi1 + zdt_2 * zu2
      zuvel3 = zuvel1 + zdt_2 * zax2
      zyj3 = zyj1 + zdt_2 * zv2
      zvvel3 = zvvel1 + zdt_2 * zay2
      CALL icb_ground(zxi3, zxi1, zu3, zyj3, zyj1, zv3, ll_bounced)
      CALL icb_accel(berg, zxi3, ze1, zuvel3, zuvel1, zax3, zyj3, ze2, zvvel3, zvvel1, zay3, zdt)
      zu3 = zuvel3 / ze1
      zv3 = zvvel3 / ze2
      zxi4 = zxi1 + zdt * zu3
      zuvel4 = zuvel1 + zdt * zax3
      zyj4 = zyj1 + zdt * zv3
      zvvel4 = zvvel1 + zdt * zay3
      CALL icb_ground(zxi4, zxi1, zu4, zyj4, zyj1, zv4, ll_bounced)
      CALL icb_accel(berg, zxi4, ze1, zuvel4, zuvel1, zax4, zyj4, ze2, zvvel4, zvvel1, zay4, zdt)
      zu4 = zuvel4 / ze1
      zv4 = zvvel4 / ze2
      zxi_n = pt % xi + zdt_6 * (zu1 + 2. * (zu2 + zu3) + zu4)
      zyj_n = pt % yj + zdt_6 * (zv1 + 2. * (zv2 + zv3) + zv4)
      zuvel_n = pt % uvel + zdt_6 * (zax1 + 2. * (zax2 + zax3) + zax4)
      zvvel_n = pt % vvel + zdt_6 * (zay1 + 2. * (zay2 + zay3) + zay4)
      CALL icb_ground(zxi_n, zxi1, zuvel_n, zyj_n, zyj1, zvvel_n, ll_bounced)
      pt % uvel = zuvel_n
      pt % vvel = zvvel_n
      pt % xi = zxi_n
      pt % yj = zyj_n
      pt % lon = icb_utl_bilin_x(glamt, pt % xi, pt % yj)
      pt % lat = icb_utl_bilin(gphit, pt % xi, pt % yj, 'T')
      berg => berg % next
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_dyn
  SUBROUTINE icb_ground(pi, pi0, pu, pj, pj0, pv, ld_bounced)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: pi, pj
    REAL(KIND = wp), INTENT(IN) :: pi0, pj0
    REAL(KIND = wp), INTENT(INOUT) :: pu, pv
    LOGICAL, INTENT(OUT) :: ld_bounced
    INTEGER :: ii, ii0
    INTEGER :: ij, ij0
    INTEGER :: ibounce_method
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('icb_ground', 'r0', 0, 0)
    ld_bounced = .FALSE.
    ii0 = INT(pi0 + 0.5)
    ij0 = INT(pj0 + 0.5)
    ii = INT(pi + 0.5)
    ij = INT(pj + 0.5)
    CALL profile_psy_data0 % PostEnd
    IF (ii == ii0 .AND. ij == ij0) RETURN
    CALL profile_psy_data1 % PreStart('icb_ground', 'r1', 0, 0)
    ii0 = mi1(ii0)
    ij0 = mj1(ij0)
    ii = mi1(ii)
    ij = mj1(ij)
    CALL profile_psy_data1 % PostEnd
    IF (tmask(ii, ij, 1) /= 0._wp) RETURN
    CALL profile_psy_data2 % PreStart('icb_ground', 'r2', 0, 0)
    ld_bounced = .TRUE.
    ibounce_method = 2
    SELECT CASE (ibounce_method)
    CASE (1)
      pi = pi0
      pj = pj0
      pu = 0._wp
      pv = 0._wp
    CASE (2)
      IF (ii0 /= ii) THEN
        pi = pi0
        pu = 0._wp
      END IF
      IF (ij0 /= ij) THEN
        pj = pj0
        pv = 0._wp
      END IF
    END SELECT
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE icb_ground
  SUBROUTINE icb_accel(berg, pxi, pe1, puvel, puvel0, pax, pyj, pe2, pvvel, pvvel0, pay, pdt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(iceberg), POINTER, INTENT(IN) :: berg
    REAL(KIND = wp), INTENT(IN) :: pxi, pyj
    REAL(KIND = wp), INTENT(IN) :: puvel, pvvel
    REAL(KIND = wp), INTENT(IN) :: puvel0, pvvel0
    REAL(KIND = wp), INTENT(OUT) :: pe1, pe2
    REAL(KIND = wp), INTENT(INOUT) :: pax, pay
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), PARAMETER :: pp_alpha = 0._wp
    REAL(KIND = wp), PARAMETER :: pp_beta = 1._wp
    REAL(KIND = wp), PARAMETER :: pp_vel_lim = 15._wp
    REAL(KIND = wp), PARAMETER :: pp_accel_lim = 1.E-2_wp
    REAL(KIND = wp), PARAMETER :: pp_Cr0 = 0.06_wp
    INTEGER :: itloop
    REAL(KIND = wp) :: zuo, zui, zua, zuwave, zssh_x, zsst, zcn, zhi
    REAL(KIND = wp) :: zvo, zvi, zva, zvwave, zssh_y
    REAL(KIND = wp) :: zff, zT, zD, zW, zL, zM, zF
    REAL(KIND = wp) :: zdrag_ocn, zdrag_atm, zdrag_ice, zwave_rad
    REAL(KIND = wp) :: z_ocn, z_atm, z_ice
    REAL(KIND = wp) :: zampl, zwmod, zCr, zLwavelength, zLcutoff, zLtop
    REAL(KIND = wp) :: zlambda, zdetA, zA11, zA12, zaxe, zaye, zD_hi
    REAL(KIND = wp) :: zuveln, zvveln, zus, zvs, zspeed, zloc_dx, zspeed_new
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('icb_accel', 'r0', 0, 0)
    nknberg = berg % number(1)
    CALL icb_utl_interp(pxi, pe1, zuo, zui, zua, zssh_x, pyj, pe2, zvo, zvi, zva, zssh_y, zsst, zcn, zhi, zff)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zM = berg % current_point % mass
    zT = berg % current_point % thickness
    zD = (rn_rho_bergs / pp_rho_seawater) * zT
    zF = zT - zD
    zW = berg % current_point % width
    zL = berg % current_point % length
    zhi = MIN(zhi, zD)
    zD_hi = MAX(0._wp, zD - zhi)
    zuwave = zua - zuo
    zvwave = zva - zvo
    zwmod = zuwave * zuwave + zvwave * zvwave
    zampl = 0.5 * 0.02025 * zwmod
    zLwavelength = 0.32 * zwmod
    zLcutoff = 0.125 * zLwavelength
    zLtop = 0.25 * zLwavelength
    zCr = pp_Cr0 * MIN(MAX(0., (zL - zLcutoff) / ((zLtop - zLcutoff) + 1.E-30)), 1.)
    zwave_rad = 0.5 * pp_rho_seawater / zM * zCr * grav * zampl * MIN(zampl, zF) * (2. * zW * zL) / (zW + zL)
    zwmod = SQRT(zua * zua + zva * zva)
    IF (zwmod /= 0._wp) THEN
      zuwave = zua / zwmod
      zvwave = zva / zwmod
    ELSE
      zuwave = 0.
      zvwave = 0.
      zwave_rad = 0.
    END IF
    z_ocn = pp_rho_seawater / zM * (0.5 * pp_Cd_wv * zW * (zD_hi) + pp_Cd_wh * zW * zL)
    z_atm = pp_rho_air / zM * (0.5 * pp_Cd_av * zW * zF + pp_Cd_ah * zW * zL)
    z_ice = pp_rho_ice / zM * (0.5 * pp_Cd_iv * zW * zhi)
    IF (ABS(zui) + ABS(zvi) == 0._wp) z_ice = 0._wp
    zuveln = puvel
    zvveln = pvvel
    DO itloop = 1, 2
      zus = 0.5 * (zuveln + puvel)
      zvs = 0.5 * (zvveln + pvvel)
      zdrag_ocn = z_ocn * SQRT((zus - zuo) * (zus - zuo) + (zvs - zvo) * (zvs - zvo))
      zdrag_atm = z_atm * SQRT((zus - zua) * (zus - zua) + (zvs - zva) * (zvs - zva))
      zdrag_ice = z_ice * SQRT((zus - zui) * (zus - zui) + (zvs - zvi) * (zvs - zvi))
      zaxe = - grav * zssh_x + zwave_rad * zuwave
      zaye = - grav * zssh_y + zwave_rad * zvwave
      IF (pp_alpha > 0._wp) THEN
        zaxe = zaxe + zff * pvvel0
        zaye = zaye - zff * puvel0
      ELSE
        zaxe = zaxe + zff * pvvel
        zaye = zaye - zff * puvel
      END IF
      IF (pp_beta > 0._wp) THEN
        zaxe = zaxe - zdrag_ocn * (puvel0 - zuo) - zdrag_atm * (puvel0 - zua) - zdrag_ice * (puvel0 - zui)
        zaye = zaye - zdrag_ocn * (pvvel0 - zvo) - zdrag_atm * (pvvel0 - zva) - zdrag_ice * (pvvel0 - zvi)
      ELSE
        zaxe = zaxe - zdrag_ocn * (puvel - zuo) - zdrag_atm * (puvel - zua) - zdrag_ice * (puvel - zui)
        zaye = zaye - zdrag_ocn * (pvvel - zvo) - zdrag_atm * (pvvel - zva) - zdrag_ice * (pvvel - zvi)
      END IF
      IF (pp_alpha + pp_beta > 0._wp) THEN
        zlambda = zdrag_ocn + zdrag_atm + zdrag_ice
        zA11 = 1._wp + pp_beta * pdt * zlambda
        zA12 = pp_alpha * pdt * zff
        zdetA = 1._wp / (zA11 * zA11 + zA12 * zA12)
        pax = zdetA * (zA11 * zaxe + zA12 * zaye)
        pay = zdetA * (zA11 * zaye - zA12 * zaxe)
      ELSE
        pax = zaxe
        pay = zaye
      END IF
      zuveln = puvel0 + pdt * pax
      zvveln = pvvel0 + pdt * pay
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('icb_accel', 'r1', 0, 0)
    IF (rn_speed_limit > 0._wp) THEN
      zspeed = SQRT(zuveln * zuveln + zvveln * zvveln)
      IF (zspeed > 0._wp) THEN
        zloc_dx = MIN(pe1, pe2)
        zspeed_new = zloc_dx / pdt * rn_speed_limit
        IF (zspeed_new < zspeed) THEN
          zuveln = zuveln * (zspeed_new / zspeed)
          zvveln = zvveln * (zspeed_new / zspeed)
          CALL icb_dia_speed
        END IF
      END IF
    END IF
    IF (nn_verbose_level > 0) THEN
      IF (ABS(zuveln) > pp_vel_lim .OR. ABS(zvveln) > pp_vel_lim) WRITE(numicb, '("pe=",i3,x,a)') narea, 'Dump triggered by excessive velocity'
      IF (ABS(pax) > pp_accel_lim .OR. ABS(pay) > pp_accel_lim) WRITE(numicb, '("pe=",i3,x,a)') narea, 'Dump triggered by excessive acceleration'
    END IF
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE icb_accel
END MODULE icbdyn