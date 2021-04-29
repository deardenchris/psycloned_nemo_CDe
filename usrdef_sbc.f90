MODULE usrdef_sbc
  USE oce
  USE dom_oce
  USE sbc_oce
  USE phycst
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: usrdef_sbc_oce
  PUBLIC :: usrdef_sbc_ice_tau
  PUBLIC :: usrdef_sbc_ice_flx
  CONTAINS
  SUBROUTINE usrdef_sbc_oce(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    INTEGER :: zyear0
    INTEGER :: zmonth0
    INTEGER :: zday0
    INTEGER :: zday_year0
    REAL(KIND = wp) :: ztau, ztau_sais
    REAL(KIND = wp) :: ztime
    REAL(KIND = wp) :: ztimemax, ztimemin
    REAL(KIND = wp) :: ztimemax1, ztimemin1
    REAL(KIND = wp) :: ztimemax2, ztimemin2
    REAL(KIND = wp) :: ztaun
    REAL(KIND = wp) :: zemp_s, zemp_n, zemp_sais, ztstar
    REAL(KIND = wp) :: zcos_sais1, zcos_sais2, ztrp, zconv, t_star
    REAL(KIND = wp) :: zsumemp, zsurf
    REAL(KIND = wp) :: zrhoa = 1.22
    REAL(KIND = wp) :: zcdrag = 1.5E-3
    REAL(KIND = wp) :: ztx, zty, zmod, zcoef
    REAL(KIND = wp) :: zyydd
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    !$ACC KERNELS
    zyydd = REAL(nyear_len(1), wp)
    zyear0 = ndate0 / 10000
    zmonth0 = (ndate0 - zyear0 * 10000) / 100
    zday0 = ndate0 - zyear0 * 10000 - zmonth0 * 100
    zday_year0 = (zmonth0 - 1) * 30. + zday0
    ztime = REAL(kt) * rdt / (rmmss * rhhmm) - (nyear - 1) * rjjhh * zyydd
    ztimemax1 = ((5. * 30.) + 21.) * 24.
    ztimemin1 = ztimemax1 + rjjhh * zyydd / 2
    ztimemax2 = ((6. * 30.) + 21.) * 24.
    ztimemin2 = ztimemax2 - rjjhh * zyydd / 2
    zemp_S = 0.7
    zemp_N = 0.8
    zemp_sais = 0.1
    zTstar = 28.3
    zcos_sais1 = COS((ztime - ztimemax1) / (ztimemin1 - ztimemax1) * rpi)
    zcos_sais2 = COS((ztime - ztimemax2) / (ztimemax2 - ztimemin2) * rpi)
    ztrp = - 40.E0
    zconv = 3.16E-5
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        t_star = zTstar * (1. + 1. / 50. * zcos_sais2) * COS(rpi * (gphit(ji, jj) - 5.) / (53.5 * (1 + 11 / 53.5 * zcos_sais2) * 2.))
        qsr(ji, jj) = 230 * COS(3.1415 * (gphit(ji, jj) - 23.5 * zcos_sais1) / (0.9 * 180))
        qns(ji, jj) = ztrp * (tsb(ji, jj, 1, jp_tem) - t_star) - qsr(ji, jj)
        IF (gphit(ji, jj) >= 14.845 .AND. 37.2 >= gphit(ji, jj)) THEN
          emp(ji, jj) = zemp_S * zconv * SIN(rpi / 2 * (gphit(ji, jj) - 37.2) / (24.6 - 37.2)) * (1 - zemp_sais / zemp_S * zcos_sais1)
        ELSE
          emp(ji, jj) = - zemp_N * zconv * SIN(rpi / 2 * (gphit(ji, jj) - 37.2) / (46.8 - 37.2)) * (1 - zemp_sais / zemp_N * zcos_sais1)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('usrdef_sbc_oce', 'r0', 0, 0)
    zsumemp = GLOB_SUM('usrdef_sbc', emp(:, :))
    zsurf = GLOB_SUM('usrdef_sbc', tmask(:, :, 1))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zsumemp = zsumemp / zsurf
    emp(:, :) = emp(:, :) - zsumemp * tmask(:, :, 1)
    sfx(:, :) = 0.0_wp
    qns(:, :) = qns(:, :) - emp(:, :) * sst_m(:, :) * rcp
    zyear0 = ndate0 / 10000
    zmonth0 = (ndate0 - zyear0 * 10000) / 100
    zday0 = ndate0 - zyear0 * 10000 - zmonth0 * 100
    zday_year0 = (zmonth0 - 1) * 30. + zday0
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('usrdef_sbc_oce', 'r1', 0, 0)
    ztime = FLOAT(kt) * rdt / (rmmss * rhhmm) - (nyear - 1) * rjjhh * zyydd
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    ztimemax = ((5. * 30.) + 21.) * 24.
    ztimemin = ztimemax + rjjhh * zyydd / 2
    ztau = 0.105 / SQRT(2.)
    ztau_sais = 0.015
    ztaun = ztau - ztau_sais * COS((ztime - ztimemax) / (ztimemin - ztimemax) * rpi)
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        utau(ji, jj) = - ztaun * SIN(rpi * (gphiu(ji, jj) - 15.) / (29. - 15.))
        vtau(ji, jj) = ztaun * SIN(rpi * (gphiv(ji, jj) - 15.) / (29. - 15.))
      END DO
    END DO
    zcoef = 1. / (zrhoa * zcdrag)
    !$ACC loop independent collapse(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ztx = utau(ji - 1, jj) + utau(ji, jj)
        zty = vtau(ji, jj - 1) + vtau(ji, jj)
        zmod = 0.5 * SQRT(ztx * ztx + zty * zty)
        taum(ji, jj) = zmod
        wndm(ji, jj) = SQRT(zmod * zcoef)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('usrdef_sbc_oce', 'r2', 0, 0)
    CALL lbc_lnk_multi('usrdef_sbc', taum(:, :), 'T', 1., wndm(:, :), 'T', 1.)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'usrdef_sbc_oce : analytical surface fluxes for GYRE configuration'
      WRITE(numout, *) '~~~~~~~~~~~ '
      WRITE(numout, *) '           nyear      = ', nyear
      WRITE(numout, *) '           nmonth     = ', nmonth
      WRITE(numout, *) '           nday       = ', nday
      WRITE(numout, *) '           nday_year  = ', nday_year
      WRITE(numout, *) '           ztime      = ', ztime
      WRITE(numout, *) '           ztimemax   = ', ztimemax
      WRITE(numout, *) '           ztimemin   = ', ztimemin
      WRITE(numout, *) '           ztimemax1  = ', ztimemax1
      WRITE(numout, *) '           ztimemin1  = ', ztimemin1
      WRITE(numout, *) '           ztimemax2  = ', ztimemax2
      WRITE(numout, *) '           ztimemin2  = ', ztimemin2
      WRITE(numout, *) '           zyear0     = ', zyear0
      WRITE(numout, *) '           zmonth0    = ', zmonth0
      WRITE(numout, *) '           zday0      = ', zday0
      WRITE(numout, *) '           zday_year0 = ', zday_year0
      WRITE(numout, *) '           zyydd      = ', zyydd
      WRITE(numout, *) '           zemp_S     = ', zemp_S
      WRITE(numout, *) '           zemp_N     = ', zemp_N
      WRITE(numout, *) '           zemp_sais  = ', zemp_sais
      WRITE(numout, *) '           zTstar     = ', zTstar
      WRITE(numout, *) '           zsumemp    = ', zsumemp
      WRITE(numout, *) '           zsurf      = ', zsurf
      WRITE(numout, *) '           ztrp       = ', ztrp
      WRITE(numout, *) '           zconv      = ', zconv
      WRITE(numout, *) '           ndastp     = ', ndastp
      WRITE(numout, *) '           adatrj     = ', adatrj
    END IF
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE usrdef_sbc_oce
  SUBROUTINE usrdef_sbc_ice_tau(kt)
    INTEGER, INTENT(IN) :: kt
  END SUBROUTINE usrdef_sbc_ice_tau
  SUBROUTINE usrdef_sbc_ice_flx(kt, phs, phi)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phs
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phi
  END SUBROUTINE usrdef_sbc_ice_flx
END MODULE usrdef_sbc