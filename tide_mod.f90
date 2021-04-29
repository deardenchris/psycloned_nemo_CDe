MODULE tide_mod
  USE dom_oce
  USE phycst
  USE daymod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tide_harmo
  PUBLIC :: tide_init_Wave
  INTEGER, PUBLIC, PARAMETER :: jpmax_harmo = 19
  TYPE, PUBLIC :: tide
    CHARACTER(LEN = 4) :: cname_tide
    REAL(KIND = wp) :: equitide
    INTEGER :: nutide
    INTEGER :: nt, ns, nh, np, np1, shift
    INTEGER :: nksi, nnu0, nnu1, nnu2, R
    INTEGER :: nformula
  END TYPE tide
  TYPE(tide), PUBLIC, DIMENSION(jpmax_harmo) :: Wave
  REAL(KIND = wp) :: sh_T, sh_s, sh_h, sh_p, sh_p1
  REAL(KIND = wp) :: sh_xi, sh_nu, sh_nuprim, sh_nusec, sh_R
  REAL(KIND = wp) :: sh_I, sh_x1ra, sh_N
  CONTAINS
  SUBROUTINE tide_init_Wave
    Wave(1) = tide('M2', 0.242297, 2, 2, - 2, 2, 0, 0, 0, 2, - 2, 0, 0, 0, 78)
    Wave(2) = tide('N2', 0.046313, 2, 2, - 3, 2, 1, 0, 0, 2, - 2, 0, 0, 0, 78)
    Wave(3) = tide('2N2', 0.006184, 2, 2, - 4, 2, 2, 0, 0, 2, - 2, 0, 0, 0, 78)
    Wave(4) = tide('S2', 0.113572, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    Wave(5) = tide('K2', 0.030875, 2, 2, 0, 2, 0, 0, 0, 0, 0, 0, - 2, 0, 235)
    Wave(6) = tide('K1', 0.142408, 1, 1, 0, 1, 0, 0, - 90, 0, 0, - 1, 0, 0, 227)
    Wave(7) = tide('O1', 0.101266, 1, 1, - 2, 1, 0, 0, + 90, 2, - 1, 0, 0, 0, 75)
    Wave(8) = tide('Q1', 0.019387, 1, 1, - 3, 1, 1, 0, + 90, 2, - 1, 0, 0, 0, 75)
    Wave(9) = tide('P1', 0.047129, 1, 1, 0, - 1, 0, 0, + 90, 0, 0, 0, 0, 0, 0)
    Wave(10) = tide('M4', 0.000000, 4, 4, - 4, 4, 0, 0, 0, 4, - 4, 0, 0, 0, 1)
    Wave(11) = tide('Mf', 0.042017, 0, 0, 2, 0, 0, 0, 0, - 2, 0, 0, 0, 0, 74)
    Wave(12) = tide('Mm', 0.022191, 0, 0, 1, 0, - 1, 0, 0, 0, 0, 0, 0, 0, 73)
    Wave(13) = tide('Msqm', 0.000667, 0, 0, 4, - 2, 0, 0, 0, - 2, 0, 0, 0, 0, 74)
    Wave(14) = tide('Mtm', 0.008049, 0, 0, 3, 0, - 1, 0, 0, - 2, 0, 0, 0, 0, 74)
    Wave(15) = tide('S1', 0.000000, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    Wave(16) = tide('MU2', 0.005841, 2, 2, - 4, 4, 0, 0, 0, 2, - 2, 0, 0, 0, 78)
    Wave(17) = tide('NU2', 0.009094, 2, 2, - 3, 4, - 1, 0, 0, 2, - 2, 0, 0, 0, 78)
    Wave(18) = tide('L2', 0.006694, 2, 2, - 1, 2, - 1, 0, + 180, 2, - 2, 0, 0, 0, 215)
    Wave(19) = tide('T2', 0.006614, 2, 2, 0, - 1, 0, 1, 0, 0, 0, 0, 0, 0, 0)
  END SUBROUTINE tide_init_Wave
  SUBROUTINE tide_harmo(pomega, pvt, put, pcor, ktide, kc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(kc), INTENT(IN) :: ktide
    INTEGER, INTENT(IN) :: kc
    REAL(KIND = wp), DIMENSION(kc), INTENT(OUT) :: pomega
    REAL(KIND = wp), DIMENSION(kc), INTENT(OUT) :: pvt, put, pcor
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('tide_harmo', 'r0', 0, 0)
    CALL astronomic_angle
    CALL tide_pulse(pomega, ktide, kc)
    CALL tide_vuf(pvt, put, pcor, ktide, kc)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE tide_harmo
  SUBROUTINE astronomic_angle
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp) :: cosI, p, q, t2, t4, sin2I, s2, tgI2, P1, sh_tgn2, at1, at2
    REAL(KIND = wp) :: zqy, zsy, zday, zdj, zhfrac
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('astronomic_angle', 'r0', 0, 0)
    zqy = AINT((nyear - 1901.) / 4.)
    zsy = nyear - 1900.
    zdj = dayjul(nyear, nmonth, nday)
    zday = zdj + zqy - 1.
    zhfrac = nsec_day / 3600.
    sh_N = (259.1560564 - 19.328185764 * zsy - .0529539336 * zday - .0022064139 * zhfrac) * rad
    sh_T = (180. + zhfrac * (360. / 24.)) * rad
    sh_h = (280.1895014 - .238724988 * zsy + .9856473288 * zday + .0410686387 * zhfrac) * rad
    sh_s = (277.0256206 + 129.38482032 * zsy + 13.176396768 * zday + .549016532 * zhfrac) * rad
    sh_p1 = (281.2208569 + .01717836 * zsy + .000047064 * zday + .000001961 * zhfrac) * rad
    sh_p = (334.3837214 + 40.66246584 * zsy + .111404016 * zday + .004641834 * zhfrac) * rad
    sh_N = MOD(sh_N, 2 * rpi)
    sh_s = MOD(sh_s, 2 * rpi)
    sh_h = MOD(sh_h, 2 * rpi)
    sh_p = MOD(sh_p, 2 * rpi)
    sh_p1 = MOD(sh_p1, 2 * rpi)
    cosI = 0.913694997 - 0.035692561 * COS(sh_N)
    sh_I = ACOS(cosI)
    sin2I = SIN(sh_I)
    sh_tgn2 = TAN(sh_N / 2.0)
    at1 = ATAN(1.01883 * sh_tgn2)
    at2 = ATAN(0.64412 * sh_tgn2)
    sh_xi = - at1 - at2 + sh_N
    IF (sh_N > rpi) sh_xi = sh_xi - 2.0 * rpi
    sh_nu = at1 - at2
    tgI2 = TAN(sh_I / 2.0)
    P1 = sh_p - sh_xi
    t2 = tgI2 * tgI2
    t4 = t2 * t2
    sh_x1ra = SQRT(1.0 - 12.0 * t2 * COS(2.0 * P1) + 36.0 * t4)
    p = SIN(2.0 * P1)
    q = 1.0 / (6.0 * t2) - COS(2.0 * P1)
    sh_R = ATAN(p / q)
    p = SIN(2.0 * sh_I) * SIN(sh_nu)
    q = SIN(2.0 * sh_I) * COS(sh_nu) + 0.3347
    sh_nuprim = ATAN(p / q)
    s2 = SIN(sh_I) * SIN(sh_I)
    p = s2 * SIN(2.0 * sh_nu)
    q = s2 * COS(2.0 * sh_nu) + 0.0727
    sh_nusec = 0.5 * ATAN(p / q)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE astronomic_angle
  SUBROUTINE tide_pulse(pomega, ktide, kc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kc
    INTEGER, DIMENSION(kc), INTENT(IN) :: ktide
    REAL(KIND = wp), DIMENSION(kc), INTENT(OUT) :: pomega
    INTEGER :: jh
    REAL(KIND = wp) :: zscale
    REAL(KIND = wp) :: zomega_T = 13149000.0_wp
    REAL(KIND = wp) :: zomega_s = 481267.892_wp
    REAL(KIND = wp) :: zomega_h = 36000.76892_wp
    REAL(KIND = wp) :: zomega_p = 4069.0322056_wp
    REAL(KIND = wp) :: zomega_n = 1934.1423972_wp
    REAL(KIND = wp) :: zomega_p1 = 1.719175_wp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('tide_pulse', 'r0', 0, 0)
    zscale = rad / (36525._wp * 86400._wp)
    DO jh = 1, kc
      pomega(jh) = (zomega_T * Wave(ktide(jh)) % nT + zomega_s * Wave(ktide(jh)) % ns + zomega_h * Wave(ktide(jh)) % nh + zomega_p * Wave(ktide(jh)) % np + zomega_p1 * Wave(ktide(jh)) % np1) * zscale
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE tide_pulse
  SUBROUTINE tide_vuf(pvt, put, pcor, ktide, kc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kc
    INTEGER, DIMENSION(kc), INTENT(IN) :: ktide
    REAL(KIND = wp), DIMENSION(kc), INTENT(OUT) :: pvt, put, pcor
    INTEGER :: jh
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('tide_vuf', 'r0', 0, 0)
    DO jh = 1, kc
      pvt(jh) = sh_T * Wave(ktide(jh)) % nT + sh_s * Wave(ktide(jh)) % ns + sh_h * Wave(ktide(jh)) % nh + sh_p * Wave(ktide(jh)) % np + sh_p1 * Wave(ktide(jh)) % np1 + Wave(ktide(jh)) % shift * rad
      put(jh) = sh_xi * Wave(ktide(jh)) % nksi + sh_nu * Wave(ktide(jh)) % nnu0 + sh_nuprim * Wave(ktide(jh)) % nnu1 + sh_nusec * Wave(ktide(jh)) % nnu2 + sh_R * Wave(ktide(jh)) % R
      pcor(jh) = nodal_factort(Wave(ktide(jh)) % nformula)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE tide_vuf
  RECURSIVE FUNCTION nodal_factort(kformula) RESULT(zf)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kformula
    REAL(KIND = wp) :: zf
    REAL(KIND = wp) :: zs, zf1, zf2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('nodal_factort', 'r0', 0, 0)
    SELECT CASE (kformula)
    CASE (0)
      zf = 1.0
    CASE (1)
      zf = nodal_factort(78)
      zf = zf * zf
    CASE (2)
      zf1 = nodal_factort(78)
      zf = nodal_factort(0)
      zf = zf1 * zf
    CASE (4)
      zf1 = nodal_factort(78)
      zf = nodal_factort(235)
      zf = zf1 * zf
    CASE (5)
      zf1 = nodal_factort(78)
      zf = nodal_factort(235)
      zf = zf * zf1 * zf1
    CASE (6)
      zf1 = nodal_factort(78)
      zf = nodal_factort(0)
      zf = zf * zf1 * zf1
    CASE (7)
      zf = nodal_factort(75)
      zf = zf * zf
    CASE (8)
      zf = nodal_factort(78)
      zf1 = nodal_factort(0)
      zf2 = nodal_factort(235)
      zf = zf * zf1 * zf2
    CASE (9)
      zf = nodal_factort(78)
      zf1 = nodal_factort(0)
      zf2 = nodal_factort(227)
      zf = zf * zf1 * zf2
    CASE (10)
      zf = nodal_factort(78)
      zf1 = nodal_factort(227)
      zf = zf * zf1
    CASE (11)
      zf = nodal_factort(75)
      zf1 = nodal_factort(0)
      zf = zf * zf1
    CASE (12)
      zf1 = nodal_factort(78)
      zf = nodal_factort(0)
      zf = zf * zf1 * zf1 * zf1
    CASE (13)
      zf1 = nodal_factort(78)
      zf = nodal_factort(75)
      zf = zf * zf1
    CASE (14)
      zf = nodal_factort(235)
      zf1 = nodal_factort(0)
      zf = zf * zf1
    CASE (15)
      zf = nodal_factort(235)
      zf1 = nodal_factort(75)
      zf = zf * zf1
    CASE (16)
      zf = nodal_factort(78)
      zf1 = nodal_factort(0)
      zf = zf * zf1 * zf1
    CASE (17)
      zf1 = nodal_factort(227)
      zf = nodal_factort(0)
      zf = zf * zf1
    CASE (18)
      zf1 = nodal_factort(78)
      zf = zf1 * zf1 * zf1
    CASE (19)
      zf = nodal_factort(78)
      zf1 = nodal_factort(0)
      zf = zf * zf1 * zf1
    CASE (73)
      zs = SIN(sh_I)
      zf = (2. / 3. - zs * zs) / 0.5021
    CASE (74)
      zs = SIN(sh_I)
      zf = zs * zs / 0.1578
    CASE (75)
      zs = COS(sh_I / 2)
      zf = SIN(sh_I) * zs * zs / 0.3800
    CASE (76)
      zf = SIN(2 * sh_I) / 0.7214
    CASE (77)
      zs = SIN(sh_I / 2)
      zf = SIN(sh_I) * zs * zs / 0.0164
    CASE (78)
      zs = COS(sh_I / 2)
      zf = zs * zs * zs * zs / 0.9154
    CASE (79)
      zs = SIN(sh_I)
      zf = zs * zs / 0.1565
    CASE (144)
      zs = SIN(sh_I / 2)
      zf = (1 - 10 * zs * zs + 15 * zs * zs * zs * zs) * COS(sh_I / 2) / 0.5873
    CASE (149)
      zs = COS(sh_I / 2)
      zf = zs * zs * zs * zs * zs * zs / 0.8758
    CASE (215)
      zs = COS(sh_I / 2)
      zf = zs * zs * zs * zs / 0.9154 * sh_x1ra
    CASE (227)
      zs = SIN(2 * sh_I)
      zf = SQRT(0.8965 * zs * zs + 0.6001 * zs * COS(sh_nu) + 0.1006)
    CASE (235)
      zs = SIN(sh_I)
      zf = SQRT(19.0444 * zs * zs * zs * zs + 2.7702 * zs * zs * COS(2 * sh_nu) + .0981)
    END SELECT
    CALL profile_psy_data0 % PostEnd
  END FUNCTION nodal_factort
  FUNCTION dayjul(kyr, kmonth, kday)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kyr, kmonth, kday
    INTEGER, DIMENSION(12) :: idayt, idays
    INTEGER :: inc, ji
    REAL(KIND = wp) :: dayjul, zyq
    DATA idayt / 0., 31., 59., 90., 120., 151., 181., 212., 243., 273., 304., 334. /
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dayjul', 'r0', 0, 0)
    idays(1) = 0.
    idays(2) = 31.
    inc = 0.
    zyq = MOD(kyr - 1900., 4.)
    IF (zyq == 0.) inc = 1.
    DO ji = 3, 12
      idays(ji) = idayt(ji) + inc
    END DO
    dayjul = idays(kmonth) + kday
    CALL profile_psy_data0 % PostEnd
  END FUNCTION dayjul
END MODULE tide_mod