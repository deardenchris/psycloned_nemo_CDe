MODULE obs_conv
  USE par_kind, ONLY: wp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: potemp, fspott, atg, theta, depth, p_to_dep, dep_to_p
  CONTAINS
  REAL(KIND = wp) FUNCTION potemp(ps, pt, pp, ppr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: ps
    REAL(KIND = wp), INTENT(IN) :: pt
    REAL(KIND = wp), INTENT(IN) :: pp
    REAL(KIND = wp), INTENT(IN) :: ppr
    REAL(KIND = wp) :: zpol
    REAL(KIND = wp), PARAMETER :: a1 = 1.067610E-05
    REAL(KIND = wp), PARAMETER :: a2 = - 1.434297E-06
    REAL(KIND = wp), PARAMETER :: a3 = - 7.566349E-09
    REAL(KIND = wp), PARAMETER :: a4 = - 8.535585E-06
    REAL(KIND = wp), PARAMETER :: a5 = 3.074672E-08
    REAL(KIND = wp), PARAMETER :: a6 = 1.918639E-08
    REAL(KIND = wp), PARAMETER :: a7 = 1.788718E-10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    zpol = a1 + a2 * ps + a3 * (pp + ppr) + a4 * pt + a5 * ps * pt + a6 * pt * pt + a7 * pt * (pp + ppr)
    potemp = pt + (pp - ppr) * zpol
    CALL profile_psy_data0 % PostEnd
  END FUNCTION potemp
  REAL(KIND = wp) FUNCTION fspott(pft, pfs, pfp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp) :: pft
    REAL(KIND = wp) :: pfs
    REAL(KIND = wp) :: pfp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    fspott = pft - pfp * ((3.6504E-4 + pft * (8.3198E-5 + pft * (- 5.4065E-7 + pft * 4.0274E-9))) + (pfs - 35.0) * (1.7439E-5 - &
&pft * 2.9778E-7) + pfp * (8.9309E-7 + pft * (- 3.1628E-8 + pft * 2.1987E-10) - (pfs - 35.0) * 4.1057E-9 + pfp * (- 1.6056E-10 + &
&pft * 5.0484E-12)))
    CALL profile_psy_data0 % PostEnd
  END FUNCTION fspott
  REAL(KIND = wp) FUNCTION atg(p_s, p_t, p_p)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: p_s
    REAL(KIND = wp), INTENT(IN) :: p_t
    REAL(KIND = wp), INTENT(IN) :: p_p
    REAL(KIND = wp) :: z_ds
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    z_ds = p_s - 35.0
    atg = (((- 2.1687E-16 * p_t + 1.8676E-14) * p_t - 4.6206E-13) * p_p + ((2.7759E-12 * p_t - 1.1351E-10) * z_ds + ((- 5.4481E-14 &
&* p_t + 8.733E-12) * p_t - 6.7795E-10) * p_t + 1.8741E-8)) * p_p + (- 4.2393E-8 * p_t + 1.8932E-6) * z_ds + ((6.6228E-10 * p_t - &
&6.836E-8) * p_t + 8.5258E-6) * p_t + 3.5803E-5
    CALL profile_psy_data0 % PostEnd
  END FUNCTION atg
  REAL(KIND = wp) FUNCTION theta(p_s, p_t0, p_p0, p_pr)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: p_s
    REAL(KIND = wp), INTENT(IN) :: p_t0
    REAL(KIND = wp), INTENT(IN) :: p_p0
    REAL(KIND = wp), INTENT(IN) :: p_pr
    REAL(KIND = wp) :: z_p
    REAL(KIND = wp) :: z_t
    REAL(KIND = wp) :: z_h
    REAL(KIND = wp) :: z_xk
    REAL(KIND = wp) :: z_q
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    z_p = p_p0
    z_t = p_t0
    z_h = p_pr - z_p
    z_xk = z_h * atg(p_s, z_t, z_p)
    Z_t = z_t + 0.5 * z_xk
    z_q = z_xk
    z_p = z_p + 0.5 * z_h
    z_xk = z_h * atg(p_s, z_t, z_p)
    z_t = z_t + 0.29289322 * (z_xk - z_q)
    z_q = 0.58578644 * z_xk + 0.121320344 * z_q
    z_xk = z_h * atg(p_s, z_t, z_p)
    z_t = z_t + 1.707106781 * (z_xk - z_q)
    z_q = 3.414213562 * z_xk - 4.121320244 * z_q
    z_p = z_p + 0.5 * z_h
    z_xk = z_h * atg(p_s, z_t, z_p)
    theta = z_t + (z_xk - 2.0 * z_q) / 6.0
    CALL profile_psy_data0 % PostEnd
  END FUNCTION theta
  REAL(KIND = wp) FUNCTION depth(p_p, p_lat)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: p_p
    REAL(KIND = wp), INTENT(IN) :: p_lat
    REAL(KIND = wp) :: z_x
    REAL(KIND = wp) :: z_gr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    z_x = SIN(p_lat / 57.29578)
    z_x = z_x * z_x
    z_gr = 9.780318 * (1.0 + (5.2788E-3 + 2.36E-5 * z_x) * z_x) + 1.092E-6 * p_p
    depth = (((- 1.82E-15 * p_p + 2.279E-10) * p_p - 2.2512E-5) * p_p + 9.72659) * p_p
    depth = depth / z_gr
    CALL profile_psy_data0 % PostEnd
  END FUNCTION depth
  REAL(KIND = wp) FUNCTION p_to_dep(p_p, p_lat)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: p_p
    REAL(KIND = wp), INTENT(IN) :: p_lat
    REAL(KIND = wp) :: z_x
    REAL(KIND = wp) :: z_c1
    REAL(KIND = wp) :: z_c2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    z_x = SIN(p_lat / 57.29578)
    z_x = z_x * z_x
    z_c1 = (5.92 + 5.25 * z_x) * 1E-3
    z_c2 = 2.21E-6
    p_to_dep = (1 - z_c1) * p_p - z_c2 * p_p * p_p
    CALL profile_psy_data0 % PostEnd
  END FUNCTION p_to_dep
  REAL(KIND = wp) FUNCTION dep_to_p(p_dep, p_lat)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: p_dep
    REAL(KIND = wp), INTENT(IN) :: p_lat
    REAL(KIND = wp) :: z_x
    REAL(KIND = wp) :: z_c1
    REAL(KIND = wp) :: z_c2
    REAL(KIND = wp) :: z_d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('wp', 'r0', 0, 0)
    z_x = SIN(p_lat / 57.29578)
    z_x = z_x * z_x
    z_c1 = (5.92 + 5.25 * z_x) * 1E-3
    z_c2 = 2.21E-6
    z_d = (z_c1 - 1) * (z_c1 - 1) - 4 * z_c2 * p_dep
    dep_to_p = ((1 - z_c1) - SQRT(z_d)) / (2 * z_c2)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION dep_to_p
END MODULE obs_conv