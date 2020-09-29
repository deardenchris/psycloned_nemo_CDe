MODULE icevar
  USE dom_oce
  USE phycst
  USE sbc_oce, ONLY: sss_m, ln_ice_embd, nn_fsbc
  USE ice
  USE ice1D
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_var_agg
  PUBLIC :: ice_var_glo2eqv
  PUBLIC :: ice_var_eqv2glo
  PUBLIC :: ice_var_salprof
  PUBLIC :: ice_var_salprof1d
  PUBLIC :: ice_var_zapsmall
  PUBLIC :: ice_var_zapneg
  PUBLIC :: ice_var_roundoff
  PUBLIC :: ice_var_itd
  PUBLIC :: ice_var_itd2
  PUBLIC :: ice_var_bv
  PUBLIC :: ice_var_enthalpy
  PUBLIC :: ice_var_sshdyn
  CONTAINS
  SUBROUTINE ice_var_agg(kn)
    INTEGER, INTENT(IN) :: kn
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z1_at_i, z1_vt_i, z1_vt_s
    !$ACC KERNELS
    vt_i(:, :) = SUM(v_i(:, :, :), dim = 3)
    vt_s(:, :) = SUM(v_s(:, :, :), dim = 3)
    at_i(:, :) = SUM(a_i(:, :, :), dim = 3)
    et_s(:, :) = SUM(SUM(e_s(:, :, :, :), dim = 4), dim = 3)
    et_i(:, :) = SUM(SUM(e_i(:, :, :, :), dim = 4), dim = 3)
    at_ip(:, :) = SUM(a_ip(:, :, :), dim = 3)
    vt_ip(:, :) = SUM(v_ip(:, :, :), dim = 3)
    ato_i(:, :) = 1._wp - at_i(:, :)
    !$ACC END KERNELS
    IF (kn > 1) THEN
      ALLOCATE(z1_at_i(jpi, jpj), z1_vt_i(jpi, jpj), z1_vt_s(jpi, jpj))
      !$ACC KERNELS
      WHERE (at_i(:, :) > epsi20)
        z1_at_i(:, :) = 1._wp / at_i(:, :)
      ELSEWHERE
        z1_at_i(:, :) = 0._wp
      END WHERE
      WHERE (vt_i(:, :) > epsi20)
        z1_vt_i(:, :) = 1._wp / vt_i(:, :)
      ELSEWHERE
        z1_vt_i(:, :) = 0._wp
      END WHERE
      WHERE (vt_s(:, :) > epsi20)
        z1_vt_s(:, :) = 1._wp / vt_s(:, :)
      ELSEWHERE
        z1_vt_s(:, :) = 0._wp
      END WHERE
      hm_i(:, :) = vt_i(:, :) * z1_at_i(:, :)
      hm_s(:, :) = vt_s(:, :) * z1_at_i(:, :)
      tm_su(:, :) = SUM(t_su(:, :, :) * a_i(:, :, :), dim = 3) * z1_at_i(:, :)
      tm_si(:, :) = SUM(t_si(:, :, :) * a_i(:, :, :), dim = 3) * z1_at_i(:, :)
      om_i(:, :) = SUM(oa_i(:, :, :), dim = 3) * z1_at_i(:, :)
      sm_i(:, :) = SUM(sv_i(:, :, :), dim = 3) * z1_vt_i(:, :)
      tm_i(:, :) = 0._wp
      tm_s(:, :) = 0._wp
      !$ACC END KERNELS
      DO jl = 1, jpl
        !$ACC KERNELS
        DO jk = 1, nlay_i
          tm_i(:, :) = tm_i(:, :) + r1_nlay_i * t_i(:, :, jk, jl) * v_i(:, :, jl) * z1_vt_i(:, :)
        END DO
        DO jk = 1, nlay_s
          tm_s(:, :) = tm_s(:, :) + r1_nlay_s * t_s(:, :, jk, jl) * v_s(:, :, jl) * z1_vt_s(:, :)
        END DO
        !$ACC END KERNELS
      END DO
      !$ACC KERNELS
      WHERE (at_i(:, :) <= epsi20)
        tm_su(:, :) = rt0
        tm_si(:, :) = rt0
        tm_i(:, :) = rt0
        tm_s(:, :) = rt0
      END WHERE
      !$ACC END KERNELS
      DEALLOCATE(z1_at_i, z1_vt_i, z1_vt_s)
    END IF
  END SUBROUTINE ice_var_agg
  SUBROUTINE ice_var_glo2eqv
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: ze_i
    REAL(KIND = wp) :: ze_s, ztmelts, zbbb, zccc
    REAL(KIND = wp) :: zhmax, z1_zhmax
    REAL(KIND = wp) :: zlay_i, zlay_s
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z1_a_i, z1_v_i
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    WHERE (a_i(:, :, :) > epsi20)
      z1_a_i(:, :, :) = 1._wp / a_i(:, :, :)
    ELSEWHERE
      z1_a_i(:, :, :) = 0._wp
    END WHERE
    WHERE (v_i(:, :, :) > epsi20)
      z1_v_i(:, :, :) = 1._wp / v_i(:, :, :)
    ELSEWHERE
      z1_v_i(:, :, :) = 0._wp
    END WHERE
    h_i(:, :, :) = v_i(:, :, :) * z1_a_i(:, :, :)
    zhmax = hi_max(jpl)
    z1_zhmax = 1._wp / hi_max(jpl)
    WHERE (h_i(:, :, jpl) > zhmax)
      h_i(:, :, jpl) = zhmax
      a_i(:, :, jpl) = v_i(:, :, jpl) * z1_zhmax
      z1_a_i(:, :, jpl) = zhmax * z1_v_i(:, :, jpl)
    END WHERE
    h_s(:, :, :) = v_s(:, :, :) * z1_a_i(:, :, :)
    o_i(:, :, :) = oa_i(:, :, :) * z1_a_i(:, :, :)
    a_ip_frac(:, :, :) = a_ip(:, :, :) * z1_a_i(:, :, :)
    WHERE (a_ip_frac(:, :, :) > epsi20)
      h_ip(:, :, :) = v_ip(:, :, :) * z1_a_i(:, :, :) / a_ip_frac(:, :, :)
    ELSEWHERE
      h_ip(:, :, :) = 0._wp
    END WHERE
    IF (nn_icesal == 2) THEN
      WHERE (v_i(:, :, :) > epsi20)
        s_i(:, :, :) = MAX(rn_simin, MIN(rn_simax, sv_i(:, :, :) * z1_v_i(:, :, :)))
      ELSEWHERE
        s_i(:, :, :) = rn_simin
      END WHERE
    END IF
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('ice_var_glo2eqv', 'r0', 0, 0)
    CALL ice_var_salprof
    zlay_i = REAL(nlay_i, wp)
    CALL profile_psy_data0 % PostEnd
    DO jl = 1, jpl
      !$ACC KERNELS
      DO jk = 1, nlay_i
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (v_i(ji, jj, jl) > epsi20) THEN
              ze_i = e_i(ji, jj, jk, jl) * z1_v_i(ji, jj, jl) * zlay_i
              ztmelts = - sz_i(ji, jj, jk, jl) * rTmlt
              zbbb = (rcp - rcpi) * ztmelts + ze_i * r1_rhoi - rLfus
              zccc = SQRT(MAX(zbbb * zbbb - 4._wp * rcpi * rLfus * ztmelts, 0._wp))
              t_i(ji, jj, jk, jl) = MAX(- 100._wp, MIN(- (zbbb + zccc) * 0.5_wp * r1_rcpi, ztmelts)) + rt0
            ELSE
              t_i(ji, jj, jk, jl) = rt0
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    CALL profile_psy_data1 % PreStart('ice_var_glo2eqv', 'r1', 0, 0)
    zlay_s = REAL(nlay_s, wp)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    DO jk = 1, nlay_s
      WHERE (v_s(:, :, :) > epsi20)
        t_s(:, :, jk, :) = rt0 + MAX(- 100._wp, MIN(r1_rcpi * (- r1_rhos * (e_s(:, :, jk, :) / v_s(:, :, :) * zlay_s) + rLfus), &
&0._wp))
      ELSEWHERE
        t_s(:, :, jk, :) = rt0
      END WHERE
    END DO
    vt_i(:, :) = SUM(v_i, dim = 3)
    vt_s(:, :) = SUM(v_s, dim = 3)
    at_i(:, :) = SUM(a_i, dim = 3)
    !$ACC END KERNELS
  END SUBROUTINE ice_var_glo2eqv
  SUBROUTINE ice_var_eqv2glo
    !$ACC KERNELS
    v_i(:, :, :) = h_i(:, :, :) * a_i(:, :, :)
    v_s(:, :, :) = h_s(:, :, :) * a_i(:, :, :)
    sv_i(:, :, :) = s_i(:, :, :) * v_i(:, :, :)
    v_ip(:, :, :) = h_ip(:, :, :) * a_ip(:, :, :)
    !$ACC END KERNELS
  END SUBROUTINE ice_var_eqv2glo
  SUBROUTINE ice_var_salprof
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: zsal, z1_dS
    REAL(KIND = wp) :: zargtemp, zs0, zs
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z_slope_s, zalpha
    REAL(KIND = wp), PARAMETER :: zsi0 = 3.5_wp
    REAL(KIND = wp), PARAMETER :: zsi1 = 4.5_wp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    SELECT CASE (nn_icesal)
    CASE (1)
      !$ACC KERNELS
      sz_i(:, :, :, :) = rn_icesal
      s_i(:, :, :) = rn_icesal
      !$ACC END KERNELS
    CASE (2)
      ALLOCATE(z_slope_s(jpi, jpj, jpl), zalpha(jpi, jpj, jpl))
      DO jl = 1, jpl
        !$ACC KERNELS
        DO jk = 1, nlay_i
          sz_i(:, :, jk, jl) = s_i(:, :, jl)
        END DO
        !$ACC END KERNELS
      END DO
      !$ACC KERNELS
      WHERE (h_i(:, :, :) > epsi20)
        z_slope_s(:, :, :) = 2._wp * s_i(:, :, :) / h_i(:, :, :)
      ELSEWHERE
        z_slope_s(:, :, :) = 0._wp
      END WHERE
      z1_dS = 1._wp / (zsi1 - zsi0)
      DO jl = 1, jpl
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zalpha(ji, jj, jl) = MAX(0._wp, MIN((zsi1 - s_i(ji, jj, jl)) * z1_dS, 1._wp))
            IF (2._wp * s_i(ji, jj, jl) >= sss_m(ji, jj)) zalpha(ji, jj, jl) = 0._wp
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data0 % PreStart('ice_var_salprof', 'r0', 0, 0)
      DO jl = 1, jpl
        DO jk = 1, nlay_i
          DO jj = 1, jpj
            DO ji = 1, jpi
              zs0 = z_slope_s(ji, jj, jl) * (REAL(jk, wp) - 0.5_wp) * h_i(ji, jj, jl) * r1_nlay_i
              zs = zalpha(ji, jj, jl) * zs0 + (1._wp - zalpha(ji, jj, jl)) * s_i(ji, jj, jl)
              sz_i(ji, jj, jk, jl) = MIN(rn_simax, MAX(zs, rn_simin))
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(z_slope_s, zalpha)
      CALL profile_psy_data0 % PostEnd
    CASE (3)
      !$ACC KERNELS
      s_i(:, :, :) = 2.30_wp
      !$ACC END KERNELS
      DO jl = 1, jpl
        DO jk = 1, nlay_i
          CALL profile_psy_data1 % PreStart('ice_var_salprof', 'r1', 0, 0)
          zargtemp = (REAL(jk, wp) - 0.5_wp) * r1_nlay_i
          CALL profile_psy_data1 % PostEnd
          !$ACC KERNELS
          sz_i(:, :, jk, jl) = 1.6_wp * (1._wp - COS(rpi * zargtemp ** (0.407_wp / (0.573_wp + zargtemp))))
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
  END SUBROUTINE ice_var_salprof
  SUBROUTINE ice_var_salprof1d
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jk
    REAL(KIND = wp) :: zargtemp, zsal, z1_dS
    REAL(KIND = wp) :: zs, zs0
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: z_slope_s, zalpha
    REAL(KIND = wp), PARAMETER :: zsi0 = 3.5_wp
    REAL(KIND = wp), PARAMETER :: zsi1 = 4.5_wp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    SELECT CASE (nn_icesal)
    CASE (1)
      !$ACC KERNELS
      sz_i_1d(1 : npti, :) = rn_icesal
      !$ACC END KERNELS
    CASE (2)
      CALL profile_psy_data0 % PreStart('ice_var_salprof1d', 'r0', 0, 0)
      ALLOCATE(z_slope_s(jpij), zalpha(jpij))
      WHERE (h_i_1d(1 : npti) > epsi20)
        z_slope_s(1 : npti) = 2._wp * s_i_1d(1 : npti) / h_i_1d(1 : npti)
      ELSEWHERE
        z_slope_s(1 : npti) = 0._wp
      END WHERE
      z1_dS = 1._wp / (zsi1 - zsi0)
      DO ji = 1, npti
        zalpha(ji) = MAX(0._wp, MIN((zsi1 - s_i_1d(ji)) * z1_dS, 1._wp))
        IF (2._wp * s_i_1d(ji) >= sss_1d(ji)) zalpha(ji) = 0._wp
      END DO
      DO jk = 1, nlay_i
        DO ji = 1, npti
          zs0 = z_slope_s(ji) * (REAL(jk, wp) - 0.5_wp) * h_i_1d(ji) * r1_nlay_i
          zs = zalpha(ji) * zs0 + (1._wp - zalpha(ji)) * s_i_1d(ji)
          sz_i_1d(ji, jk) = MIN(rn_simax, MAX(zs, rn_simin))
        END DO
      END DO
      DEALLOCATE(z_slope_s, zalpha)
      CALL profile_psy_data0 % PostEnd
    CASE (3)
      !$ACC KERNELS
      s_i_1d(1 : npti) = 2.30_wp
      !$ACC END KERNELS
      DO jk = 1, nlay_i
        CALL profile_psy_data1 % PreStart('ice_var_salprof1d', 'r1', 0, 0)
        zargtemp = (REAL(jk, wp) - 0.5_wp) * r1_nlay_i
        CALL profile_psy_data1 % PostEnd
        !$ACC KERNELS
        zsal = 1.6_wp * (1._wp - COS(rpi * zargtemp ** (0.407_wp / (0.573_wp + zargtemp))))
        DO ji = 1, npti
          sz_i_1d(ji, jk) = zsal
        END DO
        !$ACC END KERNELS
      END DO
    END SELECT
  END SUBROUTINE ice_var_salprof1d
  SUBROUTINE ice_var_zapsmall
    INTEGER :: ji, jj, jl, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zswitch
    !$ACC KERNELS
    DO jl = 1, jpl
      WHERE (a_i(:, :, jl) > epsi10)
        h_i(:, :, jl) = v_i(:, :, jl) / a_i(:, :, jl)
      ELSEWHERE
        h_i(:, :, jl) = 0._wp
      END WHERE
      WHERE (a_i(:, :, jl) < epsi10 .OR. v_i(:, :, jl) < epsi10 .OR. h_i(:, :, jl) < epsi10)
        zswitch(:, :) = 0._wp
      ELSEWHERE
        zswitch(:, :) = 1._wp
      END WHERE
      DO jk = 1, nlay_i
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            hfx_res(ji, jj) = hfx_res(ji, jj) - (1._wp - zswitch(ji, jj)) * e_i(ji, jj, jk, jl) * r1_rdtice
            e_i(ji, jj, jk, jl) = e_i(ji, jj, jk, jl) * zswitch(ji, jj)
            t_i(ji, jj, jk, jl) = t_i(ji, jj, jk, jl) * zswitch(ji, jj) + rt0 * (1._wp - zswitch(ji, jj))
          END DO
        END DO
      END DO
      DO jk = 1, nlay_s
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            hfx_res(ji, jj) = hfx_res(ji, jj) - (1._wp - zswitch(ji, jj)) * e_s(ji, jj, jk, jl) * r1_rdtice
            e_s(ji, jj, jk, jl) = e_s(ji, jj, jk, jl) * zswitch(ji, jj)
            t_s(ji, jj, jk, jl) = t_s(ji, jj, jk, jl) * zswitch(ji, jj) + rt0 * (1._wp - zswitch(ji, jj))
          END DO
        END DO
      END DO
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          sfx_res(ji, jj) = sfx_res(ji, jj) + (1._wp - zswitch(ji, jj)) * sv_i(ji, jj, jl) * rhoi * r1_rdtice
          wfx_res(ji, jj) = wfx_res(ji, jj) + (1._wp - zswitch(ji, jj)) * v_i(ji, jj, jl) * rhoi * r1_rdtice
          wfx_res(ji, jj) = wfx_res(ji, jj) + (1._wp - zswitch(ji, jj)) * v_s(ji, jj, jl) * rhos * r1_rdtice
          a_i(ji, jj, jl) = a_i(ji, jj, jl) * zswitch(ji, jj)
          v_i(ji, jj, jl) = v_i(ji, jj, jl) * zswitch(ji, jj)
          v_s(ji, jj, jl) = v_s(ji, jj, jl) * zswitch(ji, jj)
          t_su(ji, jj, jl) = t_su(ji, jj, jl) * zswitch(ji, jj) + t_bo(ji, jj) * (1._wp - zswitch(ji, jj))
          oa_i(ji, jj, jl) = oa_i(ji, jj, jl) * zswitch(ji, jj)
          sv_i(ji, jj, jl) = sv_i(ji, jj, jl) * zswitch(ji, jj)
          h_i(ji, jj, jl) = h_i(ji, jj, jl) * zswitch(ji, jj)
          h_s(ji, jj, jl) = h_s(ji, jj, jl) * zswitch(ji, jj)
          a_ip(ji, jj, jl) = a_ip(ji, jj, jl) * zswitch(ji, jj)
          v_ip(ji, jj, jl) = v_ip(ji, jj, jl) * zswitch(ji, jj)
        END DO
      END DO
    END DO
    at_i(:, :) = SUM(a_i(:, :, :), dim = 3)
    vt_i(:, :) = SUM(v_i(:, :, :), dim = 3)
    WHERE (at_i(:, :) == 0._wp) ato_i(:, :) = 1._wp
    !$ACC END KERNELS
  END SUBROUTINE ice_var_zapsmall
  SUBROUTINE ice_var_zapneg(pdt, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i)
    REAL(KIND = wp), INTENT(IN) :: pdt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pato_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_s
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: psv_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: poa_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pa_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pa_ip
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_ip
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_s
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_i
    INTEGER :: ji, jj, jl, jk
    REAL(KIND = wp) :: z1_dt
    !$ACC KERNELS
    z1_dt = 1._wp / pdt
    DO jl = 1, jpl
      WHERE (pv_i(:, :, :) <= 0._wp) pa_i(:, :, :) = 0._wp
      DO jk = 1, nlay_i
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (pe_i(ji, jj, jk, jl) < 0._wp .OR. pa_i(ji, jj, jl) <= 0._wp) THEN
              hfx_res(ji, jj) = hfx_res(ji, jj) - pe_i(ji, jj, jk, jl) * z1_dt
              pe_i(ji, jj, jk, jl) = 0._wp
            END IF
          END DO
        END DO
      END DO
      DO jk = 1, nlay_s
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (pe_s(ji, jj, jk, jl) < 0._wp .OR. pa_i(ji, jj, jl) <= 0._wp) THEN
              hfx_res(ji, jj) = hfx_res(ji, jj) - pe_s(ji, jj, jk, jl) * z1_dt
              pe_s(ji, jj, jk, jl) = 0._wp
            END IF
          END DO
        END DO
      END DO
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (pv_i(ji, jj, jl) < 0._wp .OR. pa_i(ji, jj, jl) <= 0._wp) THEN
            wfx_res(ji, jj) = wfx_res(ji, jj) + pv_i(ji, jj, jl) * rhoi * z1_dt
            pv_i(ji, jj, jl) = 0._wp
          END IF
          IF (pv_s(ji, jj, jl) < 0._wp .OR. pa_i(ji, jj, jl) <= 0._wp) THEN
            wfx_res(ji, jj) = wfx_res(ji, jj) + pv_s(ji, jj, jl) * rhos * z1_dt
            pv_s(ji, jj, jl) = 0._wp
          END IF
          IF (psv_i(ji, jj, jl) < 0._wp .OR. pa_i(ji, jj, jl) <= 0._wp) THEN
            sfx_res(ji, jj) = sfx_res(ji, jj) + psv_i(ji, jj, jl) * rhoi * z1_dt
            psv_i(ji, jj, jl) = 0._wp
          END IF
        END DO
      END DO
    END DO
    WHERE (pato_i(:, :) < 0._wp) pato_i(:, :) = 0._wp
    WHERE (poa_i(:, :, :) < 0._wp) poa_i(:, :, :) = 0._wp
    WHERE (pa_i(:, :, :) < 0._wp) pa_i(:, :, :) = 0._wp
    WHERE (pa_ip(:, :, :) < 0._wp) pa_ip(:, :, :) = 0._wp
    WHERE (pv_ip(:, :, :) < 0._wp) pv_ip(:, :, :) = 0._wp
    !$ACC END KERNELS
  END SUBROUTINE ice_var_zapneg
  SUBROUTINE ice_var_roundoff(pa_i, pv_i, pv_s, psv_i, poa_i, pa_ip, pv_ip, pe_s, pe_i)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pa_i
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pv_i
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pv_s
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: psv_i
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: poa_i
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pa_ip
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pv_ip
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pe_s
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pe_i
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !CALL profile_psy_data0 % PreStart('ice_var_roundoff', 'r0', 0, 0)
    !$ACC KERNELS ! CDe added
    WHERE (pa_i(1 : npti, :) < 0._wp .AND. pa_i(1 : npti, :) > - epsi10) pa_i(1 : npti, :) = 0._wp
    WHERE (pv_i(1 : npti, :) < 0._wp .AND. pv_i(1 : npti, :) > - epsi10) pv_i(1 : npti, :) = 0._wp
    WHERE (pv_s(1 : npti, :) < 0._wp .AND. pv_s(1 : npti, :) > - epsi10) pv_s(1 : npti, :) = 0._wp
    WHERE (psv_i(1 : npti, :) < 0._wp .AND. psv_i(1 : npti, :) > - epsi10) psv_i(1 : npti, :) = 0._wp
    WHERE (poa_i(1 : npti, :) < 0._wp .AND. poa_i(1 : npti, :) > - epsi10) poa_i(1 : npti, :) = 0._wp
    WHERE (pe_i(1 : npti, :, :) < 0._wp .AND. pe_i(1 : npti, :, :) > - epsi06) pe_i(1 : npti, :, :) = 0._wp
    WHERE (pe_s(1 : npti, :, :) < 0._wp .AND. pe_s(1 : npti, :, :) > - epsi06) pe_s(1 : npti, :, :) = 0._wp
    IF (ln_pnd_H12) THEN
      WHERE (pa_ip(1 : npti, :) < 0._wp .AND. pa_ip(1 : npti, :) > - epsi10) pa_ip(1 : npti, :) = 0._wp
      WHERE (pv_ip(1 : npti, :) < 0._wp .AND. pv_ip(1 : npti, :) > - epsi10) pv_ip(1 : npti, :) = 0._wp
    END IF
    !$ACC END KERNELS
    !CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_var_roundoff
  SUBROUTINE ice_var_itd(zhti, zhts, zati, zh_i, zh_s, za_i)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jk, jl
    INTEGER :: idim, i_fill, jl0
    REAL(KIND = wp) :: zarg, zV, zconv, zdh, zdv
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: zhti, zhts, zati
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: zh_i, zh_s, za_i
    INTEGER, DIMENSION(4) :: itest
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    idim = SIZE(zhti, 1)
    zh_i(1 : idim, 1 : jpl) = 0._wp
    zh_s(1 : idim, 1 : jpl) = 0._wp
    za_i(1 : idim, 1 : jpl) = 0._wp
    !$ACC END KERNELS
    DO ji = 1, idim
      IF (zhti(ji) > 0._wp) THEN
        CALL profile_psy_data0 % PreStart('ice_var_itd', 'r0', 0, 0)
        jl0 = jpl
        DO jl = 1, jpl
          IF ((zhti(ji) >= hi_max(jl - 1)) .AND. (zhti(ji) < hi_max(jl))) THEN
            jl0 = jl
            CYCLE
          END IF
        END DO
        CALL profile_psy_data0 % PostEnd
        !$ACC KERNELS
        itest(:) = 0
        i_fill = jpl + 1
        !$ACC END KERNELS
        DO WHILE ((SUM(itest(:)) /= 4) .AND. (i_fill >= 2))
          i_fill = i_fill - 1
          zh_i(ji, 1 : jpl) = 0._wp
          za_i(ji, 1 : jpl) = 0._wp
          itest(:) = 0
          IF (i_fill == 1) THEN
            zh_i(ji, 1) = zhti(ji)
            za_i(ji, 1) = zati(ji)
          ELSE
            DO jl = 1, i_fill - 1
              zh_i(ji, jl) = hi_mean(jl)
            END DO
            za_i(ji, jl0) = zati(ji) / SQRT(REAL(jpl))
            DO jl = 1, i_fill - 1
              IF (jl /= jl0) THEN
                zarg = (zh_i(ji, jl) - zhti(ji)) / (zhti(ji) * 0.5_wp)
                za_i(ji, jl) = za_i(ji, jl0) * EXP(- zarg ** 2)
              END IF
            END DO
            za_i(ji, i_fill) = zati(ji) - SUM(za_i(ji, 1 : i_fill - 1))
            zV = SUM(za_i(ji, 1 : i_fill - 1) * zh_i(ji, 1 : i_fill - 1))
            zh_i(ji, i_fill) = (zhti(ji) * zati(ji) - zV) / MAX(za_i(ji, i_fill), epsi10)
            IF (jl0 /= jpl) THEN
              DO jl = jpl, jl0 + 1, - 1
                IF (za_i(ji, jl) > za_i(ji, jl - 1)) THEN
                  zdv = zh_i(ji, jl) * za_i(ji, jl)
                  zh_i(ji, jl) = 0._wp
                  za_i(ji, jl) = 0._wp
                  za_i(ji, 1 : jl - 1) = za_i(ji, 1 : jl - 1) + zdv / MAX(REAL(jl - 1) * zhti(ji), epsi10)
                END IF
              END DO
            END IF
          END IF
          zconv = ABS(zati(ji) - SUM(za_i(ji, 1 : jpl)))
          IF (zconv < epsi06) itest(1) = 1
          zconv = ABS(zhti(ji) * zati(ji) - SUM(za_i(ji, 1 : jpl) * zh_i(ji, 1 : jpl)))
          IF (zconv < epsi06) itest(2) = 1
          IF (zh_i(ji, i_fill) >= hi_max(i_fill - 1)) itest(3) = 1
          itest(4) = 1
          DO jl = 1, i_fill
            IF (za_i(ji, jl) < 0._wp) itest(4) = 0
          END DO
        END DO
      END IF
    END DO
    !$ACC KERNELS
    DO jl = 1, jpl
      DO ji = 1, idim
        IF (za_i(ji, jl) > 0._wp) THEN
          zh_s(ji, jl) = zh_i(ji, jl) * (zhts(ji) / zhti(ji))
          zdh = MAX(0._wp, (rhos * zh_s(ji, jl) + (rhoi - rau0) * zh_i(ji, jl)) * r1_rau0)
          zh_i(ji, jl) = MIN(hi_max(jl), zh_i(ji, jl) + zdh)
          zh_s(ji, jl) = MAX(0._wp, zh_s(ji, jl) - zdh * rhoi * r1_rhos)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE ice_var_itd
  SUBROUTINE ice_var_itd2(zhti, zhts, zati, zh_i, zh_s, za_i)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jl, jl1, jl2
    INTEGER :: idim, icat
    REAL(KIND = wp), PARAMETER :: ztrans = 0.25_wp
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: zhti, zhts, zati
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: zh_i, zh_s, za_i
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jlfil, jlfil2
    INTEGER, DIMENSION(:), ALLOCATABLE :: jlmax, jlmin
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('ice_var_itd2', 'r0', 0, 0)
    idim = SIZE(zhti, 1)
    icat = SIZE(zhti, 2)
    ALLOCATE(jlfil(idim, jpl), jlfil2(idim, jpl))
    ALLOCATE(jlmin(idim), jlmax(idim))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zh_i(1 : idim, 1 : jpl) = 0._wp
    zh_s(1 : idim, 1 : jpl) = 0._wp
    za_i(1 : idim, 1 : jpl) = 0._wp
    jlmax(:) = 0
    jlmin(:) = 999
    jlfil(:, :) = 0
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_var_itd2', 'r1', 0, 0)
    DO jl1 = 1, jpl
      DO jl2 = 1, icat
        DO ji = 1, idim
          IF (hi_max(jl1 - 1) <= zhti(ji, jl2) .AND. hi_max(jl1) > zhti(ji, jl2)) THEN
            zh_i(ji, jl1) = zhti(ji, jl2)
            zh_s(ji, jl1) = zhts(ji, jl2)
            za_i(ji, jl1) = zati(ji, jl2)
            jlmax(ji) = MAX(jlmax(ji), jl1)
            jlmin(ji) = MIN(jlmin(ji), jl1)
            jlfil(ji, jl1) = jl1
          END IF
        END DO
      END DO
    END DO
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    DO ji = 1, idim
      jl1 = jlmin(ji)
      jl2 = jlmax(ji)
      IF (jl1 > 1) THEN
        za_i(ji, jl1 - 1) = ztrans * za_i(ji, jl1)
        zh_i(ji, jl1 - 1) = hi_mean(jl1 - 1)
        za_i(ji, jl1) = (1._wp - ztrans) * za_i(ji, jl1)
      END IF
      IF (jl2 < jpl) THEN
        za_i(ji, jl2 + 1) = ztrans * za_i(ji, jl2)
        zh_i(ji, jl2 + 1) = hi_mean(jl2 + 1)
        za_i(ji, jl2) = (1._wp - ztrans) * za_i(ji, jl2)
      END IF
    END DO
    jlfil2(:, :) = jlfil(:, :)
    DO jl = 2, jpl - 1
      DO ji = 1, idim
        IF (jlfil(ji, jl - 1) /= 0 .AND. jlfil(ji, jl) == 0) THEN
          za_i(ji, jl) = ztrans * za_i(ji, jl - 1)
          zh_i(ji, jl) = hi_mean(jl)
          jlfil(ji, jl) = jl
          za_i(ji, jl - 1) = (1._wp - ztrans) * za_i(ji, jl - 1)
        END IF
      END DO
    END DO
    DO jl = jpl - 1, 2, - 1
      DO ji = 1, idim
        IF (jlfil2(ji, jl + 1) /= 0 .AND. jlfil2(ji, jl) == 0) THEN
          za_i(ji, jl) = za_i(ji, jl) + ztrans * za_i(ji, jl + 1)
          zh_i(ji, jl) = hi_mean(jl)
          jlfil2(ji, jl) = jl
          za_i(ji, jl + 1) = (1._wp - ztrans) * za_i(ji, jl + 1)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('ice_var_itd2', 'r2', 0, 0)
    DEALLOCATE(jlfil, jlfil2)
    DEALLOCATE(jlmin, jlmax)
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE ice_var_itd2
  SUBROUTINE ice_var_bv
    INTEGER :: ji, jj, jk, jl
    !$ACC KERNELS
    bv_i(:, :, :) = 0._wp
    !$ACC END KERNELS
    DO jl = 1, jpl
      !$ACC KERNELS
      DO jk = 1, nlay_i
        WHERE (t_i(:, :, jk, jl) < rt0 - epsi10)
          bv_i(:, :, jl) = bv_i(:, :, jl) - rTmlt * sz_i(:, :, jk, jl) * r1_nlay_i / (t_i(:, :, jk, jl) - rt0)
        END WHERE
      END DO
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    WHERE (vt_i(:, :) > epsi20)
      bvm_i(:, :) = SUM(bv_i(:, :, :) * v_i(:, :, :), dim = 3) / vt_i(:, :)
    ELSEWHERE
      bvm_i(:, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
  END SUBROUTINE ice_var_bv
  SUBROUTINE ice_var_enthalpy
    INTEGER :: ji, jk
    REAL(KIND = wp) :: ztmelts
    !$ACC KERNELS
    DO jk = 1, nlay_i
      DO ji = 1, npti
        ztmelts = - rTmlt * sz_i_1d(ji, jk)
        t_i_1d(ji, jk) = MIN(t_i_1d(ji, jk), ztmelts + rt0)
        e_i_1d(ji, jk) = rhoi * (rcpi * (ztmelts - (t_i_1d(ji, jk) - rt0)) + rLfus * (1._wp - ztmelts / (t_i_1d(ji, jk) - rt0)) - &
&rcp * ztmelts)
      END DO
    END DO
    DO jk = 1, nlay_s
      DO ji = 1, npti
        e_s_1d(ji, jk) = rhos * (rcpi * (rt0 - t_s_1d(ji, jk)) + rLfus)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE ice_var_enthalpy
  FUNCTION ice_var_sshdyn(pssh, psnwice_mass, psnwice_mass_b)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pssh
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: psnwice_mass
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: psnwice_mass_b
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ice_var_sshdyn
    REAL(KIND = wp) :: zintn, zintb
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsnwiceload
    INTEGER :: ji, jj ! CDe added
    ! TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    IF (ln_ice_embd) THEN
      !$ACC KERNELS
      zintn = REAL(nn_fsbc - 1) / REAL(nn_fsbc) * 0.5_wp
      zintb = REAL(nn_fsbc + 1) / REAL(nn_fsbc) * 0.5_wp
      zsnwiceload(:, :) = (zintn * psnwice_mass(:, :) + zintb * psnwice_mass_b(:, :)) * r1_rau0
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zsnwiceload(:, :) = 0.0_wp
      !$ACC END KERNELS
    END IF
!    CALL profile_psy_data0 % PreStart('ice_var_sshdyn', 'r0', 0, 0)
    !$ACC KERNELS ! CDe added
    !$ACC LOOP INDEPENDENT COLLAPSE(2) ! CDe added
    DO jj = 1, jpj
      DO ji = 1, jpi
        ice_var_sshdyn(ji, jj) = pssh(ji, jj) + zsnwiceload(ji, jj) ! CDe re-wrote with explicit do loops
      END DO
    END DO
    !$ACC END KERNELS
!    CALL profile_psy_data0 % PostEnd
  END FUNCTION ice_var_sshdyn
END MODULE icevar
