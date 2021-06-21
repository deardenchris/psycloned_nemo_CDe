MODULE icectl
  USE phycst
  USE oce
  USE dom_oce
  USE ice
  USE ice1D
  USE sbc_oce
  USE sbc_ice
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE timing
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_cons_hsm
  PUBLIC :: ice_cons_final
  PUBLIC :: ice_ctl
  PUBLIC :: ice_prt
  PUBLIC :: ice_prt3D
  CONTAINS
  SUBROUTINE ice_cons_hsm(icount, cd_routine, pdiag_v, pdiag_s, pdiag_t, pdiag_fv, pdiag_fs, pdiag_ft)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: icount
    CHARACTER(LEN = *), INTENT(IN) :: cd_routine
    REAL(KIND = wp), INTENT(INOUT) :: pdiag_v, pdiag_s, pdiag_t, pdiag_fv, pdiag_fs, pdiag_ft
    REAL(KIND = wp) :: zv, zs, zt, zfs, zfv, zft
    REAL(KIND = wp) :: zvmin, zamin, zamax, zeimin, zesmin, zsmin
    REAL(KIND = wp) :: zvtrp, zetrp
    REAL(KIND = wp) :: zarea, zv_sill, zs_sill, zt_sill
    REAL(KIND = wp), PARAMETER :: zconv = 1.E-9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_cons_hsm', 'r0', 0, 0)
    IF (icount == 0) THEN
      pdiag_fv = glob_sum('icectl', - (wfx_bog(:, :) + wfx_bom(:, :) + wfx_sum(:, :) + wfx_sni(:, :) + wfx_opw(:, :) + wfx_res(:, :) + wfx_dyn(:, :) + wfx_lam(:, :) + wfx_pnd(:, :) + wfx_snw_sni(:, :) + wfx_snw_sum(:, :) + wfx_snw_dyn(:, :) + wfx_snw_sub(:, :) + wfx_ice_sub(:, :) + wfx_spr(:, :)) * e1e2t(:, :)) * zconv
      pdiag_fs = glob_sum('icectl', (sfx_bri(:, :) + sfx_bog(:, :) + sfx_bom(:, :) + sfx_sum(:, :) + sfx_sni(:, :) + sfx_opw(:, :) + sfx_res(:, :) + sfx_dyn(:, :) + sfx_sub(:, :) + sfx_lam(:, :)) * e1e2t(:, :)) * zconv
      pdiag_ft = glob_sum('icectl', (hfx_sum(:, :) + hfx_bom(:, :) + hfx_bog(:, :) + hfx_dif(:, :) + hfx_opw(:, :) + hfx_snw(:, :) - hfx_thd(:, :) - hfx_dyn(:, :) - hfx_res(:, :) - hfx_sub(:, :) - hfx_spr(:, :)) * e1e2t(:, :)) * zconv
      pdiag_v = glob_sum('icectl', SUM(v_i * rhoi + v_s * rhos, dim = 3) * e1e2t * zconv)
      pdiag_s = glob_sum('icectl', SUM(sv_i * rhoi, dim = 3) * e1e2t * zconv)
      pdiag_t = glob_sum('icectl', (SUM(SUM(e_i(:, :, 1 : nlay_i, :), dim = 4), dim = 3) + SUM(SUM(e_s(:, :, 1 : nlay_s, :), dim = 4), dim = 3)) * e1e2t) * zconv
    ELSE IF (icount == 1) THEN
      zfv = glob_sum('icectl', - (wfx_bog(:, :) + wfx_bom(:, :) + wfx_sum(:, :) + wfx_sni(:, :) + wfx_opw(:, :) + wfx_res(:, :) + wfx_dyn(:, :) + wfx_lam(:, :) + wfx_pnd(:, :) + wfx_snw_sni(:, :) + wfx_snw_sum(:, :) + wfx_snw_dyn(:, :) + wfx_snw_sub(:, :) + wfx_ice_sub(:, :) + wfx_spr(:, :)) * e1e2t(:, :)) * zconv - pdiag_fv
      zfs = glob_sum('icectl', (sfx_bri(:, :) + sfx_bog(:, :) + sfx_bom(:, :) + sfx_sum(:, :) + sfx_sni(:, :) + sfx_opw(:, :) + sfx_res(:, :) + sfx_dyn(:, :) + sfx_sub(:, :) + sfx_lam(:, :)) * e1e2t(:, :)) * zconv - pdiag_fs
      zft = glob_sum('icectl', (hfx_sum(:, :) + hfx_bom(:, :) + hfx_bog(:, :) + hfx_dif(:, :) + hfx_opw(:, :) + hfx_snw(:, :) - hfx_thd(:, :) - hfx_dyn(:, :) - hfx_res(:, :) - hfx_sub(:, :) - hfx_spr(:, :)) * e1e2t(:, :)) * zconv - pdiag_ft
      zv = ((glob_sum('icectl', SUM(v_i * rhoi + v_s * rhos, dim = 3) * e1e2t) * zconv - pdiag_v) * r1_rdtice - zfv) * rday
      zs = ((glob_sum('icectl', SUM(sv_i * rhoi, dim = 3) * e1e2t) * zconv - pdiag_s) * r1_rdtice + zfs) * rday
      zt = (glob_sum('icectl', (SUM(SUM(e_i(:, :, 1 : nlay_i, :), dim = 4), dim = 3) + SUM(SUM(e_s(:, :, 1 : nlay_s, :), dim = 4), dim = 3)) * e1e2t) * zconv - pdiag_t) * r1_rdtice + zft
      zvtrp = glob_sum('icectl', (diag_trp_vi * rhoi + diag_trp_vs * rhos) * e1e2t) * zconv * rday
      zetrp = glob_sum('icectl', (diag_trp_ei + diag_trp_es) * e1e2t) * zconv
      zamax = glob_max('icectl', SUM(a_i, dim = 3))
      zvmin = glob_min('icectl', v_i)
      zamin = glob_min('icectl', a_i)
      zsmin = glob_min('icectl', sv_i)
      zeimin = glob_min('icectl', SUM(e_i, dim = 3))
      zesmin = glob_min('icectl', SUM(e_s, dim = 3))
      zarea = glob_sum('icectl', SUM(a_i + epsi10, dim = 3) * e1e2t) * zconv
      zv_sill = zarea * 2.5E-5
      zs_sill = zarea * 25.E-5
      zt_sill = zarea * 10.E-5
      IF (lwp) THEN
        IF (ABS(zv) > zv_sill) WRITE(numout, *) 'violation volume [Mt/day]     (', cd_routine, ') = ', zv
        IF (ABS(zs) > zs_sill) WRITE(numout, *) 'violation saline [psu*Mt/day] (', cd_routine, ') = ', zs
        IF (ABS(zt) > zt_sill) WRITE(numout, *) 'violation enthalpy [GW]       (', cd_routine, ') = ', zt
        IF (zamax > MAX(rn_amax_n, rn_amax_s) + epsi10 .AND. cd_routine /= 'icedyn_adv' .AND. cd_routine /= 'icedyn_rdgrft') WRITE(numout, *) 'violation a_i>amax            (', cd_routine, ') = ', zamax
        IF (zvmin < 0.) WRITE(numout, *) 'violation v_i<0  [m]          (', cd_routine, ') = ', zvmin
        IF (zamin < 0.) WRITE(numout, *) 'violation a_i<0               (', cd_routine, ') = ', zamin
        IF (zsmin < 0.) WRITE(numout, *) 'violation s_i<0               (', cd_routine, ') = ', zsmin
        IF (zeimin < 0.) WRITE(numout, *) 'violation e_i<0               (', cd_routine, ') = ', zeimin
        IF (zesmin < 0.) WRITE(numout, *) 'violation e_s<0               (', cd_routine, ') = ', zesmin
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_cons_hsm
  SUBROUTINE ice_cons_final(cd_routine)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cd_routine
    REAL(KIND = wp) :: zqmass, zhfx, zsfx, zvfx
    REAL(KIND = wp) :: zarea, zv_sill, zs_sill, zt_sill
    REAL(KIND = wp), PARAMETER :: zconv = 1.E-9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_cons_final', 'r0', 0, 0)
    zvfx = glob_sum('icectl', (wfx_ice + wfx_snw + wfx_spr + wfx_sub + diag_vice + diag_vsnw) * e1e2t) * zconv * rday
    zsfx = glob_sum('icectl', (sfx + diag_sice) * e1e2t) * zconv * rday
    zarea = glob_sum('icectl', SUM(a_i + epsi10, dim = 3) * e1e2t) * zconv
    zv_sill = zarea * 2.5E-5
    zs_sill = zarea * 25.E-5
    zt_sill = zarea * 10.E-5
    IF (lwp) THEN
      IF (ABS(zvfx) > zv_sill) WRITE(numout, *) 'violation vfx  [Mt/day]       (', cd_routine, ') = ', zvfx
      IF (ABS(zsfx) > zs_sill) WRITE(numout, *) 'violation sfx  [psu*Mt/day]   (', cd_routine, ') = ', zsfx
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_cons_final
  SUBROUTINE ice_ctl(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl
    INTEGER :: inb_altests
    INTEGER :: ialert_id
    REAL(KIND = wp) :: ztmelts
    CHARACTER(LEN = 30), DIMENSION(20) :: cl_alname
    INTEGER, DIMENSION(20) :: inb_alp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    inb_altests = 10
    inb_alp(:) = 0
    ialert_id = 2
    cl_alname(ialert_id) = ' Incompat vol and con         '
    DO jl = 1, jpl
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (v_i(ji, jj, jl) /= 0._wp .AND. a_i(ji, jj, jl) == 0._wp) THEN
            inb_alp(ialert_id) = inb_alp(ialert_id) + 1
          END IF
        END DO
      END DO
    END DO
    ialert_id = 3
    cl_alname(ialert_id) = ' Very thick ice               '
    jl = jpl
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (h_i(ji, jj, jl) > 50._wp) THEN
          inb_alp(ialert_id) = inb_alp(ialert_id) + 1
        END IF
      END DO
    END DO
    ialert_id = 4
    cl_alname(ialert_id) = ' Very fast ice               '
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (MAX(ABS(u_ice(ji, jj)), ABS(v_ice(ji, jj))) > 1.5 .AND. at_i(ji, jj) > 0._wp) THEN
          inb_alp(ialert_id) = inb_alp(ialert_id) + 1
        END IF
      END DO
    END DO
    ialert_id = 6
    cl_alname(ialert_id) = ' Ice on continents           '
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (tmask(ji, jj, 1) <= 0._wp .AND. at_i(ji, jj) > 0._wp) THEN
          inb_alp(ialert_id) = inb_alp(ialert_id) + 1
        END IF
      END DO
    END DO
    ialert_id = 7
    cl_alname(ialert_id) = ' Very fresh ice               '
    DO jl = 1, jpl
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (s_i(ji, jj, jl) < 0.1 .AND. a_i(ji, jj, jl) > 0._wp) THEN
            inb_alp(ialert_id) = inb_alp(ialert_id) + 1
          END IF
        END DO
      END DO
    END DO
    ialert_id = 9
    cl_alname(ialert_id) = ' Very old   ice               '
    DO jl = 1, jpl
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (((ABS(o_i(ji, jj, jl)) > rdt_ice) .OR. (ABS(o_i(ji, jj, jl)) < 0._wp)) .AND. (a_i(ji, jj, jl) > 0._wp)) THEN
            inb_alp(ialert_id) = inb_alp(ialert_id) + 1
          END IF
        END DO
      END DO
    END DO
    ialert_id = 5
    cl_alname(ialert_id) = ' High salt flux               '
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (ABS(sfx(ji, jj)) > 1.0E-2) THEN
          inb_alp(ialert_id) = inb_alp(ialert_id) + 1
        END IF
      END DO
    END DO
    ialert_id = 8
    cl_alname(ialert_id) = ' fnsolar very big             '
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (ABS(qns(ji, jj)) > 1500._wp .AND. at_i(ji, jj) > 0._wp) THEN
          inb_alp(ialert_id) = inb_alp(ialert_id) + 1
        END IF
      END DO
    END DO
    ialert_id = 10
    cl_alname(ialert_id) = ' Very warm ice                '
    inb_alp(ialert_id) = 0
    !$ACC END KERNELS
    DO jl = 1, jpl
      !$ACC KERNELS
      DO jk = 1, nlay_i
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            ztmelts = - rTmlt * sz_i(ji, jj, jk, jl) + rt0
            IF (t_i(ji, jj, jk, jl) >= ztmelts .AND. v_i(ji, jj, jl) > 1.E-10 .AND. a_i(ji, jj, jl) > 0._wp) THEN
              inb_alp(ialert_id) = inb_alp(ialert_id) + 1
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    CALL profile_psy_data0 % PreStart('ice_ctl', 'r0', 0, 0)
    IF (lk_mpp) THEN
      DO ialert_id = 1, inb_altests
        CALL mpp_sum('icectl', inb_alp(ialert_id))
      END DO
    END IF
    IF (lwp) THEN
      ialert_id = 1
      cl_alname(ialert_id) = ' NO alerte 1      '
      WRITE(numout, *) ' time step ', kt
      WRITE(numout, *) ' All alerts at the end of ice model '
      DO ialert_id = 1, inb_altests
        WRITE(numout, *) ialert_id, cl_alname(ialert_id) // ' : ', inb_alp(ialert_id), ' times ! '
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_ctl
  SUBROUTINE ice_prt(kt, ki, kj, kn, cd1)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: ki, kj, kn
    CHARACTER(LEN = *), INTENT(IN) :: cd1
    INTEGER :: jl, ji, jj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_prt', 'r0', 0, 0)
    DO ji = mi0(ki), mi1(ki)
      DO jj = mj0(kj), mj1(kj)
        WRITE(numout, *) ' time step ', kt, ' ', cd1
        IF (kn == 1 .OR. kn == - 1) THEN
          WRITE(numout, *) ' ice_prt - Point : ', ji, jj
          WRITE(numout, *) ' ~~~~~~~~~~~~~~ '
          WRITE(numout, *) ' Simple state '
          WRITE(numout, *) ' masks s,u,v   : ', tmask(ji, jj, 1), umask(ji, jj, 1), vmask(ji, jj, 1)
          WRITE(numout, *) ' lat - long    : ', gphit(ji, jj), glamt(ji, jj)
          WRITE(numout, *) ' - Ice drift   '
          WRITE(numout, *) '   ~~~~~~~~~~~ '
          WRITE(numout, *) ' u_ice(i-1,j)  : ', u_ice(ji - 1, jj)
          WRITE(numout, *) ' u_ice(i  ,j)  : ', u_ice(ji, jj)
          WRITE(numout, *) ' v_ice(i  ,j-1): ', v_ice(ji, jj - 1)
          WRITE(numout, *) ' v_ice(i  ,j)  : ', v_ice(ji, jj)
          WRITE(numout, *) ' strength      : ', strength(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' - Cell values '
          WRITE(numout, *) '   ~~~~~~~~~~~ '
          WRITE(numout, *) ' cell area     : ', e1e2t(ji, jj)
          WRITE(numout, *) ' at_i          : ', at_i(ji, jj)
          WRITE(numout, *) ' vt_i          : ', vt_i(ji, jj)
          WRITE(numout, *) ' vt_s          : ', vt_s(ji, jj)
          DO jl = 1, jpl
            WRITE(numout, *) ' - Category (', jl, ')'
            WRITE(numout, *) ' a_i           : ', a_i(ji, jj, jl)
            WRITE(numout, *) ' h_i           : ', h_i(ji, jj, jl)
            WRITE(numout, *) ' h_s           : ', h_s(ji, jj, jl)
            WRITE(numout, *) ' v_i           : ', v_i(ji, jj, jl)
            WRITE(numout, *) ' v_s           : ', v_s(ji, jj, jl)
            WRITE(numout, *) ' e_s           : ', e_s(ji, jj, 1 : nlay_s, jl)
            WRITE(numout, *) ' e_i           : ', e_i(ji, jj, 1 : nlay_i, jl)
            WRITE(numout, *) ' t_su          : ', t_su(ji, jj, jl)
            WRITE(numout, *) ' t_snow        : ', t_s(ji, jj, 1 : nlay_s, jl)
            WRITE(numout, *) ' t_i           : ', t_i(ji, jj, 1 : nlay_i, jl)
            WRITE(numout, *) ' s_i           : ', s_i(ji, jj, jl)
            WRITE(numout, *) ' sv_i          : ', sv_i(ji, jj, jl)
            WRITE(numout, *)
          END DO
        END IF
        IF (kn == - 1) THEN
          WRITE(numout, *) ' Mechanical Check ************** '
          WRITE(numout, *) ' Check what means ice divergence '
          WRITE(numout, *) ' Total ice concentration ', at_i(ji, jj)
          WRITE(numout, *) ' Total lead fraction     ', ato_i(ji, jj)
          WRITE(numout, *) ' Sum of both             ', ato_i(ji, jj) + at_i(ji, jj)
          WRITE(numout, *) ' Sum of both minus 1     ', ato_i(ji, jj) + at_i(ji, jj) - 1.00
        END IF
        IF (kn .EQ. 2) THEN
          WRITE(numout, *) ' ice_prt - Point : ', ji, jj
          WRITE(numout, *) ' ~~~~~~~~~~~~~~ '
          WRITE(numout, *) ' Exhaustive state '
          WRITE(numout, *) ' lat - long ', gphit(ji, jj), glamt(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' - Cell values '
          WRITE(numout, *) '   ~~~~~~~~~~~ '
          WRITE(numout, *) ' cell area     : ', e1e2t(ji, jj)
          WRITE(numout, *) ' at_i          : ', at_i(ji, jj)
          WRITE(numout, *) ' vt_i          : ', vt_i(ji, jj)
          WRITE(numout, *) ' vt_s          : ', vt_s(ji, jj)
          WRITE(numout, *) ' u_ice(i-1,j)  : ', u_ice(ji - 1, jj)
          WRITE(numout, *) ' u_ice(i  ,j)  : ', u_ice(ji, jj)
          WRITE(numout, *) ' v_ice(i  ,j-1): ', v_ice(ji, jj - 1)
          WRITE(numout, *) ' v_ice(i  ,j)  : ', v_ice(ji, jj)
          WRITE(numout, *) ' strength      : ', strength(ji, jj)
          WRITE(numout, *) ' u_ice_b       : ', u_ice_b(ji, jj), ' v_ice_b       : ', v_ice_b(ji, jj)
          WRITE(numout, *)
          DO jl = 1, jpl
            WRITE(numout, *) ' - Category (', jl, ')'
            WRITE(numout, *) '   ~~~~~~~~         '
            WRITE(numout, *) ' h_i        : ', h_i(ji, jj, jl), ' h_s        : ', h_s(ji, jj, jl)
            WRITE(numout, *) ' t_i        : ', t_i(ji, jj, 1 : nlay_i, jl)
            WRITE(numout, *) ' t_su       : ', t_su(ji, jj, jl), ' t_s        : ', t_s(ji, jj, 1 : nlay_s, jl)
            WRITE(numout, *) ' s_i        : ', s_i(ji, jj, jl), ' o_i        : ', o_i(ji, jj, jl)
            WRITE(numout, *) ' a_i        : ', a_i(ji, jj, jl), ' a_i_b      : ', a_i_b(ji, jj, jl)
            WRITE(numout, *) ' v_i        : ', v_i(ji, jj, jl), ' v_i_b      : ', v_i_b(ji, jj, jl)
            WRITE(numout, *) ' v_s        : ', v_s(ji, jj, jl), ' v_s_b      : ', v_s_b(ji, jj, jl)
            WRITE(numout, *) ' e_i1       : ', e_i(ji, jj, 1, jl), ' ei1        : ', e_i_b(ji, jj, 1, jl)
            WRITE(numout, *) ' e_i2       : ', e_i(ji, jj, 2, jl), ' ei2_b      : ', e_i_b(ji, jj, 2, jl)
            WRITE(numout, *) ' e_snow     : ', e_s(ji, jj, 1, jl), ' e_snow_b   : ', e_s_b(ji, jj, 1, jl)
            WRITE(numout, *) ' sv_i       : ', sv_i(ji, jj, jl), ' sv_i_b     : ', sv_i_b(ji, jj, jl)
            WRITE(numout, *) ' oa_i       : ', oa_i(ji, jj, jl), ' oa_i_b     : ', oa_i_b(ji, jj, jl)
          END DO
          WRITE(numout, *)
          WRITE(numout, *) ' - Heat / FW fluxes '
          WRITE(numout, *) '   ~~~~~~~~~~~~~~~~ '
          WRITE(numout, *) ' - Heat fluxes in and out the ice ***'
          WRITE(numout, *) ' qsr_ini       : ', (1._wp - at_i_b(ji, jj)) * qsr(ji, jj) + SUM(a_i_b(ji, jj, :) * qsr_ice(ji, jj, :))
          WRITE(numout, *) ' qns_ini       : ', (1._wp - at_i_b(ji, jj)) * qns(ji, jj) + SUM(a_i_b(ji, jj, :) * qns_ice(ji, jj, :))
          WRITE(numout, *)
          WRITE(numout, *)
          WRITE(numout, *) ' sst        : ', sst_m(ji, jj)
          WRITE(numout, *) ' sss        : ', sss_m(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' - Stresses '
          WRITE(numout, *) '   ~~~~~~~~ '
          WRITE(numout, *) ' utau_ice   : ', utau_ice(ji, jj)
          WRITE(numout, *) ' vtau_ice   : ', vtau_ice(ji, jj)
          WRITE(numout, *) ' utau       : ', utau(ji, jj)
          WRITE(numout, *) ' vtau       : ', vtau(ji, jj)
        END IF
        IF (kn .EQ. 3) THEN
          WRITE(numout, *) ' ice_prt - Point : ', ji, jj
          WRITE(numout, *) ' ~~~~~~~~~~~~~~ '
          WRITE(numout, *) ' - Salt / Heat Fluxes '
          WRITE(numout, *) '   ~~~~~~~~~~~~~~~~ '
          WRITE(numout, *) ' lat - long ', gphit(ji, jj), glamt(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' - Heat fluxes at bottom interface ***'
          WRITE(numout, *) ' qsr       : ', qsr(ji, jj)
          WRITE(numout, *) ' qns       : ', qns(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' hfx_mass     : ', hfx_thd(ji, jj) + hfx_dyn(ji, jj) + hfx_snw(ji, jj) + hfx_res(ji, jj)
          WRITE(numout, *) ' qt_atm_oi    : ', qt_atm_oi(ji, jj)
          WRITE(numout, *) ' qt_oce_ai    : ', qt_oce_ai(ji, jj)
          WRITE(numout, *) ' dhc          : ', diag_heat(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' hfx_dyn      : ', hfx_dyn(ji, jj)
          WRITE(numout, *) ' hfx_thd      : ', hfx_thd(ji, jj)
          WRITE(numout, *) ' hfx_res      : ', hfx_res(ji, jj)
          WRITE(numout, *) ' qsb_ice_bot  : ', qsb_ice_bot(ji, jj)
          WRITE(numout, *) ' qlead        : ', qlead(ji, jj) * r1_rdtice
          WRITE(numout, *)
          WRITE(numout, *) ' - Salt fluxes at bottom interface ***'
          WRITE(numout, *) ' emp       : ', emp(ji, jj)
          WRITE(numout, *) ' sfx       : ', sfx(ji, jj)
          WRITE(numout, *) ' sfx_res   : ', sfx_res(ji, jj)
          WRITE(numout, *) ' sfx_bri   : ', sfx_bri(ji, jj)
          WRITE(numout, *) ' sfx_dyn   : ', sfx_dyn(ji, jj)
          WRITE(numout, *)
          WRITE(numout, *) ' - Momentum fluxes '
          WRITE(numout, *) ' utau      : ', utau(ji, jj)
          WRITE(numout, *) ' vtau      : ', vtau(ji, jj)
        END IF
        WRITE(numout, *) ' '
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_prt
  SUBROUTINE ice_prt3D(cd_routine)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cd_routine
    INTEGER :: jk, jl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_prt3d', 'r0', 0, 0)
    CALL prt_ctl_info(' ========== ')
    CALL prt_ctl_info(cd_routine)
    CALL prt_ctl_info(' ========== ')
    CALL prt_ctl_info(' - Cell values : ')
    CALL prt_ctl_info('   ~~~~~~~~~~~~~ ')
    CALL prt_ctl(tab2d_1 = e1e2t, clinfo1 = ' cell area   :')
    CALL prt_ctl(tab2d_1 = at_i, clinfo1 = ' at_i        :')
    CALL prt_ctl(tab2d_1 = ato_i, clinfo1 = ' ato_i       :')
    CALL prt_ctl(tab2d_1 = vt_i, clinfo1 = ' vt_i        :')
    CALL prt_ctl(tab2d_1 = vt_s, clinfo1 = ' vt_s        :')
    CALL prt_ctl(tab2d_1 = divu_i, clinfo1 = ' divu_i      :')
    CALL prt_ctl(tab2d_1 = delta_i, clinfo1 = ' delta_i     :')
    CALL prt_ctl(tab2d_1 = stress1_i, clinfo1 = ' stress1_i   :')
    CALL prt_ctl(tab2d_1 = stress2_i, clinfo1 = ' stress2_i   :')
    CALL prt_ctl(tab2d_1 = stress12_i, clinfo1 = ' stress12_i  :')
    CALL prt_ctl(tab2d_1 = strength, clinfo1 = ' strength    :')
    CALL prt_ctl(tab2d_1 = delta_i, clinfo1 = ' delta_i     :')
    CALL prt_ctl(tab2d_1 = u_ice, clinfo1 = ' u_ice       :', tab2d_2 = v_ice, clinfo2 = ' v_ice       :')
    DO jl = 1, jpl
      CALL prt_ctl_info(' ')
      CALL prt_ctl_info(' - Category : ', ivar1 = jl)
      CALL prt_ctl_info('   ~~~~~~~~~~')
      CALL prt_ctl(tab2d_1 = h_i(:, :, jl), clinfo1 = ' h_i         : ')
      CALL prt_ctl(tab2d_1 = h_s(:, :, jl), clinfo1 = ' h_s         : ')
      CALL prt_ctl(tab2d_1 = t_su(:, :, jl), clinfo1 = ' t_su        : ')
      CALL prt_ctl(tab2d_1 = t_s(:, :, 1, jl), clinfo1 = ' t_snow      : ')
      CALL prt_ctl(tab2d_1 = s_i(:, :, jl), clinfo1 = ' s_i         : ')
      CALL prt_ctl(tab2d_1 = o_i(:, :, jl), clinfo1 = ' o_i         : ')
      CALL prt_ctl(tab2d_1 = a_i(:, :, jl), clinfo1 = ' a_i         : ')
      CALL prt_ctl(tab2d_1 = v_i(:, :, jl), clinfo1 = ' v_i         : ')
      CALL prt_ctl(tab2d_1 = v_s(:, :, jl), clinfo1 = ' v_s         : ')
      CALL prt_ctl(tab2d_1 = e_i(:, :, 1, jl), clinfo1 = ' e_i1        : ')
      CALL prt_ctl(tab2d_1 = e_s(:, :, 1, jl), clinfo1 = ' e_snow      : ')
      CALL prt_ctl(tab2d_1 = sv_i(:, :, jl), clinfo1 = ' sv_i        : ')
      CALL prt_ctl(tab2d_1 = oa_i(:, :, jl), clinfo1 = ' oa_i        : ')
      DO jk = 1, nlay_i
        CALL prt_ctl_info(' - Layer : ', ivar1 = jk)
        CALL prt_ctl(tab2d_1 = t_i(:, :, jk, jl), clinfo1 = ' t_i       : ')
      END DO
    END DO
    CALL prt_ctl_info(' ')
    CALL prt_ctl_info(' - Heat / FW fluxes : ')
    CALL prt_ctl_info('   ~~~~~~~~~~~~~~~~~~ ')
    CALL prt_ctl(tab2d_1 = sst_m, clinfo1 = ' sst   : ', tab2d_2 = sss_m, clinfo2 = ' sss       : ')
    CALL prt_ctl(tab2d_1 = qsr, clinfo1 = ' qsr   : ', tab2d_2 = qns, clinfo2 = ' qns       : ')
    CALL prt_ctl(tab2d_1 = emp, clinfo1 = ' emp   : ', tab2d_2 = sfx, clinfo2 = ' sfx       : ')
    CALL prt_ctl_info(' ')
    CALL prt_ctl_info(' - Stresses : ')
    CALL prt_ctl_info('   ~~~~~~~~~~ ')
    CALL prt_ctl(tab2d_1 = utau, clinfo1 = ' utau      : ', tab2d_2 = vtau, clinfo2 = ' vtau      : ')
    CALL prt_ctl(tab2d_1 = utau_ice, clinfo1 = ' utau_ice  : ', tab2d_2 = vtau_ice, clinfo2 = ' vtau_ice  : ')
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_prt3D
END MODULE icectl