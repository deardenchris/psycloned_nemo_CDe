MODULE icedyn_rdgrft
  USE dom_oce
  USE phycst
  USE sbc_oce, ONLY: sss_m, sst_m
  USE ice1D
  USE ice
  USE icetab
  USE icevar
  USE icectl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_rdgrft
  PUBLIC :: ice_dyn_rdgrft_init
  PUBLIC :: ice_strength
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:) :: closing_net
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:) :: opning
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:) :: closing_gross
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: apartf
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: hrmin
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: hrmax
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: hraft
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: hi_hrdg
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: aridge
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: araft
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ze_i_2d
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ze_s_2d
  REAL(KIND = wp), PARAMETER :: hrdg_hi_min = 1.1_wp
  REAL(KIND = wp), PARAMETER :: hi_hrft = 0.5_wp
  LOGICAL :: ln_str_H79
  REAL(KIND = wp) :: rn_pstar
  REAL(KIND = wp) :: rn_csrdg
  LOGICAL :: ln_partf_lin
  REAL(KIND = wp) :: rn_gstar
  LOGICAL :: ln_partf_exp
  REAL(KIND = wp) :: rn_astar
  LOGICAL :: ln_ridging
  REAL(KIND = wp) :: rn_hstar
  REAL(KIND = wp) :: rn_porordg
  REAL(KIND = wp) :: rn_fsnwrdg
  REAL(KIND = wp) :: rn_fpndrdg
  LOGICAL :: ln_rafting
  REAL(KIND = wp) :: rn_hraft
  REAL(KIND = wp) :: rn_craft
  REAL(KIND = wp) :: rn_fsnwrft
  REAL(KIND = wp) :: rn_fpndrft
  CONTAINS
  INTEGER FUNCTION ice_dyn_rdgrft_alloc()
    ALLOCATE(closing_net(jpij), opning(jpij), closing_gross(jpij), apartf(jpij, 0 : jpl), hrmin(jpij, jpl), hraft(jpij, jpl), aridge(jpij, jpl), hrmax(jpij, jpl), hi_hrdg(jpij, jpl), araft(jpij, jpl), ze_i_2d(jpij, nlay_i, jpl), ze_s_2d(jpij, nlay_s, jpl), STAT = ice_dyn_rdgrft_alloc)
    CALL mpp_sum('icedyn_rdgrft', ice_dyn_rdgrft_alloc)
    IF (ice_dyn_rdgrft_alloc /= 0) CALL ctl_stop('STOP', 'ice_dyn_rdgrft_alloc: failed to allocate arrays')
  END FUNCTION ice_dyn_rdgrft_alloc
  SUBROUTINE ice_dyn_rdgrft(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl
    INTEGER :: iter, iterate_ridging
    INTEGER :: ipti
    REAL(KIND = wp) :: zfac
    INTEGER, DIMENSION(jpij) :: iptidx
    REAL(KIND = wp), DIMENSION(jpij) :: zdivu_adv
    REAL(KIND = wp), DIMENSION(jpij) :: zdivu, zdelt
    INTEGER, PARAMETER :: jp_itermax = 20
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('ice_dyn_rdgrft', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('icedyn_rdgrft')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icedyn_rdgrft', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'ice_dyn_rdgrft: ice ridging and rafting'
      IF (lwp) WRITE(numout, *) '~~~~~~~~~~~~~~'
    END IF
    at_i(:, :) = SUM(a_i, dim = 3)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    npti = 0
    nptidx(:) = 0
    ipti = 0
    iptidx(:) = 0
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (at_i(ji, jj) > epsi10) THEN
          npti = npti + 1
          nptidx(npti) = (jj - 1) * jpi + ji
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    IF (npti > 0) THEN
      CALL profile_psy_data1 % PreStart('ice_dyn_rdgrft', 'r1', 0, 0)
      CALL tab_2d_1d(npti, nptidx(1 : npti), zdivu(1 : npti), divu_i)
      CALL tab_2d_1d(npti, nptidx(1 : npti), zdelt(1 : npti), delta_i)
      CALL tab_3d_2d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i)
      CALL tab_3d_2d(npti, nptidx(1 : npti), v_i_2d(1 : npti, 1 : jpl), v_i)
      CALL tab_2d_1d(npti, nptidx(1 : npti), ato_i_1d(1 : npti), ato_i)
      DO ji = 1, npti
        closing_net(ji) = rn_csrdg * 0.5_wp * (zdelt(ji) - ABS(zdivu(ji))) - MIN(zdivu(ji), 0._wp)
        IF (ln_adv_Pra) THEN
          zdivu_adv(ji) = (1._wp - ato_i_1d(ji) - SUM(a_i_2d(ji, :))) * r1_rdtice
        ELSE IF (ln_adv_UMx) THEN
          zdivu_adv(ji) = zdivu(ji)
        END IF
        IF (zdivu_adv(ji) < 0._wp) closing_net(ji) = MAX(closing_net(ji), - zdivu_adv(ji))
        opning(ji) = closing_net(ji) + zdivu_adv(ji)
      END DO
      CALL rdgrft_prep(a_i_2d, v_i_2d, ato_i_1d, closing_net)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      DO ji = 1, npti
        IF (SUM(apartf(ji, 1 : jpl)) > 0._wp .AND. closing_gross(ji) > 0._wp) THEN
          ipti = ipti + 1
          iptidx(ipti) = nptidx(ji)
          a_i_2d(ipti, :) = a_i_2d(ji, :)
          v_i_2d(ipti, :) = v_i_2d(ji, :)
          ato_i_1d(ipti) = ato_i_1d(ji)
          closing_net(ipti) = closing_net(ji)
          zdivu_adv(ipti) = zdivu_adv(ji)
          opning(ipti) = opning(ji)
        END IF
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    nptidx(:) = iptidx(:)
    npti = ipti
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('ice_dyn_rdgrft', 'r2', 0, 0)
    IF (npti > 0) THEN
      CALL ice_dyn_1d2d(1)
      iter = 1
      iterate_ridging = 1
      DO WHILE (iterate_ridging > 0 .AND. iter < jp_itermax)
        CALL rdgrft_prep(a_i_2d, v_i_2d, ato_i_1d, closing_net)
        CALL rdgrft_shift
        iterate_ridging = 0
        DO ji = 1, npti
          zfac = 1._wp - (ato_i_1d(ji) + SUM(a_i_2d(ji, :)))
          IF (ABS(zfac) < epsi10) THEN
            closing_net(ji) = 0._wp
            opning(ji) = 0._wp
            ato_i_1d(ji) = MAX(0._wp, 1._wp - SUM(a_i_2d(ji, :)))
          ELSE
            iterate_ridging = 1
            zdivu_adv(ji) = zfac * r1_rdtice
            closing_net(ji) = MAX(0._wp, - zdivu_adv(ji))
            opning(ji) = MAX(0._wp, zdivu_adv(ji))
          END IF
        END DO
        iter = iter + 1
        IF (iter > jp_itermax) CALL ctl_stop('STOP', 'icedyn_rdgrft: non-converging ridging scheme')
      END DO
      CALL ice_dyn_1d2d(2)
    END IF
    CALL ice_var_agg(1)
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'icedyn_rdgrft', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (ln_ctl) CALL ice_prt3D('icedyn_rdgrft')
    IF (ln_timing) CALL timing_stop('icedyn_rdgrft')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE ice_dyn_rdgrft
  SUBROUTINE rdgrft_prep(pa_i, pv_i, pato_i, pclosing_net)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pato_i, pclosing_net
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pa_i, pv_i
    INTEGER :: ji, jl
    REAL(KIND = wp) :: z1_gstar, z1_astar, zhmean, zfac
    REAL(KIND = wp), DIMENSION(jpij) :: zasum, z1_asum, zaksum
    REAL(KIND = wp), DIMENSION(jpij, jpl) :: zhi
    REAL(KIND = wp), DIMENSION(jpij, - 1 : jpl) :: zGsum
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('rdgrft_prep', 'r0', 0, 0)
    z1_gstar = 1._wp / rn_gstar
    z1_astar = 1._wp / rn_astar
    WHERE (pa_i(1 : npti, :) > epsi20)
      zhi(1 : npti, :) = pv_i(1 : npti, :) / pa_i(1 : npti, :)
    ELSEWHERE
      zhi(1 : npti, :) = 0._wp
    END WHERE
    zasum(1 : npti) = pato_i(1 : npti) + SUM(pa_i(1 : npti, :), dim = 2)
    WHERE (zasum(1 : npti) > epsi20)
      z1_asum(1 : npti) = 1._wp / zasum(1 : npti)
    ELSEWHERE
      z1_asum(1 : npti) = 0._wp
    END WHERE
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zGsum(1 : npti, - 1) = 0._wp
    zGsum(1 : npti, 0) = pato_i(1 : npti) * z1_asum(1 : npti)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('rdgrft_prep', 'r1', 0, 0)
    DO jl = 1, jpl
      zGsum(1 : npti, jl) = (pato_i(1 : npti) + SUM(pa_i(1 : npti, 1 : jl), dim = 2)) * z1_asum(1 : npti)
    END DO
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    IF (ln_partf_lin) THEN
      DO jl = 0, jpl
        DO ji = 1, npti
          IF (zGsum(ji, jl) < rn_gstar) THEN
            apartf(ji, jl) = z1_gstar * (zGsum(ji, jl) - zGsum(ji, jl - 1)) * (2._wp - (zGsum(ji, jl - 1) + zGsum(ji, jl)) * z1_gstar)
          ELSE IF (zGsum(ji, jl - 1) < rn_gstar) THEN
            apartf(ji, jl) = z1_gstar * (rn_gstar - zGsum(ji, jl - 1)) * (2._wp - (zGsum(ji, jl - 1) + rn_gstar) * z1_gstar)
          ELSE
            apartf(ji, jl) = 0._wp
          END IF
        END DO
      END DO
    ELSE IF (ln_partf_exp) THEN
      zfac = 1._wp / (1._wp - EXP(- z1_astar))
      DO jl = - 1, jpl
        DO ji = 1, npti
          zGsum(ji, jl) = EXP(- zGsum(ji, jl) * z1_astar) * zfac
        END DO
      END DO
      DO jl = 0, jpl
        DO ji = 1, npti
          apartf(ji, jl) = zGsum(ji, jl - 1) - zGsum(ji, jl)
        END DO
      END DO
    END IF
    !$ACC END KERNELS
    IF (ln_rafting .AND. ln_ridging) THEN
      CALL profile_psy_data2 % PreStart('rdgrft_prep', 'r2', 0, 0)
      DO jl = 1, jpl
        DO ji = 1, npti
          aridge(ji, jl) = (1._wp + TANH(rn_craft * (zhi(ji, jl) - rn_hraft))) * 0.5_wp * apartf(ji, jl)
          araft(ji, jl) = apartf(ji, jl) - aridge(ji, jl)
        END DO
      END DO
      CALL profile_psy_data2 % PostEnd
    ELSE IF (ln_ridging .AND. .NOT. ln_rafting) THEN
      !$ACC KERNELS
      DO jl = 1, jpl
        DO ji = 1, npti
          aridge(ji, jl) = apartf(ji, jl)
          araft(ji, jl) = 0._wp
        END DO
      END DO
      !$ACC END KERNELS
      ! !$ACC KERNELS ! CDe compiler doesn't like this
    ELSE IF (ln_rafting .AND. .NOT. ln_ridging) THEN
      ! !$ACC END KERNELS
      !$ACC KERNELS ! CDe added
      DO jl = 1, jpl
        DO ji = 1, npti
          aridge(ji, jl) = 0._wp
          araft(ji, jl) = apartf(ji, jl)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS ! CDe added      
      DO jl = 1, jpl
        DO ji = 1, npti
          aridge(ji, jl) = 0._wp
          araft(ji, jl) = 0._wp
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    zfac = 1._wp / hi_hrft
    zaksum(1 : npti) = apartf(1 : npti, 0)
    DO jl = 1, jpl
      DO ji = 1, npti
        IF (apartf(ji, jl) > 0._wp) THEN
          zhmean = MAX(SQRT(rn_hstar * zhi(ji, jl)), zhi(ji, jl) * hrdg_hi_min)
          hrmin(ji, jl) = MIN(2._wp * zhi(ji, jl), 0.5_wp * (zhmean + zhi(ji, jl)))
          hrmax(ji, jl) = 2._wp * zhmean - hrmin(ji, jl)
          hraft(ji, jl) = zhi(ji, jl) * zfac
          hi_hrdg(ji, jl) = zhi(ji, jl) / MAX(zhmean, epsi20)
          zaksum(ji) = zaksum(ji) + aridge(ji, jl) * (1._wp - hi_hrdg(ji, jl)) + araft(ji, jl) * (1._wp - hi_hrft)
        ELSE
          hrmin(ji, jl) = 0._wp
          hrmax(ji, jl) = 0._wp
          hraft(ji, jl) = 0._wp
          hi_hrdg(ji, jl) = 1._wp
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (zaksum(1 : npti) > epsi20)
      closing_gross(1 : npti) = pclosing_net(1 : npti) / zaksum(1 : npti)
    ELSEWHERE
      closing_gross(1 : npti) = 0._wp
    END WHERE
    !$ACC KERNELS
    DO jl = 1, jpl
      DO ji = 1, npti
        zfac = apartf(ji, jl) * closing_gross(ji) * rdt_ice
        IF (zfac > pa_i(ji, jl) .AND. apartf(ji, jl) /= 0._wp) THEN
          closing_gross(ji) = pa_i(ji, jl) / apartf(ji, jl) * r1_rdtice
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('rdgrft_prep', 'r3', 0, 0)
    DO ji = 1, npti
      zfac = pato_i(ji) + (opning(ji) - apartf(ji, 0) * closing_gross(ji)) * rdt_ice
      IF (zfac < 0._wp) THEN
        opning(ji) = apartf(ji, 0) * closing_gross(ji) - pato_i(ji) * r1_rdtice
      ELSE IF (zfac > zasum(ji)) THEN
        opning(ji) = apartf(ji, 0) * closing_gross(ji) + (zasum(ji) - pato_i(ji)) * r1_rdtice
      END IF
    END DO
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE rdgrft_prep
  SUBROUTINE rdgrft_shift
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jl, jl1, jl2, jk
    REAL(KIND = wp) :: hL, hR, farea
    REAL(KIND = wp) :: vsw
    REAL(KIND = wp) :: afrdg, afrft
    REAL(KIND = wp) :: airdg1, oirdg1, aprdg1, virdg1, sirdg1
    REAL(KIND = wp) :: airft1, oirft1, aprft1
    REAL(KIND = wp), DIMENSION(jpij) :: airdg2, oirdg2, aprdg2, virdg2, sirdg2, vsrdg, vprdg
    REAL(KIND = wp), DIMENSION(jpij) :: airft2, oirft2, aprft2, virft, sirft, vsrft, vprft
    REAL(KIND = wp), DIMENSION(jpij) :: ersw
    REAL(KIND = wp), DIMENSION(jpij) :: zswitch, fvol
    REAL(KIND = wp), DIMENSION(jpij) :: z1_ai
    REAL(KIND = wp), DIMENSION(jpij) :: zvti
    REAL(KIND = wp), DIMENSION(jpij, nlay_s) :: esrft
    REAL(KIND = wp), DIMENSION(jpij, nlay_i) :: eirft
    REAL(KIND = wp), DIMENSION(jpij, nlay_s) :: esrdg
    REAL(KIND = wp), DIMENSION(jpij, nlay_i) :: eirdg
    INTEGER, DIMENSION(jpij) :: itest_rdg, itest_rft
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('rdgrft_shift', 'r0', 0, 0)
    zvti(1 : npti) = SUM(v_i_2d(1 : npti, :), dim = 2)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO ji = 1, npti
      ato_i_1d(ji) = MAX(0._wp, ato_i_1d(ji) + (opning(ji) - apartf(ji, 0) * closing_gross(ji)) * rdt_ice)
    END DO
    !$ACC END KERNELS
    DO jl1 = 1, jpl
      CALL profile_psy_data1 % PreStart('rdgrft_shift', 'r1', 0, 0)
      CALL tab_2d_1d(npti, nptidx(1 : npti), s_i_1d(1 : npti), s_i(:, :, jl1))
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      DO ji = 1, npti
        IF (apartf(ji, jl1) > 0._wp .AND. closing_gross(ji) > 0._wp) THEN
          IF (a_i_2d(ji, jl1) > epsi20) THEN
            z1_ai(ji) = 1._wp / a_i_2d(ji, jl1)
          ELSE
            z1_ai(ji) = 0._wp
          END IF
          airdg1 = aridge(ji, jl1) * closing_gross(ji) * rdt_ice
          airft1 = araft(ji, jl1) * closing_gross(ji) * rdt_ice
          airdg2(ji) = airdg1 * hi_hrdg(ji, jl1)
          airft2(ji) = airft1 * hi_hrft
          afrdg = airdg1 * z1_ai(ji)
          afrft = airft1 * z1_ai(ji)
          IF (zvti(ji) <= 10.) THEN
            vsw = v_i_2d(ji, jl1) * afrdg * rn_porordg
          ELSE IF (zvti(ji) >= 20.) THEN
            vsw = 0._wp
          ELSE
            vsw = v_i_2d(ji, jl1) * afrdg * rn_porordg * MAX(0._wp, 2._wp - 0.1_wp * zvti(ji))
          END IF
          ersw(ji) = - rhoi * vsw * rcp * sst_1d(ji)
          virdg1 = v_i_2d(ji, jl1) * afrdg
          virdg2(ji) = v_i_2d(ji, jl1) * afrdg + vsw
          vsrdg(ji) = v_s_2d(ji, jl1) * afrdg
          sirdg1 = sv_i_2d(ji, jl1) * afrdg
          sirdg2(ji) = sv_i_2d(ji, jl1) * afrdg + vsw * sss_1d(ji)
          oirdg1 = oa_i_2d(ji, jl1) * afrdg
          oirdg2(ji) = oa_i_2d(ji, jl1) * afrdg * hi_hrdg(ji, jl1)
          virft(ji) = v_i_2d(ji, jl1) * afrft
          vsrft(ji) = v_s_2d(ji, jl1) * afrft
          sirft(ji) = sv_i_2d(ji, jl1) * afrft
          oirft1 = oa_i_2d(ji, jl1) * afrft
          oirft2(ji) = oa_i_2d(ji, jl1) * afrft * hi_hrft
          IF (ln_pnd_H12) THEN
            aprdg1 = a_ip_2d(ji, jl1) * afrdg
            aprdg2(ji) = a_ip_2d(ji, jl1) * afrdg * hi_hrdg(ji, jl1)
            vprdg(ji) = v_ip_2d(ji, jl1) * afrdg
            aprft1 = a_ip_2d(ji, jl1) * afrft
            aprft2(ji) = a_ip_2d(ji, jl1) * afrft * hi_hrft
            vprft(ji) = v_ip_2d(ji, jl1) * afrft
          END IF
          wfx_dyn_1d(ji) = wfx_dyn_1d(ji) - vsw * rhoi * r1_rdtice
          sfx_dyn_1d(ji) = sfx_dyn_1d(ji) - vsw * sss_1d(ji) * rhoi * r1_rdtice
          hfx_dyn_1d(ji) = hfx_dyn_1d(ji) + ersw(ji) * r1_rdtice
          wfx_snw_dyn_1d(ji) = wfx_snw_dyn_1d(ji) + (rhos * vsrdg(ji) * (1._wp - rn_fsnwrdg) + rhos * vsrft(ji) * (1._wp - rn_fsnwrft)) * r1_rdtice
          IF (nn_icesal /= 2) THEN
            sirdg2(ji) = sirdg2(ji) - vsw * (sss_1d(ji) - s_i_1d(ji))
            sfx_bri_1d(ji) = sfx_bri_1d(ji) + sss_1d(ji) * vsw * rhoi * r1_rdtice - s_i_1d(ji) * vsw * rhoi * r1_rdtice
          END IF
          a_i_2d(ji, jl1) = a_i_2d(ji, jl1) - airdg1 - airft1
          v_i_2d(ji, jl1) = v_i_2d(ji, jl1) - virdg1 - virft(ji)
          v_s_2d(ji, jl1) = v_s_2d(ji, jl1) - vsrdg(ji) - vsrft(ji)
          sv_i_2d(ji, jl1) = sv_i_2d(ji, jl1) - sirdg1 - sirft(ji)
          oa_i_2d(ji, jl1) = oa_i_2d(ji, jl1) - oirdg1 - oirft1
          IF (ln_pnd_H12) THEN
            a_ip_2d(ji, jl1) = a_ip_2d(ji, jl1) - aprdg1 - aprft1
            v_ip_2d(ji, jl1) = v_ip_2d(ji, jl1) - vprdg(ji) - vprft(ji)
          END IF
        END IF
      END DO
      DO jk = 1, nlay_s
        DO ji = 1, npti
          IF (apartf(ji, jl1) > 0._wp .AND. closing_gross(ji) > 0._wp) THEN
            afrdg = aridge(ji, jl1) * closing_gross(ji) * rdt_ice * z1_ai(ji)
            afrft = araft(ji, jl1) * closing_gross(ji) * rdt_ice * z1_ai(ji)
            esrdg(ji, jk) = ze_s_2d(ji, jk, jl1) * afrdg
            esrft(ji, jk) = ze_s_2d(ji, jk, jl1) * afrft
            hfx_dyn_1d(ji) = hfx_dyn_1d(ji) + (- esrdg(ji, jk) * (1._wp - rn_fsnwrdg) - esrft(ji, jk) * (1._wp - rn_fsnwrft)) * r1_rdtice
            ze_s_2d(ji, jk, jl1) = ze_s_2d(ji, jk, jl1) * (1._wp - afrdg - afrft)
          END IF
        END DO
      END DO
      DO jk = 1, nlay_i
        DO ji = 1, npti
          IF (apartf(ji, jl1) > 0._wp .AND. closing_gross(ji) > 0._wp) THEN
            afrdg = aridge(ji, jl1) * closing_gross(ji) * rdt_ice * z1_ai(ji)
            afrft = araft(ji, jl1) * closing_gross(ji) * rdt_ice * z1_ai(ji)
            eirdg(ji, jk) = ze_i_2d(ji, jk, jl1) * afrdg + ersw(ji) * r1_nlay_i
            eirft(ji, jk) = ze_i_2d(ji, jk, jl1) * afrft
            ze_i_2d(ji, jk, jl1) = ze_i_2d(ji, jk, jl1) * (1._wp - afrdg - afrft)
          END IF
        END DO
      END DO
      itest_rdg(1 : npti) = 0
      itest_rft(1 : npti) = 0
      DO jl2 = 1, jpl
        DO ji = 1, npti
          IF (apartf(ji, jl1) > 0._wp .AND. closing_gross(ji) > 0._wp) THEN
            IF (hrmin(ji, jl1) <= hi_max(jl2) .AND. hrmax(ji, jl1) > hi_max(jl2 - 1)) THEN
              hL = MAX(hrmin(ji, jl1), hi_max(jl2 - 1))
              hR = MIN(hrmax(ji, jl1), hi_max(jl2))
              farea = (hR - hL) / (hrmax(ji, jl1) - hrmin(ji, jl1))
              fvol(ji) = (hR * hR - hL * hL) / (hrmax(ji, jl1) * hrmax(ji, jl1) - hrmin(ji, jl1) * hrmin(ji, jl1))
              itest_rdg(ji) = 1
            ELSE
              farea = 0._wp
              fvol(ji) = 0._wp
            END IF
            IF (hraft(ji, jl1) <= hi_max(jl2) .AND. hraft(ji, jl1) > hi_max(jl2 - 1)) THEN
              zswitch(ji) = 1._wp
              itest_rft(ji) = 1
            ELSE
              zswitch(ji) = 0._wp
            END IF
            IF (itest_rdg(ji) == 0 .AND. jl2 == jpl) THEN
              farea = 1._wp
              fvol(ji) = 1._wp
            END IF
            IF (itest_rft(ji) == 0 .AND. jl2 == jpl) zswitch(ji) = 1._wp
            a_i_2d(ji, jl2) = a_i_2d(ji, jl2) + (airdg2(ji) * farea + airft2(ji) * zswitch(ji))
            oa_i_2d(ji, jl2) = oa_i_2d(ji, jl2) + (oirdg2(ji) * farea + oirft2(ji) * zswitch(ji))
            v_i_2d(ji, jl2) = v_i_2d(ji, jl2) + (virdg2(ji) * fvol(ji) + virft(ji) * zswitch(ji))
            sv_i_2d(ji, jl2) = sv_i_2d(ji, jl2) + (sirdg2(ji) * fvol(ji) + sirft(ji) * zswitch(ji))
            v_s_2d(ji, jl2) = v_s_2d(ji, jl2) + (vsrdg(ji) * rn_fsnwrdg * fvol(ji) + vsrft(ji) * rn_fsnwrft * zswitch(ji))
            IF (ln_pnd_H12) THEN
              v_ip_2d(ji, jl2) = v_ip_2d(ji, jl2) + (vprdg(ji) * rn_fpndrdg * fvol(ji) + vprft(ji) * rn_fpndrft * zswitch(ji))
              a_ip_2d(ji, jl2) = a_ip_2d(ji, jl2) + (aprdg2(ji) * rn_fpndrdg * farea + aprft2(ji) * rn_fpndrft * zswitch(ji))
            END IF
          END IF
        END DO
        DO jk = 1, nlay_s
          DO ji = 1, npti
            IF (apartf(ji, jl1) > 0._wp .AND. closing_gross(ji) > 0._wp) ze_s_2d(ji, jk, jl2) = ze_s_2d(ji, jk, jl2) + (esrdg(ji, jk) * rn_fsnwrdg * fvol(ji) + esrft(ji, jk) * rn_fsnwrft * zswitch(ji))
          END DO
        END DO
        DO jk = 1, nlay_i
          DO ji = 1, npti
            IF (apartf(ji, jl1) > 0._wp .AND. closing_gross(ji) > 0._wp) ze_i_2d(ji, jk, jl2) = ze_i_2d(ji, jk, jl2) + eirdg(ji, jk) * fvol(ji) + eirft(ji, jk) * zswitch(ji)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    CALL profile_psy_data2 % PreStart('rdgrft_shift', 'r2', 0, 0)
    CALL ice_var_roundoff(a_i_2d, v_i_2d, v_s_2d, sv_i_2d, oa_i_2d, a_ip_2d, v_ip_2d, ze_s_2d, ze_i_2d)
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE rdgrft_shift
  SUBROUTINE ice_strength
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj, jl
    INTEGER :: ismooth
    INTEGER :: itframe
    REAL(KIND = wp) :: zp, z1_3
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zworka
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zstrp1, zstrp2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_str_H79) THEN
      CALL profile_psy_data0 % PreStart('ice_strength', 'r0', 0, 0)
      strength(:, :) = rn_pstar * SUM(v_i(:, :, :), dim = 3) * EXP(- rn_crhg * (1._wp - SUM(a_i(:, :, :), dim = 3)))
      ismooth = 1
      CALL profile_psy_data0 % PostEnd
    ELSE
      !$ACC KERNELS
      strength(:, :) = 0._wp
      ismooth = 0
      !$ACC END KERNELS
    END IF
    SELECT CASE (ismooth)
    CASE (1)
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          IF (SUM(a_i(ji, jj, :)) > 0._wp) THEN
            zworka(ji, jj) = (4.0 * strength(ji, jj) + strength(ji - 1, jj) * tmask(ji - 1, jj, 1) + strength(ji + 1, jj) * tmask(ji + 1, jj, 1) + strength(ji, jj - 1) * tmask(ji, jj - 1, 1) + strength(ji, jj + 1) * tmask(ji, jj + 1, 1)) / (4.0 + tmask(ji - 1, jj, 1) + tmask(ji + 1, jj, 1) + tmask(ji, jj - 1, 1) + tmask(ji, jj + 1, 1))
          ELSE
            zworka(ji, jj) = 0._wp
          END IF
        END DO
      END DO
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          strength(ji, jj) = zworka(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('ice_strength', 'r1', 0, 0)
      CALL lbc_lnk('icedyn_rdgrft', strength, 'T', 1.)
      CALL profile_psy_data1 % PostEnd
    CASE (2)
      !$ACC KERNELS
      IF (kt_ice == nit000) THEN
        zstrp1(:, :) = 0._wp
        zstrp2(:, :) = 0._wp
      END IF
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          IF (SUM(a_i(ji, jj, :)) > 0._wp) THEN
            itframe = 1
            IF (zstrp1(ji, jj) > 0._wp) itframe = itframe + 1
            IF (zstrp2(ji, jj) > 0._wp) itframe = itframe + 1
            zp = (strength(ji, jj) + zstrp1(ji, jj) + zstrp2(ji, jj)) / itframe
            zstrp2(ji, jj) = zstrp1(ji, jj)
            zstrp1(ji, jj) = strength(ji, jj)
            strength(ji, jj) = zp
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('ice_strength', 'r2', 0, 0)
      CALL lbc_lnk('icedyn_rdgrft', strength, 'T', 1.)
      CALL profile_psy_data2 % PostEnd
    END SELECT
  END SUBROUTINE ice_strength
  SUBROUTINE ice_dyn_1d2d(kn)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kn
    INTEGER :: jl, jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_dyn_1d2d', 'r0', 0, 0)
    SELECT CASE (kn)
    CASE (1)
      CALL tab_2d_1d(npti, nptidx(1 : npti), sss_1d(1 : npti), sss_m(:, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), sst_1d(1 : npti), sst_m(:, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), v_s_2d(1 : npti, 1 : jpl), v_s(:, :, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), sv_i_2d(1 : npti, 1 : jpl), sv_i(:, :, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), oa_i_2d(1 : npti, 1 : jpl), oa_i(:, :, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), a_ip_2d(1 : npti, 1 : jpl), a_ip(:, :, :))
      CALL tab_3d_2d(npti, nptidx(1 : npti), v_ip_2d(1 : npti, 1 : jpl), v_ip(:, :, :))
      DO jl = 1, jpl
        DO jk = 1, nlay_s
          CALL tab_2d_1d(npti, nptidx(1 : npti), ze_s_2d(1 : npti, jk, jl), e_s(:, :, jk, jl))
        END DO
        DO jk = 1, nlay_i
          CALL tab_2d_1d(npti, nptidx(1 : npti), ze_i_2d(1 : npti, jk, jl), e_i(:, :, jk, jl))
        END DO
      END DO
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_dyn_1d(1 : npti), sfx_dyn(:, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), sfx_bri_1d(1 : npti), sfx_bri(:, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_dyn_1d(1 : npti), wfx_dyn(:, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), hfx_dyn_1d(1 : npti), hfx_dyn(:, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_snw_dyn_1d(1 : npti), wfx_snw_dyn(:, :))
      CALL tab_2d_1d(npti, nptidx(1 : npti), wfx_pnd_1d(1 : npti), wfx_pnd(:, :))
    CASE (2)
      CALL tab_1d_2d(npti, nptidx(1 : npti), ato_i_1d(1 : npti), ato_i(:, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), a_i_2d(1 : npti, 1 : jpl), a_i(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), v_i_2d(1 : npti, 1 : jpl), v_i(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), v_s_2d(1 : npti, 1 : jpl), v_s(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), sv_i_2d(1 : npti, 1 : jpl), sv_i(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), oa_i_2d(1 : npti, 1 : jpl), oa_i(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), a_ip_2d(1 : npti, 1 : jpl), a_ip(:, :, :))
      CALL tab_2d_3d(npti, nptidx(1 : npti), v_ip_2d(1 : npti, 1 : jpl), v_ip(:, :, :))
      DO jl = 1, jpl
        DO jk = 1, nlay_s
          CALL tab_1d_2d(npti, nptidx(1 : npti), ze_s_2d(1 : npti, jk, jl), e_s(:, :, jk, jl))
        END DO
        DO jk = 1, nlay_i
          CALL tab_1d_2d(npti, nptidx(1 : npti), ze_i_2d(1 : npti, jk, jl), e_i(:, :, jk, jl))
        END DO
      END DO
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_dyn_1d(1 : npti), sfx_dyn(:, :))
      CALL tab_1d_2d(npti, nptidx(1 : npti), sfx_bri_1d(1 : npti), sfx_bri(:, :))
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_dyn_1d(1 : npti), wfx_dyn(:, :))
      CALL tab_1d_2d(npti, nptidx(1 : npti), hfx_dyn_1d(1 : npti), hfx_dyn(:, :))
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_snw_dyn_1d(1 : npti), wfx_snw_dyn(:, :))
      CALL tab_1d_2d(npti, nptidx(1 : npti), wfx_pnd_1d(1 : npti), wfx_pnd(:, :))
    END SELECT
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_dyn_1d2d
  SUBROUTINE ice_dyn_rdgrft_init
    INTEGER :: ios
    NAMELIST /namdyn_rdgrft/ ln_str_H79, rn_pstar, rn_crhg, rn_csrdg, ln_partf_lin, rn_gstar, ln_partf_exp, rn_astar, ln_ridging, rn_hstar, rn_porordg, rn_fsnwrdg, rn_fpndrdg, ln_rafting, rn_hraft, rn_craft, rn_fsnwrft, rn_fpndrft
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namdyn_rdgrft, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_rdgrft in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namdyn_rdgrft, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_rdgrft in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namdyn_rdgrft)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_dyn_rdgrft_init: ice parameters for ridging/rafting '
      WRITE(numout, *) '~~~~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namdyn_rdgrft:'
      WRITE(numout, *) '      ice strength parameterization Hibler (1979)              ln_str_H79   = ', ln_str_H79
      WRITE(numout, *) '            1st bulk-rheology parameter                        rn_pstar     = ', rn_pstar
      WRITE(numout, *) '            2nd bulk-rhelogy parameter                         rn_crhg      = ', rn_crhg
      WRITE(numout, *) '      Fraction of shear energy contributing to ridging         rn_csrdg     = ', rn_csrdg
      WRITE(numout, *) '      linear ridging participation function                    ln_partf_lin = ', ln_partf_lin
      WRITE(numout, *) '            Fraction of ice coverage contributing to ridging   rn_gstar     = ', rn_gstar
      WRITE(numout, *) '      Exponential ridging participation function               ln_partf_exp = ', ln_partf_exp
      WRITE(numout, *) '            Equivalent to G* for an exponential function       rn_astar     = ', rn_astar
      WRITE(numout, *) '      Ridging of ice sheets or not                             ln_ridging   = ', ln_ridging
      WRITE(numout, *) '            max ridged ice thickness                           rn_hstar     = ', rn_hstar
      WRITE(numout, *) '            Initial porosity of ridges                         rn_porordg   = ', rn_porordg
      WRITE(numout, *) '            Fraction of snow volume conserved during ridging   rn_fsnwrdg   = ', rn_fsnwrdg
      WRITE(numout, *) '            Fraction of pond volume conserved during ridging   rn_fpndrdg   = ', rn_fpndrdg
      WRITE(numout, *) '      Rafting of ice sheets or not                             ln_rafting   = ', ln_rafting
      WRITE(numout, *) '            Parmeter thickness (threshold between ridge-raft)  rn_hraft     = ', rn_hraft
      WRITE(numout, *) '            Rafting hyperbolic tangent coefficient             rn_craft     = ', rn_craft
      WRITE(numout, *) '            Fraction of snow volume conserved during rafting   rn_fsnwrft   = ', rn_fsnwrft
      WRITE(numout, *) '            Fraction of pond volume conserved during rafting   rn_fpndrft   = ', rn_fpndrft
    END IF
    IF ((ln_partf_lin .AND. ln_partf_exp) .OR. (.NOT. ln_partf_lin .AND. .NOT. ln_partf_exp)) THEN
      CALL ctl_stop('ice_dyn_rdgrft_init: choose one and only one participation function (ln_partf_lin or ln_partf_exp)')
    END IF
    IF (.NOT. ln_icethd) THEN
      rn_porordg = 0._wp
      rn_fsnwrdg = 1._wp
      rn_fsnwrft = 1._wp
      rn_fpndrdg = 1._wp
      rn_fpndrft = 1._wp
      IF (lwp) THEN
        WRITE(numout, *) '      ==> only ice dynamics is activated, thus some parameters must be changed'
        WRITE(numout, *) '            rn_porordg   = ', rn_porordg
        WRITE(numout, *) '            rn_fsnwrdg   = ', rn_fsnwrdg
        WRITE(numout, *) '            rn_fpndrdg   = ', rn_fpndrdg
        WRITE(numout, *) '            rn_fsnwrft   = ', rn_fsnwrft
        WRITE(numout, *) '            rn_fpndrft   = ', rn_fpndrft
      END IF
    END IF
    IF (ice_dyn_rdgrft_alloc() /= 0) CALL ctl_stop('STOP', 'ice_dyn_rdgrft_init: unable to allocate arrays')
  END SUBROUTINE ice_dyn_rdgrft_init
END MODULE icedyn_rdgrft
