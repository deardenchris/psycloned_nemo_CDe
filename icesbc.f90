MODULE icesbc
  USE oce
  USE dom_oce
  USE ice
  USE sbc_oce
  USE sbc_ice
  USE usrdef_sbc
  USE sbcblk
  USE sbccpl
  USE icealb
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_sbc_tau
  PUBLIC :: ice_sbc_flx
  PUBLIC :: ice_sbc_init
  CONTAINS
  SUBROUTINE ice_sbc_tau(kt, ksbc, utau_ice, vtau_ice)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: ksbc
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(OUT) :: utau_ice, vtau_ice
    INTEGER :: ji, jj
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zutau_ice, zvtau_ice
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('ice_sbc_tau', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ice_sbc')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_sbc_tau: Surface boundary condition for sea ice (momentum)'
      WRITE(numout, *) '~~~~~~~~~~~~~~~'
    END IF
    SELECT CASE (ksbc)
    CASE (jp_usr)
      CALL usrdef_sbc_ice_tau(kt)
    CASE (jp_blk)
      CALL blk_ice_tau
    CASE (jp_purecpl)
      CALL sbc_cpl_ice_tau(utau_ice, vtau_ice)
    END SELECT
    CALL profile_psy_data0 % PostEnd
    IF (ln_mixcpl) THEN
      CALL profile_psy_data1 % PreStart('ice_sbc_tau', 'r1', 0, 0)
      CALL sbc_cpl_ice_tau(zutau_ice, zvtau_ice)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          utau_ice(ji, jj) = utau_ice(ji, jj) * xcplmask(ji, jj, 0) + zutau_ice(ji, jj) * (1. - xcplmask(ji, jj, 0))
          vtau_ice(ji, jj) = vtau_ice(ji, jj) * xcplmask(ji, jj, 0) + zvtau_ice(ji, jj) * (1. - xcplmask(ji, jj, 0))
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('ice_sbc_tau', 'r2', 0, 0)
      CALL lbc_lnk_multi('icesbc', utau_ice, 'U', - 1., vtau_ice, 'V', - 1.)
      CALL profile_psy_data2 % PostEnd
    END IF
    CALL profile_psy_data3 % PreStart('ice_sbc_tau', 'r3', 0, 0)
    IF (ln_timing) CALL timing_stop('ice_sbc')
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE ice_sbc_tau
  SUBROUTINE ice_sbc_flx(kt, ksbc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: ksbc
    INTEGER :: ji, jj, jl
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zalb_os, zalb_cs
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zalb
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ice_sbc_flx', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ice_sbc_flx')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_sbc_flx: Surface boundary condition for sea ice (flux)'
      WRITE(numout, *) '~~~~~~~~~~~~~~~'
    END IF
    CALL ice_alb(t_su, h_i, h_s, ln_pnd_alb, a_ip_frac, h_ip, zalb_cs, zalb_os)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    alb_ice(:, :, :) = (1. - cldf_ice) * zalb_cs(:, :, :) + cldf_ice * zalb_os(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_sbc_flx', 'r1', 0, 0)
    SELECT CASE (ksbc)
    CASE (jp_usr)
      CALL usrdef_sbc_ice_flx(kt, h_s, h_i)
    CASE (jp_blk)
      CALL blk_ice_flx(t_su, h_s, h_i, alb_ice)
      IF (ln_mixcpl) CALL sbc_cpl_ice_flx(picefr = at_i_b, palbi = alb_ice, psst = sst_m, pist = t_su, phs = h_s, phi = h_i)
      IF (nn_flxdist /= - 1) CALL ice_flx_dist(t_su, alb_ice, qns_ice, qsr_ice, dqns_ice, evap_ice, devap_ice, nn_flxdist)
      IF (ln_cndflx .AND. .NOT. ln_cndemulate) CALL blk_ice_qcn(ln_virtual_itd, t_su, t_bo, h_s, h_i)
    CASE (jp_purecpl)
      CALL sbc_cpl_ice_flx(picefr = at_i_b, palbi = alb_ice, psst = sst_m, pist = t_su, phs = h_s, phi = h_i)
      IF (nn_flxdist /= - 1) CALL ice_flx_dist(t_su, alb_ice, qns_ice, qsr_ice, dqns_ice, evap_ice, devap_ice, nn_flxdist)
    END SELECT
    IF (iom_use('icealb')) THEN
      WHERE (at_i_b <= epsi06)
        zalb(:, :) = rn_alb_oce
      ELSEWHERE
        zalb(:, :) = SUM(alb_ice * a_i_b, dim = 3) / at_i_b
      END WHERE
      CALL iom_put("icealb", zalb(:, :))
    END IF
    IF (iom_use('albedo')) THEN
      zalb(:, :) = SUM(alb_ice * a_i_b, dim = 3) + rn_alb_oce * (1._wp - at_i_b)
      CALL iom_put("albedo", zalb(:, :))
    END IF
    IF (ln_timing) CALL timing_stop('ice_sbc_flx')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ice_sbc_flx
  SUBROUTINE ice_flx_dist(ptn_ice, palb_ice, pqns_ice, pqsr_ice, pdqn_ice, pevap_ice, pdevap_ice, k_flxdist)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: k_flxdist
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptn_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: palb_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pqns_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pqsr_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pdqn_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pevap_ice
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pdevap_ice
    INTEGER :: jl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z1_at_i
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z_qsr_m
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z_qns_m
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z_evap_m
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z_dqn_m
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z_devap_m
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zalb_m
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: ztem_m
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    WHERE (at_i(:, :) > 0._wp)
      z1_at_i(:, :) = 1._wp / at_i(:, :)
    ELSEWHERE
      z1_at_i(:, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    SELECT CASE (k_flxdist)
    CASE (0, 1)
      CALL profile_psy_data0 % PreStart('ice_flx_dist', 'r0', 0, 0)
      ALLOCATE(z_qns_m(jpi, jpj), z_qsr_m(jpi, jpj), z_dqn_m(jpi, jpj), z_evap_m(jpi, jpj), z_devap_m(jpi, jpj))
      z_qns_m(:, :) = SUM(a_i(:, :, :) * pqns_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      z_qsr_m(:, :) = SUM(a_i(:, :, :) * pqsr_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      z_dqn_m(:, :) = SUM(a_i(:, :, :) * pdqn_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      z_evap_m(:, :) = SUM(a_i(:, :, :) * pevap_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      z_devap_m(:, :) = SUM(a_i(:, :, :) * pdevap_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      CALL profile_psy_data0 % PostEnd
      !$ACC KERNELS
      DO jl = 1, jpl
        pqns_ice(:, :, jl) = z_qns_m(:, :)
        pqsr_ice(:, :, jl) = z_qsr_m(:, :)
        pdqn_ice(:, :, jl) = z_dqn_m(:, :)
        pevap_ice(:, :, jl) = z_evap_m(:, :)
        pdevap_ice(:, :, jl) = z_devap_m(:, :)
      END DO
      !$ACC END KERNELS
      DEALLOCATE(z_qns_m, z_qsr_m, z_dqn_m, z_evap_m, z_devap_m)
    END SELECT
    SELECT CASE (k_flxdist)
    CASE (1, 2)
      CALL profile_psy_data1 % PreStart('ice_flx_dist', 'r1', 0, 0)
      ALLOCATE(zalb_m(jpi, jpj), ztem_m(jpi, jpj))
      zalb_m(:, :) = SUM(a_i(:, :, :) * palb_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      ztem_m(:, :) = SUM(a_i(:, :, :) * ptn_ice(:, :, :), dim = 3) * z1_at_i(:, :)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      DO jl = 1, jpl
        pqns_ice(:, :, jl) = pqns_ice(:, :, jl) + pdqn_ice(:, :, jl) * (ptn_ice(:, :, jl) - ztem_m(:, :))
        pevap_ice(:, :, jl) = pevap_ice(:, :, jl) + pdevap_ice(:, :, jl) * (ptn_ice(:, :, jl) - ztem_m(:, :))
        pqsr_ice(:, :, jl) = pqsr_ice(:, :, jl) * (1._wp - palb_ice(:, :, jl)) / (1._wp - zalb_m(:, :))
      END DO
      !$ACC END KERNELS
      DEALLOCATE(zalb_m, ztem_m)
    END SELECT
  END SUBROUTINE ice_flx_dist
  SUBROUTINE ice_sbc_init
    INTEGER :: ios, ioptio
    NAMELIST /namsbc/ rn_cio, rn_blow_s, nn_flxdist, ln_cndflx, ln_cndemulate
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namsbc, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namsbc, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namsbc)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_sbc_init: ice parameters for ice dynamics '
      WRITE(numout, *) '~~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namsbc:'
      WRITE(numout, *) '      drag coefficient for oceanic stress              rn_cio        = ', rn_cio
      WRITE(numout, *) '      coefficient for ice-lead partition of snowfall   rn_blow_s     = ', rn_blow_s
      WRITE(numout, *) '      Multicategory heat flux formulation              nn_flxdist    = ', nn_flxdist
      WRITE(numout, *) '      Use conduction flux as surface condition         ln_cndflx     = ', ln_cndflx
      WRITE(numout, *) '         emulate conduction flux                       ln_cndemulate = ', ln_cndemulate
    END IF
    IF (lwp) WRITE(numout, *)
    SELECT CASE (nn_flxdist)
    CASE (- 1)
      IF (lwp) WRITE(numout, *) '   SI3: use per-category fluxes (nn_flxdist = -1) '
    CASE (0)
      IF (lwp) WRITE(numout, *) '   SI3: use average per-category fluxes (nn_flxdist = 0) '
    CASE (1)
      IF (lwp) WRITE(numout, *) '   SI3: use average then redistribute per-category fluxes (nn_flxdist = 1) '
      IF (ln_cpl) CALL ctl_stop('ice_thd_init: the chosen nn_flxdist for SI3 in coupled mode must be /=1')
    CASE (2)
      IF (lwp) WRITE(numout, *) '   SI3: Redistribute a single flux over categories (nn_flxdist = 2) '
      IF (.NOT. ln_cpl) CALL ctl_stop('ice_thd_init: the chosen nn_flxdist for SI3 in forced mode must be /=2')
    CASE DEFAULT
      CALL ctl_stop('ice_thd_init: SI3 option, nn_flxdist, should be between -1 and 2')
    END SELECT
  END SUBROUTINE ice_sbc_init
END MODULE icesbc