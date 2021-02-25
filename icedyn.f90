MODULE icedyn
  USE phycst
  USE dom_oce
  USE ice
  USE icedyn_rhg
  USE icedyn_adv
  USE icedyn_rdgrft
  USE icecor
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
  PUBLIC :: ice_dyn
  PUBLIC :: ice_dyn_init
  INTEGER :: nice_dyn
  INTEGER, PARAMETER :: np_dynALL = 1
  INTEGER, PARAMETER :: np_dynRHGADV = 2
  INTEGER, PARAMETER :: np_dynADV1D = 3
  INTEGER, PARAMETER :: np_dynADV2D = 4
  LOGICAL :: ln_dynALL
  LOGICAL :: ln_dynRHGADV
  LOGICAL :: ln_dynADV1D
  LOGICAL :: ln_dynADV2D
  REAL(KIND = wp) :: rn_uice
  REAL(KIND = wp) :: rn_vice
  CONTAINS
  SUBROUTINE ice_dyn(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zcoefu, zcoefv
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdivu_i
    IF (ln_timing) CALL timing_start('icedyn')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn: sea-ice dynamics'
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    WHERE (a_i(:, :, :) >= epsi20)
      h_i(:, :, :) = v_i(:, :, :) / a_i_b(:, :, :)
      h_s(:, :, :) = v_s(:, :, :) / a_i_b(:, :, :)
    ELSEWHERE
      h_i(:, :, :) = 0._wp
      h_s(:, :, :) = 0._wp
    END WHERE
    WHERE (a_ip(:, :, :) >= epsi20)
      h_ip(:, :, :) = v_ip(:, :, :) / a_ip(:, :, :)
    ELSEWHERE
      h_ip(:, :, :) = 0._wp
    END WHERE
    SELECT CASE (nice_dyn)
    CASE (np_dynALL)
      CALL ice_dyn_rhg(kt)
      CALL ice_dyn_adv(kt)
      CALL ice_dyn_rdgrft(kt)
      CALL ice_cor(kt, 1)
    CASE (np_dynRHGADV)
      CALL ice_dyn_rhg(kt)
      CALL ice_dyn_adv(kt)
      CALL Hpiling
      CALL ice_var_zapsmall
    CASE (np_dynADV1D)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zcoefu = (REAL(jpiglo + 1) * 0.5 - REAL(ji + nimpp - 1)) / (REAL(jpiglo + 1) * 0.5 - 1.)
          zcoefv = (REAL(jpjglo + 1) * 0.5 - REAL(jj + njmpp - 1)) / (REAL(jpjglo + 1) * 0.5 - 1.)
          u_ice(ji, jj) = rn_uice * 1.5 * SIGN(1., zcoefu) * ABS(zcoefu) * umask(ji, jj, 1)
          v_ice(ji, jj) = rn_vice * 1.5 * SIGN(1., zcoefv) * ABS(zcoefv) * vmask(ji, jj, 1)
        END DO
      END DO
      CALL ice_dyn_adv(kt)
    CASE (np_dynADV2D)
      u_ice(:, :) = rn_uice * umask(:, :, 1)
      v_ice(:, :) = rn_vice * vmask(:, :, 1)
      CALL ice_dyn_adv(kt)
    END SELECT
    IF (iom_use('icediv')) THEN
      SELECT CASE (nice_dyn)
      CASE (np_dynADV1D, np_dynADV2D)
        ALLOCATE(zdivu_i(jpi, jpj))
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdivu_i(ji, jj) = (e2u(ji, jj) * u_ice(ji, jj) - e2u(ji - 1, jj) * u_ice(ji - 1, jj) + e1v(ji, jj) * v_ice(ji, jj) - &
&e1v(ji, jj - 1) * v_ice(ji, jj - 1)) * r1_e1e2t(ji, jj)
          END DO
        END DO
        CALL lbc_lnk('icedyn', zdivu_i, 'T', 1.)
        CALL iom_put("icediv", zdivu_i(:, :))
        DEALLOCATE(zdivu_i)
      END SELECT
    END IF
    IF (ln_timing) CALL timing_stop('icedyn')
  END SUBROUTINE ice_dyn
  SUBROUTINE Hpiling
    INTEGER :: jl
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'Hpiling', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    at_i(:, :) = SUM(a_i(:, :, :), dim = 3)
    DO jl = 1, jpl
      WHERE (at_i(:, :) > epsi20)
        a_i(:, :, jl) = a_i(:, :, jl) * (1._wp + MIN(rn_amax_2d(:, :) - at_i(:, :), 0._wp) / at_i(:, :))
      END WHERE
    END DO
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'Hpiling', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
  END SUBROUTINE Hpiling
  SUBROUTINE ice_dyn_init
    INTEGER :: ios, ioptio
    NAMELIST /namdyn/ ln_dynALL, ln_dynRHGADV, ln_dynADV1D, ln_dynADV2D, rn_uice, rn_vice, rn_ishlat, ln_landfast_L16, &
&ln_landfast_home, rn_depfra, rn_icebfr, rn_lfrelax, rn_tensile
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namdyn, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namdyn, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namdyn)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn_init: ice parameters for ice dynamics '
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namdyn:'
      WRITE(numout, FMT = *) '      Full ice dynamics      (rhg + adv + ridge/raft + corr) ln_dynALL       = ', ln_dynALL
      WRITE(numout, FMT = *) '      No ridge/raft & No cor (rhg + adv)                     ln_dynRHGADV    = ', ln_dynRHGADV
      WRITE(numout, FMT = *) '      Advection 1D only      (Schar & Smolarkiewicz 1996)    ln_dynADV1D     = ', ln_dynADV1D
      WRITE(numout, FMT = *) '      Advection 2D only      (rn_uvice + adv)                ln_dynADV2D     = ', ln_dynADV2D
      WRITE(numout, FMT = *) '         with prescribed velocity given by   (u,v)_ice = (rn_uice,rn_vice)   = (', rn_uice, ',', &
&rn_vice, ')'
      WRITE(numout, FMT = *) '      lateral boundary condition for sea ice dynamics        rn_ishlat       = ', rn_ishlat
      WRITE(numout, FMT = *) '      Landfast: param from Lemieux 2016                      ln_landfast_L16 = ', ln_landfast_L16
      WRITE(numout, FMT = *) '      Landfast: param from home made                         ln_landfast_home= ', ln_landfast_home
      WRITE(numout, FMT = *) '         fraction of ocean depth that ice must reach         rn_depfra       = ', rn_depfra
      WRITE(numout, FMT = *) '         maximum bottom stress per unit area of contact      rn_icebfr       = ', rn_icebfr
      WRITE(numout, FMT = *) '         relax time scale (s-1) to reach static friction     rn_lfrelax      = ', rn_lfrelax
      WRITE(numout, FMT = *) '         isotropic tensile strength                          rn_tensile      = ', rn_tensile
      WRITE(numout, FMT = *)
    END IF
    ioptio = 0
    IF (ln_dynall) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynall
    END IF
    IF (ln_dynrhgadv) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynrhgadv
    END IF
    IF (ln_dynadv1d) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynadv1d
    END IF
    IF (ln_dynadv2d) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynadv2d
    END IF
    IF (ioptio /= 1) CALL ctl_stop('ice_dyn_init: one and only one ice dynamics option has to be defined ')
    IF (rn_ishlat == 0.) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  free-slip'
    ELSE IF (rn_ishlat == 2.) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  no-slip'
    ELSE IF (0. < rn_ishlat .AND. rn_ishlat < 2.) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  partial-slip'
    ELSE IF (2. < rn_ishlat) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  strong-slip'
    END IF
    IF (.NOT. ln_landfast_L16 .AND. .NOT. ln_landfast_home) tau_icebfr(:, :) = 0._wp
    IF (ln_landfast_L16 .AND. ln_landfast_home) THEN
      CALL ctl_stop('ice_dyn_init: choose one and only one landfast parameterization (ln_landfast_L16 or ln_landfast_home)')
    END IF
    CALL ice_dyn_rdgrft_init
    CALL ice_dyn_rhg_init
    CALL ice_dyn_adv_init
  END SUBROUTINE ice_dyn_init
END MODULE icedyn