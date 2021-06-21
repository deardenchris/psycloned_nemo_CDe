MODULE icedyn_rhg
  USE phycst
  USE dom_oce
  USE ice
  USE icedyn_rhg_evp
  USE icectl
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_rhg
  PUBLIC :: ice_dyn_rhg_init
  INTEGER :: nice_rhg
  INTEGER, PARAMETER :: np_rhgEVP = 1
  LOGICAL :: ln_rhg_EVP
  CONTAINS
  SUBROUTINE ice_dyn_rhg(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ice_dyn_rhg', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('icedyn_rhg')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_dyn_rhg: sea-ice rheology'
      WRITE(numout, *) '~~~~~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_landfast_home) THEN
      !$ACC KERNELS
      tau_icebfr(:, :) = 0._wp
      DO jl = 1, jpl
        WHERE (h_i(:, :, jl) > ht_n(:, :) * rn_depfra) tau_icebfr(:, :) = tau_icebfr(:, :) + a_i(:, :, jl) * rn_icebfr
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data1 % PreStart('ice_dyn_rhg', 'r1', 0, 0)
    SELECT CASE (nice_rhg)
    CASE (np_rhgEVP)
      CALL ice_dyn_rhg_evp(kt, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i)
    END SELECT
    IF (lrst_ice) THEN
      IF (ln_rhg_EVP) CALL rhg_evp_rst('WRITE', kt)
    END IF
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (ln_ctl) CALL ice_prt3D('icedyn_rhg')
    IF (ln_timing) CALL timing_stop('icedyn_rhg')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ice_dyn_rhg
  SUBROUTINE ice_dyn_rhg_init
    INTEGER :: ios, ioptio
    NAMELIST /namdyn_rhg/ ln_rhg_EVP, ln_aEVP, rn_creepl, rn_ecc, nn_nevp, rn_relast
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namdyn_rhg, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_rhg in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namdyn_rhg, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_rhg in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namdyn_rhg)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_dyn_rhg_init: ice parameters for ice dynamics '
      WRITE(numout, *) '~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist : namdyn_rhg:'
      WRITE(numout, *) '      rheology EVP (icedyn_rhg_evp)                        ln_rhg_EVP = ', ln_rhg_EVP
      WRITE(numout, *) '         use adaptive EVP (aEVP)                           ln_aEVP    = ', ln_aEVP
      WRITE(numout, *) '         creep limit                                       rn_creepl  = ', rn_creepl
      WRITE(numout, *) '         eccentricity of the elliptical yield curve        rn_ecc     = ', rn_ecc
      WRITE(numout, *) '         number of iterations for subcycling               nn_nevp    = ', nn_nevp
      WRITE(numout, *) '         ratio of elastic timescale over ice time step     rn_relast  = ', rn_relast
    END IF
    ioptio = 0
    IF (ln_rhg_evp) THEN
      ioptio = ioptio + 1
      nice_rhg = np_rhgevp
    END IF
    IF (ioptio /= 1) CALL ctl_stop('ice_dyn_rhg_init: choose one and only one ice rheology')
    IF (ln_rhg_EVP) CALL rhg_evp_rst('READ')
  END SUBROUTINE ice_dyn_rhg_init
END MODULE icedyn_rhg