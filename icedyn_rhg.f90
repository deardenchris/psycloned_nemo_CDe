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
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jl
    IF (ln_timing) CALL timing_start('icedyn_rhg')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn_rhg: sea-ice rheology'
      WRITE(numout, FMT = *) '~~~~~~~~~~~'
    END IF
    IF (ln_landfast_home) THEN
      tau_icebfr(:, :) = 0._wp
      DO jl = 1, jpl
        WHERE (h_i(:, :, jl) > ht_n(:, :) * rn_depfra) tau_icebfr(:, :) = tau_icebfr(:, :) + a_i(:, :, jl) * rn_icebfr
      END DO
    END IF
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
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn_rhg_init: ice parameters for ice dynamics '
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist : namdyn_rhg:'
      WRITE(numout, FMT = *) '      rheology EVP (icedyn_rhg_evp)                        ln_rhg_EVP = ', ln_rhg_EVP
      WRITE(numout, FMT = *) '         use adaptive EVP (aEVP)                           ln_aEVP    = ', ln_aEVP
      WRITE(numout, FMT = *) '         creep limit                                       rn_creepl  = ', rn_creepl
      WRITE(numout, FMT = *) '         eccentricity of the elliptical yield curve        rn_ecc     = ', rn_ecc
      WRITE(numout, FMT = *) '         number of iterations for subcycling               nn_nevp    = ', nn_nevp
      WRITE(numout, FMT = *) '         ratio of elastic timescale over ice time step     rn_relast  = ', rn_relast
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