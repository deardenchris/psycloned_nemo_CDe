MODULE domhgr
  USE dom_oce
  USE par_oce
  USE phycst
  USE usrdef_hgr
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_hgr
  CONTAINS
  SUBROUTINE dom_hgr
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj
    INTEGER :: ie1e2u_v
    INTEGER :: iff
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('dom_hgr', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dom_hgr')
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dom_hgr : define the horizontal mesh from ithe following par_oce parameters '
      WRITE(numout, *) '~~~~~~~   '
      WRITE(numout, *) '   namcfg : read (=T) or user defined (=F) configuration    ln_read_cfg  = ', ln_read_cfg
    END IF
    IF (ln_read_cfg) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   read horizontal mesh in ', TRIM(cn_domcfg), ' file'
      CALL hgr_read(glamt, glamu, glamv, glamf, gphit, gphiu, gphiv, gphif, iff, ff_f, ff_t, e1t, e1u, e1v, e1f, e2t, e2u, e2v, e2f, ie1e2u_v, e1e2u, e1e2v)
    ELSE
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '          User defined horizontal mesh (usr_def_hgr)'
      CALL usr_def_hgr(glamt, glamu, glamv, glamf, gphit, gphiu, gphiv, gphif, iff, ff_f, ff_t, e1t, e1u, e1v, e1f, e2t, e2u, e2v, e2f, ie1e2u_v, e1e2u, e1e2v)
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (iff == 0) THEN
      IF (lwp) WRITE(numout, *) '          Coriolis parameter calculated on the sphere from gphif & gphit'
      !$ACC KERNELS
      ff_f(:, :) = 2. * omega * SIN(rad * gphif(:, :))
      ff_t(:, :) = 2. * omega * SIN(rad * gphit(:, :))
      !$ACC END KERNELS
    ELSE
      CALL profile_psy_data1 % PreStart('dom_hgr', 'r1', 0, 0)
      IF (ln_read_cfg) THEN
        IF (lwp) WRITE(numout, *) '          Coriolis parameter have been read in ', TRIM(cn_domcfg), ' file'
      ELSE
        IF (lwp) WRITE(numout, *) '          Coriolis parameter have been set in usr_def_hgr routine'
      END IF
      CALL profile_psy_data1 % PostEnd
    END IF
    !$ACC KERNELS
    r1_e1t(:, :) = 1._wp / e1t(:, :)
    r1_e2t(:, :) = 1._wp / e2t(:, :)
    r1_e1u(:, :) = 1._wp / e1u(:, :)
    r1_e2u(:, :) = 1._wp / e2u(:, :)
    r1_e1v(:, :) = 1._wp / e1v(:, :)
    r1_e2v(:, :) = 1._wp / e2v(:, :)
    r1_e1f(:, :) = 1._wp / e1f(:, :)
    r1_e2f(:, :) = 1._wp / e2f(:, :)
    e1e2t(:, :) = e1t(:, :) * e2t(:, :)
    r1_e1e2t(:, :) = 1._wp / e1e2t(:, :)
    e1e2f(:, :) = e1f(:, :) * e2f(:, :)
    r1_e1e2f(:, :) = 1._wp / e1e2f(:, :)
    !$ACC END KERNELS
    IF (ie1e2u_v == 0) THEN
      IF (lwp) WRITE(numout, *) '          u- & v-surfaces calculated as e1 e2 product'
      !$ACC KERNELS
      e1e2u(:, :) = e1u(:, :) * e2u(:, :)
      e1e2v(:, :) = e1v(:, :) * e2v(:, :)
      !$ACC END KERNELS
    ELSE
      CALL profile_psy_data2 % PreStart('dom_hgr', 'r2', 0, 0)
      IF (lwp) WRITE(numout, *) '          u- & v-surfaces have been read in "mesh_mask" file:'
      IF (lwp) WRITE(numout, *) '                     grid size reduction in strait(s) is used'
      CALL profile_psy_data2 % PostEnd
    END IF
    !$ACC KERNELS
    r1_e1e2u(:, :) = 1._wp / e1e2u(:, :)
    r1_e1e2v(:, :) = 1._wp / e1e2v(:, :)
    e2_e1u(:, :) = e2u(:, :) / e1u(:, :)
    e1_e2v(:, :) = e1v(:, :) / e2v(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('dom_hgr', 'r3', 0, 0)
    IF (ln_timing) CALL timing_stop('dom_hgr')
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE dom_hgr
  SUBROUTINE hgr_read(plamt, plamu, plamv, plamf, pphit, pphiu, pphiv, pphif, kff, pff_f, pff_t, pe1t, pe1u, pe1v, pe1f, pe2t, pe2u, pe2v, pe2f, ke1e2u_v, pe1e2u, pe1e2v)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: plamt, plamu, plamv, plamf
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pphit, pphiu, pphiv, pphif
    INTEGER, INTENT(OUT) :: kff
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pff_f, pff_t
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pe1t, pe1u, pe1v, pe1f
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pe2t, pe2u, pe2v, pe2f
    INTEGER, INTENT(OUT) :: ke1e2u_v
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: pe1e2u, pe1e2v
    INTEGER :: inum
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('hgr_read', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   hgr_read : read the horizontal coordinates in mesh_mask'
      WRITE(numout, *) '   ~~~~~~~~      jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpk = ', jpk
    END IF
    CALL iom_open(cn_domcfg, inum)
    CALL iom_get(inum, jpdom_data, 'glamt', plamt, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'glamu', plamu, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'glamv', plamv, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'glamf', plamf, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'gphit', pphit, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'gphiu', pphiu, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'gphiv', pphiv, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'gphif', pphif, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e1t', pe1t, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e1u', pe1u, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e1v', pe1v, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e1f', pe1f, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e2t', pe2t, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e2u', pe2u, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e2v', pe2v, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e2f', pe2f, lrowattr = ln_use_jattr)
    IF (iom_varid(inum, 'ff_f', ldstop = .FALSE.) > 0 .AND. iom_varid(inum, 'ff_t', ldstop = .FALSE.) > 0) THEN
      IF (lwp) WRITE(numout, *) '           Coriolis factor at f- and t-points read in ', TRIM(cn_domcfg), ' file'
      CALL iom_get(inum, jpdom_data, 'ff_f', pff_f, lrowattr = ln_use_jattr)
      CALL iom_get(inum, jpdom_data, 'ff_t', pff_t, lrowattr = ln_use_jattr)
      kff = 1
    ELSE
      kff = 0
    END IF
    IF (iom_varid(inum, 'e1e2u', ldstop = .FALSE.) > 0) THEN
      IF (lwp) WRITE(numout, *) '           e1e2u & e1e2v read in ', TRIM(cn_domcfg), ' file'
      CALL iom_get(inum, jpdom_data, 'e1e2u', pe1e2u, lrowattr = ln_use_jattr)
      CALL iom_get(inum, jpdom_data, 'e1e2v', pe1e2v, lrowattr = ln_use_jattr)
      ke1e2u_v = 1
    ELSE
      ke1e2u_v = 0
    END IF
    CALL iom_close(inum)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE hgr_read
END MODULE domhgr