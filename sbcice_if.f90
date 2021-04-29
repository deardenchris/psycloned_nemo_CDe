MODULE sbcice_if
  USE oce
  USE dom_oce
  USE phycst
  USE eosbn2
  USE sbc_oce
  USE sbc_ice, ONLY: a_i
  USE in_out_manager
  USE iom
  USE fldread
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_ice_if
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ice
  CONTAINS
  SUBROUTINE sbc_ice_if(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    INTEGER :: ierror
    INTEGER :: ios
    REAL(KIND = wp) :: ztrp, zsice, zt_fzp, zfr_obs
    REAL(KIND = wp) :: zqri, zqrj, zqrp, zqi
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N) :: sn_ice
    NAMELIST /namsbc_iif/ cn_dir, sn_ice
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('sbc_ice_if', 'r0', 0, 0)
    IF (kt == nit000) THEN
      REWIND(UNIT = numnam_ref)
      READ(numnam_ref, namsbc_iif, IOSTAT = ios, ERR = 901)
901   IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_iif in reference namelist', lwp)
      REWIND(UNIT = numnam_cfg)
      READ(numnam_cfg, namsbc_iif, IOSTAT = ios, ERR = 902)
902   IF (ios > 0) CALL ctl_nam(ios, 'namsbc_iif in configuration namelist', lwp)
      IF (lwm) WRITE(numond, namsbc_iif)
      ALLOCATE(sf_ice(1), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ice_if: unable to allocate sf_ice structure')
      ALLOCATE(sf_ice(1) % fnow(jpi, jpj, 1))
      IF (sn_ice % ln_tint) ALLOCATE(sf_ice(1) % fdta(jpi, jpj, 1, 2))
      CALL fld_fill(sf_ice, (/sn_ice/), cn_dir, 'sbc_ice_if', 'ice-if sea-ice model', 'namsbc_iif')
    END IF
    CALL fld_read(kt, nn_fsbc, sf_ice)
    CALL profile_psy_data0 % PostEnd
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      CALL profile_psy_data1 % PreStart('sbc_ice_if', 'r1', 0, 0)
      ztrp = - 40.
      zsice = - 0.04 / 0.8
      CALL eos_fzp(sss_m(:, :), fr_i(:, :))
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      fr_i(:, :) = fr_i(:, :) * tmask(:, :, 1)
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('sbc_ice_if', 'r2', 0, 0)
      IF (ln_cpl) a_i(:, :, 1) = fr_i(:, :)
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zt_fzp = fr_i(ji, jj)
          zfr_obs = sf_ice(1) % fnow(ji, jj, 1)
          IF (sst_m(ji, jj) <= zt_fzp) THEN
            fr_i(ji, jj) = 1.E0
          ELSE
            fr_i(ji, jj) = 0.E0
          END IF
          tsn(ji, jj, 1, jp_tem) = MAX(tsn(ji, jj, 1, jp_tem), zt_fzp)
          qsr(ji, jj) = (1. - zfr_obs) * qsr(ji, jj)
          zqri = ztrp * (tsb(ji, jj, 1, jp_tem) - (zt_fzp - 1.))
          zqrj = ztrp * MIN(0., tsb(ji, jj, 1, jp_tem) - zt_fzp)
          zqrp = (zfr_obs * ((1. - fr_i(ji, jj)) * zqri + fr_i(ji, jj) * zqrj)) * tmask(ji, jj, 1)
          zqi = - 3. + SIGN(1._wp, ff_f(ji, jj))
          qns(ji, jj) = ((1. - zfr_obs) * qns(ji, jj) + zfr_obs * fr_i(ji, jj) * zqi) * tmask(ji, jj, 1) + zqrp
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END SUBROUTINE sbc_ice_if
END MODULE sbcice_if