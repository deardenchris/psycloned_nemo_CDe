MODULE sbcssr
  USE oce
  USE dom_oce
  USE sbc_oce
  USE phycst
  USE sbcrnf
  USE fldread
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_ssr
  PUBLIC :: sbc_ssr_init
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: erp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: qrp
  INTEGER, PUBLIC :: nn_sstr
  INTEGER, PUBLIC :: nn_sssr
  REAL(KIND = wp) :: rn_dqdt
  REAL(KIND = wp) :: rn_deds
  LOGICAL :: ln_sssr_bnd
  REAL(KIND = wp) :: rn_sssr_bnd
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: buffer
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_sst
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_sss
  CONTAINS
  SUBROUTINE sbc_ssr(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zerp
    REAL(KIND = wp) :: zqrp
    REAL(KIND = wp) :: zsrp
    REAL(KIND = wp) :: zerp_bnd
    INTEGER :: ierror
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N) :: sn_sst, sn_sss
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    IF (nn_sstr + nn_sssr /= 0) THEN
      CALL profile_psy_data0 % PreStart('sbc_ssr', 'r0', 0, 0)
      IF (nn_sstr == 1) CALL fld_read(kt, nn_fsbc, sf_sst)
      IF (nn_sssr >= 1) CALL fld_read(kt, nn_fsbc, sf_sss)
      CALL profile_psy_data0 % PostEnd
      IF (MOD(kt - 1, nn_fsbc) == 0) THEN
        IF (nn_sstr == 1) THEN
          !$ACC KERNELS
          !$ACC loop independent collapse(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              zqrp = rn_dqdt * (sst_m(ji, jj) - sf_sst(1) % fnow(ji, jj, 1)) * tmask(ji, jj, 1)
              qns(ji, jj) = qns(ji, jj) + zqrp
              qrp(ji, jj) = zqrp
            END DO
          END DO
          !$ACC END KERNELS
          CALL profile_psy_data1 % PreStart('sbc_ssr', 'r1', 0, 0)
          CALL iom_put("qrp", qrp)
          CALL profile_psy_data1 % PostEnd
        END IF
        IF (nn_sssr == 1) THEN
          !$ACC KERNELS
          zsrp = rn_deds / rday
          !$ACC loop independent collapse(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              zerp = zsrp * (1. - 2. * rnfmsk(ji, jj)) * (sss_m(ji, jj) - sf_sss(1) % fnow(ji, jj, 1)) * tmask(ji, jj, 1)
              sfx(ji, jj) = sfx(ji, jj) + zerp
              erp(ji, jj) = zerp / MAX(sss_m(ji, jj), 1.E-20)
            END DO
          END DO
          !$ACC END KERNELS
          CALL profile_psy_data2 % PreStart('sbc_ssr', 'r2', 0, 0)
          CALL iom_put("erp", erp)
          CALL profile_psy_data2 % PostEnd
        ELSE IF (nn_sssr == 2) THEN
          CALL profile_psy_data3 % PreStart('sbc_ssr', 'r3', 0, 0)
          zsrp = rn_deds / rday
          zerp_bnd = rn_sssr_bnd / rday
          DO jj = 1, jpj
            DO ji = 1, jpi
              zerp = zsrp * (1. - 2. * rnfmsk(ji, jj)) * (sss_m(ji, jj) - sf_sss(1) % fnow(ji, jj, 1)) / MAX(sss_m(ji, jj), 1.E-20) * tmask(ji, jj, 1)
              IF (ln_sssr_bnd) zerp = SIGN(1., zerp) * MIN(zerp_bnd, ABS(zerp))
              emp(ji, jj) = emp(ji, jj) + zerp
              qns(ji, jj) = qns(ji, jj) - zerp * rcp * sst_m(ji, jj)
              erp(ji, jj) = zerp
            END DO
          END DO
          CALL iom_put("erp", erp)
          CALL profile_psy_data3 % PostEnd
        END IF
      END IF
    END IF
  END SUBROUTINE sbc_ssr
  SUBROUTINE sbc_ssr_init
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zerp
    REAL(KIND = wp) :: zqrp
    REAL(KIND = wp) :: zsrp
    REAL(KIND = wp) :: zerp_bnd
    INTEGER :: ierror
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N) :: sn_sst, sn_sss
    NAMELIST /namsbc_ssr/ cn_dir, nn_sstr, nn_sssr, rn_dqdt, rn_deds, sn_sst, sn_sss, ln_sssr_bnd, rn_sssr_bnd
    INTEGER :: ios
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'sbc_ssr : SST and/or SSS damping term '
      WRITE(numout, *) '~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_ssr, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_ssr in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_ssr, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_ssr in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_ssr)
    IF (lwp) THEN
      WRITE(numout, *) '   Namelist namsbc_ssr :'
      WRITE(numout, *) '      SST restoring term (Yes=1)             nn_sstr        = ', nn_sstr
      WRITE(numout, *) '         dQ/dT (restoring magnitude on SST)     rn_dqdt     = ', rn_dqdt, ' W/m2/K'
      WRITE(numout, *) '      SSS damping term (Yes=1, salt   flux)  nn_sssr        = ', nn_sssr
      WRITE(numout, *) '                       (Yes=2, volume flux) '
      WRITE(numout, *) '         dE/dS (restoring magnitude on SST)     rn_deds     = ', rn_deds, ' mm/day'
      WRITE(numout, *) '         flag to bound erp term                 ln_sssr_bnd = ', ln_sssr_bnd
      WRITE(numout, *) '         ABS(Max./Min.) erp threshold           rn_sssr_bnd = ', rn_sssr_bnd, ' mm/day'
    END IF
    ALLOCATE(qrp(jpi, jpj), erp(jpi, jpj), STAT = ierror)
    IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate erp and qrp array')
    IF (nn_sstr == 1) THEN
      ALLOCATE(sf_sst(1), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate sf_sst structure')
      ALLOCATE(sf_sst(1) % fnow(jpi, jpj, 1), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate sf_sst now array')
      CALL fld_fill(sf_sst, (/sn_sst/), cn_dir, 'sbc_ssr', 'SST restoring term toward SST data', 'namsbc_ssr', no_print)
      IF (sf_sst(1) % ln_tint) ALLOCATE(sf_sst(1) % fdta(jpi, jpj, 1, 2), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate sf_sst data array')
    END IF
    IF (nn_sssr >= 1) THEN
      ALLOCATE(sf_sss(1), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate sf_sss structure')
      ALLOCATE(sf_sss(1) % fnow(jpi, jpj, 1), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate sf_sss now array')
      CALL fld_fill(sf_sss, (/sn_sss/), cn_dir, 'sbc_ssr', 'SSS restoring term toward SSS data', 'namsbc_ssr', no_print)
      IF (sf_sss(1) % ln_tint) ALLOCATE(sf_sss(1) % fdta(jpi, jpj, 1, 2), STAT = ierror)
      IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_ssr: unable to allocate sf_sss data array')
    END IF
    !$ACC KERNELS
    IF (nn_sstr /= 1) qrp(:, :) = 0._wp
    IF (nn_sssr /= 1 .OR. nn_sssr /= 2) erp(:, :) = 0._wp
    !$ACC END KERNELS
  END SUBROUTINE sbc_ssr_init
END MODULE sbcssr