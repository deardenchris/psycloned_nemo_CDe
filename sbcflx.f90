MODULE sbcflx
  USE oce
  USE dom_oce
  USE sbc_oce
  USE sbcdcy
  USE phycst
  USE fldread
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_flx
  INTEGER, PARAMETER :: jpfld = 5
  INTEGER, PARAMETER :: jp_utau = 1
  INTEGER, PARAMETER :: jp_vtau = 2
  INTEGER, PARAMETER :: jp_qtot = 3
  INTEGER, PARAMETER :: jp_qsr = 4
  INTEGER, PARAMETER :: jp_emp = 5
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf
  CONTAINS
  SUBROUTINE sbc_flx(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jf
    INTEGER :: ierror
    INTEGER :: ios
    REAL(KIND = wp) :: zfact
    REAL(KIND = wp) :: zrhoa = 1.22
    REAL(KIND = wp) :: zcdrag = 1.5E-3
    REAL(KIND = wp) :: ztx, zty, zmod, zcoef
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N), DIMENSION(jpfld) :: slf_i
    TYPE(FLD_N) :: sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp
    NAMELIST /namsbc_flx/ cn_dir, sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    IF (kt == nit000) THEN
      CALL profile_psy_data0 % PreStart('sbc_flx', 'r0', 0, 0)
      REWIND(UNIT = numnam_ref)
      READ(numnam_ref, namsbc_flx, IOSTAT = ios, ERR = 901)
901   IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_flx in reference namelist', lwp)
      REWIND(UNIT = numnam_cfg)
      READ(numnam_cfg, namsbc_flx, IOSTAT = ios, ERR = 902)
902   IF (ios > 0) CALL ctl_nam(ios, 'namsbc_flx in configuration namelist', lwp)
      IF (lwm) WRITE(numond, namsbc_flx)
      IF (ln_dm2dc .AND. sn_qsr % nfreqh /= 24) CALL ctl_stop('sbc_blk_core: ln_dm2dc can be activated only with daily short-wave &
&forcing')
      slf_i(jp_utau) = sn_utau
      slf_i(jp_vtau) = sn_vtau
      slf_i(jp_qtot) = sn_qtot
      slf_i(jp_qsr) = sn_qsr
      slf_i(jp_emp) = sn_emp
      ALLOCATE(sf(jpfld), STAT = ierror)
      CALL profile_psy_data0 % PostEnd
      IF (ierror > 0) THEN
        CALL ctl_stop('sbc_flx: unable to allocate sf structure')
        RETURN
      END IF
      CALL profile_psy_data1 % PreStart('sbc_flx', 'r1', 0, 0)
      DO ji = 1, jpfld
        ALLOCATE(sf(ji) % fnow(jpi, jpj, 1))
        IF (slf_i(ji) % ln_tint) ALLOCATE(sf(ji) % fdta(jpi, jpj, 1, 2))
      END DO
      CALL fld_fill(sf, slf_i, cn_dir, 'sbc_flx', 'flux formulation for ocean surface boundary condition', 'namsbc_flx')
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      sfx(:, :) = 0.0_wp
      !$ACC END KERNELS
    END IF
    CALL fld_read(kt, nn_fsbc, sf)
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      CALL profile_psy_data2 % PreStart('sbc_flx', 'r2', 0, 0)
      IF (ln_dm2dc) THEN
        qsr(:, :) = sbc_dcy(sf(jp_qsr) % fnow(:, :, 1))
      ELSE
        qsr(:, :) = sf(jp_qsr) % fnow(:, :, 1)
      END IF
      DO jj = 1, jpj
        DO ji = 1, jpi
          utau(ji, jj) = sf(jp_utau) % fnow(ji, jj, 1)
          vtau(ji, jj) = sf(jp_vtau) % fnow(ji, jj, 1)
          qns(ji, jj) = sf(jp_qtot) % fnow(ji, jj, 1) - sf(jp_qsr) % fnow(ji, jj, 1)
          emp(ji, jj) = sf(jp_emp) % fnow(ji, jj, 1)
        END DO
      END DO
      CALL profile_psy_data2 % PostEnd
      !$ACC KERNELS
      qns(:, :) = qns(:, :) - emp(:, :) * sst_m(:, :) * rcp
      qns(:, :) = qns(:, :) * tmask(:, :, 1)
      emp(:, :) = emp(:, :) * tmask(:, :, 1)
      zcoef = 1. / (zrhoa * zcdrag)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ztx = utau(ji - 1, jj) + utau(ji, jj)
          zty = vtau(ji, jj - 1) + vtau(ji, jj)
          zmod = 0.5 * SQRT(ztx * ztx + zty * zty)
          taum(ji, jj) = zmod
          wndm(ji, jj) = SQRT(zmod * zcoef)
        END DO
      END DO
      taum(:, :) = taum(:, :) * tmask(:, :, 1)
      wndm(:, :) = wndm(:, :) * tmask(:, :, 1)
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('sbc_flx', 'r3', 0, 0)
      CALL lbc_lnk('sbcflx', taum(:, :), 'T', 1.)
      CALL lbc_lnk('sbcflx', wndm(:, :), 'T', 1.)
      IF (nitend - nit000 <= 100 .AND. lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '        read daily momentum, heat and freshwater fluxes OK'
        DO jf = 1, jpfld
          IF (jf == jp_utau .OR. jf == jp_vtau) zfact = 1.
          IF (jf == jp_qtot .OR. jf == jp_qsr) zfact = 0.1
          IF (jf == jp_emp) zfact = 86400.
          WRITE(numout, FMT = *)
          WRITE(numout, FMT = *) ' day: ', ndastp, TRIM(sf(jf) % clvar), ' * ', zfact
        END DO
      END IF
      CALL profile_psy_data3 % PostEnd
    END IF
  END SUBROUTINE sbc_flx
END MODULE sbcflx