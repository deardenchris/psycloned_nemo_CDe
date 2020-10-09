MODULE step_diu
  USE diurnal_bulk
  USE cool_skin
  USE iom
  USE sbc_oce
  USE sbcmod
  USE diaobs
  USE oce
  USE daymod
  USE restart
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stp_diurnal
  CONTAINS
  SUBROUTINE stp_diurnal(kstp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kstp
    INTEGER :: jk
    INTEGER :: indic
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z_fvel_bkginc, z_hflux_bkginc
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('stp_diurnal', 'r0', 0, 0)
    IF (ln_diurnal_only) THEN
      indic = 0
      IF (kstp /= nit000) CALL day(kstp)
      CALL iom_setkt(kstp - nit000 + 1, cxios_context)
      IF (ln_crs) THEN
        CALL iom_setkt(kstp - nit000 + 1, TRIM(cxios_context) // "_crs")
      END IF
      CALL sbc(kstp)
    END IF
    IF (.NOT. ln_diurnal) CALL ctl_stop("stp_diurnal: ln_diurnal not set")
    IF (.NOT. ln_blk) CALL ctl_stop("stp_diurnal: diurnal flux processing only implemented for bulk forcing")
    CALL diurnal_sst_coolskin_step(qns, taum, rhop(:, :, 1), rdt)
    CALL iom_put("sst_wl", x_dsst)
    CALL iom_put("sst_cs", x_csdsst)
    CALL diurnal_sst_takaya_step(kstp, qsr, qns, taum, rhop(:, :, 1), rdt)
    IF (ln_diurnal_only) THEN
      IF (ln_diaobs) CALL dia_obs(kstp)
      IF (kstp == nit000) CALL iom_close(numror)
      IF (lrst_oce) CALL rst_write(kstp)
      IF (ln_timing .AND. kstp == nit000) CALL timing_reset
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE stp_diurnal
END MODULE step_diu