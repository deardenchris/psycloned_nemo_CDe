MODULE trdmxl_rst
  USE dom_oce
  USE trd_oce
  USE in_out_manager
  USE iom
  USE restart
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_mxl_rst_read
  PUBLIC :: trd_mxl_rst_write
  INTEGER :: nummxlw
  CONTAINS
  SUBROUTINE trd_mxl_rst_write(kt)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = 35) :: charout
    INTEGER :: jk
    CHARACTER(LEN = 20) :: clkt
    CHARACTER(LEN = 50) :: clname
    CHARACTER(LEN = 256) :: clpath
    IF (kt == nitrst - 1 .OR. nstock == 1 .OR. (kt == nitend .AND. MOD(nitend - 1, nstock) == 0)) THEN
      IF (nitrst > 999999999) THEN
        WRITE(clkt, *) nitrst
      ELSE
        WRITE(clkt, '(i8.8)') nitrst
      END IF
      clname = TRIM(cexper) // "_" // TRIM(ADJUSTL(clkt)) // "_" // TRIM(cn_trdrst_out)
      clpath = TRIM(cn_ocerst_outdir)
      IF (clpath(LEN_TRIM(clpath) :) /= '/') clpath = TRIM(clpath) // '/'
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '             open ocean restart_mxl NetCDF file: ' // clname
        IF (kt == nitrst - 1) THEN
          WRITE(numout, *) '             kt = nitrst - 1 = ', kt, ' date= ', ndastp
        ELSE
          WRITE(numout, *) '             kt = ', kt, ' date= ', ndastp
        END IF
      END IF
      CALL iom_open(TRIM(clpath) // TRIM(clname), nummxlw, ldwrt = .TRUE.)
    END IF
    IF (kt == nitrst .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'trdmxl_rst: output for ML diags. restart, with trd_mxl_rst_write routine kt =', kt
      WRITE(numout, *) '~~~~~~~~~~'
      WRITE(numout, *)
    END IF
    IF (ln_trdmxl_instant) THEN
      CALL iom_rstput(kt, nitrst, nummxlw, 'tmlbb', tmlbb)
      CALL iom_rstput(kt, nitrst, nummxlw, 'tmlbn', tmlbn)
      CALL iom_rstput(kt, nitrst, nummxlw, 'tmlatfb', tmlatfb)
      CALL iom_rstput(kt, nitrst, nummxlw, 'smlbb', smlbb)
      CALL iom_rstput(kt, nitrst, nummxlw, 'smlbn', smlbn)
      CALL iom_rstput(kt, nitrst, nummxlw, 'smlatfb', smlatfb)
    ELSE
      CALL iom_rstput(kt, nitrst, nummxlw, 'hmxlbn', hmxlbn)
      CALL iom_rstput(kt, nitrst, nummxlw, 'tmlbn', tmlbn)
      CALL iom_rstput(kt, nitrst, nummxlw, 'tml_sumb', tml_sumb)
      DO jk = 1, jpltrd
        IF (jk < 10) THEN
          WRITE(charout, FMT = "('tmltrd_csum_ub_', I1)") jk
        ELSE
          WRITE(charout, FMT = "('tmltrd_csum_ub_', I2)") jk
        END IF
        CALL iom_rstput(kt, nitrst, nummxlw, charout, tmltrd_csum_ub(:, :, jk))
      END DO
      CALL iom_rstput(kt, nitrst, nummxlw, 'tmltrd_atf_sumb', tmltrd_atf_sumb)
      CALL iom_rstput(kt, nitrst, nummxlw, 'smlbn', smlbn)
      CALL iom_rstput(kt, nitrst, nummxlw, 'sml_sumb', sml_sumb)
      DO jk = 1, jpltrd
        IF (jk < 10) THEN
          WRITE(charout, FMT = "('smltrd_csum_ub_', I1)") jk
        ELSE
          WRITE(charout, FMT = "('smltrd_csum_ub_', I2)") jk
        END IF
        CALL iom_rstput(kt, nitrst, nummxlw, charout, smltrd_csum_ub(:, :, jk))
      END DO
      CALL iom_rstput(kt, nitrst, nummxlw, 'smltrd_atf_sumb', smltrd_atf_sumb)
    END IF
    IF (kt == nitrst) THEN
      CALL iom_close(nummxlw)
      lrst_oce = .FALSE.
    END IF
  END SUBROUTINE trd_mxl_rst_write
  SUBROUTINE trd_mxl_rst_read
    INTEGER :: inum
    CHARACTER(LEN = 35) :: charout
    INTEGER :: jk
    LOGICAL :: llok
    CHARACTER(LEN = 256) :: clpath
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) ' trd_mxl_rst_read : read the NetCDF mixed layer trend restart file'
      WRITE(numout, *) ' ~~~~~~~~~~~~~~~~'
    END IF
    clpath = TRIM(cn_ocerst_indir)
    IF (clpath(LEN_TRIM(clpath) :) /= '/') clpath = TRIM(clpath) // '/'
    CALL iom_open(TRIM(clpath) // TRIM(cn_trdrst_in), inum)
    IF (ln_trdmxl_instant) THEN
      CALL iom_get(inum, jpdom_autoglo, 'tmlbb', tmlbb)
      CALL iom_get(inum, jpdom_autoglo, 'tmlbn', tmlbn)
      CALL iom_get(inum, jpdom_autoglo, 'tmlatfb', tmlatfb)
      CALL iom_get(inum, jpdom_autoglo, 'smlbb', smlbb)
      CALL iom_get(inum, jpdom_autoglo, 'smlbn', smlbn)
      CALL iom_get(inum, jpdom_autoglo, 'smlatfb', smlatfb)
    ELSE
      CALL iom_get(inum, jpdom_autoglo, 'hmxlbn', hmxlbn)
      CALL iom_get(inum, jpdom_autoglo, 'tmlbn', tmlbn)
      CALL iom_get(inum, jpdom_autoglo, 'tml_sumb', tml_sumb)
      DO jk = 1, jpltrd
        IF (jk < 10) THEN
          WRITE(charout, FMT = "('tmltrd_csum_ub_', I1)") jk
        ELSE
          WRITE(charout, FMT = "('tmltrd_csum_ub_', I2)") jk
        END IF
        CALL iom_get(inum, jpdom_autoglo, charout, tmltrd_csum_ub(:, :, jk))
      END DO
      CALL iom_get(inum, jpdom_autoglo, 'tmltrd_atf_sumb', tmltrd_atf_sumb)
      CALL iom_get(inum, jpdom_autoglo, 'smlbn', smlbn)
      CALL iom_get(inum, jpdom_autoglo, 'sml_sumb', sml_sumb)
      DO jk = 1, jpltrd
        IF (jk < 10) THEN
          WRITE(charout, FMT = "('smltrd_csum_ub_', I1)") jk
        ELSE
          WRITE(charout, FMT = "('smltrd_csum_ub_', I2)") jk
        END IF
        CALL iom_get(inum, jpdom_autoglo, charout, smltrd_csum_ub(:, :, jk))
      END DO
      CALL iom_get(inum, jpdom_autoglo, 'smltrd_atf_sumb', smltrd_atf_sumb)
      CALL iom_close(inum)
    END IF
  END SUBROUTINE trd_mxl_rst_read
END MODULE trdmxl_rst