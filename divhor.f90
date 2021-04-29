MODULE divhor
  USE oce
  USE dom_oce
  USE sbc_oce, ONLY: ln_rnf, ln_isf
  USE sbcrnf
  USE sbcisf
  USE iscplhsb
  USE iscplini
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: div_hor
  CONTAINS
  SUBROUTINE div_hor(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zraur, zdep
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('div_hor', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('div_hor')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'div_hor : horizontal velocity divergence '
      IF (lwp) WRITE(numout, *) '~~~~~~~   '
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          hdivn(ji, jj, jk) = (e2u(ji, jj) * e3u_n(ji, jj, jk) * un(ji, jj, jk) - e2u(ji - 1, jj) * e3u_n(ji - 1, jj, jk) * un(ji - 1, jj, jk) + e1v(ji, jj) * e3v_n(ji, jj, jk) * vn(ji, jj, jk) - e1v(ji, jj - 1) * e3v_n(ji, jj - 1, jk) * vn(ji, jj - 1, jk)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('div_hor', 'r1', 0, 0)
    IF (ln_rnf) CALL sbc_rnf_div(hdivn)
    IF (ln_isf) CALL sbc_isf_div(hdivn)
    IF (ln_iscpl .AND. ln_hsb) CALL iscpl_div(hdivn)
    CALL lbc_lnk('divhor', hdivn, 'T', 1.)
    IF (ln_timing) CALL timing_stop('div_hor')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE div_hor
END MODULE divhor