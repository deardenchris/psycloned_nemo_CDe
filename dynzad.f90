MODULE dynzad
  USE oce
  USE dom_oce
  USE sbc_oce
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE lib_mpp
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_zad
  CONTAINS
  SUBROUTINE dyn_zad(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zua, zva
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zww
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwuw, zwvw
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdu, ztrdv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('dyn_zad', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('dyn_zad')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_zad : 2nd order vertical advection scheme'
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    DO jk = 2, jpkm1
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          zww(ji, jj) = 0.25_wp * e1e2t(ji, jj) * wn(ji, jj, jk)
        END DO
      END DO
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwuw(ji, jj, jk) = (zww(ji + 1, jj) + zww(ji, jj)) * (un(ji, jj, jk - 1) - un(ji, jj, jk))
          zwvw(ji, jj, jk) = (zww(ji, jj + 1) + zww(ji, jj)) * (vn(ji, jj, jk - 1) - vn(ji, jj, jk))
        END DO
      END DO
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zwuw(ji, jj, 1) = 0._wp
        zwvw(ji, jj, 1) = 0._wp
        zwuw(ji, jj, jpk) = 0._wp
        zwvw(ji, jj, jpk) = 0._wp
      END DO
    END DO
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zwuw(ji, jj, jk) + zwuw(ji, jj, jk + 1)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (zwvw(ji, jj, jk) + zwvw(ji, jj, jk + 1)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (l_trddyn) THEN
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dyn_zad', 'r1', 0, 0)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_zad, kt)
      DEALLOCATE(ztrdu, ztrdv)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('dyn_zad', 'r2', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' zad  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_zad')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE dyn_zad
END MODULE dynzad