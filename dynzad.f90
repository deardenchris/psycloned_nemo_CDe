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
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zua, zva
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zww
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwuw, zwvw
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdu, ztrdv
    IF (ln_timing) CALL timing_start('dyn_zad')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_zad : 2nd order vertical advection scheme'
    END IF
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
    END IF
    ! !$OMP parallel default(shared), private(ji,jj,jk) ! CDe causes runtime crash
    ! !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpj
        DO ji = 2, jpi
          zww(ji, jj) = 0.25_wp * e1e2t(ji, jj) * wn(ji, jj, jk)
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zwuw(ji, jj, jk) = (zww(ji + 1, jj) + zww(ji, jj)) * (un(ji, jj, jk - 1) - un(ji, jj, jk))
          zwvw(ji, jj, jk) = (zww(ji, jj + 1) + zww(ji, jj)) * (vn(ji, jj, jk - 1) - vn(ji, jj, jk))
        END DO
      END DO
    END DO
    ! !$OMP end do
    ! !$OMP end parallel
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zwuw(ji, jj, 1) = 0._wp
        zwvw(ji, jj, 1) = 0._wp
        zwuw(ji, jj, jpk) = 0._wp
        zwvw(ji, jj, jpk) = 0._wp
      END DO
    END DO
    !$OMP parallel default(shared), private(ji,jj,jk) ! CDe
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zwuw(ji, jj, jk) + zwuw(ji, jj, jk + 1)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (zwvw(ji, jj, jk) + zwvw(ji, jj, jk + 1)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (l_trddyn) THEN
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_zad, kt)
      DEALLOCATE(ztrdu, ztrdv)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' zad  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_zad')
  END SUBROUTINE dyn_zad
END MODULE dynzad
