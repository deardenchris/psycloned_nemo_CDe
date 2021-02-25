MODULE dynspg_exp
  USE oce
  USE dom_oce
  USE sbc_oce
  USE phycst
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE prtctl
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_spg_exp
  CONTAINS
  SUBROUTINE dyn_spg_exp(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_spg_exp : surface pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   (explicit free surface)'
      spgu(:, :) = 0._wp
      spgv(:, :) = 0._wp
      IF (.NOT. ln_linssh .AND. lwp) WRITE(numout, FMT = *) '      non linear free surface: spg is included in dynhpg'
    END IF
    IF (ln_linssh) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          spgu(ji, jj) = - grav * (sshn(ji + 1, jj) - sshn(ji, jj)) * r1_e1u(ji, jj)
          spgv(ji, jj) = - grav * (sshn(ji, jj + 1) - sshn(ji, jj)) * r1_e2v(ji, jj)
        END DO
      END DO
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ua(ji, jj, jk) = ua(ji, jj, jk) + spgu(ji, jj)
            va(ji, jj, jk) = va(ji, jj, jk) + spgv(ji, jj)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
  END SUBROUTINE dyn_spg_exp
END MODULE dynspg_exp