MODULE dynldf_lap_blp
  USE oce
  USE dom_oce
  USE ldfdyn
  USE ldfslp
  USE zdf_oce
  USE in_out_manager
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_ldf_lap
  PUBLIC :: dyn_ldf_blp
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zulap, zvlap
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zcur, zdiv
  CONTAINS
  SUBROUTINE dyn_ldf_lap(kt, pub, pvb, pua, pva, kpass)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kpass
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pub, pvb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zsign
    REAL(KIND = wp) :: zua, zva
    IF (kt == nit000) THEN
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) 'dyn_ldf : iso-level harmonic (laplacian) operator, pass=', kpass
        WRITE(numout, FMT = *) '~~~~~~~ '
      END IF
      IF (kpass == 1) ALLOCATE(zcur(jpi, jpj), zdiv(jpi, jpj))
    END IF
    IF (kpass == 1) THEN
      zsign = 1._wp
    ELSE
      zsign = - 1._wp
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      DO jj = 2, jpj
        DO ji = 2, jpi
          zcur(ji - 1, jj - 1) = ahmf(ji - 1, jj - 1, jk) * e3f_n(ji - 1, jj - 1, jk) * r1_e1e2f(ji - 1, jj - 1) * (e2v(ji, jj - &
&1) * pvb(ji, jj - 1, jk) - e2v(ji - 1, jj - 1) * pvb(ji - 1, jj - 1, jk) - e1u(ji - 1, jj) * pub(ji - 1, jj, jk) + e1u(ji - 1, jj &
&- 1) * pub(ji - 1, jj - 1, jk))
          zdiv(ji, jj) = ahmt(ji, jj, jk) * r1_e1e2t(ji, jj) / e3t_b(ji, jj, jk) * (e2u(ji, jj) * e3u_b(ji, jj, jk) * pub(ji, jj, &
&jk) - e2u(ji - 1, jj) * e3u_b(ji - 1, jj, jk) * pub(ji - 1, jj, jk) + e1v(ji, jj) * e3v_b(ji, jj, jk) * pvb(ji, jj, jk) - e1v(ji, &
&jj - 1) * e3v_b(ji, jj - 1, jk) * pvb(ji, jj - 1, jk))
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pua(ji, jj, jk) = pua(ji, jj, jk) + zsign * (- (zcur(ji, jj) - zcur(ji, jj - 1)) * r1_e2u(ji, jj) / e3u_n(ji, jj, jk) + &
&(zdiv(ji + 1, jj) - zdiv(ji, jj)) * r1_e1u(ji, jj))
          pva(ji, jj, jk) = pva(ji, jj, jk) + zsign * ((zcur(ji, jj) - zcur(ji - 1, jj)) * r1_e1v(ji, jj) / e3v_n(ji, jj, jk) + &
&(zdiv(ji, jj + 1) - zdiv(ji, jj)) * r1_e2v(ji, jj))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
  END SUBROUTINE dyn_ldf_lap
  SUBROUTINE dyn_ldf_blp(kt, pub, pvb, pua, pva)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pub, pvb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pua, pva
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_ldf_blp : bilaplacian operator momentum '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      ALLOCATE(zulap(jpi, jpj, jpk), zvlap(jpi, jpj, jpk))
    END IF
    zulap(:, :, :) = 0._wp
    zvlap(:, :, :) = 0._wp
    CALL dyn_ldf_lap(kt, pub, pvb, zulap, zvlap, 1)
    CALL lbc_lnk_multi('dynldf_lap_blp', zulap, 'U', - 1., zvlap, 'V', - 1.)
    CALL dyn_ldf_lap(kt, zulap, zvlap, pua, pva, 2)
  END SUBROUTINE dyn_ldf_blp
END MODULE dynldf_lap_blp