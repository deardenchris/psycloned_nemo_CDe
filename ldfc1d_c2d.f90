MODULE ldfc1d_c2d
  USE oce
  USE dom_oce
  USE phycst
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ldf_c1d
  PUBLIC :: ldf_c2d
  REAL(KIND = wp) :: r1_2 = 0.5_wp
  REAL(KIND = wp) :: r1_4 = 0.25_wp
  REAL(KIND = wp) :: r1_12 = 1._wp / 12._wp
  CONTAINS
  SUBROUTINE ldf_c1d(cd_type, pahs1, pahs2, pah1, pah2)
    CHARACTER(LEN = 3), INTENT(IN) :: cd_type
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pahs1, pahs2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(INOUT) :: pah1, pah2
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zh, zc, zdep1
    REAL(KIND = wp) :: zw, zdep2
    REAL(KIND = wp) :: zratio
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) '   ldf_c1d : set a given profile to eddy mixing coefficients'
    zratio = 0.25_wp
    zh = 500._wp
    zw = 1._wp / 200._wp
    zc = (1._wp - zratio) / (1._wp + TANH(zh * zw))
    SELECT CASE (cd_type)
    CASE ('DYN')
      DO jk = jpkm1, 1, - 1
        pah1(:, :, jk) = pahs1(:, :) * (zratio + zc * (1._wp + TANH(- (gdept_0(:, :, jk) - zh) * zw)))
      END DO
      DO jk = jpkm1, 1, - 1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdep2 = (gdept_0(ji, jj + 1, jk) + gdept_0(ji + 1, jj + 1, jk) + gdept_0(ji, jj, jk) + gdept_0(ji + 1, jj, jk)) * r1_4
            pah2(ji, jj, jk) = pahs2(ji, jj) * (zratio + zc * (1._wp + TANH(- (zdep2 - zh) * zw)))
          END DO
        END DO
      END DO
      CALL lbc_lnk('ldfc1d_c2d', pah2, 'F', 1.)
    CASE ('TRA')
      DO jk = jpkm1, 1, - 1
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdep1 = (gdept_0(ji, jj, jk) + gdept_0(ji + 1, jj, jk)) * 0.5_wp
            zdep2 = (gdept_0(ji, jj, jk) + gdept_0(ji, jj + 1, jk)) * 0.5_wp
            pah1(ji, jj, jk) = pahs1(ji, jj) * (zratio + zc * (1._wp + TANH(- (zdep1 - zh) * zw)))
            pah2(ji, jj, jk) = pahs2(ji, jj) * (zratio + zc * (1._wp + TANH(- (zdep2 - zh) * zw)))
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('ldfc1d_c2d', pah1, 'U', 1., pah2, 'V', 1.)
    CASE DEFAULT
      CALL ctl_stop('ldf_c1d: ', cd_type, ' Unknown, i.e. /= DYN or TRA')
    END SELECT
  END SUBROUTINE ldf_c1d
  SUBROUTINE ldf_c2d(cd_type, pUfac, knn, pah1, pah2)
    CHARACTER(LEN = 3), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: pUfac
    INTEGER, INTENT(IN) :: knn
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pah1, pah2
    INTEGER :: ji, jj, jk
    INTEGER :: inn
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) '   ldf_c2d :   aht = Ufac * max(e1,e2)   with Ufac = ', pUfac, ' m/s'
    SELECT CASE (cd_type)
    CASE ('DYN')
      DO jj = 1, jpj
        DO ji = 1, jpi
          pah1(ji, jj, 1) = pUfac * MAX(e1t(ji, jj), e2t(ji, jj)) ** knn
          pah2(ji, jj, 1) = pUfac * MAX(e1f(ji, jj), e2f(ji, jj)) ** knn
        END DO
      END DO
    CASE ('TRA')
      DO jj = 1, jpj
        DO ji = 1, jpi
          pah1(ji, jj, 1) = pUfac * MAX(e1u(ji, jj), e2u(ji, jj)) ** knn
          pah2(ji, jj, 1) = pUfac * MAX(e1v(ji, jj), e2v(ji, jj)) ** knn
        END DO
      END DO
    CASE DEFAULT
      CALL ctl_stop('ldf_c2d: ', cd_type, ' Unknown, i.e. /= DYN or TRA')
    END SELECT
    DO jk = 2, jpkm1
      pah1(:, :, jk) = pah1(:, :, 1)
      pah2(:, :, jk) = pah2(:, :, 1)
    END DO
  END SUBROUTINE ldf_c2d
END MODULE ldfc1d_c2d