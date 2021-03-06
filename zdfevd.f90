MODULE zdfevd
  USE oce
  USE dom_oce
  USE zdf_oce
  USE trd_oce
  USE trdtra
  USE in_out_manager
  USE iom
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_evd
  CONTAINS
  SUBROUTINE zdf_evd(kt, p_avm, p_avt)
    INTEGER, INTENT(IN   ) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm, p_avt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zavt_evd, zavm_evd
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'zdf_evd : Enhanced Vertical Diffusion (evd)'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
      IF (lwp) WRITE(numout, FMT = *)
    END IF
    !$ACC KERNELS
    zavt_evd(:, :, :) = p_avt(:, :, :)
    !$ACC END KERNELS
    SELECT CASE (nn_evdm)
    CASE (1)
      !$ACC KERNELS
      zavm_evd(:, :, :) = p_avm(:, :, :)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (MIN(rn2(ji, jj, jk), rn2b(ji, jj, jk)) <= - 1.E-12) THEN
              p_avt(ji, jj, jk) = rn_evd * wmask(ji, jj, jk)
              p_avm(ji, jj, jk) = rn_evd * wmask(ji, jj, jk)
            END IF
          END DO
        END DO
      END DO
      zavm_evd(:, :, :) = p_avm(:, :, :) - zavm_evd(:, :, :)
      !$ACC END KERNELS
      CALL iom_put("avm_evd", zavm_evd)
    CASE DEFAULT
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF (MIN(rn2(ji, jj, jk), rn2b(ji, jj, jk)) <= - 1.E-12) p_avt(ji, jj, jk) = rn_evd * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END SELECT
    !$ACC KERNELS
    zavt_evd(:, :, :) = p_avt(:, :, :) - zavt_evd(:, :, :)
    !$ACC END KERNELS
    CALL iom_put("avt_evd", zavt_evd)
    IF (l_trdtra) CALL trd_tra(kt, 'TRA', jp_tem, jptra_evd, zavt_evd)
  END SUBROUTINE zdf_evd
END MODULE zdfevd