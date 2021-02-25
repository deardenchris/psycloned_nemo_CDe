MODULE zdfswm
  USE dom_oce
  USE zdf_oce
  USE sbc_oce
  USE sbcwave
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_swm
  PUBLIC :: zdf_swm_init
  CONTAINS
  SUBROUTINE zdf_swm(kt, p_avm, p_avt, p_avs)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avt, p_avs
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoef, zqb
    zcoef = 1._wp * 0.353553_wp
    !$OMP parallel default(shared), private(ji,jj,jk,zqb)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zqb = zcoef * hsw(ji, jj) * tsd2d(ji, jj) * EXP(- 3. * wnum(ji, jj) * gdepw_n(ji, jj, jk)) * wmask(ji, jj, jk)
          p_avt(ji, jj, jk) = p_avt(ji, jj, jk) + zqb
          p_avs(ji, jj, jk) = p_avs(ji, jj, jk) + zqb
          p_avm(ji, jj, jk) = p_avm(ji, jj, jk) + zqb
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
  END SUBROUTINE zdf_swm
  SUBROUTINE zdf_swm_init
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'zdf_swm_init : surface wave-driven mixing'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    END IF
    IF (.NOT. ln_wave .OR. .NOT. ln_sdw) CALL ctl_stop('zdf_swm_init: ln_zdfswm=T but ln_wave and ln_sdw /= T')
  END SUBROUTINE zdf_swm_init
END MODULE zdfswm