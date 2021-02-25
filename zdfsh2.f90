MODULE zdfsh2
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_sh2
  CONTAINS
  SUBROUTINE zdf_sh2(pub, pvb, pun, pvn, p_avm, p_sh2)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pub, pvb, pun, pvn
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: p_avm
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: p_sh2
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsh2u, zsh2v
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zsh2u(ji, jj) = (p_avm(ji + 1, jj, jk) + p_avm(ji, jj, jk)) * (pun(ji, jj, jk - 1) - pun(ji, jj, jk)) * (pub(ji, jj, jk &
&- 1) - pub(ji, jj, jk)) / (e3uw_n(ji, jj, jk) * e3uw_b(ji, jj, jk)) * wumask(ji, jj, jk)
          zsh2v(ji, jj) = (p_avm(ji, jj + 1, jk) + p_avm(ji, jj, jk)) * (pvn(ji, jj, jk - 1) - pvn(ji, jj, jk)) * (pvb(ji, jj, jk &
&- 1) - pvb(ji, jj, jk)) / (e3vw_n(ji, jj, jk) * e3vw_b(ji, jj, jk)) * wvmask(ji, jj, jk)
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          p_sh2(ji, jj, jk) = 0.25 * ((zsh2u(ji - 1, jj) + zsh2u(ji, jj)) * (2. - umask(ji - 1, jj, jk) * umask(ji, jj, jk)) + &
&(zsh2v(ji, jj - 1) + zsh2v(ji, jj)) * (2. - vmask(ji, jj - 1, jk) * vmask(ji, jj, jk)))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
  END SUBROUTINE zdf_sh2
END MODULE zdfsh2