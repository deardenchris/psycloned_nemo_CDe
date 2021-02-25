MODULE obs_inter_z1d
  USE par_kind, ONLY: wp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_int_z1d, obs_int_z1d_spl
  CONTAINS
  SUBROUTINE obs_int_z1d(kpk, kkco, k1dint, kdep, pobsdep, pobsk, pobs2k, pobs, pdep, pobsmask)
    INTEGER, INTENT(IN) :: kpk
    INTEGER, INTENT(IN) :: k1dint
    INTEGER, INTENT(IN) :: kdep
    INTEGER, INTENT(IN), DIMENSION(kdep) :: kkco
    REAL(KIND = wp), INTENT(IN), DIMENSION(kdep) :: pobsdep
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpk) :: pobsk, pobs2k, pdep, pobsmask
    REAL(KIND = wp), INTENT(OUT), DIMENSION(kdep) :: pobs
    REAL(KIND = wp) :: z1dm
    REAL(KIND = wp) :: z1dp
    REAL(KIND = wp) :: zsum
    REAL(KIND = wp) :: zsum2
    INTEGER :: jdep
    DO jdep = 1, kdep
      z1dm = (pdep(kkco(jdep)) - pobsdep(jdep))
      z1dp = (pobsdep(jdep) - pdep(kkco(jdep) - 1))
      IF (pobsmask(kkco(jdep)) == 0.0_wp) z1dp = 0.0_wp
      zsum = z1dm + z1dp
      IF (k1dint == 0) THEN
        pobs(jdep) = (z1dm * pobsk(kkco(jdep) - 1) + z1dp * pobsk(kkco(jdep))) / zsum
      ELSE IF (k1dint == 1) THEN
        zsum2 = zsum * zsum
        pobs(jdep) = (z1dm * pobsk(kkco(jdep) - 1) + z1dp * pobsk(kkco(jdep)) + (z1dm * (z1dm * z1dm - zsum2) * pobs2k(kkco(jdep) &
&- 1) + z1dp * (z1dp * z1dp - zsum2) * pobs2k(kkco(jdep))) / 6.0_wp) / zsum
      END IF
    END DO
  END SUBROUTINE obs_int_z1d
  SUBROUTINE obs_int_z1d_spl(kpk, pobsk, pobs2k, pdep, pobsmask)
    INTEGER, INTENT(IN) :: kpk
    REAL(KIND = wp), INTENT(IN), DIMENSION(kpk) :: pobsk, pdep, pobsmask
    REAL(KIND = wp), INTENT(OUT), DIMENSION(kpk) :: pobs2k
    INTEGER :: jk
    REAL(KIND = wp) :: za
    REAL(KIND = wp) :: zb
    REAL(KIND = wp) :: zc
    REAL(KIND = wp) :: zpa
    REAL(KIND = wp) :: zkm
    REAL(KIND = wp) :: zkp
    REAL(KIND = wp) :: zk
    REAL(KIND = wp), DIMENSION(kpk - 1) :: zs, zp, zu, zv
    zs(1) = 0.0_wp
    zp(1) = 0.0_wp
    zv(1) = - 0.5_wp
    DO jk = 2, kpk - 1
      zs(jk) = (pdep(jk) - pdep(jk - 1)) / (pdep(jk + 1) - pdep(jk - 1))
      zp(jk) = zs(jk) * zv(jk - 1) + 2.0_wp
      zv(jk) = (zs(jk) - 1.0_wp) / zp(jk)
    END DO
    zu(1) = 0.0_wp
    DO jk = 2, kpk - 1
      za = pdep(jk + 1) - pdep(jk - 1)
      zb = pdep(jk + 1) - pdep(jk)
      zc = pdep(jk) - pdep(jk - 1)
      zpa = 6.0_wp / (zp(jk) * za)
      zkm = zpa / zc
      zkp = zpa / zb
      zk = - (zkm + zkp)
      zu(jk) = pobsk(jk + 1) * zkp + pobsk(jk) * zk + pobsk(jk - 1) * zkm + zu(jk - 1) * (- zs(jk) / zp(jk))
    END DO
    pobs2k(kpk) = 0.0_wp
    DO jk = kpk - 1, 1, - 1
      pobs2k(jk) = zv(jk) * pobs2k(jk + 1) + zu(jk)
      IF (pobsmask(jk + 1) == 0.0_wp) pobs2k(jk) = 0.0_wp
    END DO
  END SUBROUTINE obs_int_z1d_spl
END MODULE obs_inter_z1d