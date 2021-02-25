MODULE dynadv_cen2
  USE oce
  USE dom_oce
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE lib_mpp
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_adv_cen2
  CONTAINS
  SUBROUTINE dyn_adv_cen2(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zfu_t, zfu_f, zfu_uw, zfu
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zfv_t, zfv_f, zfv_vw, zfv, zfw
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dyn_adv_cen2 : 2nd order flux form momentum advection'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
    END IF
    IF (l_trddyn) THEN
      zfu_uw(:, :, :) = ua(:, :, :)
      zfv_vw(:, :, :) = va(:, :, :)
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      zfu(:, :, jk) = 0.25_wp * e2u(:, :) * e3u_n(:, :, jk) * un(:, :, jk)
      zfv(:, :, jk) = 0.25_wp * e1v(:, :) * e3v_n(:, :, jk) * vn(:, :, jk)
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zfu_t(ji + 1, jj, jk) = (zfu(ji, jj, jk) + zfu(ji + 1, jj, jk)) * (un(ji, jj, jk) + un(ji + 1, jj, jk))
          zfv_f(ji, jj, jk) = (zfv(ji, jj, jk) + zfv(ji + 1, jj, jk)) * (un(ji, jj, jk) + un(ji, jj + 1, jk))
          zfu_f(ji, jj, jk) = (zfu(ji, jj, jk) + zfu(ji, jj + 1, jk)) * (vn(ji, jj, jk) + vn(ji + 1, jj, jk))
          zfv_t(ji, jj + 1, jk) = (zfv(ji, jj, jk) + zfv(ji, jj + 1, jk)) * (vn(ji, jj, jk) + vn(ji, jj + 1, jk))
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zfu_t(ji + 1, jj, jk) - zfu_t(ji, jj, jk) + zfv_f(ji, jj, jk) - zfv_f(ji, jj - 1, &
&jk)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (zfu_f(ji, jj, jk) - zfu_f(ji - 1, jj, jk) + zfv_t(ji, jj + 1, jk) - zfv_t(ji, jj, &
&jk)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (l_trddyn) THEN
      zfu_uw(:, :, :) = ua(:, :, :) - zfu_uw(:, :, :)
      zfv_vw(:, :, :) = va(:, :, :) - zfv_vw(:, :, :)
      CALL trd_dyn(zfu_uw, zfv_vw, jpdyn_keg, kt)
      zfu_t(:, :, :) = ua(:, :, :)
      zfv_t(:, :, :) = va(:, :, :)
    END IF
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zfu_uw(ji, jj, jpk) = 0._wp
        zfv_vw(ji, jj, jpk) = 0._wp
        zfu_uw(ji, jj, 1) = 0._wp
        zfv_vw(ji, jj, 1) = 0._wp
      END DO
    END DO
    IF (ln_linssh) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfu_uw(ji, jj, 1) = 0.5_wp * (e1e2t(ji, jj) * wn(ji, jj, 1) + e1e2t(ji + 1, jj) * wn(ji + 1, jj, 1)) * un(ji, jj, 1)
          zfv_vw(ji, jj, 1) = 0.5_wp * (e1e2t(ji, jj) * wn(ji, jj, 1) + e1e2t(ji, jj + 1) * wn(ji, jj + 1, 1)) * vn(ji, jj, 1)
        END DO
      END DO
    END IF
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      DO jj = 2, jpj
        DO ji = 2, jpi
          zfw(ji, jj, jk) = 0.25_wp * e1e2t(ji, jj) * wn(ji, jj, jk)
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfu_uw(ji, jj, jk) = (zfw(ji, jj, jk) + zfw(ji + 1, jj, jk)) * (un(ji, jj, jk) + un(ji, jj, jk - 1))
          zfv_vw(ji, jj, jk) = (zfw(ji, jj, jk) + zfw(ji, jj + 1, jk)) * (vn(ji, jj, jk) + vn(ji, jj, jk - 1))
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zfu_uw(ji, jj, jk) - zfu_uw(ji, jj, jk + 1)) * r1_e1e2u(ji, jj) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (zfv_vw(ji, jj, jk) - zfv_vw(ji, jj, jk + 1)) * r1_e1e2v(ji, jj) / e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (l_trddyn) THEN
      zfu_t(:, :, :) = ua(:, :, :) - zfu_t(:, :, :)
      zfv_t(:, :, :) = va(:, :, :) - zfv_t(:, :, :)
      CALL trd_dyn(zfu_t, zfv_t, jpdyn_zad, kt)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' cen2 adv - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
  END SUBROUTINE dyn_adv_cen2
END MODULE dynadv_cen2