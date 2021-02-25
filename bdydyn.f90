MODULE bdydyn
  USE oce
  USE dom_oce
  USE bdy_oce
  USE bdydyn2d
  USE bdydyn3d
  USE lbclnk
  USE in_out_manager
  USE domvvl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_dyn
  CONTAINS
  SUBROUTINE bdy_dyn(kt, dyn3d_only)
    INTEGER, INTENT(IN) :: kt
    LOGICAL, INTENT(IN), OPTIONAL :: dyn3d_only
    INTEGER :: jk, ii, ij, ib_bdy, ib, igrd
    LOGICAL :: ll_dyn2d, ll_dyn3d, ll_orlanski
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: pua2d, pva2d
    ll_dyn2d = .TRUE.
    ll_dyn3d = .TRUE.
    IF (PRESENT(dyn3d_only)) THEN
      IF (dyn3d_only) ll_dyn2d = .FALSE.
    END IF
    ll_orlanski = .FALSE.
    DO ib_bdy = 1, nb_bdy
      IF (cn_dyn2d(ib_bdy) == 'orlanski' .OR. cn_dyn2d(ib_bdy) == 'orlanski_npo' .OR. cn_dyn3d(ib_bdy) == 'orlanski' .OR. &
&cn_dyn3d(ib_bdy) == 'orlanski_npo') ll_orlanski = .TRUE.
    END DO
    pua2d(:, :) = 0._wp
    pva2d(:, :) = 0._wp
    ! !$OMP parallel default(shared), private(jk) ! CDe race condition??
    ! !$OMP do schedule(static)
    DO jk = 1, jpkm1
      pua2d(:, :) = pua2d(:, :) + e3u_a(:, :, jk) * ua(:, :, jk) * umask(:, :, jk)
      pva2d(:, :) = pva2d(:, :) + e3v_a(:, :, jk) * va(:, :, jk) * vmask(:, :, jk)
    END DO
    ! !$OMP end do
    ! !$OMP end parallel
    pua2d(:, :) = pua2d(:, :) * r1_hu_a(:, :)
    pva2d(:, :) = pva2d(:, :) * r1_hv_a(:, :)
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      ua(:, :, jk) = (ua(:, :, jk) - pua2d(:, :)) * umask(:, :, jk)
      va(:, :, jk) = (va(:, :, jk) - pva2d(:, :)) * vmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (ll_orlanski) THEN
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ub(:, :, jk) = (ub(:, :, jk) - ub_b(:, :)) * umask(:, :, jk)
        vb(:, :, jk) = (vb(:, :, jk) - vb_b(:, :)) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
    IF (ll_dyn2d) CALL bdy_dyn2d(kt, pua2d, pva2d, ub_b, vb_b, r1_hu_a(:, :), r1_hv_a(:, :), ssha)
    IF (ll_dyn3d) CALL bdy_dyn3d(kt)
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      ua(:, :, jk) = (ua(:, :, jk) + pua2d(:, :)) * umask(:, :, jk)
      va(:, :, jk) = (va(:, :, jk) + pva2d(:, :)) * vmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (ll_orlanski) THEN
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ub(:, :, jk) = (ub(:, :, jk) + ub_b(:, :)) * umask(:, :, jk)
        vb(:, :, jk) = (vb(:, :, jk) + vb_b(:, :)) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
  END SUBROUTINE bdy_dyn
END MODULE bdydyn
