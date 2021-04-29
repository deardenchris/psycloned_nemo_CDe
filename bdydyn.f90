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
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    LOGICAL, INTENT(IN), OPTIONAL :: dyn3d_only
    INTEGER :: jk, ii, ij, ib_bdy, ib, igrd
    LOGICAL :: ll_dyn2d, ll_dyn3d, ll_orlanski
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: pua2d, pva2d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('bdy_dyn', 'r0', 0, 0)
    ll_dyn2d = .TRUE.
    ll_dyn3d = .TRUE.
    IF (PRESENT(dyn3d_only)) THEN
      IF (dyn3d_only) ll_dyn2d = .FALSE.
    END IF
    ll_orlanski = .FALSE.
    DO ib_bdy = 1, nb_bdy
      IF (cn_dyn2d(ib_bdy) == 'orlanski' .OR. cn_dyn2d(ib_bdy) == 'orlanski_npo' .OR. cn_dyn3d(ib_bdy) == 'orlanski' .OR. cn_dyn3d(ib_bdy) == 'orlanski_npo') ll_orlanski = .TRUE.
    END DO
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    pua2d(:, :) = 0._wp
    pva2d(:, :) = 0._wp
    DO jk = 1, jpkm1
      pua2d(:, :) = pua2d(:, :) + e3u_a(:, :, jk) * ua(:, :, jk) * umask(:, :, jk)
      pva2d(:, :) = pva2d(:, :) + e3v_a(:, :, jk) * va(:, :, jk) * vmask(:, :, jk)
    END DO
    pua2d(:, :) = pua2d(:, :) * r1_hu_a(:, :)
    pva2d(:, :) = pva2d(:, :) * r1_hv_a(:, :)
    DO jk = 1, jpkm1
      ua(:, :, jk) = (ua(:, :, jk) - pua2d(:, :)) * umask(:, :, jk)
      va(:, :, jk) = (va(:, :, jk) - pva2d(:, :)) * vmask(:, :, jk)
    END DO
    !$ACC END KERNELS
    IF (ll_orlanski) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        ub(:, :, jk) = (ub(:, :, jk) - ub_b(:, :)) * umask(:, :, jk)
        vb(:, :, jk) = (vb(:, :, jk) - vb_b(:, :)) * vmask(:, :, jk)
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data1 % PreStart('bdy_dyn', 'r1', 0, 0)
    IF (ll_dyn2d) CALL bdy_dyn2d(kt, pua2d, pva2d, ub_b, vb_b, r1_hu_a(:, :), r1_hv_a(:, :), ssha)
    IF (ll_dyn3d) CALL bdy_dyn3d(kt)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    DO jk = 1, jpkm1
      ua(:, :, jk) = (ua(:, :, jk) + pua2d(:, :)) * umask(:, :, jk)
      va(:, :, jk) = (va(:, :, jk) + pva2d(:, :)) * vmask(:, :, jk)
    END DO
    !$ACC END KERNELS
    IF (ll_orlanski) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        ub(:, :, jk) = (ub(:, :, jk) + ub_b(:, :)) * umask(:, :, jk)
        vb(:, :, jk) = (vb(:, :, jk) + vb_b(:, :)) * vmask(:, :, jk)
      END DO
      !$ACC END KERNELS
    END IF
  END SUBROUTINE bdy_dyn
END MODULE bdydyn