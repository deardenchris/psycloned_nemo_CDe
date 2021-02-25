MODULE depth_e3
  USE oce
  USE dom_oce
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  INTERFACE depth_to_e3
    MODULE PROCEDURE depth_to_e3_1d, depth_to_e3_3d
  END INTERFACE
  INTERFACE e3_to_depth
    MODULE PROCEDURE e3_to_depth_1d, e3_to_depth_3d
  END INTERFACE
  PUBLIC :: depth_to_e3
  PUBLIC :: e3_to_depth
  CONTAINS
  SUBROUTINE depth_to_e3_1d(pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d)
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pdept_1d, pdepw_1d
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pe3t_1d, pe3w_1d
    INTEGER :: jk
    pe3w_1d(1) = 2._wp * (pdept_1d(1) - pdepw_1d(1))
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      pe3w_1d(jk + 1) = pdept_1d(jk + 1) - pdept_1d(jk)
      pe3t_1d(jk) = pdepw_1d(jk + 1) - pdepw_1d(jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    pe3t_1d(jpk) = 2._wp * (pdept_1d(jpk) - pdepw_1d(jpk))
  END SUBROUTINE depth_to_e3_1d
  SUBROUTINE depth_to_e3_3d(pdept_3d, pdepw_3d, pe3t_3d, pe3w_3d)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pdept_3d, pdepw_3d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3t_3d, pe3w_3d
    INTEGER :: jk
    pe3w_3d(:, :, 1) = 2._wp * (pdept_3d(:, :, 1) - pdepw_3d(:, :, 1))
    DO jk = 1, jpkm1
      pe3w_3d(:, :, jk + 1) = pdept_3d(:, :, jk + 1) - pdept_3d(:, :, jk)
      pe3t_3d(:, :, jk) = pdepw_3d(:, :, jk + 1) - pdepw_3d(:, :, jk)
    END DO
    pe3t_3d(:, :, jpk) = 2._wp * (pdept_3d(:, :, jpk) - pdepw_3d(:, :, jpk))
  END SUBROUTINE depth_to_e3_3d
  SUBROUTINE e3_to_depth_1d(pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d)
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pe3t_1d, pe3w_1d
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pdept_1d, pdepw_1d
    INTEGER :: jk
    pdepw_1d(1) = 0.0_wp
    pdept_1d(1) = 0.5_wp * pe3w_1d(1)
    DO jk = 2, jpk
      pdepw_1d(jk) = pdepw_1d(jk - 1) + pe3t_1d(jk - 1)
      pdept_1d(jk) = pdept_1d(jk - 1) + pe3w_1d(jk)
    END DO
  END SUBROUTINE e3_to_depth_1d
  SUBROUTINE e3_to_depth_3d(pe3t_3d, pe3w_3d, pdept_3d, pdepw_3d)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pe3t_3d, pe3w_3d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pdept_3d, pdepw_3d
    INTEGER :: jk
    pdepw_3d(:, :, 1) = 0.0_wp
    pdept_3d(:, :, 1) = 0.5_wp * pe3w_3d(:, :, 1)
    DO jk = 2, jpk
      pdepw_3d(:, :, jk) = pdepw_3d(:, :, jk - 1) + pe3t_3d(:, :, jk - 1)
      pdept_3d(:, :, jk) = pdept_3d(:, :, jk - 1) + pe3w_3d(:, :, jk)
    END DO
  END SUBROUTINE e3_to_depth_3d
END MODULE depth_e3