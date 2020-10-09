MODULE ocealb
  USE phycst
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: oce_alb
  CONTAINS
  SUBROUTINE oce_alb(palb_os, palb_cs)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: palb_os
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT) :: palb_cs
    REAL(KIND = wp) :: zcoef
    REAL(KIND = wp) :: rmue = 0.40
    !$ACC KERNELS
    zcoef = 0.05 / (1.1 * rmue ** 1.4 + 0.15)
    palb_cs(:, :) = zcoef
    palb_os(:, :) = 0.06
    !$ACC END KERNELS
  END SUBROUTINE oce_alb
END MODULE ocealb