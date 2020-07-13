MODULE lib_fortran
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: glob_sum
  PUBLIC :: glob_sum_full
  PUBLIC :: DDPDD
  PUBLIC :: glob_min, glob_max
  PUBLIC :: SIGN
  INTERFACE glob_sum
    MODULE PROCEDURE glob_sum_1d, glob_sum_2d, glob_sum_3d, glob_sum_2d_a, glob_sum_3d_a
  END INTERFACE
  INTERFACE glob_sum_full
    MODULE PROCEDURE glob_sum_full_2d, glob_sum_full_3d
  END INTERFACE
  INTERFACE glob_min
    MODULE PROCEDURE glob_min_2d, glob_min_3d, glob_min_2d_a, glob_min_3d_a
  END INTERFACE
  INTERFACE glob_max
    MODULE PROCEDURE glob_max_2d, glob_max_3d, glob_max_2d_a, glob_max_3d_a
  END INTERFACE
  INTERFACE SIGN
    MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D, SIGN_ARRAY_1D_A, SIGN_ARRAY_2D_A, SIGN_ARRAY_3D_A, SIGN_ARRAY_1D_B, SIGN_ARRAY_2D_B, SIGN_ARRAY_3D_B
  END INTERFACE
  CONTAINS
  FUNCTION glob_sum_1d(ptab, kdim)
    INTEGER, INTENT(IN) :: kdim
    REAL(KIND = wp), INTENT(IN), DIMENSION(kdim) :: ptab
    REAL(KIND = wp) :: glob_sum_1d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji

    ctmp = CMPLX(0.E0, 0.E0, wp)
    !$ACC LOOP SEQ
    DO ji = 1, kdim
      ztmp = ptab(ji)
      CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
    END DO
    !$ACC END LOOP
    IF (lk_mpp) CALL mpp_sum(ctmp)
    glob_sum_1d = REAL(ctmp, wp)
  END FUNCTION glob_sum_1d

  FUNCTION glob_sum_c1d(ptab, kdim)
    INTEGER, INTENT(IN) :: kdim
    COMPLEX(KIND = wp), INTENT(IN), DIMENSION(kdim) :: ptab
    REAL(KIND = wp) :: glob_sum_c1d
    COMPLEX(KIND = wp) :: ctmp
    INTEGER :: ji

    ctmp = CMPLX(0.E0, 0.E0, wp)
    !$ACC LOOP SEQ
    DO ji = 1, kdim
      CALL DDPDD(ptab(ji), ctmp)
    END DO
    !$ACC END LOOP
    IF (lk_mpp) CALL mpp_sum(ctmp)
    glob_sum_c1d = REAL(ctmp, wp)
  END FUNCTION glob_sum_c1d

  FUNCTION glob_sum_2d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab
    REAL(KIND = wp) :: glob_sum_2d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj
    COMPLEX(KIND = wp) :: hsum(jpj)
    
    !$ACC KERNELS
    DO jj = 1, jpj
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO ji = 1, jpi
        ztmp = ptab(ji, jj) * tmask_i(ji, jj)
        CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
      END DO
      hsum(jj) = ctmp
    END DO
    !$ACC END KERNELS
    glob_sum_2d =  glob_sum_c1d(hsum, jpj)
  END FUNCTION glob_sum_2d
  FUNCTION glob_sum_3d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab
    REAL(KIND = wp) :: glob_sum_3d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ijpk
    COMPLEX(KIND = wp), allocatable :: hsum(:)

    ijpk = SIZE(ptab, 3)
    allocate(hsum(ijpk))

    !$ACC KERNELS
    DO jk = 1, ijpk
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ztmp = ptab(ji, jj, jk) * tmask_i(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
      hsum(jk) = ctmp
    END DO
    !$ACC END KERNELS
    glob_sum_3d = glob_sum_c1d(hsum, ijpk)
    deallocate(hsum)
  END FUNCTION glob_sum_3d
  FUNCTION glob_sum_2d_a(ptab1, ptab2)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab1, ptab2
    REAL(KIND = wp) :: glob_sum_2d_a
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj
    COMPLEX(KIND = wp) :: hsum(jpj)

    !$ACC KERNELS
    DO jj = 1, jpj
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO ji = 1, jpi
        ztmp = ptab1(ji, jj) * tmask_i(ji, jj)
        CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        ztmp = ptab2(ji, jj) * tmask_i(ji, jj)
        CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
      END DO
      hsum(jj) = ctmp
    END DO
    !$ACC END KERNELS
    glob_sum_2d_a = glob_sum_c1d(hsum, jpj)
  END FUNCTION glob_sum_2d_a
  FUNCTION glob_sum_3d_a(ptab1, ptab2)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab1, ptab2
    REAL(KIND = wp) :: glob_sum_3d_a
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ijpk
    COMPLEX(KIND = wp), allocatable :: hsum(:)

    ijpk = SIZE(ptab1, 3)
    allocate(hsum(ijpk))
    ztmp = 0.E0
    !$ACC KERNELS
    DO jk = 1, ijpk
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ztmp = ptab1(ji, jj, jk) * tmask_i(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
          ztmp = ptab2(ji, jj, jk) * tmask_i(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
      hsum(jk) = ctmp
    END DO
    !$ACC END KERNELS

    glob_sum_3d_a = glob_sum_c1d(hsum, ijpk)
    deallocate(hsum)
  END FUNCTION glob_sum_3d_a
  FUNCTION glob_sum_full_2d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab
    REAL(KIND = wp) :: glob_sum_full_2d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj
    COMPLEX(KIND = wp) :: hsum(jpj)

    !$ACC KERNELS
    DO jj = 1, jpj
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO ji = 1, jpi
        ztmp = ptab(ji, jj) * tmask_h(ji, jj)
        CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
      END DO
      hsum(jj) = ctmp
    END DO
    !$ACC END KERNELS
    glob_sum_full_2d = glob_sum_c1d(hsum, jpj)
  END FUNCTION glob_sum_full_2d
  FUNCTION glob_sum_full_3d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab
    REAL(KIND = wp) :: glob_sum_full_3d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ijpk
    COMPLEX(KIND = wp), allocatable :: hsum(:)

    ijpk = SIZE(ptab, 3)
    allocate(hsum(ijpk))
    !$ACC KERNELS
    DO jk = 1, ijpk
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ztmp = ptab(ji, jj, jk) * tmask_h(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
      hsum(jk) = ctmp
    END DO
    !$ACC END KERNELS
    glob_sum_full_3d = glob_sum_c1d(hsum, ijpk)
    deallocate(hsum)
  END FUNCTION glob_sum_full_3d
  FUNCTION glob_min_2d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab
    REAL(KIND = wp) :: glob_min_2d
    !$ACC KERNELS ! CDe added
    glob_min_2d = MINVAL(ptab(:, :) * tmask_i(:, :))
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_min(glob_min_2d)
  END FUNCTION glob_min_2d
  FUNCTION glob_min_3d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab
    REAL(KIND = wp) :: glob_min_3d
    INTEGER :: jk
    INTEGER :: ijpk
    ijpk = SIZE(ptab, 3)
    !$ACC KERNELS ! CDe added
    glob_min_3d = MINVAL(ptab(:, :, 1) * tmask_i(:, :))
    DO jk = 2, ijpk
      glob_min_3d = MIN(glob_min_3d, MINVAL(ptab(:, :, jk) * tmask_i(:, :)))
    END DO
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_min(glob_min_3d)
  END FUNCTION glob_min_3d
  FUNCTION glob_min_2d_a(ptab1, ptab2)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab1, ptab2
    REAL(KIND = wp), DIMENSION(2) :: glob_min_2d_a
    !$ACC KERNELS ! CDe added
    glob_min_2d_a(1) = MINVAL(ptab1(:, :) * tmask_i(:, :))
    glob_min_2d_a(2) = MINVAL(ptab2(:, :) * tmask_i(:, :))
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_min(glob_min_2d_a, 2)
  END FUNCTION glob_min_2d_a
  FUNCTION glob_min_3d_a(ptab1, ptab2)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab1, ptab2
    REAL(KIND = wp), DIMENSION(2) :: glob_min_3d_a
    INTEGER :: jk
    INTEGER :: ijpk
    ijpk = SIZE(ptab1, 3)
    !$ACC KERNELS ! CDe added
    glob_min_3d_a(1) = MINVAL(ptab1(:, :, 1) * tmask_i(:, :))
    glob_min_3d_a(2) = MINVAL(ptab2(:, :, 1) * tmask_i(:, :))
    DO jk = 2, ijpk
      glob_min_3d_a(1) = MIN(glob_min_3d_a(1), MINVAL(ptab1(:, :, jk) * tmask_i(:, :)))
      glob_min_3d_a(2) = MIN(glob_min_3d_a(2), MINVAL(ptab2(:, :, jk) * tmask_i(:, :)))
    END DO
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_min(glob_min_3d_a, 2)
  END FUNCTION glob_min_3d_a
  FUNCTION glob_max_2d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab
    REAL(KIND = wp) :: glob_max_2d
    !$ACC KERNELS ! CDe added
    glob_max_2d = MAXVAL(ptab(:, :) * tmask_i(:, :))
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_max(glob_max_2d)
  END FUNCTION glob_max_2d
  FUNCTION glob_max_3d(ptab)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab
    REAL(KIND = wp) :: glob_max_3d
    INTEGER :: jk
    INTEGER :: ijpk
    ijpk = SIZE(ptab, 3)
    !$ACC KERNELS ! CDe added
    glob_max_3d = MAXVAL(ptab(:, :, 1) * tmask_i(:, :))
    DO jk = 2, ijpk
      glob_max_3d = MAX(glob_max_3d, MAXVAL(ptab(:, :, jk) * tmask_i(:, :)))
    END DO
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_max(glob_max_3d)
  END FUNCTION glob_max_3d
  FUNCTION glob_max_2d_a(ptab1, ptab2)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: ptab1, ptab2
    REAL(KIND = wp), DIMENSION(2) :: glob_max_2d_a
    !$ACC KERNELS ! CDe added
    glob_max_2d_a(1) = MAXVAL(ptab1(:, :) * tmask_i(:, :))
    glob_max_2d_a(2) = MAXVAL(ptab2(:, :) * tmask_i(:, :))
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_max(glob_max_2d_a, 2)
  END FUNCTION glob_max_2d_a
  FUNCTION glob_max_3d_a(ptab1, ptab2)
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: ptab1, ptab2
    REAL(KIND = wp), DIMENSION(2) :: glob_max_3d_a
    INTEGER :: jk
    INTEGER :: ijpk
    ijpk = SIZE(ptab1, 3)
    !$ACC KERNELS ! CDE added
    glob_max_3d_a(1) = MAXVAL(ptab1(:, :, 1) * tmask_i(:, :))
    glob_max_3d_a(2) = MAXVAL(ptab2(:, :, 1) * tmask_i(:, :))
    DO jk = 2, ijpk
      glob_max_3d_a(1) = MAX(glob_max_3d_a(1), MAXVAL(ptab1(:, :, jk) * tmask_i(:, :)))
      glob_max_3d_a(2) = MAX(glob_max_3d_a(2), MAXVAL(ptab2(:, :, jk) * tmask_i(:, :)))
    END DO
    !$ACC END KERNELS
    IF (lk_mpp) CALL mpp_max(glob_max_3d_a, 2)
  END FUNCTION glob_max_3d_a
  SUBROUTINE DDPDD(ydda, yddb)
    COMPLEX(KIND = wp), INTENT(IN   ) :: ydda
    COMPLEX(KIND = wp), INTENT(INOUT) :: yddb
    REAL(KIND = wp) :: zerr, zt1, zt2
    zt1 = REAL(ydda) + REAL(yddb)
    zerr = zt1 - REAL(ydda)
    zt2 = ((REAL(yddb) - zerr) + (REAL(ydda) - (zt1 - zerr))) + AIMAG(ydda) + AIMAG(yddb)
    yddb = CMPLX(zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), wp)
  END SUBROUTINE DDPDD
  FUNCTION SIGN_SCALAR(pa, pb)
    REAL(KIND = wp) :: pa, pb
    !$ACC routine
    REAL(KIND = wp) :: SIGN_SCALAR
    IF (pb >= 0.E0) THEN
      sign_scalar = abs(pa)
    ELSE
      sign_scalar = - abs(pa)
    END IF
  END FUNCTION SIGN_SCALAR
  FUNCTION SIGN_ARRAY_1D(pa, pb)
    REAL(KIND = wp) :: pa, pb(:)
    REAL(KIND = wp) :: SIGN_ARRAY_1D(SIZE(pb, 1))
    !$ACC routine
    WHERE (pb >= 0.E0)
      sign_array_1d = abs(pa)
    ELSEWHERE
      sign_array_1d = - abs(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_1D
  FUNCTION SIGN_ARRAY_2D(pa, pb)
    REAL(KIND = wp) :: pa, pb(:, :)
    REAL(KIND = wp) :: SIGN_ARRAY_2D(SIZE(pb, 1), SIZE(pb, 2))
    !$ACC routine
    WHERE (pb >= 0.E0)
      sign_array_2d = abs(pa)
    ELSEWHERE
      sign_array_2d = - abs(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_2D
  FUNCTION SIGN_ARRAY_3D(pa, pb)
    REAL(KIND = wp) :: pa, pb(:, :, :)
    REAL(KIND = wp) :: SIGN_ARRAY_3D(SIZE(pb, 1), SIZE(pb, 2), SIZE(pb, 3))
    !$ACC routine
    WHERE (pb >= 0.E0)
      sign_array_3d = abs(pa)
    ELSEWHERE
      sign_array_3d = - abs(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_3D
  FUNCTION SIGN_ARRAY_1D_A(pa, pb)
    REAL(KIND = wp) :: pa(:), pb(:)
    REAL(KIND = wp) :: SIGN_ARRAY_1D_A(SIZE(pb, 1))
    WHERE (pb >= 0.E0)
      sign_array_1d_a = abs(pa)
    ELSEWHERE
      sign_array_1d_a = - abs(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_1D_A
  FUNCTION SIGN_ARRAY_2D_A(pa, pb)
    REAL(KIND = wp) :: pa(:, :), pb(:, :)
    REAL(KIND = wp) :: SIGN_ARRAY_2D_A(SIZE(pb, 1), SIZE(pb, 2))
    !$ACC routine ! CDe added
    WHERE (pb >= 0.E0)
      sign_array_2d_a = abs(pa)
    ELSEWHERE
      sign_array_2d_a = - abs(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_2D_A
  FUNCTION SIGN_ARRAY_3D_A(pa, pb)
    REAL(KIND = wp) :: pa(:, :, :), pb(:, :, :)
    REAL(KIND = wp) :: SIGN_ARRAY_3D_A(SIZE(pb, 1), SIZE(pb, 2), SIZE(pb, 3))
    WHERE (pb >= 0.E0)
      sign_array_3d_a = abs(pa)
    ELSEWHERE
      sign_array_3d_a = - abs(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_3D_A
  FUNCTION SIGN_ARRAY_1D_B(pa, pb)
    REAL(KIND = wp) :: pa(:), pb
    REAL(KIND = wp) :: SIGN_ARRAY_1D_B(SIZE(pa, 1))
    IF (pb >= 0.E0) THEN
      sign_array_1d_b = abs(pa)
    ELSE
      sign_array_1d_b = - abs(pa)
    END IF
  END FUNCTION SIGN_ARRAY_1D_B
  FUNCTION SIGN_ARRAY_2D_B(pa, pb)
    REAL(KIND = wp) :: pa(:, :), pb
    REAL(KIND = wp) :: SIGN_ARRAY_2D_B(SIZE(pa, 1), SIZE(pa, 2))
    IF (pb >= 0.E0) THEN
      sign_array_2d_b = abs(pa)
    ELSE
      sign_array_2d_b = - abs(pa)
    END IF
  END FUNCTION SIGN_ARRAY_2D_B
  FUNCTION SIGN_ARRAY_3D_B(pa, pb)
    REAL(KIND = wp) :: pa(:, :, :), pb
    REAL(KIND = wp) :: SIGN_ARRAY_3D_B(SIZE(pa, 1), SIZE(pa, 2), SIZE(pa, 3))
    IF (pb >= 0.E0) THEN
      sign_array_3d_b = abs(pa)
    ELSE
      sign_array_3d_b = - abs(pa)
    END IF
  END FUNCTION SIGN_ARRAY_3D_B
END MODULE lib_fortran
