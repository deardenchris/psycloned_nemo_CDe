MODULE lib_fortran
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: glob_sum
  PUBLIC :: glob_sum_full
  PUBLIC :: local_sum
  PUBLIC :: sum3x3
  PUBLIC :: DDPDD
  PUBLIC :: glob_min, glob_max
  PUBLIC :: SIGN
  INTERFACE glob_sum
    MODULE PROCEDURE glob_sum_1d, glob_sum_2d, glob_sum_3d
  END INTERFACE
  INTERFACE glob_sum_full
    MODULE PROCEDURE glob_sum_full_2d, glob_sum_full_3d
  END INTERFACE
  INTERFACE local_sum
    MODULE PROCEDURE local_sum_2d, local_sum_3d
  END INTERFACE
  INTERFACE sum3x3
    MODULE PROCEDURE sum3x3_2d, sum3x3_3d
  END INTERFACE
  INTERFACE glob_min
    MODULE PROCEDURE glob_min_2d, glob_min_3d
  END INTERFACE
  INTERFACE glob_max
    MODULE PROCEDURE glob_max_2d, glob_max_3d
  END INTERFACE
  INTERFACE SIGN
    MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D, SIGN_ARRAY_1D_A, SIGN_ARRAY_2D_A, SIGN_ARRAY_3D_A, &
&SIGN_ARRAY_1D_B, SIGN_ARRAY_2D_B, SIGN_ARRAY_3D_B
  END INTERFACE
  CONTAINS
  FUNCTION glob_sum_c1d(ptab, kdim, ldcom, cdname)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(IN) :: kdim
    COMPLEX(KIND = wp), INTENT(IN), DIMENSION(kdim) :: ptab
    LOGICAL, INTENT(IN) :: ldcom
    REAL(KIND = wp) :: glob_sum_c1d
    COMPLEX(KIND = wp) :: ctmp
    INTEGER :: ji
    ctmp = CMPLX(0.E0, 0.E0, wp)
    DO ji = 1, kdim
      CALL DDPDD(ptab(ji), ctmp)
    END DO
    IF (ldcom) CALL mpp_sum(cdname, ctmp)
    glob_sum_c1d = REAL(ctmp, wp)
  END FUNCTION glob_sum_c1d
  FUNCTION glob_sum_1d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:)
    REAL(KIND = wp) :: glob_sum_1d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ipi, ipj, ipk
    ipi = SIZE(ptab, 1)
    ipj = 1
    ipk = 1
    ctmp = CMPLX(0.E0, 0.E0, wp)
    DO jk = 1, ipk
      DO jj = 1, ipj
        DO ji = 1, ipi
          ztmp = ptab(ji) * 1.
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
    END DO
    CALL mpp_sum(cdname, ctmp)
    glob_sum_1d = REAL(ctmp, wp)
  END FUNCTION glob_sum_1d
  FUNCTION glob_sum_2d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp) :: glob_sum_2d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ipi, ipj, ipk
    COMPLEX(KIND = wp), ALLOCATABLE :: hsum(:)
    ipi = SIZE(ptab, 1)
    ipj = SIZE(ptab, 2)
    ipk = 1
    ALLOCATE(hsum(ipj))
    DO jk = 1, ipk
      DO jj = 1, ipj
        ctmp = CMPLX(0.E0, 0.E0, wp)
        DO ji = 1, ipi
          ztmp = ptab(ji, jj) * tmask_i(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
        hsum(jj) = ctmp
      END DO
    END DO
    glob_sum_2d = glob_sum_c1d(hsum, ipj, .TRUE. .AND. lk_mpp, cdname)
    DEALLOCATE(hsum)
  END FUNCTION glob_sum_2d
  FUNCTION glob_sum_full_2d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp) :: glob_sum_full_2d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ipi, ipj, ipk
    COMPLEX(KIND = wp), ALLOCATABLE :: hsum(:)
    ipi = SIZE(ptab, 1)
    ipj = SIZE(ptab, 2)
    ipk = 1
    ALLOCATE(hsum(ipj))
    DO jk = 1, ipk
      DO jj = 1, ipj
        ctmp = CMPLX(0.E0, 0.E0, wp)
        DO ji = 1, ipi
          ztmp = ptab(ji, jj) * tmask_h(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
        hsum(jj) = ctmp
      END DO
    END DO
    glob_sum_full_2d = glob_sum_c1d(hsum, ipj, .TRUE. .AND. lk_mpp, cdname)
    DEALLOCATE(hsum)
  END FUNCTION glob_sum_full_2d
  FUNCTION glob_sum_3d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp) :: glob_sum_3d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ipi, ipj, ipk
    COMPLEX(KIND = wp), ALLOCATABLE :: hsum(:)
    ipi = SIZE(ptab, 1)
    ipj = SIZE(ptab, 2)
    ipk = SIZE(ptab, 3)
    ALLOCATE(hsum(ipk))
    DO jk = 1, ipk
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO jj = 1, ipj
        DO ji = 1, ipi
          ztmp = ptab(ji, jj, jk) * tmask_i(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
      hsum(jk) = ctmp
    END DO
    glob_sum_3d = glob_sum_c1d(hsum, ipk, .TRUE. .AND. lk_mpp, cdname)
    DEALLOCATE(hsum)
  END FUNCTION glob_sum_3d
  FUNCTION glob_sum_full_3d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp) :: glob_sum_full_3d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ipi, ipj, ipk
    COMPLEX(KIND = wp), ALLOCATABLE :: hsum(:)
    ipi = SIZE(ptab, 1)
    ipj = SIZE(ptab, 2)
    ipk = SIZE(ptab, 3)
    ALLOCATE(hsum(ipk))
    DO jk = 1, ipk
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO jj = 1, ipj
        DO ji = 1, ipi
          ztmp = ptab(ji, jj, jk) * tmask_h(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
      hsum(jk) = ctmp
    END DO
    glob_sum_full_3d = glob_sum_c1d(hsum, ipk, .TRUE. .AND. lk_mpp, cdname)
    DEALLOCATE(hsum)
  END FUNCTION glob_sum_full_3d
  FUNCTION glob_min_2d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp) :: glob_min_2d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: jk
    INTEGER :: ipk
    ipk = 1
    ztmp = MINVAL(ptab(:, :) * tmask_i(:, :))
    CALL mpp_min(cdname, ztmp)
    glob_min_2d = ztmp
  END FUNCTION glob_min_2d
  FUNCTION glob_max_2d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    REAL(KIND = wp) :: glob_max_2d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: jk
    INTEGER :: ipk
    ipk = 1
    ztmp = MAXVAL(ptab(:, :) * tmask_i(:, :))
    CALL mpp_max(cdname, ztmp)
    glob_max_2d = ztmp
  END FUNCTION glob_max_2d
  FUNCTION glob_min_3d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp) :: glob_min_3d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: jk
    INTEGER :: ipk
    ipk = SIZE(ptab, 3)
    ztmp = MINVAL(ptab(:, :, 1) * tmask_i(:, :))
    DO jk = 2, ipk
      ztmp = MIN(ztmp, MINVAL(ptab(:, :, jk) * tmask_i(:, :)))
    END DO
    CALL mpp_min(cdname, ztmp)
    glob_min_3d = ztmp
  END FUNCTION glob_min_3d
  FUNCTION glob_max_3d(cdname, ptab)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    REAL(KIND = wp) :: glob_max_3d
    REAL(KIND = wp) :: FUNCTION_GLOB_OP
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: jk
    INTEGER :: ipk
    ipk = SIZE(ptab, 3)
    ztmp = MAXVAL(ptab(:, :, 1) * tmask_i(:, :))
    DO jk = 2, ipk
      ztmp = MAX(ztmp, MAXVAL(ptab(:, :, jk) * tmask_i(:, :)))
    END DO
    CALL mpp_max(cdname, ztmp)
    glob_max_3d = ztmp
  END FUNCTION glob_max_3d
  FUNCTION local_sum_2d(ptab)
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :)
    COMPLEX(KIND = wp) :: local_sum_2d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj
    INTEGER :: ipi, ipj
    COMPLEX(KIND = wp), ALLOCATABLE :: hsum(:)
    ipi = SIZE(ptab, 1)
    ipj = SIZE(ptab, 2)
    ALLOCATE(hsum(ipj))
    DO jj = 1, ipj
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO ji = 1, ipi
        ztmp = ptab(ji, jj) * tmask_i(ji, jj)
        CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
      END DO
      hsum(jj) = ctmp
    END DO
    local_sum_2d = glob_sum_c1d(hsum, ipj, .FALSE., 'NONE')
    DEALLOCATE(hsum)
  END FUNCTION local_sum_2d
  FUNCTION local_sum_3d(ptab)
    REAL(KIND = wp), INTENT(IN) :: ptab(:, :, :)
    COMPLEX(KIND = wp) :: local_sum_3d
    COMPLEX(KIND = wp) :: ctmp
    REAL(KIND = wp) :: ztmp
    INTEGER :: ji, jj, jk
    INTEGER :: ipi, ipj, ipk
    COMPLEX(KIND = wp), ALLOCATABLE :: hsum(:)
    ipi = SIZE(ptab, 1)
    ipj = SIZE(ptab, 2)
    ipk = SIZE(ptab, 3)
    ALLOCATE(hsum(ipk))
    ctmp = CMPLX(0.E0, 0.E0, wp)
    DO jk = 1, ipk
      ctmp = CMPLX(0.E0, 0.E0, wp)
      DO jj = 1, ipj
        DO ji = 1, ipi
          ztmp = ptab(ji, jj, jk) * tmask_i(ji, jj)
          CALL DDPDD(CMPLX(ztmp, 0.E0, wp), ctmp)
        END DO
      END DO
      hsum(jk) = ctmp
    END DO
    local_sum_3d = glob_sum_c1d(hsum, ipk, .FALSE., 'NONE')
    DEALLOCATE(hsum)
  END FUNCTION local_sum_3d
  SUBROUTINE sum3x3_2d(p2d)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: p2d
    INTEGER :: ji, ji2, jj, jj2
    IF (SIZE(p2d, 1) /= jpi) CALL ctl_stop('STOP', 'wrong call of sum3x3_2d, the first dimension is not equal to jpi')
    IF (SIZE(p2d, 2) /= jpj) CALL ctl_stop('STOP', 'wrong call of sum3x3_2d, the second dimension is not equal to jpj')
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (MOD(mig(ji), 3) == 1 .AND. MOD(mjg(jj), 3) == 1) THEN
          ji2 = MIN(mig(ji) + 2, jpiglo) - nimpp + 1
          jj2 = MIN(mjg(jj) + 2, jpjglo) - njmpp + 1
          IF (ji2 <= jpi .AND. jj2 <= jpj) THEN
            p2d(ji : ji2, jj : jj2) = SUM(p2d(ji : ji2, jj : jj2))
          END IF
        END IF
      END DO
    END DO
    CALL lbc_lnk('lib_fortran', p2d, 'T', 1.)
    IF (nbondi /= - 1) THEN
      IF (MOD(mig(1), 3) == 1) p2d(1, :) = p2d(2, :)
      IF (MOD(mig(1), 3) == 2) p2d(2, :) = p2d(1, :)
    END IF
    IF (nbondi /= 1) THEN
      IF (MOD(mig(jpi - 2), 3) == 1) p2d(jpi, :) = p2d(jpi - 1, :)
      IF (MOD(mig(jpi - 2), 3) == 0) p2d(jpi - 1, :) = p2d(jpi, :)
    END IF
    IF (nbondj /= - 1) THEN
      IF (MOD(mjg(1), 3) == 1) p2d(:, 1) = p2d(:, 2)
      IF (MOD(mjg(1), 3) == 2) p2d(:, 2) = p2d(:, 1)
    END IF
    IF (nbondj /= 1) THEN
      IF (MOD(mjg(jpj - 2), 3) == 1) p2d(:, jpj) = p2d(:, jpj - 1)
      IF (MOD(mjg(jpj - 2), 3) == 0) p2d(:, jpj - 1) = p2d(:, jpj)
    END IF
    CALL lbc_lnk('lib_fortran', p2d, 'T', 1.)
  END SUBROUTINE sum3x3_2d
  SUBROUTINE sum3x3_3d(p3d)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p3d
    INTEGER :: ji, ji2, jj, jj2, jn
    INTEGER :: ipn
    IF (SIZE(p3d, 1) /= jpi) CALL ctl_stop('STOP', 'wrong call of sum3x3_3d, the first dimension is not equal to jpi')
    IF (SIZE(p3d, 2) /= jpj) CALL ctl_stop('STOP', 'wrong call of sum3x3_3d, the second dimension is not equal to jpj')
    ipn = SIZE(p3d, 3)
    DO jn = 1, ipn
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (MOD(mig(ji), 3) == 1 .AND. MOD(mjg(jj), 3) == 1) THEN
            ji2 = MIN(mig(ji) + 2, jpiglo) - nimpp + 1
            jj2 = MIN(mjg(jj) + 2, jpjglo) - njmpp + 1
            IF (ji2 <= jpi .AND. jj2 <= jpj) THEN
              p3d(ji : ji2, jj : jj2, jn) = SUM(p3d(ji : ji2, jj : jj2, jn))
            END IF
          END IF
        END DO
      END DO
    END DO
    CALL lbc_lnk('lib_fortran', p3d, 'T', 1.)
    IF (nbondi /= - 1) THEN
      IF (MOD(mig(1), 3) == 1) p3d(1, :, :) = p3d(2, :, :)
      IF (MOD(mig(1), 3) == 2) p3d(2, :, :) = p3d(1, :, :)
    END IF
    IF (nbondi /= 1) THEN
      IF (MOD(mig(jpi - 2), 3) == 1) p3d(jpi, :, :) = p3d(jpi - 1, :, :)
      IF (MOD(mig(jpi - 2), 3) == 0) p3d(jpi - 1, :, :) = p3d(jpi, :, :)
    END IF
    IF (nbondj /= - 1) THEN
      IF (MOD(mjg(1), 3) == 1) p3d(:, 1, :) = p3d(:, 2, :)
      IF (MOD(mjg(1), 3) == 2) p3d(:, 2, :) = p3d(:, 1, :)
    END IF
    IF (nbondj /= 1) THEN
      IF (MOD(mjg(jpj - 2), 3) == 1) p3d(:, jpj, :) = p3d(:, jpj - 1, :)
      IF (MOD(mjg(jpj - 2), 3) == 0) p3d(:, jpj - 1, :) = p3d(:, jpj, :)
    END IF
    CALL lbc_lnk('lib_fortran', p3d, 'T', 1.)
  END SUBROUTINE sum3x3_3d
  SUBROUTINE DDPDD(ydda, yddb)
    COMPLEX(KIND = wp), INTENT(IN) :: ydda
    COMPLEX(KIND = wp), INTENT(INOUT) :: yddb
    REAL(KIND = wp) :: zerr, zt1, zt2
    zt1 = REAL(ydda) + REAL(yddb)
    zerr = zt1 - REAL(ydda)
    zt2 = ((REAL(yddb) - zerr) + (REAL(ydda) - (zt1 - zerr))) + AIMAG(ydda) + AIMAG(yddb)
    yddb = CMPLX(zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), wp)
  END SUBROUTINE DDPDD
  FUNCTION SIGN_SCALAR(pa, pb)
    REAL(KIND = wp) :: pa, pb
    REAL(KIND = wp) :: SIGN_SCALAR
    IF (pb >= 0.E0) THEN
      sign_scalar = ABS(pa)
    ELSE
      sign_scalar = - ABS(pa)
    END IF
  END FUNCTION SIGN_SCALAR
  FUNCTION SIGN_ARRAY_1D(pa, pb)
    REAL(KIND = wp) :: pa, pb(:)
    REAL(KIND = wp) :: SIGN_ARRAY_1D(SIZE(pb, 1))
    WHERE (pb >= 0.E0)
      sign_array_1d = ABS(pa)
    ELSEWHERE
      sign_array_1d = - ABS(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_1D
  FUNCTION SIGN_ARRAY_2D(pa, pb)
    REAL(KIND = wp) :: pa, pb(:, :)
    REAL(KIND = wp) :: SIGN_ARRAY_2D(SIZE(pb, 1), SIZE(pb, 2))
    WHERE (pb >= 0.E0)
      sign_array_2d = ABS(pa)
    ELSEWHERE
      sign_array_2d = - ABS(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_2D
  FUNCTION SIGN_ARRAY_3D(pa, pb)
    REAL(KIND = wp) :: pa, pb(:, :, :)
    REAL(KIND = wp) :: SIGN_ARRAY_3D(SIZE(pb, 1), SIZE(pb, 2), SIZE(pb, 3))
    WHERE (pb >= 0.E0)
      sign_array_3d = ABS(pa)
    ELSEWHERE
      sign_array_3d = - ABS(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_3D
  FUNCTION SIGN_ARRAY_1D_A(pa, pb)
    REAL(KIND = wp) :: pa(:), pb(:)
    REAL(KIND = wp) :: SIGN_ARRAY_1D_A(SIZE(pb, 1))
    WHERE (pb >= 0.E0)
      sign_array_1d_a = ABS(pa)
    ELSEWHERE
      sign_array_1d_a = - ABS(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_1D_A
  FUNCTION SIGN_ARRAY_2D_A(pa, pb)
    REAL(KIND = wp) :: pa(:, :), pb(:, :)
    REAL(KIND = wp) :: SIGN_ARRAY_2D_A(SIZE(pb, 1), SIZE(pb, 2))
    WHERE (pb >= 0.E0)
      sign_array_2d_a = ABS(pa)
    ELSEWHERE
      sign_array_2d_a = - ABS(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_2D_A
  FUNCTION SIGN_ARRAY_3D_A(pa, pb)
    REAL(KIND = wp) :: pa(:, :, :), pb(:, :, :)
    REAL(KIND = wp) :: SIGN_ARRAY_3D_A(SIZE(pb, 1), SIZE(pb, 2), SIZE(pb, 3))
    WHERE (pb >= 0.E0)
      sign_array_3d_a = ABS(pa)
    ELSEWHERE
      sign_array_3d_a = - ABS(pa)
    END WHERE
  END FUNCTION SIGN_ARRAY_3D_A
  FUNCTION SIGN_ARRAY_1D_B(pa, pb)
    REAL(KIND = wp) :: pa(:), pb
    REAL(KIND = wp) :: SIGN_ARRAY_1D_B(SIZE(pa, 1))
    IF (pb >= 0.E0) THEN
      sign_array_1d_b = ABS(pa)
    ELSE
      sign_array_1d_b = - ABS(pa)
    END IF
  END FUNCTION SIGN_ARRAY_1D_B
  FUNCTION SIGN_ARRAY_2D_B(pa, pb)
    REAL(KIND = wp) :: pa(:, :), pb
    REAL(KIND = wp) :: SIGN_ARRAY_2D_B(SIZE(pa, 1), SIZE(pa, 2))
    IF (pb >= 0.E0) THEN
      sign_array_2d_b = ABS(pa)
    ELSE
      sign_array_2d_b = - ABS(pa)
    END IF
  END FUNCTION SIGN_ARRAY_2D_B
  FUNCTION SIGN_ARRAY_3D_B(pa, pb)
    REAL(KIND = wp) :: pa(:, :, :), pb
    REAL(KIND = wp) :: SIGN_ARRAY_3D_B(SIZE(pa, 1), SIZE(pa, 2), SIZE(pa, 3))
    IF (pb >= 0.E0) THEN
      sign_array_3d_b = ABS(pa)
    ELSE
      sign_array_3d_b = - ABS(pa)
    END IF
  END FUNCTION SIGN_ARRAY_3D_B
END MODULE lib_fortran