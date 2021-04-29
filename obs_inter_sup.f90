MODULE obs_inter_sup
  USE par_kind
  USE dom_oce
  USE mpp_map
  USE lib_mpp
  USE obs_mpp
  USE obs_grid
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_int_comm_3d, obs_int_comm_2d
  CONTAINS
  SUBROUTINE obs_int_comm_3d(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval, kproc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kptsi
    INTEGER, INTENT(IN) :: kptsj
    INTEGER, INTENT(IN) :: kobs
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kpk
    INTEGER, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kgrdi, kgrdj
    INTEGER, OPTIONAL, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kproc
    REAL(KIND = wp), DIMENSION(kpi, kpj, kpk), INTENT(IN) :: pval
    REAL(KIND = wp), DIMENSION(kptsi, kptsj, kpk, kobs), INTENT(OUT) :: pgval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_comm_3d', 'r0', 0, 0)
    IF (ln_grid_global) THEN
      IF (PRESENT(kproc)) THEN
        CALL obs_int_comm_3d_global(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval, kproc = kproc)
      ELSE
        CALL obs_int_comm_3d_global(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval)
      END IF
    ELSE
      CALL obs_int_comm_3d_local(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_int_comm_3d
  SUBROUTINE obs_int_comm_2d(kptsi, kptsj, kobs, kpi, kpj, kgrdi, kgrdj, pval, pgval, kproc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kptsi
    INTEGER, INTENT(IN) :: kptsj
    INTEGER, INTENT(IN) :: kobs
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kgrdi, kgrdj
    INTEGER, OPTIONAL, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kproc
    REAL(KIND = wp), DIMENSION(kpi, kpj), INTENT(IN) :: pval
    REAL(KIND = wp), DIMENSION(kptsi, kptsj, kobs), INTENT(OUT) :: pgval
    REAL(KIND = wp), DIMENSION(jpi, jpj, 1) :: zval
    REAL(KIND = wp), DIMENSION(kptsi, kptsj, 1, kobs) :: zgval
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    zval(:, :, 1) = pval(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('obs_int_comm_2d', 'r0', 0, 0)
    IF (PRESENT(kproc)) THEN
      CALL obs_int_comm_3d(kptsi, kptsj, kobs, kpi, kpj, 1, kgrdi, kgrdj, zval, zgval, kproc = kproc)
    ELSE
      CALL obs_int_comm_3d(kptsi, kptsj, kobs, kpi, kpj, 1, kgrdi, kgrdj, zval, zgval)
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    pgval(:, :, :) = zgval(:, :, 1, :)
    !$ACC END KERNELS
  END SUBROUTINE obs_int_comm_2d
  SUBROUTINE obs_int_comm_3d_global(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval, kproc)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kptsi
    INTEGER, INTENT(IN) :: kptsj
    INTEGER, INTENT(IN) :: kobs
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kpk
    INTEGER, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kgrdi, kgrdj
    INTEGER, OPTIONAL, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kproc
    REAL(KIND = wp), DIMENSION(kpi, kpj, kpk), INTENT(IN) :: pval
    REAL(KIND = wp), DIMENSION(kptsi, kptsj, kpk, kobs), INTENT(OUT) :: pgval
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zsend, zrecv
    INTEGER, DIMENSION(:), ALLOCATABLE :: igrdij_send, igrdij_recv
    INTEGER, DIMENSION(kptsi, kptsj, kobs) :: iorder, iproc
    INTEGER :: nplocal(jpnij)
    INTEGER :: npglobal(jpnij)
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jp
    INTEGER :: jobs
    INTEGER :: it
    INTEGER :: itot
    INTEGER :: ii
    INTEGER :: ij
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('obs_int_comm_3d_global', 'r0', 0, 0)
    IF ((MAXVAL(kgrdi) > jpiglo) .OR. (MINVAL(kgrdi) < 1) .OR. (MAXVAL(kgrdj) > jpjglo) .OR. (MINVAL(kgrdj) < 1)) THEN
      CALL ctl_stop('Error in obs_int_comm_3d_global', 'Point outside global domain')
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    nplocal(:) = 0
    !$ACC END KERNELS
    IF (PRESENT(kproc)) THEN
      !$ACC KERNELS
      iproc(:, :, :) = kproc(:, :, :)
      DO jobs = 1, kobs
        !$ACC loop independent collapse(2)
        DO jj = 1, kptsj
          DO ji = 1, kptsi
            nplocal(iproc(ji, jj, jobs)) = nplocal(iproc(ji, jj, jobs)) + 1
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jobs = 1, kobs
        !$ACC loop independent collapse(2)
        DO jj = 1, kptsj
          DO ji = 1, kptsi
            iproc(ji, jj, jobs) = mppmap(kgrdi(ji, jj, jobs), kgrdj(ji, jj, jobs))
            nplocal(iproc(ji, jj, jobs)) = nplocal(iproc(ji, jj, jobs)) + 1
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data1 % PreStart('obs_int_comm_3d_global', 'r1', 0, 0)
    CALL mpp_alltoall_int(1, nplocal, npglobal)
    itot = SUM(npglobal)
    ALLOCATE(igrdij_send(kptsi * kptsj * kobs * 2), igrdij_recv(itot * 2), zsend(kpk, itot), zrecv(kpk, kptsi * kptsj * kobs))
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    it = 0
    DO jp = 1, jpnij
      DO jobs = 1, kobs
        !$ACC loop independent collapse(2)
        DO jj = 1, kptsj
          DO ji = 1, kptsi
            IF (iproc(ji, jj, jobs) == jp) THEN
              it = it + 1
              iorder(ji, jj, jobs) = it
              igrdij_send(2 * it - 1) = kgrdi(ji, jj, jobs)
              igrdij_send(2 * it) = kgrdj(ji, jj, jobs)
            END IF
          END DO
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('obs_int_comm_3d_global', 'r2', 0, 0)
    CALL mpp_alltoallv_int(igrdij_send, kptsi * kptsj * kobs * 2, nplocal(:) * 2, igrdij_recv, itot * 2, npglobal(:) * 2)
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    DO ji = 1, itot
      ii = mi1(igrdij_recv(2 * ji - 1))
      ij = mj1(igrdij_recv(2 * ji))
      DO jk = 1, kpk
        zsend(jk, ji) = pval(ii, ij, jk)
      END DO
    END DO
    nplocal(:) = kpk * nplocal(:)
    npglobal(:) = kpk * npglobal(:)
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('obs_int_comm_3d_global', 'r3', 0, 0)
    CALL mpp_alltoallv_real(zsend, kpk * itot, npglobal, zrecv, kpk * kptsi * kptsj * kobs, nplocal)
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    DO jobs = 1, kobs
      !$ACC loop independent collapse(2)
      DO jj = 1, kptsj
        DO ji = 1, kptsi
          it = iorder(ji, jj, jobs)
          DO jk = 1, kpk
            pgval(ji, jj, jk, jobs) = zrecv(jk, it)
          END DO
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    DEALLOCATE(igrdij_send, igrdij_recv, zsend, zrecv)
  END SUBROUTINE obs_int_comm_3d_global
  SUBROUTINE obs_int_comm_3d_local(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kptsi
    INTEGER, INTENT(IN) :: kptsj
    INTEGER, INTENT(IN) :: kobs
    INTEGER, INTENT(IN) :: kpi
    INTEGER, INTENT(IN) :: kpj
    INTEGER, INTENT(IN) :: kpk
    INTEGER, DIMENSION(kptsi, kptsj, kobs), INTENT(IN) :: kgrdi, kgrdj
    REAL(KIND = wp), DIMENSION(kpi, kpj, kpk), INTENT(IN) :: pval
    REAL(KIND = wp), DIMENSION(kptsi, kptsj, kpk, kobs), INTENT(OUT) :: pgval
    INTEGER :: ji
    INTEGER :: jj
    INTEGER :: jk
    INTEGER :: jobs
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_int_comm_3d_local', 'r0', 0, 0)
    IF ((MAXVAL(kgrdi) > jpi) .OR. (MINVAL(kgrdi) < 1) .OR. (MAXVAL(kgrdj) > jpj) .OR. (MINVAL(kgrdj) < 1)) THEN
      CALL ctl_stop('Error in obs_int_comm_3d_local', 'Point outside local domain')
    END IF
    CALL profile_psy_data0 % PostEnd
    DO jobs = 1, kobs
      DO jj = 1, kptsj
        DO ji = 1, kptsi
          !$ACC KERNELS
          DO jk = 1, kpk
            pgval(ji, jj, jk, jobs) = pval(kgrdi(ji, jj, jobs), kgrdj(ji, jj, jobs), jk)
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END DO
  END SUBROUTINE obs_int_comm_3d_local
END MODULE obs_inter_sup