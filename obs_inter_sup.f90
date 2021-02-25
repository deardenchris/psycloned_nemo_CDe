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
    IF (ln_grid_global) THEN
      IF (PRESENT(kproc)) THEN
        CALL obs_int_comm_3d_global(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval, kproc = kproc)
      ELSE
        CALL obs_int_comm_3d_global(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval)
      END IF
    ELSE
      CALL obs_int_comm_3d_local(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval)
    END IF
  END SUBROUTINE obs_int_comm_3d
  SUBROUTINE obs_int_comm_2d(kptsi, kptsj, kobs, kpi, kpj, kgrdi, kgrdj, pval, pgval, kproc)
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
    zval(:, :, 1) = pval(:, :)
    IF (PRESENT(kproc)) THEN
      CALL obs_int_comm_3d(kptsi, kptsj, kobs, kpi, kpj, 1, kgrdi, kgrdj, zval, zgval, kproc = kproc)
    ELSE
      CALL obs_int_comm_3d(kptsi, kptsj, kobs, kpi, kpj, 1, kgrdi, kgrdj, zval, zgval)
    END IF
    pgval(:, :, :) = zgval(:, :, 1, :)
  END SUBROUTINE obs_int_comm_2d
  SUBROUTINE obs_int_comm_3d_global(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval, kproc)
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
    IF ((MAXVAL(kgrdi) > jpiglo) .OR. (MINVAL(kgrdi) < 1) .OR. (MAXVAL(kgrdj) > jpjglo) .OR. (MINVAL(kgrdj) < 1)) THEN
      CALL ctl_stop('Error in obs_int_comm_3d_global', 'Point outside global domain')
    END IF
    nplocal(:) = 0
    IF (PRESENT(kproc)) THEN
      iproc(:, :, :) = kproc(:, :, :)
      DO jobs = 1, kobs
        DO jj = 1, kptsj
          DO ji = 1, kptsi
            nplocal(iproc(ji, jj, jobs)) = nplocal(iproc(ji, jj, jobs)) + 1
          END DO
        END DO
      END DO
    ELSE
      DO jobs = 1, kobs
        DO jj = 1, kptsj
          DO ji = 1, kptsi
            iproc(ji, jj, jobs) = mppmap(kgrdi(ji, jj, jobs), kgrdj(ji, jj, jobs))
            nplocal(iproc(ji, jj, jobs)) = nplocal(iproc(ji, jj, jobs)) + 1
          END DO
        END DO
      END DO
    END IF
    CALL mpp_alltoall_int(1, nplocal, npglobal)
    itot = SUM(npglobal)
    ALLOCATE(igrdij_send(kptsi * kptsj * kobs * 2), igrdij_recv(itot * 2), zsend(kpk, itot), zrecv(kpk, kptsi * kptsj * kobs))
    it = 0
    DO jp = 1, jpnij
      DO jobs = 1, kobs
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
    CALL mpp_alltoallv_int(igrdij_send, kptsi * kptsj * kobs * 2, nplocal(:) * 2, igrdij_recv, itot * 2, npglobal(:) * 2)
    DO ji = 1, itot
      ii = mi1(igrdij_recv(2 * ji - 1))
      ij = mj1(igrdij_recv(2 * ji))
      DO jk = 1, kpk
        zsend(jk, ji) = pval(ii, ij, jk)
      END DO
    END DO
    nplocal(:) = kpk * nplocal(:)
    npglobal(:) = kpk * npglobal(:)
    CALL mpp_alltoallv_real(zsend, kpk * itot, npglobal, zrecv, kpk * kptsi * kptsj * kobs, nplocal)
    DO jobs = 1, kobs
      DO jj = 1, kptsj
        DO ji = 1, kptsi
          it = iorder(ji, jj, jobs)
          DO jk = 1, kpk
            pgval(ji, jj, jk, jobs) = zrecv(jk, it)
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(igrdij_send, igrdij_recv, zsend, zrecv)
  END SUBROUTINE obs_int_comm_3d_global
  SUBROUTINE obs_int_comm_3d_local(kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, pval, pgval)
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
    IF ((MAXVAL(kgrdi) > jpi) .OR. (MINVAL(kgrdi) < 1) .OR. (MAXVAL(kgrdj) > jpj) .OR. (MINVAL(kgrdj) < 1)) THEN
      CALL ctl_stop('Error in obs_int_comm_3d_local', 'Point outside local domain')
    END IF
    DO jobs = 1, kobs
      DO jj = 1, kptsj
        DO ji = 1, kptsi
          DO jk = 1, kpk
            pgval(ji, jj, jk, jobs) = pval(kgrdi(ji, jj, jobs), kgrdj(ji, jj, jobs), jk)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE obs_int_comm_3d_local
END MODULE obs_inter_sup