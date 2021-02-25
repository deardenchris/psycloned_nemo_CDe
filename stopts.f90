MODULE stopts
  USE dom_oce
  USE lbclnk
  USE phycst
  USE stopar
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sto_pts
  PUBLIC :: sto_pts_init
  REAL(KIND = wp), PUBLIC, DIMENSION(:, :, :, :, :), ALLOCATABLE :: pts_ran
  CONTAINS
  SUBROUTINE sto_pts(pts)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts), INTENT(INOUT) :: pts
    INTEGER :: ji, jj, jk, jts, jdof
    INTEGER :: jim1, jjm1, jkm1
    INTEGER :: jip1, jjp1, jkp1
    REAL(KIND = wp) :: zdtsim, zdtsjm, zdtskm
    REAL(KIND = wp) :: zdtsip, zdtsjp, zdtskp, zdts
    DO jts = 1, jpts
      CALL lbc_lnk('stopts', pts(:, :, :, jts), 'T', 1._wp)
    END DO
    DO jdof = 1, nn_sto_eos
      DO jts = 1, jpts
        !$OMP parallel default(shared),  &
!$omp& private(ji,jim1,jip1,jj,jjm1,jjp1,jk,jkm1,jkp1,zdts,zdtsim,zdtsip,zdtsjm,zdtsjp,zdtskm,zdtskp)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          jkm1 = MAX(jk - 1, 1)
          jkp1 = MIN(jk + 1, jpkm1)
          DO jj = 1, jpj
            jjm1 = MAX(jj - 1, 1)
            jjp1 = MIN(jj + 1, jpj)
            DO ji = 1, jpi
              jim1 = MAX(ji - 1, 1)
              jip1 = MIN(ji + 1, jpi)
              zdtsip = (pts(jip1, jj, jk, jts) - pts(ji, jj, jk, jts)) * tmask(jip1, jj, jk)
              zdtsim = (pts(ji, jj, jk, jts) - pts(jim1, jj, jk, jts)) * tmask(jim1, jj, jk)
              zdtsjp = (pts(ji, jjp1, jk, jts) - pts(ji, jj, jk, jts)) * tmask(ji, jjp1, jk)
              zdtsjm = (pts(ji, jj, jk, jts) - pts(ji, jjm1, jk, jts)) * tmask(ji, jjm1, jk)
              zdtskp = (pts(ji, jj, jkp1, jts) - pts(ji, jj, jk, jts)) * tmask(ji, jj, jkp1)
              zdtskm = (pts(ji, jj, jk, jts) - pts(ji, jj, jkm1, jts)) * tmask(ji, jj, jkm1)
              zdts = (zdtsip + zdtsim) * sto2d(ji, jj, jsto_eosi(jdof)) + (zdtsjp + zdtsjm) * sto2d(ji, jj, jsto_eosj(jdof)) + &
&(zdtskp + zdtskm) * sto2d(ji, jj, jsto_eosk(jdof))
              zdts = zdts * tmask(ji, jj, jk) * SIN(gphit(ji, jj) * rad)
              pts_ran(ji, jj, jk, jts, jdof) = zdts * 0.5_wp
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      END DO
    END DO
    DO jdof = 1, nn_sto_eos
      !$OMP parallel default(shared), private(ji,jj,jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            pts_ran(ji, jj, jk, jp_sal, jdof) = MIN(ABS(pts_ran(ji, jj, jk, jp_sal, jdof)), MAX(pts(ji, jj, jk, jp_sal), 0._wp)) * &
&SIGN(1._wp, pts_ran(ji, jj, jk, jp_sal, jdof))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    END DO
    DO jdof = 1, nn_sto_eos
      DO jts = 1, jpts
        CALL lbc_lnk('stopts', pts_ran(:, :, :, jts, jdof), 'T', 1._wp)
      END DO
    END DO
  END SUBROUTINE sto_pts
  SUBROUTINE sto_pts_init
    ALLOCATE(pts_ran(jpi, jpj, jpk, jpts, nn_sto_eos))
  END SUBROUTINE sto_pts_init
END MODULE stopts