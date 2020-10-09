MODULE obs_rot_vel
  USE par_kind
  USE par_oce
  USE in_out_manager
  USE dom_oce
  USE obs_grid
  USE obs_utils
  USE obs_profiles_def
  USE obs_inter_h2d
  USE obs_inter_sup
  USE geo2ocean
  USE obs_fbm
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: obs_rotvel
  CONTAINS
  SUBROUTINE obs_rotvel(profdata, k2dint, pu, pv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    TYPE(obs_prof), INTENT(INOUT) :: profdata
    INTEGER, INTENT(IN) :: k2dint
    REAL(KIND = wp), DIMENSION(*) :: pu, pv
    REAL(KIND = wp), DIMENSION(2, 2, 1) :: zweig
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: zmasku, zmaskv, zcoslu, zsinlu, zcoslv, zsinlv, zglamu, zgphiu, zglamv, &
&zgphiv
    REAL(KIND = wp), DIMENSION(1) :: zsinu, zcosu, zsinv, zcosv
    REAL(KIND = wp) :: zsin
    REAL(KIND = wp) :: zcos
    REAL(KIND = wp), DIMENSION(1) :: zobsmask
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsingu, zcosgu, zsingv, zcosgv
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: igrdiu, igrdju, igrdiv, igrdjv
    INTEGER :: ji
    INTEGER :: jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('obs_rotvel', 'r0', 0, 0)
    ALLOCATE(igrdiu(2, 2, profdata % nprof), igrdju(2, 2, profdata % nprof), zglamu(2, 2, profdata % nprof), zgphiu(2, 2, profdata &
&% nprof), zmasku(2, 2, profdata % nprof), zcoslu(2, 2, profdata % nprof), zsinlu(2, 2, profdata % nprof), igrdiv(2, 2, profdata % &
&nprof), igrdjv(2, 2, profdata % nprof), zglamv(2, 2, profdata % nprof), zgphiv(2, 2, profdata % nprof), zmaskv(2, 2, profdata % &
&nprof), zcoslv(2, 2, profdata % nprof), zsinlv(2, 2, profdata % nprof))
    CALL obs_rot(zsingu, zcosgu, zsingv, zcosgv)
    DO ji = 1, profdata % nprof
      igrdiu(1, 1, ji) = profdata % mi(ji, 1) - 1
      igrdju(1, 1, ji) = profdata % mj(ji, 1) - 1
      igrdiu(1, 2, ji) = profdata % mi(ji, 1) - 1
      igrdju(1, 2, ji) = profdata % mj(ji, 1)
      igrdiu(2, 1, ji) = profdata % mi(ji, 1)
      igrdju(2, 1, ji) = profdata % mj(ji, 1) - 1
      igrdiu(2, 2, ji) = profdata % mi(ji, 1)
      igrdju(2, 2, ji) = profdata % mj(ji, 1)
      igrdiv(1, 1, ji) = profdata % mi(ji, 2) - 1
      igrdjv(1, 1, ji) = profdata % mj(ji, 2) - 1
      igrdiv(1, 2, ji) = profdata % mi(ji, 2) - 1
      igrdjv(1, 2, ji) = profdata % mj(ji, 2)
      igrdiv(2, 1, ji) = profdata % mi(ji, 2)
      igrdjv(2, 1, ji) = profdata % mj(ji, 2) - 1
      igrdiv(2, 2, ji) = profdata % mi(ji, 2)
      igrdjv(2, 2, ji) = profdata % mj(ji, 2)
    END DO
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiu, igrdju, glamu, zglamu)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiu, igrdju, gphiu, zgphiu)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiu, igrdju, umask(:, :, 1), zmasku)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiu, igrdju, zsingu, zsinlu)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiu, igrdju, zcosgu, zcoslu)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiv, igrdjv, glamv, zglamv)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiv, igrdjv, gphiv, zgphiv)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiv, igrdjv, vmask(:, :, 1), zmaskv)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiv, igrdjv, zsingv, zsinlv)
    CALL obs_int_comm_2d(2, 2, profdata % nprof, jpi, jpj, igrdiv, igrdjv, zcosgv, zcoslv)
    DO ji = 1, profdata % nprof
      CALL obs_int_h2d_init(1, 1, k2dint, profdata % rlam(ji), profdata % rphi(ji), zglamu(:, :, ji), zgphiu(:, :, ji), zmasku(:, &
&:, ji), zweig, zobsmask)
      CALL obs_int_h2d(1, 1, zweig, zsinlu(:, :, ji), zsinu)
      CALL obs_int_h2d(1, 1, zweig, zcoslu(:, :, ji), zcosu)
      CALL obs_int_h2d_init(1, 1, k2dint, profdata % rlam(ji), profdata % rphi(ji), zglamv(:, :, ji), zgphiv(:, :, ji), zmaskv(:, &
&:, ji), zweig, zobsmask)
      CALL obs_int_h2d(1, 1, zweig, zsinlv(:, :, ji), zsinv)
      CALL obs_int_h2d(1, 1, zweig, zcoslv(:, :, ji), zcosv)
      zcos = 0.5_wp * (zcosu(1) + zcosv(1))
      zsin = 0.5_wp * (zsinu(1) + zsinv(1))
      IF ((profdata % npvsta(ji, 1) /= profdata % npvsta(ji, 2)) .OR. (profdata % npvend(ji, 1) /= profdata % npvend(ji, 2))) THEN
        CALL fatal_error('Different number of U and V observations ' // 'in a profile in obs_rotvel', 190)
      END IF
      DO jk = profdata % npvsta(ji, 1), profdata % npvend(ji, 1)
        IF ((profdata % var(1) % vmod(jk) /= fbrmdi) .AND. (profdata % var(2) % vmod(jk) /= fbrmdi)) THEN
          pu(jk) = profdata % var(1) % vmod(jk) * zcos - profdata % var(2) % vmod(jk) * zsin
          pv(jk) = profdata % var(2) % vmod(jk) * zcos + profdata % var(1) % vmod(jk) * zsin
        ELSE
          pu(jk) = fbrmdi
          pv(jk) = fbrmdi
        END IF
      END DO
    END DO
    DEALLOCATE(igrdiu, igrdju, zglamu, zgphiu, zmasku, zcoslu, zsinlu, igrdiv, igrdjv, zglamv, zgphiv, zmaskv, zcoslv, zsinlv)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE obs_rotvel
END MODULE obs_rot_vel