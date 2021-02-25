MODULE trddyn
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE zdf_oce
  USE trd_oce
  USE trdken
  USE trdglo
  USE trdvor
  USE trdmxl
  USE in_out_manager
  USE lbclnk
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_dyn
  CONTAINS
  SUBROUTINE trd_dyn(putrd, pvtrd, ktrd, kt)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    putrd(:, :, :) = putrd(:, :, :) * umask(:, :, :)
    pvtrd(:, :, :) = pvtrd(:, :, :) * vmask(:, :, :)
    IF (ln_dyn_trd) CALL trd_dyn_iom(putrd, pvtrd, ktrd, kt)
    IF (ln_glo_trd) CALL trd_glo(putrd, pvtrd, ktrd, 'DYN', kt)
    IF (ln_KE_trd) CALL trd_ken(putrd, pvtrd, ktrd, kt)
    IF (ln_vor_trd) CALL trd_vor(putrd, pvtrd, ktrd, kt)
  END SUBROUTINE trd_dyn
  SUBROUTINE trd_dyn_iom(putrd, pvtrd, ktrd, kt)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2dx, z2dy
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z3dx, z3dy
    SELECT CASE (ktrd)
    CASE (jpdyn_hpg)
      CALL iom_put("utrd_hpg", putrd)
      CALL iom_put("vtrd_hpg", pvtrd)
    CASE (jpdyn_spg)
      CALL iom_put("utrd_spg", putrd)
      CALL iom_put("vtrd_spg", pvtrd)
    CASE (jpdyn_pvo)
      CALL iom_put("utrd_pvo", putrd)
      CALL iom_put("vtrd_pvo", pvtrd)
    CASE (jpdyn_rvo)
      CALL iom_put("utrd_rvo", putrd)
      CALL iom_put("vtrd_rvo", pvtrd)
    CASE (jpdyn_keg)
      CALL iom_put("utrd_keg", putrd)
      CALL iom_put("vtrd_keg", pvtrd)
      ALLOCATE(z3dx(jpi, jpj, jpk), z3dy(jpi, jpj, jpk))
      z3dx(:, :, :) = 0._wp
      z3dy(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z3dx(ji, jj, jk) = un(ji, jj, jk) * (un(ji + 1, jj, jk) - un(ji - 1, jj, jk)) / (2._wp * e1u(ji, jj))
            z3dy(ji, jj, jk) = vn(ji, jj, jk) * (vn(ji, jj + 1, jk) - vn(ji, jj - 1, jk)) / (2._wp * e2v(ji, jj))
          END DO
        END DO
      END DO
      CALL lbc_lnk_multi('trddyn', z3dx, 'U', - 1., z3dy, 'V', - 1.)
      CALL iom_put("utrd_udx", z3dx)
      CALL iom_put("vtrd_vdy", z3dy)
      DEALLOCATE(z3dx, z3dy)
    CASE (jpdyn_zad)
      CALL iom_put("utrd_zad", putrd)
      CALL iom_put("vtrd_zad", pvtrd)
    CASE (jpdyn_ldf)
      CALL iom_put("utrd_ldf", putrd)
      CALL iom_put("vtrd_ldf", pvtrd)
    CASE (jpdyn_zdf)
      CALL iom_put("utrd_zdf", putrd)
      CALL iom_put("vtrd_zdf", pvtrd)
      ALLOCATE(z2dx(jpi, jpj), z2dy(jpi, jpj))
      z2dx(:, :) = (utau_b(:, :) + utau(:, :)) / (e3u_n(:, :, 1) * rau0)
      z2dy(:, :) = (vtau_b(:, :) + vtau(:, :)) / (e3v_n(:, :, 1) * rau0)
      CALL iom_put("utrd_tau", z2dx)
      CALL iom_put("vtrd_tau", z2dy)
      DEALLOCATE(z2dx, z2dy)
    CASE (jpdyn_bfr)
      CALL iom_put("utrd_bfr", putrd)
      CALL iom_put("vtrd_bfr", pvtrd)
    CASE (jpdyn_atf)
      CALL iom_put("utrd_atf", putrd)
      CALL iom_put("vtrd_atf", pvtrd)
    END SELECT
  END SUBROUTINE trd_dyn_iom
END MODULE trddyn