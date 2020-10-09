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
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    putrd(:, :, :) = putrd(:, :, :) * umask(:, :, :)
    pvtrd(:, :, :) = pvtrd(:, :, :) * vmask(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('trd_dyn', 'r0', 0, 0)
    IF (ln_dyn_trd) CALL trd_dyn_iom(putrd, pvtrd, ktrd, kt)
    IF (ln_glo_trd) CALL trd_glo(putrd, pvtrd, ktrd, 'DYN', kt)
    IF (ln_KE_trd) CALL trd_ken(putrd, pvtrd, ktrd, kt)
    IF (ln_vor_trd) CALL trd_vor(putrd, pvtrd, ktrd, kt)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE trd_dyn
  SUBROUTINE trd_dyn_iom(putrd, pvtrd, ktrd, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2dx, z2dy
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z3dx, z3dy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data9
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data10
    SELECT CASE (ktrd)
    CASE (jpdyn_hpg)
      CALL profile_psy_data0 % PreStart('trd_dyn_iom', 'r0', 0, 0)
      CALL iom_put("utrd_hpg", putrd)
      CALL iom_put("vtrd_hpg", pvtrd)
      CALL profile_psy_data0 % PostEnd
    CASE (jpdyn_spg)
      CALL profile_psy_data1 % PreStart('trd_dyn_iom', 'r1', 0, 0)
      CALL iom_put("utrd_spg", putrd)
      CALL iom_put("vtrd_spg", pvtrd)
      CALL profile_psy_data1 % PostEnd
    CASE (jpdyn_pvo)
      CALL profile_psy_data2 % PreStart('trd_dyn_iom', 'r2', 0, 0)
      CALL iom_put("utrd_pvo", putrd)
      CALL iom_put("vtrd_pvo", pvtrd)
      CALL profile_psy_data2 % PostEnd
    CASE (jpdyn_rvo)
      CALL profile_psy_data3 % PreStart('trd_dyn_iom', 'r3', 0, 0)
      CALL iom_put("utrd_rvo", putrd)
      CALL iom_put("vtrd_rvo", pvtrd)
      CALL profile_psy_data3 % PostEnd
    CASE (jpdyn_keg)
      CALL profile_psy_data4 % PreStart('trd_dyn_iom', 'r4', 0, 0)
      CALL iom_put("utrd_keg", putrd)
      CALL iom_put("vtrd_keg", pvtrd)
      ALLOCATE(z3dx(jpi, jpj, jpk), z3dy(jpi, jpj, jpk))
      CALL profile_psy_data4 % PostEnd
      !$ACC KERNELS
      z3dx(:, :, :) = 0._wp
      z3dy(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z3dx(ji, jj, jk) = un(ji, jj, jk) * (un(ji + 1, jj, jk) - un(ji - 1, jj, jk)) / (2._wp * e1u(ji, jj))
            z3dy(ji, jj, jk) = vn(ji, jj, jk) * (vn(ji, jj + 1, jk) - vn(ji, jj - 1, jk)) / (2._wp * e2v(ji, jj))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('trd_dyn_iom', 'r5', 0, 0)
      CALL lbc_lnk_multi('trddyn', z3dx, 'U', - 1., z3dy, 'V', - 1.)
      CALL iom_put("utrd_udx", z3dx)
      CALL iom_put("vtrd_vdy", z3dy)
      DEALLOCATE(z3dx, z3dy)
      CALL profile_psy_data5 % PostEnd
    CASE (jpdyn_zad)
      CALL profile_psy_data6 % PreStart('trd_dyn_iom', 'r6', 0, 0)
      CALL iom_put("utrd_zad", putrd)
      CALL iom_put("vtrd_zad", pvtrd)
      CALL profile_psy_data6 % PostEnd
    CASE (jpdyn_ldf)
      CALL profile_psy_data7 % PreStart('trd_dyn_iom', 'r7', 0, 0)
      CALL iom_put("utrd_ldf", putrd)
      CALL iom_put("vtrd_ldf", pvtrd)
      CALL profile_psy_data7 % PostEnd
    CASE (jpdyn_zdf)
      CALL profile_psy_data8 % PreStart('trd_dyn_iom', 'r8', 0, 0)
      CALL iom_put("utrd_zdf", putrd)
      CALL iom_put("vtrd_zdf", pvtrd)
      ALLOCATE(z2dx(jpi, jpj), z2dy(jpi, jpj))
      CALL profile_psy_data8 % PostEnd
      !$ACC KERNELS
      z2dx(:, :) = (utau_b(:, :) + utau(:, :)) / (e3u_n(:, :, 1) * rau0)
      z2dy(:, :) = (vtau_b(:, :) + vtau(:, :)) / (e3v_n(:, :, 1) * rau0)
      !$ACC END KERNELS
      CALL profile_psy_data9 % PreStart('trd_dyn_iom', 'r9', 0, 0)
      CALL iom_put("utrd_tau", z2dx)
      CALL iom_put("vtrd_tau", z2dy)
      DEALLOCATE(z2dx, z2dy)
      CALL profile_psy_data9 % PostEnd
    CASE (jpdyn_bfr)
      CALL profile_psy_data10 % PreStart('trd_dyn_iom', 'r10', 0, 0)
      CALL iom_put("utrd_bfr", putrd)
      CALL iom_put("vtrd_bfr", pvtrd)
      CALL profile_psy_data10 % PostEnd
    CASE (jpdyn_atf)
      CALL iom_put("utrd_atf", putrd)
      CALL iom_put("vtrd_atf", pvtrd)
    END SELECT
  END SUBROUTINE trd_dyn_iom
END MODULE trddyn