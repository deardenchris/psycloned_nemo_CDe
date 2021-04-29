MODULE trdken
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE zdf_oce
  USE ldftra
  USE trd_oce
  USE trdvor
  USE trdglo
  USE trdmxl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE ldfslp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_ken
  PUBLIC :: trd_ken_init
  INTEGER :: nkstp
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: bu, bv
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: r1_bt
  CONTAINS
  INTEGER FUNCTION trd_ken_alloc()
    ALLOCATE(bu(jpi, jpj, jpk), bv(jpi, jpj, jpk), r1_bt(jpi, jpj, jpk), STAT = trd_ken_alloc)
    CALL mpp_sum('trdken', trd_ken_alloc)
    IF (trd_ken_alloc /= 0) CALL ctl_stop('STOP', 'trd_ken_alloc: failed to allocate arrays')
  END FUNCTION trd_ken_alloc
  SUBROUTINE trd_ken(putrd, pvtrd, ktrd, kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbu, ikbv
    INTEGER :: ikbum1, ikbvm1
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: z2dx, z2dy, zke2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zke
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
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data11
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data12
    CALL profile_psy_data0 % PreStart('trd_ken', 'r0', 0, 0)
    CALL lbc_lnk_multi('trdken', putrd, 'U', - 1., pvtrd, 'V', - 1.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    nkstp = kt
    DO jk = 1, jpkm1
      bu(:, :, jk) = e1e2u(:, :) * e3u_n(:, :, jk)
      bv(:, :, jk) = e1e2v(:, :) * e3v_n(:, :, jk)
      r1_bt(:, :, jk) = r1_e1e2t(:, :) / e3t_n(:, :, jk) * tmask(:, :, jk)
    END DO
    zke(:, :, jpk) = 0._wp
    zke(1, :, :) = 0._wp
    zke(:, 1, :) = 0._wp
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          zke(ji, jj, jk) = 0.5_wp * rau0 * (un(ji, jj, jk) * putrd(ji, jj, jk) * bu(ji, jj, jk) + un(ji - 1, jj, jk) * putrd(ji - 1, jj, jk) * bu(ji - 1, jj, jk) + vn(ji, jj, jk) * pvtrd(ji, jj, jk) * bv(ji, jj, jk) + vn(ji, jj - 1, jk) * pvtrd(ji, jj - 1, jk) * bv(ji, jj - 1, jk)) * r1_bt(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    SELECT CASE (ktrd)
    CASE (jpdyn_hpg)
      CALL profile_psy_data1 % PreStart('trd_ken', 'r1', 0, 0)
      CALL iom_put("ketrd_hpg", zke)
      CALL profile_psy_data1 % PostEnd
    CASE (jpdyn_spg)
      CALL profile_psy_data2 % PreStart('trd_ken', 'r2', 0, 0)
      CALL iom_put("ketrd_spg", zke)
      CALL profile_psy_data2 % PostEnd
    CASE (jpdyn_pvo)
      CALL profile_psy_data3 % PreStart('trd_ken', 'r3', 0, 0)
      CALL iom_put("ketrd_pvo", zke)
      CALL profile_psy_data3 % PostEnd
    CASE (jpdyn_rvo)
      CALL profile_psy_data4 % PreStart('trd_ken', 'r4', 0, 0)
      CALL iom_put("ketrd_rvo", zke)
      CALL profile_psy_data4 % PostEnd
    CASE (jpdyn_keg)
      CALL profile_psy_data5 % PreStart('trd_ken', 'r5', 0, 0)
      CALL iom_put("ketrd_keg", zke)
      CALL profile_psy_data5 % PostEnd
    CASE (jpdyn_zad)
      CALL profile_psy_data6 % PreStart('trd_ken', 'r6', 0, 0)
      CALL iom_put("ketrd_zad", zke)
      CALL profile_psy_data6 % PostEnd
    CASE (jpdyn_ldf)
      CALL profile_psy_data7 % PreStart('trd_ken', 'r7', 0, 0)
      CALL iom_put("ketrd_ldf", zke)
      CALL profile_psy_data7 % PostEnd
    CASE (jpdyn_zdf)
      CALL profile_psy_data8 % PreStart('trd_ken', 'r8', 0, 0)
      CALL iom_put("ketrd_zdf", zke)
      ALLOCATE(z2dx(jpi, jpj), z2dy(jpi, jpj), zke2d(jpi, jpj))
      CALL profile_psy_data8 % PostEnd
      !$ACC KERNELS
      z2dx(:, :) = un(:, :, 1) * (utau_b(:, :) + utau(:, :)) * e1e2u(:, :) * umask(:, :, 1)
      z2dy(:, :) = vn(:, :, 1) * (vtau_b(:, :) + vtau(:, :)) * e1e2v(:, :) * vmask(:, :, 1)
      zke2d(1, :) = 0._wp
      zke2d(:, 1) = 0._wp
      !$ACC loop independent collapse(2)
      DO jj = 2, jpj
        DO ji = 2, jpi
          zke2d(ji, jj) = r1_rau0 * 0.5_wp * (z2dx(ji, jj) + z2dx(ji - 1, jj) + z2dy(ji, jj) + z2dy(ji, jj - 1)) * r1_bt(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data9 % PreStart('trd_ken', 'r9', 0, 0)
      CALL iom_put("ketrd_tau", zke2d)
      DEALLOCATE(z2dx, z2dy, zke2d)
      CALL profile_psy_data9 % PostEnd
    CASE (jpdyn_bfr)
      CALL profile_psy_data10 % PreStart('trd_ken', 'r10', 0, 0)
      CALL iom_put("ketrd_bfr", zke)
      CALL profile_psy_data10 % PostEnd
    CASE (jpdyn_atf)
      CALL profile_psy_data11 % PreStart('trd_ken', 'r11', 0, 0)
      CALL iom_put("ketrd_atf", zke)
      CALL profile_psy_data11 % PostEnd
    CASE (jpdyn_ken)
      !$ACC KERNELS
      zke(:, :, :) = 0.5_wp * zke(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data12 % PreStart('trd_ken', 'r12', 0, 0)
      CALL iom_put("KE", zke)
      CALL ken_p2k(kt, zke)
      CALL iom_put("ketrd_convP2K", zke)
      CALL profile_psy_data12 % PostEnd
    END SELECT
  END SUBROUTINE trd_ken
  SUBROUTINE ken_p2k(kt, pconv)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pconv
    INTEGER :: ji, jj, jk
    INTEGER :: iku, ikv
    REAL(KIND = wp) :: zcoef
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zconv
    !$ACC KERNELS
    zcoef = - rau0 * grav * 0.5_wp
    zconv(:, :, 1) = zcoef * (2._wp * rhd(:, :, 1)) * wn(:, :, 1) * e3w_n(:, :, 1)
    DO jk = 2, jpk
      zconv(:, :, jk) = zcoef * (rhd(:, :, jk) + rhd(:, :, jk - 1)) * wn(:, :, jk) * e3w_n(:, :, jk)
    END DO
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zcoef = 0.5_wp / e3t_n(ji, jj, jk)
          pconv(ji, jj, jk) = zcoef * (zconv(ji, jj, jk) + zconv(ji, jj, jk + 1)) * tmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE ken_p2k
  SUBROUTINE trd_ken_init
    INTEGER :: ji, jj, jk
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'trd_ken_init : 3D Kinetic Energy trends'
      WRITE(numout, *) '~~~~~~~~~~~~~'
    END IF
    IF (trd_ken_alloc() /= 0) CALL ctl_stop('trd_ken_alloc: failed to allocate arrays')
  END SUBROUTINE trd_ken_init
END MODULE trdken