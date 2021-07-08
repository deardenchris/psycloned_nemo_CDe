MODULE icethd_ent
  USE dom_oce
  USE domain
  USE phycst
  USE ice
  USE ice1D
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_ent
  CONTAINS
  SUBROUTINE ice_thd_ent(qnew)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: qnew
    INTEGER :: ji
    INTEGER :: jk0, jk1
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i + 2) :: zeh_cum0, zh_cum0
    REAL(KIND = wp), DIMENSION(jpij, 0 : nlay_i) :: zeh_cum1, zh_cum1
    REAL(KIND = wp), DIMENSION(jpij) :: zhnew
!    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
!    CALL profile_psy_data0 % PreStart('ice_thd_ent', 'r0', 0, 0)
    !$ACC KERNELS ! CDe added
    zeh_cum0(1 : npti, 0) = 0._wp
    zh_cum0(1 : npti, 0) = 0._wp
    DO jk0 = 1, nlay_i + 2
      DO ji = 1, npti
        zeh_cum0(ji, jk0) = zeh_cum0(ji, jk0 - 1) + eh_i_old(ji, jk0 - 1)
        zh_cum0(ji, jk0) = zh_cum0(ji, jk0 - 1) + h_i_old(ji, jk0 - 1)
      END DO
    END DO
    DO ji = 1, npti
      zhnew(ji) = SUM(h_i_old(ji, 0 : nlay_i + 1)) * r1_nlay_i
    END DO
    zh_cum1(1 : npti, 0) = 0._wp
    DO jk1 = 1, nlay_i
      DO ji = 1, npti
        zh_cum1(ji, jk1) = zh_cum1(ji, jk1 - 1) + zhnew(ji)
      END DO
    END DO
    zeh_cum1(1 : npti, 0 : nlay_i) = 0._wp
    DO jk0 = 1, nlay_i + 2
      DO jk1 = 1, nlay_i - 1
        DO ji = 1, npti
          IF (zh_cum1(ji, jk1) <= zh_cum0(ji, jk0) .AND. zh_cum1(ji, jk1) > zh_cum0(ji, jk0 - 1)) THEN
            zeh_cum1(ji, jk1) = (zeh_cum0(ji, jk0 - 1) * (zh_cum0(ji, jk0) - zh_cum1(ji, jk1)) + zeh_cum0(ji, jk0) * (zh_cum1(ji, jk1) - zh_cum0(ji, jk0 - 1))) / (zh_cum0(ji, jk0) - zh_cum0(ji, jk0 - 1))
          END IF
        END DO
      END DO
    END DO
    zeh_cum1(1 : npti, nlay_i) = zeh_cum0(1 : npti, nlay_i + 2)
    DO jk1 = 1, nlay_i
      DO ji = 1, npti
        rswitch = MAX(0._wp, SIGN(1._wp, zhnew(ji) - epsi20))
        qnew(ji, jk1) = rswitch * (zeh_cum1(ji, jk1) - zeh_cum1(ji, jk1 - 1)) / MAX(zhnew(ji), epsi20)
      END DO
    END DO
    DO ji = 1, npti
      hfx_err_rem_1d(ji) = hfx_err_rem_1d(ji) + a_i_1d(ji) * r1_rdtice * (SUM(qnew(ji, 1 : nlay_i)) * zhnew(ji) - SUM(eh_i_old(ji, 0 : nlay_i + 1)))
    END DO
!    CALL profile_psy_data0 % PostEnd
    !$ACC END KERNELS
  END SUBROUTINE ice_thd_ent
END MODULE icethd_ent
