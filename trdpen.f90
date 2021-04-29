MODULE trdpen
  USE oce
  USE dom_oce
  USE sbc_oce
  USE zdf_oce
  USE trd_oce
  USE eosbn2
  USE ldftra
  USE zdfddm
  USE phycst
  USE in_out_manager
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_pen
  PUBLIC :: trd_pen_init
  INTEGER :: nkstp
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :, :) :: rab_pe
  CONTAINS
  INTEGER FUNCTION trd_pen_alloc()
    ALLOCATE(rab_pe(jpi, jpj, jpk, jpts), STAT = trd_pen_alloc)
    CALL mpp_sum('trdpen', trd_pen_alloc)
    IF (trd_pen_alloc /= 0) CALL ctl_stop('STOP', 'trd_pen_alloc: failed to allocate arrays')
  END FUNCTION trd_pen_alloc
  SUBROUTINE trd_pen(ptrdx, ptrdy, ktrd, kt, pdt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptrdx, ptrdy
    INTEGER, INTENT(IN) :: ktrd
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), INTENT(IN) :: pdt
    INTEGER :: jk
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zpe
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    !$ACC KERNELS
    zpe(:, :, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('trd_pen', 'r0', 0, 0)
    IF (kt /= nkstp) THEN
      nkstp = kt
      CALL eos_pen(tsn, rab_PE, zpe)
      CALL iom_put("alphaPE", rab_pe(:, :, :, jp_tem))
      CALL iom_put("betaPE", rab_pe(:, :, :, jp_sal))
      CALL iom_put("PEanom", zpe)
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zpe(:, :, jpk) = 0._wp
    DO jk = 1, jpkm1
      zpe(:, :, jk) = (- (rab_n(:, :, jk, jp_tem) + rab_pe(:, :, jk, jp_tem)) * ptrdx(:, :, jk) + (rab_n(:, :, jk, jp_sal) + rab_pe(:, :, jk, jp_sal)) * ptrdy(:, :, jk))
    END DO
    !$ACC END KERNELS
    SELECT CASE (ktrd)
    CASE (jptra_xad)
      CALL profile_psy_data1 % PreStart('trd_pen', 'r1', 0, 0)
      CALL iom_put("petrd_xad", zpe)
      CALL profile_psy_data1 % PostEnd
    CASE (jptra_yad)
      CALL profile_psy_data2 % PreStart('trd_pen', 'r2', 0, 0)
      CALL iom_put("petrd_yad", zpe)
      CALL profile_psy_data2 % PostEnd
    CASE (jptra_zad)
      CALL profile_psy_data3 % PreStart('trd_pen', 'r3', 0, 0)
      CALL iom_put("petrd_zad", zpe)
      CALL profile_psy_data3 % PostEnd
      IF (ln_linssh) THEN
        ALLOCATE(z2d(jpi, jpj))
        !$ACC KERNELS
        z2d(:, :) = wn(:, :, 1) * (- (rab_n(:, :, 1, jp_tem) + rab_pe(:, :, 1, jp_tem)) * tsn(:, :, 1, jp_tem) + (rab_n(:, :, 1, jp_sal) + rab_pe(:, :, 1, jp_sal)) * tsn(:, :, 1, jp_sal)) / e3t_n(:, :, 1)
        !$ACC END KERNELS
        CALL profile_psy_data4 % PreStart('trd_pen', 'r4', 0, 0)
        CALL iom_put("petrd_sad", z2d)
        DEALLOCATE(z2d)
        CALL profile_psy_data4 % PostEnd
      END IF
    CASE (jptra_ldf)
      CALL profile_psy_data5 % PreStart('trd_pen', 'r5', 0, 0)
      CALL iom_put("petrd_ldf", zpe)
      CALL profile_psy_data5 % PostEnd
    CASE (jptra_zdf)
      CALL iom_put("petrd_zdf", zpe)
    CASE (jptra_zdfp)
      CALL iom_put("petrd_zdfp", zpe)
    CASE (jptra_dmp)
      CALL iom_put("petrd_dmp", zpe)
    CASE (jptra_bbl)
      CALL iom_put("petrd_bbl", zpe)
    CASE (jptra_npc)
      CALL iom_put("petrd_npc", zpe)
    CASE (jptra_nsr)
      CALL iom_put("petrd_nsr", zpe)
    CASE (jptra_qsr)
      CALL iom_put("petrd_qsr", zpe)
    CASE (jptra_bbc)
      CALL iom_put("petrd_bbc", zpe)
    CASE (jptra_atf)
      CALL iom_put("petrd_atf", zpe)
    END SELECT
  END SUBROUTINE trd_pen
  SUBROUTINE trd_pen_init
    INTEGER :: ji, jj, jk
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'trd_pen_init : 3D Potential ENergy trends'
      WRITE(numout, *) '~~~~~~~~~~~~~'
    END IF
    IF (trd_pen_alloc() /= 0) CALL ctl_stop('trd_pen_alloc: failed to allocate arrays')
    !$ACC KERNELS
    rab_pe(:, :, :, :) = 0._wp
    !$ACC END KERNELS
    IF (.NOT. ln_linssh) CALL ctl_stop('trd_pen_init : PE trends not coded for variable volume')
    nkstp = nit000 - 1
  END SUBROUTINE trd_pen_init
END MODULE trdpen