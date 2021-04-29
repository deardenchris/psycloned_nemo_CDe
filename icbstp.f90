MODULE icbstp
  USE par_oce
  USE dom_oce
  USE sbc_oce
  USE phycst
  USE icb_oce
  USE icbini
  USE icbutl
  USE icbrst
  USE icbdyn
  USE icbclv
  USE icbthm
  USE icblbc
  USE icbtrj
  USE icbdia
  USE in_out_manager
  USE lib_mpp
  USE iom
  USE fldread
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_stp
  PUBLIC :: icb_end
  CONTAINS
  SUBROUTINE icb_stp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    LOGICAL :: ll_sample_traj, ll_budget, ll_verbose
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('icb_stp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('icb_stp')
    nktberg = kt
    CALL profile_psy_data0 % PostEnd
    IF (nn_test_icebergs < 0 .OR. ln_use_calving) THEN
      CALL profile_psy_data1 % PreStart('icb_stp', 'r1', 0, 0)
      CALL fld_read(kt, 1, sf_icb)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      src_calving(:, :) = sf_icb(1) % fnow(:, :, 1)
      src_calving_hflx(:, :) = 0._wp
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    berg_grid % floating_melt(:, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('icb_stp', 'r2', 0, 0)
    CALL icb_dia_step
    ll_verbose = .FALSE.
    IF (nn_verbose_write > 0 .AND. MOD(kt - 1, nn_verbose_write) == 0) ll_verbose = (nn_verbose_level >= 0)
    IF (ll_verbose) WRITE(numicb, 9100) nktberg, ndastp, nsec_day
9100 FORMAT('kt= ', I8, ' day= ', I8, ' secs=', I8)
    CALL icb_utl_copy
    CALL icb_clv_flx(kt)
    CALL icb_clv(kt)
    IF (ASSOCIATED(first_berg)) CALL icb_dyn(kt)
    IF (lk_mpp) THEN
      CALL icb_lbc_mpp
    ELSE
      CALL icb_lbc
    END IF
    IF (ASSOCIATED(first_berg)) CALL icb_thm(kt)
    ll_sample_traj = .FALSE.
    IF (nn_sample_rate > 0 .AND. MOD(kt - 1, nn_sample_rate) == 0) ll_sample_traj = .TRUE.
    IF (ll_sample_traj .AND. ASSOCIATED(first_berg)) CALL icb_trj_write(kt)
    CALL iom_put("calving", berg_grid % calving(:, :))
    CALL iom_put("berg_floating_melt", berg_grid % floating_melt(:, :))
    CALL iom_put("berg_stored_ice", berg_grid % stored_ice(:, :, :))
    CALL icb_dia_put
    IF (nn_verbose_level >= 2) CALL icb_utl_print('icb_stp, status', kt)
    ll_budget = .FALSE.
    IF (nn_verbose_write > 0 .AND. MOD(kt - 1, nn_verbose_write) == 0) ll_budget = ln_bergdia
    CALL icb_dia(ll_budget)
    IF (lrst_oce) THEN
      CALL icb_rst_write(kt)
      IF (nn_sample_rate > 0) CALL icb_trj_sync
    END IF
    IF (ln_timing) CALL timing_stop('icb_stp')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE icb_stp
  SUBROUTINE icb_end(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_end', 'r0', 0, 0)
    IF (nn_sample_rate > 0) CALL icb_trj_end
    IF (lwp) WRITE(numout, '(a,i6)') 'icebergs: icb_end complete', narea
    IF (nn_verbose_level > 0) THEN
      CALL flush(numicb)
      CLOSE(UNIT = numicb)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_end
END MODULE icbstp