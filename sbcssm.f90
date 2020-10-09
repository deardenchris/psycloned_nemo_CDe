MODULE sbcssm
  USE oce
  USE dom_oce
  USE sbc_oce
  USE sbcapr
  USE eosbn2
  USE traqsr, ONLY: ln_traqsr
  USE in_out_manager
  USE prtctl
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_ssm
  PUBLIC :: sbc_ssm_init
  LOGICAL, SAVE :: l_ssm_mean = .FALSE.
  CONTAINS
  SUBROUTINE sbc_ssm(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zcoef, zf_sbc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts) :: zts
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zts(ji, jj, jp_tem) = tsn(ji, jj, mikt(ji, jj), jp_tem)
        zts(ji, jj, jp_sal) = tsn(ji, jj, mikt(ji, jj), jp_sal)
      END DO
    END DO
    !$ACC END KERNELS
    IF (nn_fsbc == 1) THEN
      !$ACC KERNELS
      ssu_m(:, :) = ub(:, :, 1)
      ssv_m(:, :) = vb(:, :, 1)
      !$ACC END KERNELS
      IF (l_usect) THEN
        CALL profile_psy_data0 % PreStart('sbc_ssm', 'r0', 0, 0)
        sst_m(:, :) = eos_pt_from_ct(zts(:, :, jp_tem), zts(:, :, jp_sal))
        CALL profile_psy_data0 % PostEnd
      ELSE
        !$ACC KERNELS
        sst_m(:, :) = zts(:, :, jp_tem)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      sss_m(:, :) = zts(:, :, jp_sal)
      !$ACC END KERNELS
      IF (ln_apr_dyn) THEN
        !$ACC KERNELS
        ssh_m(:, :) = sshn(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :))
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        ssh_m(:, :) = sshn(:, :)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      e3t_m(:, :) = e3t_n(:, :, 1)
      frq_m(:, :) = fraqsr_1lev(:, :)
      !$ACC END KERNELS
    ELSE
      IF (kt == nit000 .AND. .NOT. l_ssm_mean) THEN
        CALL profile_psy_data1 % PreStart('sbc_ssm', 'r1', 0, 0)
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'sbc_ssm : mean fields initialised to instantaneous values'
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~   '
        CALL profile_psy_data1 % PostEnd
        !$ACC KERNELS
        zcoef = REAL(nn_fsbc - 1, wp)
        ssu_m(:, :) = zcoef * ub(:, :, 1)
        ssv_m(:, :) = zcoef * vb(:, :, 1)
        !$ACC END KERNELS
        IF (l_usect) THEN
          CALL profile_psy_data2 % PreStart('sbc_ssm', 'r2', 0, 0)
          sst_m(:, :) = zcoef * eos_pt_from_ct(zts(:, :, jp_tem), zts(:, :, jp_sal))
          CALL profile_psy_data2 % PostEnd
        ELSE
          !$ACC KERNELS
          sst_m(:, :) = zcoef * zts(:, :, jp_tem)
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        sss_m(:, :) = zcoef * zts(:, :, jp_sal)
        !$ACC END KERNELS
        IF (ln_apr_dyn) THEN
          !$ACC KERNELS
          ssh_m(:, :) = zcoef * (sshn(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :)))
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          ssh_m(:, :) = zcoef * sshn(:, :)
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        e3t_m(:, :) = zcoef * e3t_n(:, :, 1)
        frq_m(:, :) = zcoef * fraqsr_1lev(:, :)
        !$ACC END KERNELS
      ELSE IF (MOD(kt - 2, nn_fsbc) == 0) THEN
        !$ACC KERNELS
        ssu_m(:, :) = 0._wp
        ssv_m(:, :) = 0._wp
        sst_m(:, :) = 0._wp
        sss_m(:, :) = 0._wp
        ssh_m(:, :) = 0._wp
        e3t_m(:, :) = 0._wp
        frq_m(:, :) = 0._wp
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      ssu_m(:, :) = ssu_m(:, :) + ub(:, :, 1)
      ssv_m(:, :) = ssv_m(:, :) + vb(:, :, 1)
      !$ACC END KERNELS
      IF (l_usect) THEN
        CALL profile_psy_data3 % PreStart('sbc_ssm', 'r3', 0, 0)
        sst_m(:, :) = sst_m(:, :) + eos_pt_from_ct(zts(:, :, jp_tem), zts(:, :, jp_sal))
        CALL profile_psy_data3 % PostEnd
      ELSE
        !$ACC KERNELS
        sst_m(:, :) = sst_m(:, :) + zts(:, :, jp_tem)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      sss_m(:, :) = sss_m(:, :) + zts(:, :, jp_sal)
      !$ACC END KERNELS
      IF (ln_apr_dyn) THEN
        !$ACC KERNELS
        ssh_m(:, :) = ssh_m(:, :) + sshn(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :))
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        ssh_m(:, :) = ssh_m(:, :) + sshn(:, :)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      e3t_m(:, :) = e3t_m(:, :) + e3t_n(:, :, 1)
      frq_m(:, :) = frq_m(:, :) + fraqsr_1lev(:, :)
      IF (MOD(kt - 1, nn_fsbc) == 0) THEN
        zcoef = 1. / REAL(nn_fsbc, wp)
        sst_m(:, :) = sst_m(:, :) * zcoef
        sss_m(:, :) = sss_m(:, :) * zcoef
        ssu_m(:, :) = ssu_m(:, :) * zcoef
        ssv_m(:, :) = ssv_m(:, :) * zcoef
        ssh_m(:, :) = ssh_m(:, :) * zcoef
        e3t_m(:, :) = e3t_m(:, :) * zcoef
        frq_m(:, :) = frq_m(:, :) * zcoef
      END IF
      !$ACC END KERNELS
      CALL profile_psy_data4 % PreStart('sbc_ssm', 'r4', 0, 0)
      IF (lrst_oce) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'sbc_ssm : sea surface mean fields written in ocean restart file ', 'at it= ', kt, ' date= &
&', ndastp
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
        zf_sbc = REAL(nn_fsbc, wp)
        IF (lwxios) CALL iom_swap(cwxios_context)
        CALL iom_rstput(kt, nitrst, numrow, 'nn_fsbc', zf_sbc, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'ssu_m', ssu_m, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'ssv_m', ssv_m, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'sst_m', sst_m, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'sss_m', sss_m, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'ssh_m', ssh_m, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'e3t_m', e3t_m, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'frq_m', frq_m, ldxios = lwxios)
        IF (lwxios) CALL iom_swap(cxios_context)
      END IF
      CALL profile_psy_data4 % PostEnd
    END IF
    CALL profile_psy_data5 % PreStart('sbc_ssm', 'r5', 0, 0)
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      CALL iom_put('ssu_m', ssu_m)
      CALL iom_put('ssv_m', ssv_m)
      CALL iom_put('sst_m', sst_m)
      CALL iom_put('sss_m', sss_m)
      CALL iom_put('ssh_m', ssh_m)
      CALL iom_put('e3t_m', e3t_m)
      CALL iom_put('frq_m', frq_m)
    END IF
    CALL profile_psy_data5 % PostEnd
  END SUBROUTINE sbc_ssm
  SUBROUTINE sbc_ssm_init
    REAL(KIND = wp) :: zcoef, zf_sbc
    IF (nn_fsbc == 1) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc_ssm_init : sea surface mean fields, nn_fsbc=1 : instantaneous values'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~ '
    ELSE
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc_ssm_init : sea surface mean fields'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~ '
      IF (ln_rstart .AND. iom_varid(numror, 'nn_fsbc', ldstop = .FALSE.) > 0) THEN
        l_ssm_mean = .TRUE.
        CALL iom_get(numror, 'nn_fsbc', zf_sbc, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'ssu_m', ssu_m, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'ssv_m', ssv_m, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'sst_m', sst_m, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'sss_m', sss_m, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'ssh_m', ssh_m, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'e3t_m', e3t_m, ldxios = lrxios)
        IF (iom_varid(numror, 'frq_m', ldstop = .FALSE.) > 0) THEN
          CALL iom_get(numror, jpdom_autoglo, 'frq_m', frq_m, ldxios = lrxios)
        ELSE
          !$ACC KERNELS
          frq_m(:, :) = 1._wp
          !$ACC END KERNELS
        END IF
        IF (zf_sbc /= REAL(nn_fsbc, wp)) THEN
          IF (lwp) WRITE(numout, FMT = *) '   restart with a change in the frequency of mean from ', zf_sbc, ' to ', nn_fsbc
          !$ACC KERNELS
          zcoef = REAL(nn_fsbc - 1, wp) / zf_sbc
          ssu_m(:, :) = zcoef * ssu_m(:, :)
          ssv_m(:, :) = zcoef * ssv_m(:, :)
          sst_m(:, :) = zcoef * sst_m(:, :)
          sss_m(:, :) = zcoef * sss_m(:, :)
          ssh_m(:, :) = zcoef * ssh_m(:, :)
          e3t_m(:, :) = zcoef * e3t_m(:, :)
          frq_m(:, :) = zcoef * frq_m(:, :)
          !$ACC END KERNELS
        ELSE
          IF (lwp) WRITE(numout, FMT = *) '   mean fields read in the ocean restart file'
        END IF
      END IF
    END IF
    IF (.NOT. l_ssm_mean) THEN
      IF (lwp) WRITE(numout, FMT = *) '   default initialisation of ss._m arrays'
      !$ACC KERNELS
      ssu_m(:, :) = ub(:, :, 1)
      ssv_m(:, :) = vb(:, :, 1)
      !$ACC END KERNELS
      IF (l_usect) THEN
        sst_m(:, :) = eos_pt_from_ct(tsn(:, :, 1, jp_tem), tsn(:, :, 1, jp_sal))
      ELSE
        !$ACC KERNELS
        sst_m(:, :) = tsn(:, :, 1, jp_tem)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      sss_m(:, :) = tsn(:, :, 1, jp_sal)
      ssh_m(:, :) = sshn(:, :)
      e3t_m(:, :) = e3t_n(:, :, 1)
      frq_m(:, :) = 1._wp
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    IF (.NOT. ln_traqsr) fraqsr_1lev(:, :) = 1._wp
    !$ACC END KERNELS
    IF (lwxios .AND. nn_fsbc > 1) THEN
      CALL iom_set_rstw_var_active('nn_fsbc')
      CALL iom_set_rstw_var_active('ssu_m')
      CALL iom_set_rstw_var_active('ssv_m')
      CALL iom_set_rstw_var_active('sst_m')
      CALL iom_set_rstw_var_active('sss_m')
      CALL iom_set_rstw_var_active('ssh_m')
      CALL iom_set_rstw_var_active('e3t_m')
      CALL iom_set_rstw_var_active('frq_m')
    END IF
  END SUBROUTINE sbc_ssm_init
END MODULE sbcssm