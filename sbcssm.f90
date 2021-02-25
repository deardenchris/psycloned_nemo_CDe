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
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zcoef, zf_sbc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts) :: zts
    DO jj = 1, jpj
      DO ji = 1, jpi
        zts(ji, jj, jp_tem) = tsn(ji, jj, mikt(ji, jj), jp_tem)
        zts(ji, jj, jp_sal) = tsn(ji, jj, mikt(ji, jj), jp_sal)
      END DO
    END DO
    IF (nn_fsbc == 1) THEN
      ssu_m(:, :) = ub(:, :, 1)
      ssv_m(:, :) = vb(:, :, 1)
      IF (l_usect) THEN
        sst_m(:, :) = eos_pt_from_ct(zts(:, :, jp_tem), zts(:, :, jp_sal))
      ELSE
        sst_m(:, :) = zts(:, :, jp_tem)
      END IF
      sss_m(:, :) = zts(:, :, jp_sal)
      IF (ln_apr_dyn) THEN
        ssh_m(:, :) = sshn(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :))
      ELSE
        ssh_m(:, :) = sshn(:, :)
      END IF
      e3t_m(:, :) = e3t_n(:, :, 1)
      frq_m(:, :) = fraqsr_1lev(:, :)
    ELSE
      IF (kt == nit000 .AND. .NOT. l_ssm_mean) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'sbc_ssm : mean fields initialised to instantaneous values'
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~   '
        zcoef = REAL(nn_fsbc - 1, wp)
        ssu_m(:, :) = zcoef * ub(:, :, 1)
        ssv_m(:, :) = zcoef * vb(:, :, 1)
        IF (l_usect) THEN
          sst_m(:, :) = zcoef * eos_pt_from_ct(zts(:, :, jp_tem), zts(:, :, jp_sal))
        ELSE
          sst_m(:, :) = zcoef * zts(:, :, jp_tem)
        END IF
        sss_m(:, :) = zcoef * zts(:, :, jp_sal)
        IF (ln_apr_dyn) THEN
          ssh_m(:, :) = zcoef * (sshn(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :)))
        ELSE
          ssh_m(:, :) = zcoef * sshn(:, :)
        END IF
        e3t_m(:, :) = zcoef * e3t_n(:, :, 1)
        frq_m(:, :) = zcoef * fraqsr_1lev(:, :)
      ELSE IF (MOD(kt - 2, nn_fsbc) == 0) THEN
        ssu_m(:, :) = 0._wp
        ssv_m(:, :) = 0._wp
        sst_m(:, :) = 0._wp
        sss_m(:, :) = 0._wp
        ssh_m(:, :) = 0._wp
        e3t_m(:, :) = 0._wp
        frq_m(:, :) = 0._wp
      END IF
      ssu_m(:, :) = ssu_m(:, :) + ub(:, :, 1)
      ssv_m(:, :) = ssv_m(:, :) + vb(:, :, 1)
      IF (l_usect) THEN
        sst_m(:, :) = sst_m(:, :) + eos_pt_from_ct(zts(:, :, jp_tem), zts(:, :, jp_sal))
      ELSE
        sst_m(:, :) = sst_m(:, :) + zts(:, :, jp_tem)
      END IF
      sss_m(:, :) = sss_m(:, :) + zts(:, :, jp_sal)
      IF (ln_apr_dyn) THEN
        ssh_m(:, :) = ssh_m(:, :) + sshn(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :))
      ELSE
        ssh_m(:, :) = ssh_m(:, :) + sshn(:, :)
      END IF
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
    END IF
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      CALL iom_put('ssu_m', ssu_m)
      CALL iom_put('ssv_m', ssv_m)
      CALL iom_put('sst_m', sst_m)
      CALL iom_put('sss_m', sss_m)
      CALL iom_put('ssh_m', ssh_m)
      CALL iom_put('e3t_m', e3t_m)
      CALL iom_put('frq_m', frq_m)
    END IF
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
          frq_m(:, :) = 1._wp
        END IF
        IF (zf_sbc /= REAL(nn_fsbc, wp)) THEN
          IF (lwp) WRITE(numout, FMT = *) '   restart with a change in the frequency of mean from ', zf_sbc, ' to ', nn_fsbc
          zcoef = REAL(nn_fsbc - 1, wp) / zf_sbc
          ssu_m(:, :) = zcoef * ssu_m(:, :)
          ssv_m(:, :) = zcoef * ssv_m(:, :)
          sst_m(:, :) = zcoef * sst_m(:, :)
          sss_m(:, :) = zcoef * sss_m(:, :)
          ssh_m(:, :) = zcoef * ssh_m(:, :)
          e3t_m(:, :) = zcoef * e3t_m(:, :)
          frq_m(:, :) = zcoef * frq_m(:, :)
        ELSE
          IF (lwp) WRITE(numout, FMT = *) '   mean fields read in the ocean restart file'
        END IF
      END IF
    END IF
    IF (.NOT. l_ssm_mean) THEN
      IF (lwp) WRITE(numout, FMT = *) '   default initialisation of ss._m arrays'
      ssu_m(:, :) = ub(:, :, 1)
      ssv_m(:, :) = vb(:, :, 1)
      IF (l_usect) THEN
        sst_m(:, :) = eos_pt_from_ct(tsn(:, :, 1, jp_tem), tsn(:, :, 1, jp_sal))
      ELSE
        sst_m(:, :) = tsn(:, :, 1, jp_tem)
      END IF
      sss_m(:, :) = tsn(:, :, 1, jp_sal)
      ssh_m(:, :) = sshn(:, :)
      e3t_m(:, :) = e3t_n(:, :, 1)
      frq_m(:, :) = 1._wp
    END IF
    IF (.NOT. ln_traqsr) fraqsr_1lev(:, :) = 1._wp
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