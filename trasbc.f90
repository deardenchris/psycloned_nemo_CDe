MODULE trasbc
  USE oce
  USE sbc_oce
  USE dom_oce
  USE phycst
  USE eosbn2
  USE sbcmod
  USE sbcrnf
  USE sbcisf
  USE iscplini
  USE traqsr
  USE trd_oce
  USE trdtra
  USE in_out_manager
  USE prtctl
  USE iom
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_sbc
  CONTAINS
  SUBROUTINE tra_sbc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jn
    INTEGER :: ikt, ikb
    REAL(KIND = wp) :: zfact, z1_e3t, zdep, ztim
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdt, ztrds
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    CALL profile_psy_data0 % PreStart('tra_sbc', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('tra_sbc')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_sbc : TRAcer Surface Boundary Condition'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal)
      !$ACC END KERNELS
    END IF
    IF (.NOT. ln_traqsr) THEN
      !$ACC KERNELS
      qns(:, :) = qns(:, :) + qsr(:, :)
      qsr(:, :) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (kt == nit000) THEN
      IF (ln_rstart .AND. iom_varid(numror, 'sbc_hc_b', ldstop = .FALSE.) > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) '          nit000-1 sbc tracer content field read in the restart file'
        !$ACC KERNELS
        zfact = 0.5_wp
        sbc_tsc(:, :, :) = 0._wp
        !$ACC END KERNELS
        CALL profile_psy_data1 % PreStart('tra_sbc', 'r1', 0, 0)
        CALL iom_get(numror, jpdom_autoglo, 'sbc_hc_b', sbc_tsc_b(:, :, jp_tem), ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'sbc_sc_b', sbc_tsc_b(:, :, jp_sal), ldxios = lrxios)
        CALL profile_psy_data1 % PostEnd
      ELSE
        !$ACC KERNELS
        zfact = 1._wp
        sbc_tsc(:, :, :) = 0._wp
        sbc_tsc_b(:, :, :) = 0._wp
        !$ACC END KERNELS
      END IF
    ELSE
      !$ACC KERNELS
      zfact = 0.5_wp
      sbc_tsc_b(:, :, :) = sbc_tsc(:, :, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpj
      DO ji = 2, jpim1
        sbc_tsc(ji, jj, jp_tem) = r1_rau0_rcp * qns(ji, jj)
        sbc_tsc(ji, jj, jp_sal) = r1_rau0 * sfx(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_linssh) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpj
        DO ji = 2, jpim1
          sbc_tsc(ji, jj, jp_tem) = sbc_tsc(ji, jj, jp_tem) + r1_rau0 * emp(ji, jj) * tsn(ji, jj, 1, jp_tem)
          sbc_tsc(ji, jj, jp_sal) = sbc_tsc(ji, jj, jp_sal) + r1_rau0 * emp(ji, jj) * tsn(ji, jj, 1, jp_sal)
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data2 % PreStart('tra_sbc', 'r2', 0, 0)
      IF (iom_use('emp_x_sst')) CALL iom_put("emp_x_sst", emp(:, :) * tsn(:, :, 1, jp_tem))
      IF (iom_use('emp_x_sss')) CALL iom_put("emp_x_sss", emp(:, :) * tsn(:, :, 1, jp_sal))
      CALL profile_psy_data2 % PostEnd
    END IF
    !$ACC KERNELS
    DO jn = 1, jpts
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpj
        DO ji = 2, jpim1
          tsa(ji, jj, 1, jn) = tsa(ji, jj, 1, jn) + zfact * (sbc_tsc_b(ji, jj, jn) + sbc_tsc(ji, jj, jn)) / e3t_n(ji, jj, 1)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('tra_sbc', 'r3', 0, 0)
    IF (lrst_oce) THEN
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'sbc_hc_b', sbc_tsc(:, :, jp_tem), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sbc_sc_b', sbc_tsc(:, :, jp_sal), ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
    CALL profile_psy_data3 % PostEnd
    IF (ln_isf) THEN
      !$ACC KERNELS
      zfact = 0.5_wp
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpj
        DO ji = 2, jpim1
          ikt = misfkt(ji, jj)
          ikb = misfkb(ji, jj)
          DO jk = ikt, ikb - 1
            tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) + zfact * (risf_tsc_b(ji, jj, jp_tem) + risf_tsc(ji, jj, jp_tem)) * &
&r1_hisf_tbl(ji, jj)
          END DO
          tsa(ji, jj, ikb, jp_tem) = tsa(ji, jj, ikb, jp_tem) + zfact * (risf_tsc_b(ji, jj, jp_tem) + risf_tsc(ji, jj, jp_tem)) * &
&r1_hisf_tbl(ji, jj) * ralpha(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data4 % PreStart('tra_sbc', 'r4', 0, 0)
    IF (ln_rnf) THEN
      zfact = 0.5_wp
      DO jj = 2, jpj
        DO ji = 2, jpim1
          IF (rnf(ji, jj) /= 0._wp) THEN
            zdep = zfact / h_rnf(ji, jj)
            DO jk = 1, nk_rnf(ji, jj)
              tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) + (rnf_tsc_b(ji, jj, jp_tem) + rnf_tsc(ji, jj, jp_tem)) * zdep
              IF (ln_rnf_sal) tsa(ji, jj, jk, jp_sal) = tsa(ji, jj, jk, jp_sal) + (rnf_tsc_b(ji, jj, jp_sal) + rnf_tsc(ji, jj, &
&jp_sal)) * zdep
            END DO
          END IF
        END DO
      END DO
    END IF
    IF (iom_use('rnf_x_sst')) CALL iom_put("rnf_x_sst", rnf * tsn(:, :, 1, jp_tem))
    IF (iom_use('rnf_x_sss')) CALL iom_put("rnf_x_sss", rnf * tsn(:, :, 1, jp_sal))
    CALL profile_psy_data4 % PostEnd
    !$ACC KERNELS
    IF (ln_iscpl .AND. ln_hsb) THEN
      DO jk = 1, jpk
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpj
          DO ji = 2, jpim1
            zdep = 1._wp / e3t_n(ji, jj, jk)
            tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) - htsc_iscpl(ji, jj, jk, jp_tem) * zdep
            tsa(ji, jj, jk, jp_sal) = tsa(ji, jj, jk, jp_sal) - htsc_iscpl(ji, jj, jk, jp_sal) * zdep
          END DO
        END DO
      END DO
    END IF
    !$ACC END KERNELS
    IF (l_trdtra) THEN
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem) - ztrdt(:, :, :)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal) - ztrds(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data5 % PreStart('tra_sbc', 'r5', 0, 0)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_nsr, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_nsr, ztrds)
      DEALLOCATE(ztrdt, ztrds)
      CALL profile_psy_data5 % PostEnd
    END IF
    CALL profile_psy_data6 % PreStart('tra_sbc', 'r6', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' sbc  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, &
&jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
    IF (ln_timing) CALL timing_stop('tra_sbc')
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE tra_sbc
END MODULE trasbc