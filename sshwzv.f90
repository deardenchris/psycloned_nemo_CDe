MODULE sshwzv
  USE oce
  USE dom_oce
  USE sbc_oce
  USE domvvl
  USE divhor
  USE phycst
  USE bdy_oce, ONLY: ln_bdy, bdytmask
  USE bdydyn2d
  USE iom
  USE in_out_manager
  USE restart
  USE prtctl
  USE lbclnk
  USE lib_mpp
  USE timing
  USE wet_dry
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ssh_nxt
  PUBLIC :: wzv
  PUBLIC :: wAimp
  PUBLIC :: ssh_swp
  CONTAINS
  SUBROUTINE ssh_nxt(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    REAL(KIND = wp) :: z2dt, zcoef
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhdiv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ssh_nxt', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ssh_nxt')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'ssh_nxt : after sea surface height'
      IF (lwp) WRITE(numout, *) '~~~~~~~ '
    END IF
    z2dt = 2._wp * rdt
    IF (neuler == 0 .AND. kt == nit000) z2dt = rdt
    zcoef = 0.5_wp * r1_rau0
    IF (ln_wd_il) THEN
      CALL wad_lmt(sshb, zcoef * (emp_b(:, :) + emp(:, :)), z2dt)
    END IF
    CALL div_hor(kt)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zhdiv(:, :) = 0._wp
    DO jk = 1, jpkm1
      zhdiv(:, :) = zhdiv(:, :) + e3t_n(:, :, jk) * hdivn(:, :, jk)
    END DO
    ssha(:, :) = (sshb(:, :) - z2dt * (zcoef * (emp_b(:, :) + emp(:, :)) + zhdiv(:, :))) * ssmask(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ssh_nxt', 'r1', 0, 0)
    IF (.NOT. ln_dynspg_ts) THEN
      IF (ln_bdy) THEN
        CALL lbc_lnk('sshwzv', ssha, 'T', 1.)
        CALL bdy_ssh(ssha)
      END IF
    END IF
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = ssha, clinfo1 = ' ssha  - : ', mask1 = tmask)
    IF (ln_timing) CALL timing_stop('ssh_nxt')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ssh_nxt
  SUBROUTINE wzv(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: z1_2dt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zhdiv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('wzv', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('wzv')
    CALL profile_psy_data0 % PostEnd
    IF (kt == nit000) THEN
      CALL profile_psy_data1 % PreStart('wzv', 'r1', 0, 0)
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'wzv : now vertical velocity '
      IF (lwp) WRITE(numout, *) '~~~~~ '
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      wn(:, :, jpk) = 0._wp
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data2 % PreStart('wzv', 'r2', 0, 0)
    z1_2dt = 1. / (2. * rdt)
    IF (neuler == 0 .AND. kt == nit000) z1_2dt = 1. / rdt
    CALL profile_psy_data2 % PostEnd
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      ALLOCATE(zhdiv(jpi, jpj, jpk))
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zhdiv(ji, jj, jk) = r1_e1e2t(ji, jj) * (un_td(ji, jj, jk) - un_td(ji - 1, jj, jk) + vn_td(ji, jj, jk) - vn_td(ji, jj - 1, jk))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('wzv', 'r3', 0, 0)
      CALL lbc_lnk('sshwzv', zhdiv, 'T', 1.)
      CALL profile_psy_data3 % PostEnd
      !$ACC KERNELS
      DO jk = jpkm1, 1, - 1
        wn(:, :, jk) = wn(:, :, jk + 1) - (e3t_n(:, :, jk) * hdivn(:, :, jk) + zhdiv(:, :, jk) + z1_2dt * (e3t_a(:, :, jk) - e3t_b(:, :, jk))) * tmask(:, :, jk)
      END DO
      !$ACC END KERNELS
      DEALLOCATE(zhdiv)
    ELSE
      !$ACC KERNELS
      DO jk = jpkm1, 1, - 1
        wn(:, :, jk) = wn(:, :, jk + 1) - (e3t_n(:, :, jk) * hdivn(:, :, jk) + z1_2dt * (e3t_a(:, :, jk) - e3t_b(:, :, jk))) * tmask(:, :, jk)
      END DO
      !$ACC END KERNELS
    END IF
    IF (ln_bdy) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        wn(:, :, jk) = wn(:, :, jk) * bdytmask(:, :)
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data4 % PreStart('wzv', 'r4', 0, 0)
    IF (ln_timing) CALL timing_stop('wzv')
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE wzv
  SUBROUTINE ssh_swp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp) :: zcoef
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ssh_swp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ssh_swp')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'ssh_swp : Asselin time filter and swap of sea surface height'
      IF (lwp) WRITE(numout, *) '~~~~~~~ '
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (neuler == 0 .AND. kt == nit000) THEN
      !$ACC KERNELS
      sshn(:, :) = ssha(:, :)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      sshb(:, :) = sshn(:, :) + atfp * (sshb(:, :) - 2 * sshn(:, :) + ssha(:, :))
      !$ACC END KERNELS
      IF (.NOT. ln_linssh) THEN
        !$ACC KERNELS
        zcoef = atfp * rdt * r1_rau0
        sshb(:, :) = sshb(:, :) - zcoef * (emp_b(:, :) - emp(:, :) - rnf_b(:, :) + rnf(:, :) + fwfisf_b(:, :) - fwfisf(:, :)) * ssmask(:, :)
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      sshn(:, :) = ssha(:, :)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data1 % PreStart('ssh_swp', 'r1', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = sshb, clinfo1 = ' sshb  - : ', mask1 = tmask)
    IF (ln_timing) CALL timing_stop('ssh_swp')
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ssh_swp
  SUBROUTINE wAimp(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zCu, zcff, z1_e3w
    REAL(KIND = wp), PARAMETER :: Cu_min = 0.15_wp
    REAL(KIND = wp), PARAMETER :: Cu_max = 0.27
    REAL(KIND = wp), PARAMETER :: Cu_cut = 2._wp * Cu_max - Cu_min
    REAL(KIND = wp), PARAMETER :: Fcu = 4._wp * Cu_max * (Cu_max - Cu_min)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('waimp', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('wAimp')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'wAimp : Courant number-based partitioning of now vertical velocity '
      IF (lwp) WRITE(numout, *) '~~~~~ '
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO jk = 1, jpkm1
      !$ACC loop independent collapse(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          z1_e3w = 1._wp / e3w_n(ji, jj, jk)
          Cu_adv(ji, jj, jk) = 2._wp * rdt * ((MAX(wn(ji, jj, jk), 0._wp) - MIN(wn(ji, jj, jk + 1), 0._wp)) + (MAX(e2u(ji, jj) * e3uw_n(ji, jj, jk) * un(ji, jj, jk), 0._wp) - MIN(e2u(ji - 1, jj) * e3uw_n(ji - 1, jj, jk) * un(ji - 1, jj, jk), 0._wp)) * r1_e1e2t(ji, jj) + (MAX(e1v(ji, jj) * e3vw_n(ji, jj, jk) * vn(ji, jj, jk), 0._wp) - MIN(e1v(ji, jj - 1) * e3vw_n(ji, jj - 1, jk) * vn(ji, jj - 1, jk), 0._wp)) * r1_e1e2t(ji, jj)) * z1_e3w
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('waimp', 'r1', 0, 0)
    CALL lbc_lnk('sshwzv', Cu_adv, 'T', 1.)
    CALL iom_put("Courant", Cu_adv)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    IF (MAXVAL(Cu_adv(:, :, :)) > Cu_min) THEN
      DO jk = 1, jpkm1
        !$ACC loop independent collapse(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zCu = MAX(Cu_adv(ji, jj, jk), Cu_adv(ji, jj, jk + 1))
            IF (zCu <= Cu_min) THEN
              zcff = 0._wp
            ELSE IF (zCu < Cu_cut) THEN
              zcff = (zCu - Cu_min) ** 2
              zcff = zcff / (Fcu + zcff)
            ELSE
              zcff = (zCu - Cu_max) / zCu
            END IF
            zcff = MIN(1._wp, zcff)
            wi(ji, jj, jk) = zcff * wn(ji, jj, jk)
            wn(ji, jj, jk) = (1._wp - zcff) * wn(ji, jj, jk)
            Cu_adv(ji, jj, jk) = zcff
          END DO
        END DO
      END DO
    ELSE
      Cu_adv(:, :, :) = 0._wp
      wi(:, :, :) = 0._wp
    END IF
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('waimp', 'r2', 0, 0)
    CALL iom_put("wimp", wi)
    CALL iom_put("wi_cff", Cu_adv)
    CALL iom_put("wexp", wn)
    IF (ln_timing) CALL timing_stop('wAimp')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE wAimp
END MODULE sshwzv