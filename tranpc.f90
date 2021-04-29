MODULE tranpc
  USE oce
  USE dom_oce
  USE phycst
  USE zdf_oce
  USE trd_oce
  USE trdtra
  USE eosbn2
  USE lbclnk
  USE in_out_manager
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_npc
  CONTAINS
  SUBROUTINE tra_npc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: inpcc
    INTEGER :: jiter, ikbot, ikp, ikup, ikdown, ilayer, ik_low
    LOGICAL :: l_bottom_reached, l_column_treated
    REAL(KIND = wp) :: zta, zalfa, zsum_temp, zsum_alfa, zaw, zdz, zsum_z
    REAL(KIND = wp) :: zsa, zbeta, zsum_sali, zsum_beta, zbw, zrw, z1_r2dt
    REAL(KIND = wp), PARAMETER :: zn2_zero = 1.E-14_wp
    REAL(KIND = wp), DIMENSION(jpk) :: zvn2
    REAL(KIND = wp), DIMENSION(jpk, jpts) :: zvts, zvab
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zn2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts) :: zab
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdt, ztrds
    LOGICAL, PARAMETER :: l_LB_debug = .FALSE.
    INTEGER :: ilc1, jlc1, klc1, nncpu
    LOGICAL :: lp_monitor_point = .FALSE.
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    CALL profile_psy_data0 % PreStart('tra_npc', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('tra_npc')
    CALL profile_psy_data0 % PostEnd
    IF (MOD(kt, nn_npc) == 0) THEN
      IF (l_trdtra) THEN
        ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
        !$ACC KERNELS
        ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
        ztrds(:, :, :) = tsa(:, :, :, jp_sal)
        !$ACC END KERNELS
      END IF
      CALL profile_psy_data1 % PreStart('tra_npc', 'r1', 0, 0)
      IF (l_LB_debug) THEN
        ilc1 = 45
        jlc1 = 3
        nncpu = 1
        klc1 = mbkt(ilc1, jlc1)
      END IF
      CALL eos_rab(tsa, zab)
      CALL bn2(tsa, zab, zn2)
      inpcc = 0
      CALL profile_psy_data1 % PostEnd
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          IF (tmask(ji, jj, 2) == 1) THEN
            !$ACC KERNELS
            zvts(:, jp_tem) = tsa(ji, jj, :, jp_tem)
            zvts(:, jp_sal) = tsa(ji, jj, :, jp_sal)
            zvab(:, jp_tem) = zab(ji, jj, :, jp_tem)
            zvab(:, jp_sal) = zab(ji, jj, :, jp_sal)
            zvn2(:) = zn2(ji, jj, :)
            !$ACC END KERNELS
            CALL profile_psy_data2 % PreStart('tra_npc', 'r2', 0, 0)
            IF (l_LB_debug) THEN
              lp_monitor_point = .FALSE.
              IF ((ji == ilc1) .AND. (jj == jlc1)) lp_monitor_point = .TRUE.
              lp_monitor_point = (narea == nncpu) .AND. lp_monitor_point
            END IF
            ikbot = mbkt(ji, jj)
            ikp = 1
            ilayer = 0
            jiter = 0
            l_column_treated = .FALSE.
            DO WHILE (.NOT. l_column_treated)
              jiter = jiter + 1
              IF (jiter >= 400) EXIT
              l_bottom_reached = .FALSE.
              DO WHILE (.NOT. l_bottom_reached)
                ikp = ikp + 1
                IF (zvn2(ikp) < - zn2_zero) THEN
                  ilayer = ilayer + 1
                  IF (lp_monitor_point) THEN
                    WRITE(numout, *)
                    IF (ilayer == 1 .AND. jiter == 1) THEN
                      WRITE(numout, *)
                      WRITE(numout, *) 'Time step = ', kt, ' !!!'
                    END IF
                    WRITE(numout, *) ' * Iteration #', jiter, ': found instable portion #', ilayer, ' in column! Starting at ikp =', ikp
                    WRITE(numout, *) ' *** N2 for point (i,j) = ', ji, ' , ', jj
                    DO jk = 1, klc1
                      WRITE(numout, *) jk, zvn2(jk)
                    END DO
                    WRITE(numout, *)
                  END IF
                  IF (jiter == 1) inpcc = inpcc + 1
                  IF (lp_monitor_point) WRITE(numout, *) 'Negative N2 at ikp =', ikp, ' for layer #', ilayer
                  ikup = ikp - 1
                  IF (ikp > 2) THEN
                    DO jk = ikp - 1, 2, - 1
                      IF (ABS(zvn2(jk)) < zn2_zero) THEN
                        ikup = ikup - 1
                      ELSE
                        EXIT
                      END IF
                    END DO
                  END IF
                  IF (ikup < 1) CALL ctl_stop('tra_npc :  PROBLEM #1')
                  zsum_temp = 0._wp
                  zsum_sali = 0._wp
                  zsum_alfa = 0._wp
                  zsum_beta = 0._wp
                  zsum_z = 0._wp
                  DO jk = ikup, ikbot
                    zdz = e3t_n(ji, jj, jk)
                    zsum_temp = zsum_temp + zvts(jk, jp_tem) * zdz
                    zsum_sali = zsum_sali + zvts(jk, jp_sal) * zdz
                    zsum_alfa = zsum_alfa + zvab(jk, jp_tem) * zdz
                    zsum_beta = zsum_beta + zvab(jk, jp_sal) * zdz
                    zsum_z = zsum_z + zdz
                    IF (jk == ikbot) EXIT
                    IF (zvn2(jk + 1) > zn2_zero) EXIT
                  END DO
                  ikdown = jk
                  IF (ikup == ikdown) CALL ctl_stop('tra_npc :  PROBLEM #2')
                  zta = zsum_temp / zsum_z
                  zsa = zsum_sali / zsum_z
                  zalfa = zsum_alfa / zsum_z
                  zbeta = zsum_beta / zsum_z
                  IF (lp_monitor_point) THEN
                    WRITE(numout, *) 'MIXED T, S, alfa and beta between ikup =', ikup, ' and ikdown =', ikdown, ', in layer #', ilayer
                    WRITE(numout, *) '  => Mean temp. in that portion =', zta
                    WRITE(numout, *) '  => Mean sali. in that portion =', zsa
                    WRITE(numout, *) '  => Mean Alfa  in that portion =', zalfa
                    WRITE(numout, *) '  => Mean Beta  in that portion =', zbeta
                  END IF
                  DO jk = ikup, ikdown
                    zvts(jk, jp_tem) = zta
                    zvts(jk, jp_sal) = zsa
                    zvab(jk, jp_tem) = zalfa
                    zvab(jk, jp_sal) = zbeta
                  END DO
                  ikup = MAX(2, ikup)
                  ik_low = MIN(ikdown + 1, ikbot)
                  DO jk = ikup, ik_low
                    zrw = (gdepw_n(ji, jj, jk) - gdept_n(ji, jj, jk)) / (gdept_n(ji, jj, jk - 1) - gdept_n(ji, jj, jk))
                    zaw = zvab(jk, jp_tem) * (1._wp - zrw) + zvab(jk - 1, jp_tem) * zrw
                    zbw = zvab(jk, jp_sal) * (1._wp - zrw) + zvab(jk - 1, jp_sal) * zrw
                    zvn2(jk) = grav * (zaw * (zvts(jk - 1, jp_tem) - zvts(jk, jp_tem)) - zbw * (zvts(jk - 1, jp_sal) - zvts(jk, jp_sal))) / e3w_n(ji, jj, jk) * tmask(ji, jj, jk)
                  END DO
                  ikp = MIN(ikdown + 1, ikbot)
                END IF
                IF (ikp == ikbot) l_bottom_reached = .TRUE.
              END DO
              IF (ikp /= ikbot) CALL ctl_stop('tra_npc :  PROBLEM #3')
              IF (ilayer > 0) THEN
                IF (lp_monitor_point) THEN
                  WRITE(numout, *)
                  WRITE(numout, *) 'After ', jiter, ' iteration(s), we neutralized ', ilayer, ' instable layer(s)'
                  WRITE(numout, *) '   ==> N2 at i,j=', ji, ',', jj, ' now looks like this:'
                  DO jk = 1, klc1
                    WRITE(numout, *) jk, zvn2(jk)
                  END DO
                  WRITE(numout, *)
                END IF
                ikp = 1
                ilayer = 0
              END IF
              IF (ikp >= ikbot) l_column_treated = .TRUE.
            END DO
            CALL profile_psy_data2 % PostEnd
            !$ACC KERNELS
            tsa(ji, jj, :, jp_tem) = zvts(:, jp_tem)
            tsa(ji, jj, :, jp_sal) = zvts(:, jp_sal)
            !$ACC END KERNELS
            IF (lp_monitor_point) WRITE(numout, *)
          END IF
        END DO
      END DO
      IF (l_trdtra) THEN
        !$ACC KERNELS
        z1_r2dt = 1._wp / (2._wp * rdt)
        ztrdt(:, :, :) = (tsa(:, :, :, jp_tem) - ztrdt(:, :, :)) * z1_r2dt
        ztrds(:, :, :) = (tsa(:, :, :, jp_sal) - ztrds(:, :, :)) * z1_r2dt
        !$ACC END KERNELS
        CALL profile_psy_data3 % PreStart('tra_npc', 'r3', 0, 0)
        CALL trd_tra(kt, 'TRA', jp_tem, jptra_npc, ztrdt)
        CALL trd_tra(kt, 'TRA', jp_sal, jptra_npc, ztrds)
        DEALLOCATE(ztrdt, ztrds)
        CALL profile_psy_data3 % PostEnd
      END IF
      CALL profile_psy_data4 % PreStart('tra_npc', 'r4', 0, 0)
      CALL lbc_lnk_multi('tranpc', tsa(:, :, :, jp_tem), 'T', 1., tsa(:, :, :, jp_sal), 'T', 1.)
      IF (lwp .AND. l_LB_debug) THEN
        WRITE(numout, *) 'Exiting tra_npc , kt = ', kt, ', => numb. of statically instable water-columns: ', inpcc
        WRITE(numout, *)
      END IF
      CALL profile_psy_data4 % PostEnd
    END IF
    CALL profile_psy_data5 % PreStart('tra_npc', 'r5', 0, 0)
    IF (ln_timing) CALL timing_stop('tra_npc')
    CALL profile_psy_data5 % PostEnd
  END SUBROUTINE tra_npc
END MODULE tranpc