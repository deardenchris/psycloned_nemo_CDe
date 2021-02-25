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
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    REAL(KIND = wp) :: z2dt, zcoef
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhdiv
    IF (ln_timing) CALL timing_start('ssh_nxt')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'ssh_nxt : after sea surface height'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
    END IF
    z2dt = 2._wp * rdt
    IF (neuler == 0 .AND. kt == nit000) z2dt = rdt
    zcoef = 0.5_wp * r1_rau0
    IF (ln_wd_il) THEN
      CALL wad_lmt(sshb, zcoef * (emp_b(:, :) + emp(:, :)), z2dt)
    END IF
    CALL div_hor(kt)
    zhdiv(:, :) = 0._wp
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      zhdiv(:, :) = zhdiv(:, :) + e3t_n(:, :, jk) * hdivn(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    ssha(:, :) = (sshb(:, :) - z2dt * (zcoef * (emp_b(:, :) + emp(:, :)) + zhdiv(:, :))) * ssmask(:, :)
    IF (.NOT. ln_dynspg_ts) THEN
      IF (ln_bdy) THEN
        CALL lbc_lnk('sshwzv', ssha, 'T', 1.)
        CALL bdy_ssh(ssha)
      END IF
    END IF
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = ssha, clinfo1 = ' ssha  - : ', mask1 = tmask)
    IF (ln_timing) CALL timing_stop('ssh_nxt')
  END SUBROUTINE ssh_nxt
  SUBROUTINE wzv(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: z1_2dt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zhdiv
    IF (ln_timing) CALL timing_start('wzv')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'wzv : now vertical velocity '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~ '
      wn(:, :, jpk) = 0._wp
    END IF
    z1_2dt = 1. / (2. * rdt)
    IF (neuler == 0 .AND. kt == nit000) z1_2dt = 1. / rdt
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      ALLOCATE(zhdiv(jpi, jpj, jpk))
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zhdiv(ji, jj, jk) = r1_e1e2t(ji, jj) * (un_td(ji, jj, jk) - un_td(ji - 1, jj, jk) + vn_td(ji, jj, jk) - vn_td(ji, jj - &
&1, jk))
          END DO
        END DO
      END DO
      CALL lbc_lnk('sshwzv', zhdiv, 'T', 1.)
      DO jk = jpkm1, 1, - 1
        wn(:, :, jk) = wn(:, :, jk + 1) - (e3t_n(:, :, jk) * hdivn(:, :, jk) + zhdiv(:, :, jk) + z1_2dt * (e3t_a(:, :, jk) - &
&e3t_b(:, :, jk))) * tmask(:, :, jk)
      END DO
      DEALLOCATE(zhdiv)
    ELSE
      DO jk = jpkm1, 1, - 1
        wn(:, :, jk) = wn(:, :, jk + 1) - (e3t_n(:, :, jk) * hdivn(:, :, jk) + z1_2dt * (e3t_a(:, :, jk) - e3t_b(:, :, jk))) * &
&tmask(:, :, jk)
      END DO
    END IF
    IF (ln_bdy) THEN
      DO jk = 1, jpkm1
        wn(:, :, jk) = wn(:, :, jk) * bdytmask(:, :)
      END DO
    END IF
    IF (ln_timing) CALL timing_stop('wzv')
  END SUBROUTINE wzv
  SUBROUTINE ssh_swp(kt)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp) :: zcoef
    IF (ln_timing) CALL timing_start('ssh_swp')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'ssh_swp : Asselin time filter and swap of sea surface height'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
    END IF
    IF (neuler == 0 .AND. kt == nit000) THEN
      sshn(:, :) = ssha(:, :)
    ELSE
      sshb(:, :) = sshn(:, :) + atfp * (sshb(:, :) - 2 * sshn(:, :) + ssha(:, :))
      IF (.NOT. ln_linssh) THEN
        zcoef = atfp * rdt * r1_rau0
        sshb(:, :) = sshb(:, :) - zcoef * (emp_b(:, :) - emp(:, :) - rnf_b(:, :) + rnf(:, :) + fwfisf_b(:, :) - fwfisf(:, :)) * &
&ssmask(:, :)
      END IF
      sshn(:, :) = ssha(:, :)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = sshb, clinfo1 = ' sshb  - : ', mask1 = tmask)
    IF (ln_timing) CALL timing_stop('ssh_swp')
  END SUBROUTINE ssh_swp
  SUBROUTINE wAimp(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zCu, zcff, z1_e3w
    REAL(KIND = wp), PARAMETER :: Cu_min = 0.15_wp
    REAL(KIND = wp), PARAMETER :: Cu_max = 0.27
    REAL(KIND = wp), PARAMETER :: Cu_cut = 2._wp * Cu_max - Cu_min
    REAL(KIND = wp), PARAMETER :: Fcu = 4._wp * Cu_max * (Cu_max - Cu_min)
    IF (ln_timing) CALL timing_start('wAimp')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'wAimp : Courant number-based partitioning of now vertical velocity '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~ '
    END IF
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          z1_e3w = 1._wp / e3w_n(ji, jj, jk)
          Cu_adv(ji, jj, jk) = 2._wp * rdt * ((MAX(wn(ji, jj, jk), 0._wp) - MIN(wn(ji, jj, jk + 1), 0._wp)) + (MAX(e2u(ji, jj) * &
&e3uw_n(ji, jj, jk) * un(ji, jj, jk), 0._wp) - MIN(e2u(ji - 1, jj) * e3uw_n(ji - 1, jj, jk) * un(ji - 1, jj, jk), 0._wp)) * &
&r1_e1e2t(ji, jj) + (MAX(e1v(ji, jj) * e3vw_n(ji, jj, jk) * vn(ji, jj, jk), 0._wp) - MIN(e1v(ji, jj - 1) * e3vw_n(ji, jj - 1, jk) &
&* vn(ji, jj - 1, jk), 0._wp)) * r1_e1e2t(ji, jj)) * z1_e3w
        END DO
      END DO
    END DO
    CALL lbc_lnk('sshwzv', Cu_adv, 'T', 1.)
    CALL iom_put("Courant", Cu_adv)
    IF (MAXVAL(Cu_adv(:, :, :)) > Cu_min) THEN
      DO jk = 1, jpkm1
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
    CALL iom_put("wimp", wi)
    CALL iom_put("wi_cff", Cu_adv)
    CALL iom_put("wexp", wn)
    IF (ln_timing) CALL timing_stop('wAimp')
  END SUBROUTINE wAimp
END MODULE sshwzv