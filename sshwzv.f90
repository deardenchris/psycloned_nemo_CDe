MODULE sshwzv
  USE oce
  USE dom_oce
  USE sbc_oce
  USE domvvl
  USE divhor
  USE phycst
  USE bdy_oce, ONLY: ln_bdy, bdytmask
  USE bdydyn2d
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
    !$ACC KERNELS
    zhdiv(:, :) = 0._wp
    DO jk = 1, jpkm1
      zhdiv(:, :) = zhdiv(:, :) + e3t_n(:, :, jk) * hdivn(:, :, jk)
    END DO
    ssha(:, :) = (sshb(:, :) - z2dt * (zcoef * (emp_b(:, :) + emp(:, :)) + zhdiv(:, :))) * ssmask(:, :)
    !$ACC END KERNELS
    IF (.NOT. ln_dynspg_ts) THEN
      IF (ln_bdy) THEN
        CALL lbc_lnk(ssha, 'T', 1.)
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
      !$ACC KERNELS
      wn(:, :, jpk) = 0._wp
      !$ACC END KERNELS
    END IF
    z1_2dt = 1. / (2. * rdt)
    IF (neuler == 0 .AND. kt == nit000) z1_2dt = 1. / rdt
    IF (ln_vvl_ztilde .OR. ln_vvl_layer) THEN
      ALLOCATE(zhdiv(jpi, jpj, jpk))
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zhdiv(ji, jj, jk) = r1_e1e2t(ji, jj) * (un_td(ji, jj, jk) - un_td(ji - 1, jj, jk) + vn_td(ji, jj, jk) - vn_td(ji, jj - 1, jk))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(zhdiv, 'T', 1.)
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
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = sshb, clinfo1 = ' sshb  - : ', mask1 = tmask)
    IF (ln_timing) CALL timing_stop('ssh_swp')
  END SUBROUTINE ssh_swp
END MODULE sshwzv