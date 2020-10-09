MODULE istate
  USE oce
  USE dom_oce
  USE daymod
  USE divhor
  USE dtatsd
  USE dtauvd
  USE domvvl
  USE iscplrst
  USE wet_dry
  USE usrdef_istate
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE restart
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: istate_init
  CONTAINS
  SUBROUTINE istate_init
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: zuvd
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'istate_init : Initialization of the dynamics and tracers'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
    CALL dta_tsd_init
    !$ACC KERNELS
    rhd(:, :, :) = 0._wp
    rhop(:, :, :) = 0._wp
    rn2b(:, :, :) = 0._wp
    rn2(:, :, :) = 0._wp
    tsa(:, :, :, :) = 0._wp
    rab_b(:, :, :, :) = 0._wp
    rab_n(:, :, :, :) = 0._wp
    !$ACC END KERNELS
    IF (ln_rstart) THEN
      CALL rst_read
      IF (ln_iscpl) CALL iscpl_stp
      CALL day_init
    ELSE
      numror = 0
      neuler = 0
      CALL day_init
      IF (ln_tsd_init) THEN
        CALL dta_tsd(nit000, tsb)
        !$ACC KERNELS
        sshb(:, :) = 0._wp
        !$ACC END KERNELS
        IF (ll_wd) THEN
          !$ACC KERNELS
          sshb(:, :) = - ssh_ref
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              IF (ht_0(ji, jj) + sshb(ji, jj) < rn_wdmin1) THEN
                sshb(ji, jj) = tmask(ji, jj, 1) * (rn_wdmin1 - (ht_0(ji, jj)))
              END IF
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        ub(:, :, :) = 0._wp
        vb(:, :, :) = 0._wp
        !$ACC END KERNELS
      ELSE
        CALL usr_def_istate(gdept_b, tmask, tsb, ub, vb, sshb)
      END IF
      !$ACC KERNELS
      tsn(:, :, :, :) = tsb(:, :, :, :)
      sshn(:, :) = sshb(:, :)
      un(:, :, :) = ub(:, :, :)
      vn(:, :, :) = vb(:, :, :)
      hdivn(:, :, jpk) = 0._wp
      !$ACC END KERNELS
      CALL div_hor(0)
    END IF
    !$ACC KERNELS
    un_b(:, :) = 0._wp
    vn_b(:, :) = 0._wp
    ub_b(:, :) = 0._wp
    vb_b(:, :) = 0._wp
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          un_b(ji, jj) = un_b(ji, jj) + e3u_n(ji, jj, jk) * un(ji, jj, jk) * umask(ji, jj, jk)
          vn_b(ji, jj) = vn_b(ji, jj) + e3v_n(ji, jj, jk) * vn(ji, jj, jk) * vmask(ji, jj, jk)
          ub_b(ji, jj) = ub_b(ji, jj) + e3u_b(ji, jj, jk) * ub(ji, jj, jk) * umask(ji, jj, jk)
          vb_b(ji, jj) = vb_b(ji, jj) + e3v_b(ji, jj, jk) * vb(ji, jj, jk) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    un_b(:, :) = un_b(:, :) * r1_hu_n(:, :)
    vn_b(:, :) = vn_b(:, :) * r1_hv_n(:, :)
    ub_b(:, :) = ub_b(:, :) * r1_hu_b(:, :)
    vb_b(:, :) = vb_b(:, :) * r1_hv_b(:, :)
    !$ACC END KERNELS
  END SUBROUTINE istate_init
END MODULE istate