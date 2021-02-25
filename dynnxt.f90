MODULE dynnxt
  USE oce
  USE dom_oce
  USE sbc_oce
  USE sbcrnf
  USE sbcisf
  USE phycst
  USE dynadv
  USE dynspg_ts
  USE domvvl
  USE bdy_oce, ONLY: ln_bdy
  USE bdydta
  USE bdydyn
  USE bdyvol
  USE trd_oce
  USE trddyn
  USE trdken
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  USE prtctl
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_nxt
  CONTAINS
  SUBROUTINE dyn_nxt(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikt
    REAL(KIND = wp) :: zue3a, zue3n, zue3b, zuf, zcoef
    REAL(KIND = wp) :: zve3a, zve3n, zve3b, zvf, z1_2dt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zue, zve
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ze3u_f, ze3v_f, zua, zva
    IF (ln_timing) CALL timing_start('dyn_nxt')
    IF (ln_dynspg_ts) ALLOCATE(zue(jpi, jpj), zve(jpi, jpj))
    IF (l_trddyn) ALLOCATE(zua(jpi, jpj, jpk), zva(jpi, jpj, jpk))
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_nxt : time stepping'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    IF (ln_dynspg_ts) THEN
      zue(:, :) = e3u_a(:, :, 1) * ua(:, :, 1) * umask(:, :, 1)
      zve(:, :) = e3v_a(:, :, 1) * va(:, :, 1) * vmask(:, :, 1)
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        zue(:, :) = zue(:, :) + e3u_a(:, :, jk) * ua(:, :, jk) * umask(:, :, jk)
        zve(:, :) = zve(:, :) + e3v_a(:, :, jk) * va(:, :, jk) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        ua(:, :, jk) = (ua(:, :, jk) - zue(:, :) * r1_hu_a(:, :) + ua_b(:, :)) * umask(:, :, jk)
        va(:, :, jk) = (va(:, :, jk) - zve(:, :) * r1_hv_a(:, :) + va_b(:, :)) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (.NOT. ln_bt_fw) THEN
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          un(:, :, jk) = (un(:, :, jk) - un_adv(:, :) * r1_hu_n(:, :) + un_b(:, :)) * umask(:, :, jk)
          vn(:, :, jk) = (vn(:, :, jk) - vn_adv(:, :) * r1_hv_n(:, :) + vn_b(:, :)) * vmask(:, :, jk)
        END DO
        !$OMP end do
        !$OMP end parallel
      END IF
    END IF
    CALL lbc_lnk_multi('dynnxt', ua, 'U', - 1., va, 'V', - 1.)
    IF (ln_bdy .AND. ln_dynspg_exp) CALL bdy_dyn(kt)
    IF (ln_bdy .AND. ln_dynspg_ts) CALL bdy_dyn(kt, dyn3d_only = .TRUE.)
    IF (l_trddyn) THEN
      z1_2dt = 1._wp / (2. * rdt)
      IF (neuler == 0 .AND. kt == nit000) z1_2dt = 1._wp / rdt
      IF (ln_KE_trd) CALL trd_dyn(ua, va, jpdyn_ken, kt)
      IF (ln_dyn_trd) THEN
        zua(:, :, :) = (ua(:, :, :) - ub(:, :, :)) * z1_2dt
        zva(:, :, :) = (va(:, :, :) - vb(:, :, :)) * z1_2dt
        CALL iom_put("utrd_tot", zua)
        CALL iom_put("vtrd_tot", zva)
      END IF
      zua(:, :, :) = un(:, :, :)
      zva(:, :, :) = vn(:, :, :)
    END IF
    IF (neuler == 0 .AND. kt == nit000) THEN
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        un(:, :, jk) = ua(:, :, jk)
        vn(:, :, jk) = va(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
      IF (.NOT. ln_linssh) THEN
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          e3t_n(:, :, jk) = e3t_a(:, :, jk)
          e3u_n(:, :, jk) = e3u_a(:, :, jk)
          e3v_n(:, :, jk) = e3v_a(:, :, jk)
        END DO
        !$OMP end do
        !$OMP end parallel
      END IF
    ELSE
      IF (ln_linssh) THEN
        !$OMP parallel default(shared), private(ji,jj,jk,zuf,zvf)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              zuf = un(ji, jj, jk) + atfp * (ub(ji, jj, jk) - 2._wp * un(ji, jj, jk) + ua(ji, jj, jk))
              zvf = vn(ji, jj, jk) + atfp * (vb(ji, jj, jk) - 2._wp * vn(ji, jj, jk) + va(ji, jj, jk))
              ub(ji, jj, jk) = zuf
              vb(ji, jj, jk) = zvf
              un(ji, jj, jk) = ua(ji, jj, jk)
              vn(ji, jj, jk) = va(ji, jj, jk)
            END DO
          END DO
        END DO
        !$OMP end do
        !$OMP end parallel
      ELSE
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          e3t_b(:, :, jk) = e3t_n(:, :, jk) + atfp * (e3t_b(:, :, jk) - 2._wp * e3t_n(:, :, jk) + e3t_a(:, :, jk))
        END DO
        !$OMP end do
        !$OMP end parallel
        zcoef = atfp * rdt * r1_rau0
        e3t_b(:, :, 1) = e3t_b(:, :, 1) - zcoef * (emp_b(:, :) - emp(:, :)) * tmask(:, :, 1)
        IF (ln_rnf) THEN
          IF (ln_rnf_depth) THEN
            !$OMP parallel default(shared), private(ji,jj,jk)
            !$OMP do schedule(static)
            DO jk = 1, jpkm1
              DO jj = 1, jpj
                DO ji = 1, jpi
                  IF (jk <= nk_rnf(ji, jj)) THEN
                    e3t_b(ji, jj, jk) = e3t_b(ji, jj, jk) - zcoef * (- rnf_b(ji, jj) + rnf(ji, jj)) * (e3t_n(ji, jj, jk) / &
&h_rnf(ji, jj)) * tmask(ji, jj, jk)
                  END IF
                END DO
              END DO
            END DO
            !$OMP end do
            !$OMP end parallel
          ELSE
            e3t_b(:, :, 1) = e3t_b(:, :, 1) - zcoef * (- rnf_b(:, :) + rnf(:, :)) * tmask(:, :, 1)
          END IF
        END IF
        IF (ln_isf) THEN
          !$OMP parallel default(shared), private(ji,jj,jk)
          !$OMP do schedule(static)
          DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji = 1, jpi
                IF (misfkt(ji, jj) <= jk .AND. jk < misfkb(ji, jj)) THEN
                  e3t_b(ji, jj, jk) = e3t_b(ji, jj, jk) - zcoef * (fwfisf_b(ji, jj) - fwfisf(ji, jj)) * (e3t_n(ji, jj, jk) * &
&r1_hisf_tbl(ji, jj)) * tmask(ji, jj, jk)
                ELSE IF (jk == misfkb(ji, jj)) THEN
                  e3t_b(ji, jj, jk) = e3t_b(ji, jj, jk) - zcoef * (fwfisf_b(ji, jj) - fwfisf(ji, jj)) * (e3t_n(ji, jj, jk) * &
&r1_hisf_tbl(ji, jj)) * ralpha(ji, jj) * tmask(ji, jj, jk)
                END IF
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
        END IF
        IF (ln_dynadv_vec) THEN
          CALL dom_vvl_interpol(e3t_b(:, :, :), e3u_b(:, :, :), 'U')
          CALL dom_vvl_interpol(e3t_b(:, :, :), e3v_b(:, :, :), 'V')
          !$OMP parallel default(shared), private(ji,jj,jk,zuf,zvf)
          !$OMP do schedule(static)
          DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji = 1, jpi
                zuf = un(ji, jj, jk) + atfp * (ub(ji, jj, jk) - 2._wp * un(ji, jj, jk) + ua(ji, jj, jk))
                zvf = vn(ji, jj, jk) + atfp * (vb(ji, jj, jk) - 2._wp * vn(ji, jj, jk) + va(ji, jj, jk))
                ub(ji, jj, jk) = zuf
                vb(ji, jj, jk) = zvf
                un(ji, jj, jk) = ua(ji, jj, jk)
                vn(ji, jj, jk) = va(ji, jj, jk)
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
        ELSE
          ALLOCATE(ze3u_f(jpi, jpj, jpk), ze3v_f(jpi, jpj, jpk))
          CALL dom_vvl_interpol(e3t_b(:, :, :), ze3u_f, 'U')
          CALL dom_vvl_interpol(e3t_b(:, :, :), ze3v_f, 'V')
          !$OMP parallel default(shared), private(ji,jj,jk,zue3a,zue3b,zue3n,zuf,zve3a,zve3b,zve3n,zvf)
          !$OMP do schedule(static)
          DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji = 1, jpi
                zue3a = e3u_a(ji, jj, jk) * ua(ji, jj, jk)
                zve3a = e3v_a(ji, jj, jk) * va(ji, jj, jk)
                zue3n = e3u_n(ji, jj, jk) * un(ji, jj, jk)
                zve3n = e3v_n(ji, jj, jk) * vn(ji, jj, jk)
                zue3b = e3u_b(ji, jj, jk) * ub(ji, jj, jk)
                zve3b = e3v_b(ji, jj, jk) * vb(ji, jj, jk)
                zuf = (zue3n + atfp * (zue3b - 2._wp * zue3n + zue3a)) / ze3u_f(ji, jj, jk)
                zvf = (zve3n + atfp * (zve3b - 2._wp * zve3n + zve3a)) / ze3v_f(ji, jj, jk)
                ub(ji, jj, jk) = zuf
                vb(ji, jj, jk) = zvf
                un(ji, jj, jk) = ua(ji, jj, jk)
                vn(ji, jj, jk) = va(ji, jj, jk)
              END DO
            END DO
          END DO
          !$OMP end do
          !$OMP end parallel
          e3u_b(:, :, 1 : jpkm1) = ze3u_f(:, :, 1 : jpkm1)
          e3v_b(:, :, 1 : jpkm1) = ze3v_f(:, :, 1 : jpkm1)
          DEALLOCATE(ze3u_f, ze3v_f)
        END IF
      END IF
      IF (ln_dynspg_ts .AND. ln_bt_fw) THEN
        zue(:, :) = e3u_b(:, :, 1) * ub(:, :, 1) * umask(:, :, 1)
        zve(:, :) = e3v_b(:, :, 1) * vb(:, :, 1) * vmask(:, :, 1)
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 2, jpkm1
          zue(:, :) = zue(:, :) + e3u_b(:, :, jk) * ub(:, :, jk) * umask(:, :, jk)
          zve(:, :) = zve(:, :) + e3v_b(:, :, jk) * vb(:, :, jk) * vmask(:, :, jk)
        END DO
        !$OMP end do
        !$OMP do schedule(static)
        DO jk = 1, jpkm1
          ub(:, :, jk) = ub(:, :, jk) - (zue(:, :) * r1_hu_n(:, :) - un_b(:, :)) * umask(:, :, jk)
          vb(:, :, jk) = vb(:, :, jk) - (zve(:, :) * r1_hv_n(:, :) - vn_b(:, :)) * vmask(:, :, jk)
        END DO
        !$OMP end do
        !$OMP end parallel
      END IF
    END IF
    IF (.NOT. ln_linssh) THEN
      hu_b(:, :) = e3u_b(:, :, 1) * umask(:, :, 1)
      hv_b(:, :) = e3v_b(:, :, 1) * vmask(:, :, 1)
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 2, jpkm1
        hu_b(:, :) = hu_b(:, :) + e3u_b(:, :, jk) * umask(:, :, jk)
        hv_b(:, :) = hv_b(:, :) + e3v_b(:, :, jk) * vmask(:, :, jk)
      END DO
      !$OMP end do
      !$OMP end parallel
      r1_hu_b(:, :) = ssumask(:, :) / (hu_b(:, :) + 1._wp - ssumask(:, :))
      r1_hv_b(:, :) = ssvmask(:, :) / (hv_b(:, :) + 1._wp - ssvmask(:, :))
    END IF
    un_b(:, :) = e3u_a(:, :, 1) * un(:, :, 1) * umask(:, :, 1)
    ub_b(:, :) = e3u_b(:, :, 1) * ub(:, :, 1) * umask(:, :, 1)
    vn_b(:, :) = e3v_a(:, :, 1) * vn(:, :, 1) * vmask(:, :, 1)
    vb_b(:, :) = e3v_b(:, :, 1) * vb(:, :, 1) * vmask(:, :, 1)
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 2, jpkm1
      un_b(:, :) = un_b(:, :) + e3u_a(:, :, jk) * un(:, :, jk) * umask(:, :, jk)
      ub_b(:, :) = ub_b(:, :) + e3u_b(:, :, jk) * ub(:, :, jk) * umask(:, :, jk)
      vn_b(:, :) = vn_b(:, :) + e3v_a(:, :, jk) * vn(:, :, jk) * vmask(:, :, jk)
      vb_b(:, :) = vb_b(:, :) + e3v_b(:, :, jk) * vb(:, :, jk) * vmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    un_b(:, :) = un_b(:, :) * r1_hu_a(:, :)
    vn_b(:, :) = vn_b(:, :) * r1_hv_a(:, :)
    ub_b(:, :) = ub_b(:, :) * r1_hu_b(:, :)
    vb_b(:, :) = vb_b(:, :) * r1_hv_b(:, :)
    IF (.NOT. ln_dynspg_ts) THEN
      CALL iom_put("ubar", un_b(:, :))
      CALL iom_put("vbar", vn_b(:, :))
    END IF
    IF (l_trddyn) THEN
      zua(:, :, :) = (ub(:, :, :) - zua(:, :, :)) * z1_2dt
      zva(:, :, :) = (vb(:, :, :) - zva(:, :, :)) * z1_2dt
      CALL trd_dyn(zua, zva, jpdyn_atf, kt)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = un, clinfo1 = ' nxt  - Un: ', mask1 = umask, tab3d_2 = vn, clinfo2 = ' Vn: ', mask2 = vmask)
    IF (ln_dynspg_ts) DEALLOCATE(zue, zve)
    IF (l_trddyn) DEALLOCATE(zua, zva)
    IF (ln_timing) CALL timing_stop('dyn_nxt')
  END SUBROUTINE dyn_nxt
END MODULE dynnxt