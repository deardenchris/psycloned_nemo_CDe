MODULE traldf_triad
  USE oce
  USE dom_oce
  USE phycst
  USE trc_oce
  USE zdf_oce
  USE ldftra
  USE ldfslp
  USE traldf_iso
  USE diaptr
  USE diaar5
  USE zpshde
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_ldf_triad
  REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE, SAVE :: zdkt3d
  LOGICAL :: l_ptr
  LOGICAL :: l_hst
  CONTAINS
  SUBROUTINE tra_ldf_triad(kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, ptbb, pta, kjpt, kpass)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kit000
    CHARACTER(LEN = 3), INTENT(IN) :: cdtype
    INTEGER, INTENT(IN) :: kjpt
    INTEGER, INTENT(IN) :: kpass
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pahu, pahv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgu, pgv
    REAL(KIND = wp), DIMENSION(jpi, jpj, kjpt), INTENT(IN) :: pgui, pgvi
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN) :: ptbb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    INTEGER :: ip, jp, kp
    INTEGER :: ierr
    REAL(KIND = wp) :: zmsku, zabe1, zcof1, zcoef3
    REAL(KIND = wp) :: zmskv, zabe2, zcof2, zcoef4
    REAL(KIND = wp) :: zcoef0, ze3w_2, zsign, z2dt, z1_2dt
    REAL(KIND = wp) :: zslope_skew, zslope_iso, zslope2, zbu, zbv
    REAL(KIND = wp) :: ze1ur, ze2vr, ze3wr, zdxt, zdyt, zdzt
    REAL(KIND = wp) :: zah, zah_slp, zaei_slp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdit, zdjt, zftu, zftv, ztfw, zpsi_uw, zpsi_vw
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('tra_ldf_triad', 'r0', 0, 0)
    IF (.NOT. ALLOCATED(zdkt3d)) THEN
      ALLOCATE(zdkt3d(jpi, jpj, 0 : 1), STAT = ierr)
      CALL mpp_sum('traldf_triad', ierr)
      IF (ierr > 0) CALL ctl_stop('STOP', 'tra_ldf_triad: unable to allocate arrays')
    END IF
    IF (kpass == 1 .AND. kt == kit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_ldf_triad : rotated laplacian diffusion operator on ', cdtype
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~'
    END IF
    l_hst = .FALSE.
    l_ptr = .FALSE.
    IF (cdtype == 'TRA' .AND. ln_diaptr) l_ptr = .TRUE.
    IF (cdtype == 'TRA' .AND. (iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. iom_use("uadv_salttr") .OR. &
&iom_use("vadv_salttr"))) l_hst = .TRUE.
    IF (neuler == 0 .AND. kt == kit000) THEN
      z2dt = rdt
    ELSE
      z2dt = 2. * rdt
    END IF
    z1_2dt = 1._wp / z2dt
    IF (kpass == 1) THEN
      zsign = 1._wp
    ELSE
      zsign = - 1._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (kpass == 1) THEN
      !$ACC KERNELS
      akz(:, :, :) = 0._wp
      ah_wslp2(:, :, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_ldfeiv_dia) THEN
        !$ACC KERNELS
        zpsi_uw(:, :, :) = 0._wp
        zpsi_vw(:, :, :) = 0._wp
        !$ACC END KERNELS
      END IF
      CALL profile_psy_data1 % PreStart('tra_ldf_triad', 'r1', 0, 0)
      DO ip = 0, 1
        DO kp = 0, 1
          DO jk = 1, jpkm1
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                ze3wr = 1._wp / e3w_n(ji + ip, jj, jk + kp)
                zbu = e1e2u(ji, jj) * e3u_n(ji, jj, jk)
                zah = 0.25_wp * pahu(ji, jj, jk)
                zslope_skew = triadi_g(ji + ip, jj, jk, 1 - ip, kp)
                zslope2 = zslope_skew + (gdept_n(ji + 1, jj, jk) - gdept_n(ji, jj, jk)) * r1_e1u(ji, jj) * umask(ji, jj, jk + kp)
                zslope2 = zslope2 * zslope2
                ah_wslp2(ji + ip, jj, jk + kp) = ah_wslp2(ji + ip, jj, jk + kp) + zah * zbu * ze3wr * r1_e1e2t(ji + ip, jj) * &
&zslope2
                akz(ji + ip, jj, jk + kp) = akz(ji + ip, jj, jk + kp) + zah * r1_e1u(ji, jj) * r1_e1u(ji, jj) * umask(ji, jj, jk + &
&kp)
                IF (ln_ldfeiv_dia) zpsi_uw(ji, jj, jk + kp) = zpsi_uw(ji, jj, jk + kp) + 0.25_wp * aeiu(ji, jj, jk) * e2u(ji, jj) &
&* zslope_skew
              END DO
            END DO
          END DO
        END DO
      END DO
      DO jp = 0, 1
        DO kp = 0, 1
          DO jk = 1, jpkm1
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                ze3wr = 1.0_wp / e3w_n(ji, jj + jp, jk + kp)
                zbv = e1e2v(ji, jj) * e3v_n(ji, jj, jk)
                zah = 0.25_wp * pahv(ji, jj, jk)
                zslope_skew = triadj_g(ji, jj + jp, jk, 1 - jp, kp)
                zslope2 = zslope_skew + (gdept_n(ji, jj + 1, jk) - gdept_n(ji, jj, jk)) * r1_e2v(ji, jj) * vmask(ji, jj, jk + kp)
                zslope2 = zslope2 * zslope2
                ah_wslp2(ji, jj + jp, jk + kp) = ah_wslp2(ji, jj + jp, jk + kp) + zah * zbv * ze3wr * r1_e1e2t(ji, jj + jp) * &
&zslope2
                akz(ji, jj + jp, jk + kp) = akz(ji, jj + jp, jk + kp) + zah * r1_e2v(ji, jj) * r1_e2v(ji, jj) * vmask(ji, jj, jk + &
&kp)
                IF (ln_ldfeiv_dia) zpsi_vw(ji, jj, jk + kp) = zpsi_vw(ji, jj, jk + kp) + 0.25 * aeiv(ji, jj, jk) * e1v(ji, jj) * &
&zslope_skew
              END DO
            END DO
          END DO
        END DO
      END DO
      CALL profile_psy_data1 % PostEnd
      IF (ln_traldf_msc) THEN
        IF (ln_traldf_blp) THEN
          !$ACC KERNELS
          DO jk = 2, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                akz(ji, jj, jk) = 16._wp * ah_wslp2(ji, jj, jk) * (akz(ji, jj, jk) + ah_wslp2(ji, jj, jk) / (e3w_n(ji, jj, jk) * &
&e3w_n(ji, jj, jk)))
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        ELSE IF (ln_traldf_lap) THEN
          !$ACC KERNELS
          DO jk = 2, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpjm1
              DO ji = 1, jpim1
                ze3w_2 = e3w_n(ji, jj, jk) * e3w_n(ji, jj, jk)
                zcoef0 = z2dt * (akz(ji, jj, jk) + ah_wslp2(ji, jj, jk) / ze3w_2)
                akz(ji, jj, jk) = MAX(zcoef0 - 0.5_wp, 0._wp) * ze3w_2 * z1_2dt
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END IF
      ELSE
        !$ACC KERNELS
        akz(:, :, :) = ah_wslp2(:, :, :)
        !$ACC END KERNELS
      END IF
      IF (ln_ldfeiv_dia .AND. cdtype == 'TRA') CALL ldf_eiv_dia(zpsi_uw, zpsi_vw)
    END IF
    DO jn = 1, kjpt
      !$ACC KERNELS
      ztfw(:, :, :) = 0._wp
      zftu(:, :, :) = 0._wp
      zftv(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdit(ji, jj, jk) = (ptb(ji + 1, jj, jk, jn) - ptb(ji, jj, jk, jn)) * umask(ji, jj, jk)
            zdjt(ji, jj, jk) = (ptb(ji, jj + 1, jk, jn) - ptb(ji, jj, jk, jn)) * vmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      IF (ln_zps .AND. l_grad_zps) THEN
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zdit(ji, jj, mbku(ji, jj)) = pgu(ji, jj, jn)
            zdjt(ji, jj, mbkv(ji, jj)) = pgv(ji, jj, jn)
          END DO
        END DO
        !$ACC END KERNELS
        IF (ln_isfcav) THEN
          !$ACC KERNELS
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 1, jpim1
              IF (miku(ji, jj) > 1) zdit(ji, jj, miku(ji, jj)) = pgui(ji, jj, jn)
              IF (mikv(ji, jj) > 1) zdjt(ji, jj, mikv(ji, jj)) = pgvi(ji, jj, jn)
            END DO
          END DO
          !$ACC END KERNELS
        END IF
      END IF
      DO jk = 1, jpkm1
        !$ACC KERNELS
        zdkt3d(:, :, 1) = (ptb(:, :, jk, jn) - ptb(:, :, jk + 1, jn)) * tmask(:, :, jk + 1)
        IF (jk == 1) THEN
          zdkt3d(:, :, 0) = zdkt3d(:, :, 1)
        ELSE
          zdkt3d(:, :, 0) = (ptb(:, :, jk - 1, jn) - ptb(:, :, jk, jn)) * tmask(:, :, jk)
        END IF
        zaei_slp = 0._wp
        !$ACC END KERNELS
        CALL profile_psy_data2 % PreStart('tra_ldf_triad', 'r2', 0, 0)
        IF (ln_botmix_triad) THEN
          DO ip = 0, 1
            DO kp = 0, 1
              DO jj = 1, jpjm1
                DO ji = 1, jpim1
                  ze1ur = r1_e1u(ji, jj)
                  zdxt = zdit(ji, jj, jk) * ze1ur
                  ze3wr = 1._wp / e3w_n(ji + ip, jj, jk + kp)
                  zdzt = zdkt3d(ji + ip, jj, kp) * ze3wr
                  zslope_skew = triadi_g(ji + ip, jj, jk, 1 - ip, kp)
                  zslope_iso = triadi(ji + ip, jj, jk, 1 - ip, kp)
                  zbu = 0.25_wp * e1e2u(ji, jj) * e3u_n(ji, jj, jk)
                  zah = pahu(ji, jj, jk)
                  zah_slp = zah * zslope_iso
                  IF (ln_ldfeiv) zaei_slp = aeiu(ji, jj, jk) * zslope_skew
                  zftu(ji, jj, jk) = zftu(ji, jj, jk) - (zah * zdxt + (zah_slp - zaei_slp) * zdzt) * zbu * ze1ur
                  ztfw(ji + ip, jj, jk + kp) = ztfw(ji + ip, jj, jk + kp) - (zah_slp + zaei_slp) * zdxt * zbu * ze3wr
                END DO
              END DO
            END DO
          END DO
          DO jp = 0, 1
            DO kp = 0, 1
              DO jj = 1, jpjm1
                DO ji = 1, jpim1
                  ze2vr = r1_e2v(ji, jj)
                  zdyt = zdjt(ji, jj, jk) * ze2vr
                  ze3wr = 1._wp / e3w_n(ji, jj + jp, jk + kp)
                  zdzt = zdkt3d(ji, jj + jp, kp) * ze3wr
                  zslope_skew = triadj_g(ji, jj + jp, jk, 1 - jp, kp)
                  zslope_iso = triadj(ji, jj + jp, jk, 1 - jp, kp)
                  zbv = 0.25_wp * e1e2v(ji, jj) * e3v_n(ji, jj, jk)
                  zah = pahv(ji, jj, jk)
                  zah_slp = zah * zslope_iso
                  IF (ln_ldfeiv) zaei_slp = aeiv(ji, jj, jk) * zslope_skew
                  zftv(ji, jj, jk) = zftv(ji, jj, jk) - (zah * zdyt + (zah_slp - zaei_slp) * zdzt) * zbv * ze2vr
                  ztfw(ji, jj + jp, jk + kp) = ztfw(ji, jj + jp, jk + kp) - (zah_slp + zaei_slp) * zdyt * zbv * ze3wr
                END DO
              END DO
            END DO
          END DO
        ELSE
          DO ip = 0, 1
            DO kp = 0, 1
              DO jj = 1, jpjm1
                DO ji = 1, jpim1
                  ze1ur = r1_e1u(ji, jj)
                  zdxt = zdit(ji, jj, jk) * ze1ur
                  ze3wr = 1._wp / e3w_n(ji + ip, jj, jk + kp)
                  zdzt = zdkt3d(ji + ip, jj, kp) * ze3wr
                  zslope_skew = triadi_g(ji + ip, jj, jk, 1 - ip, kp)
                  zslope_iso = triadi(ji + ip, jj, jk, 1 - ip, kp)
                  zbu = 0.25_wp * e1e2u(ji, jj) * e3u_n(ji, jj, jk)
                  zah = pahu(ji, jj, jk) * umask(ji, jj, jk + kp)
                  zah_slp = zah * zslope_iso
                  IF (ln_ldfeiv) zaei_slp = aeiu(ji, jj, jk) * zslope_skew
                  zftu(ji, jj, jk) = zftu(ji, jj, jk) - (zah * zdxt + (zah_slp - zaei_slp) * zdzt) * zbu * ze1ur
                  ztfw(ji + ip, jj, jk + kp) = ztfw(ji + ip, jj, jk + kp) - (zah_slp + zaei_slp) * zdxt * zbu * ze3wr
                END DO
              END DO
            END DO
          END DO
          DO jp = 0, 1
            DO kp = 0, 1
              DO jj = 1, jpjm1
                DO ji = 1, jpim1
                  ze2vr = r1_e2v(ji, jj)
                  zdyt = zdjt(ji, jj, jk) * ze2vr
                  ze3wr = 1._wp / e3w_n(ji, jj + jp, jk + kp)
                  zdzt = zdkt3d(ji, jj + jp, kp) * ze3wr
                  zslope_skew = triadj_g(ji, jj + jp, jk, 1 - jp, kp)
                  zslope_iso = triadj(ji, jj + jp, jk, 1 - jp, kp)
                  zbv = 0.25_wp * e1e2v(ji, jj) * e3v_n(ji, jj, jk)
                  zah = pahv(ji, jj, jk) * vmask(ji, jj, jk + kp)
                  zah_slp = zah * zslope_iso
                  IF (ln_ldfeiv) zaei_slp = aeiv(ji, jj, jk) * zslope_skew
                  zftv(ji, jj, jk) = zftv(ji, jj, jk) - (zah * zdyt + (zah_slp - zaei_slp) * zdzt) * zbv * ze2vr
                  ztfw(ji, jj + jp, jk + kp) = ztfw(ji, jj + jp, jk + kp) - (zah_slp + zaei_slp) * zdyt * zbv * ze3wr
                END DO
              END DO
            END DO
          END DO
        END IF
        CALL profile_psy_data2 % PostEnd
        !$ACC KERNELS
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + zsign * (zftu(ji - 1, jj, jk) - zftu(ji, jj, jk) + zftv(ji, jj - 1, jk) - &
&zftv(ji, jj, jk)) / (e1e2t(ji, jj) * e3t_n(ji, jj, jk))
          END DO
        END DO
        !$ACC END KERNELS
      END DO
      IF (ln_traldf_lap) THEN
        !$ACC KERNELS
        DO jk = 2, jpkm1
          !$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO jj = 1, jpjm1
            DO ji = 2, jpim1
              ztfw(ji, jj, jk) = ztfw(ji, jj, jk) - e1e2t(ji, jj) / e3w_n(ji, jj, jk) * tmask(ji, jj, jk) * (ah_wslp2(ji, jj, jk) &
&- akz(ji, jj, jk)) * (ptb(ji, jj, jk - 1, jn) - ptb(ji, jj, jk, jn))
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        SELECT CASE (kpass)
        CASE (1)
          DO jk = 2, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpjm1
              DO ji = 2, jpim1
                ztfw(ji, jj, jk) = ztfw(ji, jj, jk) - e1e2t(ji, jj) / e3w_n(ji, jj, jk) * tmask(ji, jj, jk) * ah_wslp2(ji, jj, jk) &
&* (ptb(ji, jj, jk - 1, jn) - ptb(ji, jj, jk, jn))
              END DO
            END DO
          END DO
        CASE (2)
          DO jk = 2, jpkm1
            !$ACC LOOP INDEPENDENT COLLAPSE(2)
            DO jj = 1, jpjm1
              DO ji = 2, jpim1
                ztfw(ji, jj, jk) = ztfw(ji, jj, jk) - e1e2t(ji, jj) / e3w_n(ji, jj, jk) * tmask(ji, jj, jk) * (ah_wslp2(ji, jj, &
&jk) * (ptb(ji, jj, jk - 1, jn) - ptb(ji, jj, jk, jn)) + akz(ji, jj, jk) * (ptbb(ji, jj, jk - 1, jn) - ptbb(ji, jj, jk, jn)))
              END DO
            END DO
          END DO
        END SELECT
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = pta(ji, jj, jk, jn) + zsign * (ztfw(ji, jj, jk + 1) - ztfw(ji, jj, jk)) / (e1e2t(ji, jj) * &
&e3t_n(ji, jj, jk))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL profile_psy_data3 % PreStart('tra_ldf_triad', 'r3', 0, 0)
      IF ((kpass == 1 .AND. ln_traldf_lap) .OR. (kpass == 2 .AND. ln_traldf_blp)) THEN
        IF (l_ptr) CALL dia_ptr_hst(jn, 'ldf', zftv(:, :, :))
        IF (l_hst) CALL dia_ar5_hst(jn, 'ldf', zftu(:, :, :), zftv(:, :, :))
      END IF
      CALL profile_psy_data3 % PostEnd
    END DO
  END SUBROUTINE tra_ldf_triad
END MODULE traldf_triad