MODULE iscplrst
  USE oce
  USE dom_oce
  USE domwri
  USE domvvl, ONLY: dom_vvl_interpol
  USE phycst
  USE sbc_oce
  USE iscplini
  USE iscplhsb
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iscpl_stp
  CONTAINS
  SUBROUTINE iscpl_stp
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: inum0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsmask_b
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztmask_b, zumask_b, zvmask_b
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ze3t_b, ze3u_b, ze3v_b
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdepw_b
    CHARACTER(LEN = 20) :: cfile
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iscpl_stp', 'r0', 0, 0)
    CALL iom_get(numror, jpdom_autoglo, 'tmask', ztmask_b, ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'umask', zumask_b, ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'vmask', zvmask_b, ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'smask', zsmask_b, ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'e3t_n', ze3t_b(:, :, :), ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'e3u_n', ze3u_b(:, :, :), ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'e3v_n', ze3v_b(:, :, :), ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'gdepw_n', zdepw_b(:, :, :), ldxios = lrxios)
    CALL iscpl_init
    CALL iscpl_rst_interpol(ztmask_b, zumask_b, zvmask_b, zsmask_b, ze3t_b, ze3u_b, ze3v_b, zdepw_b)
    IF (ln_hsb) THEN
      IF (iscpl_alloc() /= 0) CALL ctl_stop('STOP', 'rst_iscpl : unable to allocate rst_iscpl arrays')
      CALL iscpl_cons(ztmask_b, zsmask_b, ze3t_b, htsc_iscpl, hdiv_iscpl, rdt_iscpl)
    END IF
    IF (ln_meshmask .AND. ln_iscpl) CALL dom_wri
    IF (ln_hsb) THEN
      cfile = 'correction'
      cfile = TRIM(cfile)
      CALL iom_open(cfile, inum0, ldwrt = .TRUE.)
      CALL iom_rstput(0, 0, inum0, 'vol_cor', hdiv_iscpl(:, :, :))
      CALL iom_rstput(0, 0, inum0, 'tem_cor', htsc_iscpl(:, :, :, jp_tem))
      CALL iom_rstput(0, 0, inum0, 'sal_cor', htsc_iscpl(:, :, :, jp_sal))
      CALL iom_close(inum0)
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    neuler = 0
    tsb(:, :, :, :) = tsn(:, :, :, :)
    ub(:, :, :) = un(:, :, :)
    vb(:, :, :) = vn(:, :, :)
    sshb(:, :) = sshn(:, :)
    e3t_b(:, :, :) = e3t_n(:, :, :)
    e3u_b(:, :, :) = e3u_n(:, :, :)
    e3v_b(:, :, :) = e3v_n(:, :, :)
    e3uw_b(:, :, :) = e3uw_n(:, :, :)
    e3vw_b(:, :, :) = e3vw_n(:, :, :)
    gdept_b(:, :, :) = gdept_n(:, :, :)
    gdepw_b(:, :, :) = gdepw_n(:, :, :)
    hu_b(:, :) = hu_n(:, :)
    hv_b(:, :) = hv_n(:, :)
    r1_hu_b(:, :) = r1_hu_n(:, :)
    r1_hv_b(:, :) = r1_hv_n(:, :)
    !$ACC END KERNELS
  END SUBROUTINE iscpl_stp
  SUBROUTINE iscpl_rst_interpol(ptmask_b, pumask_b, pvmask_b, psmask_b, pe3t_b, pe3u_b, pe3v_b, pdepw_b)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptmask_b, pumask_b, pvmask_b
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pe3t_b, pe3u_b, pe3v_b
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pdepw_b
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: psmask_b
    INTEGER :: ji, jj, jk, iz
    INTEGER :: jip1, jim1, jjp1, jjm1, jkp1, jkm1
    REAL(KIND = wp) :: summsk, zsum, zsum1, zarea, zsumn, zsumb
    REAL(KIND = wp) :: zdz, zdzm1, zdzp1
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdmask, zsmask0, zucorr, zbub, zbun, zssh0, zhu1, zde3t
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdsmask, zsmask1, zvcorr, zbvb, zbvn, zssh1, zhv1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztmask0, zwmaskn, ztrp
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ztmask1, zwmaskb, ztmp3d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts) :: zts0
    !$ACC KERNELS
    tsn(:, :, :, jp_tem) = tsn(:, :, :, jp_tem) * ptmask_b(:, :, :)
    tsn(:, :, :, jp_sal) = tsn(:, :, :, jp_sal) * ptmask_b(:, :, :)
    zwmaskn(:, :, 1) = tmask(:, :, 1)
    zwmaskb(:, :, 1) = ptmask_b(:, :, 1)
    !$ACC END KERNELS
    DO jk = 2, jpk
      !$ACC KERNELS
      zwmaskn(:, :, jk) = tmask(:, :, jk) * tmask(:, :, jk - 1)
      zwmaskb(:, :, jk) = ptmask_b(:, :, jk) * ptmask_b(:, :, jk - 1)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    sshb(:, :) = sshn(:, :)
    zssh0(:, :) = sshn(:, :)
    zsmask0(:, :) = psmask_b(:, :)
    zsmask1(:, :) = psmask_b(:, :)
    !$ACC END KERNELS
    DO iz = 1, 10
      !$ACC KERNELS
      zdsmask(:, :) = ssmask(:, :) - zsmask0(:, :)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpj - 1
        DO ji = 2, jpim1
          jip1 = ji + 1
          jim1 = ji - 1
          jjp1 = jj + 1
          jjm1 = jj - 1
          summsk = (zsmask0(jip1, jj) + zsmask0(jim1, jj) + zsmask0(ji, jjp1) + zsmask0(ji, jjm1))
          IF (zdsmask(ji, jj) == 1._wp .AND. summsk /= 0._wp) THEN
            sshn(ji, jj) = (zssh0(jip1, jj) * zsmask0(jip1, jj) + zssh0(jim1, jj) * zsmask0(jim1, jj) + zssh0(ji, jjp1) * &
&zsmask0(ji, jjp1) + zssh0(ji, jjm1) * zsmask0(ji, jjm1)) / summsk
            zsmask1(ji, jj) = 1._wp
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi('iscplrst', sshn, 'T', 1., zsmask1, 'T', 1.)
      zssh0 = sshn
      zsmask0 = zsmask1
    END DO
    !$ACC KERNELS
    sshn(:, :) = sshn(:, :) * ssmask(:, :)
    !$ACC END KERNELS
    IF (.NOT. ln_linssh) THEN
      !$ACC KERNELS
      DO jk = 1, jpk
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (tmask(ji, jj, 1) == 0._wp .OR. ptmask_b(ji, jj, 1) == 0._wp) THEN
              e3t_n(ji, jj, jk) = e3t_0(ji, jj, jk) * (1._wp + sshn(ji, jj) / (ht_0(ji, jj) + 1._wp - ssmask(ji, jj)) * tmask(ji, &
&jj, jk))
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL dom_vvl_interpol(e3t_n(:, :, :), e3u_n(:, :, :), 'U')
      CALL dom_vvl_interpol(e3t_n(:, :, :), e3v_n(:, :, :), 'V')
      CALL dom_vvl_interpol(e3u_n(:, :, :), e3f_n(:, :, :), 'F')
      CALL dom_vvl_interpol(e3t_n(:, :, :), e3w_n(:, :, :), 'W')
      CALL dom_vvl_interpol(e3u_n(:, :, :), e3uw_n(:, :, :), 'UW')
      CALL dom_vvl_interpol(e3v_n(:, :, :), e3vw_n(:, :, :), 'VW')
      !$ACC KERNELS
      gdept_n(:, :, 1) = 0.5_wp * e3w_n(:, :, 1)
      gdepw_n(:, :, 1) = 0.0_wp
      gde3w_n(:, :, 1) = gdept_n(:, :, 1) - sshn(:, :)
      !$ACC END KERNELS
      DO jj = 1, jpj
        DO ji = 1, jpi
          !$ACC KERNELS
          DO jk = 2, mikt(ji, jj) - 1
            gdept_n(ji, jj, jk) = gdept_0(ji, jj, jk)
            gdepw_n(ji, jj, jk) = gdepw_0(ji, jj, jk)
            gde3w_n(ji, jj, jk) = gdept_0(ji, jj, jk) - sshn(ji, jj)
          END DO
          IF (mikt(ji, jj) > 1) THEN
            jk = mikt(ji, jj)
            gdept_n(ji, jj, jk) = gdepw_0(ji, jj, jk) + 0.5_wp * e3w_n(ji, jj, jk)
            gdepw_n(ji, jj, jk) = gdepw_0(ji, jj, jk)
            gde3w_n(ji, jj, jk) = gdept_n(ji, jj, jk) - sshn(ji, jj)
          END IF
          DO jk = mikt(ji, jj) + 1, jpk
            gdept_n(ji, jj, jk) = gdept_n(ji, jj, jk - 1) + e3w_n(ji, jj, jk)
            gdepw_n(ji, jj, jk) = gdepw_n(ji, jj, jk - 1) + e3t_n(ji, jj, jk - 1)
            gde3w_n(ji, jj, jk) = gdept_n(ji, jj, jk) - sshn(ji, jj)
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      !$ACC KERNELS
      ht_n(:, :) = 0._wp
      hu_n(:, :) = 0._wp
      hv_n(:, :) = 0._wp
      !$ACC END KERNELS
      DO jk = 1, jpkm1
        !$ACC KERNELS
        hu_n(:, :) = hu_n(:, :) + e3u_n(:, :, jk) * umask(:, :, jk)
        hv_n(:, :) = hv_n(:, :) + e3v_n(:, :, jk) * vmask(:, :, jk)
        ht_n(:, :) = ht_n(:, :) + e3t_n(:, :, jk) * tmask(:, :, jk)
        !$ACC END KERNELS
      END DO
      !$ACC KERNELS
      r1_hu_n(:, :) = 1._wp / (hu_n(:, :) + 1._wp - ssumask(:, :)) * ssumask(:, :)
      r1_hv_n(:, :) = 1._wp / (hv_n(:, :) + 1._wp - ssvmask(:, :)) * ssvmask(:, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    ub(:, :, :) = un(:, :, :)
    vb(:, :, :) = vn(:, :, :)
    DO jk = 1, jpk
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          un(ji, jj, jk) = ub(ji, jj, jk) * pe3u_b(ji, jj, jk) * pumask_b(ji, jj, jk) / e3u_n(ji, jj, jk) * umask(ji, jj, jk)
          vn(ji, jj, jk) = vb(ji, jj, jk) * pe3v_b(ji, jj, jk) * pvmask_b(ji, jj, jk) / e3v_n(ji, jj, jk) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    ztrp(:, :, :) = ub(:, :, :) * pe3u_b(:, :, :)
    zbub(:, :) = SUM(ztrp, DIM = 3)
    ztrp(:, :, :) = vb(:, :, :) * pe3v_b(:, :, :)
    zbvb(:, :) = SUM(ztrp, DIM = 3)
    ztrp(:, :, :) = un(:, :, :) * e3u_n(:, :, :)
    zbun(:, :) = SUM(ztrp, DIM = 3)
    ztrp(:, :, :) = vn(:, :, :) * e3v_n(:, :, :)
    zbvn(:, :) = SUM(ztrp, DIM = 3)
    zhu1 = 0.0_wp
    zhv1 = 0.0_wp
    !$ACC END KERNELS
    DO jk = 1, jpk
      !$ACC KERNELS
      zhu1(:, :) = zhu1(:, :) + e3u_n(:, :, jk) * umask(:, :, jk)
      zhv1(:, :) = zhv1(:, :) + e3v_n(:, :, jk) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    zucorr = 0._wp
    zvcorr = 0._wp
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (zbun(ji, jj) /= zbub(ji, jj) .AND. zhu1(ji, jj) /= 0._wp) THEN
          zucorr(ji, jj) = (zbun(ji, jj) - zbub(ji, jj)) / zhu1(ji, jj)
        END IF
        IF (zbvn(ji, jj) /= zbvb(ji, jj) .AND. zhv1(ji, jj) /= 0._wp) THEN
          zvcorr(ji, jj) = (zbvn(ji, jj) - zbvb(ji, jj)) / zhv1(ji, jj)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    DO jk = 1, jpk
      !$ACC KERNELS
      un(:, :, jk) = (un(:, :, jk) - zucorr(:, :)) * umask(:, :, jk)
      vn(:, :, jk) = (vn(:, :, jk) - zvcorr(:, :)) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    tsb(:, :, :, :) = tsn(:, :, :, :)
    zts0(:, :, :, :) = tsn(:, :, :, :)
    ztmask1(:, :, :) = ptmask_b(:, :, :)
    ztmask0(:, :, :) = ptmask_b(:, :, :)
    !$ACC END KERNELS
    DO iz = 1, nn_drown
      !$ACC KERNELS
      DO jk = 1, jpk - 1
        zdmask = tmask(:, :, jk) - ztmask0(:, :, jk)
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpj - 1
          DO ji = 2, jpim1
            jip1 = ji + 1
            jim1 = ji - 1
            jjp1 = jj + 1
            jjm1 = jj - 1
            summsk = (ztmask0(jip1, jj, jk) + ztmask0(jim1, jj, jk) + ztmask0(ji, jjp1, jk) + ztmask0(ji, jjm1, jk))
            IF (zdmask(ji, jj) == 1._wp .AND. summsk /= 0._wp) THEN
              tsn(ji, jj, jk, 1) = (zts0(jip1, jj, jk, 1) * ztmask0(jip1, jj, jk) + zts0(jim1, jj, jk, 1) * ztmask0(jim1, jj, jk) &
&+ zts0(ji, jjp1, jk, 1) * ztmask0(ji, jjp1, jk) + zts0(ji, jjm1, jk, 1) * ztmask0(ji, jjm1, jk)) / summsk
              tsn(ji, jj, jk, 2) = (zts0(jip1, jj, jk, 2) * ztmask0(jip1, jj, jk) + zts0(jim1, jj, jk, 2) * ztmask0(jim1, jj, jk) &
&+ zts0(ji, jjp1, jk, 2) * ztmask0(ji, jjp1, jk) + zts0(ji, jjm1, jk, 2) * ztmask0(ji, jjm1, jk)) / summsk
              ztmask1(ji, jj, jk) = 1
            ELSE IF (zdmask(ji, jj) == 1._wp .AND. summsk == 0._wp) THEN
              jkm1 = MAX(1, jk - 1)
              jkp1 = MIN(jpk, jk + 1)
              summsk = (ztmask0(ji, jj, jkm1) + ztmask0(ji, jj, jkp1))
              IF (zdmask(ji, jj) == 1._wp .AND. summsk /= 0._wp) THEN
                tsn(ji, jj, jk, 1) = (zts0(ji, jj, jkp1, 1) * ztmask0(ji, jj, jkp1) + zts0(ji, jj, jkm1, 1) * ztmask0(ji, jj, &
&jkm1)) / summsk
                tsn(ji, jj, jk, 2) = (zts0(ji, jj, jkp1, 2) * ztmask0(ji, jj, jkp1) + zts0(ji, jj, jkm1, 2) * ztmask0(ji, jj, &
&jkm1)) / summsk
                ztmask1(ji, jj, jk) = 1._wp
              END IF
            END IF
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi('iscplrst', tsn(:, :, :, jp_tem), 'T', 1., tsn(:, :, :, jp_sal), 'T', 1., ztmask1, 'T', 1.)
      !$ACC KERNELS
      zts0(:, :, :, :) = tsn(:, :, :, :)
      ztmask0 = ztmask1
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    tsn(:, :, :, jp_tem) = tsn(:, :, :, jp_tem) * tmask(:, :, :)
    tsn(:, :, :, jp_sal) = tsn(:, :, :, jp_sal) * tmask(:, :, :)
    !$ACC END KERNELS
    IF (.NOT. ln_linssh) THEN
      DO jk = 2, jpk - 1
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (zwmaskn(ji, jj, jk) * zwmaskb(ji, jj, jk) == 1._wp .AND. (tmask(ji, jj, 1) == 0._wp .OR. ptmask_b(ji, jj, 1) == &
&0._wp)) THEN
              zdzp1 = MAX(0._wp, gdepw_n(ji, jj, jk + 1) - pdepw_b(ji, jj, jk + 1))
              zdz = gdepw_n(ji, jj, jk + 1) - pdepw_b(ji, jj, jk)
              zdzm1 = MAX(0._wp, pdepw_b(ji, jj, jk) - gdepw_n(ji, jj, jk))
              IF (zdz .LT. 0._wp) THEN
                CALL ctl_stop('STOP', 'rst_iscpl : unable to compute the interpolation')
              END IF
              tsn(ji, jj, jk, jp_tem) = (zdzp1 * tsb(ji, jj, jk + 1, jp_tem) + zdz * tsb(ji, jj, jk, jp_tem) + zdzm1 * tsb(ji, jj, &
&jk - 1, jp_tem)) / e3t_n(ji, jj, jk)
              tsn(ji, jj, jk, jp_sal) = (zdzp1 * tsb(ji, jj, jk + 1, jp_sal) + zdz * tsb(ji, jj, jk, jp_sal) + zdzm1 * tsb(ji, jj, &
&jk - 1, jp_sal)) / e3t_n(ji, jj, jk)
            END IF
          END DO
        END DO
      END DO
    END IF
    !$ACC KERNELS
    WHERE (tmask(:, :, :) == 1._wp .AND. tsn(:, :, :, 2) == 0._wp)
      tsn(:, :, :, 2) = - 99._wp
      tmask(:, :, :) = 0._wp
      umask(:, :, :) = 0._wp
      vmask(:, :, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    WHERE (SUM(tmask, dim = 3) == 0)
      mbkt(:, :) = 1
      mbku(:, :) = 1
      mbkv(:, :) = 1
      mikt(:, :) = 1
      miku(:, :) = 1
      mikv(:, :) = 1
    END WHERE
  END SUBROUTINE iscpl_rst_interpol
END MODULE iscplrst