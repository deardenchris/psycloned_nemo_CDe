MODULE domzgr
  USE oce
  USE dom_oce
  USE usrdef_zgr
  USE closea
  USE depth_e3
  USE wet_dry, ONLY: ll_wd, ssh_ref
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_zgr
  CONTAINS
  SUBROUTINE dom_zgr(k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(:, :), INTENT(OUT) :: k_top, k_bot
    INTEGER :: jk
    INTEGER :: ioptio, ibat, ios
    REAL(KIND = wp) :: zrefdep
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('dom_zgr', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dom_zgr : vertical coordinate'
      WRITE(numout, *) '~~~~~~~'
    END IF
    IF (ln_linssh .AND. lwp) WRITE(numout, *) '   linear free surface: the vertical mesh does not change in time'
    IF (ln_read_cfg) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   ==>>>   Read vertical mesh in ', TRIM(cn_domcfg), ' file'
      CALL zgr_read(ln_zco, ln_zps, ln_sco, ln_isfcav, gdept_1d, gdepw_1d, e3t_1d, e3w_1d, gdept_0, gdepw_0, e3t_0, e3u_0, e3v_0, e3f_0, e3w_0, e3uw_0, e3vw_0, k_top, k_bot)
    ELSE
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '          User defined vertical mesh (usr_def_zgr)'
      CALL usr_def_zgr(ln_zco, ln_zps, ln_sco, ln_isfcav, gdept_1d, gdepw_1d, e3t_1d, e3w_1d, gdept_0, gdepw_0, e3t_0, e3u_0, e3v_0, e3f_0, e3w_0, e3uw_0, e3vw_0, k_top, k_bot)
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    gde3w_0(:, :, 1) = 0.5_wp * e3w_0(:, :, 1)
    DO jk = 2, jpk
      gde3w_0(:, :, jk) = gde3w_0(:, :, jk - 1) + e3w_0(:, :, jk)
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dom_zgr', 'r1', 0, 0)
    IF (.NOT. ln_closea) CALL clo_bat(k_top, k_bot)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   Type of vertical coordinate (read in ', TRIM(cn_domcfg), ' file or set in userdef_nam) :'
      WRITE(numout, *) '      z-coordinate - full steps      ln_zco    = ', ln_zco
      WRITE(numout, *) '      z-coordinate - partial steps   ln_zps    = ', ln_zps
      WRITE(numout, *) '      s- or hybrid z-s-coordinate    ln_sco    = ', ln_sco
      WRITE(numout, *) '      ice shelf cavities             ln_isfcav = ', ln_isfcav
    END IF
    ioptio = 0
    IF (ln_zco) ioptio = ioptio + 1
    IF (ln_zps) ioptio = ioptio + 1
    IF (ln_sco) ioptio = ioptio + 1
    IF (ioptio /= 1) CALL ctl_stop(' none or several vertical coordinate options used')
    CALL zgr_top_bot(k_top, k_bot)
    zrefdep = 10._wp - 0.1_wp * MINVAL(e3w_1d)
    nlb10 = MINLOC(gdepw_1d, mask = gdepw_1d > zrefdep, dim = 1)
    nla10 = nlb10 - 1
    IF (nprint == 1 .AND. lwp) THEN
      WRITE(numout, *) ' MIN val k_top   ', MINVAL(k_top(:, :)), ' MAX ', MAXVAL(k_top(:, :))
      WRITE(numout, *) ' MIN val k_bot   ', MINVAL(k_bot(:, :)), ' MAX ', MAXVAL(k_bot(:, :))
      WRITE(numout, *) ' MIN val depth t ', MINVAL(gdept_0(:, :, :)), ' w ', MINVAL(gdepw_0(:, :, :)), '3w ', MINVAL(gde3w_0(:, :, :))
      WRITE(numout, *) ' MIN val e3    t ', MINVAL(e3t_0(:, :, :)), ' f ', MINVAL(e3f_0(:, :, :)), ' u ', MINVAL(e3u_0(:, :, :)), ' u ', MINVAL(e3v_0(:, :, :)), ' uw', MINVAL(e3uw_0(:, :, :)), ' vw', MINVAL(e3vw_0(:, :, :)), ' w ', MINVAL(e3w_0(:, :, :))
      WRITE(numout, *) ' MAX val depth t ', MAXVAL(gdept_0(:, :, :)), ' w ', MAXVAL(gdepw_0(:, :, :)), '3w ', MAXVAL(gde3w_0(:, :, :))
      WRITE(numout, *) ' MAX val e3    t ', MAXVAL(e3t_0(:, :, :)), ' f ', MAXVAL(e3f_0(:, :, :)), ' u ', MAXVAL(e3u_0(:, :, :)), ' u ', MAXVAL(e3v_0(:, :, :)), ' uw', MAXVAL(e3uw_0(:, :, :)), ' vw', MAXVAL(e3vw_0(:, :, :)), ' w ', MAXVAL(e3w_0(:, :, :))
    END IF
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE dom_zgr
  SUBROUTINE zgr_read(ld_zco, ld_zps, ld_sco, ld_isfcav, pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d, pdept, pdepw, pe3t, pe3u, pe3v, pe3f, pe3w, pe3uw, pe3vw, k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, INTENT(OUT) :: ld_zco, ld_zps, ld_sco
    LOGICAL, INTENT(OUT) :: ld_isfcav
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pdept_1d, pdepw_1d
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: pe3t_1d, pe3w_1d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pdept, pdepw
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3t, pe3u, pe3v, pe3f
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: pe3w, pe3uw, pe3vw
    INTEGER, DIMENSION(:, :), INTENT(OUT) :: k_top, k_bot
    INTEGER :: jk
    INTEGER :: inum
    REAL(KIND = WP) :: z_zco, z_zps, z_sco, z_cav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('zgr_read', 'r0', 0, 0)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   zgr_read : read the vertical coordinates in ', TRIM(cn_domcfg), ' file'
      WRITE(numout, *) '   ~~~~~~~~'
    END IF
    CALL iom_open(cn_domcfg, inum)
    CALL iom_get(inum, 'ln_zco', z_zco)
    CALL iom_get(inum, 'ln_zps', z_zps)
    CALL iom_get(inum, 'ln_sco', z_sco)
    IF (z_zco == 0._wp) THEN
      ld_zco = .FALSE.
    ELSE
      ld_zco = .TRUE.
    END IF
    IF (z_zps == 0._wp) THEN
      ld_zps = .FALSE.
    ELSE
      ld_zps = .TRUE.
    END IF
    IF (z_sco == 0._wp) THEN
      ld_sco = .FALSE.
    ELSE
      ld_sco = .TRUE.
    END IF
    CALL iom_get(inum, 'ln_isfcav', z_cav)
    IF (z_cav == 0._wp) THEN
      ld_isfcav = .FALSE.
    ELSE
      ld_isfcav = .TRUE.
    END IF
    CALL iom_get(inum, jpdom_unknown, 'e3t_1d', pe3t_1d)
    CALL iom_get(inum, jpdom_unknown, 'e3w_1d', pe3w_1d)
    CALL iom_get(inum, jpdom_data, 'e3t_0', pe3t, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e3u_0', pe3u, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e3v_0', pe3v, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e3f_0', pe3f, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e3w_0', pe3w, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e3uw_0', pe3uw, lrowattr = ln_use_jattr)
    CALL iom_get(inum, jpdom_data, 'e3vw_0', pe3vw, lrowattr = ln_use_jattr)
    IF (iom_varid(inum, 'gdept_1d', ldstop = .FALSE.) > 0 .AND. iom_varid(inum, 'gdepw_1d', ldstop = .FALSE.) > 0 .AND. iom_varid(inum, 'gdept_0', ldstop = .FALSE.) > 0 .AND. iom_varid(inum, 'gdepw_0', ldstop = .FALSE.) > 0) THEN
      CALL ctl_warn('zgr_read : old definition of depths and scale factors used ', '           depths at t- and w-points read in the domain configuration file')
      CALL iom_get(inum, jpdom_unknown, 'gdept_1d', pdept_1d)
      CALL iom_get(inum, jpdom_unknown, 'gdepw_1d', pdepw_1d)
      CALL iom_get(inum, jpdom_data, 'gdept_0', pdept, lrowattr = ln_use_jattr)
      CALL iom_get(inum, jpdom_data, 'gdepw_0', pdepw, lrowattr = ln_use_jattr)
    ELSE
      CALL e3_to_depth(pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d)
      CALL e3_to_depth(pe3t, pe3w, pdept, pdepw)
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '              Reference 1D z-coordinate depth and scale factors:'
        WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')")
        WRITE(numout, "(10x, i4, 4f9.2)") (jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk)
      END IF
    END IF
    CALL iom_get(inum, jpdom_data, 'top_level', z2d, lrowattr = ln_use_jattr)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    k_top(:, :) = NINT(z2d(:, :))
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('zgr_read', 'r1', 0, 0)
    CALL iom_get(inum, jpdom_data, 'bottom_level', z2d, lrowattr = ln_use_jattr)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    k_bot(:, :) = NINT(z2d(:, :))
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('zgr_read', 'r2', 0, 0)
    IF (ll_wd) CALL iom_get(inum, 'rn_wd_ref_depth', ssh_ref)
    CALL iom_close(inum)
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE zgr_read
  SUBROUTINE zgr_top_bot(k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(:, :), INTENT(IN) :: k_top, k_bot
    INTEGER :: ji, jj
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    CALL profile_psy_data0 % PreStart('zgr_top_bot', 'r0', 0, 0)
    IF (lwp) WRITE(numout, *)
    IF (lwp) WRITE(numout, *) '    zgr_top_bot : ocean top and bottom k-index of T-, U-, V- and W-levels '
    IF (lwp) WRITE(numout, *) '    ~~~~~~~~~~~'
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    mikt(:, :) = MAX(k_top(:, :), 1)
    mbkt(:, :) = MAX(k_bot(:, :), 1)
    !$ACC loop independent collapse(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        miku(ji, jj) = MAX(mikt(ji + 1, jj), mikt(ji, jj))
        mikv(ji, jj) = MAX(mikt(ji, jj + 1), mikt(ji, jj))
        mikf(ji, jj) = MAX(mikt(ji, jj + 1), mikt(ji, jj), mikt(ji + 1, jj), mikt(ji + 1, jj + 1))
        mbku(ji, jj) = MIN(mbkt(ji + 1, jj), mbkt(ji, jj))
        mbkv(ji, jj) = MIN(mbkt(ji, jj + 1), mbkt(ji, jj))
      END DO
    END DO
    zk(:, :) = REAL(miku(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('zgr_top_bot', 'r1', 0, 0)
    CALL lbc_lnk('domzgr', zk, 'U', 1.)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    miku(:, :) = MAX(NINT(zk(:, :)), 1)
    zk(:, :) = REAL(mikv(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('zgr_top_bot', 'r2', 0, 0)
    CALL lbc_lnk('domzgr', zk, 'V', 1.)
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    mikv(:, :) = MAX(NINT(zk(:, :)), 1)
    zk(:, :) = REAL(mikf(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('zgr_top_bot', 'r3', 0, 0)
    CALL lbc_lnk('domzgr', zk, 'F', 1.)
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    mikf(:, :) = MAX(NINT(zk(:, :)), 1)
    zk(:, :) = REAL(mbku(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data4 % PreStart('zgr_top_bot', 'r4', 0, 0)
    CALL lbc_lnk('domzgr', zk, 'U', 1.)
    CALL profile_psy_data4 % PostEnd
    !$ACC KERNELS
    mbku(:, :) = MAX(NINT(zk(:, :)), 1)
    zk(:, :) = REAL(mbkv(:, :), wp)
    !$ACC END KERNELS
    CALL profile_psy_data5 % PreStart('zgr_top_bot', 'r5', 0, 0)
    CALL lbc_lnk('domzgr', zk, 'V', 1.)
    CALL profile_psy_data5 % PostEnd
    !$ACC KERNELS
    mbkv(:, :) = MAX(NINT(zk(:, :)), 1)
    !$ACC END KERNELS
  END SUBROUTINE zgr_top_bot
END MODULE domzgr