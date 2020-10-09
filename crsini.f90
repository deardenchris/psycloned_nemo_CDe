MODULE crsini
  USE par_kind, ONLY: wp
  USE par_oce
  USE dom_oce
  USE crs
  USE phycst, ONLY: omega, rad
  USE crsdom
  USE crsdomwri
  USE crslbclnk
  USE iom
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: crs_init
  CONTAINS
  SUBROUTINE crs_init
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    INTEGER :: ios
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: ze3t, ze3u, ze3v, ze3w
    NAMELIST /namcrs/ nn_factx, nn_facty, nn_binref, ln_msh_crs, nn_crs_kz, ln_crs_wn
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namcrs, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namcrs in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namcrs, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namcrs in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namcrs)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'crs_init : Initializing the grid coarsening module'
      WRITE(numout, FMT = *) '~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namcrs '
      WRITE(numout, FMT = *) '      coarsening factor in i-direction      nn_factx   = ', nn_factx
      WRITE(numout, FMT = *) '      coarsening factor in j-direction      nn_facty   = ', nn_facty
      WRITE(numout, FMT = *) '      bin centering preference              nn_binref  = ', nn_binref
      WRITE(numout, FMT = *) '      create a mesh file (=T)               ln_msh_crs = ', ln_msh_crs
      WRITE(numout, FMT = *) '      type of Kz coarsening (0,1,2)         nn_crs_kz  = ', nn_crs_kz
      WRITE(numout, FMT = *) '      wn coarsened or computed using hdivn  ln_crs_wn  = ', ln_crs_wn
    END IF
    rfactx_r = 1. / nn_factx
    rfacty_r = 1. / nn_facty
    CALL crs_dom_def
    CALL crs_dom_msk
    IF (nresty /= 0 .AND. nrestx /= 0) THEN
      CALL crs_dom_coordinates(gphit, glamt, 'T', gphit_crs, glamt_crs)
      CALL crs_dom_coordinates(gphiu, glamu, 'U', gphiu_crs, glamu_crs)
      CALL crs_dom_coordinates(gphiv, glamv, 'V', gphiv_crs, glamv_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'F', gphif_crs, glamf_crs)
    ELSE IF (nresty /= 0 .AND. nrestx == 0) THEN
      CALL crs_dom_coordinates(gphiu, glamu, 'T', gphit_crs, glamt_crs)
      CALL crs_dom_coordinates(gphiu, glamu, 'U', gphiu_crs, glamu_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'V', gphiv_crs, glamv_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'F', gphif_crs, glamf_crs)
    ELSE IF (nresty == 0 .AND. nrestx /= 0) THEN
      CALL crs_dom_coordinates(gphiv, glamv, 'T', gphit_crs, glamt_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'U', gphiu_crs, glamu_crs)
      CALL crs_dom_coordinates(gphiv, glamv, 'V', gphiv_crs, glamv_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'F', gphif_crs, glamf_crs)
    ELSE
      CALL crs_dom_coordinates(gphif, glamf, 'T', gphit_crs, glamt_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'U', gphiu_crs, glamu_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'V', gphiv_crs, glamv_crs)
      CALL crs_dom_coordinates(gphif, glamf, 'F', gphif_crs, glamf_crs)
    END IF
    CALL crs_dom_hgr(e1t, e2t, 'T', e1t_crs, e2t_crs)
    CALL crs_dom_hgr(e1u, e2u, 'U', e1u_crs, e2u_crs)
    CALL crs_dom_hgr(e1v, e2v, 'V', e1v_crs, e2v_crs)
    CALL crs_dom_hgr(e1f, e2f, 'F', e1f_crs, e2f_crs)
    !$ACC KERNELS
    e1e2t_crs(:, :) = e1t_crs(:, :) * e2t_crs(:, :)
    ff_crs(:, :) = 2. * omega * SIN(rad * gphif_crs(:, :))
    !$ACC END KERNELS
    CALL ctl_warn('crsini: CAUTION, CRS only designed for Coriolis defined on the sphere')
    CALL crs_dom_bat
    !$ACC KERNELS
    ze3t(:, :, :) = e3t_n(:, :, :)
    ze3u(:, :, :) = e3u_n(:, :, :)
    ze3v(:, :, :) = e3v_n(:, :, :)
    ze3w(:, :, :) = e3w_n(:, :, :)
    !$ACC END KERNELS
    CALL crs_dom_sfc(tmask, 'W', e1e2w_crs, e1e2w_msk, p_e1 = e1t, p_e2 = e2t)
    CALL crs_dom_sfc(umask, 'U', e2e3u_crs, e2e3u_msk, p_e2 = e2u, p_e3 = ze3u)
    CALL crs_dom_sfc(vmask, 'V', e1e3v_crs, e1e3v_msk, p_e1 = e1v, p_e3 = ze3v)
    !$ACC KERNELS
    facsurfu(:, :, :) = umask_crs(:, :, :) * e2e3u_msk(:, :, :) / e2e3u_crs(:, :, :)
    facsurfv(:, :, :) = vmask_crs(:, :, :) * e1e3v_msk(:, :, :) / e1e3v_crs(:, :, :)
    !$ACC END KERNELS
    CALL crs_dom_e3(e1t, e2t, ze3t, e1e2w_crs, 'T', tmask, e3t_crs, e3t_max_crs)
    CALL crs_dom_e3(e1u, e2u, ze3u, e2e3u_crs, 'U', umask, e3u_crs, e3u_max_crs)
    CALL crs_dom_e3(e1v, e2v, ze3v, e1e3v_crs, 'V', vmask, e3v_crs, e3v_max_crs)
    CALL crs_dom_e3(e1t, e2t, ze3w, e1e2w_crs, 'W', tmask, e3w_crs, e3w_max_crs)
    !$ACC KERNELS
    DO jk = 1, jpk
      DO ji = 1, jpi_crs
        DO jj = 1, jpj_crs
          IF (e3t_crs(ji, jj, jk) == 0._wp) e3t_crs(ji, jj, jk) = e3t_1d(jk)
          IF (e3w_crs(ji, jj, jk) == 0._wp) e3w_crs(ji, jj, jk) = e3w_1d(jk)
          IF (e3u_crs(ji, jj, jk) == 0._wp) e3u_crs(ji, jj, jk) = e3t_1d(jk)
          IF (e3v_crs(ji, jj, jk) == 0._wp) e3v_crs(ji, jj, jk) = e3t_1d(jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL crs_dom_ope(gdept_0, 'MAX', 'T', tmask, gdept_crs, p_e3 = ze3t, psgn = 1.0)
    CALL crs_dom_ope(gdepw_0, 'MAX', 'W', tmask, gdepw_crs, p_e3 = ze3w, psgn = 1.0)
    CALL crs_dom_facvol(tmask, 'T', e1t, e2t, ze3t, ocean_volume_crs_t, facvol_t)
    !$ACC KERNELS
    bt_crs(:, :, :) = ocean_volume_crs_t(:, :, :) * facvol_t(:, :, :)
    r1_bt_crs(:, :, :) = 0._wp
    !$ACC END KERNELS
    WHERE (bt_crs /= 0._wp) r1_bt_crs(:, :, :) = 1._wp / bt_crs(:, :, :)
    CALL crs_dom_facvol(tmask, 'W', e1t, e2t, ze3w, ocean_volume_crs_w, facvol_w)
    IF (ln_msh_crs) THEN
      CALL dom_grid_crs
      CALL crs_dom_wri
      CALL dom_grid_glo
    END IF
  END SUBROUTINE crs_init
END MODULE crsini