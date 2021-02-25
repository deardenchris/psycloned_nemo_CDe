MODULE crsdomwri
  USE timing
  USE dom_oce
  USE in_out_manager
  USE par_kind, ONLY: wp
  USE lib_mpp
  USE iom_def
  USE iom
  USE crs
  USE crsdom
  USE crslbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: crs_dom_wri
  CONTAINS
  SUBROUTINE crs_dom_wri
    INTEGER :: ji, jj, jk
    INTEGER :: inum
    INTEGER :: iif, iil, ijf, ijl
    CHARACTER(LEN = 21) :: clnam
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs) :: zprt, zprw
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk) :: zdepu, zdepv
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'crs_dom_wri : create NetCDF mesh and mask file'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~'
    clnam = 'mesh_mask_crs'
    CALL iom_open(TRIM(clnam), inum, ldwrt = .TRUE.)
    CALL iom_rstput(0, 0, inum, 'tmask', tmask_crs, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'umask', umask_crs, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'vmask', vmask_crs, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'fmask', fmask_crs, ktype = jp_i1)
    tmask_i_crs(:, :) = tmask_crs(:, :, 1)
    iif = nn_hls
    iil = nlci_crs - nn_hls + 1
    ijf = nn_hls
    ijl = nlcj_crs - nn_hls + 1
    tmask_i_crs(1 : iif, :) = 0._wp
    tmask_i_crs(iil : jpi_crs, :) = 0._wp
    tmask_i_crs(:, 1 : ijf) = 0._wp
    tmask_i_crs(:, ijl : jpj_crs) = 0._wp
    tpol_crs(1 : jpiglo_crs, :) = 1._wp
    fpol_crs(1 : jpiglo_crs, :) = 1._wp
    IF (jperio == 3 .OR. jperio == 4) THEN
      tpol_crs(jpiglo_crs / 2 + 1 : jpiglo_crs, :) = 0._wp
      fpol_crs(1 : jpiglo_crs, :) = 0._wp
      IF (mjg_crs(nlej_crs) == jpiglo_crs) THEN
        DO ji = iif + 1, iil - 1
          tmask_i_crs(ji, nlej_crs - 1) = tmask_i_crs(ji, nlej_crs - 1) * tpol_crs(mig_crs(ji), 1)
        END DO
      END IF
    END IF
    IF (jperio == 5 .OR. jperio == 6) THEN
      tpol_crs(1 : jpiglo_crs, :) = 0._wp
      fpol_crs(jpiglo_crs / 2 + 1 : jpiglo_crs, :) = 0._wp
    END IF
    CALL iom_rstput(0, 0, inum, 'tmaskutil', tmask_i_crs, ktype = jp_i1)
    CALL dom_uniq_crs(zprw, 'U')
    zprt = umask_crs(:, :, 1) * zprw
    CALL iom_rstput(0, 0, inum, 'umaskutil', zprt, ktype = jp_i1)
    CALL dom_uniq_crs(zprw, 'V')
    zprt = vmask_crs(:, :, 1) * zprw
    CALL iom_rstput(0, 0, inum, 'vmaskutil', zprt, ktype = jp_i1)
    CALL dom_uniq_crs(zprw, 'F')
    zprt = fmask_crs(:, :, 1) * zprw
    CALL iom_rstput(0, 0, inum, 'fmaskutil', zprt, ktype = jp_i1)
    CALL iom_rstput(0, 0, inum, 'glamt', glamt_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'glamu', glamu_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'glamv', glamv_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'glamf', glamf_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gphit', gphit_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gphiu', gphiu_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gphiv', gphiv_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gphif', gphif_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'e1t', e1t_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1u', e1u_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1v', e1v_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1f', e1f_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2t', e2t_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2u', e2u_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2v', e2v_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2f', e2f_crs, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'ff', ff_crs, ktype = jp_r8)
    zprt(:, :) = tmask_crs(:, :, 1) * REAL(mbkt_crs(:, :), wp)
    CALL iom_rstput(0, 0, inum, 'mbathy', zprt, ktype = jp_i2)
    CALL iom_rstput(0, 0, inum, 'e3t', e3t_crs)
    CALL iom_rstput(0, 0, inum, 'e3w', e3w_crs)
    CALL iom_rstput(0, 0, inum, 'e3u', e3u_crs)
    CALL iom_rstput(0, 0, inum, 'e3v', e3v_crs)
    CALL iom_rstput(0, 0, inum, 'gdept', gdept_crs, ktype = jp_r4)
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpk
      DO jj = 1, jpj_crsm1
        DO ji = 1, jpi_crsm1
          zdepu(ji, jj, jk) = MIN(gdept_crs(ji, jj, jk), gdept_crs(ji + 1, jj, jk)) * umask_crs(ji, jj, jk)
          zdepv(ji, jj, jk) = MIN(gdept_crs(ji, jj, jk), gdept_crs(ji, jj + 1, jk)) * vmask_crs(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL crs_lbc_lnk(zdepu, 'U', 1.)
    CALL crs_lbc_lnk(zdepv, 'V', 1.)
    CALL iom_rstput(0, 0, inum, 'gdepu', zdepu, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gdepv', zdepv, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gdepw', gdepw_crs, ktype = jp_r4)
    CALL iom_rstput(0, 0, inum, 'gdept_1d', gdept_1d)
    CALL iom_rstput(0, 0, inum, 'gdepw_1d', gdepw_1d)
    CALL iom_rstput(0, 0, inum, 'e3t_1d', e3t_1d)
    CALL iom_rstput(0, 0, inum, 'e3w_1d', e3w_1d)
    CALL iom_rstput(0, 0, inum, 'ocean_volume_t', ocean_volume_crs_t)
    CALL iom_rstput(0, 0, inum, 'facvol_t', facvol_t)
    CALL iom_rstput(0, 0, inum, 'facvol_w', facvol_w)
    CALL iom_rstput(0, 0, inum, 'facsurfu', facsurfu)
    CALL iom_rstput(0, 0, inum, 'facsurfv', facsurfv)
    CALL iom_rstput(0, 0, inum, 'e1e2w_msk', e1e2w_msk)
    CALL iom_rstput(0, 0, inum, 'e2e3u_msk', e2e3u_msk)
    CALL iom_rstput(0, 0, inum, 'e1e3v_msk', e1e3v_msk)
    CALL iom_rstput(0, 0, inum, 'e1e2w', e1e2w_crs)
    CALL iom_rstput(0, 0, inum, 'e2e3u', e2e3u_crs)
    CALL iom_rstput(0, 0, inum, 'e1e3v', e1e3v_crs)
    CALL iom_rstput(0, 0, inum, 'bt', bt_crs)
    CALL iom_rstput(0, 0, inum, 'r1_bt', r1_bt_crs)
    CALL iom_rstput(0, 0, inum, 'crs_surfu_wgt', crs_surfu_wgt)
    CALL iom_rstput(0, 0, inum, 'crs_surfv_wgt', crs_surfv_wgt)
    CALL iom_rstput(0, 0, inum, 'crs_volt_wgt', crs_volt_wgt)
    CALL iom_close(inum)
  END SUBROUTINE crs_dom_wri
  SUBROUTINE dom_uniq_crs(puniq, cdgrd)
    CHARACTER(LEN = 1), INTENT(IN) :: cdgrd
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: puniq
    REAL(KIND = wp) :: zshift
    INTEGER :: ji
    LOGICAL, DIMENSION(SIZE(puniq, 1), SIZE(puniq, 2), 1) :: lldbl
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs) :: ztstref
    zshift = jpi_crs * jpj_crs * (narea - 1)
    ztstref(:, :) = RESHAPE((/(zshift + REAL(ji, wp), ji = 1, jpi_crs * jpj_crs)/), (/jpi_crs, jpj_crs/))
    puniq(:, :) = ztstref(:, :)
    CALL crs_lbc_lnk(puniq, cdgrd, 1.)
    lldbl(:, :, 1) = puniq(:, :) == ztstref(:, :)
    puniq(:, :) = 1.
    puniq(nldi_crs : nlei_crs, nldj_crs : nlej_crs) = REAL(COUNT(lldbl(nldi_crs : nlei_crs, nldj_crs : nlej_crs, :), dim = 3), wp)
  END SUBROUTINE dom_uniq_crs
END MODULE crsdomwri