MODULE icbrst
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  USE netcdf
  USE iom
  USE ioipsl, ONLY: ju2ymds
  USE icb_oce
  USE icbutl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_rst_read
  PUBLIC :: icb_rst_write
  INTEGER :: nlonid, nlatid, nxid, nyid, nuvelid, nvvelid
  INTEGER :: nmassid, nthicknessid, nwidthid, nlengthid
  INTEGER :: nyearid, ndayid
  INTEGER :: nscaling_id, nmass_of_bits_id, nheat_density_id, numberid
  INTEGER :: nsiceid, nsheatid, ncalvid, ncalvhid, nkountid
  INTEGER :: nret, ncid, nc_dim
  INTEGER, DIMENSION(3) :: nstrt3, nlngth3
  CONTAINS
  SUBROUTINE icb_rst_read
    INTEGER :: idim, ivar, iatt
    INTEGER :: jn, iunlim_dim, ibergs_in_file
    INTEGER :: ii, ij, iclass, ibase_err, imax_icb
    REAL(KIND = wp), DIMENSION(nkounts) :: zdata
    LOGICAL :: ll_found_restart
    CHARACTER(LEN = 256) :: cl_path
    CHARACTER(LEN = 256) :: cl_filename
    CHARACTER(LEN = NF90_MAX_NAME) :: cl_dname
    TYPE(iceberg) :: localberg
    TYPE(point) :: localpt
    cl_path = TRIM(cn_ocerst_indir)
    IF (cl_path(LEN_TRIM(cl_path) :) /= '/') cl_path = TRIM(cl_path) // '/'
    cl_filename = TRIM(cn_ocerst_in) // '_icebergs'
    CALL iom_open(TRIM(cl_path) // cl_filename, ncid)
    imax_icb = 0
    IF (iom_file(ncid) % iduld .GE. 0) THEN
      ibergs_in_file = iom_file(ncid) % lenuld
      DO jn = 1, ibergs_in_file
        CALL iom_get(ncid, 'xi', localpt % xi, ktime = jn)
        CALL iom_get(ncid, 'yj', localpt % yj, ktime = jn)
        ii = INT(localpt % xi + 0.5)
        ij = INT(localpt % yj + 0.5)
        IF (ii .GE. nldi + nimpp - 1 .AND. ii .LE. nlei + nimpp - 1 .AND. ij .GE. nldj + njmpp - 1 .AND. ij .LE. nlej + njmpp - 1) &
&THEN
          CALL iom_get(ncid, jpdom_unknown, 'number', zdata(:), ktime = jn, kstart = (/1/), kcount = (/nkounts/))
          localberg % number(:) = INT(zdata(:))
          imax_icb = MAX(imax_icb, INT(zdata(1)))
          CALL iom_get(ncid, 'mass_scaling', localberg % mass_scaling, ktime = jn)
          CALL iom_get(ncid, 'lon', localpt % lon, ktime = jn)
          CALL iom_get(ncid, 'lat', localpt % lat, ktime = jn)
          CALL iom_get(ncid, 'uvel', localpt % uvel, ktime = jn)
          CALL iom_get(ncid, 'vvel', localpt % vvel, ktime = jn)
          CALL iom_get(ncid, 'mass', localpt % mass, ktime = jn)
          CALL iom_get(ncid, 'thickness', localpt % thickness, ktime = jn)
          CALL iom_get(ncid, 'width', localpt % width, ktime = jn)
          CALL iom_get(ncid, 'length', localpt % length, ktime = jn)
          CALL iom_get(ncid, 'year', zdata(1), ktime = jn)
          localpt % year = INT(zdata(1))
          CALL iom_get(ncid, 'day', localpt % day, ktime = jn)
          CALL iom_get(ncid, 'mass_of_bits', localpt % mass_of_bits, ktime = jn)
          CALL iom_get(ncid, 'heat_density', localpt % heat_density, ktime = jn)
          CALL icb_utl_add(localberg, localpt)
        END IF
      END DO
    ELSE
      ibergs_in_file = 0
    END IF
    CALL iom_get(ncid, jpdom_autoglo, 'calving', src_calving)
    CALL iom_get(ncid, jpdom_autoglo, 'calving_hflx', src_calving_hflx)
    CALL iom_get(ncid, jpdom_autoglo, 'stored_heat', berg_grid(1) % stored_heat)
    CALL iom_get(ncid, jpdom_autoglo_xy, 'stored_ice', berg_grid(1) % stored_ice, kstart = (/1, 1, 1/), kcount = (/1, 1, nclasses/))
    CALL iom_get(ncid, jpdom_unknown, 'kount', zdata(:))
    num_bergs(:) = INT(zdata(:))
    CALL iom_close(ncid)
    jn = icb_utl_count()
    IF (lwp .AND. nn_verbose_level >= 0) WRITE(numout, FMT = '(2(a,i5))') 'icebergs, read_restart_bergs: # bergs =', jn, ' on PE', &
&narea - 1
    IF (lk_mpp) THEN
      IF (INDEX(iom_file(ncid) % name, 'icebergs.nc') .EQ. 0) CALL mpp_sum('icbrst', ibergs_in_file)
      CALL mpp_sum('icbrst', jn)
    END IF
    IF (lwp) WRITE(numout, FMT = '(a,i5,a,i5,a)') 'icebergs, icb_rst_read: there were', ibergs_in_file, ' bergs in the restart &
&file and', jn, ' bergs have been read'
    ibase_err = 0
    IF (num_bergs(1) < 0 .AND. num_bergs(1) /= narea - jpnij) THEN
      ibase_err = 1
    ELSE IF (MOD(num_bergs(1) - narea, jpnij) /= 0) THEN
      ibase_err = 1
    END IF
    IF (lk_mpp) THEN
      CALL mpp_sum('icbrst', ibase_err)
    END IF
    IF (ibase_err > 0) THEN
      IF (lk_mpp) THEN
        CALL mpp_max('icbrst', imax_icb)
      END IF
      num_bergs(1) = imax_icb - jpnij + narea
    END IF
    IF (lwp .AND. nn_verbose_level >= 0) WRITE(numout, FMT = '(a)') 'icebergs, icb_rst_read: completed'
  END SUBROUTINE icb_rst_read
  SUBROUTINE icb_rst_write(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jn
    INTEGER :: ix_dim, iy_dim, ik_dim, in_dim
    INTEGER :: iyear, imonth, iday
    REAL(KIND = wp) :: zsec
    REAL(KIND = wp) :: zfjulday
    CHARACTER(LEN = 256) :: cl_path
    CHARACTER(LEN = 256) :: cl_filename
    CHARACTER(LEN = 20) :: clkt
    TYPE(iceberg), POINTER :: this
    TYPE(point), POINTER :: pt
    RETURN
    IF (kt == nitrst) THEN
      cl_path = TRIM(cn_ocerst_outdir)
      IF (cl_path(LEN_TRIM(cl_path) :) /= '/') cl_path = TRIM(cl_path) // '/'
      IF (ln_rstdate) THEN
        zfjulday = fjulday + rdt / rday
        IF (ABS(zfjulday - REAL(NINT(zfjulday), wp)) < 0.1 / rday) zfjulday = REAL(NINT(zfjulday), wp)
        CALL ju2ymds(zfjulday, iyear, imonth, iday, zsec)
        WRITE(clkt, FMT = '(i4.4,2i2.2)') iyear, imonth, iday
      ELSE
        IF (kt > 999999999) THEN
          WRITE(clkt, FMT = *) kt
        ELSE
          WRITE(clkt, FMT = '(i8.8)') kt
        END IF
      END IF
      IF (lk_mpp) THEN
        WRITE(cl_filename, FMT = '(A,"_icebergs_",A,"_restart_",I6.6,".nc")') TRIM(cexper), TRIM(ADJUSTL(clkt)), narea - 1
      ELSE
        WRITE(cl_filename, FMT = '(A,"_icebergs_",A,"_restart.nc")') TRIM(cexper), TRIM(ADJUSTL(clkt))
      END IF
      IF (lwp .AND. nn_verbose_level >= 0) WRITE(numout, FMT = '(2a)') 'icebergs, write_restart: creating ', TRIM(cl_path) // &
&TRIM(cl_filename)
      nret = NF90_CREATE(TRIM(cl_path) // TRIM(cl_filename), NF90_CLOBBER, ncid)
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_create failed')
      nret = NF90_DEF_DIM(ncid, 'x', jpi, ix_dim)
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim x failed')
      nret = NF90_DEF_DIM(ncid, 'y', jpj, iy_dim)
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim y failed')
      nret = NF90_DEF_DIM(ncid, 'c', nclasses, nc_dim)
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim c failed')
      nret = NF90_DEF_DIM(ncid, 'k', nkounts, ik_dim)
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim k failed')
      IF (lk_mpp) THEN
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_number_total', jpnij)
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_number', narea - 1)
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', (/1, 2/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_size_global', (/jpiglo, jpjglo/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_size_local', (/jpi, jpj/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_position_first', (/nimpp, njmpp/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_position_last', (/nimpp + jpi - 1, njmpp + jpj - 1/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', (/nldi - 1, nldj - 1/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', (/jpi - nlei, jpj - nlej/))
        nret = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'DOMAIN_type', 'BOX')
      END IF
      IF (ASSOCIATED(first_berg)) THEN
        nret = NF90_DEF_DIM(ncid, 'n', NF90_UNLIMITED, in_dim)
        IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_def_dim n failed')
      END IF
      nret = NF90_DEF_VAR(ncid, 'kount', NF90_INT, (/ik_dim/), nkountid)
      nret = NF90_DEF_VAR(ncid, 'calving', NF90_DOUBLE, (/ix_dim, iy_dim/), ncalvid)
      nret = NF90_DEF_VAR(ncid, 'calving_hflx', NF90_DOUBLE, (/ix_dim, iy_dim/), ncalvhid)
      nret = NF90_DEF_VAR(ncid, 'stored_ice', NF90_DOUBLE, (/ix_dim, iy_dim, nc_dim/), nsiceid)
      nret = NF90_DEF_VAR(ncid, 'stored_heat', NF90_DOUBLE, (/ix_dim, iy_dim/), nsheatid)
      nret = NF90_PUT_ATT(ncid, ncalvid, 'long_name', 'iceberg calving')
      nret = NF90_PUT_ATT(ncid, ncalvid, 'units', 'some')
      nret = NF90_PUT_ATT(ncid, ncalvhid, 'long_name', 'heat flux associated with iceberg calving')
      nret = NF90_PUT_ATT(ncid, ncalvhid, 'units', 'some')
      nret = NF90_PUT_ATT(ncid, nsiceid, 'long_name', 'stored ice used to calve icebergs')
      nret = NF90_PUT_ATT(ncid, nsiceid, 'units', 'kg/s')
      nret = NF90_PUT_ATT(ncid, nsheatid, 'long_name', 'heat in stored ice used to calve icebergs')
      nret = NF90_PUT_ATT(ncid, nsheatid, 'units', 'J/kg/s')
      IF (ASSOCIATED(first_berg)) THEN
        nret = NF90_DEF_VAR(ncid, 'lon', NF90_DOUBLE, in_dim, nlonid)
        nret = NF90_DEF_VAR(ncid, 'lat', NF90_DOUBLE, in_dim, nlatid)
        nret = NF90_DEF_VAR(ncid, 'xi', NF90_DOUBLE, in_dim, nxid)
        nret = NF90_DEF_VAR(ncid, 'yj', NF90_DOUBLE, in_dim, nyid)
        nret = NF90_DEF_VAR(ncid, 'uvel', NF90_DOUBLE, in_dim, nuvelid)
        nret = NF90_DEF_VAR(ncid, 'vvel', NF90_DOUBLE, in_dim, nvvelid)
        nret = NF90_DEF_VAR(ncid, 'mass', NF90_DOUBLE, in_dim, nmassid)
        nret = NF90_DEF_VAR(ncid, 'thickness', NF90_DOUBLE, in_dim, nthicknessid)
        nret = NF90_DEF_VAR(ncid, 'width', NF90_DOUBLE, in_dim, nwidthid)
        nret = NF90_DEF_VAR(ncid, 'length', NF90_DOUBLE, in_dim, nlengthid)
        nret = NF90_DEF_VAR(ncid, 'number', NF90_INT, (/ik_dim, in_dim/), numberid)
        nret = NF90_DEF_VAR(ncid, 'year', NF90_INT, in_dim, nyearid)
        nret = NF90_DEF_VAR(ncid, 'day', NF90_DOUBLE, in_dim, ndayid)
        nret = NF90_DEF_VAR(ncid, 'mass_scaling', NF90_DOUBLE, in_dim, nscaling_id)
        nret = NF90_DEF_VAR(ncid, 'mass_of_bits', NF90_DOUBLE, in_dim, nmass_of_bits_id)
        nret = NF90_DEF_VAR(ncid, 'heat_density', NF90_DOUBLE, in_dim, nheat_density_id)
        nret = NF90_PUT_ATT(ncid, nlonid, 'long_name', 'longitude')
        nret = NF90_PUT_ATT(ncid, nlonid, 'units', 'degrees_E')
        nret = NF90_PUT_ATT(ncid, nlatid, 'long_name', 'latitude')
        nret = NF90_PUT_ATT(ncid, nlatid, 'units', 'degrees_N')
        nret = NF90_PUT_ATT(ncid, nxid, 'long_name', 'x grid box position')
        nret = NF90_PUT_ATT(ncid, nxid, 'units', 'fractional')
        nret = NF90_PUT_ATT(ncid, nyid, 'long_name', 'y grid box position')
        nret = NF90_PUT_ATT(ncid, nyid, 'units', 'fractional')
        nret = NF90_PUT_ATT(ncid, nuvelid, 'long_name', 'zonal velocity')
        nret = NF90_PUT_ATT(ncid, nuvelid, 'units', 'm/s')
        nret = NF90_PUT_ATT(ncid, nvvelid, 'long_name', 'meridional velocity')
        nret = NF90_PUT_ATT(ncid, nvvelid, 'units', 'm/s')
        nret = NF90_PUT_ATT(ncid, nmassid, 'long_name', 'mass')
        nret = NF90_PUT_ATT(ncid, nmassid, 'units', 'kg')
        nret = NF90_PUT_ATT(ncid, nthicknessid, 'long_name', 'thickness')
        nret = NF90_PUT_ATT(ncid, nthicknessid, 'units', 'm')
        nret = NF90_PUT_ATT(ncid, nwidthid, 'long_name', 'width')
        nret = NF90_PUT_ATT(ncid, nwidthid, 'units', 'm')
        nret = NF90_PUT_ATT(ncid, nlengthid, 'long_name', 'length')
        nret = NF90_PUT_ATT(ncid, nlengthid, 'units', 'm')
        nret = NF90_PUT_ATT(ncid, numberid, 'long_name', 'iceberg number on this processor')
        nret = NF90_PUT_ATT(ncid, numberid, 'units', 'count')
        nret = NF90_PUT_ATT(ncid, nyearid, 'long_name', 'calendar year of calving event')
        nret = NF90_PUT_ATT(ncid, nyearid, 'units', 'years')
        nret = NF90_PUT_ATT(ncid, ndayid, 'long_name', 'year day of calving event')
        nret = NF90_PUT_ATT(ncid, ndayid, 'units', 'days')
        nret = NF90_PUT_ATT(ncid, nscaling_id, 'long_name', 'scaling factor for mass of calving berg')
        nret = NF90_PUT_ATT(ncid, nscaling_id, 'units', 'none')
        nret = NF90_PUT_ATT(ncid, nmass_of_bits_id, 'long_name', 'mass of bergy bits')
        nret = NF90_PUT_ATT(ncid, nmass_of_bits_id, 'units', 'kg')
        nret = NF90_PUT_ATT(ncid, nheat_density_id, 'long_name', 'heat density')
        nret = NF90_PUT_ATT(ncid, nheat_density_id, 'units', 'J/kg')
      END IF
      nret = NF90_ENDDEF(ncid)
      nstrt3(1) = 1
      nstrt3(2) = 1
      nlngth3(1) = jpi
      nlngth3(2) = jpj
      nlngth3(3) = 1
      DO jn = 1, nclasses
        griddata(:, :, 1) = berg_grid(1) % stored_ice(:, :, jn)
        nstrt3(3) = jn
        nret = NF90_PUT_VAR(ncid, nsiceid, griddata, nstrt3, nlngth3)
        IF (nret .NE. NF90_NOERR) THEN
          IF (lwp) WRITE(numout, FMT = *) TRIM(NF90_STRERROR(nret))
          CALL ctl_stop('icebergs, write_restart: nf_put_var stored_ice failed')
        END IF
      END DO
      IF (lwp) WRITE(numout, FMT = *) 'file: ', TRIM(cl_path) // TRIM(cl_filename), ' var: stored_ice  written'
      nret = NF90_PUT_VAR(ncid, nkountid, num_bergs(:))
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var kount failed')
      nret = NF90_PUT_VAR(ncid, nsheatid, berg_grid(1) % stored_heat(:, :))
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var stored_heat failed')
      IF (lwp) WRITE(numout, FMT = *) 'file: ', TRIM(cl_path) // TRIM(cl_filename), ' var: stored_heat written'
      nret = NF90_PUT_VAR(ncid, ncalvid, src_calving(:, :))
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var calving failed')
      nret = NF90_PUT_VAR(ncid, ncalvhid, src_calving_hflx(:, :))
      IF (nret .NE. NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_put_var calving_hflx failed')
      IF (lwp) WRITE(numout, FMT = *) 'file: ', TRIM(cl_path) // TRIM(cl_filename), ' var: calving written'
      IF (ASSOCIATED(first_berg)) THEN
        this => first_berg
        jn = 0
        DO WHILE (ASSOCIATED(this))
          pt => this % current_point
          jn = jn + 1
          nret = NF90_PUT_VAR(ncid, numberid, this % number, (/1, jn/), (/nkounts, 1/))
          nret = NF90_PUT_VAR(ncid, nscaling_id, this % mass_scaling, (/jn/))
          nret = NF90_PUT_VAR(ncid, nlonid, pt % lon, (/jn/))
          nret = NF90_PUT_VAR(ncid, nlatid, pt % lat, (/jn/))
          nret = NF90_PUT_VAR(ncid, nxid, pt % xi, (/jn/))
          nret = NF90_PUT_VAR(ncid, nyid, pt % yj, (/jn/))
          nret = NF90_PUT_VAR(ncid, nuvelid, pt % uvel, (/jn/))
          nret = NF90_PUT_VAR(ncid, nvvelid, pt % vvel, (/jn/))
          nret = NF90_PUT_VAR(ncid, nmassid, pt % mass, (/jn/))
          nret = NF90_PUT_VAR(ncid, nthicknessid, pt % thickness, (/jn/))
          nret = NF90_PUT_VAR(ncid, nwidthid, pt % width, (/jn/))
          nret = NF90_PUT_VAR(ncid, nlengthid, pt % length, (/jn/))
          nret = NF90_PUT_VAR(ncid, nyearid, pt % year, (/jn/))
          nret = NF90_PUT_VAR(ncid, ndayid, pt % day, (/jn/))
          nret = NF90_PUT_VAR(ncid, nmass_of_bits_id, pt % mass_of_bits, (/jn/))
          nret = NF90_PUT_VAR(ncid, nheat_density_id, pt % heat_density, (/jn/))
          this => this % next
        END DO
      END IF
      nret = NF90_CLOSE(ncid)
      IF (nret /= NF90_NOERR) CALL ctl_stop('icebergs, write_restart: nf_close failed')
      jn = icb_utl_count()
      IF (lwp .AND. nn_verbose_level >= 0) WRITE(numout, FMT = '(2(a,i5))') 'icebergs, icb_rst_write: # bergs =', jn, ' on PE', &
&narea - 1
      IF (lk_mpp) THEN
        CALL mpp_sum('icbrst', jn)
      END IF
      IF (lwp) WRITE(numout, FMT = '(a,i5,a,i5,a)') 'icebergs, icb_rst_write: ', jn, ' bergs in total have been written at &
&timestep ', kt
    END IF
  END SUBROUTINE icb_rst_write
END MODULE icbrst
