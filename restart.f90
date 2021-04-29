MODULE restart
  USE oce
  USE dom_oce
  USE sbc_ice
  USE phycst
  USE eosbn2
  USE trdmxl_oce
  USE in_out_manager
  USE iom
  USE ioipsl, ONLY: ju2ymds
  USE diurnal_bulk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rst_opn
  PUBLIC :: rst_write
  PUBLIC :: rst_read
  PUBLIC :: rst_read_open
  CONTAINS
  SUBROUTINE rst_opn(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: iyear, imonth, iday
    REAL(KIND = wp) :: zsec
    REAL(KIND = wp) :: zfjulday
    CHARACTER(LEN = 20) :: clkt
    CHARACTER(LEN = 50) :: clname
    CHARACTER(LEN = lc) :: clpath
    CHARACTER(LEN = 52) :: clpname
    CHARACTER(LEN = 256) :: clinfo
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('rst_opn', 'r0', 0, 0)
    IF (kt == nit000) THEN
      lrst_oce = .FALSE.
      IF (ln_rst_list) THEN
        nrst_lst = 1
        nitrst = nstocklist(nrst_lst)
      ELSE
        nitrst = nitend
      END IF
    END IF
    IF (.NOT. ln_rst_list .AND. MOD(kt - 1, nstock) == 0) THEN
      nitrst = kt + nstock - 1
      IF (nitrst > nitend) nitrst = nitend
    END IF
    IF (kt == nitrst - 1 .OR. nstock == 1 .OR. (kt == nitend .AND. .NOT. lrst_oce)) THEN
      IF (nitrst <= nitend .AND. nitrst > 0) THEN
        IF (ln_rstdate) THEN
          zfjulday = fjulday + rdt / rday
          IF (ABS(zfjulday - REAL(NINT(zfjulday), wp)) < 0.1 / rday) zfjulday = REAL(NINT(zfjulday), wp)
          CALL ju2ymds(zfjulday, iyear, imonth, iday, zsec)
          WRITE(clkt, '(i4.4,2i2.2)') iyear, imonth, iday
        ELSE
          IF (nitrst > 999999999) THEN
            WRITE(clkt, *) nitrst
          ELSE
            WRITE(clkt, '(i8.8)') nitrst
          END IF
        END IF
        clname = TRIM(cexper) // "_" // TRIM(ADJUSTL(clkt)) // "_" // TRIM(cn_ocerst_out)
        clpath = TRIM(cn_ocerst_outdir)
        IF (clpath(LEN_TRIM(clpath) :) /= '/') clpath = TRIM(clpath) // '/'
        IF (lwp) THEN
          WRITE(numout, *)
          IF (.NOT. lwxios) THEN
            WRITE(numout, *) '             open ocean restart NetCDF file: ', TRIM(clpath) // TRIM(clname)
            IF (snc4set % luse) WRITE(numout, *) '             opened for NetCDF4 chunking and compression'
            IF (kt == nitrst - 1) THEN
              WRITE(numout, *) '             kt = nitrst - 1 = ', kt
            ELSE
              WRITE(numout, *) '             kt = ', kt
            END IF
          END IF
        END IF
        IF (.NOT. lwxios) THEN
          CALL iom_open(TRIM(clpath) // TRIM(clname), numrow, ldwrt = .TRUE.)
        ELSE
          clinfo = 'Can not use XIOS in rst_opn'
          CALL ctl_stop(TRIM(clinfo))
        END IF
        lrst_oce = .TRUE.
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE rst_opn
  SUBROUTINE rst_write(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('rst_write', 'r0', 0, 0)
    IF (lwxios) CALL iom_swap(cwxios_context)
    CALL iom_rstput(kt, nitrst, numrow, 'rdt', rdt, ldxios = lwxios)
    CALL iom_delay_rst('WRITE', 'OCE', numrow)
    IF (.NOT. ln_diurnal_only) THEN
      CALL iom_rstput(kt, nitrst, numrow, 'ub', ub, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'vb', vb, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'tb', tsb(:, :, :, jp_tem), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sb', tsb(:, :, :, jp_sal), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sshb', sshb, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'un', un, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'vn', vn, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'tn', tsn(:, :, :, jp_tem), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sn', tsn(:, :, :, jp_sal), ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sshn', sshn, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'rhop', rhop, ldxios = lwxios)
      IF (ln_iscpl) THEN
        CALL iom_rstput(kt, nitrst, numrow, 'tmask', tmask, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'umask', umask, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'vmask', vmask, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'smask', ssmask, ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'e3t_n', e3t_n(:, :, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'e3u_n', e3u_n(:, :, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'e3v_n', e3v_n(:, :, :), ldxios = lwxios)
        CALL iom_rstput(kt, nitrst, numrow, 'gdepw_n', gdepw_n(:, :, :), ldxios = lwxios)
      END IF
    END IF
    IF (ln_diurnal) CALL iom_rstput(kt, nitrst, numrow, 'Dsst', x_dsst, ldxios = lwxios)
    IF (lwxios) CALL iom_swap(cxios_context)
    IF (kt == nitrst) THEN
      IF (.NOT. lwxios) THEN
        CALL iom_close(numrow)
      ELSE
        CALL iom_context_finalize(cwxios_context)
      END IF
      lrst_oce = .FALSE.
      IF (ln_rst_list) THEN
        nrst_lst = MIN(nrst_lst + 1, SIZE(nstocklist, 1))
        nitrst = nstocklist(nrst_lst)
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE rst_write
  SUBROUTINE rst_read_open
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: llok
    CHARACTER(LEN = lc) :: clpath
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('rst_read_open', 'r0', 0, 0)
    IF (numror <= 0) THEN
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) 'rst_read : read oce NetCDF restart file'
        IF (snc4set % luse) WRITE(numout, *) 'rst_read : configured with NetCDF4 support'
        WRITE(numout, *) '~~~~~~~~'
      END IF
      lxios_sini = .FALSE.
      clpath = TRIM(cn_ocerst_indir)
      IF (clpath(LEN_TRIM(clpath) :) /= '/') clpath = TRIM(clpath) // '/'
      CALL iom_open(TRIM(clpath) // cn_ocerst_in, numror)
      IF (.NOT. lxios_set) lrxios = lrxios .AND. lxios_sini
      IF (lrxios) THEN
        crxios_context = 'nemo_rst'
        IF (.NOT. lxios_set) THEN
          IF (lwp) WRITE(numout, *) 'Enable restart reading by XIOS'
          CALL iom_init(crxios_context, ld_tmppatch = .FALSE.)
          lxios_set = .TRUE.
        END IF
      END IF
      IF (TRIM(Agrif_CFixed()) /= '0' .AND. lrxios) THEN
        CALL iom_init(crxios_context, ld_tmppatch = .FALSE.)
        IF (lwp) WRITE(numout, *) 'Enable restart reading by XIOS for AGRIF'
        lxios_set = .TRUE.
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE rst_read_open
  SUBROUTINE rst_read
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp) :: zrdt
    INTEGER :: jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: w3d
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('rst_read', 'r0', 0, 0)
    CALL rst_read_open
    IF (iom_varid(numror, 'rdt', ldstop = .FALSE.) > 0) THEN
      CALL iom_get(numror, 'rdt', zrdt, ldxios = lrxios)
      IF (zrdt /= rdt) neuler = 0
    END IF
    CALL iom_delay_rst('READ', 'OCE', numror)
    IF (ln_diurnal) CALL iom_get(numror, jpdom_autoglo, 'Dsst', x_dsst, ldxios = lrxios)
    CALL profile_psy_data0 % PostEnd
    IF (ln_diurnal_only) THEN
      CALL profile_psy_data1 % PreStart('rst_read', 'r1', 0, 0)
      IF (lwp) WRITE(numout, *) "rst_read:- ln_diurnal_only set, setting rhop=rau0"
      rhop = rau0
      CALL iom_get(numror, jpdom_autoglo, 'tn', w3d, ldxios = lrxios)
      CALL profile_psy_data1 % PostEnd
      !$ACC KERNELS
      tsn(:, :, 1, jp_tem) = w3d(:, :, 1)
      !$ACC END KERNELS
      RETURN
    END IF
    CALL profile_psy_data2 % PreStart('rst_read', 'r2', 0, 0)
    IF (iom_varid(numror, 'ub', ldstop = .FALSE.) > 0) THEN
      CALL iom_get(numror, jpdom_autoglo, 'ub', ub, ldxios = lrxios)
      CALL iom_get(numror, jpdom_autoglo, 'vb', vb, ldxios = lrxios)
      CALL iom_get(numror, jpdom_autoglo, 'tb', tsb(:, :, :, jp_tem), ldxios = lrxios)
      CALL iom_get(numror, jpdom_autoglo, 'sb', tsb(:, :, :, jp_sal), ldxios = lrxios)
      CALL iom_get(numror, jpdom_autoglo, 'sshb', sshb, ldxios = lrxios)
    ELSE
      neuler = 0
    END IF
    CALL iom_get(numror, jpdom_autoglo, 'un', un, ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'vn', vn, ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'tn', tsn(:, :, :, jp_tem), ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'sn', tsn(:, :, :, jp_sal), ldxios = lrxios)
    CALL iom_get(numror, jpdom_autoglo, 'sshn', sshn, ldxios = lrxios)
    IF (iom_varid(numror, 'rhop', ldstop = .FALSE.) > 0) THEN
      CALL iom_get(numror, jpdom_autoglo, 'rhop', rhop, ldxios = lrxios)
    ELSE
      CALL eos(tsn, rhd, rhop, gdept_n(:, :, :))
    END IF
    CALL profile_psy_data2 % PostEnd
    IF (neuler == 0) THEN
      !$ACC KERNELS
      tsb(:, :, :, :) = tsn(:, :, :, :)
      ub(:, :, :) = un(:, :, :)
      vb(:, :, :) = vn(:, :, :)
      sshb(:, :) = sshn(:, :)
      !$ACC END KERNELS
      IF (.NOT. ln_linssh) THEN
        !$ACC KERNELS
        DO jk = 1, jpk
          e3t_b(:, :, jk) = e3t_n(:, :, jk)
        END DO
        !$ACC END KERNELS
      END IF
    END IF
  END SUBROUTINE rst_read
END MODULE restart