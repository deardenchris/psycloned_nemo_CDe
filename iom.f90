MODULE iom
  USE dom_oce
  USE c1d
  USE flo_oce
  USE lbclnk
  USE iom_def
  USE iom_nf90
  USE in_out_manager
  USE lib_mpp
  USE domngb
  USE phycst
  USE ioipsl, ONLY: ju2ymds
  USE crs
  USE lib_fortran
  USE diurnal_bulk, ONLY: ln_diurnal_only, ln_diurnal
  IMPLICIT NONE
  PUBLIC
  LOGICAL, PUBLIC, PARAMETER :: lk_iomput = .FALSE.
  PUBLIC :: iom_init, iom_swap, iom_open, iom_close, iom_setkt, iom_varid, iom_get
  PUBLIC :: iom_chkatt, iom_getatt, iom_putatt, iom_getszuld, iom_rstput, iom_delay_rst, iom_put
  PUBLIC :: iom_use, iom_context_finalize
  PRIVATE :: iom_rp0d, iom_rp1d, iom_rp2d, iom_rp3d
  PRIVATE :: iom_g0d, iom_g1d, iom_g2d, iom_g3d, iom_get_123d
  PRIVATE :: iom_p1d, iom_p2d, iom_p3d
  PUBLIC :: iom_set_rstw_var_active, iom_set_rstw_core, iom_set_rst_vars
  INTERFACE iom_get
    MODULE PROCEDURE iom_g0d, iom_g1d, iom_g2d, iom_g3d
  END INTERFACE
  INTERFACE iom_getatt
    MODULE PROCEDURE iom_g0d_iatt, iom_g1d_iatt, iom_g0d_ratt, iom_g1d_ratt, iom_g0d_catt
  END INTERFACE
  INTERFACE iom_putatt
    MODULE PROCEDURE iom_p0d_iatt, iom_p1d_iatt, iom_p0d_ratt, iom_p1d_ratt, iom_p0d_catt
  END INTERFACE
  INTERFACE iom_rstput
    MODULE PROCEDURE iom_rp0d, iom_rp1d, iom_rp2d, iom_rp3d
  END INTERFACE
  INTERFACE iom_put
    MODULE PROCEDURE iom_p0d, iom_p1d, iom_p2d, iom_p3d
  END INTERFACE iom_put
  CONTAINS
  SUBROUTINE iom_init(cdname, fname, ld_tmppatch)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    CHARACTER(LEN = *), OPTIONAL, INTENT(IN) :: fname
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_tmppatch
  END SUBROUTINE iom_init
  SUBROUTINE iom_set_rstw_var_active(field)
    CHARACTER(LEN = *), INTENT(IN) :: field
    INTEGER :: i
    LOGICAL :: llis_set
    CHARACTER(LEN = 256) :: clinfo
    clinfo = 'iom_set_rstw_var_active: key_iomput is needed to use XIOS restart read/write functionality'
    CALL ctl_stop('STOP', TRIM(clinfo))
  END SUBROUTINE iom_set_rstw_var_active
  SUBROUTINE iom_set_rstr_active
    CHARACTER(LEN = 256) :: rst_file
  END SUBROUTINE iom_set_rstr_active
  SUBROUTINE iom_set_rstw_core(cdmdl)
    CHARACTER(LEN = *), INTENT(IN) :: cdmdl
    CHARACTER(LEN = 256) :: clinfo
    clinfo = 'iom_set_rstw_core: key_iomput is needed to use XIOS restart read/write functionality'
    CALL ctl_stop('STOP', TRIM(clinfo))
  END SUBROUTINE iom_set_rstw_core
  SUBROUTINE iom_set_rst_vars(fields)
    TYPE(RST_FIELD), INTENT(INOUT) :: fields(max_rst_fields)
    INTEGER :: i
    i = 0
    i = i + 1
    fields(i) % vname = "rdt"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "un"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "ub"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "vn"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "vb"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "tn"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "tb"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "sn"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "sb"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "sshn"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sshb"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "rhop"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "kt"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "ndastp"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "adatrj"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "utau_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "vtau_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "qns_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "emp_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sfx_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "en"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "avt_k"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "avm_k"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "dissl"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "sbc_hc_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sbc_sc_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "qsr_hc_b"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "fraqsr_1lev"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "greenland_icesheet_mass"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "greenland_icesheet_timelapsed"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "greenland_icesheet_mass_roc"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "antarctica_icesheet_mass"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "antarctica_icesheet_timelapsed"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "antarctica_icesheet_mass_roc"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "frc_v"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "frc_t"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "frc_s"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "frc_wn_t"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "frc_wn_s"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "ssh_ini"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "e3t_ini"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "hc_loc_ini"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "sc_loc_ini"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "ssh_hc_loc_ini"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ssh_sc_loc_ini"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "tilde_e3t_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "tilde_e3t_n"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "hdiv_lf"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ub2_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "vb2_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sshbb_e"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ubb_e"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "vbb_e"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sshb_e"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ub_e"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "vb_e"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "fwf_isf_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "isf_sc_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "isf_hc_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ssh_ibb"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "rnf_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "rnf_hc_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "rnf_sc_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "nn_fsbc"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "ssu_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ssv_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sst_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "sss_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ssh_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "e3t_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "frq_m"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "avmb"
    fields(i) % grid = "grid_vector"
    i = i + 1
    fields(i) % vname = "avtb"
    fields(i) % grid = "grid_vector"
    i = i + 1
    fields(i) % vname = "ub2_i_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "vb2_i_b"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "ntime"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "Dsst"
    fields(i) % grid = "grid_scalar"
    i = i + 1
    fields(i) % vname = "tmask"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "umask"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "vmask"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "smask"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "gdepw_n"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "e3t_n"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "e3u_n"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "e3v_n"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "surf_ini"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "e3t_b"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "hmxl_n"
    fields(i) % grid = "grid_N_3D"
    i = i + 1
    fields(i) % vname = "un_bf"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "vn_bf"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "hbl"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "hbli"
    fields(i) % grid = "grid_N"
    i = i + 1
    fields(i) % vname = "wn"
    fields(i) % grid = "grid_N_3D"
    IF (i - 1 > max_rst_fields) THEN
      WRITE(ctmp1, FMT = *) 'E R R O R : iom_set_rst_vars SIZE of RST_FIELD array is too small'
      CALL ctl_stop('iom_set_rst_vars:', ctmp1)
    END IF
  END SUBROUTINE iom_set_rst_vars
  SUBROUTINE iom_set_rstw_active(cdrst_file)
    CHARACTER(LEN = *) :: cdrst_file
  END SUBROUTINE iom_set_rstw_active
  SUBROUTINE iom_set_rst_context(ld_rstr)
    LOGICAL, INTENT(IN) :: ld_rstr
  END SUBROUTINE iom_set_rst_context
  SUBROUTINE iom_swap(cdname)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
  END SUBROUTINE iom_swap
  SUBROUTINE iom_open(cdname, kiomid, ldwrt, kdom, ldstop, ldiof, kdlev)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    INTEGER, INTENT(OUT) :: kiomid
    LOGICAL, INTENT(IN), OPTIONAL :: ldwrt
    INTEGER, INTENT(IN), OPTIONAL :: kdom
    LOGICAL, INTENT(IN), OPTIONAL :: ldstop
    LOGICAL, INTENT(IN), OPTIONAL :: ldiof
    INTEGER, INTENT(IN), OPTIONAL :: kdlev
    CHARACTER(LEN = 256) :: clname
    CHARACTER(LEN = 256) :: cltmpn
    CHARACTER(LEN = 10) :: clsuffix
    CHARACTER(LEN = 15) :: clcpu
    CHARACTER(LEN = 256) :: clinfo
    LOGICAL :: llok
    LOGICAL :: llwrt
    LOGICAL :: llnoov
    LOGICAL :: llstop
    LOGICAL :: lliof
    INTEGER :: icnt
    INTEGER :: iln, ils
    INTEGER :: idom
    INTEGER :: istop
    INTEGER, DIMENSION(2, 5) :: idompar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_open', 'r0', 0, 0)
    kiomid = - 1
    clinfo = '                    iom_open ~~~  '
    istop = nstop
    IF (Agrif_Root()) THEN
      IF (iom_open_init == 0) THEN
        iom_file(:) % nfid = 0
        iom_open_init = 1
      END IF
    END IF
    IF (PRESENT(ldwrt)) THEN
      llwrt = ldwrt
    ELSE
      llwrt = .FALSE.
    END IF
    IF (PRESENT(ldstop)) THEN
      llstop = ldstop
    ELSE
      llstop = .TRUE.
    END IF
    IF (PRESENT(ldiof)) THEN
      lliof = ldiof
    ELSE
      lliof = .FALSE.
    END IF
    llnoov = (jpni * jpnj) == jpnij .AND. .NOT. lk_agrif
    clname = TRIM(cdname)
    IF (.NOT. Agrif_Root() .AND. .NOT. lliof) THEN
      iln = INDEX(clname, '/')
      cltmpn = clname(1 : iln)
      clname = clname(iln + 1 : LEN_TRIM(clname))
      clname = TRIM(cltmpn) // TRIM(Agrif_CFixed()) // '_' // TRIM(clname)
    END IF
    clsuffix = '.nc'
    iln = LEN_TRIM(clname)
    ils = LEN_TRIM(clsuffix)
    IF (iln <= ils .OR. INDEX(TRIM(clname), TRIM(clsuffix), back = .TRUE.) /= iln - ils + 1) clname = TRIM(clname) // TRIM(clsuffix)
    cltmpn = clname
    INQUIRE(FILE = clname, EXIST = llok)
    IF (.NOT. llok) THEN
      WRITE(clcpu, FMT = *) narea - 1
      clcpu = TRIM(ADJUSTL(clcpu))
      iln = INDEX(clname, TRIM(clsuffix), back = .TRUE.)
      clname = clname(1 : iln - 1) // '_' // TRIM(clcpu) // TRIM(clsuffix)
      icnt = 0
      INQUIRE(FILE = clname, EXIST = llok)
      DO WHILE (.NOT. llok .AND. icnt < jpmax_digits)
        clcpu = "0" // TRIM(clcpu)
        clname = clname(1 : iln - 1) // '_' // TRIM(clcpu) // TRIM(clsuffix)
        INQUIRE(FILE = clname, EXIST = llok)
        icnt = icnt + 1
      END DO
    ELSE
      lxios_sini = .TRUE.
    END IF
    IF (llwrt) THEN
      IF (llnoov) THEN
        idom = jpdom_local_noovlap
      ELSE
        idom = jpdom_local_full
      END IF
      IF (PRESENT(kdom)) idom = kdom
      SELECT CASE (idom)
      CASE (jpdom_local_full)
        idompar(:, 1) = (/jpi, jpj/)
        idompar(:, 2) = (/nimpp, njmpp/)
        idompar(:, 3) = (/nimpp + jpi - 1, njmpp + jpj - 1/)
        idompar(:, 4) = (/nldi - 1, nldj - 1/)
        idompar(:, 5) = (/jpi - nlei, jpj - nlej/)
      CASE (jpdom_local_noextra)
        idompar(:, 1) = (/nlci, nlcj/)
        idompar(:, 2) = (/nimpp, njmpp/)
        idompar(:, 3) = (/nimpp + nlci - 1, njmpp + nlcj - 1/)
        idompar(:, 4) = (/nldi - 1, nldj - 1/)
        idompar(:, 5) = (/nlci - nlei, nlcj - nlej/)
      CASE (jpdom_local_noovlap)
        idompar(:, 1) = (/nlei - nldi + 1, nlej - nldj + 1/)
        idompar(:, 2) = (/nimpp + nldi - 1, njmpp + nldj - 1/)
        idompar(:, 3) = (/nimpp + nlei - 1, njmpp + nlej - 1/)
        idompar(:, 4) = (/0, 0/)
        idompar(:, 5) = (/0, 0/)
      CASE DEFAULT
        CALL ctl_stop(TRIM(clinfo), 'wrong value of kdom, only jpdom_local* cases are accepted')
      END SELECT
    END IF
    IF (MINVAL(iom_file(:) % nfid) /= 0) CALL ctl_stop(TRIM(clinfo), 'No more free file identifier', 'increase jpmax_files in &
&iom_def')
    IF (.NOT. llok) THEN
      IF (.NOT. llwrt) THEN
        IF (llstop) THEN
          CALL ctl_stop(TRIM(clinfo), 'File ' // TRIM(cltmpn) // '* not found')
        ELSE
          istop = nstop + 1
        END IF
      ELSE
        clname = cltmpn
      END IF
    ELSE
      IF (llwrt .AND. .NOT. ln_clobber) THEN
        CALL ctl_stop(TRIM(clinfo), 'We want to write in a new file but ' // TRIM(clname) // ' already exists...')
        istop = nstop + 1
      ELSE IF (llwrt) THEN
        clname = cltmpn
      END IF
    END IF
    IF (istop == nstop) THEN
      CALL iom_nf90_open(clname, kiomid, llwrt, llok, idompar, kdlev = kdlev)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_open
  SUBROUTINE iom_close(kiomid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT), OPTIONAL :: kiomid
    INTEGER :: jf
    INTEGER :: i_s, i_e
    CHARACTER(LEN = 100) :: clinfo
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_close', 'r0', 0, 0)
    clinfo = '                    iom_close ~~~  '
    IF (PRESENT(kiomid)) THEN
      i_s = kiomid
      i_e = kiomid
    ELSE
      i_s = 1
      i_e = jpmax_files
    END IF
    IF (i_s > 0) THEN
      DO jf = i_s, i_e
        IF (iom_file(jf) % nfid > 0) THEN
          CALL iom_nf90_close(jf)
          iom_file(jf) % nfid = 0
          IF (PRESENT(kiomid)) kiomid = 0
          IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' close file: ' // TRIM(iom_file(jf) % name) // ' ok'
        ELSE IF (PRESENT(kiomid)) THEN
          WRITE(ctmp1, FMT = *) '--->', kiomid
          CALL ctl_stop(TRIM(clinfo) // ' Invalid file identifier', ctmp1)
        END IF
      END DO
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_close
  FUNCTION iom_varid(kiomid, cdvar, kdimsz, kndims, ldstop)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: kdimsz
    INTEGER, INTENT(OUT), OPTIONAL :: kndims
    LOGICAL, INTENT(IN), OPTIONAL :: ldstop
    INTEGER :: iom_varid, iiv, i_nvd
    LOGICAL :: ll_fnd
    CHARACTER(LEN = 100) :: clinfo
    LOGICAL :: llstop
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_varid', 'r0', 0, 0)
    iom_varid = 0
    IF (PRESENT(ldstop)) THEN
      llstop = ldstop
    ELSE
      llstop = .TRUE.
    END IF
    IF (kiomid > 0) THEN
      clinfo = 'iom_varid, file: ' // TRIM(iom_file(kiomid) % name) // ', var: ' // TRIM(cdvar)
      IF (iom_file(kiomid) % nfid == 0) THEN
        CALL ctl_stop(TRIM(clinfo), 'the file is not open')
      ELSE
        ll_fnd = .FALSE.
        iiv = 0
        DO WHILE (.NOT. ll_fnd .AND. iiv < iom_file(kiomid) % nvars)
          iiv = iiv + 1
          ll_fnd = (TRIM(cdvar) == TRIM(iom_file(kiomid) % cn_var(iiv)))
        END DO
        IF (.NOT. ll_fnd) THEN
          iiv = iiv + 1
          IF (iiv <= jpmax_vars) THEN
            iom_varid = iom_nf90_varid(kiomid, cdvar, iiv, kdimsz, kndims)
          ELSE
            CALL ctl_stop(TRIM(clinfo), 'Too many variables in the file ' // iom_file(kiomid) % name, &
&'increase the parameter jpmax_vars')
          END IF
          IF (llstop .AND. iom_varid == - 1) CALL ctl_stop(TRIM(clinfo) // ' not found')
        ELSE
          iom_varid = iiv
          IF (PRESENT(kdimsz)) THEN
            i_nvd = iom_file(kiomid) % ndims(iiv)
            IF (i_nvd <= SIZE(kdimsz)) THEN
              kdimsz(1 : i_nvd) = iom_file(kiomid) % dimsz(1 : i_nvd, iiv)
            ELSE
              WRITE(ctmp1, FMT = *) i_nvd, SIZE(kdimsz)
              CALL ctl_stop(TRIM(clinfo), 'error in kdimsz size' // TRIM(ctmp1))
            END IF
          END IF
          IF (PRESENT(kndims)) kndims = iom_file(kiomid) % ndims(iiv)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION iom_varid
  SUBROUTINE iom_g0d(kiomid, cdvar, pvar, ktime, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(OUT) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktime
    LOGICAL, INTENT(IN), OPTIONAL :: ldxios
    INTEGER :: idvar
    INTEGER :: idmspc
    INTEGER, DIMENSION(1) :: itime
    CHARACTER(LEN = 100) :: clinfo
    CHARACTER(LEN = 100) :: clname
    CHARACTER(LEN = 1) :: cldmspc
    LOGICAL :: llxios
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g0d', 'r0', 0, 0)
    llxios = .FALSE.
    IF (PRESENT(ldxios)) llxios = ldxios
    IF (.NOT. llxios) THEN
      itime = 1
      IF (PRESENT(ktime)) itime = ktime
      clname = iom_file(kiomid) % name
      clinfo = '          iom_g0d, file: ' // TRIM(clname) // ', var: ' // TRIM(cdvar)
      IF (kiomid > 0) THEN
        idvar = iom_varid(kiomid, cdvar)
        IF (iom_file(kiomid) % nfid > 0 .AND. idvar > 0) THEN
          idmspc = iom_file(kiomid) % ndims(idvar)
          IF (iom_file(kiomid) % luld(idvar)) idmspc = idmspc - 1
          WRITE(cldmspc, FMT = '(i1)') idmspc
          IF (idmspc > 0) CALL ctl_stop(TRIM(clinfo), 'When reading to a 0D array, we do not accept data', 'with 1 or more spatial &
&dimensions: ' // cldmspc // ' were found.', 'Use ncwa -a to suppress the unnecessary dimensions')
          CALL iom_nf90_get(kiomid, idvar, pvar, itime)
        END IF
      END IF
    ELSE
      WRITE(ctmp1, FMT = *) 'Can not use XIOS in iom_g0d, file: ' // TRIM(clname) // ', var:' // TRIM(cdvar)
      CALL ctl_stop('iom_g0d', ctmp1)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g0d
  SUBROUTINE iom_g1d(kiomid, kdom, cdvar, pvar, ktime, kstart, kcount, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    INTEGER, INTENT(IN) :: kdom
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktime
    INTEGER, INTENT(IN), DIMENSION(1), OPTIONAL :: kstart
    INTEGER, INTENT(IN), DIMENSION(1), OPTIONAL :: kcount
    LOGICAL, INTENT(IN), OPTIONAL :: ldxios
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g1d', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_get_123d(kiomid, kdom, cdvar, pv_r1d = pvar, ktime = ktime, kstart = kstart, &
&kcount = kcount, ldxios = ldxios)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g1d
  SUBROUTINE iom_g2d(kiomid, kdom, cdvar, pvar, ktime, kstart, kcount, lrowattr, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    INTEGER, INTENT(IN) :: kdom
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:, :) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktime
    INTEGER, INTENT(IN), DIMENSION(2), OPTIONAL :: kstart
    INTEGER, INTENT(IN), DIMENSION(2), OPTIONAL :: kcount
    LOGICAL, INTENT(IN), OPTIONAL :: lrowattr
    LOGICAL, INTENT(IN), OPTIONAL :: ldxios
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g2d', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_get_123d(kiomid, kdom, cdvar, pv_r2d = pvar, ktime = ktime, kstart = kstart, &
&kcount = kcount, lrowattr = lrowattr, ldxios = ldxios)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g2d
  SUBROUTINE iom_g3d(kiomid, kdom, cdvar, pvar, ktime, kstart, kcount, lrowattr, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    INTEGER, INTENT(IN) :: kdom
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:, :, :) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktime
    INTEGER, INTENT(IN), DIMENSION(3), OPTIONAL :: kstart
    INTEGER, INTENT(IN), DIMENSION(3), OPTIONAL :: kcount
    LOGICAL, INTENT(IN), OPTIONAL :: lrowattr
    LOGICAL, INTENT(IN), OPTIONAL :: ldxios
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g3d', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_get_123d(kiomid, kdom, cdvar, pv_r3d = pvar, ktime = ktime, kstart = kstart, &
&kcount = kcount, lrowattr = lrowattr, ldxios = ldxios)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g3d
  SUBROUTINE iom_get_123d(kiomid, kdom, cdvar, pv_r1d, pv_r2d, pv_r3d, ktime, kstart, kcount, lrowattr, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    INTEGER, INTENT(IN) :: kdom
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT), OPTIONAL :: pv_r1d
    REAL(KIND = wp), DIMENSION(:, :), INTENT(OUT), OPTIONAL :: pv_r2d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT), OPTIONAL :: pv_r3d
    INTEGER, INTENT(IN), OPTIONAL :: ktime
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: kstart
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: kcount
    LOGICAL, INTENT(IN), OPTIONAL :: lrowattr
    LOGICAL, INTENT(IN), OPTIONAL :: ldxios
    LOGICAL :: llxios
    LOGICAL :: llnoov
    LOGICAL :: luse_jattr
    INTEGER :: jstartrow
    INTEGER :: jl
    INTEGER :: idom
    INTEGER :: idvar
    INTEGER :: inbdim
    INTEGER :: idmspc
    INTEGER :: itime
    INTEGER :: istop
    INTEGER :: ix1, ix2, iy1, iy2
    INTEGER :: ji, jj
    INTEGER :: irankpv
    INTEGER :: ind1, ind2
    INTEGER, DIMENSION(jpmax_dims) :: istart
    INTEGER, DIMENSION(jpmax_dims) :: icnt
    INTEGER, DIMENSION(jpmax_dims) :: idimsz
    INTEGER, DIMENSION(jpmax_dims) :: ishape
    REAL(KIND = wp) :: zscf, zofs
    INTEGER :: itmp
    CHARACTER(LEN = 256) :: clinfo
    CHARACTER(LEN = 256) :: clname
    CHARACTER(LEN = 1) :: clrankpv, cldmspc
    LOGICAL :: ll_depth_spec
    INTEGER :: inlev
    REAL(KIND = wp) :: gma, gmi
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    CALL profile_psy_data0 % PreStart('iom_get_123d', 'r0', 0, 0)
    inlev = - 1
    IF (PRESENT(pv_r3d)) inlev = SIZE(pv_r3d, 3)
    llxios = .FALSE.
    IF (PRESENT(ldxios)) llxios = ldxios
    idvar = iom_varid(kiomid, cdvar)
    idom = kdom
    CALL profile_psy_data0 % PostEnd
    IF (.NOT. llxios) THEN
      CALL profile_psy_data1 % PreStart('iom_get_123d', 'r1', 0, 0)
      clname = iom_file(kiomid) % name
      clinfo = '          iom_get_123d, file: ' // TRIM(clname) // ', var: ' // TRIM(cdvar)
      llnoov = (jpni * jpnj) == jpnij .AND. .NOT. lk_agrif
      IF (PRESENT(kcount) .AND. (.NOT. PRESENT(kstart))) CALL ctl_stop(TRIM(clinfo), 'kcount present needs kstart present')
      IF (PRESENT(kstart) .AND. (.NOT. PRESENT(kcount))) CALL ctl_stop(TRIM(clinfo), 'kstart present needs kcount present')
      IF (PRESENT(kstart) .AND. idom /= jpdom_unknown .AND. idom /= jpdom_autoglo_xy) CALL ctl_stop(TRIM(clinfo), 'kstart present &
&needs kdom = jpdom_unknown or kdom = jpdom_autoglo_xy')
      luse_jattr = .FALSE.
      IF (PRESENT(lrowattr)) THEN
        IF (lrowattr .AND. idom /= jpdom_data) CALL ctl_stop(TRIM(clinfo), 'lrowattr present and true needs kdom = jpdom_data')
        IF (lrowattr .AND. idom == jpdom_data) luse_jattr = .TRUE.
      END IF
      istop = nstop
      IF (idvar > 0) THEN
        idimsz(:) = iom_file(kiomid) % dimsz(:, idvar)
        inbdim = iom_file(kiomid) % ndims(idvar)
        idmspc = inbdim
        IF (iom_file(kiomid) % luld(idvar)) idmspc = inbdim - 1
        IF (idmspc > 3) CALL ctl_stop(TRIM(clinfo), 'the file has more than 3 spatial dimensions this case is not coded...')
        IF (idom == jpdom_autoglo_xy) THEN
          ll_depth_spec = .TRUE.
          idom = jpdom_autoglo
        ELSE
          ll_depth_spec = .FALSE.
        END IF
        IF (idom == jpdom_autoglo .OR. idom == jpdom_autodta) THEN
          IF (idom == jpdom_autoglo) THEN
            idom = jpdom_global
          ELSE
            idom = jpdom_data
          END IF
          ind1 = INDEX(clname, '_', back = .TRUE.) + 1
          ind2 = INDEX(clname, '.', back = .TRUE.) - 1
          IF (ind2 > ind1) THEN
            IF (VERIFY(clname(ind1 : ind2), '0123456789') == 0) idom = jpdom_local
          END IF
        END IF
        IF (idom == jpdom_local) THEN
          IF (idimsz(1) == jpi .AND. idimsz(2) == jpj) THEN
            idom = jpdom_local_full
          ELSE IF (idimsz(1) == nlci .AND. idimsz(2) == nlcj) THEN
            idom = jpdom_local_noextra
          ELSE IF (idimsz(1) == (nlei - nldi + 1) .AND. idimsz(2) == (nlej - nldj + 1)) THEN
            idom = jpdom_local_noovlap
          ELSE
            CALL ctl_stop(TRIM(clinfo), 'impossible to identify the local domain')
          END IF
        END IF
        itime = 1
        IF (PRESENT(ktime)) itime = ktime
        irankpv = 1 * COUNT((/PRESENT(pv_r1d)/)) + 2 * COUNT((/PRESENT(pv_r2d)/)) + 3 * COUNT((/PRESENT(pv_r3d)/))
        WRITE(clrankpv, FMT = '(i1)') irankpv
        WRITE(cldmspc, FMT = '(i1)') idmspc
        IF (idmspc < irankpv) THEN
          CALL ctl_stop(TRIM(clinfo), 'The file has only ' // cldmspc // ' spatial dimension', &
&'it is impossible to read a ' // clrankpv // 'D array from this file...')
        ELSE IF (idmspc == irankpv) THEN
          IF (PRESENT(pv_r1d) .AND. idom /= jpdom_unknown) CALL ctl_stop(TRIM(clinfo), 'case not coded...You must use &
&jpdom_unknown')
        ELSE IF (idmspc > irankpv) THEN
          IF (PRESENT(pv_r2d) .AND. itime == 1 .AND. idimsz(3) == 1 .AND. idmspc == 3) THEN
            CALL ctl_warn(TRIM(clinfo), '2D array but 3 spatial dimensions for the data...', &
&'As the size of the z dimension is 1 and as we try to read the first record, ', 'we accept this case, &
&even if there is a possible mix-up between z and time dimension')
            idmspc = idmspc - 1
          ELSE
            CALL ctl_stop(TRIM(clinfo), 'To keep iom lisibility, when reading a ' // clrankpv // 'D array,', &
&'we do not accept data with ' // cldmspc // ' spatial dimensions', 'Use ncwa -a to suppress the unnecessary dimensions')
          END IF
        END IF
        icnt(:) = 1
        istart(:) = 1
        istart(idmspc + 1) = itime
        IF (PRESENT(kstart) .AND. .NOT. ll_depth_spec) THEN
          istart(1 : idmspc) = kstart(1 : idmspc)
          icnt(1 : idmspc) = kcount(1 : idmspc)
        ELSE
          IF (idom == jpdom_unknown) THEN
            icnt(1 : idmspc) = idimsz(1 : idmspc)
          ELSE
            IF (.NOT. PRESENT(pv_r1d)) THEN
              IF (idom == jpdom_data) THEN
                jstartrow = 1
                IF (luse_jattr) THEN
                  CALL iom_getatt(kiomid, 'open_ocean_jstart', jstartrow)
                  jstartrow = MAX(1, jstartrow)
                END IF
                istart(1 : 2) = (/mig(1), mjg(1) + jstartrow - 1/)
              ELSE IF (idom == jpdom_global) THEN
                istart(1 : 2) = (/nimpp, njmpp/)
              END IF
              IF (llnoov .AND. idom /= jpdom_local_noovlap) istart(1 : 2) = istart(1 : 2) + (/nldi - 1, nldj - 1/)
              IF (llnoov) THEN
                icnt(1 : 2) = (/nlei - nldi + 1, nlej - nldj + 1/)
              ELSE
                icnt(1 : 2) = (/nlci, nlcj/)
              END IF
              IF (PRESENT(pv_r3d)) THEN
                IF (idom == jpdom_data) THEN
                  icnt(3) = inlev
                ELSE IF (ll_depth_spec .AND. PRESENT(kstart)) THEN
                  istart(3) = kstart(3)
                  icnt(3) = kcount(3)
                ELSE
                  icnt(3) = inlev
                END IF
              END IF
            END IF
          END IF
        END IF
        DO jl = 1, jpmax_dims
          itmp = istart(jl) + icnt(jl) - 1
          IF (itmp > idimsz(jl) .AND. idimsz(jl) /= 0) THEN
            WRITE(ctmp1, FMT = "('(istart(', i1, ') + icnt(', i1, ') - 1) = ', i5)") jl, jl, itmp
            WRITE(ctmp2, FMT = "(' is larger than idimsz(', i1,') = ', i5)") jl, idimsz(jl)
            CALL ctl_stop(TRIM(clinfo), 'start and count too big regarding to the size of the data, ', ctmp1, ctmp2)
          END IF
        END DO
        IF (idom == jpdom_unknown) THEN
          IF (irankpv == 1) ishape(1 : 1) = SHAPE(pv_r1d)
          IF (irankpv == 2) ishape(1 : 2) = SHAPE(pv_r2d)
          IF (irankpv == 3) ishape(1 : 3) = SHAPE(pv_r3d)
          ctmp1 = 'd'
        ELSE
          IF (irankpv == 2) THEN
            IF (llnoov) THEN
              ishape(1 : 2) = SHAPE(pv_r2d(nldi : nlei, nldj : nlej))
              ctmp1 = 'd(nldi:nlei,nldj:nlej)'
            ELSE
              ishape(1 : 2) = SHAPE(pv_r2d(1 : nlci, 1 : nlcj))
              ctmp1 = 'd(1:nlci,1:nlcj)'
            END IF
          END IF
          IF (irankpv == 3) THEN
            IF (llnoov) THEN
              ishape(1 : 3) = SHAPE(pv_r3d(nldi : nlei, nldj : nlej, :))
              ctmp1 = 'd(nldi:nlei,nldj:nlej,:)'
            ELSE
              ishape(1 : 3) = SHAPE(pv_r3d(1 : nlci, 1 : nlcj, :))
              ctmp1 = 'd(1:nlci,1:nlcj,:)'
            END IF
          END IF
        END IF
        DO jl = 1, irankpv
          WRITE(ctmp2, FMT = "(', ', i1,'): ', i5,' /= icnt(', i1,'):', i5)") jl, ishape(jl), jl, icnt(jl)
          IF (ishape(jl) /= icnt(jl)) CALL ctl_stop(TRIM(clinfo), 'size(pv_r' // clrankpv // TRIM(ctmp1) // TRIM(ctmp2))
        END DO
      END IF
      CALL profile_psy_data1 % PostEnd
      IF (idvar > 0 .AND. istop == nstop) THEN
        CALL profile_psy_data2 % PreStart('iom_get_123d', 'r2', 0, 0)
        IF (llnoov) THEN
          IF (idom /= jpdom_unknown) THEN
            ix1 = nldi
            ix2 = nlei
            iy1 = nldj
            iy2 = nlej
          ELSE
            ix1 = 1
            ix2 = icnt(1)
            iy1 = 1
            iy2 = icnt(2)
          END IF
        ELSE
          IF (idom /= jpdom_unknown) THEN
            ix1 = 1
            ix2 = nlci
            iy1 = 1
            iy2 = nlcj
          ELSE
            ix1 = 1
            ix2 = icnt(1)
            iy1 = 1
            iy2 = icnt(2)
          END IF
        END IF
        CALL iom_nf90_get(kiomid, idvar, inbdim, istart, icnt, ix1, ix2, iy1, iy2, pv_r1d, pv_r2d, pv_r3d)
        CALL profile_psy_data2 % PostEnd
        IF (istop == nstop) THEN
          CALL profile_psy_data3 % PreStart('iom_get_123d', 'r3', 0, 0)
          IF (lwp) WRITE(numout, FMT = "(10x,' read ',a,' (rec: ',i6,') in ',a,' ok')") TRIM(cdvar), itime, TRIM(iom_file(kiomid) &
&% name)
          IF (PRESENT(pv_r2d) .AND. idom /= jpdom_unknown) THEN
            CALL lbc_lnk('iom', pv_r2d, 'Z', - 999., 'no0')
          ELSE IF (PRESENT(pv_r3d) .AND. idom /= jpdom_unknown) THEN
            IF (icnt(3) == inlev) THEN
              CALL lbc_lnk('iom', pv_r3d, 'Z', - 999., 'no0')
            ELSE
              DO jj = nlcj + 1, jpj
                pv_r3d(1 : nlci, jj, :) = pv_r3d(1 : nlci, nlej, :)
              END DO
              DO ji = nlci + 1, jpi
                pv_r3d(ji, :, :) = pv_r3d(nlei, :, :)
              END DO
            END IF
          END IF
          CALL profile_psy_data3 % PostEnd
        ELSE
          RETURN
        END IF
      ELSE
        RETURN
      END IF
    ELSE
      CALL profile_psy_data4 % PreStart('iom_get_123d', 'r4', 0, 0)
      istop = istop + 1
      clinfo = 'Can not use XIOS in iom_get_123d, file: ' // TRIM(clname) // ', var:' // TRIM(cdvar)
      CALL profile_psy_data4 % PostEnd
    END IF
    CALL profile_psy_data5 % PreStart('iom_get_123d', 'r5', 0, 0)
    IF (lk_c1d .AND. PRESENT(pv_r2d)) CALL lbc_lnk('iom', pv_r2d, 'Z', 1.)
    IF (lk_c1d .AND. PRESENT(pv_r3d)) CALL lbc_lnk('iom', pv_r3d, 'Z', 1.)
    zscf = iom_file(kiomid) % scf(idvar)
    zofs = iom_file(kiomid) % ofs(idvar)
    IF (PRESENT(pv_r1d)) THEN
      IF (zscf /= 1.) pv_r1d(:) = pv_r1d(:) * zscf
      IF (zofs /= 0.) pv_r1d(:) = pv_r1d(:) + zofs
    ELSE IF (PRESENT(pv_r2d)) THEN
      IF (zscf /= 1.) pv_r2d(:, :) = pv_r2d(:, :) * zscf
      IF (zofs /= 0.) pv_r2d(:, :) = pv_r2d(:, :) + zofs
    ELSE IF (PRESENT(pv_r3d)) THEN
      IF (zscf /= 1.) pv_r3d(:, :, :) = pv_r3d(:, :, :) * zscf
      IF (zofs /= 0.) pv_r3d(:, :, :) = pv_r3d(:, :, :) + zofs
    END IF
    CALL profile_psy_data5 % PostEnd
  END SUBROUTINE iom_get_123d
  FUNCTION iom_getszuld(kiomid)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    INTEGER :: iom_getszuld
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_getszuld', 'r0', 0, 0)
    iom_getszuld = - 1
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % iduld > 0) iom_getszuld = iom_file(kiomid) % lenuld
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION iom_getszuld
  SUBROUTINE iom_chkatt(kiomid, cdatt, llok, ksize, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    LOGICAL, INTENT(OUT) :: llok
    INTEGER, INTENT(OUT), OPTIONAL :: ksize
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_chkatt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_chkatt(kiomid, cdatt, llok, ksize = ksize, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_chkatt
  SUBROUTINE iom_g0d_iatt(kiomid, cdatt, katt0d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    INTEGER, INTENT(OUT) :: katt0d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g0d_iatt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_getatt(kiomid, cdatt, katt0d = katt0d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g0d_iatt
  SUBROUTINE iom_g1d_iatt(kiomid, cdatt, katt1d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    INTEGER, DIMENSION(:), INTENT(OUT) :: katt1d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g1d_iatt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_getatt(kiomid, cdatt, katt1d = katt1d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g1d_iatt
  SUBROUTINE iom_g0d_ratt(kiomid, cdatt, patt0d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    REAL(KIND = wp), INTENT(OUT) :: patt0d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g0d_ratt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_getatt(kiomid, cdatt, patt0d = patt0d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g0d_ratt
  SUBROUTINE iom_g1d_ratt(kiomid, cdatt, patt1d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    REAL(KIND = wp), DIMENSION(:), INTENT(OUT) :: patt1d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g1d_ratt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_getatt(kiomid, cdatt, patt1d = patt1d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g1d_ratt
  SUBROUTINE iom_g0d_catt(kiomid, cdatt, cdatt0d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    CHARACTER(LEN = *), INTENT(OUT) :: cdatt0d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_g0d_catt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_getatt(kiomid, cdatt, cdatt0d = cdatt0d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_g0d_catt
  SUBROUTINE iom_p0d_iatt(kiomid, cdatt, katt0d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    INTEGER, INTENT(IN) :: katt0d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_p0d_iatt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_putatt(kiomid, cdatt, katt0d = katt0d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_p0d_iatt
  SUBROUTINE iom_p1d_iatt(kiomid, cdatt, katt1d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    INTEGER, DIMENSION(:), INTENT(IN) :: katt1d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_p1d_iatt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_putatt(kiomid, cdatt, katt1d = katt1d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_p1d_iatt
  SUBROUTINE iom_p0d_ratt(kiomid, cdatt, patt0d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    REAL(KIND = wp), INTENT(IN) :: patt0d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_p0d_ratt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_putatt(kiomid, cdatt, patt0d = patt0d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_p0d_ratt
  SUBROUTINE iom_p1d_ratt(kiomid, cdatt, patt1d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: patt1d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_p1d_ratt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_putatt(kiomid, cdatt, patt1d = patt1d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_p1d_ratt
  SUBROUTINE iom_p0d_catt(kiomid, cdatt, cdatt0d, cdvar)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdatt
    CHARACTER(LEN = *), INTENT(IN) :: cdatt0d
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cdvar
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_p0d_catt', 'r0', 0, 0)
    IF (kiomid > 0) THEN
      IF (iom_file(kiomid) % nfid > 0) CALL iom_nf90_putatt(kiomid, cdatt, cdatt0d = cdatt0d, cdvar = cdvar)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_p0d_catt
  SUBROUTINE iom_rp0d(kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kwrite
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(IN) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktype
    LOGICAL, OPTIONAL :: ldxios
    LOGICAL :: llx
    INTEGER :: ivid
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    RETURN
    CALL profile_psy_data0 % PreStart('iom_rp0d', 'r0', 0, 0)
    llx = .FALSE.
    IF (PRESENT(ldxios)) llx = ldxios
    IF (llx) THEN
    ELSE
      IF (kiomid > 0) THEN
        IF (iom_file(kiomid) % nfid > 0) THEN
          ivid = iom_varid(kiomid, cdvar, ldstop = .FALSE.)
          CALL iom_nf90_rstput(kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r0d = pvar)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_rp0d
  SUBROUTINE iom_rp1d(kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kwrite
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(IN), DIMENSION(:) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktype
    LOGICAL, OPTIONAL :: ldxios
    LOGICAL :: llx
    INTEGER :: ivid
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    RETURN
    CALL profile_psy_data0 % PreStart('iom_rp1d', 'r0', 0, 0)
    llx = .FALSE.
    IF (PRESENT(ldxios)) llx = ldxios
    IF (llx) THEN
    ELSE
      IF (kiomid > 0) THEN
        IF (iom_file(kiomid) % nfid > 0) THEN
          ivid = iom_varid(kiomid, cdvar, ldstop = .FALSE.)
          CALL iom_nf90_rstput(kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r1d = pvar)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_rp1d
  SUBROUTINE iom_rp2d(kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kwrite
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktype
    LOGICAL, OPTIONAL :: ldxios
    LOGICAL :: llx
    INTEGER :: ivid
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    RETURN
    CALL profile_psy_data0 % PreStart('iom_rp2d', 'r0', 0, 0)
    llx = .FALSE.
    IF (PRESENT(ldxios)) llx = ldxios
    IF (llx) THEN
    ELSE
      IF (kiomid > 0) THEN
        IF (iom_file(kiomid) % nfid > 0) THEN
          ivid = iom_varid(kiomid, cdvar, ldstop = .FALSE.)
          CALL iom_nf90_rstput(kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r2d = pvar)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_rp2d
  SUBROUTINE iom_rp3d(kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kwrite
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :) :: pvar
    INTEGER, INTENT(IN), OPTIONAL :: ktype
    LOGICAL, OPTIONAL :: ldxios
    LOGICAL :: llx
    INTEGER :: ivid
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    RETURN
    CALL profile_psy_data0 % PreStart('iom_rp3d', 'r0', 0, 0)
    llx = .FALSE.
    IF (PRESENT(ldxios)) llx = ldxios
    IF (llx) THEN
    ELSE
      IF (kiomid > 0) THEN
        IF (iom_file(kiomid) % nfid > 0) THEN
          ivid = iom_varid(kiomid, cdvar, ldstop = .FALSE.)
          CALL iom_nf90_rstput(kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r3d = pvar)
        END IF
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE iom_rp3d
  SUBROUTINE iom_delay_rst(cdaction, cdcpnt, kncid)
    CHARACTER(LEN = *), INTENT(IN) :: cdaction
    CHARACTER(LEN = *), INTENT(IN) :: cdcpnt
    INTEGER, INTENT(IN) :: kncid
    INTEGER :: ji
    INTEGER :: indim
    LOGICAL :: llattexist
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: zreal1d
    IF (TRIM(cdaction) == 'READ') THEN
      DO ji = 1, nbdelay
        IF (c_delaycpnt(ji) == cdcpnt) THEN
          CALL iom_chkatt(kncid, 'DELAY_' // c_delaylist(ji), llattexist, indim)
          IF (llattexist) THEN
            ALLOCATE(todelay(ji) % z1d(indim))
            CALL iom_getatt(kncid, 'DELAY_' // c_delaylist(ji), todelay(ji) % z1d(:))
            ndelayid(ji) = 0
          END IF
        END IF
      END DO
    ELSE
      DO ji = 1, nbdelay
        IF (c_delaycpnt(ji) == cdcpnt) THEN
          IF (ASSOCIATED(todelay(ji) % z1d)) THEN
            CALL mpp_delay_rcv(ji)
            CALL iom_putatt(kncid, 'DELAY_' // c_delaylist(ji), todelay(ji) % z1d(:))
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE iom_delay_rst
  SUBROUTINE iom_p0d(cdname, pfield0d)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), INTENT(IN) :: pfield0d
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zz
    IF (.FALSE.) WRITE(numout, FMT = *) cdname, pfield0d
  END SUBROUTINE iom_p0d
  SUBROUTINE iom_p1d(cdname, pfield1d)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:), INTENT(IN) :: pfield1d
    IF (.FALSE.) WRITE(numout, FMT = *) cdname, pfield1d
  END SUBROUTINE iom_p1d
  SUBROUTINE iom_p2d(cdname, pfield2d)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: pfield2d
    IF (.FALSE.) WRITE(numout, FMT = *) cdname, pfield2d
  END SUBROUTINE iom_p2d
  SUBROUTINE iom_p3d(cdname, pfield3d)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pfield3d
    IF (.FALSE.) WRITE(numout, FMT = *) cdname, pfield3d
  END SUBROUTINE iom_p3d
  SUBROUTINE iom_setkt(kt, cdname)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    IF (.FALSE.) WRITE(numout, FMT = *) kt, cdname
  END SUBROUTINE iom_setkt
  SUBROUTINE iom_context_finalize(cdname)
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    IF (.FALSE.) WRITE(numout, FMT = *) cdname
  END SUBROUTINE iom_context_finalize
  LOGICAL FUNCTION iom_use(cdname)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN) :: cdname
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('iom_use', 'r0', 0, 0)
    iom_use = .FALSE.
    CALL profile_psy_data0 % PostEnd
  END FUNCTION iom_use
END MODULE iom