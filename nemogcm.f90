MODULE nemogcm
  USE step_oce
  USE phycst
  USE domain
  USE closea
  USE usrdef_nam
  USE tideini
  USE bdy_oce, ONLY: ln_bdy
  USE bdyini
  USE istate
  USE ldfdyn
  USE ldftra
  USE trdini
  USE asminc
  USE asmbkg
  USE diaptr
  USE diadct
  USE diaobs
  USE diacfl
  USE step
  USE icbini
  USE icbstp
  USE cpl_oasis3
  USE c1d
  USE step_c1d
  USE dyndmp
  USE stopar
  USE stopts
  USE diurnal_bulk
  USE step_diu
  USE crsini
  USE diatmb
  USE dia25h
  USE sbc_oce, ONLY: lk_oasis
  USE wet_dry
  USE lib_mpp
  USE mppini
  USE lbcnfd, ONLY: isendto, nsndto, nfsloop, nfeloop
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nemo_gcm
  PUBLIC :: nemo_init
  PUBLIC :: nemo_alloc
  CHARACTER(LEN = lc) :: cform_aaa = "( /, 'AAAAAAAA', / ) "
  CONTAINS
  SUBROUTINE nemo_gcm
    INTEGER :: istp
    CALL nemo_init
    CALL mpp_max('nemogcm', nstop)
    IF (lwp) WRITE(numout, cform_aaa)
    istp = nit000
    IF (.NOT. ln_diurnal_only) THEN
      DO WHILE (istp <= nitend .AND. nstop == 0)
        CALL stp(istp)
        istp = istp + 1
      END DO
    ELSE
      DO WHILE (istp <= nitend .AND. nstop == 0)
        CALL stp_diurnal(istp)
        istp = istp + 1
      END DO
    END IF
    IF (ln_diaobs) CALL dia_obs_wri
    IF (ln_icebergs) CALL icb_end(nitend)
    IF (lwp) WRITE(numout, cform_aaa)
    IF (nstop /= 0 .AND. lwp) THEN
      WRITE(numout, cform_err)
      WRITE(numout, *) '   ==>>>   nemo_gcm: a total of ', nstop, ' errors have been found'
      WRITE(numout, *)
    END IF
    IF (ln_timing) CALL timing_finalize
    CALL nemo_closefile
    IF (lk_oasis) THEN
      CALL cpl_finalize
    ELSE IF (lk_mpp) THEN
      CALL mppstop(ldfinal = .TRUE.)
    END IF
    IF (lwm) THEN
      IF (nstop == 0) THEN
        STOP 0
      ELSE
        STOP 999
      END IF
    END IF
  END SUBROUTINE nemo_gcm
  SUBROUTINE nemo_init
    INTEGER :: ji
    INTEGER :: ios, ilocal_comm
    CHARACTER(LEN = 120), DIMENSION(60) :: cltxt, cltxt2, clnam
    NAMELIST /namctl/ ln_ctl, sn_cfctl, nn_print, nn_ictls, nn_ictle, nn_isplt, nn_jsplt, nn_jctls, nn_jctle, ln_timing, ln_diacfl
    NAMELIST /namcfg/ ln_read_cfg, cn_domcfg, ln_closea, ln_write_cfg, cn_domcfg_out, ln_use_jattr
    cltxt = ''
    cltxt2 = ''
    clnam = ''
    cxios_context = 'nemo'
    CALL ctl_opn(numnam_ref, 'namelist_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE.)
    CALL ctl_opn(numnam_cfg, 'namelist_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE.)
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namctl, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namctl in reference namelist', .TRUE.)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namctl, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namctl in configuration namelist', .TRUE.)
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namcfg, IOSTAT = ios, ERR = 903)
903 IF (ios /= 0) CALL ctl_nam(ios, 'namcfg in reference namelist', .TRUE.)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namcfg, IOSTAT = ios, ERR = 904)
904 IF (ios > 0) CALL ctl_nam(ios, 'namcfg in configuration namelist', .TRUE.)
    IF (ln_read_cfg) THEN
      CALL domain_cfg(cltxt2, cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio)
    ELSE
      CALL usr_def_nam(cltxt2, clnam, cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio)
    END IF
    IF (lk_oasis) THEN
      IF (Agrif_Root()) THEN
        CALL cpl_init("oceanx", ilocal_comm)
      END IF
      narea = mynode(cltxt, 'output.namelist.dyn', numnam_ref, numnam_cfg, numond, nstop, ilocal_comm)
    ELSE
      ilocal_comm = 0
      narea = mynode(cltxt, 'output.namelist.dyn', numnam_ref, numnam_cfg, numond, nstop)
    END IF
    narea = narea + 1
    IF (sn_cfctl % l_config) THEN
      IF ((narea < sn_cfctl % procmin .OR. narea > sn_cfctl % procmax) .OR. (MOD(narea - sn_cfctl % procmin, sn_cfctl % procincr) /= 0)) CALL nemo_set_cfctl(sn_cfctl, .FALSE., .FALSE.)
    ELSE
      CALL nemo_set_cfctl(sn_cfctl, ln_ctl, .TRUE.)
    END IF
    lwm = (narea == 1)
    lwp = (narea == 1) .OR. ln_ctl
    IF (lwm) THEN
      WRITE(numond, namctl)
      WRITE(numond, namcfg)
      IF (.NOT. ln_read_cfg) THEN
        DO ji = 1, SIZE(clnam)
          IF (TRIM(clnam(ji)) /= '') WRITE(numond, *) clnam(ji)
        END DO
      END IF
    END IF
    IF (lwp) THEN
      CALL ctl_opn(numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE., narea)
      WRITE(numout, *)
      WRITE(numout, *) '   CNRS - NERC - Met OFFICE - MERCATOR-ocean - INGV - CMCC'
      WRITE(numout, *) '                       NEMO team'
      WRITE(numout, *) '            Ocean General Circulation Model'
      WRITE(numout, *) '                NEMO version 4.0  (2019) '
      WRITE(numout, *)
      WRITE(numout, *) "           ._      ._      ._      ._      ._    "
      WRITE(numout, *) "       _.-._)`\_.-._)`\_.-._)`\_.-._)`\_.-._)`\_ "
      WRITE(numout, *)
      WRITE(numout, *) "           o         _,           _,             "
      WRITE(numout, *) "            o      .' (        .-' /             "
      WRITE(numout, *) "           o     _/..._'.    .'   /              "
      WRITE(numout, *) "      (    o .-'`      ` '-./  _.'               "
      WRITE(numout, *) "       )    ( o)           ;= <_         (       "
      WRITE(numout, *) "      (      '-.,\\__ __.-;`\   '.        )      "
      WRITE(numout, *) "       )  )       \) |`\ \)  '.   \      (   (   "
      WRITE(numout, *) "      (  (           \_/       '-._\      )   )  "
      WRITE(numout, *) "       )  ) jgs                     `    (   (   "
      WRITE(numout, *) "     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "
      WRITE(numout, *)
      DO ji = 1, SIZE(cltxt)
        IF (TRIM(cltxt(ji)) /= '') WRITE(numout, *) TRIM(cltxt(ji))
      END DO
      WRITE(numout, *)
      WRITE(numout, *)
      DO ji = 1, SIZE(cltxt2)
        IF (TRIM(cltxt2(ji)) /= '') WRITE(numout, *) TRIM(cltxt2(ji))
      END DO
      WRITE(numout, cform_aaa)
    END IF
    CALL ctl_opn(numnul, '/dev/null', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE.)
    CALL mpp_init
    CALL nemo_alloc
    CALL nemo_ctl
    IF (ln_timing) CALL timing_init
    IF (ln_timing) CALL timing_start('nemo_init')
    CALL phy_cst
    CALL eos_init
    IF (lk_c1d) CALL c1d_init
    CALL wad_init
    CALL dom_init("OPA")
    IF (ln_crs) CALL crs_init
    IF (ln_ctl) CALL prt_ctl_init
    CALL diurnal_sst_bulk_init
    IF (ln_diurnal) CALL diurnal_sst_coolskin_init
    IF (ln_diurnal_only) THEN
      CALL istate_init
      CALL sbc_init
      CALL tra_qsr_init
      IF (ln_diaobs) THEN
        CALL dia_obs_init
        CALL dia_obs(nit000 - 1)
      END IF
      IF (lk_asminc) CALL asm_inc_init
      RETURN
    END IF
    CALL istate_init
    CALL tide_init
    CALL sbc_init
    CALL bdy_init
    CALL zdf_phy_init
    CALL ldf_tra_init
    CALL ldf_eiv_init
    CALL ldf_dyn_init
    IF (ln_traqsr) CALL tra_qsr_init
    CALL tra_bbc_init
    CALL tra_bbl_init
    CALL tra_dmp_init
    CALL tra_adv_init
    CALL tra_ldf_init
    IF (lk_c1d) CALL dyn_dmp_init
    CALL dyn_adv_init
    CALL dyn_vor_init
    CALL dyn_ldf_init
    CALL dyn_hpg_init
    CALL dyn_spg_init
    IF (l_ldfslp) CALL ldf_slp_init
    CALL icb_init(rdt, nit000)
    CALL sto_par_init
    IF (ln_sto_eos) CALL sto_pts_init
    IF (lk_floats) CALL flo_init
    IF (ln_diacfl) CALL dia_cfl_init
    CALL dia_ptr_init
    IF (lk_diadct) CALL dia_dct_init
    CALL dia_hsb_init
    CALL trd_init
    CALL dia_obs_init
    CALL dia_tmb_init
    CALL dia_25h_init
    IF (ln_diaobs) CALL dia_obs(nit000 - 1)
    IF (lk_asminc) CALL asm_inc_init
    IF (lwp) WRITE(numout, cform_aaa)
    IF (ln_timing) CALL timing_stop('nemo_init')
  END SUBROUTINE nemo_init
  SUBROUTINE nemo_ctl
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'nemo_ctl: Control prints'
      WRITE(numout, *) '~~~~~~~~'
      WRITE(numout, *) '   Namelist namctl'
      WRITE(numout, *) '      run control (for debugging)     ln_ctl     = ', ln_ctl
      WRITE(numout, *) '       finer control over o/p sn_cfctl%l_config  = ', sn_cfctl % l_config
      WRITE(numout, *) '                              sn_cfctl%l_runstat = ', sn_cfctl % l_runstat
      WRITE(numout, *) '                              sn_cfctl%l_trcstat = ', sn_cfctl % l_trcstat
      WRITE(numout, *) '                              sn_cfctl%l_oceout  = ', sn_cfctl % l_oceout
      WRITE(numout, *) '                              sn_cfctl%l_layout  = ', sn_cfctl % l_layout
      WRITE(numout, *) '                              sn_cfctl%l_mppout  = ', sn_cfctl % l_mppout
      WRITE(numout, *) '                              sn_cfctl%l_mpptop  = ', sn_cfctl % l_mpptop
      WRITE(numout, *) '                              sn_cfctl%procmin   = ', sn_cfctl % procmin
      WRITE(numout, *) '                              sn_cfctl%procmax   = ', sn_cfctl % procmax
      WRITE(numout, *) '                              sn_cfctl%procincr  = ', sn_cfctl % procincr
      WRITE(numout, *) '                              sn_cfctl%ptimincr  = ', sn_cfctl % ptimincr
      WRITE(numout, *) '      level of print                  nn_print   = ', nn_print
      WRITE(numout, *) '      Start i indice for SUM control  nn_ictls   = ', nn_ictls
      WRITE(numout, *) '      End i indice for SUM control    nn_ictle   = ', nn_ictle
      WRITE(numout, *) '      Start j indice for SUM control  nn_jctls   = ', nn_jctls
      WRITE(numout, *) '      End j indice for SUM control    nn_jctle   = ', nn_jctle
      WRITE(numout, *) '      number of proc. following i     nn_isplt   = ', nn_isplt
      WRITE(numout, *) '      number of proc. following j     nn_jsplt   = ', nn_jsplt
      WRITE(numout, *) '      timing by routine               ln_timing  = ', ln_timing
      WRITE(numout, *) '      CFL diagnostics                 ln_diacfl  = ', ln_diacfl
    END IF
    nprint = nn_print
    nictls = nn_ictls
    nictle = nn_ictle
    njctls = nn_jctls
    njctle = nn_jctle
    isplt = nn_isplt
    jsplt = nn_jsplt
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) '   Namelist namcfg'
      WRITE(numout, *) '      read domain configuration file                ln_read_cfg      = ', ln_read_cfg
      WRITE(numout, *) '         filename to be read                           cn_domcfg     = ', TRIM(cn_domcfg)
      WRITE(numout, *) '         keep closed seas in the domain (if exist)     ln_closea     = ', ln_closea
      WRITE(numout, *) '      create a configuration definition file        ln_write_cfg     = ', ln_write_cfg
      WRITE(numout, *) '         filename to be written                        cn_domcfg_out = ', TRIM(cn_domcfg_out)
      WRITE(numout, *) '      use file attribute if exists as i/p j-start   ln_use_jattr     = ', ln_use_jattr
    END IF
    IF (.NOT. ln_read_cfg) ln_closea = .FALSE.
    IF (ln_ctl) THEN
      IF (lk_mpp .AND. jpnij > 1) THEN
        isplt = jpni
        jsplt = jpnj
        ijsplt = jpni * jpnj
      ELSE
        IF (isplt == 1 .AND. jsplt == 1) THEN
          CALL ctl_warn(' - isplt & jsplt are equal to 1', ' - the print control will be done over the whole domain')
        END IF
        ijsplt = isplt * jsplt
      END IF
      IF (lwp) WRITE(numout, *) '          - The total number of processors over which the'
      IF (lwp) WRITE(numout, *) '            print control will be done is ijsplt : ', ijsplt
      IF (nictls + nictle + njctls + njctle == 0) THEN
        lsp_area = .FALSE.
      ELSE
        lsp_area = .TRUE.
        IF (nictls < 1 .OR. nictls > jpiglo) THEN
          CALL ctl_warn('          - nictls must be 1<=nictls>=jpiglo, it is forced to 1')
          nictls = 1
        END IF
        IF (nictle < 1 .OR. nictle > jpiglo) THEN
          CALL ctl_warn('          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo')
          nictle = jpiglo
        END IF
        IF (njctls < 1 .OR. njctls > jpjglo) THEN
          CALL ctl_warn('          - njctls must be 1<=njctls>=jpjglo, it is forced to 1')
          njctls = 1
        END IF
        IF (njctle < 1 .OR. njctle > jpjglo) THEN
          CALL ctl_warn('          - njctle must be 1<=njctle>=jpjglo, it is forced to jpjglo')
          njctle = jpjglo
        END IF
      END IF
    END IF
    IF (1._wp /= SIGN(1._wp, - 0._wp)) CALL ctl_stop('nemo_ctl: The intrinsec SIGN function follows f2003 standard.', 'Compile with key_nosignedzero enabled:', '--> add -Dkey_nosignedzero to the definition of %CPP in your arch file')
  END SUBROUTINE nemo_ctl
  SUBROUTINE nemo_closefile
    IF (lk_mpp) CALL mppsync
    CALL iom_close
    IF (numstp /= - 1) CLOSE(UNIT = numstp)
    IF (numrun /= - 1) CLOSE(UNIT = numrun)
    IF (numnam_ref /= - 1) CLOSE(UNIT = numnam_ref)
    IF (numnam_cfg /= - 1) CLOSE(UNIT = numnam_cfg)
    IF (lwm .AND. numond /= - 1) CLOSE(UNIT = numond)
    IF (numnam_ice_ref /= - 1) CLOSE(UNIT = numnam_ice_ref)
    IF (numnam_ice_cfg /= - 1) CLOSE(UNIT = numnam_ice_cfg)
    IF (lwm .AND. numoni /= - 1) CLOSE(UNIT = numoni)
    IF (numevo_ice /= - 1) CLOSE(UNIT = numevo_ice)
    IF (numout /= 6) CLOSE(UNIT = numout)
    IF (numdct_vol /= - 1) CLOSE(UNIT = numdct_vol)
    IF (numdct_heat /= - 1) CLOSE(UNIT = numdct_heat)
    IF (numdct_salt /= - 1) CLOSE(UNIT = numdct_salt)
    numout = 6
  END SUBROUTINE nemo_closefile
  SUBROUTINE nemo_alloc
    USE diawri, ONLY: dia_wri_alloc
    USE dom_oce, ONLY: dom_oce_alloc
    USE trc_oce, ONLY: trc_oce_alloc
    USE bdy_oce, ONLY: bdy_oce_alloc
    INTEGER :: ierr
    ierr = oce_alloc()
    ierr = ierr + dia_wri_alloc()
    ierr = ierr + dom_oce_alloc()
    ierr = ierr + zdf_oce_alloc()
    ierr = ierr + trc_oce_alloc()
    ierr = ierr + bdy_oce_alloc()
    CALL mpp_sum('nemogcm', ierr)
    IF (ierr /= 0) CALL ctl_stop('STOP', 'nemo_alloc: unable to allocate standard ocean arrays')
  END SUBROUTINE nemo_alloc
  SUBROUTINE nemo_set_cfctl(sn_cfctl, setto, for_all)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: setto, for_all
    TYPE(sn_ctl) :: sn_cfctl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('nemo_set_cfctl', 'r0', 0, 0)
    IF (for_all) THEN
      sn_cfctl % l_runstat = setto
      sn_cfctl % l_trcstat = setto
    END IF
    sn_cfctl % l_oceout = setto
    sn_cfctl % l_layout = setto
    sn_cfctl % l_mppout = setto
    sn_cfctl % l_mpptop = setto
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE nemo_set_cfctl
END MODULE nemogcm