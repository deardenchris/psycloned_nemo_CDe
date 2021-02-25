MODULE domain
  USE oce
  USE dom_oce
  USE sbc_oce
  USE trc_oce
  USE phycst
  USE closea
  USE domhgr
  USE domzgr
  USE dommsk
  USE domwri
  USE domvvl
  USE c1d
  USE dyncor_c1d
  USE wet_dry, ONLY: ll_wd
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_init
  PUBLIC :: domain_cfg
  CONTAINS
  SUBROUTINE dom_init(cdstr)
    INTEGER :: ji, jj, jk, ik
    INTEGER :: iconf = 0
    CHARACTER(LEN = 64) :: cform = "(A12, 3(A13, I7))"
    CHARACTER(LEN = *), INTENT(IN) :: cdstr
    INTEGER, DIMENSION(jpi, jpj) :: ik_top, ik_bot
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z1_hu_0, z1_hv_0
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dom_init : domain initialization'
      WRITE(numout, FMT = *) '~~~~~~~~'
      WRITE(numout, FMT = *) '   Domain info'
      WRITE(numout, FMT = *) '      dimension of model:'
      WRITE(numout, FMT = *) '             Local domain      Global domain       Data domain '
      WRITE(numout, cform) '        ', '   jpi     : ', jpi, '   jpiglo  : ', jpiglo
      WRITE(numout, cform) '        ', '   jpj     : ', jpj, '   jpjglo  : ', jpjglo
      WRITE(numout, cform) '        ', '   jpk     : ', jpk, '   jpkglo  : ', jpkglo
      WRITE(numout, cform) '       ', '   jpij    : ', jpij
      WRITE(numout, FMT = *) '      mpp local domain info (mpp):'
      WRITE(numout, FMT = *) '              jpni    : ', jpni, '   nn_hls  : ', nn_hls
      WRITE(numout, FMT = *) '              jpnj    : ', jpnj, '   nn_hls  : ', nn_hls
      WRITE(numout, FMT = *) '              jpnij   : ', jpnij
      WRITE(numout, FMT = *) '      lateral boundary of the Global domain : jperio  = ', jperio
      SELECT CASE (jperio)
      CASE (0)
        WRITE(numout, FMT = *) '         (i.e. closed)'
      CASE (1)
        WRITE(numout, FMT = *) '         (i.e. cyclic east-west)'
      CASE (2)
        WRITE(numout, FMT = *) '         (i.e. equatorial symmetric)'
      CASE (3)
        WRITE(numout, FMT = *) '         (i.e. north fold with T-point pivot)'
      CASE (4)
        WRITE(numout, FMT = *) '         (i.e. cyclic east-west and north fold with T-point pivot)'
      CASE (5)
        WRITE(numout, FMT = *) '         (i.e. north fold with F-point pivot)'
      CASE (6)
        WRITE(numout, FMT = *) '         (i.e. cyclic east-west and north fold with F-point pivot)'
      CASE (7)
        WRITE(numout, FMT = *) '         (i.e. cyclic east-west and north-south)'
      CASE DEFAULT
        CALL ctl_stop('jperio is out of range')
      END SELECT
      WRITE(numout, FMT = *) '      Ocean model configuration used:'
      WRITE(numout, FMT = *) '         cn_cfg = ', TRIM(cn_cfg), '   nn_cfg = ', nn_cfg
    END IF
    lwxios = .FALSE.
    ln_xios_read = .FALSE.
    CALL dom_glo
    CALL dom_nam
    IF (lwxios) THEN
      CALL iom_set_rst_vars(rst_wfields)
      CALL iom_set_rstw_core(cdstr)
    END IF
    IF (cdstr == 'SAS') THEN
      IF (lrxios) THEN
        IF (lwp) WRITE(numout, FMT = *) 'Disable reading restart file using XIOS for SAS'
        lrxios = .FALSE.
      END IF
    END IF
    CALL dom_hgr
    CALL dom_zgr(ik_top, ik_bot)
    CALL dom_msk(ik_top, ik_bot)
    IF (ln_closea) CALL dom_clo
    DO jj = 1, jpj
      DO ji = 1, jpi
        ik = mikt(ji, jj)
        risfdep(ji, jj) = gdepw_0(ji, jj, ik)
      END DO
    END DO
    ht_0(:, :) = 0._wp
    hu_0(:, :) = 0._wp
    hv_0(:, :) = 0._wp
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpk
      ht_0(:, :) = ht_0(:, :) + e3t_0(:, :, jk) * tmask(:, :, jk)
      hu_0(:, :) = hu_0(:, :) + e3u_0(:, :, jk) * umask(:, :, jk)
      hv_0(:, :) = hv_0(:, :) + e3v_0(:, :, jk) * vmask(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (ln_linssh) THEN
      gdept_b = gdept_0
      gdept_n = gdept_0
      gdepw_b = gdepw_0
      gdepw_n = gdepw_0
      gde3w_n = gde3w_0
      e3t_b = e3t_0
      e3t_n = e3t_0
      e3t_a = e3t_0
      e3u_b = e3u_0
      e3u_n = e3u_0
      e3u_a = e3u_0
      e3v_b = e3v_0
      e3v_n = e3v_0
      e3v_a = e3v_0
      e3f_n = e3f_0
      e3w_b = e3w_0
      e3w_n = e3w_0
      e3uw_b = e3uw_0
      e3uw_n = e3uw_0
      e3vw_b = e3vw_0
      e3vw_n = e3vw_0
      z1_hu_0(:, :) = ssumask(:, :) / (hu_0(:, :) + 1._wp - ssumask(:, :))
      z1_hv_0(:, :) = ssvmask(:, :) / (hv_0(:, :) + 1._wp - ssvmask(:, :))
      ht_n = ht_0
      hu_b = hu_0
      hu_n = hu_0
      hu_a = hu_0
      hv_b = hv_0
      hv_n = hv_0
      hv_a = hv_0
      r1_hu_b = z1_hu_0
      r1_hu_n = z1_hu_0
      r1_hu_a = z1_hu_0
      r1_hv_b = z1_hv_0
      r1_hv_n = z1_hv_0
      r1_hv_a = z1_hv_0
    ELSE
      IF (.NOT. l_offline) CALL dom_vvl_init
    END IF
    IF (lk_c1d) CALL cor_c1d
    IF (ln_meshmask .AND. .NOT. ln_iscpl) CALL dom_wri
    IF (ln_meshmask .AND. ln_iscpl .AND. .NOT. ln_rstart) CALL dom_wri
    IF (.NOT. ln_rstart) CALL dom_ctl
    IF (ln_write_cfg) CALL cfg_write
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dom_init :   ==>>>   END of domain initialization'
      WRITE(numout, FMT = *) '~~~~~~~~'
      WRITE(numout, FMT = *)
    END IF
  END SUBROUTINE dom_init
  SUBROUTINE dom_glo
    INTEGER :: ji, jj
    DO ji = 1, jpi
      mig(ji) = ji + nimpp - 1
    END DO
    DO jj = 1, jpj
      mjg(jj) = jj + njmpp - 1
    END DO
    DO ji = 1, jpiglo
      mi0(ji) = MAX(1, MIN(ji - nimpp + 1, jpi + 1))
      mi1(ji) = MAX(0, MIN(ji - nimpp + 1, jpi))
    END DO
    DO jj = 1, jpjglo
      mj0(jj) = MAX(1, MIN(jj - njmpp + 1, jpj + 1))
      mj1(jj) = MAX(0, MIN(jj - njmpp + 1, jpj))
    END DO
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dom_glo : domain: global <<==>> local '
      WRITE(numout, FMT = *) '~~~~~~~ '
      WRITE(numout, FMT = *) '   global domain:   jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpkglo = ', jpkglo
      WRITE(numout, FMT = *) '   local  domain:   jpi    = ', jpi, ' jpj    = ', jpj, ' jpk    = ', jpk
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   conversion from local to global domain indices (and vise versa) done'
      IF (nn_print >= 1) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '          conversion local  ==> global i-index domain (mig)'
        WRITE(numout, 25) (mig(ji), ji = 1, jpi)
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '          conversion global ==> local  i-index domain'
        WRITE(numout, FMT = *) '             starting index (mi0)'
        WRITE(numout, 25) (mi0(ji), ji = 1, jpiglo)
        WRITE(numout, FMT = *) '             ending index (mi1)'
        WRITE(numout, 25) (mi1(ji), ji = 1, jpiglo)
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '          conversion local  ==> global j-index domain (mjg)'
        WRITE(numout, 25) (mjg(jj), jj = 1, jpj)
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '          conversion global ==> local  j-index domain'
        WRITE(numout, FMT = *) '             starting index (mj0)'
        WRITE(numout, 25) (mj0(jj), jj = 1, jpjglo)
        WRITE(numout, FMT = *) '             ending index (mj1)'
        WRITE(numout, 25) (mj1(jj), jj = 1, jpjglo)
      END IF
    END IF
25  FORMAT(100(10X, 19I4, /))
  END SUBROUTINE dom_glo
  SUBROUTINE dom_nam
    USE ioipsl
    INTEGER :: ios
    NAMELIST /namrun/ cn_ocerst_indir, cn_ocerst_outdir, nn_stocklist, ln_rst_list, nn_no, cn_exp, cn_ocerst_in, cn_ocerst_out, &
&ln_rstart, nn_rstctl, nn_it000, nn_itend, nn_date0, nn_time0, nn_leapy, nn_istate, nn_stock, nn_write, ln_mskland, ln_clobber, &
&nn_chunksz, nn_euler, ln_cfmeta, ln_iscpl, ln_xios_read, nn_wxios, ln_rstdate
    NAMELIST /namdom/ ln_linssh, rn_isfhmin, rn_rdt, rn_atfp, ln_crs, ln_meshmask
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dom_nam : domain initialization through namelist read'
      WRITE(numout, FMT = *) '~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namrun, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namrun in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namrun, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namrun in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namrun)
    IF (lwp) THEN
      WRITE(numout, FMT = *) '   Namelist : namrun   ---   run parameters'
      WRITE(numout, FMT = *) '      Assimilation cycle              nn_no           = ', nn_no
      WRITE(numout, FMT = *) '      experiment name for output      cn_exp          = ', TRIM(cn_exp)
      WRITE(numout, FMT = *) '      file prefix restart input       cn_ocerst_in    = ', TRIM(cn_ocerst_in)
      WRITE(numout, FMT = *) '      restart input directory         cn_ocerst_indir = ', TRIM(cn_ocerst_indir)
      WRITE(numout, FMT = *) '      file prefix restart output      cn_ocerst_out   = ', TRIM(cn_ocerst_out)
      WRITE(numout, FMT = *) '      restart output directory        cn_ocerst_outdir= ', TRIM(cn_ocerst_outdir)
      WRITE(numout, FMT = *) '      restart logical                 ln_rstart       = ', ln_rstart
      WRITE(numout, FMT = *) '      start with forward time step    nn_euler        = ', nn_euler
      WRITE(numout, FMT = *) '      control of time step            nn_rstctl       = ', nn_rstctl
      WRITE(numout, FMT = *) '      number of the first time step   nn_it000        = ', nn_it000
      WRITE(numout, FMT = *) '      number of the last time step    nn_itend        = ', nn_itend
      WRITE(numout, FMT = *) '      initial calendar date aammjj    nn_date0        = ', nn_date0
      WRITE(numout, FMT = *) '      initial time of day in hhmm     nn_time0        = ', nn_time0
      WRITE(numout, FMT = *) '      leap year calendar (0/1)        nn_leapy        = ', nn_leapy
      WRITE(numout, FMT = *) '      initial state output            nn_istate       = ', nn_istate
      IF (ln_rst_list) THEN
        WRITE(numout, FMT = *) '      list of restart dump times      nn_stocklist    =', nn_stocklist
      ELSE
        WRITE(numout, FMT = *) '      frequency of restart file       nn_stock        = ', nn_stock
      END IF
      WRITE(numout, FMT = *) '      frequency of output file        nn_write        = ', nn_write
      WRITE(numout, FMT = *) '      mask land points                ln_mskland      = ', ln_mskland
      WRITE(numout, FMT = *) '      date-stamp restart files        ln_rstdate = ', ln_rstdate
      WRITE(numout, FMT = *) '      additional CF standard metadata ln_cfmeta       = ', ln_cfmeta
      WRITE(numout, FMT = *) '      overwrite an existing file      ln_clobber      = ', ln_clobber
      WRITE(numout, FMT = *) '      NetCDF chunksize (bytes)        nn_chunksz      = ', nn_chunksz
      WRITE(numout, FMT = *) '      IS coupling at the restart step ln_iscpl        = ', ln_iscpl
      IF (TRIM(Agrif_CFixed()) == '0') THEN
        WRITE(numout, FMT = *) '      READ restart for a single file using XIOS ln_xios_read =', ln_xios_read
        WRITE(numout, FMT = *) '      Write restart using XIOS        nn_wxios   = ', nn_wxios
      ELSE
        WRITE(numout, FMT = *) "      AGRIF: nn_wxios will be ingored. See setting for parent"
        WRITE(numout, FMT = *) "      AGRIF: ln_xios_read will be ingored. See setting for parent"
      END IF
    END IF
    cexper = cn_exp
    nrstdt = nn_rstctl
    nit000 = nn_it000
    nitend = nn_itend
    ndate0 = nn_date0
    nleapy = nn_leapy
    ninist = nn_istate
    nstock = nn_stock
    nstocklist = nn_stocklist
    nwrite = nn_write
    neuler = nn_euler
    IF (neuler == 1 .AND. .NOT. ln_rstart) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Start from rest (ln_rstart=F)'
      IF (lwp) WRITE(numout, FMT = *) '           an Euler initial time step is used : nn_euler is forced to 0 '
      neuler = 0
    END IF
    IF (nstock == 0 .OR. nstock > nitend) THEN
      WRITE(ctmp1, FMT = *) 'nstock = ', nstock, ' it is forced to ', nitend
      CALL ctl_warn(ctmp1)
      nstock = nitend
    END IF
    IF (nwrite == 0) THEN
      WRITE(ctmp1, FMT = *) 'nwrite = ', nwrite, ' it is forced to ', nitend
      CALL ctl_warn(ctmp1)
      nwrite = nitend
    END IF
    IF (lwp) WRITE(numout, FMT = *)
    SELECT CASE (nleapy)
    CASE (1)
      CALL ioconf_calendar('gregorian')
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   The IOIPSL calendar is "gregorian", i.e. leap year'
    CASE (0)
      CALL ioconf_calendar('noleap')
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   The IOIPSL calendar is "noleap", i.e. no leap year'
    CASE (30)
      CALL ioconf_calendar('360d')
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   The IOIPSL calendar is "360d", i.e. 360 days in a year'
    END SELECT
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdom, IOSTAT = ios, ERR = 903)
903 IF (ios /= 0) CALL ctl_nam(ios, 'namdom in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdom, IOSTAT = ios, ERR = 904)
904 IF (ios > 0) CALL ctl_nam(ios, 'namdom in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdom)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Namelist : namdom   ---   space & time domain'
      WRITE(numout, FMT = *) '      linear free surface (=T)                ln_linssh   = ', ln_linssh
      WRITE(numout, FMT = *) '      create mesh/mask file                   ln_meshmask = ', ln_meshmask
      WRITE(numout, FMT = *) '      treshold to open the isf cavity         rn_isfhmin  = ', rn_isfhmin, ' [m]'
      WRITE(numout, FMT = *) '      ocean time step                         rn_rdt      = ', rn_rdt
      WRITE(numout, FMT = *) '      asselin time filter parameter           rn_atfp     = ', rn_atfp
      WRITE(numout, FMT = *) '      online coarsening of dynamical fields   ln_crs      = ', ln_crs
    END IF
    atfp = rn_atfp
    rdt = rn_rdt
    IF (TRIM(Agrif_CFixed()) == '0') THEN
      lrxios = ln_xios_read .AND. ln_rstart
      IF (nn_wxios > 0) lwxios = .TRUE.
      nxioso = nn_wxios
    END IF
    snc4set % luse = .FALSE.
  END SUBROUTINE dom_nam
  SUBROUTINE dom_ctl
    INTEGER, DIMENSION(2) :: imi1, imi2, ima1, ima2
    INTEGER, DIMENSION(2) :: iloc
    REAL(KIND = wp) :: ze1min, ze1max, ze2min, ze2max
    IF (lk_mpp) THEN
      CALL mpp_minloc('domain', e1t(:, :), tmask_i(:, :), ze1min, imi1)
      CALL mpp_minloc('domain', e2t(:, :), tmask_i(:, :), ze2min, imi2)
      CALL mpp_maxloc('domain', e1t(:, :), tmask_i(:, :), ze1max, ima1)
      CALL mpp_maxloc('domain', e2t(:, :), tmask_i(:, :), ze2max, ima2)
    ELSE
      ze1min = MINVAL(e1t(:, :), mask = tmask_i(:, :) == 1._wp)
      ze2min = MINVAL(e2t(:, :), mask = tmask_i(:, :) == 1._wp)
      ze1max = MAXVAL(e1t(:, :), mask = tmask_i(:, :) == 1._wp)
      ze2max = MAXVAL(e2t(:, :), mask = tmask_i(:, :) == 1._wp)
      iloc = MINLOC(e1t(:, :), mask = tmask_i(:, :) == 1._wp)
      imi1(1) = iloc(1) + nimpp - 1
      imi1(2) = iloc(2) + njmpp - 1
      iloc = MINLOC(e2t(:, :), mask = tmask_i(:, :) == 1._wp)
      imi2(1) = iloc(1) + nimpp - 1
      imi2(2) = iloc(2) + njmpp - 1
      iloc = MAXLOC(e1t(:, :), mask = tmask_i(:, :) == 1._wp)
      ima1(1) = iloc(1) + nimpp - 1
      ima1(2) = iloc(2) + njmpp - 1
      iloc = MAXLOC(e2t(:, :), mask = tmask_i(:, :) == 1._wp)
      ima2(1) = iloc(1) + nimpp - 1
      ima2(2) = iloc(2) + njmpp - 1
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dom_ctl : extrema of the masked scale factors'
      WRITE(numout, FMT = *) '~~~~~~~'
      WRITE(numout, "(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, ima1(1), ima1(2)
      WRITE(numout, "(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, imi1(1), imi1(2)
      WRITE(numout, "(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, ima2(1), ima2(2)
      WRITE(numout, "(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, imi2(1), imi2(2)
    END IF
  END SUBROUTINE dom_ctl
  SUBROUTINE domain_cfg(ldtxt, cd_cfg, kk_cfg, kpi, kpj, kpk, kperio)
    CHARACTER(LEN = *), DIMENSION(:), INTENT(OUT) :: ldtxt
    CHARACTER(LEN = *), INTENT(OUT) :: cd_cfg
    INTEGER, INTENT(OUT) :: kk_cfg
    INTEGER, INTENT(OUT) :: kpi, kpj, kpk
    INTEGER, INTENT(OUT) :: kperio
    INTEGER :: inum, ii
    REAL(KIND = wp) :: zorca_res
    REAL(KIND = wp) :: ziglo, zjglo, zkglo, zperio
    ii = 1
    WRITE(ldtxt(ii), FMT = *) '           '
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) 'domain_cfg : domain size read in ', TRIM(cn_domcfg), ' file'
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '~~~~~~~~~~ '
    ii = ii + 1
    CALL iom_open(cn_domcfg, inum)
    IF (iom_varid(inum, 'ORCA', ldstop = .FALSE.) > 0 .AND. iom_varid(inum, 'ORCA_index', ldstop = .FALSE.) > 0) THEN
      cd_cfg = 'ORCA'
      CALL iom_get(inum, 'ORCA_index', zorca_res)
      kk_cfg = NINT(zorca_res)
      WRITE(ldtxt(ii), FMT = *) '   .'
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '   ==>>>   ORCA configuration '
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '   .'
      ii = ii + 1
    ELSE
      cd_cfg = 'UNKNOWN'
      kk_cfg = - 9999999
      CALL iom_getatt(inum, 'cn_cfg', cd_cfg)
      CALL iom_getatt(inum, 'nn_cfg', kk_cfg)
      IF (TRIM(cd_cfg) == '!') cd_cfg = 'UNKNOWN'
      IF (kk_cfg == - 999) kk_cfg = - 9999999
    END IF
    CALL iom_get(inum, 'jpiglo', ziglo)
    kpi = NINT(ziglo)
    CALL iom_get(inum, 'jpjglo', zjglo)
    kpj = NINT(zjglo)
    CALL iom_get(inum, 'jpkglo', zkglo)
    kpk = NINT(zkglo)
    CALL iom_get(inum, 'jperio', zperio)
    kperio = NINT(zperio)
    CALL iom_close(inum)
    WRITE(ldtxt(ii), FMT = *) '      cn_cfg = ', TRIM(cd_cfg), '   nn_cfg = ', kk_cfg
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      jpiglo = ', kpi
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      jpjglo = ', kpj
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      jpkglo = ', kpk
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      type of global domain lateral boundary   jperio = ', kperio
    ii = ii + 1
  END SUBROUTINE domain_cfg
  SUBROUTINE cfg_write
    INTEGER :: ji, jj, jk
    INTEGER :: izco, izps, isco, icav
    INTEGER :: inum
    CHARACTER(LEN = 21) :: clnam
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'cfg_write : create the domain configuration file (', TRIM(cn_domcfg_out), '.nc)'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~'
    clnam = cn_domcfg_out
    CALL iom_open(TRIM(clnam), inum, ldwrt = .TRUE.)
    IF (cn_cfg == "ORCA") THEN
      CALL iom_rstput(0, 0, inum, 'ORCA', 1._wp, ktype = jp_i4)
      CALL iom_rstput(0, 0, inum, 'ORCA_index', REAL(nn_cfg, wp), ktype = jp_i4)
    END IF
    CALL iom_rstput(0, 0, inum, 'jpiglo', REAL(jpiglo, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'jpjglo', REAL(jpjglo, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'jpkglo', REAL(jpk, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'jperio', REAL(jperio, wp), ktype = jp_i4)
    IF (ln_zco) THEN
      izco = 1
    ELSE
      izco = 0
    END IF
    IF (ln_zps) THEN
      izps = 1
    ELSE
      izps = 0
    END IF
    IF (ln_sco) THEN
      isco = 1
    ELSE
      isco = 0
    END IF
    CALL iom_rstput(0, 0, inum, 'ln_zco', REAL(izco, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'ln_zps', REAL(izps, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'ln_sco', REAL(isco, wp), ktype = jp_i4)
    IF (ln_isfcav) THEN
      icav = 1
    ELSE
      icav = 0
    END IF
    CALL iom_rstput(0, 0, inum, 'ln_isfcav', REAL(icav, wp), ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'glamt', glamt, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'glamu', glamu, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'glamv', glamv, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'glamf', glamf, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphit', gphit, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphiu', gphiu, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphiv', gphiv, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'gphif', gphif, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1t', e1t, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1u', e1u, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1v', e1v, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e1f', e1f, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2t', e2t, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2u', e2u, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2v', e2v, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e2f', e2f, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'ff_f', ff_f, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'ff_t', ff_t, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3t_1d', e3t_1d, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3w_1d', e3w_1d, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3t_0', e3t_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3u_0', e3u_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3v_0', e3v_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3f_0', e3f_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3w_0', e3w_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3uw_0', e3uw_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'e3vw_0', e3vw_0, ktype = jp_r8)
    CALL iom_rstput(0, 0, inum, 'top_level', REAL(mikt, wp) * ssmask, ktype = jp_i4)
    CALL iom_rstput(0, 0, inum, 'bottom_level', REAL(mbkt, wp) * ssmask, ktype = jp_i4)
    IF (ln_sco) THEN
      CALL dom_stiff(z2d)
      CALL iom_rstput(0, 0, inum, 'stiffness', z2d)
    END IF
    IF (ll_wd) THEN
      CALL iom_rstput(0, 0, inum, 'ht_0', ht_0, ktype = jp_r8)
    END IF
    CALL iom_putatt(inum, 'nn_cfg', nn_cfg)
    CALL iom_putatt(inum, 'cn_cfg', TRIM(cn_cfg))
    CALL iom_close(inum)
  END SUBROUTINE cfg_write
END MODULE domain