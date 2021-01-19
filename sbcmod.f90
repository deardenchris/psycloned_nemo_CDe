MODULE sbcmod
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE trc_oce
  USE sbc_ice
  USE sbcdcy
  USE sbcssm
  USE sbcflx
  USE sbcblk
  USE sbcice_if
  USE sbcice_cice
  USE sbcisf
  USE sbccpl
  USE cpl_oasis3
  USE sbcssr
  USE sbcrnf
  USE sbcapr
  USE sbcisf
  USE sbcfwb
  USE icbstp
  USE icb_oce, ONLY: ln_passive_mode
  USE traqsr
  USE sbcwave
  USE bdy_oce, ONLY: ln_bdy
  USE usrdef_sbc
  USE closea
  USE prtctl
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE timing
  USE wet_dry
  USE diurnal_bulk, ONLY: ln_diurnal_only
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc
  PUBLIC :: sbc_init
  INTEGER :: nsbc
  CONTAINS
  SUBROUTINE sbc_init
    INTEGER :: ios, icpt
    LOGICAL :: ll_purecpl, ll_opa, ll_not_nemo
    NAMELIST /namsbc/ nn_fsbc, ln_usr, ln_flx, ln_blk, ln_cpl, ln_mixcpl, nn_components, nn_ice, ln_ice_embd, ln_traqsr, ln_dm2dc, &
&ln_rnf, nn_fwb, ln_ssr, ln_isf, ln_apr_dyn, ln_wave, ln_cdgw, ln_sdw, ln_tauwoc, ln_stcor, ln_tauw, nn_lsm, nn_sdrift
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'sbc_init : surface boundary condition setting'
      WRITE(numout, FMT = *) '~~~~~~~~ '
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc)
    ncom_fsbc = nn_fsbc
    IF (lk_si3) nn_ice = 2
    IF (lk_cice) nn_ice = 3
    IF (lwp) THEN
      WRITE(numout, FMT = *) '   Namelist namsbc (partly overwritten with CPP key setting)'
      WRITE(numout, FMT = *) '      frequency update of sbc (and ice)             nn_fsbc       = ', nn_fsbc
      WRITE(numout, FMT = *) '      Type of air-sea fluxes : '
      WRITE(numout, FMT = *) '         user defined formulation                   ln_usr        = ', ln_usr
      WRITE(numout, FMT = *) '         flux         formulation                   ln_flx        = ', ln_flx
      WRITE(numout, FMT = *) '         bulk         formulation                   ln_blk        = ', ln_blk
      WRITE(numout, FMT = *) '      Type of coupling (Ocean/Ice/Atmosphere) : '
      WRITE(numout, FMT = *) '         ocean-atmosphere coupled formulation       ln_cpl        = ', ln_cpl
      WRITE(numout, FMT = *) '         mixed forced-coupled     formulation       ln_mixcpl     = ', ln_mixcpl
      WRITE(numout, FMT = *) '         OASIS coupling (with atm or sas)           lk_oasis      = ', lk_oasis
      WRITE(numout, FMT = *) '         components of your executable              nn_components = ', nn_components
      WRITE(numout, FMT = *) '      Sea-ice : '
      WRITE(numout, FMT = *) '         ice management in the sbc (=0/1/2/3)       nn_ice        = ', nn_ice
      WRITE(numout, FMT = *) '         ice embedded into ocean                    ln_ice_embd   = ', ln_ice_embd
      WRITE(numout, FMT = *) '      Misc. options of sbc : '
      WRITE(numout, FMT = *) '         Light penetration in temperature Eq.       ln_traqsr     = ', ln_traqsr
      WRITE(numout, FMT = *) '            daily mean to diurnal cycle qsr            ln_dm2dc   = ', ln_dm2dc
      WRITE(numout, FMT = *) '         Sea Surface Restoring on SST and/or SSS    ln_ssr        = ', ln_ssr
      WRITE(numout, FMT = *) '         FreshWater Budget control  (=0/1/2)        nn_fwb        = ', nn_fwb
      WRITE(numout, FMT = *) '         Patm gradient added in ocean & ice Eqs.    ln_apr_dyn    = ', ln_apr_dyn
      WRITE(numout, FMT = *) '         runoff / runoff mouths                     ln_rnf        = ', ln_rnf
      WRITE(numout, FMT = *) '         iceshelf formulation                       ln_isf        = ', ln_isf
      WRITE(numout, FMT = *) '         nb of iterations if land-sea-mask applied  nn_lsm        = ', nn_lsm
      WRITE(numout, FMT = *) '         surface wave                               ln_wave       = ', ln_wave
      WRITE(numout, FMT = *) '               Stokes drift corr. to vert. velocity ln_sdw        = ', ln_sdw
      WRITE(numout, FMT = *) '                  vertical parametrization          nn_sdrift     = ', nn_sdrift
      WRITE(numout, FMT = *) '               wave modified ocean stress           ln_tauwoc     = ', ln_tauwoc
      WRITE(numout, FMT = *) '               wave modified ocean stress component ln_tauw       = ', ln_tauw
      WRITE(numout, FMT = *) '               Stokes coriolis term                 ln_stcor      = ', ln_stcor
      WRITE(numout, FMT = *) '               neutral drag coefficient (CORE,NCAR) ln_cdgw       = ', ln_cdgw
    END IF
    IF (.NOT. ln_wave) THEN
      ln_sdw = .FALSE.
      ln_cdgw = .FALSE.
      ln_tauwoc = .FALSE.
      ln_tauw = .FALSE.
      ln_stcor = .FALSE.
    END IF
    IF (ln_sdw) THEN
      IF (.NOT. (nn_sdrift == jp_breivik_2014 .OR. nn_sdrift == jp_li_2017 .OR. nn_sdrift == jp_peakfr)) CALL ctl_stop('The chosen &
&nn_sdrift for Stokes drift vertical velocity must be 0, 1, or 2')
    END IF
    ll_st_bv2014 = (nn_sdrift == jp_breivik_2014)
    ll_st_li2017 = (nn_sdrift == jp_li_2017)
    ll_st_bv_li = (ll_st_bv2014 .OR. ll_st_li2017)
    ll_st_peakfr = (nn_sdrift == jp_peakfr)
    IF (ln_tauwoc .AND. ln_tauw) CALL ctl_stop('More than one method for modifying the ocean stress has been selected ', &
&'(ln_tauwoc=.true. and ln_tauw=.true.)')
    IF (ln_tauwoc) CALL ctl_warn('You are subtracting the wave stress to the ocean (ln_tauwoc=.true.)')
    IF (ln_tauw) CALL ctl_warn('The wave modified ocean stress components are used (ln_tauw=.true.) ', 'This will override any &
&other specification of the ocean stress')
    IF (.NOT. ln_usr) THEN
      IF (MOD(rday, rdt) /= 0.) CALL ctl_stop('the time step must devide the number of second of in a day')
      IF (MOD(rday, 2.) /= 0.) CALL ctl_stop('the number of second of in a day must be an even number')
      IF (MOD(rdt, 2.) /= 0.) CALL ctl_stop('the time step (in second) must be an even number')
    END IF
    IF (lwp) WRITE(numout, FMT = *)
    SELECT CASE (nn_components)
    CASE (jp_iam_nemo)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   NEMO configured as a single executable (i.e. including both OPA and Surface &
&module)'
    CASE (jp_iam_opa)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Multi executable configuration. Here, OPA component'
      IF (.NOT. lk_oasis) CALL ctl_stop('sbc_init : OPA-SAS coupled via OASIS, but key_oasis3 disabled')
      IF (ln_cpl) CALL ctl_stop('sbc_init : OPA-SAS coupled via OASIS, but ln_cpl = T in OPA')
      IF (ln_mixcpl) CALL ctl_stop('sbc_init : OPA-SAS coupled via OASIS, but ln_mixcpl = T in OPA')
    CASE (jp_iam_sas)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Multi executable configuration. Here, SAS component'
      IF (.NOT. lk_oasis) CALL ctl_stop('sbc_init : OPA-SAS coupled via OASIS, but key_oasis3 disabled')
      IF (ln_mixcpl) CALL ctl_stop('sbc_init : OPA-SAS coupled via OASIS, but ln_mixcpl = T in OPA')
    CASE DEFAULT
      CALL ctl_stop('sbc_init : unsupported value for nn_components')
    END SELECT
    IF (ln_cpl) THEN
      IF (.NOT. lk_oasis) CALL ctl_stop('sbc_init : coupled mode with an atmosphere model (ln_cpl=T)', '           required to &
&defined key_oasis3')
    END IF
    IF (ln_mixcpl) THEN
      IF (.NOT. lk_oasis) CALL ctl_stop('sbc_init : mixed forced-coupled mode (ln_mixcpl=T) ', '           required to defined &
&key_oasis3')
      IF (.NOT. ln_cpl) CALL ctl_stop('sbc_init : mixed forced-coupled mode (ln_mixcpl=T) requires ln_cpl = T')
      IF (nn_components /= jp_iam_nemo) CALL ctl_stop('sbc_init : the mixed forced-coupled mode (ln_mixcpl=T) ', '          not &
&yet working with sas-opa coupling via oasis')
    END IF
    SELECT CASE (nn_ice)
    CASE (0)
    CASE (1)
    CASE (2)
    CASE (3)
      IF (.NOT. (ln_blk .OR. ln_cpl)) CALL ctl_stop('sbc_init : CICE sea-ice model requires ln_blk or ln_cpl = T')
      IF (lk_agrif) CALL ctl_stop('sbc_init : CICE sea-ice model not currently available with AGRIF')
    CASE DEFAULT
    END SELECT
    IF (sbc_oce_alloc() /= 0) CALL ctl_stop('sbc_init : unable to allocate sbc_oce arrays')
    IF (sbc_ice_alloc() /= 0) CALL ctl_stop('sbc_init : unable to allocate sbc_ice arrays')
    IF (.NOT. ln_isf) THEN
      IF (sbc_isf_alloc() /= 0) CALL ctl_stop('STOP', 'sbc_init : unable to allocate sbc_isf arrays')
      !$ACC KERNELS
      fwfisf(:, :) = 0._wp
      risf_tsc(:, :, :) = 0._wp
      fwfisf_b(:, :) = 0._wp
      risf_tsc_b(:, :, :) = 0._wp
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    IF (nn_ice == 0) THEN
      IF (nn_components /= jp_iam_opa) fr_i(:, :) = 0._wp
    END IF
    sfx(:, :) = 0._wp
    fmmflx(:, :) = 0._wp
    taum(:, :) = 0._wp
    !$ACC END KERNELS
    IF (ln_dm2dc) THEN
      nday_qsr = - 1
      IF (.NOT. (ln_flx .OR. ln_blk) .AND. nn_components /= jp_iam_opa) CALL ctl_stop('qsr diurnal cycle from daily values &
&requires a flux or bulk formulation')
    END IF
    ll_purecpl = ln_cpl .AND. .NOT. ln_mixcpl
    ll_opa = nn_components == jp_iam_opa
    ll_not_nemo = nn_components /= jp_iam_nemo
    icpt = 0
    IF (ln_usr) THEN
      nsbc = jp_usr
      icpt = icpt + 1
    END IF
    IF (ln_flx) THEN
      nsbc = jp_flx
      icpt = icpt + 1
    END IF
    IF (ln_blk) THEN
      nsbc = jp_blk
      icpt = icpt + 1
    END IF
    IF (ll_purecpl) THEN
      nsbc = jp_purecpl
      icpt = icpt + 1
    END IF
    IF (ll_opa) THEN
      nsbc = jp_none
      icpt = icpt + 1
    END IF
    IF (icpt /= 1) CALL ctl_stop('sbc_init : choose ONE and only ONE sbc option')
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      SELECT CASE (nsbc)
      CASE (jp_usr)
        WRITE(numout, FMT = *) '   ==>>>   user defined forcing formulation'
      CASE (jp_flx)
        WRITE(numout, FMT = *) '   ==>>>   flux formulation'
      CASE (jp_blk)
        WRITE(numout, FMT = *) '   ==>>>   bulk formulation'
      CASE (jp_purecpl)
        WRITE(numout, FMT = *) '   ==>>>   pure coupled formulation'
      CASE (jp_none)
        WRITE(numout, FMT = *) '   ==>>>   OPA coupled to SAS via oasis'
        IF (ln_mixcpl) WRITE(numout, FMT = *) '               + forced-coupled mixed formulation'
      END SELECT
      IF (ll_not_nemo) WRITE(numout, FMT = *) '               + OASIS coupled SAS'
    END IF
    IF (lk_oasis) CALL sbc_cpl_init(nn_ice)
    IF (nn_components /= jp_iam_nemo) THEN
      IF (nn_components == jp_iam_opa) nn_fsbc = cpl_freq('O_SFLX') / NINT(rdt)
      IF (nn_components == jp_iam_sas) nn_fsbc = cpl_freq('I_SFLX') / NINT(rdt)
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) "   OPA-SAS coupled via OASIS : nn_fsbc re-defined from OASIS namcouple ", nn_fsbc
        WRITE(numout, FMT = *)
      END IF
    END IF
    IF (MOD(nitend - nit000 + 1, nn_fsbc) /= 0 .OR. MOD(nstock, nn_fsbc) /= 0) THEN
      WRITE(ctmp1, FMT = *) 'sbc_init : experiment length (', nitend - nit000 + 1, ') or nstock (', nstock, ' is NOT a multiple of &
&nn_fsbc (', nn_fsbc, ')'
      CALL ctl_stop(ctmp1, 'Impossible to properly do model restart')
    END IF
    IF (MOD(rday, REAL(nn_fsbc, wp) * rdt) /= 0) CALL ctl_warn('sbc_init : nn_fsbc is NOT a multiple of the number of time steps &
&in a day')
    IF (ln_dm2dc .AND. NINT(rday) / (nn_fsbc * NINT(rdt)) < 8) CALL ctl_warn('sbc_init : diurnal cycle for qsr: the sampling of &
&the diurnal cycle is too small...')
    CALL sbc_ssm_init
    IF (ln_blk) CALL sbc_blk_init
    IF (ln_ssr) CALL sbc_ssr_init
    IF (ln_isf) CALL sbc_isf_init
    CALL sbc_rnf_init
    IF (ln_apr_dyn) CALL sbc_apr_init
    IF (nn_ice == 3) CALL cice_sbc_init(nsbc)
    IF (ln_wave) CALL sbc_wave_init
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('utau_b')
      CALL iom_set_rstw_var_active('vtau_b')
      CALL iom_set_rstw_var_active('qns_b')
      CALL iom_set_rstw_var_active('emp_b')
      CALL iom_set_rstw_var_active('sfx_b')
    END IF
  END SUBROUTINE sbc_init
  SUBROUTINE sbc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    LOGICAL :: ll_sas, ll_opa
    REAL(KIND = wp) :: zthscl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwdht, zwght
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_timing) CALL timing_start('sbc')
    IF (kt /= nit000) THEN
      !$ACC KERNELS
      utau_b(:, :) = utau(:, :)
      vtau_b(:, :) = vtau(:, :)
      qns_b(:, :) = qns(:, :)
      emp_b(:, :) = emp(:, :)
      sfx_b(:, :) = sfx(:, :)
      !$ACC END KERNELS
      IF (ln_rnf) THEN
        !$ACC KERNELS
        rnf_b(:, :) = rnf(:, :)
        rnf_tsc_b(:, :, :) = rnf_tsc(:, :, :)
        !$ACC END KERNELS
      END IF
      IF (ln_isf) THEN
        !$ACC KERNELS
        fwfisf_b(:, :) = fwfisf(:, :)
        risf_tsc_b(:, :, :) = risf_tsc(:, :, :)
        !$ACC END KERNELS
      END IF
    END IF
    CALL profile_psy_data0 % PreStart('sbc', 'r0', 0, 0)
    ll_sas = nn_components == jp_iam_sas
    ll_opa = nn_components == jp_iam_opa
    IF (.NOT. ll_sas) CALL sbc_ssm(kt)
    IF (ln_wave) CALL sbc_wave(kt)
    SELECT CASE (nsbc)
    CASE (jp_usr)
      CALL usrdef_sbc_oce(kt)
    CASE (jp_flx)
      CALL sbc_flx(kt)
    CASE (jp_blk)
      IF (ll_sas) THEN
        IF (ln_timing) CALL timing_start('sbc_sas_rcv')
        CALL sbc_cpl_rcv(kt, nn_fsbc, nn_ice)
        IF (ln_timing) CALL timing_stop('sbc_sas_rcv')
      END IF
      CALL sbc_blk(kt)
    CASE (jp_purecpl)
      CALL sbc_cpl_rcv(kt, nn_fsbc, nn_ice)
    CASE (jp_none)
      IF (ll_opa) THEN
        IF (ln_timing) CALL timing_start('sbc_opa_rcv')
        CALL sbc_cpl_rcv(kt, nn_fsbc, nn_ice)
        IF (ln_timing) CALL timing_stop('sbc_opa_rcv')
      END IF
    END SELECT
    IF (ln_mixcpl) CALL sbc_cpl_rcv(kt, nn_fsbc, nn_ice)
    IF (ln_wave .AND. (ln_tauwoc .OR. ln_tauw)) CALL sbc_wstress
    SELECT CASE (nn_ice)
    CASE (1)
      CALL sbc_ice_if(kt)
    CASE (3)
      CALL sbc_ice_cice(kt, nsbc)
    END SELECT
    IF (ln_icebergs) THEN
      CALL icb_stp(kt)
      IF (.NOT. ln_passive_mode) CALL lbc_lnk('sbcmod', emp, 'T', 1.)
    END IF
    IF (ln_isf) CALL sbc_isf(kt)
    IF (ln_rnf) CALL sbc_rnf(kt)
    IF (ln_ssr) CALL sbc_ssr(kt)
    IF (nn_fwb /= 0) CALL sbc_fwb(kt, nn_fwb, nn_fsbc)
    IF (l_sbc_clo .AND. (.NOT. ln_diurnal_only)) CALL sbc_clo(kt)
    CALL profile_psy_data0 % PostEnd
    IF (ll_wd) THEN
      !$ACC KERNELS
      zthscl = atanh(rn_wd_sbcfra)
      zwdht(:, :) = sshn(:, :) + ht_0(:, :) - rn_wdmin1
      WHERE (zwdht(:, :) <= 0.0)
        taum(:, :) = 0.0
        utau(:, :) = 0.0
        vtau(:, :) = 0.0
        qns(:, :) = 0.0
        qsr(:, :) = 0.0
        emp(:, :) = MIN(emp(:, :), 0.0)
        sfx(:, :) = 0.0
      END WHERE
      zwght(:, :) = TANH(zthscl * zwdht(:, :))
      WHERE (zwdht(:, :) > 0.0 .AND. zwdht(:, :) < rn_wd_sbcdep)
        qsr(:, :) = qsr(:, :) * zwght(:, :)
        qns(:, :) = qns(:, :) * zwght(:, :)
        taum(:, :) = taum(:, :) * zwght(:, :)
        utau(:, :) = utau(:, :) * zwght(:, :)
        vtau(:, :) = vtau(:, :) * zwght(:, :)
        sfx(:, :) = sfx(:, :) * zwght(:, :)
        emp(:, :) = emp(:, :) * zwght(:, :)
      END WHERE
      !$ACC END KERNELS
    END IF
    IF (kt == nit000) THEN
      IF (ln_rstart .AND. iom_varid(numror, 'utau_b', ldstop = .FALSE.) > 0) THEN
        CALL profile_psy_data1 % PreStart('sbc', 'r1', 0, 0)
        IF (lwp) WRITE(numout, FMT = *) '          nit000-1 surface forcing fields red in the restart file'
        CALL iom_get(numror, jpdom_autoglo, 'utau_b', utau_b, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'vtau_b', vtau_b, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'qns_b', qns_b, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'emp_b', emp_b, ldxios = lrxios)
        CALL profile_psy_data1 % PostEnd
        IF (iom_varid(numror, 'sfx_b', ldstop = .FALSE.) > 0) THEN
          CALL iom_get(numror, jpdom_autoglo, 'sfx_b', sfx_b, ldxios = lrxios)
        ELSE
          !$ACC KERNELS
          sfx_b(:, :) = sfx(:, :)
          !$ACC END KERNELS
        END IF
      ELSE
        IF (lwp) WRITE(numout, FMT = *) '          nit000-1 surface forcing fields set to nit000'
        !$ACC KERNELS
        utau_b(:, :) = utau(:, :)
        vtau_b(:, :) = vtau(:, :)
        qns_b(:, :) = qns(:, :)
        emp_b(:, :) = emp(:, :)
        sfx_b(:, :) = sfx(:, :)
        !$ACC END KERNELS
      END IF
    END IF
    CALL profile_psy_data2 % PreStart('sbc', 'r2', 0, 0)
    IF (lrst_oce) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc : ocean surface forcing fields written in ocean restart file ', 'at it= ', kt, ' date= &
&', ndastp
      IF (lwp) WRITE(numout, FMT = *) '~~~~'
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'utau_b', utau, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'vtau_b', vtau, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'qns_b', qns, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'emp_b', emp, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'sfx_b', sfx, ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
    IF (MOD(kt - 1, nn_fsbc) == 0) THEN
      CALL iom_put("empmr", emp - rnf)
      CALL iom_put("empbmr", emp_b - rnf)
      CALL iom_put("saltflx", sfx)
      CALL iom_put("fmmflx", fmmflx)
      CALL iom_put("qt", qns + qsr)
      CALL iom_put("qns", qns)
      CALL iom_put("qsr", qsr)
      IF (nn_ice > 0 .OR. ll_opa) CALL iom_put("ice_cover", fr_i)
      CALL iom_put("taum", taum)
      CALL iom_put("wspd", wndm)
    END IF
    CALL iom_put("utau", utau)
    CALL iom_put("vtau", vtau)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab2d_1 = fr_i, clinfo1 = ' fr_i    - : ', mask1 = tmask)
      CALL prt_ctl(tab2d_1 = (emp - rnf + fwfisf), clinfo1 = ' emp-rnf - : ', mask1 = tmask)
      CALL prt_ctl(tab2d_1 = (sfx - rnf + fwfisf), clinfo1 = ' sfx-rnf - : ', mask1 = tmask)
      CALL prt_ctl(tab2d_1 = qns, clinfo1 = ' qns      - : ', mask1 = tmask)
      CALL prt_ctl(tab2d_1 = qsr, clinfo1 = ' qsr      - : ', mask1 = tmask)
      CALL prt_ctl(tab3d_1 = tmask, clinfo1 = ' tmask    - : ', mask1 = tmask, kdim = jpk)
      CALL prt_ctl(tab3d_1 = tsn(:, :, :, jp_tem), clinfo1 = ' sst      - : ', mask1 = tmask, kdim = 1)
      CALL prt_ctl(tab3d_1 = tsn(:, :, :, jp_sal), clinfo1 = ' sss      - : ', mask1 = tmask, kdim = 1)
      CALL prt_ctl(tab2d_1 = utau, clinfo1 = ' utau     - : ', mask1 = umask, tab2d_2 = vtau, clinfo2 = ' vtau     - : ', &
&mask2 = vmask)
    END IF
    IF (kt == nitend) CALL sbc_final
    IF (ln_timing) CALL timing_stop('sbc')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE sbc
  SUBROUTINE sbc_final
    IF (nn_ice == 3) CALL cice_sbc_final
  END SUBROUTINE sbc_final
END MODULE sbcmod