MODULE icedia
  USE dom_oce
  USE phycst
  USE daymod
  USE sbc_oce, ONLY: sfx, nn_fsbc
  USE ice
  USE icerst
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dia
  PUBLIC :: ice_dia_init
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: vol_loc_ini, sal_loc_ini, tem_loc_ini
  REAL(KIND = wp) :: frc_sal, frc_voltop, frc_volbot, frc_temtop, frc_tembot
  CONTAINS
  INTEGER FUNCTION ice_dia_alloc()
    ALLOCATE(vol_loc_ini(jpi, jpj), sal_loc_ini(jpi, jpj), tem_loc_ini(jpi, jpj), STAT = ice_dia_alloc)
    CALL mpp_sum('icedia', ice_dia_alloc)
    IF (ice_dia_alloc /= 0) CALL ctl_stop('STOP', 'ice_dia_alloc: failed to allocate arrays')
  END FUNCTION ice_dia_alloc
  SUBROUTINE ice_dia(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp) :: zbg_ivol, zbg_item, zbg_area, zbg_isal
    REAL(KIND = wp) :: zbg_svol, zbg_stem
    REAL(KIND = wp) :: z_frc_voltop, z_frc_temtop, z_frc_sal
    REAL(KIND = wp) :: z_frc_volbot, z_frc_tembot
    REAL(KIND = wp) :: zdiff_vol, zdiff_sal, zdiff_tem
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_dia', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('ice_dia')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'icedia: output ice diagnostics (integrated over the domain)'
      WRITE(numout, *) '~~~~~~'
    END IF
    zbg_ivol = glob_sum('icedia', vt_i(:, :) * e1e2t(:, :)) * 1.E-9
    zbg_svol = glob_sum('icedia', vt_s(:, :) * e1e2t(:, :)) * 1.E-9
    zbg_area = glob_sum('icedia', at_i(:, :) * e1e2t(:, :)) * 1.E-6
    zbg_isal = glob_sum('icedia', SUM(sv_i(:, :, :), dim = 3) * e1e2t(:, :)) * 1.E-9
    zbg_item = glob_sum('icedia', et_i * e1e2t(:, :)) * 1.E-20
    zbg_stem = glob_sum('icedia', et_s * e1e2t(:, :)) * 1.E-20
    z_frc_volbot = r1_rau0 * glob_sum('icedia', - (wfx_ice(:, :) + wfx_snw(:, :) + wfx_err_sub(:, :)) * e1e2t(:, :)) * 1.E-9
    z_frc_voltop = r1_rau0 * glob_sum('icedia', - (wfx_sub(:, :) + wfx_spr(:, :)) * e1e2t(:, :)) * 1.E-9
    z_frc_sal = r1_rau0 * glob_sum('icedia', - sfx(:, :) * e1e2t(:, :)) * 1.E-9
    z_frc_tembot = glob_sum('icedia', qt_oce_ai(:, :) * e1e2t(:, :)) * 1.E-20
    z_frc_temtop = glob_sum('icedia', qt_atm_oi(:, :) * e1e2t(:, :)) * 1.E-20
    frc_voltop = frc_voltop + z_frc_voltop * rdt_ice
    frc_volbot = frc_volbot + z_frc_volbot * rdt_ice
    frc_sal = frc_sal + z_frc_sal * rdt_ice
    frc_temtop = frc_temtop + z_frc_temtop * rdt_ice
    frc_tembot = frc_tembot + z_frc_tembot * rdt_ice
    zdiff_vol = r1_rau0 * glob_sum('icedia', (rhoi * vt_i(:, :) + rhos * vt_s(:, :) - vol_loc_ini(:, :)) * e1e2t(:, :)) * 1.E-9
    zdiff_sal = r1_rau0 * glob_sum('icedia', (rhoi * SUM(sv_i(:, :, :), dim = 3) - sal_loc_ini(:, :)) * e1e2t(:, :)) * 1.E-9
    zdiff_tem = glob_sum('icedia', (et_i(:, :) + et_s(:, :) - tem_loc_ini(:, :)) * e1e2t(:, :)) * 1.E-20
    zdiff_vol = zdiff_vol - (frc_voltop + frc_volbot)
    zdiff_sal = zdiff_sal - frc_sal
    zdiff_tem = zdiff_tem - (frc_tembot - frc_temtop)
    IF (iom_use('ibgvolume')) CALL iom_put('ibgvolume', zdiff_vol)
    IF (iom_use('ibgsaltco')) CALL iom_put('ibgsaltco', zdiff_sal)
    IF (iom_use('ibgheatco')) CALL iom_put('ibgheatco', zdiff_tem)
    IF (iom_use('ibgheatfx')) CALL iom_put('ibgheatfx', zdiff_tem / glob_sum('icedia', e1e2t(:, :) * 1.E-20 * kt * rdt))
    IF (iom_use('ibgfrcvoltop')) CALL iom_put('ibgfrcvoltop', frc_voltop)
    IF (iom_use('ibgfrcvolbot')) CALL iom_put('ibgfrcvolbot', frc_volbot)
    IF (iom_use('ibgfrcsal')) CALL iom_put('ibgfrcsal', frc_sal)
    IF (iom_use('ibgfrctemtop')) CALL iom_put('ibgfrctemtop', frc_temtop)
    IF (iom_use('ibgfrctembot')) CALL iom_put('ibgfrctembot', frc_tembot)
    IF (iom_use('ibgfrchfxtop')) CALL iom_put('ibgfrchfxtop', frc_temtop / glob_sum('icedia', e1e2t(:, :)) * 1.E-20 * kt * rdt)
    IF (iom_use('ibgfrchfxbot')) CALL iom_put('ibgfrchfxbot', frc_tembot / glob_sum('icedia', e1e2t(:, :)) * 1.E-20 * kt * rdt)
    IF (iom_use('ibgvol_tot')) CALL iom_put('ibgvol_tot', zbg_ivol)
    IF (iom_use('sbgvol_tot')) CALL iom_put('sbgvol_tot', zbg_svol)
    IF (iom_use('ibgarea_tot')) CALL iom_put('ibgarea_tot', zbg_area)
    IF (iom_use('ibgsalt_tot')) CALL iom_put('ibgsalt_tot', zbg_isal)
    IF (iom_use('ibgheat_tot')) CALL iom_put('ibgheat_tot', zbg_item)
    IF (iom_use('sbgheat_tot')) CALL iom_put('sbgheat_tot', zbg_stem)
    IF (lrst_ice) CALL ice_dia_rst('WRITE', kt_ice)
    IF (ln_timing) CALL timing_stop('ice_dia')
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_dia
  SUBROUTINE ice_dia_init
    INTEGER :: ios, ierror
    NAMELIST /namdia/ ln_icediachk, ln_icediahsb, ln_icectl, iiceprt, jiceprt
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namdia, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdia in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namdia, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdia in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namdia)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_dia_init: ice diagnostics'
      WRITE(numout, *) ' ~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namdia:'
      WRITE(numout, *) '      Diagnose online heat/mass/salt budget      ln_icediachk = ', ln_icediachk
      WRITE(numout, *) '      Output          heat/mass/salt budget      ln_icediahsb = ', ln_icediahsb
      WRITE(numout, *) '      control prints for a given grid point      ln_icectl    = ', ln_icectl
      WRITE(numout, *) '         chosen grid point position         (iiceprt,jiceprt) = (', iiceprt, ',', jiceprt, ')'
    END IF
    IF (ln_icediahsb) THEN
      IF (ice_dia_alloc() /= 0) CALL ctl_stop('STOP', 'ice_dia_init : unable to allocate arrays')
      CALL ice_dia_rst('READ')
    END IF
  END SUBROUTINE ice_dia_init
  SUBROUTINE ice_dia_rst(cdrw, kt)
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER, OPTIONAL, INTENT(IN) :: kt
    INTEGER :: iter
    REAL(KIND = wp) :: ziter
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        CALL iom_get(numrir, 'kt_ice', ziter)
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) 'ice_dia_rst read at time step = ', ziter
        IF (lwp) WRITE(numout, *) '~~~~~~~~~~'
        CALL iom_get(numrir, 'frc_voltop', frc_voltop)
        CALL iom_get(numrir, 'frc_volbot', frc_volbot)
        CALL iom_get(numrir, 'frc_temtop', frc_temtop)
        CALL iom_get(numrir, 'frc_tembot', frc_tembot)
        CALL iom_get(numrir, 'frc_sal', frc_sal)
        CALL iom_get(numrir, jpdom_autoglo, 'vol_loc_ini', vol_loc_ini)
        CALL iom_get(numrir, jpdom_autoglo, 'tem_loc_ini', tem_loc_ini)
        CALL iom_get(numrir, jpdom_autoglo, 'sal_loc_ini', sal_loc_ini)
      ELSE
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) ' ice_dia at initial state '
        IF (lwp) WRITE(numout, *) '~~~~~~~'
        !$ACC KERNELS
        frc_voltop = 0._wp
        frc_volbot = 0._wp
        frc_temtop = 0._wp
        frc_tembot = 0._wp
        frc_sal = 0._wp
        vol_loc_ini(:, :) = rhoi * vt_i(:, :) + rhos * vt_s(:, :)
        tem_loc_ini(:, :) = et_i(:, :) + et_s(:, :)
        !$ACC END KERNELS
        sal_loc_ini(:, :) = rhoi * SUM(sv_i(:, :, :), dim = 3)
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      iter = kt + nn_fsbc - 1
      IF (iter == nitrst) THEN
        IF (lwp) WRITE(numout, *)
        IF (lwp) WRITE(numout, *) 'ice_dia_rst write at time step = ', kt
        IF (lwp) WRITE(numout, *) '~~~~~~~~~~~'
      END IF
      CALL iom_rstput(iter, nitrst, numriw, 'frc_voltop', frc_voltop)
      CALL iom_rstput(iter, nitrst, numriw, 'frc_volbot', frc_volbot)
      CALL iom_rstput(iter, nitrst, numriw, 'frc_temtop', frc_temtop)
      CALL iom_rstput(iter, nitrst, numriw, 'frc_tembot', frc_tembot)
      CALL iom_rstput(iter, nitrst, numriw, 'frc_sal', frc_sal)
      CALL iom_rstput(iter, nitrst, numriw, 'vol_loc_ini', vol_loc_ini)
      CALL iom_rstput(iter, nitrst, numriw, 'tem_loc_ini', tem_loc_ini)
      CALL iom_rstput(iter, nitrst, numriw, 'sal_loc_ini', sal_loc_ini)
    END IF
  END SUBROUTINE ice_dia_rst
END MODULE icedia