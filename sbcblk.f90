MODULE sbcblk
  USE oce
  USE dom_oce
  USE phycst
  USE fldread
  USE sbc_oce
  USE cyclone
  USE sbcdcy
  USE sbcwave, ONLY: cdn_wave
  USE sbc_ice
  USE lib_fortran
  USE ice, ONLY: u_ice, v_ice, jpl, a_i_b, at_i_b, t_su, rn_cnd_s, hfx_err_dif
  USE icethd_dh
  USE sbcblk_algo_ncar
  USE sbcblk_algo_coare
  USE sbcblk_algo_coare3p5
  USE sbcblk_algo_ecmwf
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_blk_init
  PUBLIC :: sbc_blk
  PUBLIC :: blk_ice_tau
  PUBLIC :: blk_ice_flx
  PUBLIC :: blk_ice_qcn
  REAL(KIND = wp), PARAMETER :: Cp_dry = 1005.0
  REAL(KIND = wp), PARAMETER :: Cp_vap = 1860.0
  REAL(KIND = wp), PARAMETER :: R_dry = 287.05_wp
  REAL(KIND = wp), PARAMETER :: R_vap = 461.495_wp
  REAL(KIND = wp), PARAMETER :: reps0 = R_dry / R_vap
  REAL(KIND = wp), PARAMETER :: rctv0 = R_vap / R_dry
  INTEGER, PARAMETER :: jpfld = 10
  INTEGER, PARAMETER :: jp_wndi = 1
  INTEGER, PARAMETER :: jp_wndj = 2
  INTEGER, PARAMETER :: jp_tair = 3
  INTEGER, PARAMETER :: jp_humi = 4
  INTEGER, PARAMETER :: jp_qsr = 5
  INTEGER, PARAMETER :: jp_qlw = 6
  INTEGER, PARAMETER :: jp_prec = 7
  INTEGER, PARAMETER :: jp_snow = 8
  INTEGER, PARAMETER :: jp_slp = 9
  INTEGER, PARAMETER :: jp_tdif = 10
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf
  REAL(KIND = wp), PARAMETER :: cpa = 1000.5
  REAL(KIND = wp), PARAMETER :: Ls = 2.839E6
  REAL(KIND = wp), PARAMETER :: Stef = 5.67E-8
  REAL(KIND = wp), PARAMETER :: Cd_ice = 1.4E-3
  REAL(KIND = wp), PARAMETER :: albo = 0.066
  LOGICAL :: ln_NCAR
  LOGICAL :: ln_COARE_3p0
  LOGICAL :: ln_COARE_3p5
  LOGICAL :: ln_ECMWF
  LOGICAL :: ln_taudif
  REAL(KIND = wp) :: rn_pfac
  REAL(KIND = wp) :: rn_efac
  REAL(KIND = wp) :: rn_vfac
  REAL(KIND = wp) :: rn_zqt
  REAL(KIND = wp) :: rn_zu
  LOGICAL :: ln_Cd_L12 = .FALSE.
  LOGICAL :: ln_Cd_L15 = .FALSE.
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: Cd_atm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: Ch_atm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: Ce_atm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: t_zu
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: q_zu
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: cdn_oce, chn_oce, cen_oce
  INTEGER :: nblk
  INTEGER, PARAMETER :: np_NCAR = 1
  INTEGER, PARAMETER :: np_COARE_3p0 = 2
  INTEGER, PARAMETER :: np_COARE_3p5 = 3
  INTEGER, PARAMETER :: np_ECMWF = 4
  CONTAINS
  INTEGER FUNCTION sbc_blk_alloc()
    ALLOCATE(Cd_atm(jpi, jpj), Ch_atm(jpi, jpj), Ce_atm(jpi, jpj), t_zu(jpi, jpj), q_zu(jpi, jpj), cdn_oce(jpi, jpj), chn_oce(jpi, &
&jpj), cen_oce(jpi, jpj), STAT = sbc_blk_alloc)
    CALL mpp_sum('sbcblk', sbc_blk_alloc)
    IF (sbc_blk_alloc /= 0) CALL ctl_stop('STOP', 'sbc_blk_alloc: failed to allocate arrays')
  END FUNCTION sbc_blk_alloc
  SUBROUTINE sbc_blk_init
    INTEGER :: ifpr, jfld
    INTEGER :: ios, ierror, ioptio
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N), DIMENSION(jpfld) :: slf_i
    TYPE(FLD_N) :: sn_wndi, sn_wndj, sn_humi, sn_qsr
    TYPE(FLD_N) :: sn_qlw, sn_tair, sn_prec, sn_snow
    TYPE(FLD_N) :: sn_slp, sn_tdif
    NAMELIST /namsbc_blk/ sn_wndi, sn_wndj, sn_humi, sn_qsr, sn_qlw, sn_tair, sn_prec, sn_snow, sn_slp, sn_tdif, ln_NCAR, &
&ln_COARE_3p0, ln_COARE_3p5, ln_ECMWF, cn_dir, ln_taudif, rn_zqt, rn_zu, rn_pfac, rn_efac, rn_vfac, ln_Cd_L12, ln_Cd_L15
    IF (sbc_blk_alloc() /= 0) CALL ctl_stop('STOP', 'sbc_blk : unable to allocate standard arrays')
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_blk, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_blk in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_blk, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_blk in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_blk)
    ioptio = 0
    IF (ln_ncar) THEN
      nblk = np_ncar
      ioptio = ioptio + 1
    END IF
    IF (ln_coare_3p0) THEN
      nblk = np_coare_3p0
      ioptio = ioptio + 1
    END IF
    IF (ln_coare_3p5) THEN
      nblk = np_coare_3p5
      ioptio = ioptio + 1
    END IF
    IF (ln_ecmwf) THEN
      nblk = np_ecmwf
      ioptio = ioptio + 1
    END IF
    IF (ioptio /= 1) CALL ctl_stop('sbc_blk_init: Choose one and only one bulk algorithm')
    IF (ln_dm2dc) THEN
      IF (sn_qsr % nfreqh /= 24) CALL ctl_stop('sbc_blk_init: ln_dm2dc=T only with daily short-wave input')
      IF (sn_qsr % ln_tint) THEN
        CALL ctl_warn('sbc_blk_init: ln_dm2dc=T daily qsr time interpolation done by sbcdcy module', &
&'              ==> We force time interpolation = .false. for qsr')
        sn_qsr % ln_tint = .FALSE.
      END IF
    END IF
    slf_i(jp_wndi) = sn_wndi
    slf_i(jp_wndj) = sn_wndj
    slf_i(jp_qsr) = sn_qsr
    slf_i(jp_qlw) = sn_qlw
    slf_i(jp_tair) = sn_tair
    slf_i(jp_humi) = sn_humi
    slf_i(jp_prec) = sn_prec
    slf_i(jp_snow) = sn_snow
    slf_i(jp_slp) = sn_slp
    slf_i(jp_tdif) = sn_tdif
    lhftau = ln_taudif
    jfld = jpfld - COUNT((/.NOT. lhftau/))
    ALLOCATE(sf(jfld), STAT = ierror)
    IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_blk_init: unable to allocate sf structure')
    DO ifpr = 1, jfld
      ALLOCATE(sf(ifpr) % fnow(jpi, jpj, 1))
      IF (slf_i(ifpr) % ln_tint) ALLOCATE(sf(ifpr) % fdta(jpi, jpj, 1, 2))
      IF (slf_i(ifpr) % nfreqh > 0. .AND. MOD(3600. * slf_i(ifpr) % nfreqh, REAL(nn_fsbc) * rdt) /= 0.) CALL &
&ctl_warn('sbc_blk_init: sbcmod timestep rdt*nn_fsbc is NOT a submultiple of atmospheric forcing frequency.', '               This &
&is not ideal. You should consider changing either rdt or nn_fsbc value...')
    END DO
    CALL fld_fill(sf, slf_i, cn_dir, 'sbc_blk_init', 'surface boundary condition -- bulk formulae', 'namsbc_blk')
    IF (ln_wave) THEN
      IF (.NOT. (ln_cdgw .OR. ln_sdw .OR. ln_tauwoc .OR. ln_stcor)) THEN
        CALL ctl_stop('STOP', 'Ask for wave coupling but ln_cdgw=F, ln_sdw=F, ln_tauwoc=F, ln_stcor=F')
      ELSE IF (ln_cdgw .AND. .NOT. ln_NCAR) THEN
        CALL ctl_stop('drag coefficient read from wave model definable only with NCAR and CORE bulk formulae')
      ELSE IF (ln_stcor .AND. .NOT. ln_sdw) THEN
        CALL ctl_stop('Stokes-Coriolis term calculated only if activated Stokes Drift ln_sdw=T')
      END IF
    ELSE
      IF (ln_cdgw .OR. ln_sdw .OR. ln_tauwoc .OR. ln_stcor) CALL ctl_stop('Not Activated Wave Module (ln_wave=F) but asked &
&coupling ', 'with drag coefficient (ln_cdgw =T) ', 'or Stokes Drift (ln_sdw=T) ', 'or ocean stress modification due to waves &
&(ln_tauwoc=T) ', 'or Stokes-Coriolis term (ln_stcori=T)')
    END IF
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Namelist namsbc_blk (other than data information):'
      WRITE(numout, FMT = *) '      "NCAR"      algorithm   (Large and Yeager 2008)     ln_NCAR      = ', ln_NCAR
      WRITE(numout, FMT = *) '      "COARE 3.0" algorithm   (Fairall et al. 2003)       ln_COARE_3p0 = ', ln_COARE_3p0
      WRITE(numout, FMT = *) '      "COARE 3.5" algorithm   (Edson et al. 2013)         ln_COARE_3p5 = ', ln_COARE_3p0
      WRITE(numout, FMT = *) '      "ECMWF"     algorithm   (IFS cycle 31)              ln_ECMWF     = ', ln_ECMWF
      WRITE(numout, FMT = *) '      add High freq.contribution to the stress module     ln_taudif    = ', ln_taudif
      WRITE(numout, FMT = *) '      Air temperature and humidity reference height (m)   rn_zqt       = ', rn_zqt
      WRITE(numout, FMT = *) '      Wind vector reference height (m)                    rn_zu        = ', rn_zu
      WRITE(numout, FMT = *) '      factor applied on precipitation (total & snow)      rn_pfac      = ', rn_pfac
      WRITE(numout, FMT = *) '      factor applied on evaporation                       rn_efac      = ', rn_efac
      WRITE(numout, FMT = *) '      factor applied on ocean/ice velocity                rn_vfac      = ', rn_vfac
      WRITE(numout, FMT = *) '         (form absolute (=0) to relative winds(=1))'
      WRITE(numout, FMT = *) '      use ice-atm drag from Lupkes2012                    ln_Cd_L12    = ', ln_Cd_L12
      WRITE(numout, FMT = *) '      use ice-atm drag from Lupkes2015                    ln_Cd_L15    = ', ln_Cd_L15
      WRITE(numout, FMT = *)
      SELECT CASE (nblk)
      CASE (np_ncar)
        WRITE(numout, FMT = *) '   ==>>>   "NCAR" algorithm        (Large and Yeager 2008)'
      CASE (np_coare_3p0)
        WRITE(numout, FMT = *) '   ==>>>   "COARE 3.0" algorithm   (Fairall et al. 2003)'
      CASE (np_coare_3p5)
        WRITE(numout, FMT = *) '   ==>>>   "COARE 3.5" algorithm   (Edson et al. 2013)'
      CASE (np_ecmwf)
        WRITE(numout, FMT = *) '   ==>>>   "ECMWF" algorithm       (IFS cycle 31)'
      END SELECT
    END IF
  END SUBROUTINE sbc_blk_init
  SUBROUTINE sbc_blk(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_blk', 'r0', 0, 0)
    CALL fld_read(kt, nn_fsbc, sf)
    IF (MOD(kt - 1, nn_fsbc) == 0) CALL blk_oce(kt, sf, sst_m, ssu_m, ssv_m)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sbc_blk
  SUBROUTINE blk_oce(kt, sf, pst, pu, pv)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    TYPE(fld), INTENT(INOUT), DIMENSION(:) :: sf
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: pst
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: pu
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: pv
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zztmp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwnd_i, zwnd_j
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsq
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zqlw, zqsb
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zqla, zevap
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zst
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zU_zu
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztpot
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zrhoa
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    !$ACC KERNELS
    zst(:, :) = pst(:, :) + rt0
    zwnd_i(:, :) = 0._wp
    zwnd_j(:, :) = 0._wp
    ! CALL profile_psy_data0 % PreStart('blk_oce', 'r0', 0, 0)
    !$ACC LOOP COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zwnd_i(ji, jj) = (sf(jp_wndi) % fnow(ji, jj, 1) - rn_vfac * 0.5 * (pu(ji - 1, jj) + pu(ji, jj)))
        zwnd_j(ji, jj) = (sf(jp_wndj) % fnow(ji, jj, 1) - rn_vfac * 0.5 * (pv(ji, jj - 1) + pv(ji, jj)))
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('blk_oce', 'r0', 0, 0)
    CALL lbc_lnk_multi('sbcblk', zwnd_i, 'T', - 1., zwnd_j, 'T', - 1.)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    wndm(:, :) = SQRT(zwnd_i(:, :) * zwnd_i(:, :) + zwnd_j(:, :) * zwnd_j(:, :)) * tmask(:, :, 1)
    zztmp = 1. - albo
    !$ACC END KERNELS
    !CALL profile_psy_data1 % PreStart('blk_oce', 'r1', 0, 0)
    IF (ln_dm2dc) THEN
      qsr(:, :) = sbc_dcy(sf(jp_qsr) % fnow(:, :, 1))
      !$ACC KERNELS ! CDe added
      qsr(:, :) = zztmp * qsr(:,:) * tmask(:, :, 1)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS ! CDe added      
      qsr(:, :) = zztmp * sf(jp_qsr) % fnow(:, :, 1) * tmask(:, :, 1)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS ! CDe added
    zqlw(:, :) = (sf(jp_qlw) % fnow(:, :, 1) - Stef * zst(:, :) * zst(:, :) * zst(:, :) * zst(:, :)) * tmask(:, :, 1)
    !$ACC END KERNELS
    zsq(:, :) = q_sat(zst(:, :), sf(jp_slp) % fnow(:, :, 1))
    !$ACC KERNELS ! CDe added
    zsq(:, :) = 0.98 * zsq(:,:)
    !$ACC END KERNELS
    ztpot(:,:) = gamma_moist(sf(jp_tair) % fnow(:, :, 1), sf(jp_humi) % fnow(:, :, 1))
    !$ACC KERNELS ! CDe added
    ztpot(:,:) = sf(jp_tair) % fnow(:, :, 1) + ztpot(:,:) * rn_zqt
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('blk_oce', 'r1', 0, 0)
    SELECT CASE (nblk)
    CASE (np_ncar)
      CALL turb_ncar(rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi) % fnow, wndm, cd_atm, ch_atm, ce_atm, t_zu, q_zu, zu_zu, cdn_oce, &
&chn_oce, cen_oce)
    CASE (np_coare_3p0)
      CALL turb_coare(rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi) % fnow, wndm, cd_atm, ch_atm, ce_atm, t_zu, q_zu, zu_zu, &
&cdn_oce, chn_oce, cen_oce)
    CASE (np_coare_3p5)
      CALL turb_coare3p5(rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi) % fnow, wndm, cd_atm, ch_atm, ce_atm, t_zu, q_zu, zu_zu, &
&cdn_oce, chn_oce, cen_oce)
    CASE (np_ecmwf)
      CALL turb_ecmwf(rn_zqt, rn_zu, zst, ztpot, zsq, sf(jp_humi) % fnow, wndm, cd_atm, ch_atm, ce_atm, t_zu, q_zu, zu_zu, &
&cdn_oce, chn_oce, cen_oce)
    CASE DEFAULT
      CALL ctl_stop('STOP', 'sbc_oce: non-existing bulk formula selected')
    END SELECT
    IF (ABS(rn_zu - rn_zqt) > 0.01) THEN
      zrhoa(:, :) = rho_air(t_zu(:, :), q_zu(:, :), sf(jp_slp) % fnow(:, :, 1))
    ELSE
      zrhoa(:, :) = rho_air(sf(jp_tair) % fnow(:, :, 1), sf(jp_humi) % fnow(:, :, 1), sf(jp_slp) % fnow(:, :, 1))
    END IF
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zztmp = zrhoa(ji, jj) * zU_zu(ji, jj) * Cd_atm(ji, jj)
        taum(ji, jj) = zztmp * wndm(ji, jj)
        zwnd_i(ji, jj) = zztmp * zwnd_i(ji, jj)
        zwnd_j(ji, jj) = zztmp * zwnd_j(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('blk_oce', 'r2', 0, 0)
    IF (lhftau) taum(:, :) = taum(:, :) + sf(jp_tdif) % fnow(:, :, 1)
    CALL iom_put("taum_oce", taum)
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        utau(ji, jj) = 0.5 * (2. - umask(ji, jj, 1)) * (zwnd_i(ji, jj) + zwnd_i(ji + 1, jj)) * MAX(tmask(ji, jj, 1), tmask(ji + 1, &
&jj, 1))
        vtau(ji, jj) = 0.5 * (2. - vmask(ji, jj, 1)) * (zwnd_j(ji, jj) + zwnd_j(ji, jj + 1)) * MAX(tmask(ji, jj, 1), tmask(ji, jj &
&+ 1, 1))
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('sbcblk', utau, 'U', - 1., vtau, 'V', - 1.)
    !$ACC KERNELS
    zqla(:, :) = zrhoa(:, :) * zU_zu(:, :) * tmask(:, :, 1)
    !$ACC END KERNELS
    IF (ABS(rn_zu - rn_zqt) < 0.01_wp) THEN
      !$ACC KERNELS ! CDe added
      zevap(:, :) = rn_efac * MAX(0._wp, zqla(:, :) * Ce_atm(:, :) * (zsq(:, :) - sf(jp_humi) % fnow(:, :, 1)))
      !CALL profile_psy_data3 % PreStart('blk_oce', 'r3', 0, 0)
      !$ACC END KERNELS
      zqsb(:, :) = cp_air(sf(jp_humi) % fnow(:, :, 1))
      !$ACC KERNELS ! CDe added
      zqsb(:, :) = zqsb(:,:) * zqla(:, :) * Ch_atm(:, :) * (zst(:, :) - ztpot(:, :))
      !$ACC END KERNELS
      !CALL profile_psy_data3 % PostEnd
    ELSE
      !$ACC KERNELS
      zevap(:, :) = rn_efac * MAX(0._wp, zqla(:, :) * Ce_atm(:, :) * (zsq(:, :) - q_zu(:, :)))
      !$ACC END KERNELS
      !CALL profile_psy_data4 % PreStart('blk_oce', 'r4', 0, 0)
      zqsb(:, :) = cp_air(sf(jp_humi) % fnow(:, :, 1))
      !$ACC KERNELS ! CDe added
      zqsb(:, :) = zqsb(:,:) * zqla(:, :) * Ch_atm(:, :) * (zst(:, :) - t_zu(:, :))
      !$ACC END KERNELS
      !CALL profile_psy_data4 % PostEnd
    END IF
    !CALL profile_psy_data5 % PreStart('blk_oce', 'r5', 0, 0)
    zqla(:, :) = L_vap(zst(:, :))
    !$ACC KERNELS ! CDe added
    zqla(:, :) = zqla(:,:) * zevap(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data5 % PreStart('blk_oce', 'r5', 0, 0)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab2d_1 = zqla, clinfo1 = ' blk_oce: zqla   : ', tab2d_2 = Ce_atm, clinfo2 = ' Ce_oce  : ')
      CALL prt_ctl(tab2d_1 = zqsb, clinfo1 = ' blk_oce: zqsb   : ', tab2d_2 = Ch_atm, clinfo2 = ' Ch_oce  : ')
      CALL prt_ctl(tab2d_1 = zqlw, clinfo1 = ' blk_oce: zqlw   : ', tab2d_2 = qsr, clinfo2 = ' qsr : ')
      CALL prt_ctl(tab2d_1 = zsq, clinfo1 = ' blk_oce: zsq    : ', tab2d_2 = zst, clinfo2 = ' zst : ')
      CALL prt_ctl(tab2d_1 = utau, clinfo1 = ' blk_oce: utau   : ', mask1 = umask, tab2d_2 = vtau, clinfo2 = ' vtau : ', &
&mask2 = vmask)
      CALL prt_ctl(tab2d_1 = wndm, clinfo1 = ' blk_oce: wndm   : ')
      CALL prt_ctl(tab2d_1 = zst, clinfo1 = ' blk_oce: zst    : ')
    END IF
    CALL profile_psy_data5 % PostEnd
    !$ACC KERNELS ! CDe added
    emp(:, :) = (zevap(:, :) - sf(jp_prec) % fnow(:, :, 1) * rn_pfac) * tmask(:, :, 1)
    qns(:, :) = zqlw(:, :) - zqsb(:, :) - zqla(:, :) - sf(jp_snow) % fnow(:, :, 1) * rn_pfac * rLfus - zevap(:, :) * pst(:, :) * &
&rcp + (sf(jp_prec) % fnow(:, :, 1) - sf(jp_snow) % fnow(:, :, 1)) * rn_pfac * (sf(jp_tair) % fnow(:, :, 1) - rt0) * rcp + &
&sf(jp_snow) % fnow(:, :, 1) * rn_pfac * (MIN(sf(jp_tair) % fnow(:, :, 1), rt0) - rt0) * rcpi
    qns(:, :) = qns(:, :) * tmask(:, :, 1)
    qns_oce(:, :) = zqlw(:, :) - zqsb(:, :) - zqla(:, :)
    qsr_oce(:, :) = qsr(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data6 % PreStart('blk_oce', 'r6', 0, 0)
    IF (nn_ice == 0) THEN
      CALL iom_put("qlw_oce", zqlw)
      CALL iom_put("qsb_oce", - zqsb)
      CALL iom_put("qla_oce", - zqla)
      CALL iom_put("qemp_oce", qns - zqlw + zqsb + zqla)
      CALL iom_put("qns_oce", qns)
      CALL iom_put("qsr_oce", qsr)
      CALL iom_put("qt_oce", qns + qsr)
      tprecip(:, :) = sf(jp_prec) % fnow(:, :, 1) * rn_pfac * tmask(:, :, 1)
      sprecip(:, :) = sf(jp_snow) % fnow(:, :, 1) * rn_pfac * tmask(:, :, 1)
      CALL iom_put('snowpre', sprecip)
      CALL iom_put('precip', tprecip)
    END IF
    IF (ln_ctl) THEN
      CALL prt_ctl(tab2d_1 = zqsb, clinfo1 = ' blk_oce: zqsb   : ', tab2d_2 = zqlw, clinfo2 = ' zqlw  : ')
      CALL prt_ctl(tab2d_1 = zqla, clinfo1 = ' blk_oce: zqla   : ', tab2d_2 = qsr, clinfo2 = ' qsr   : ')
      CALL prt_ctl(tab2d_1 = pst, clinfo1 = ' blk_oce: pst    : ', tab2d_2 = emp, clinfo2 = ' emp   : ')
      CALL prt_ctl(tab2d_1 = utau, clinfo1 = ' blk_oce: utau   : ', mask1 = umask, tab2d_2 = vtau, clinfo2 = ' vtau  : ', &
&mask2 = vmask)
    END IF
    CALL profile_psy_data6 % PostEnd
  END SUBROUTINE blk_oce
  FUNCTION rho_air(ptak, pqa, pslp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptak
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqa
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pslp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: rho_air
    !$ACC KERNELS ! CDe added
    rho_air(:,:) = pslp(:,:) / (R_dry * ptak(:,:) * (1._wp + rctv0 * pqa(:,:)))
    !$ACC END KERNELS
  END FUNCTION rho_air
  FUNCTION cp_air(pqa)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqa
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: cp_air
    !$ACC KERNELS ! CDe added
    Cp_air(:,:) = Cp_dry + Cp_vap * pqa(:,:)
    !$ACC END KERNELS
  END FUNCTION cp_air
  FUNCTION q_sat(ptak, pslp)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptak
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pslp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: q_sat
    INTEGER :: ji, jj
    REAL(KIND = wp) :: ze_sat, ztmp
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        ztmp = rt0 / ptak(ji, jj)
        ze_sat = 10. ** (10.79574 * (1. - ztmp) - 5.028 * LOG10(ptak(ji, jj) / rt0) + 1.50475 * 10. ** (- 4) * (1. - 10. ** (- &
&8.2969 * (ptak(ji, jj) / rt0 - 1.))) + 0.42873 * 10. ** (- 3) * (10. ** (4.76955 * (1. - ztmp)) - 1.) + 0.78614)
        q_sat(ji, jj) = reps0 * ze_sat / (0.01_wp * pslp(ji, jj) - (1._wp - reps0) * ze_sat)
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION q_sat
  FUNCTION gamma_moist(ptak, pqa)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptak
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqa
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: gamma_moist
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zrv, ziRT
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zrv = pqa(ji, jj) / (1. - pqa(ji, jj))
        ziRT = 1. / (R_dry * ptak(ji, jj))
        gamma_moist(ji, jj) = grav * (1. + rLevap * zrv * ziRT) / (Cp_dry + rLevap * rLevap * zrv * reps0 * ziRT / ptak(ji, jj))
      END DO
    END DO
    !$ACC END KERNELS
  END FUNCTION gamma_moist
  FUNCTION L_vap(psst)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: L_vap
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: psst
    !$ACC KERNELS ! CDe added
    L_vap(:,:) = (2.501 - 0.00237 * (psst(:, :) - rt0)) * 1.E6
    !$ACC END KERNELS
  END FUNCTION L_vap
  SUBROUTINE blk_ice_tau
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zwndi_f, zwndj_f, zwnorm_f
    REAL(KIND = wp) :: zwndi_t, zwndj_t
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zrhoa
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    !$ACC KERNELS
    Cd_atm(:, :) = Cd_ice
    Ch_atm(:, :) = Cd_ice
    Ce_atm(:, :) = Cd_ice
    wndm_ice(:, :) = 0._wp
    !$ACC LOOP COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zwndi_t = (sf(jp_wndi) % fnow(ji, jj, 1) - rn_vfac * 0.5 * (u_ice(ji - 1, jj) + u_ice(ji, jj)))
        zwndj_t = (sf(jp_wndj) % fnow(ji, jj, 1) - rn_vfac * 0.5 * (v_ice(ji, jj - 1) + v_ice(ji, jj)))
        wndm_ice(ji, jj) = SQRT(zwndi_t * zwndi_t + zwndj_t * zwndj_t) * tmask(ji, jj, 1)
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('blk_ice_tau', 'r0', 0, 0)
    CALL lbc_lnk('sbcblk', wndm_ice, 'T', 1.)
    CALL profile_psy_data0 % PostEnd
    IF (ln_Cd_L12) THEN
      CALL Cdn10_Lupkes2012(Cd_atm)
      !$ACC KERNELS
      Ch_atm(:, :) = Cd_atm(:, :)
      !$ACC END KERNELS
    ELSE IF (ln_Cd_L15) THEN
      CALL profile_psy_data1 % PreStart('blk_ice_tau', 'r1', 0, 0)
      CALL Cdn10_Lupkes2015(Cd_atm, Ch_atm)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('blk_ice_tau', 'r2', 0, 0)
    zrhoa(:, :) = rho_air(sf(jp_tair) % fnow(:, :, 1), sf(jp_humi) % fnow(:, :, 1), sf(jp_slp) % fnow(:, :, 1))
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    utau_ice(:, :) = 0._wp
    vtau_ice(:, :) = 0._wp
    !$ACC LOOP COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        utau_ice(ji, jj) = 0.5 * zrhoa(ji, jj) * Cd_atm(ji, jj) * (wndm_ice(ji + 1, jj) + wndm_ice(ji, jj)) * (0.5 * (sf(jp_wndi) &
&% fnow(ji + 1, jj, 1) + sf(jp_wndi) % fnow(ji, jj, 1)) - rn_vfac * u_ice(ji, jj))
        vtau_ice(ji, jj) = 0.5 * zrhoa(ji, jj) * Cd_atm(ji, jj) * (wndm_ice(ji, jj + 1) + wndm_ice(ji, jj)) * (0.5 * (sf(jp_wndj) &
&% fnow(ji, jj + 1, 1) + sf(jp_wndj) % fnow(ji, jj, 1)) - rn_vfac * v_ice(ji, jj))
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('blk_ice_tau', 'r3', 0, 0)
    CALL lbc_lnk_multi('sbcblk', utau_ice, 'U', - 1., vtau_ice, 'V', - 1.)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab2d_1 = utau_ice, clinfo1 = ' blk_ice: utau_ice : ', tab2d_2 = vtau_ice, clinfo2 = ' vtau_ice : ')
      CALL prt_ctl(tab2d_1 = wndm_ice, clinfo1 = ' blk_ice: wndm_ice : ')
    END IF
    CALL profile_psy_data3 % PostEnd
  END SUBROUTINE blk_ice_tau
  SUBROUTINE blk_ice_flx(ptsu, phs, phi, palb)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptsu
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phs
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phi
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: palb
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: zst3
    REAL(KIND = wp) :: zcoef_dqlw, zcoef_dqla
    REAL(KIND = wp) :: zztmp, z1_rLsub
    REAL(KIND = wp) :: zfr1, zfr2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z1_st
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z_qlw
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z_qsb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z_dqlw
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: z_dqsb
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zevap, zsnw
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zrhoa
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('blk_ice_flx', 'r0', 0, 0)
    zcoef_dqlw = 4.0 * 0.95 * Stef
    zcoef_dqla = - Ls * 11637800. * (- 5897.8)
    zrhoa(:, :) = rho_air(sf(jp_tair) % fnow(:, :, 1), sf(jp_humi) % fnow(:, :, 1), sf(jp_slp) % fnow(:, :, 1))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zztmp = 1. / (1. - albo)
    WHERE (ptsu(:, :, :) /= 0._wp)
      z1_st(:, :, :) = 1._wp / ptsu(:, :, :)
    ELSEWHERE
      z1_st(:, :, :) = 0._wp
    END WHERE
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('blk_ice_flx', 'r1', 0, 0)
    DO jl = 1, jpl
      DO jj = 1, jpj
        DO ji = 1, jpi
          zst3 = ptsu(ji, jj, jl) * ptsu(ji, jj, jl) * ptsu(ji, jj, jl)
          qsr_ice(ji, jj, jl) = zztmp * (1. - palb(ji, jj, jl)) * qsr(ji, jj)
          z_qlw(ji, jj, jl) = 0.95 * (sf(jp_qlw) % fnow(ji, jj, 1) - Stef * ptsu(ji, jj, jl) * zst3) * tmask(ji, jj, 1)
          z_dqlw(ji, jj, jl) = zcoef_dqlw * zst3
          z_qsb(ji, jj, jl) = zrhoa(ji, jj) * cpa * Ch_atm(ji, jj) * wndm_ice(ji, jj) * (ptsu(ji, jj, jl) - sf(jp_tair) % fnow(ji, &
&jj, 1))
          qla_ice(ji, jj, jl) = rn_efac * MAX(0.E0, zrhoa(ji, jj) * Ls * Ch_atm(ji, jj) * wndm_ice(ji, jj) * (11637800. * EXP(- &
&5897.8 * z1_st(ji, jj, jl)) / zrhoa(ji, jj) - sf(jp_humi) % fnow(ji, jj, 1)))
          IF (qla_ice(ji, jj, jl) > 0._wp) THEN
            dqla_ice(ji, jj, jl) = rn_efac * zcoef_dqla * Ch_atm(ji, jj) * wndm_ice(ji, jj) * z1_st(ji, jj, jl) * z1_st(ji, jj, &
&jl) * EXP(- 5897.8 * z1_st(ji, jj, jl))
          ELSE
            dqla_ice(ji, jj, jl) = 0._wp
          END IF
          z_dqsb(ji, jj, jl) = zrhoa(ji, jj) * cpa * Ch_atm(ji, jj) * wndm_ice(ji, jj)
          qns_ice(ji, jj, jl) = z_qlw(ji, jj, jl) - z_qsb(ji, jj, jl) - qla_ice(ji, jj, jl)
          dqns_ice(ji, jj, jl) = - (z_dqlw(ji, jj, jl) + z_dqsb(ji, jj, jl) + dqla_ice(ji, jj, jl))
        END DO
      END DO
    END DO
    tprecip(:, :) = sf(jp_prec) % fnow(:, :, 1) * rn_pfac * tmask(:, :, 1)
    sprecip(:, :) = sf(jp_snow) % fnow(:, :, 1) * rn_pfac * tmask(:, :, 1)
    CALL iom_put('snowpre', sprecip)
    CALL iom_put('precip', tprecip)
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    z1_rLsub = 1._wp / rLsub
    evap_ice(:, :, :) = rn_efac * qla_ice(:, :, :) * z1_rLsub
    devap_ice(:, :, :) = rn_efac * dqla_ice(:, :, :) * z1_rLsub
    zevap(:, :) = rn_efac * (emp(:, :) + tprecip(:, :))
    zsnw(:, :) = 0._wp
    !$ACC END KERNELS
    CALL ice_thd_snwblow((1. - at_i_b(:, :)), zsnw)
    !$ACC KERNELS
    emp_oce(:, :) = (1._wp - at_i_b(:, :)) * zevap(:, :) - (tprecip(:, :) - sprecip(:, :)) - sprecip(:, :) * (1._wp - zsnw)
    emp_ice(:, :) = SUM(a_i_b(:, :, :) * evap_ice(:, :, :), dim = 3) - sprecip(:, :) * zsnw
    emp_tot(:, :) = emp_oce(:, :) + emp_ice(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('blk_ice_flx', 'r2', 0, 0)
    qemp_oce(:, :) = - (1._wp - at_i_b(:, :)) * zevap(:, :) * sst_m(:, :) * rcp + (tprecip(:, :) - sprecip(:, :)) * (sf(jp_tair) % &
&fnow(:, :, 1) - rt0) * rcp + sprecip(:, :) * (1._wp - zsnw) * ((MIN(sf(jp_tair) % fnow(:, :, 1), rt0) - rt0) * rcpi * tmask(:, :, &
&1) - rLfus)
    qemp_ice(:, :) = sprecip(:, :) * zsnw * ((MIN(sf(jp_tair) % fnow(:, :, 1), rt0) - rt0) * rcpi * tmask(:, :, 1) - rLfus)
    CALL profile_psy_data2 % PostEnd
    !$ACC KERNELS
    qns_tot(:, :) = (1._wp - at_i_b(:, :)) * qns_oce(:, :) + SUM(a_i_b(:, :, :) * qns_ice(:, :, :), dim = 3) + qemp_ice(:, :) + &
&qemp_oce(:, :)
    qsr_tot(:, :) = (1._wp - at_i_b(:, :)) * qsr_oce(:, :) + SUM(a_i_b(:, :, :) * qsr_ice(:, :, :), dim = 3)
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('blk_ice_flx', 'r3', 0, 0)
    qprec_ice(:, :) = rhos * ((MIN(sf(jp_tair) % fnow(:, :, 1), rt0) - rt0) * rcpi * tmask(:, :, 1) - rLfus)
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    DO jl = 1, jpl
      qevap_ice(:, :, jl) = 0._wp
    END DO
    zfr1 = (0.18 * (1.0 - cldf_ice) + 0.35 * cldf_ice)
    zfr2 = (0.82 * (1.0 - cldf_ice) + 0.65 * cldf_ice)
    !$ACC END KERNELS
    CALL profile_psy_data4 % PreStart('blk_ice_flx', 'r4', 0, 0)
    WHERE (phs(:, :, :) <= 0._wp .AND. phi(:, :, :) < 0.1_wp)
      qtr_ice_top(:, :, :) = qsr_ice(:, :, :) * (zfr1 + zfr2 * (1._wp - phi(:, :, :) * 10._wp))
    ELSEWHERE(phs(:, :, :) <= 0._wp .AND. phi(:, :, :) >= 0.1_wp)
      qtr_ice_top(:, :, :) = qsr_ice(:, :, :) * zfr1
    ELSEWHERE
      qtr_ice_top(:, :, :) = 0._wp
    END WHERE
    IF (ln_ctl) THEN
      CALL prt_ctl(tab3d_1 = qla_ice, clinfo1 = ' blk_ice: qla_ice  : ', tab3d_2 = z_qsb, clinfo2 = ' z_qsb    : ', kdim = jpl)
      CALL prt_ctl(tab3d_1 = z_qlw, clinfo1 = ' blk_ice: z_qlw    : ', tab3d_2 = dqla_ice, clinfo2 = ' dqla_ice : ', kdim = jpl)
      CALL prt_ctl(tab3d_1 = z_dqsb, clinfo1 = ' blk_ice: z_dqsb   : ', tab3d_2 = z_dqlw, clinfo2 = ' z_dqlw   : ', kdim = jpl)
      CALL prt_ctl(tab3d_1 = dqns_ice, clinfo1 = ' blk_ice: dqns_ice : ', tab3d_2 = qsr_ice, clinfo2 = ' qsr_ice  : ', kdim = jpl)
      CALL prt_ctl(tab3d_1 = ptsu, clinfo1 = ' blk_ice: ptsu     : ', tab3d_2 = qns_ice, clinfo2 = ' qns_ice  : ', kdim = jpl)
      CALL prt_ctl(tab2d_1 = tprecip, clinfo1 = ' blk_ice: tprecip  : ', tab2d_2 = sprecip, clinfo2 = ' sprecip  : ')
    END IF
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE blk_ice_flx
  SUBROUTINE blk_ice_qcn(ld_virtual_itd, ptsu, ptb, phs, phi)
    LOGICAL, INTENT(IN) :: ld_virtual_itd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: ptsu
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: ptb
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phs
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phi
    INTEGER, PARAMETER :: nit = 10
    REAL(KIND = wp), PARAMETER :: zepsilon = 0.1_wp
    INTEGER :: ji, jj, jl
    INTEGER :: iter
    REAL(KIND = wp) :: zfac, zfac2, zfac3
    REAL(KIND = wp) :: zkeff_h, ztsu, ztsu0
    REAL(KIND = wp) :: zqc, zqnet
    REAL(KIND = wp) :: zhe, zqa0
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zgfac
    !$ACC KERNELS
    zgfac(:, :, :) = 1._wp
    !$ACC END KERNELS
    IF (ld_virtual_itd) THEN
      !$ACC KERNELS
      zfac = 1._wp / (rn_cnd_s + rcnd_i)
      zfac2 = EXP(1._wp) * 0.5_wp * zepsilon
      zfac3 = 2._wp / zepsilon
      DO jl = 1, jpl
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 1, jpj
          DO ji = 1, jpi
            zhe = (rn_cnd_s * phi(ji, jj, jl) + rcnd_i * phs(ji, jj, jl)) * zfac
            IF (zhe >= zfac2) zgfac(ji, jj, jl) = MIN(2._wp, 0.5_wp * (1._wp + LOG(zhe * zfac3)))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    zfac = rcnd_i * rn_cnd_s
    DO jl = 1, jpl
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zkeff_h = zfac * zgfac(ji, jj, jl) / (rcnd_i * phs(ji, jj, jl) + rn_cnd_s * MAX(0.01, phi(ji, jj, jl)))
          ztsu = ptsu(ji, jj, jl)
          ztsu0 = ptsu(ji, jj, jl)
          zqa0 = qsr_ice(ji, jj, jl) - qtr_ice_top(ji, jj, jl) + qns_ice(ji, jj, jl)
          DO iter = 1, nit
            zqc = zkeff_h * (ztsu - ptb(ji, jj))
            zqnet = zqa0 + dqns_ice(ji, jj, jl) * (ztsu - ptsu(ji, jj, jl)) - zqc
            ztsu = ztsu - zqnet / (dqns_ice(ji, jj, jl) - zkeff_h)
          END DO
          ptsu(ji, jj, jl) = MIN(rt0, ztsu)
          qcn_ice(ji, jj, jl) = zkeff_h * (ptsu(ji, jj, jl) - ptb(ji, jj))
          qns_ice(ji, jj, jl) = qns_ice(ji, jj, jl) + dqns_ice(ji, jj, jl) * (ptsu(ji, jj, jl) - ztsu0)
          qml_ice(ji, jj, jl) = (qsr_ice(ji, jj, jl) - qtr_ice_top(ji, jj, jl) + qns_ice(ji, jj, jl) - qcn_ice(ji, jj, jl)) * &
&MAX(0._wp, SIGN(1._wp, ptsu(ji, jj, jl) - rt0))
          hfx_err_dif(ji, jj) = hfx_err_dif(ji, jj) - (dqns_ice(ji, jj, jl) * (ptsu(ji, jj, jl) - ztsu0)) * a_i_b(ji, jj, jl)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE blk_ice_qcn
  SUBROUTINE Cdn10_Lupkes2012(Cd)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: Cd
    REAL(KIND = wp), PARAMETER :: zCe = 2.23E-03_wp
    REAL(KIND = wp), PARAMETER :: znu = 1._wp
    REAL(KIND = wp), PARAMETER :: zmu = 1._wp
    REAL(KIND = wp), PARAMETER :: zbeta = 1._wp
    REAL(KIND = wp) :: zcoef
    !$ACC KERNELS
    zcoef = znu + 1._wp / (10._wp * zbeta)
    Cd(:, :) = Cd_ice + zCe * (1._wp - at_i_b(:, :)) ** zcoef * at_i_b(:, :) ** (zmu - 1._wp)
    !$ACC END KERNELS
  END SUBROUTINE Cdn10_Lupkes2012
  SUBROUTINE Cdn10_Lupkes2015(Cd, Ch)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: Cd
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: Ch
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztm_su, zst, zqo_sat, zqi_sat
    REAL(KIND = wp), PARAMETER :: z0_skin_ice = 0.69E-3_wp
    REAL(KIND = wp), PARAMETER :: z0_form_ice = 0.57E-3_wp
    REAL(KIND = wp), PARAMETER :: z0_ice = 1.00E-3_wp
    REAL(KIND = wp), PARAMETER :: zce10 = 2.80E-3_wp
    REAL(KIND = wp), PARAMETER :: zbeta = 1.1_wp
    REAL(KIND = wp), PARAMETER :: zc = 5._wp
    REAL(KIND = wp), PARAMETER :: zc2 = zc * zc
    REAL(KIND = wp), PARAMETER :: zam = 2. * zc
    REAL(KIND = wp), PARAMETER :: zah = 3. * zc
    REAL(KIND = wp), PARAMETER :: z1_alpha = 1._wp / 0.2_wp
    REAL(KIND = wp), PARAMETER :: z1_alphaf = z1_alpha
    REAL(KIND = wp), PARAMETER :: zbetah = 1.E-3_wp
    REAL(KIND = wp), PARAMETER :: zgamma = 1.25_wp
    REAL(KIND = wp), PARAMETER :: z1_gamma = 1._wp / zgamma
    REAL(KIND = wp), PARAMETER :: r1_3 = 1._wp / 3._wp
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zthetav_os, zthetav_is, zthetav_zu
    REAL(KIND = wp) :: zrib_o, zrib_i
    REAL(KIND = wp) :: zCdn_skin_ice, zCdn_form_ice, zCdn_ice
    REAL(KIND = wp) :: zChn_skin_ice, zChn_form_ice
    REAL(KIND = wp) :: z0w, z0i, zfmi, zfmw, zfhi, zfhw
    REAL(KIND = wp) :: zCdn_form_tmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    WHERE (at_i_b(:, :) > 1.E-20)
      ztm_su(:, :) = SUM(t_su(:, :, :) * a_i_b(:, :, :), dim = 3) / at_i_b(:, :)
    ELSEWHERE
      ztm_su(:, :) = rt0
    END WHERE
    zCdn_form_tmp = zce10 * (LOG(10._wp / z0_form_ice + 1._wp) / LOG(rn_zu / z0_form_ice + 1._wp)) ** 2
    zCdn_skin_ice = (vkarmn / LOG(rn_zu / z0_skin_ice + 1._wp)) ** 2
    zCdn_ice = zCdn_skin_ice
    zChn_skin_ice = vkarmn ** 2 / (LOG(rn_zu / z0_ice + 1._wp) * LOG(rn_zu * z1_alpha / z0_skin_ice + 1._wp))
    zst(:, :) = sst_m(:, :) + rt0
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('cdn10_lupkes2015', 'r0', 0, 0)
    zqo_sat(:, :) = 0.98_wp * q_sat(zst(:, :), sf(jp_slp) % fnow(:, :, 1))
    zqi_sat(:, :) = 0.98_wp * q_sat(ztm_su(:, :), sf(jp_slp) % fnow(:, :, 1))
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zthetav_os = zst(ji, jj) * (1._wp + rctv0 * zqo_sat(ji, jj))
        zthetav_is = ztm_su(ji, jj) * (1._wp + rctv0 * zqi_sat(ji, jj))
        zthetav_zu = t_zu(ji, jj) * (1._wp + rctv0 * q_zu(ji, jj))
        zrib_o = grav / zthetav_os * (zthetav_zu - zthetav_os) * rn_zu / MAX(0.5, wndm(ji, jj)) ** 2
        zrib_i = grav / zthetav_is * (zthetav_zu - zthetav_is) * rn_zu / MAX(0.5, wndm_ice(ji, jj)) ** 2
        zCdn_form_ice = zCdn_form_tmp * at_i_b(ji, jj) * (1._wp - at_i_b(ji, jj)) ** zbeta
        zChn_form_ice = zCdn_form_ice / (1._wp + (LOG(z1_alphaf) / vkarmn) * SQRT(zCdn_form_ice))
        z0w = rn_zu * EXP(- 1._wp * vkarmn / SQRT(Cdn_oce(ji, jj)))
        z0i = z0_skin_ice
        IF (zrib_o <= 0._wp) THEN
          zfmw = 1._wp - zam * zrib_o / (1._wp + 3._wp * zc2 * Cdn_oce(ji, jj) * SQRT(- zrib_o * (rn_zu / z0w + 1._wp)))
          zfhw = (1._wp + (zbetah * (zthetav_os - zthetav_zu) ** r1_3 / (Chn_oce(ji, jj) * MAX(0.01, wndm(ji, jj)))) ** zgamma) ** &
&z1_gamma
        ELSE
          zfmw = 1._wp / (1._wp + zam * zrib_o / SQRT(1._wp + zrib_o))
          zfhw = 1._wp / (1._wp + zah * zrib_o / SQRT(1._wp + zrib_o))
        END IF
        IF (zrib_i <= 0._wp) THEN
          zfmi = 1._wp - zam * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT(- zrib_i * (rn_zu / z0i + 1._wp)))
          zfhi = 1._wp - zah * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT(- zrib_i * (rn_zu / z0i + 1._wp)))
        ELSE
          zfmi = 1._wp / (1._wp + zam * zrib_i / SQRT(1._wp + zrib_i))
          zfhi = 1._wp / (1._wp + zah * zrib_i / SQRT(1._wp + zrib_i))
        END IF
        Cd(ji, jj) = zCdn_skin_ice * zfmi + zCdn_form_ice * (zfmi * at_i_b(ji, jj) + zfmw * (1._wp - at_i_b(ji, jj))) / &
&MAX(1.E-06, at_i_b(ji, jj))
        Ch(ji, jj) = zChn_skin_ice * zfhi + zChn_form_ice * (zfhi * at_i_b(ji, jj) + zfhw * (1._wp - at_i_b(ji, jj))) / &
&MAX(1.E-06, at_i_b(ji, jj))
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('sbcblk', Cd, 'T', 1., Ch, 'T', 1.)
  END SUBROUTINE Cdn10_Lupkes2015
END MODULE sbcblk
