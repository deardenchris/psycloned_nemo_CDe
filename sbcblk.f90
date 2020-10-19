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
    !  qsr(:, :) = zztmp * sbc_dcy(sf(jp_qsr) % fnow(:, :, 1)) * tmask(:, :, 1)
    qsr(:, :) = sbc_dcy(sf(jp_qsr) % fnow(:, :, 1))
    !$ACC KERNELS ! CDe
    qsr(:, :) = zztmp * qsr(:,:)  * tmask(:, :, 1)
    !$ACC END KERNELS
    ELSE
      !$ACC KERNELS ! CDe      
      qsr(:, :) = zztmp * sf(jp_qsr) % fnow(:, :, 1) * tmask(:, :, 1)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS ! CDe
    zqlw(:, :) = (sf(jp_qlw) % fnow(:, :, 1) - Stef * zst(:, :) * zst(:, :) * zst(:, :) * zst(:, :)) * tmask(:, :, 1)
!    zsq(:, :) = 0.98 * q_sat(zst(:, :), sf(jp_slp) % fnow(:, :, 1))
!    ztpot = sf(jp_tair) % fnow(:, :, 1) + gamma_moist(sf(jp_tair) % fnow(:, :, 1), sf(jp_humi) % fnow(:, :, 1)) * rn_zqt
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
    !$ACC KERNELS ! CDe
    IF (lhftau) taum(:, :) = taum(:, :) + sf(jp_tdif) % fnow(:, :, 1)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('blk_oce', 'r2', 0, 0)
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
      !$ACC KERNELS ! CDe
      zevap(:, :) = rn_efac * MAX(0._wp, zqla(:, :) * Ce_atm(:, :) * (zsq(:, :) - sf(jp_humi) % fnow(:, :, 1)))
      !CALL profile_psy_data3 % PreStart('blk_oce', 'r3', 0, 0)
      !$ACC END KERNELS
      !zqsb(:, :) = cp_air(sf(jp_humi) % fnow(:, :, 1)) * zqla(:, :) * Ch_atm(:, :) * (zst(:, :) - ztpot(:, :))
      !CALL profile_psy_data3 % PostEnd
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
    !zqla(:, :) = L_vap(zst(:, :)) * zevap(:, :)
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
    !TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !CALL profile_psy_data0 % PreStart('rho_air', 'r0', 0, 0)
    !$ACC KERNELS ! CDe
    rho_air = pslp / (R_dry * ptak * (1._wp + rctv0 * pqa))
    !$ACC END KERNELS
    !CALL profile_psy_data0 % PostEnd
  END FUNCTION rho_air
  FUNCTION cp_air(pqa)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pqa
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: cp_air
    !TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !CALL profile_psy_data0 % PreStart('cp_air', 'r0', 0, 0)
    !$ACC KERNELS
    Cp_air = Cp_dry + Cp_vap * pqa
    !$ACC END KERNELS
    !CALL profile_psy_data0 % PostEnd
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
    !TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !CALL profile_psy_data0 % PreStart('l_vap', 'r0', 0, 0)
    !$ACC KERNELS ! CDe
    L_vap = (2.501 - 0.00237 * (psst(:, :) - rt0)) * 1.E6
    !$ACC END KERNELS
    !CALL profile_psy_data0 % PostEnd
  END FUNCTION L_vap
END MODULE sbcblk
