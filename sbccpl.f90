MODULE sbccpl
  USE dom_oce
  USE sbc_oce
  USE trc_oce
  USE sbc_ice
  USE sbcapr
  USE sbcdcy
  USE sbcwave
  USE phycst
  USE cpl_oasis3
  USE geo2ocean
  USE oce, ONLY: tsn, un, vn, sshn, ub, vb, sshb, fraqsr_1lev
  USE ocealb
  USE eosbn2
  USE sbcrnf, ONLY: l_rnfcpl
  USE sbcisf, ONLY: l_isfcpl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_cpl_init
  PUBLIC :: sbc_cpl_rcv
  PUBLIC :: sbc_cpl_snd
  PUBLIC :: sbc_cpl_ice_tau
  PUBLIC :: sbc_cpl_ice_flx
  PUBLIC :: sbc_cpl_alloc
  INTEGER, PARAMETER :: jpr_otx1 = 1
  INTEGER, PARAMETER :: jpr_oty1 = 2
  INTEGER, PARAMETER :: jpr_otz1 = 3
  INTEGER, PARAMETER :: jpr_otx2 = 4
  INTEGER, PARAMETER :: jpr_oty2 = 5
  INTEGER, PARAMETER :: jpr_otz2 = 6
  INTEGER, PARAMETER :: jpr_itx1 = 7
  INTEGER, PARAMETER :: jpr_ity1 = 8
  INTEGER, PARAMETER :: jpr_itz1 = 9
  INTEGER, PARAMETER :: jpr_itx2 = 10
  INTEGER, PARAMETER :: jpr_ity2 = 11
  INTEGER, PARAMETER :: jpr_itz2 = 12
  INTEGER, PARAMETER :: jpr_qsroce = 13
  INTEGER, PARAMETER :: jpr_qsrice = 14
  INTEGER, PARAMETER :: jpr_qsrmix = 15
  INTEGER, PARAMETER :: jpr_qnsoce = 16
  INTEGER, PARAMETER :: jpr_qnsice = 17
  INTEGER, PARAMETER :: jpr_qnsmix = 18
  INTEGER, PARAMETER :: jpr_rain = 19
  INTEGER, PARAMETER :: jpr_snow = 20
  INTEGER, PARAMETER :: jpr_tevp = 21
  INTEGER, PARAMETER :: jpr_ievp = 22
  INTEGER, PARAMETER :: jpr_sbpr = 23
  INTEGER, PARAMETER :: jpr_semp = 24
  INTEGER, PARAMETER :: jpr_oemp = 25
  INTEGER, PARAMETER :: jpr_w10m = 26
  INTEGER, PARAMETER :: jpr_dqnsdt = 27
  INTEGER, PARAMETER :: jpr_rnf = 28
  INTEGER, PARAMETER :: jpr_cal = 29
  INTEGER, PARAMETER :: jpr_taum = 30
  INTEGER, PARAMETER :: jpr_co2 = 31
  INTEGER, PARAMETER :: jpr_topm = 32
  INTEGER, PARAMETER :: jpr_botm = 33
  INTEGER, PARAMETER :: jpr_sflx = 34
  INTEGER, PARAMETER :: jpr_toce = 35
  INTEGER, PARAMETER :: jpr_soce = 36
  INTEGER, PARAMETER :: jpr_ocx1 = 37
  INTEGER, PARAMETER :: jpr_ocy1 = 38
  INTEGER, PARAMETER :: jpr_ssh = 39
  INTEGER, PARAMETER :: jpr_fice = 40
  INTEGER, PARAMETER :: jpr_e3t1st = 41
  INTEGER, PARAMETER :: jpr_fraqsr = 42
  INTEGER, PARAMETER :: jpr_mslp = 43
  INTEGER, PARAMETER :: jpr_hsig = 44
  INTEGER, PARAMETER :: jpr_phioc = 45
  INTEGER, PARAMETER :: jpr_sdrftx = 46
  INTEGER, PARAMETER :: jpr_sdrfty = 47
  INTEGER, PARAMETER :: jpr_wper = 48
  INTEGER, PARAMETER :: jpr_wnum = 49
  INTEGER, PARAMETER :: jpr_tauwoc = 50
  INTEGER, PARAMETER :: jpr_wdrag = 51
  INTEGER, PARAMETER :: jpr_isf = 52
  INTEGER, PARAMETER :: jpr_icb = 53
  INTEGER, PARAMETER :: jpr_wfreq = 54
  INTEGER, PARAMETER :: jpr_tauwx = 55
  INTEGER, PARAMETER :: jpr_tauwy = 56
  INTEGER, PARAMETER :: jpr_ts_ice = 57
  INTEGER, PARAMETER :: jprcv = 57
  INTEGER, PARAMETER :: jps_fice = 1
  INTEGER, PARAMETER :: jps_toce = 2
  INTEGER, PARAMETER :: jps_tice = 3
  INTEGER, PARAMETER :: jps_tmix = 4
  INTEGER, PARAMETER :: jps_albice = 5
  INTEGER, PARAMETER :: jps_albmix = 6
  INTEGER, PARAMETER :: jps_hice = 7
  INTEGER, PARAMETER :: jps_hsnw = 8
  INTEGER, PARAMETER :: jps_ocx1 = 9
  INTEGER, PARAMETER :: jps_ocy1 = 10
  INTEGER, PARAMETER :: jps_ocz1 = 11
  INTEGER, PARAMETER :: jps_ivx1 = 12
  INTEGER, PARAMETER :: jps_ivy1 = 13
  INTEGER, PARAMETER :: jps_ivz1 = 14
  INTEGER, PARAMETER :: jps_co2 = 15
  INTEGER, PARAMETER :: jps_soce = 16
  INTEGER, PARAMETER :: jps_ssh = 17
  INTEGER, PARAMETER :: jps_qsroce = 18
  INTEGER, PARAMETER :: jps_qnsoce = 19
  INTEGER, PARAMETER :: jps_oemp = 20
  INTEGER, PARAMETER :: jps_sflx = 21
  INTEGER, PARAMETER :: jps_otx1 = 22
  INTEGER, PARAMETER :: jps_oty1 = 23
  INTEGER, PARAMETER :: jps_rnf = 24
  INTEGER, PARAMETER :: jps_taum = 25
  INTEGER, PARAMETER :: jps_fice2 = 26
  INTEGER, PARAMETER :: jps_e3t1st = 27
  INTEGER, PARAMETER :: jps_fraqsr = 28
  INTEGER, PARAMETER :: jps_ficet = 29
  INTEGER, PARAMETER :: jps_ocxw = 30
  INTEGER, PARAMETER :: jps_ocyw = 31
  INTEGER, PARAMETER :: jps_wlev = 32
  INTEGER, PARAMETER :: jps_fice1 = 33
  INTEGER, PARAMETER :: jps_a_p = 34
  INTEGER, PARAMETER :: jps_ht_p = 35
  INTEGER, PARAMETER :: jps_kice = 36
  INTEGER, PARAMETER :: jps_sstfrz = 37
  INTEGER, PARAMETER :: jps_ttilyr = 38
  INTEGER, PARAMETER :: jpsnd = 38
  TYPE :: FLD_C
    CHARACTER(LEN = 32) :: cldes
    CHARACTER(LEN = 32) :: clcat
    CHARACTER(LEN = 32) :: clvref
    CHARACTER(LEN = 32) :: clvor
    CHARACTER(LEN = 32) :: clvgrd
  END TYPE FLD_C
  TYPE(FLD_C) :: sn_snd_temp, sn_snd_alb, sn_snd_thick, sn_snd_crt, sn_snd_co2, sn_snd_thick1, sn_snd_cond, sn_snd_mpnd, sn_snd_sstfrz, sn_snd_ttilyr
  TYPE(FLD_C) :: sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau, sn_rcv_tauw, sn_rcv_dqnsdt, sn_rcv_qsr, sn_rcv_qns, sn_rcv_emp, sn_rcv_rnf, sn_rcv_ts_ice
  TYPE(FLD_C) :: sn_rcv_cal, sn_rcv_iceflx, sn_rcv_co2, sn_rcv_mslp, sn_rcv_icb, sn_rcv_isf
  TYPE(FLD_C) :: sn_snd_ifrac, sn_snd_crtw, sn_snd_wlev
  TYPE(FLD_C) :: sn_rcv_hsig, sn_rcv_phioc, sn_rcv_sdrfx, sn_rcv_sdrfy, sn_rcv_wper, sn_rcv_wnum, sn_rcv_tauwoc, sn_rcv_wdrag, sn_rcv_wfreq
  INTEGER :: nn_cplmodel
  LOGICAL :: ln_usecplmask
  TYPE :: DYNARR
    REAL(KIND = wp), POINTER, DIMENSION(:, :, :) :: z3
  END TYPE DYNARR
  TYPE(DYNARR), SAVE, DIMENSION(jprcv) :: frcv
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: alb_oce_mix
  REAL(KIND = wp) :: rpref = 101000._wp
  REAL(KIND = wp) :: r1_grau
  INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: nrcvinfo
  CONTAINS
  INTEGER FUNCTION sbc_cpl_alloc()
    INTEGER :: ierr(4)
    ierr(:) = 0
    ALLOCATE(alb_oce_mix(jpi, jpj), nrcvinfo(jprcv), STAT = ierr(1))
    ALLOCATE(a_i(jpi, jpj, 1), STAT = ierr(2))
    ALLOCATE(xcplmask(jpi, jpj, 0 : nn_cplmodel), STAT = ierr(3))
    IF (.NOT. ln_apr_dyn) ALLOCATE(ssh_ib(jpi, jpj), ssh_ibb(jpi, jpj), apr(jpi, jpj), STAT = ierr(4))
    sbc_cpl_alloc = MAXVAL(ierr)
    CALL mpp_sum('sbccpl', sbc_cpl_alloc)
    IF (sbc_cpl_alloc > 0) CALL ctl_warn('sbc_cpl_alloc: allocation of arrays failed')
  END FUNCTION sbc_cpl_alloc
  SUBROUTINE sbc_cpl_init(k_ice)
    INTEGER, INTENT(IN) :: k_ice
    INTEGER :: jn
    INTEGER :: ios, inum
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zacs, zaos
    NAMELIST /namsbc_cpl/ sn_snd_temp, sn_snd_alb, sn_snd_thick, sn_snd_crt, sn_snd_co2, sn_snd_ttilyr, sn_snd_cond, sn_snd_mpnd, sn_snd_sstfrz, sn_snd_thick1, sn_snd_ifrac, sn_snd_crtw, sn_snd_wlev, sn_rcv_hsig, sn_rcv_phioc, sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau, sn_rcv_dqnsdt, sn_rcv_qsr, sn_rcv_sdrfx, sn_rcv_sdrfy, sn_rcv_wper, sn_rcv_wnum, sn_rcv_tauwoc, sn_rcv_wdrag, sn_rcv_qns, sn_rcv_emp, sn_rcv_rnf, sn_rcv_cal, sn_rcv_iceflx, sn_rcv_co2, nn_cplmodel, ln_usecplmask, sn_rcv_mslp, sn_rcv_icb, sn_rcv_isf, sn_rcv_wfreq, sn_rcv_tauw, nn_cats_cpl, sn_rcv_ts_ice
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_cpl, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_cpl in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_cpl, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_cpl in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_cpl)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'sbc_cpl_init : namsbc_cpl namelist '
      WRITE(numout, *) '~~~~~~~~~~~~'
    END IF
    IF (lwp .AND. ln_cpl) THEN
      WRITE(numout, *) '  received fields (mutiple ice categogies)'
      WRITE(numout, *) '      10m wind module                 = ', TRIM(sn_rcv_w10m % cldes), ' (', TRIM(sn_rcv_w10m % clcat), ')'
      WRITE(numout, *) '      stress module                   = ', TRIM(sn_rcv_taumod % cldes), ' (', TRIM(sn_rcv_taumod % clcat), ')'
      WRITE(numout, *) '      surface stress                  = ', TRIM(sn_rcv_tau % cldes), ' (', TRIM(sn_rcv_tau % clcat), ')'
      WRITE(numout, *) '                     - referential    = ', sn_rcv_tau % clvref
      WRITE(numout, *) '                     - orientation    = ', sn_rcv_tau % clvor
      WRITE(numout, *) '                     - mesh           = ', sn_rcv_tau % clvgrd
      WRITE(numout, *) '      non-solar heat flux sensitivity = ', TRIM(sn_rcv_dqnsdt % cldes), ' (', TRIM(sn_rcv_dqnsdt % clcat), ')'
      WRITE(numout, *) '      solar heat flux                 = ', TRIM(sn_rcv_qsr % cldes), ' (', TRIM(sn_rcv_qsr % clcat), ')'
      WRITE(numout, *) '      non-solar heat flux             = ', TRIM(sn_rcv_qns % cldes), ' (', TRIM(sn_rcv_qns % clcat), ')'
      WRITE(numout, *) '      freshwater budget               = ', TRIM(sn_rcv_emp % cldes), ' (', TRIM(sn_rcv_emp % clcat), ')'
      WRITE(numout, *) '      runoffs                         = ', TRIM(sn_rcv_rnf % cldes), ' (', TRIM(sn_rcv_rnf % clcat), ')'
      WRITE(numout, *) '      calving                         = ', TRIM(sn_rcv_cal % cldes), ' (', TRIM(sn_rcv_cal % clcat), ')'
      WRITE(numout, *) '      iceberg                         = ', TRIM(sn_rcv_icb % cldes), ' (', TRIM(sn_rcv_icb % clcat), ')'
      WRITE(numout, *) '      ice shelf                       = ', TRIM(sn_rcv_isf % cldes), ' (', TRIM(sn_rcv_isf % clcat), ')'
      WRITE(numout, *) '      sea ice heat fluxes             = ', TRIM(sn_rcv_iceflx % cldes), ' (', TRIM(sn_rcv_iceflx % clcat), ')'
      WRITE(numout, *) '      atm co2                         = ', TRIM(sn_rcv_co2 % cldes), ' (', TRIM(sn_rcv_co2 % clcat), ')'
      WRITE(numout, *) '      significant wave heigth         = ', TRIM(sn_rcv_hsig % cldes), ' (', TRIM(sn_rcv_hsig % clcat), ')'
      WRITE(numout, *) '      wave to oce energy flux         = ', TRIM(sn_rcv_phioc % cldes), ' (', TRIM(sn_rcv_phioc % clcat), ')'
      WRITE(numout, *) '      Surface Stokes drift grid u     = ', TRIM(sn_rcv_sdrfx % cldes), ' (', TRIM(sn_rcv_sdrfx % clcat), ')'
      WRITE(numout, *) '      Surface Stokes drift grid v     = ', TRIM(sn_rcv_sdrfy % cldes), ' (', TRIM(sn_rcv_sdrfy % clcat), ')'
      WRITE(numout, *) '      Mean wave period                = ', TRIM(sn_rcv_wper % cldes), ' (', TRIM(sn_rcv_wper % clcat), ')'
      WRITE(numout, *) '      Mean wave number                = ', TRIM(sn_rcv_wnum % cldes), ' (', TRIM(sn_rcv_wnum % clcat), ')'
      WRITE(numout, *) '      Wave peak frequency             = ', TRIM(sn_rcv_wfreq % cldes), ' (', TRIM(sn_rcv_wfreq % clcat), ')'
      WRITE(numout, *) '      Stress frac adsorbed by waves   = ', TRIM(sn_rcv_tauwoc % cldes), ' (', TRIM(sn_rcv_tauwoc % clcat), ')'
      WRITE(numout, *) '      Stress components by waves      = ', TRIM(sn_rcv_tauw % cldes), ' (', TRIM(sn_rcv_tauw % clcat), ')'
      WRITE(numout, *) '      Neutral surf drag coefficient   = ', TRIM(sn_rcv_wdrag % cldes), ' (', TRIM(sn_rcv_wdrag % clcat), ')'
      WRITE(numout, *) '      Sea ice surface skin temperature= ', TRIM(sn_rcv_ts_ice % cldes), ' (', TRIM(sn_rcv_ts_ice % clcat), ')'
      WRITE(numout, *) '  sent fields (multiple ice categories)'
      WRITE(numout, *) '      surface temperature             = ', TRIM(sn_snd_temp % cldes), ' (', TRIM(sn_snd_temp % clcat), ')'
      WRITE(numout, *) '      top ice layer temperature       = ', TRIM(sn_snd_ttilyr % cldes), ' (', TRIM(sn_snd_ttilyr % clcat), ')'
      WRITE(numout, *) '      albedo                          = ', TRIM(sn_snd_alb % cldes), ' (', TRIM(sn_snd_alb % clcat), ')'
      WRITE(numout, *) '      ice/snow thickness              = ', TRIM(sn_snd_thick % cldes), ' (', TRIM(sn_snd_thick % clcat), ')'
      WRITE(numout, *) '      total ice fraction              = ', TRIM(sn_snd_ifrac % cldes), ' (', TRIM(sn_snd_ifrac % clcat), ')'
      WRITE(numout, *) '      surface current                 = ', TRIM(sn_snd_crt % cldes), ' (', TRIM(sn_snd_crt % clcat), ')'
      WRITE(numout, *) '                      - referential   = ', sn_snd_crt % clvref
      WRITE(numout, *) '                      - orientation   = ', sn_snd_crt % clvor
      WRITE(numout, *) '                      - mesh          = ', sn_snd_crt % clvgrd
      WRITE(numout, *) '      oce co2 flux                    = ', TRIM(sn_snd_co2 % cldes), ' (', TRIM(sn_snd_co2 % clcat), ')'
      WRITE(numout, *) '      ice effective conductivity      = ', TRIM(sn_snd_cond % cldes), ' (', TRIM(sn_snd_cond % clcat), ')'
      WRITE(numout, *) '      meltponds fraction and depth    = ', TRIM(sn_snd_mpnd % cldes), ' (', TRIM(sn_snd_mpnd % clcat), ')'
      WRITE(numout, *) '      sea surface freezing temp       = ', TRIM(sn_snd_sstfrz % cldes), ' (', TRIM(sn_snd_sstfrz % clcat), ')'
      WRITE(numout, *) '      water level                     = ', TRIM(sn_snd_wlev % cldes), ' (', TRIM(sn_snd_wlev % clcat), ')'
      WRITE(numout, *) '      mean sea level pressure         = ', TRIM(sn_rcv_mslp % cldes), ' (', TRIM(sn_rcv_mslp % clcat), ')'
      WRITE(numout, *) '      surface current to waves        = ', TRIM(sn_snd_crtw % cldes), ' (', TRIM(sn_snd_crtw % clcat), ')'
      WRITE(numout, *) '                      - referential   = ', sn_snd_crtw % clvref
      WRITE(numout, *) '                      - orientation   = ', sn_snd_crtw % clvor
      WRITE(numout, *) '                      - mesh          = ', sn_snd_crtw % clvgrd
      WRITE(numout, *) '  nn_cplmodel                         = ', nn_cplmodel
      WRITE(numout, *) '  ln_usecplmask                       = ', ln_usecplmask
      WRITE(numout, *) '  nn_cats_cpl                         = ', nn_cats_cpl
    END IF
    IF (sbc_cpl_alloc() /= 0) CALL ctl_stop('STOP', 'sbc_cpl_alloc : unable to allocate arrays')
    nrcvinfo(:) = OASIS_idle
    srcv(:) % laction = .FALSE.
    srcv(:) % clgrid = 'T'
    srcv(:) % nsgn = 1.
    srcv(:) % nct = 1
    srcv(jpr_otx1) % clname = 'O_OTaux1'
    srcv(jpr_oty1) % clname = 'O_OTauy1'
    srcv(jpr_otz1) % clname = 'O_OTauz1'
    srcv(jpr_otx2) % clname = 'O_OTaux2'
    srcv(jpr_oty2) % clname = 'O_OTauy2'
    srcv(jpr_otz2) % clname = 'O_OTauz2'
    srcv(jpr_itx1) % clname = 'O_ITaux1'
    srcv(jpr_ity1) % clname = 'O_ITauy1'
    srcv(jpr_itz1) % clname = 'O_ITauz1'
    srcv(jpr_itx2) % clname = 'O_ITaux2'
    srcv(jpr_ity2) % clname = 'O_ITauy2'
    srcv(jpr_itz2) % clname = 'O_ITauz2'
    IF (TRIM(sn_rcv_tau % cldes) == 'oce only' .OR. TRIM(sn_rcv_tau % cldes) == 'oce and ice') THEN
      IF (TRIM(sn_rcv_tau % clvor) == 'local grid') srcv(jpr_otx1 : jpr_itz2) % nsgn = - 1.
      SELECT CASE (TRIM(sn_rcv_tau % clvgrd))
      CASE ('T')
        srcv(jpr_otx1 : jpr_itz2) % clgrid = 'T'
        srcv(jpr_otx1 : jpr_otz1) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz1) % laction = .TRUE.
      CASE ('U,V')
        srcv(jpr_otx1 : jpr_otz1) % clgrid = 'U'
        srcv(jpr_otx2 : jpr_otz2) % clgrid = 'V'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'U'
        srcv(jpr_itx2 : jpr_itz2) % clgrid = 'V'
        srcv(jpr_otx1 : jpr_itz2) % laction = .TRUE.
      CASE ('U,V,T')
        srcv(jpr_otx1 : jpr_otz1) % clgrid = 'U'
        srcv(jpr_otx2 : jpr_otz2) % clgrid = 'V'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'T'
        srcv(jpr_otx1 : jpr_otz2) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz1) % laction = .TRUE.
      CASE ('U,V,I')
        srcv(jpr_otx1 : jpr_otz1) % clgrid = 'U'
        srcv(jpr_otx2 : jpr_otz2) % clgrid = 'V'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'I'
        srcv(jpr_otx1 : jpr_otz2) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz1) % laction = .TRUE.
      CASE ('U,V,F')
        srcv(jpr_otx1 : jpr_otz1) % clgrid = 'U'
        srcv(jpr_otx2 : jpr_otz2) % clgrid = 'V'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'F'
        srcv(jpr_otx1 : jpr_otz2) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz1) % laction = .TRUE.
      CASE ('T,I')
        srcv(jpr_otx1 : jpr_itz2) % clgrid = 'T'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'I'
        srcv(jpr_otx1 : jpr_otz1) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz1) % laction = .TRUE.
      CASE ('T,F')
        srcv(jpr_otx1 : jpr_itz2) % clgrid = 'T'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'F'
        srcv(jpr_otx1 : jpr_otz1) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz1) % laction = .TRUE.
      CASE ('T,U,V')
        srcv(jpr_otx1 : jpr_otz1) % clgrid = 'T'
        srcv(jpr_itx1 : jpr_itz1) % clgrid = 'U'
        srcv(jpr_itx2 : jpr_itz2) % clgrid = 'V'
        srcv(jpr_otx1 : jpr_otz1) % laction = .TRUE.
        srcv(jpr_itx1 : jpr_itz2) % laction = .TRUE.
      CASE DEFAULT
        CALL ctl_stop('sbc_cpl_init: wrong definition of sn_rcv_tau%clvgrd')
      END SELECT
      IF (TRIM(sn_rcv_tau % clvref) == 'spherical') srcv((/jpr_otz1, jpr_otz2, jpr_itz1, jpr_itz2/)) % laction = .FALSE.
      IF (TRIM(sn_rcv_tau % clvor) == 'local grid') THEN
        srcv(jpr_otx2 : jpr_otz2) % laction = .FALSE.
        srcv(jpr_itx2 : jpr_itz2) % laction = .FALSE.
        srcv(jpr_oty1) % clgrid = srcv(jpr_oty2) % clgrid
        srcv(jpr_ity1) % clgrid = srcv(jpr_ity2) % clgrid
      END IF
      IF (TRIM(sn_rcv_tau % cldes) /= 'oce and ice') THEN
        srcv(jpr_itx1 : jpr_itz2) % laction = .FALSE.
        srcv(jpr_itx1) % clgrid = 'U'
        srcv(jpr_ity1) % clgrid = 'V'
      END IF
    END IF
    srcv(jpr_rain) % clname = 'OTotRain'
    srcv(jpr_snow) % clname = 'OTotSnow'
    srcv(jpr_tevp) % clname = 'OTotEvap'
    srcv(jpr_ievp) % clname = 'OIceEvap'
    srcv(jpr_sbpr) % clname = 'OSubMPre'
    srcv(jpr_semp) % clname = 'OISubMSn'
    srcv(jpr_oemp) % clname = 'OOEvaMPr'
    SELECT CASE (TRIM(sn_rcv_emp % cldes))
    CASE ('none')
    CASE ('oce only')
      srcv(jpr_oemp) % laction = .TRUE.
    CASE ('conservative')
      srcv((/jpr_rain, jpr_snow, jpr_ievp, jpr_tevp/)) % laction = .TRUE.
      IF (k_ice <= 1) srcv(jpr_ievp) % laction = .FALSE.
    CASE ('oce and ice')
      srcv((/jpr_ievp, jpr_sbpr, jpr_semp, jpr_oemp/)) % laction = .TRUE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_rcv_emp%cldes')
    END SELECT
    srcv(jpr_rnf) % clname = 'O_Runoff'
    IF (TRIM(sn_rcv_rnf % cldes) == 'coupled') THEN
      srcv(jpr_rnf) % laction = .TRUE.
      l_rnfcpl = .TRUE.
      ln_rnf = nn_components /= jp_iam_sas
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   runoffs received from oasis -> force ln_rnf = ', ln_rnf
    END IF
    srcv(jpr_cal) % clname = 'OCalving'
    IF (TRIM(sn_rcv_cal % cldes) == 'coupled') srcv(jpr_cal) % laction = .TRUE.
    srcv(jpr_isf) % clname = 'OIcshelf'
    IF (TRIM(sn_rcv_isf % cldes) == 'coupled') srcv(jpr_isf) % laction = .TRUE.
    srcv(jpr_icb) % clname = 'OIceberg'
    IF (TRIM(sn_rcv_icb % cldes) == 'coupled') srcv(jpr_icb) % laction = .TRUE.
    IF (srcv(jpr_isf) % laction .AND. ln_isf) THEN
      l_isfcpl = .TRUE.
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   iceshelf received from oasis '
    END IF
    srcv(jpr_qnsoce) % clname = 'O_QnsOce'
    srcv(jpr_qnsice) % clname = 'O_QnsIce'
    srcv(jpr_qnsmix) % clname = 'O_QnsMix'
    SELECT CASE (TRIM(sn_rcv_qns % cldes))
    CASE ('none')
    CASE ('oce only')
      srcv(jpr_qnsoce) % laction = .TRUE.
    CASE ('conservative')
      srcv((/jpr_qnsice, jpr_qnsmix/)) % laction = .TRUE.
    CASE ('oce and ice')
      srcv((/jpr_qnsice, jpr_qnsoce/)) % laction = .TRUE.
    CASE ('mixed oce-ice')
      srcv(jpr_qnsmix) % laction = .TRUE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_rcv_qns%cldes')
    END SELECT
    IF (TRIM(sn_rcv_qns % cldes) == 'mixed oce-ice' .AND. nn_cats_cpl > 1) CALL ctl_stop('sbc_cpl_init: sn_rcv_qns%cldes not currently allowed to be mixed oce-ice for multi-category ice')
    srcv(jpr_qsroce) % clname = 'O_QsrOce'
    srcv(jpr_qsrice) % clname = 'O_QsrIce'
    srcv(jpr_qsrmix) % clname = 'O_QsrMix'
    SELECT CASE (TRIM(sn_rcv_qsr % cldes))
    CASE ('none')
    CASE ('oce only')
      srcv(jpr_qsroce) % laction = .TRUE.
    CASE ('conservative')
      srcv((/jpr_qsrice, jpr_qsrmix/)) % laction = .TRUE.
    CASE ('oce and ice')
      srcv((/jpr_qsrice, jpr_qsroce/)) % laction = .TRUE.
    CASE ('mixed oce-ice')
      srcv(jpr_qsrmix) % laction = .TRUE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_rcv_qsr%cldes')
    END SELECT
    IF (TRIM(sn_rcv_qsr % cldes) == 'mixed oce-ice' .AND. nn_cats_cpl > 1) CALL ctl_stop('sbc_cpl_init: sn_rcv_qsr%cldes not currently allowed to be mixed oce-ice for multi-category ice')
    srcv(jpr_dqnsdt) % clname = 'O_dQnsdT'
    IF (TRIM(sn_rcv_dqnsdt % cldes) == 'coupled') srcv(jpr_dqnsdt) % laction = .TRUE.
    IF (TRIM(sn_rcv_dqnsdt % cldes) == 'none' .AND. TRIM(sn_rcv_qns % cldes) == 'mixed oce-ice') CALL ctl_stop('sbc_cpl_init: namsbc_cpl namelist mismatch between sn_rcv_qns%cldes and sn_rcv_dqnsdt%cldes')
    srcv(jpr_w10m) % clname = 'O_Wind10'
    IF (TRIM(sn_rcv_w10m % cldes) == 'coupled') srcv(jpr_w10m) % laction = .TRUE.
    srcv(jpr_taum) % clname = 'O_TauMod'
    IF (TRIM(sn_rcv_taumod % cldes) == 'coupled') srcv(jpr_taum) % laction = .TRUE.
    lhftau = srcv(jpr_taum) % laction
    srcv(jpr_co2) % clname = 'O_AtmCO2'
    IF (TRIM(sn_rcv_co2 % cldes) == 'coupled') THEN
      srcv(jpr_co2) % laction = .TRUE.
      l_co2cpl = .TRUE.
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) '   Atmospheric pco2 received from oasis '
      IF (lwp) WRITE(numout, *)
    END IF
    srcv(jpr_mslp) % clname = 'O_MSLP'
    IF (TRIM(sn_rcv_mslp % cldes) == 'coupled') srcv(jpr_mslp) % laction = .TRUE.
    srcv(jpr_topm) % clname = 'OTopMlt'
    srcv(jpr_botm) % clname = 'OBotMlt'
    IF (TRIM(sn_rcv_iceflx % cldes) == 'coupled') THEN
      IF (TRIM(sn_rcv_iceflx % clcat) == 'yes') THEN
        srcv(jpr_topm : jpr_botm) % nct = nn_cats_cpl
      ELSE
        CALL ctl_stop('sbc_cpl_init: sn_rcv_iceflx%clcat should always be set to yes currently')
      END IF
      srcv(jpr_topm : jpr_botm) % laction = .TRUE.
    END IF
    srcv(jpr_ts_ice) % clname = 'OTsfIce'
    IF (TRIM(sn_rcv_ts_ice % cldes) == 'ice') srcv(jpr_ts_ice) % laction = .TRUE.
    IF (TRIM(sn_rcv_ts_ice % clcat) == 'yes') srcv(jpr_ts_ice) % nct = nn_cats_cpl
    IF (TRIM(sn_rcv_emp % clcat) == 'yes') srcv(jpr_ievp) % nct = nn_cats_cpl
    srcv(jpr_hsig) % clname = 'O_Hsigwa'
    IF (TRIM(sn_rcv_hsig % cldes) == 'coupled') THEN
      srcv(jpr_hsig) % laction = .TRUE.
      cpl_hsig = .TRUE.
    END IF
    srcv(jpr_phioc) % clname = 'O_PhiOce'
    IF (TRIM(sn_rcv_phioc % cldes) == 'coupled') THEN
      srcv(jpr_phioc) % laction = .TRUE.
      cpl_phioc = .TRUE.
    END IF
    srcv(jpr_sdrftx) % clname = 'O_Sdrfx'
    IF (TRIM(sn_rcv_sdrfx % cldes) == 'coupled') THEN
      srcv(jpr_sdrftx) % laction = .TRUE.
      cpl_sdrftx = .TRUE.
    END IF
    srcv(jpr_sdrfty) % clname = 'O_Sdrfy'
    IF (TRIM(sn_rcv_sdrfy % cldes) == 'coupled') THEN
      srcv(jpr_sdrfty) % laction = .TRUE.
      cpl_sdrfty = .TRUE.
    END IF
    srcv(jpr_wper) % clname = 'O_WPer'
    IF (TRIM(sn_rcv_wper % cldes) == 'coupled') THEN
      srcv(jpr_wper) % laction = .TRUE.
      cpl_wper = .TRUE.
    END IF
    srcv(jpr_wfreq) % clname = 'O_WFreq'
    IF (TRIM(sn_rcv_wfreq % cldes) == 'coupled') THEN
      srcv(jpr_wfreq) % laction = .TRUE.
      cpl_wfreq = .TRUE.
    END IF
    srcv(jpr_wnum) % clname = 'O_WNum'
    IF (TRIM(sn_rcv_wnum % cldes) == 'coupled') THEN
      srcv(jpr_wnum) % laction = .TRUE.
      cpl_wnum = .TRUE.
    END IF
    srcv(jpr_tauwoc) % clname = 'O_TauOce'
    IF (TRIM(sn_rcv_tauwoc % cldes) == 'coupled') THEN
      srcv(jpr_tauwoc) % laction = .TRUE.
      cpl_tauwoc = .TRUE.
    END IF
    srcv(jpr_tauwx) % clname = 'O_Tauwx'
    srcv(jpr_tauwy) % clname = 'O_Tauwy'
    IF (TRIM(sn_rcv_tauw % cldes) == 'coupled') THEN
      srcv(jpr_tauwx) % laction = .TRUE.
      srcv(jpr_tauwy) % laction = .TRUE.
      cpl_tauw = .TRUE.
    END IF
    srcv(jpr_wdrag) % clname = 'O_WDrag'
    IF (TRIM(sn_rcv_wdrag % cldes) == 'coupled') THEN
      srcv(jpr_wdrag) % laction = .TRUE.
      cpl_wdrag = .TRUE.
    END IF
    IF (srcv(jpr_tauwoc) % laction .AND. srcv(jpr_tauwx) % laction .AND. srcv(jpr_tauwy) % laction) CALL ctl_stop('More than one method for modifying the ocean stress has been selected ', '(sn_rcv_tauwoc=coupled and sn_rcv_tauw=coupled)')
    srcv(jpr_sflx) % clname = 'O_SFLX'
    srcv(jpr_fice) % clname = 'RIceFrc'
    IF (nn_components == jp_iam_opa) THEN
      srcv(:) % laction = .FALSE.
      srcv(:) % clgrid = 'T'
      srcv(:) % nsgn = 1.
      srcv((/jpr_qsroce, jpr_qnsoce, jpr_oemp, jpr_sflx, jpr_fice, jpr_otx1, jpr_oty1, jpr_taum/)) % laction = .TRUE.
      srcv(jpr_otx1) % clgrid = 'U'
      srcv(jpr_oty1) % clgrid = 'V'
      srcv((/jpr_otx1, jpr_oty1/)) % nsgn = - 1.
      sn_rcv_tau % clvgrd = 'U,V'
      sn_rcv_tau % clvor = 'local grid'
      sn_rcv_tau % clvref = 'spherical'
      sn_rcv_emp % cldes = 'oce only'
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '               Special conditions for SAS-OPA coupling  '
        WRITE(numout, *) '               OPA component  '
        WRITE(numout, *)
        WRITE(numout, *) '  received fields from SAS component '
        WRITE(numout, *) '                  ice cover '
        WRITE(numout, *) '                  oce only EMP  '
        WRITE(numout, *) '                  salt flux  '
        WRITE(numout, *) '                  mixed oce-ice solar flux  '
        WRITE(numout, *) '                  mixed oce-ice non solar flux  '
        WRITE(numout, *) '                  wind stress U,V on local grid and sperical coordinates '
        WRITE(numout, *) '                  wind stress module'
        WRITE(numout, *)
      END IF
    END IF
    srcv(jpr_toce) % clname = 'I_SSTSST'
    srcv(jpr_soce) % clname = 'I_SSSal'
    srcv(jpr_ocx1) % clname = 'I_OCurx1'
    srcv(jpr_ocy1) % clname = 'I_OCury1'
    srcv(jpr_ssh) % clname = 'I_SSHght'
    srcv(jpr_e3t1st) % clname = 'I_E3T1st'
    srcv(jpr_fraqsr) % clname = 'I_FraQsr'
    IF (nn_components == jp_iam_sas) THEN
      IF (.NOT. ln_cpl) srcv(:) % laction = .FALSE.
      IF (.NOT. ln_cpl) srcv(:) % clgrid = 'T'
      IF (.NOT. ln_cpl) srcv(:) % nsgn = 1.
      srcv((/jpr_toce, jpr_soce, jpr_ssh, jpr_fraqsr, jpr_ocx1, jpr_ocy1/)) % laction = .TRUE.
      srcv(jpr_e3t1st) % laction = .NOT. ln_linssh
      srcv(jpr_ocx1) % clgrid = 'U'
      srcv(jpr_ocy1) % clgrid = 'V'
      srcv(jpr_ocx1 : jpr_ocy1) % nsgn = - 1.
      DO jn = 1, jprcv
        IF (srcv(jn) % clname(1 : 1) == "O") srcv(jn) % clname = "S" // srcv(jn) % clname(2 : LEN(srcv(jn) % clname))
      END DO
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '               Special conditions for SAS-OPA coupling  '
        WRITE(numout, *) '               SAS component  '
        WRITE(numout, *)
        IF (.NOT. ln_cpl) THEN
          WRITE(numout, *) '  received fields from OPA component '
        ELSE
          WRITE(numout, *) '  Additional received fields from OPA component : '
        END IF
        WRITE(numout, *) '               sea surface temperature (Celsius) '
        WRITE(numout, *) '               sea surface salinity '
        WRITE(numout, *) '               surface currents '
        WRITE(numout, *) '               sea surface height '
        WRITE(numout, *) '               thickness of first ocean T level '
        WRITE(numout, *) '               fraction of solar net radiation absorbed in the first ocean level'
        WRITE(numout, *)
      END IF
    END IF
    DO jn = 1, jprcv
      IF (srcv(jn) % laction) ALLOCATE(frcv(jn) % z3(jpi, jpj, srcv(jn) % nct))
    END DO
    IF (.NOT. srcv(jpr_taum) % laction) ALLOCATE(frcv(jpr_taum) % z3(jpi, jpj, srcv(jpr_taum) % nct))
    IF (.NOT. srcv(jpr_w10m) % laction) ALLOCATE(frcv(jpr_w10m) % z3(jpi, jpj, srcv(jpr_w10m) % nct))
    IF (.NOT. srcv(jpr_otx1) % laction) ALLOCATE(frcv(jpr_otx1) % z3(jpi, jpj, srcv(jpr_otx1) % nct))
    IF (.NOT. srcv(jpr_oty1) % laction) ALLOCATE(frcv(jpr_oty1) % z3(jpi, jpj, srcv(jpr_oty1) % nct))
    IF (k_ice /= 0) THEN
      IF (.NOT. srcv(jpr_itx1) % laction) ALLOCATE(frcv(jpr_itx1) % z3(jpi, jpj, srcv(jpr_itx1) % nct))
      IF (.NOT. srcv(jpr_ity1) % laction) ALLOCATE(frcv(jpr_ity1) % z3(jpi, jpj, srcv(jpr_ity1) % nct))
    END IF
    ssnd(:) % laction = .FALSE.
    ssnd(:) % clgrid = 'T'
    ssnd(:) % nsgn = 1.
    ssnd(:) % nct = 1
    ssnd(jps_toce) % clname = 'O_SSTSST'
    ssnd(jps_tice) % clname = 'O_TepIce'
    ssnd(jps_ttilyr) % clname = 'O_TtiLyr'
    ssnd(jps_tmix) % clname = 'O_TepMix'
    SELECT CASE (TRIM(sn_snd_temp % cldes))
    CASE ('none')
    CASE ('oce only')
      ssnd(jps_toce) % laction = .TRUE.
    CASE ('oce and ice', 'weighted oce and ice', 'oce and weighted ice')
      ssnd((/jps_toce, jps_tice/)) % laction = .TRUE.
      IF (TRIM(sn_snd_temp % clcat) == 'yes') ssnd(jps_tice) % nct = nn_cats_cpl
    CASE ('mixed oce-ice')
      ssnd(jps_tmix) % laction = .TRUE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_temp%cldes')
    END SELECT
    ssnd(jps_albice) % clname = 'O_AlbIce'
    ssnd(jps_albmix) % clname = 'O_AlbMix'
    SELECT CASE (TRIM(sn_snd_alb % cldes))
    CASE ('none')
    CASE ('ice', 'weighted ice')
      ssnd(jps_albice) % laction = .TRUE.
    CASE ('mixed oce-ice')
      ssnd(jps_albmix) % laction = .TRUE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_alb%cldes')
    END SELECT
    IF (TRIM(sn_snd_alb % cldes) == 'mixed oce-ice' .OR. TRIM(sn_rcv_qsr % cldes) == 'mixed oce-ice') THEN
      CALL oce_alb(zaos, zacs)
      alb_oce_mix(:, :) = (zacs(:, :) + zaos(:, :)) * 0.5
    END IF
    ssnd(jps_fice) % clname = 'OIceFrc'
    ssnd(jps_ficet) % clname = 'OIceFrcT'
    ssnd(jps_hice) % clname = 'OIceTck'
    ssnd(jps_a_p) % clname = 'OPndFrc'
    ssnd(jps_ht_p) % clname = 'OPndTck'
    ssnd(jps_hsnw) % clname = 'OSnwTck'
    ssnd(jps_fice1) % clname = 'OIceFrd'
    IF (k_ice /= 0) THEN
      ssnd(jps_fice) % laction = .TRUE.
      ssnd(jps_fice1) % laction = .TRUE.
      IF (TRIM(sn_snd_thick % clcat) == 'yes') ssnd(jps_fice) % nct = nn_cats_cpl
      IF (TRIM(sn_snd_thick1 % clcat) == 'yes') ssnd(jps_fice1) % nct = nn_cats_cpl
    END IF
    IF (TRIM(sn_snd_ifrac % cldes) == 'coupled') ssnd(jps_ficet) % laction = .TRUE.
    SELECT CASE (TRIM(sn_snd_thick % cldes))
    CASE ('none')
    CASE ('ice and snow')
      ssnd(jps_hice : jps_hsnw) % laction = .TRUE.
      IF (TRIM(sn_snd_thick % clcat) == 'yes') THEN
        ssnd(jps_hice : jps_hsnw) % nct = nn_cats_cpl
      END IF
    CASE ('weighted ice and snow')
      ssnd(jps_hice : jps_hsnw) % laction = .TRUE.
      IF (TRIM(sn_snd_thick % clcat) == 'yes') ssnd(jps_hice : jps_hsnw) % nct = nn_cats_cpl
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_thick%cldes')
    END SELECT
    ssnd(jps_a_p) % clname = 'OPndFrc'
    ssnd(jps_ht_p) % clname = 'OPndTck'
    SELECT CASE (TRIM(sn_snd_mpnd % cldes))
    CASE ('none')
      ssnd(jps_a_p) % laction = .FALSE.
      ssnd(jps_ht_p) % laction = .FALSE.
    CASE ('ice only')
      ssnd(jps_a_p) % laction = .TRUE.
      ssnd(jps_ht_p) % laction = .TRUE.
      IF (TRIM(sn_snd_mpnd % clcat) == 'yes') THEN
        ssnd(jps_a_p) % nct = nn_cats_cpl
        ssnd(jps_ht_p) % nct = nn_cats_cpl
      ELSE
        IF (nn_cats_cpl > 1) THEN
          CALL ctl_stop('sbc_cpl_init: use weighted ice option for sn_snd_mpnd%cldes if not exchanging category fields')
        END IF
      END IF
    CASE ('weighted ice')
      ssnd(jps_a_p) % laction = .TRUE.
      ssnd(jps_ht_p) % laction = .TRUE.
      IF (TRIM(sn_snd_mpnd % clcat) == 'yes') THEN
        ssnd(jps_a_p) % nct = nn_cats_cpl
        ssnd(jps_ht_p) % nct = nn_cats_cpl
      END IF
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_mpnd%cldes; ' // sn_snd_mpnd % cldes)
    END SELECT
    ssnd(jps_ocx1) % clname = 'O_OCurx1'
    ssnd(jps_ivx1) % clname = 'O_IVelx1'
    ssnd(jps_ocy1) % clname = 'O_OCury1'
    ssnd(jps_ivy1) % clname = 'O_IVely1'
    ssnd(jps_ocz1) % clname = 'O_OCurz1'
    ssnd(jps_ivz1) % clname = 'O_IVelz1'
    ssnd(jps_ocxw) % clname = 'O_OCurxw'
    ssnd(jps_ocyw) % clname = 'O_OCuryw'
    ssnd(jps_ocx1 : jps_ivz1) % nsgn = - 1.
    IF (sn_snd_crt % clvgrd == 'U,V') THEN
      ssnd(jps_ocx1) % clgrid = 'U'
      ssnd(jps_ocy1) % clgrid = 'V'
    ELSE IF (sn_snd_crt % clvgrd /= 'T') THEN
      CALL ctl_stop('sn_snd_crt%clvgrd must be equal to T')
      ssnd(jps_ocx1 : jps_ivz1) % clgrid = 'T'
    END IF
    ssnd(jps_ocx1 : jps_ivz1) % laction = .TRUE.
    IF (TRIM(sn_snd_crt % clvref) == 'spherical') ssnd((/jps_ocz1, jps_ivz1/)) % laction = .FALSE.
    IF (TRIM(sn_snd_crt % clvor) == 'eastward-northward') ssnd(jps_ocx1 : jps_ivz1) % nsgn = 1.
    SELECT CASE (TRIM(sn_snd_crt % cldes))
    CASE ('none')
      ssnd(jps_ocx1 : jps_ivz1) % laction = .FALSE.
    CASE ('oce only')
      ssnd(jps_ivx1 : jps_ivz1) % laction = .FALSE.
    CASE ('weighted oce and ice')
    CASE ('mixed oce-ice')
      ssnd(jps_ivx1 : jps_ivz1) % laction = .FALSE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_crt%cldes')
    END SELECT
    ssnd(jps_ocxw : jps_ocyw) % nsgn = - 1.
    IF (sn_snd_crtw % clvgrd == 'U,V') THEN
      ssnd(jps_ocxw) % clgrid = 'U'
      ssnd(jps_ocyw) % clgrid = 'V'
    ELSE IF (sn_snd_crtw % clvgrd /= 'T') THEN
      CALL ctl_stop('sn_snd_crtw%clvgrd must be equal to T')
    END IF
    IF (TRIM(sn_snd_crtw % clvor) == 'eastward-northward') ssnd(jps_ocxw : jps_ocyw) % nsgn = 1.
    SELECT CASE (TRIM(sn_snd_crtw % cldes))
    CASE ('none')
      ssnd(jps_ocxw : jps_ocyw) % laction = .FALSE.
    CASE ('oce only')
      ssnd(jps_ocxw : jps_ocyw) % laction = .TRUE.
    CASE ('weighted oce and ice')
    CASE ('mixed oce-ice')
      ssnd(jps_ivx1 : jps_ivz1) % laction = .FALSE.
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_crtw%cldes')
    END SELECT
    ssnd(jps_co2) % clname = 'O_CO2FLX'
    IF (TRIM(sn_snd_co2 % cldes) == 'coupled') ssnd(jps_co2) % laction = .TRUE.
    ssnd(jps_sstfrz) % clname = 'O_SSTFrz'
    IF (TRIM(sn_snd_sstfrz % cldes) == 'coupled') ssnd(jps_sstfrz) % laction = .TRUE.
    ssnd(jps_ttilyr) % clname = 'O_TtiLyr'
    SELECT CASE (TRIM(sn_snd_ttilyr % cldes))
    CASE ('none')
      ssnd(jps_ttilyr) % laction = .FALSE.
    CASE ('ice only')
      ssnd(jps_ttilyr) % laction = .TRUE.
      IF (TRIM(sn_snd_ttilyr % clcat) == 'yes') THEN
        ssnd(jps_ttilyr) % nct = nn_cats_cpl
      ELSE
        IF (nn_cats_cpl > 1) THEN
          CALL ctl_stop('sbc_cpl_init: use weighted ice option for sn_snd_ttilyr%cldes if not exchanging category fields')
        END IF
      END IF
    CASE ('weighted ice')
      ssnd(jps_ttilyr) % laction = .TRUE.
      IF (TRIM(sn_snd_ttilyr % clcat) == 'yes') ssnd(jps_ttilyr) % nct = nn_cats_cpl
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_ttilyr%cldes;' // sn_snd_ttilyr % cldes)
    END SELECT
    ssnd(jps_kice) % clname = 'OIceKn'
    SELECT CASE (TRIM(sn_snd_cond % cldes))
    CASE ('none')
      ssnd(jps_kice) % laction = .FALSE.
    CASE ('ice only')
      ssnd(jps_kice) % laction = .TRUE.
      IF (TRIM(sn_snd_cond % clcat) == 'yes') THEN
        ssnd(jps_kice) % nct = nn_cats_cpl
      ELSE
        IF (nn_cats_cpl > 1) THEN
          CALL ctl_stop('sbc_cpl_init: use weighted ice option for sn_snd_cond%cldes if not exchanging category fields')
        END IF
      END IF
    CASE ('weighted ice')
      ssnd(jps_kice) % laction = .TRUE.
      IF (TRIM(sn_snd_cond % clcat) == 'yes') ssnd(jps_kice) % nct = nn_cats_cpl
    CASE DEFAULT
      CALL ctl_stop('sbc_cpl_init: wrong definition of sn_snd_cond%cldes;' // sn_snd_cond % cldes)
    END SELECT
    ssnd(jps_wlev) % clname = 'O_Wlevel'
    IF (TRIM(sn_snd_wlev % cldes) == 'coupled') ssnd(jps_wlev) % laction = .TRUE.
    ssnd(jps_ssh) % clname = 'O_SSHght'
    ssnd(jps_soce) % clname = 'O_SSSal'
    ssnd(jps_e3t1st) % clname = 'O_E3T1st'
    ssnd(jps_fraqsr) % clname = 'O_FraQsr'
    IF (nn_components == jp_iam_opa) THEN
      ssnd(:) % laction = .FALSE.
      ssnd((/jps_toce, jps_soce, jps_ssh, jps_fraqsr, jps_ocx1, jps_ocy1/)) % laction = .TRUE.
      ssnd(jps_e3t1st) % laction = .NOT. ln_linssh
      ssnd(jps_ocx1) % clgrid = 'U'
      ssnd(jps_ocy1) % clgrid = 'V'
      sn_snd_crt % clvgrd = 'U,V'
      sn_snd_crt % clvor = 'local grid'
      sn_snd_crt % clvref = 'spherical'
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '  sent fields to SAS component '
        WRITE(numout, *) '               sea surface temperature (T before, Celsius) '
        WRITE(numout, *) '               sea surface salinity '
        WRITE(numout, *) '               surface currents U,V on local grid and spherical coordinates'
        WRITE(numout, *) '               sea surface height '
        WRITE(numout, *) '               thickness of first ocean T level '
        WRITE(numout, *) '               fraction of solar net radiation absorbed in the first ocean level'
        WRITE(numout, *)
      END IF
    END IF
    ssnd(jps_sflx) % clname = 'I_SFLX'
    ssnd(jps_fice2) % clname = 'IIceFrc'
    ssnd(jps_qsroce) % clname = 'I_QsrOce'
    ssnd(jps_qnsoce) % clname = 'I_QnsOce'
    ssnd(jps_oemp) % clname = 'IOEvaMPr'
    ssnd(jps_otx1) % clname = 'I_OTaux1'
    ssnd(jps_oty1) % clname = 'I_OTauy1'
    ssnd(jps_rnf) % clname = 'I_Runoff'
    ssnd(jps_taum) % clname = 'I_TauMod'
    IF (nn_components == jp_iam_sas) THEN
      IF (.NOT. ln_cpl) ssnd(:) % laction = .FALSE.
      ssnd((/jps_qsroce, jps_qnsoce, jps_oemp, jps_fice2, jps_sflx, jps_otx1, jps_oty1, jps_taum/)) % laction = .TRUE.
      DO jn = 1, jpsnd
        IF (ssnd(jn) % clname(1 : 1) == "O") ssnd(jn) % clname = "S" // ssnd(jn) % clname(2 : LEN(ssnd(jn) % clname))
      END DO
      IF (lwp) THEN
        WRITE(numout, *)
        IF (.NOT. ln_cpl) THEN
          WRITE(numout, *) '  sent fields to OPA component '
        ELSE
          WRITE(numout, *) '  Additional sent fields to OPA component : '
        END IF
        WRITE(numout, *) '                  ice cover '
        WRITE(numout, *) '                  oce only EMP  '
        WRITE(numout, *) '                  salt flux  '
        WRITE(numout, *) '                  mixed oce-ice solar flux  '
        WRITE(numout, *) '                  mixed oce-ice non solar flux  '
        WRITE(numout, *) '                  wind stress U,V components'
        WRITE(numout, *) '                  wind stress module'
      END IF
    END IF
    CALL cpl_define(jprcv, jpsnd, nn_cplmodel)
    IF (ln_usecplmask) THEN
      xcplmask(:, :, :) = 0.
      CALL iom_open('cplmask', inum)
      CALL iom_get(inum, jpdom_unknown, 'cplmask', xcplmask(1 : nlci, 1 : nlcj, 1 : nn_cplmodel), kstart = (/mig(1), mjg(1), 1/), kcount = (/nlci, nlcj, nn_cplmodel/))
      CALL iom_close(inum)
    ELSE
      xcplmask(:, :, :) = 1.
    END IF
    xcplmask(:, :, 0) = 1. - SUM(xcplmask(:, :, 1 : nn_cplmodel), dim = 3)
    ncpl_qsr_freq = cpl_freq('O_QsrOce') + cpl_freq('O_QsrMix') + cpl_freq('I_QsrOce') + cpl_freq('I_QsrMix')
    IF (ln_dm2dc .AND. ln_cpl .AND. ncpl_qsr_freq /= 86400) CALL ctl_stop('sbc_cpl_init: diurnal cycle reconstruction (ln_dm2dc) needs daily couping for solar radiation')
    IF (ln_dm2dc .AND. ln_cpl) ncpl_qsr_freq = 86400 / ncpl_qsr_freq
  END SUBROUTINE sbc_cpl_init
  SUBROUTINE sbc_cpl_rcv(kt, k_fsbc, k_ice)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE zdf_oce, ONLY: ln_zdfswm
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: k_fsbc
    INTEGER, INTENT(IN) :: k_ice
    LOGICAL :: llnewtx, llnewtau
    INTEGER :: ji, jj, jn
    INTEGER :: isec
    REAL(KIND = wp) :: zcumulneg, zcumulpos
    REAL(KIND = wp) :: zcoef
    REAL(KIND = wp) :: zrhoa = 1.22
    REAL(KIND = wp) :: zcdrag = 1.5E-3
    REAL(KIND = wp) :: zzx, zzy
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztx, zty, zmsk, zemp, zqns, zqsr
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_cpl_rcv', 'r0', 0, 0)
    IF (ln_mixcpl) zmsk(:, :) = 1. - xcplmask(:, :, 0)
    isec = (kt - nit000) * NINT(rdt)
    DO jn = 1, jprcv
      IF (srcv(jn) % laction) CALL cpl_rcv(jn, isec, frcv(jn) % z3, xcplmask(:, :, 1 : nn_cplmodel), nrcvinfo(jn))
    END DO
    IF (srcv(jpr_otx1) % laction) THEN
      IF (nrcvinfo(jpr_otx1) == OASIS_Rcv) THEN
        IF (TRIM(sn_rcv_tau % clvref) == 'cartesian') THEN
          CALL geo2oce(frcv(jpr_otx1) % z3(:, :, 1), frcv(jpr_oty1) % z3(:, :, 1), frcv(jpr_otz1) % z3(:, :, 1), srcv(jpr_otx1) % clgrid, ztx, zty)
          frcv(jpr_otx1) % z3(:, :, 1) = ztx(:, :)
          frcv(jpr_oty1) % z3(:, :, 1) = zty(:, :)
          IF (srcv(jpr_otx2) % laction) THEN
            CALL geo2oce(frcv(jpr_otx2) % z3(:, :, 1), frcv(jpr_oty2) % z3(:, :, 1), frcv(jpr_otz2) % z3(:, :, 1), srcv(jpr_otx2) % clgrid, ztx, zty)
            frcv(jpr_otx2) % z3(:, :, 1) = ztx(:, :)
            frcv(jpr_oty2) % z3(:, :, 1) = zty(:, :)
          END IF
        END IF
        IF (TRIM(sn_rcv_tau % clvor) == 'eastward-northward') THEN
          CALL rot_rep(frcv(jpr_otx1) % z3(:, :, 1), frcv(jpr_oty1) % z3(:, :, 1), srcv(jpr_otx1) % clgrid, 'en->i', ztx)
          IF (srcv(jpr_otx2) % laction) THEN
            CALL rot_rep(frcv(jpr_otx2) % z3(:, :, 1), frcv(jpr_oty2) % z3(:, :, 1), srcv(jpr_otx2) % clgrid, 'en->j', zty)
          ELSE
            CALL rot_rep(frcv(jpr_otx1) % z3(:, :, 1), frcv(jpr_oty1) % z3(:, :, 1), srcv(jpr_otx1) % clgrid, 'en->j', zty)
          END IF
          frcv(jpr_otx1) % z3(:, :, 1) = ztx(:, :)
          frcv(jpr_oty1) % z3(:, :, 1) = zty(:, :)
        END IF
        IF (srcv(jpr_otx1) % clgrid == 'T') THEN
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              frcv(jpr_otx1) % z3(ji, jj, 1) = 0.5 * (frcv(jpr_otx1) % z3(ji + 1, jj, 1) + frcv(jpr_otx1) % z3(ji, jj, 1))
              frcv(jpr_oty1) % z3(ji, jj, 1) = 0.5 * (frcv(jpr_oty1) % z3(ji, jj + 1, 1) + frcv(jpr_oty1) % z3(ji, jj, 1))
            END DO
          END DO
          CALL lbc_lnk_multi('sbccpl', frcv(jpr_otx1) % z3(:, :, 1), 'U', - 1., frcv(jpr_oty1) % z3(:, :, 1), 'V', - 1.)
        END IF
        llnewtx = .TRUE.
      ELSE
        llnewtx = .FALSE.
      END IF
    ELSE
      frcv(jpr_otx1) % z3(:, :, 1) = 0.E0
      frcv(jpr_oty1) % z3(:, :, 1) = 0.E0
      llnewtx = .TRUE.
    END IF
    IF (.NOT. srcv(jpr_taum) % laction) THEN
      IF (llnewtx) THEN
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zzx = frcv(jpr_otx1) % z3(ji - 1, jj, 1) + frcv(jpr_otx1) % z3(ji, jj, 1)
            zzy = frcv(jpr_oty1) % z3(ji, jj - 1, 1) + frcv(jpr_oty1) % z3(ji, jj, 1)
            frcv(jpr_taum) % z3(ji, jj, 1) = 0.5 * SQRT(zzx * zzx + zzy * zzy)
          END DO
        END DO
        CALL lbc_lnk('sbccpl', frcv(jpr_taum) % z3(:, :, 1), 'T', 1.)
        llnewtau = .TRUE.
      ELSE
        llnewtau = .FALSE.
      END IF
    ELSE
      llnewtau = nrcvinfo(jpr_taum) == OASIS_Rcv
      IF (llnewtau) THEN
        frcv(jpr_taum) % z3(:, :, 1) = MAX(0._wp, frcv(jpr_taum) % z3(:, :, 1))
      END IF
    END IF
    IF (.NOT. srcv(jpr_w10m) % laction) THEN
      IF (llnewtau) THEN
        zcoef = 1. / (zrhoa * zcdrag)
        DO jj = 1, jpj
          DO ji = 1, jpi
            frcv(jpr_w10m) % z3(ji, jj, 1) = SQRT(frcv(jpr_taum) % z3(ji, jj, 1) * zcoef)
          END DO
        END DO
      END IF
    END IF
    IF (MOD(kt - 1, k_fsbc) == 0) THEN
      IF (ln_mixcpl) THEN
        utau(:, :) = utau(:, :) * xcplmask(:, :, 0) + frcv(jpr_otx1) % z3(:, :, 1) * zmsk(:, :)
        vtau(:, :) = vtau(:, :) * xcplmask(:, :, 0) + frcv(jpr_oty1) % z3(:, :, 1) * zmsk(:, :)
        taum(:, :) = taum(:, :) * xcplmask(:, :, 0) + frcv(jpr_taum) % z3(:, :, 1) * zmsk(:, :)
        wndm(:, :) = wndm(:, :) * xcplmask(:, :, 0) + frcv(jpr_w10m) % z3(:, :, 1) * zmsk(:, :)
      ELSE
        utau(:, :) = frcv(jpr_otx1) % z3(:, :, 1)
        vtau(:, :) = frcv(jpr_oty1) % z3(:, :, 1)
        taum(:, :) = frcv(jpr_taum) % z3(:, :, 1)
        wndm(:, :) = frcv(jpr_w10m) % z3(:, :, 1)
      END IF
      CALL iom_put("taum_oce", taum)
    END IF
    IF (srcv(jpr_co2) % laction) atm_co2(:, :) = frcv(jpr_co2) % z3(:, :, 1)
    IF (srcv(jpr_mslp) % laction) THEN
      IF (kt /= nit000) ssh_ibb(:, :) = ssh_ib(:, :)
      r1_grau = 1.E0 / (grav * rau0)
      ssh_ib(:, :) = - (frcv(jpr_mslp) % z3(:, :, 1) - rpref) * r1_grau
      apr(:, :) = frcv(jpr_mslp) % z3(:, :, 1)
      IF (kt == nit000) ssh_ibb(:, :) = ssh_ib(:, :)
    END IF
    IF (ln_sdw) THEN
      IF (srcv(jpr_sdrftx) % laction) ut0sd(:, :) = frcv(jpr_sdrftx) % z3(:, :, 1)
      IF (srcv(jpr_sdrfty) % laction) vt0sd(:, :) = frcv(jpr_sdrfty) % z3(:, :, 1)
      IF (srcv(jpr_wper) % laction) wmp(:, :) = frcv(jpr_wper) % z3(:, :, 1)
      IF (srcv(jpr_hsig) % laction) hsw(:, :) = frcv(jpr_hsig) % z3(:, :, 1)
      IF (srcv(jpr_wfreq) % laction) wfreq(:, :) = frcv(jpr_wfreq) % z3(:, :, 1)
      IF (srcv(jpr_wnum) % laction .AND. ln_zdfswm) wnum(:, :) = frcv(jpr_wnum) % z3(:, :, 1)
      IF (srcv(jpr_sdrftx) % laction .OR. srcv(jpr_sdrfty) % laction .OR. srcv(jpr_wper) % laction .OR. srcv(jpr_hsig) % laction .OR. srcv(jpr_wfreq) % laction) THEN
        CALL sbc_stokes
      END IF
    END IF
    IF (srcv(jpr_tauwoc) % laction .AND. ln_tauwoc) tauoc_wave(:, :) = frcv(jpr_tauwoc) % z3(:, :, 1)
    IF (srcv(jpr_tauwx) % laction .AND. srcv(jpr_tauwy) % laction .AND. ln_tauw) THEN
      tauw_x(:, :) = frcv(jpr_tauwx) % z3(:, :, 1)
      tauw_y(:, :) = frcv(jpr_tauwy) % z3(:, :, 1)
    END IF
    IF (srcv(jpr_wdrag) % laction .AND. ln_cdgw) cdn_wave(:, :) = frcv(jpr_wdrag) % z3(:, :, 1)
    IF (srcv(jpr_soce) % laction) THEN
      sss_m(:, :) = frcv(jpr_soce) % z3(:, :, 1)
      CALL iom_put('sss_m', sss_m)
    END IF
    IF (srcv(jpr_toce) % laction) THEN
      sst_m(:, :) = frcv(jpr_toce) % z3(:, :, 1)
      IF (srcv(jpr_soce) % laction .AND. l_useCT) THEN
        sst_m(:, :) = eos_pt_from_ct(sst_m(:, :), sss_m(:, :))
      END IF
    END IF
    IF (srcv(jpr_ssh) % laction) THEN
      ssh_m(:, :) = frcv(jpr_ssh) % z3(:, :, 1)
      CALL iom_put('ssh_m', ssh_m)
    END IF
    IF (srcv(jpr_ocx1) % laction) THEN
      ssu_m(:, :) = frcv(jpr_ocx1) % z3(:, :, 1)
      ub(:, :, 1) = ssu_m(:, :)
      un(:, :, 1) = ssu_m(:, :)
      CALL iom_put('ssu_m', ssu_m)
    END IF
    IF (srcv(jpr_ocy1) % laction) THEN
      ssv_m(:, :) = frcv(jpr_ocy1) % z3(:, :, 1)
      vb(:, :, 1) = ssv_m(:, :)
      vn(:, :, 1) = ssv_m(:, :)
      CALL iom_put('ssv_m', ssv_m)
    END IF
    IF (srcv(jpr_e3t1st) % laction) THEN
      e3t_m(:, :) = frcv(jpr_e3t1st) % z3(:, :, 1)
      CALL iom_put('e3t_m', e3t_m(:, :))
    END IF
    IF (srcv(jpr_fraqsr) % laction) THEN
      frq_m(:, :) = frcv(jpr_fraqsr) % z3(:, :, 1)
      CALL iom_put('frq_m', frq_m)
    END IF
    IF (k_ice <= 1 .AND. MOD(kt - 1, k_fsbc) == 0) THEN
      IF (srcv(jpr_oemp) % laction .OR. srcv(jpr_rain) % laction) THEN
        SELECT CASE (TRIM(sn_rcv_emp % cldes))
        CASE ('conservative')
          zemp(:, :) = frcv(jpr_tevp) % z3(:, :, 1) - (frcv(jpr_rain) % z3(:, :, 1) + frcv(jpr_snow) % z3(:, :, 1))
        CASE ('oce only', 'oce and ice')
          zemp(:, :) = frcv(jpr_oemp) % z3(:, :, 1)
        CASE DEFAULT
          CALL ctl_stop('sbc_cpl_rcv: wrong definition of sn_rcv_emp%cldes')
        END SELECT
      ELSE
        zemp(:, :) = 0._wp
      END IF
      IF (srcv(jpr_rnf) % laction) rnf(:, :) = frcv(jpr_rnf) % z3(:, :, 1)
      IF (srcv(jpr_cal) % laction) zemp(:, :) = zemp(:, :) - frcv(jpr_cal) % z3(:, :, 1)
      IF (srcv(jpr_icb) % laction) THEN
        fwficb(:, :) = frcv(jpr_icb) % z3(:, :, 1)
        rnf(:, :) = rnf(:, :) + fwficb(:, :)
      END IF
      IF (srcv(jpr_isf) % laction) fwfisf(:, :) = - frcv(jpr_isf) % z3(:, :, 1)
      IF (ln_mixcpl) THEN
        emp(:, :) = emp(:, :) * xcplmask(:, :, 0) + zemp(:, :) * zmsk(:, :)
      ELSE
        emp(:, :) = zemp(:, :)
      END IF
      IF (srcv(jpr_qnsoce) % laction) THEN
        zqns(:, :) = frcv(jpr_qnsoce) % z3(:, :, 1)
      ELSE IF (srcv(jpr_qnsmix) % laction) THEN
        zqns(:, :) = frcv(jpr_qnsmix) % z3(:, :, 1)
      ELSE
        zqns(:, :) = 0._wp
      END IF
      IF (nn_components /= jp_iam_opa) THEN
        zqns(:, :) = zqns(:, :) - zemp(:, :) * sst_m(:, :) * rcp
        IF (srcv(jpr_snow) % laction) THEN
          zqns(:, :) = zqns(:, :) - frcv(jpr_snow) % z3(:, :, 1) * rLfus
        END IF
      END IF
      IF (srcv(jpr_icb) % laction) zqns(:, :) = zqns(:, :) - frcv(jpr_icb) % z3(:, :, 1) * rLfus
      IF (ln_mixcpl) THEN
        qns(:, :) = qns(:, :) * xcplmask(:, :, 0) + zqns(:, :) * zmsk(:, :)
      ELSE
        qns(:, :) = zqns(:, :)
      END IF
      IF (srcv(jpr_qsroce) % laction) THEN
        zqsr(:, :) = frcv(jpr_qsroce) % z3(:, :, 1)
      ELSE IF (srcv(jpr_qsrmix) % laction) THEN
        zqsr(:, :) = frcv(jpr_qsrmix) % z3(:, :, 1)
      ELSE
        zqsr(:, :) = 0._wp
      END IF
      IF (ln_dm2dc .AND. ln_cpl) zqsr(:, :) = sbc_dcy(zqsr)
      IF (ln_mixcpl) THEN
        qsr(:, :) = qsr(:, :) * xcplmask(:, :, 0) + zqsr(:, :) * zmsk(:, :)
      ELSE
        qsr(:, :) = zqsr(:, :)
      END IF
      IF (srcv(jpr_sflx) % laction) sfx(:, :) = frcv(jpr_sflx) % z3(:, :, 1)
      IF (srcv(jpr_fice) % laction) fr_i(:, :) = frcv(jpr_fice) % z3(:, :, 1)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sbc_cpl_rcv
  SUBROUTINE sbc_cpl_ice_tau(p_taui, p_tauj)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:, :) :: p_taui
    REAL(KIND = wp), INTENT(OUT), DIMENSION(:, :) :: p_tauj
    INTEGER :: ji, jj
    INTEGER :: itx
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztx, zty
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_cpl_ice_tau', 'r0', 0, 0)
    IF (srcv(jpr_itx1) % laction) THEN
      itx = jpr_itx1
    ELSE
      itx = jpr_otx1
    END IF
    IF (nrcvinfo(itx) == OASIS_Rcv) THEN
      IF (srcv(jpr_itx1) % laction) THEN
        IF (TRIM(sn_rcv_tau % clvref) == 'cartesian') THEN
          CALL geo2oce(frcv(jpr_itx1) % z3(:, :, 1), frcv(jpr_ity1) % z3(:, :, 1), frcv(jpr_itz1) % z3(:, :, 1), srcv(jpr_itx1) % clgrid, ztx, zty)
          frcv(jpr_itx1) % z3(:, :, 1) = ztx(:, :)
          frcv(jpr_ity1) % z3(:, :, 1) = zty(:, :)
          IF (srcv(jpr_itx2) % laction) THEN
            CALL geo2oce(frcv(jpr_itx2) % z3(:, :, 1), frcv(jpr_ity2) % z3(:, :, 1), frcv(jpr_itz2) % z3(:, :, 1), srcv(jpr_itx2) % clgrid, ztx, zty)
            frcv(jpr_itx2) % z3(:, :, 1) = ztx(:, :)
            frcv(jpr_ity2) % z3(:, :, 1) = zty(:, :)
          END IF
        END IF
        IF (TRIM(sn_rcv_tau % clvor) == 'eastward-northward') THEN
          CALL rot_rep(frcv(jpr_itx1) % z3(:, :, 1), frcv(jpr_ity1) % z3(:, :, 1), srcv(jpr_itx1) % clgrid, 'en->i', ztx)
          IF (srcv(jpr_itx2) % laction) THEN
            CALL rot_rep(frcv(jpr_itx2) % z3(:, :, 1), frcv(jpr_ity2) % z3(:, :, 1), srcv(jpr_itx2) % clgrid, 'en->j', zty)
          ELSE
            CALL rot_rep(frcv(jpr_itx1) % z3(:, :, 1), frcv(jpr_ity1) % z3(:, :, 1), srcv(jpr_itx1) % clgrid, 'en->j', zty)
          END IF
          frcv(jpr_itx1) % z3(:, :, 1) = ztx(:, :)
          frcv(jpr_ity1) % z3(:, :, 1) = zty(:, :)
        END IF
      ELSE
        frcv(jpr_itx1) % z3(:, :, 1) = frcv(jpr_otx1) % z3(:, :, 1)
        frcv(jpr_ity1) % z3(:, :, 1) = frcv(jpr_oty1) % z3(:, :, 1)
      END IF
      SELECT CASE (srcv(jpr_itx1) % clgrid)
      CASE ('U')
        p_taui(:, :) = frcv(jpr_itx1) % z3(:, :, 1)
        p_tauj(:, :) = frcv(jpr_ity1) % z3(:, :, 1)
      CASE ('F')
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_taui(ji, jj) = 0.5 * (frcv(jpr_itx1) % z3(ji, jj, 1) + frcv(jpr_itx1) % z3(ji, jj - 1, 1))
            p_tauj(ji, jj) = 0.5 * (frcv(jpr_ity1) % z3(ji, jj, 1) + frcv(jpr_ity1) % z3(ji - 1, jj, 1))
          END DO
        END DO
      CASE ('T')
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_taui(ji, jj) = 0.5 * (frcv(jpr_itx1) % z3(ji + 1, jj, 1) + frcv(jpr_itx1) % z3(ji, jj, 1))
            p_tauj(ji, jj) = 0.5 * (frcv(jpr_ity1) % z3(ji, jj + 1, 1) + frcv(jpr_ity1) % z3(ji, jj, 1))
          END DO
        END DO
      CASE ('I')
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_taui(ji, jj) = 0.5 * (frcv(jpr_itx1) % z3(ji + 1, jj + 1, 1) + frcv(jpr_itx1) % z3(ji + 1, jj, 1))
            p_tauj(ji, jj) = 0.5 * (frcv(jpr_ity1) % z3(ji + 1, jj + 1, 1) + frcv(jpr_ity1) % z3(ji, jj + 1, 1))
          END DO
        END DO
      END SELECT
      IF (srcv(jpr_itx1) % clgrid /= 'U') THEN
        CALL lbc_lnk_multi('sbccpl', p_taui, 'U', - 1., p_tauj, 'V', - 1.)
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sbc_cpl_ice_tau
  SUBROUTINE sbc_cpl_ice_flx(picefr, palbi, psst, pist, phs, phi)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :) :: picefr
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :), OPTIONAL :: palbi
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :), OPTIONAL :: psst
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :), OPTIONAL :: pist
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :), OPTIONAL :: phs
    REAL(KIND = wp), INTENT(IN), DIMENSION(:, :, :), OPTIONAL :: phi
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: ztri
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zcptn, zcptrain, zcptsnw, ziceld, zmsk, zsnw
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zemp_tot, zemp_ice, zemp_oce, ztprecip, zsprecip, zevap_oce, zdevap_ice
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zqns_tot, zqns_oce, zqsr_tot, zqsr_oce, zqprec_ice, zqemp_oce, zqemp_ice
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zqns_ice, zqsr_ice, zdqns_ice, zqevap_ice, zevap_ice
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_cpl_ice_flx', 'r0', 0, 0)
    IF (ln_mixcpl) zmsk(:, :) = 1. - xcplmask(:, :, 0)
    ziceld(:, :) = 1._wp - picefr(:, :)
    zcptn(:, :) = rcp * sst_m(:, :)
    SELECT CASE (TRIM(sn_rcv_emp % cldes))
    CASE ('conservative')
      zsprecip(:, :) = frcv(jpr_snow) % z3(:, :, 1)
      ztprecip(:, :) = frcv(jpr_rain) % z3(:, :, 1) + zsprecip(:, :)
      zemp_tot(:, :) = frcv(jpr_tevp) % z3(:, :, 1) - ztprecip(:, :)
      zemp_ice(:, :) = (frcv(jpr_ievp) % z3(:, :, 1) - frcv(jpr_snow) % z3(:, :, 1)) * picefr(:, :)
    CASE ('oce and ice')
      zemp_tot(:, :) = ziceld(:, :) * frcv(jpr_oemp) % z3(:, :, 1) + picefr(:, :) * frcv(jpr_sbpr) % z3(:, :, 1)
      zemp_ice(:, :) = frcv(jpr_semp) % z3(:, :, 1) * picefr(:, :)
      zsprecip(:, :) = frcv(jpr_ievp) % z3(:, :, 1) - frcv(jpr_semp) % z3(:, :, 1)
      ztprecip(:, :) = frcv(jpr_semp) % z3(:, :, 1) - frcv(jpr_sbpr) % z3(:, :, 1) + zsprecip(:, :)
    CASE ('none')
      CALL ctl_stop('STOP', 'sbccpl/sbc_cpl_ice_flx: some fields are not defined. Change sn_rcv_emp value in namelist namsbc_cpl')
    END SELECT
    zsnw(:, :) = picefr(:, :)
    IF (srcv(jpr_rnf) % laction) THEN
      rnf(:, :) = frcv(jpr_rnf) % z3(:, :, 1)
    END IF
    IF (srcv(jpr_cal) % laction) THEN
      zemp_tot(:, :) = zemp_tot(:, :) - frcv(jpr_cal) % z3(:, :, 1)
    END IF
    IF (srcv(jpr_icb) % laction) THEN
      fwficb(:, :) = frcv(jpr_icb) % z3(:, :, 1)
      rnf(:, :) = rnf(:, :) + fwficb(:, :)
    END IF
    IF (srcv(jpr_isf) % laction) THEN
      fwfisf(:, :) = - frcv(jpr_isf) % z3(:, :, 1)
    END IF
    IF (ln_mixcpl) THEN
      emp_tot(:, :) = emp_tot(:, :) * xcplmask(:, :, 0) + zemp_tot(:, :) * zmsk(:, :)
      emp_ice(:, :) = emp_ice(:, :) * xcplmask(:, :, 0) + zemp_ice(:, :) * zmsk(:, :)
      sprecip(:, :) = sprecip(:, :) * xcplmask(:, :, 0) + zsprecip(:, :) * zmsk(:, :)
      tprecip(:, :) = tprecip(:, :) * xcplmask(:, :, 0) + ztprecip(:, :) * zmsk(:, :)
    ELSE
      emp_tot(:, :) = zemp_tot(:, :)
      emp_ice(:, :) = zemp_ice(:, :)
      sprecip(:, :) = zsprecip(:, :)
      tprecip(:, :) = ztprecip(:, :)
    END IF
    IF (srcv(jpr_cal) % laction) CALL iom_put('calving_cea', frcv(jpr_cal) % z3(:, :, 1) * tmask(:, :, 1))
    IF (srcv(jpr_icb) % laction) CALL iom_put('iceberg_cea', frcv(jpr_icb) % z3(:, :, 1) * tmask(:, :, 1))
    IF (iom_use('snowpre')) CALL iom_put('snowpre', sprecip(:, :))
    IF (iom_use('precip')) CALL iom_put('precip', tprecip(:, :))
    IF (iom_use('rain')) CALL iom_put('rain', tprecip(:, :) - sprecip(:, :))
    IF (iom_use('snow_ao_cea')) CALL iom_put('snow_ao_cea', sprecip(:, :) * (1._wp - zsnw(:, :)))
    IF (iom_use('snow_ai_cea')) CALL iom_put('snow_ai_cea', sprecip(:, :) * zsnw(:, :))
    IF (iom_use('subl_ai_cea')) CALL iom_put('subl_ai_cea', frcv(jpr_ievp) % z3(:, :, 1) * picefr(:, :) * tmask(:, :, 1))
    IF (iom_use('evap_ao_cea')) CALL iom_put('evap_ao_cea', (frcv(jpr_tevp) % z3(:, :, 1) - frcv(jpr_ievp) % z3(:, :, 1) * picefr(:, :)) * tmask(:, :, 1))
    SELECT CASE (TRIM(sn_rcv_qns % cldes))
    CASE ('oce only')
      zqns_tot(:, :) = frcv(jpr_qnsoce) % z3(:, :, 1)
    CASE ('conservative')
      zqns_tot(:, :) = frcv(jpr_qnsmix) % z3(:, :, 1)
      IF (TRIM(sn_rcv_qns % clcat) == 'yes') THEN
        zqns_ice(:, :, 1 : jpl) = frcv(jpr_qnsice) % z3(:, :, 1 : jpl)
      ELSE
        DO jl = 1, jpl
          zqns_ice(:, :, jl) = frcv(jpr_qnsice) % z3(:, :, 1)
        END DO
      END IF
    CASE ('oce and ice')
      zqns_tot(:, :) = ziceld(:, :) * frcv(jpr_qnsoce) % z3(:, :, 1)
      IF (TRIM(sn_rcv_qns % clcat) == 'yes') THEN
        DO jl = 1, jpl
          zqns_tot(:, :) = zqns_tot(:, :) + a_i(:, :, jl) * frcv(jpr_qnsice) % z3(:, :, jl)
          zqns_ice(:, :, jl) = frcv(jpr_qnsice) % z3(:, :, jl)
        END DO
      ELSE
        qns_tot(:, :) = qns_tot(:, :) + picefr(:, :) * frcv(jpr_qnsice) % z3(:, :, 1)
        DO jl = 1, jpl
          zqns_tot(:, :) = zqns_tot(:, :) + picefr(:, :) * frcv(jpr_qnsice) % z3(:, :, 1)
          zqns_ice(:, :, jl) = frcv(jpr_qnsice) % z3(:, :, 1)
        END DO
      END IF
    CASE ('mixed oce-ice')
      zqns_tot(:, :) = frcv(jpr_qnsmix) % z3(:, :, 1)
      zqns_ice(:, :, 1) = frcv(jpr_qnsmix) % z3(:, :, 1) + frcv(jpr_dqnsdt) % z3(:, :, 1) * (pist(:, :, 1) - ((rt0 + psst(:, :)) * ziceld(:, :) + pist(:, :, 1) * picefr(:, :)))
    END SELECT
    IF (srcv(jpr_cal) % laction) zqns_tot(:, :) = zqns_tot(:, :) - frcv(jpr_cal) % z3(:, :, 1) * rLfus
    IF (srcv(jpr_icb) % laction) zqns_tot(:, :) = zqns_tot(:, :) - frcv(jpr_icb) % z3(:, :, 1) * rLfus
    zcptsnw(:, :) = zcptn(:, :)
    zcptrain(:, :) = zcptn(:, :)
    zqns_tot(:, :) = zqns_tot(:, :) - (ziceld(:, :) * zsprecip(:, :) * rLfus) - (zemp_tot(:, :) - zemp_ice(:, :)) * zcptn(:, :)
    IF (ln_mixcpl) THEN
      qns_tot(:, :) = qns(:, :) * ziceld(:, :) + SUM(qns_ice(:, :, :) * a_i(:, :, :), dim = 3)
      qns_tot(:, :) = qns_tot(:, :) * xcplmask(:, :, 0) + zqns_tot(:, :) * zmsk(:, :)
      DO jl = 1, jpl
        qns_ice(:, :, jl) = qns_ice(:, :, jl) * xcplmask(:, :, 0) + zqns_ice(:, :, jl) * zmsk(:, :)
      END DO
    ELSE
      qns_tot(:, :) = zqns_tot(:, :)
      qns_ice(:, :, :) = zqns_ice(:, :, :)
    END IF
    IF (srcv(jpr_cal) % laction) CALL iom_put('hflx_cal_cea', - frcv(jpr_cal) % z3(:, :, 1) * rLfus)
    IF (srcv(jpr_icb) % laction) CALL iom_put('hflx_icb_cea', - frcv(jpr_icb) % z3(:, :, 1) * rLfus)
    IF (iom_use('hflx_rain_cea')) CALL iom_put('hflx_rain_cea', (tprecip(:, :) - sprecip(:, :)) * zcptrain(:, :))
    IF (iom_use('hflx_evap_cea')) CALL iom_put('hflx_evap_cea', (frcv(jpr_tevp) % z3(:, :, 1) - frcv(jpr_ievp) % z3(:, :, 1) * picefr(:, :)) * zcptn(:, :) * tmask(:, :, 1))
    IF (iom_use('hflx_snow_cea')) CALL iom_put('hflx_snow_cea', sprecip(:, :) * (zcptsnw(:, :) - rLfus))
    IF (iom_use('hflx_snow_ao_cea')) CALL iom_put('hflx_snow_ao_cea', sprecip(:, :) * (zcptsnw(:, :) - rLfus) * (1._wp - zsnw(:, :)))
    IF (iom_use('hflx_snow_ai_cea')) CALL iom_put('hflx_snow_ai_cea', sprecip(:, :) * (zcptsnw(:, :) - rLfus) * zsnw(:, :))
    SELECT CASE (TRIM(sn_rcv_qsr % cldes))
    CASE ('oce only')
      zqsr_tot(:, :) = MAX(0._wp, frcv(jpr_qsroce) % z3(:, :, 1))
    CASE ('conservative')
      zqsr_tot(:, :) = frcv(jpr_qsrmix) % z3(:, :, 1)
      IF (TRIM(sn_rcv_qsr % clcat) == 'yes') THEN
        zqsr_ice(:, :, 1 : jpl) = frcv(jpr_qsrice) % z3(:, :, 1 : jpl)
      ELSE
        DO jl = 1, jpl
          zqsr_ice(:, :, jl) = frcv(jpr_qsrice) % z3(:, :, 1)
        END DO
      END IF
      zqsr_tot(:, :) = frcv(jpr_qsrmix) % z3(:, :, 1)
      zqsr_ice(:, :, 1) = frcv(jpr_qsrice) % z3(:, :, 1)
    CASE ('oce and ice')
      zqsr_tot(:, :) = ziceld(:, :) * frcv(jpr_qsroce) % z3(:, :, 1)
      IF (TRIM(sn_rcv_qsr % clcat) == 'yes') THEN
        DO jl = 1, jpl
          zqsr_tot(:, :) = zqsr_tot(:, :) + a_i(:, :, jl) * frcv(jpr_qsrice) % z3(:, :, jl)
          zqsr_ice(:, :, jl) = frcv(jpr_qsrice) % z3(:, :, jl)
        END DO
      ELSE
        qsr_tot(:, :) = qsr_tot(:, :) + picefr(:, :) * frcv(jpr_qsrice) % z3(:, :, 1)
        DO jl = 1, jpl
          zqsr_tot(:, :) = zqsr_tot(:, :) + picefr(:, :) * frcv(jpr_qsrice) % z3(:, :, 1)
          zqsr_ice(:, :, jl) = frcv(jpr_qsrice) % z3(:, :, 1)
        END DO
      END IF
    CASE ('mixed oce-ice')
      zqsr_tot(:, :) = frcv(jpr_qsrmix) % z3(:, :, 1)
      zqsr_ice(:, :, 1) = frcv(jpr_qsrmix) % z3(:, :, 1) * (1. - palbi(:, :, 1)) / (1. - (alb_oce_mix(:, :) * ziceld(:, :) + palbi(:, :, 1) * picefr(:, :)))
    CASE ('none')
      CALL ctl_stop('STOP', 'sbccpl/sbc_cpl_ice_flx: some fields are not defined. Change sn_rcv_qsr value in namelist namsbc_cpl')
    END SELECT
    IF (ln_dm2dc .AND. ln_cpl) THEN
      zqsr_tot(:, :) = sbc_dcy(zqsr_tot(:, :))
      DO jl = 1, jpl
        zqsr_ice(:, :, jl) = sbc_dcy(zqsr_ice(:, :, jl))
      END DO
    END IF
    IF (ln_mixcpl) THEN
      qsr_tot(:, :) = qsr(:, :) * ziceld(:, :) + SUM(qsr_ice(:, :, :) * a_i(:, :, :), dim = 3)
      qsr_tot(:, :) = qsr_tot(:, :) * xcplmask(:, :, 0) + zqsr_tot(:, :) * zmsk(:, :)
      DO jl = 1, jpl
        qsr_ice(:, :, jl) = qsr_ice(:, :, jl) * xcplmask(:, :, 0) + zqsr_ice(:, :, jl) * zmsk(:, :)
      END DO
    ELSE
      qsr_tot(:, :) = zqsr_tot(:, :)
      qsr_ice(:, :, :) = zqsr_ice(:, :, :)
    END IF
    SELECT CASE (TRIM(sn_rcv_dqnsdt % cldes))
    CASE ('coupled')
      IF (TRIM(sn_rcv_dqnsdt % clcat) == 'yes') THEN
        zdqns_ice(:, :, 1 : jpl) = frcv(jpr_dqnsdt) % z3(:, :, 1 : jpl)
      ELSE
        DO jl = 1, jpl
          zdqns_ice(:, :, jl) = frcv(jpr_dqnsdt) % z3(:, :, 1)
        END DO
      END IF
    END SELECT
    IF (ln_mixcpl) THEN
      DO jl = 1, jpl
        dqns_ice(:, :, jl) = dqns_ice(:, :, jl) * xcplmask(:, :, 0) + zdqns_ice(:, :, jl) * zmsk(:, :)
      END DO
    ELSE
      dqns_ice(:, :, :) = zdqns_ice(:, :, :)
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sbc_cpl_ice_flx
  SUBROUTINE sbc_cpl_snd(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jl
    INTEGER :: isec, info
    REAL(KIND = wp) :: zumax, zvmax
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfr_l, ztmp1, ztmp2, zotx1, zoty1, zotz1, zitx1, zity1, zitz1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: ztmp3, ztmp4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('sbc_cpl_snd', 'r0', 0, 0)
    isec = (kt - nit000) * NINT(rdt)
    zfr_l(:, :) = 1. - fr_i(:, :)
    IF (ssnd(jps_toce) % laction .OR. ssnd(jps_tice) % laction .OR. ssnd(jps_tmix) % laction) THEN
      IF (nn_components == jp_iam_opa) THEN
        ztmp1(:, :) = tsn(:, :, 1, jp_tem)
      ELSE
        IF (l_usect) THEN
          ztmp1(:, :) = eos_pt_from_ct(tsn(:, :, 1, jp_tem), tsn(:, :, 1, jp_sal))
        ELSE
          ztmp1(:, :) = tsn(:, :, 1, jp_tem)
        END IF
        SELECT CASE (sn_snd_temp % cldes)
        CASE ('oce only')
          ztmp1(:, :) = ztmp1(:, :) + rt0
        CASE ('oce and ice')
          ztmp1(:, :) = ztmp1(:, :) + rt0
          SELECT CASE (sn_snd_temp % clcat)
          CASE ('yes')
            ztmp3(:, :, 1 : jpl) = tn_ice(:, :, 1 : jpl)
          CASE ('no')
            WHERE (SUM(a_i, dim = 3) /= 0.)
              ztmp3(:, :, 1) = SUM(tn_ice * a_i, dim = 3) / SUM(a_i, dim = 3)
            ELSEWHERE
              ztmp3(:, :, 1) = rt0
            END WHERE
          CASE DEFAULT
            CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_temp%clcat')
          END SELECT
        CASE ('weighted oce and ice')
          ztmp1(:, :) = (ztmp1(:, :) + rt0) * zfr_l(:, :)
          SELECT CASE (sn_snd_temp % clcat)
          CASE ('yes')
            ztmp3(:, :, 1 : jpl) = tn_ice(:, :, 1 : jpl) * a_i(:, :, 1 : jpl)
          CASE ('no')
            ztmp3(:, :, :) = 0.0
            DO jl = 1, jpl
              ztmp3(:, :, 1) = ztmp3(:, :, 1) + tn_ice(:, :, jl) * a_i(:, :, jl)
            END DO
          CASE DEFAULT
            CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_temp%clcat')
          END SELECT
        CASE ('oce and weighted ice')
          ztmp1(:, :) = tsn(:, :, 1, jp_tem) + rt0
          SELECT CASE (sn_snd_temp % clcat)
          CASE ('yes')
            ztmp3(:, :, 1 : jpl) = tn_ice(:, :, 1 : jpl) * a_i(:, :, 1 : jpl)
          CASE ('no')
            ztmp3(:, :, :) = 0.0
            DO jl = 1, jpl
              ztmp3(:, :, 1) = ztmp3(:, :, 1) + tn_ice(:, :, jl) * a_i(:, :, jl)
            END DO
          CASE DEFAULT
            CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_temp%clcat')
          END SELECT
        CASE ('mixed oce-ice')
          ztmp1(:, :) = (ztmp1(:, :) + rt0) * zfr_l(:, :)
          DO jl = 1, jpl
            ztmp1(:, :) = ztmp1(:, :) + tn_ice(:, :, jl) * a_i(:, :, jl)
          END DO
        CASE DEFAULT
          CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_temp%cldes')
        END SELECT
      END IF
      IF (ssnd(jps_toce) % laction) CALL cpl_snd(jps_toce, isec, RESHAPE(ztmp1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_tice) % laction) CALL cpl_snd(jps_tice, isec, ztmp3, info)
      IF (ssnd(jps_tmix) % laction) CALL cpl_snd(jps_tmix, isec, RESHAPE(ztmp1, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_albice) % laction) THEN
      SELECT CASE (sn_snd_alb % cldes)
      CASE ('ice')
        SELECT CASE (sn_snd_alb % clcat)
        CASE ('yes')
          ztmp3(:, :, 1 : jpl) = alb_ice(:, :, 1 : jpl)
        CASE ('no')
          WHERE (SUM(a_i, dim = 3) /= 0.)
            ztmp1(:, :) = SUM(alb_ice(:, :, 1 : jpl) * a_i(:, :, 1 : jpl), dim = 3) / SUM(a_i(:, :, 1 : jpl), dim = 3)
          ELSEWHERE
            ztmp1(:, :) = alb_oce_mix(:, :)
          END WHERE
        CASE DEFAULT
          CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_alb%clcat')
        END SELECT
      CASE ('weighted ice')
        SELECT CASE (sn_snd_alb % clcat)
        CASE ('yes')
          ztmp3(:, :, 1 : jpl) = alb_ice(:, :, 1 : jpl) * a_i(:, :, 1 : jpl)
        CASE ('no')
          WHERE (fr_i(:, :) > 0.)
            ztmp1(:, :) = SUM(alb_ice(:, :, 1 : jpl) * a_i(:, :, 1 : jpl), dim = 3)
          ELSEWHERE
            ztmp1(:, :) = 0.
          END WHERE
        CASE DEFAULT
          CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_ice%clcat')
        END SELECT
      CASE DEFAULT
        CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_alb%cldes')
      END SELECT
      SELECT CASE (sn_snd_alb % clcat)
      CASE ('yes')
        CALL cpl_snd(jps_albice, isec, ztmp3, info)
      CASE ('no')
        CALL cpl_snd(jps_albice, isec, RESHAPE(ztmp1, (/jpi, jpj, 1/)), info)
      END SELECT
    END IF
    IF (ssnd(jps_albmix) % laction) THEN
      ztmp1(:, :) = alb_oce_mix(:, :) * zfr_l(:, :)
      DO jl = 1, jpl
        ztmp1(:, :) = ztmp1(:, :) + alb_ice(:, :, jl) * a_i(:, :, jl)
      END DO
      CALL cpl_snd(jps_albmix, isec, RESHAPE(ztmp1, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_fice) % laction) THEN
      SELECT CASE (sn_snd_thick % clcat)
      CASE ('yes')
        ztmp3(:, :, 1 : jpl) = a_i(:, :, 1 : jpl)
      CASE ('no')
        ztmp3(:, :, 1) = fr_i(:, :)
      CASE DEFAULT
        CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_thick%clcat')
      END SELECT
      IF (ssnd(jps_fice) % laction) CALL cpl_snd(jps_fice, isec, ztmp3, info)
    END IF
    IF (ssnd(jps_fice1) % laction) THEN
      SELECT CASE (sn_snd_thick1 % clcat)
      CASE ('yes')
        ztmp3(:, :, 1 : jpl) = a_i(:, :, 1 : jpl)
      CASE ('no')
        ztmp3(:, :, 1) = fr_i(:, :)
      CASE DEFAULT
        CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_thick1%clcat')
      END SELECT
      CALL cpl_snd(jps_fice1, isec, ztmp3, info)
    END IF
    IF (ssnd(jps_fice2) % laction) THEN
      ztmp3(:, :, 1) = fr_i(:, :)
      IF (ssnd(jps_fice2) % laction) CALL cpl_snd(jps_fice2, isec, ztmp3, info)
    END IF
    IF (ssnd(jps_hice) % laction .OR. ssnd(jps_hsnw) % laction) THEN
      SELECT CASE (sn_snd_thick % cldes)
      CASE ('none')
      CASE ('weighted ice and snow')
        SELECT CASE (sn_snd_thick % clcat)
        CASE ('yes')
          ztmp3(:, :, 1 : jpl) = h_i(:, :, 1 : jpl) * a_i(:, :, 1 : jpl)
          ztmp4(:, :, 1 : jpl) = h_s(:, :, 1 : jpl) * a_i(:, :, 1 : jpl)
        CASE ('no')
          ztmp3(:, :, :) = 0.0
          ztmp4(:, :, :) = 0.0
          DO jl = 1, jpl
            ztmp3(:, :, 1) = ztmp3(:, :, 1) + h_i(:, :, jl) * a_i(:, :, jl)
            ztmp4(:, :, 1) = ztmp4(:, :, 1) + h_s(:, :, jl) * a_i(:, :, jl)
          END DO
        CASE DEFAULT
          CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_thick%clcat')
        END SELECT
      CASE ('ice and snow')
        SELECT CASE (sn_snd_thick % clcat)
        CASE ('yes')
          ztmp3(:, :, 1 : jpl) = h_i(:, :, 1 : jpl)
          ztmp4(:, :, 1 : jpl) = h_s(:, :, 1 : jpl)
        CASE ('no')
          WHERE (SUM(a_i, dim = 3) /= 0.)
            ztmp3(:, :, 1) = SUM(h_i * a_i, dim = 3) / SUM(a_i, dim = 3)
            ztmp4(:, :, 1) = SUM(h_s * a_i, dim = 3) / SUM(a_i, dim = 3)
          ELSEWHERE
            ztmp3(:, :, 1) = 0.
            ztmp4(:, :, 1) = 0.
          END WHERE
        CASE DEFAULT
          CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_thick%clcat')
        END SELECT
      CASE DEFAULT
        CALL ctl_stop('sbc_cpl_snd: wrong definition of sn_snd_thick%cldes')
      END SELECT
      IF (ssnd(jps_hice) % laction) CALL cpl_snd(jps_hice, isec, ztmp3, info)
      IF (ssnd(jps_hsnw) % laction) CALL cpl_snd(jps_hsnw, isec, ztmp4, info)
    END IF
    IF (ssnd(jps_co2) % laction .AND. l_co2cpl) CALL cpl_snd(jps_co2, isec, RESHAPE(oce_co2, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_ocx1) % laction) THEN
      IF (nn_components == jp_iam_opa) THEN
        zotx1(:, :) = un(:, :, 1)
        zoty1(:, :) = vn(:, :, 1)
      ELSE
        SELECT CASE (TRIM(sn_snd_crt % cldes))
        CASE ('oce only')
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zotx1(ji, jj) = 0.5 * (un(ji, jj, 1) + un(ji - 1, jj, 1))
              zoty1(ji, jj) = 0.5 * (vn(ji, jj, 1) + vn(ji, jj - 1, 1))
            END DO
          END DO
        CASE ('weighted oce and ice')
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zotx1(ji, jj) = 0.5 * (un(ji, jj, 1) + un(ji - 1, jj, 1)) * zfr_l(ji, jj)
              zoty1(ji, jj) = 0.5 * (vn(ji, jj, 1) + vn(ji, jj - 1, 1)) * zfr_l(ji, jj)
              zitx1(ji, jj) = 0.5 * (u_ice(ji, jj) + u_ice(ji - 1, jj)) * fr_i(ji, jj)
              zity1(ji, jj) = 0.5 * (v_ice(ji, jj) + v_ice(ji, jj - 1)) * fr_i(ji, jj)
            END DO
          END DO
          CALL lbc_lnk_multi('sbccpl', zitx1, 'T', - 1., zity1, 'T', - 1.)
        CASE ('mixed oce-ice')
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zotx1(ji, jj) = 0.5 * (un(ji, jj, 1) + un(ji - 1, jj, 1)) * zfr_l(ji, jj) + 0.5 * (u_ice(ji, jj) + u_ice(ji - 1, jj)) * fr_i(ji, jj)
              zoty1(ji, jj) = 0.5 * (vn(ji, jj, 1) + vn(ji, jj - 1, 1)) * zfr_l(ji, jj) + 0.5 * (v_ice(ji, jj) + v_ice(ji, jj - 1)) * fr_i(ji, jj)
            END DO
          END DO
        END SELECT
        CALL lbc_lnk_multi('sbccpl', zotx1, ssnd(jps_ocx1) % clgrid, - 1., zoty1, ssnd(jps_ocy1) % clgrid, - 1.)
      END IF
      IF (TRIM(sn_snd_crt % clvor) == 'eastward-northward') THEN
        CALL rot_rep(zotx1, zoty1, ssnd(jps_ocx1) % clgrid, 'ij->e', ztmp1)
        CALL rot_rep(zotx1, zoty1, ssnd(jps_ocx1) % clgrid, 'ij->n', ztmp2)
        zotx1(:, :) = ztmp1(:, :)
        zoty1(:, :) = ztmp2(:, :)
        IF (ssnd(jps_ivx1) % laction) THEN
          CALL rot_rep(zitx1, zity1, ssnd(jps_ivx1) % clgrid, 'ij->e', ztmp1)
          CALL rot_rep(zitx1, zity1, ssnd(jps_ivx1) % clgrid, 'ij->n', ztmp2)
          zitx1(:, :) = ztmp1(:, :)
          zity1(:, :) = ztmp2(:, :)
        END IF
      END IF
      IF (TRIM(sn_snd_crt % clvref) == 'cartesian') THEN
        ztmp1(:, :) = zotx1(:, :)
        ztmp2(:, :) = zoty1(:, :)
        CALL oce2geo(ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1)
        IF (ssnd(jps_ivx1) % laction) THEN
          ztmp1(:, :) = zitx1(:, :)
          ztmp1(:, :) = zity1(:, :)
          CALL oce2geo(ztmp1, ztmp2, 'T', zitx1, zity1, zitz1)
        END IF
      END IF
      IF (ssnd(jps_ocx1) % laction) CALL cpl_snd(jps_ocx1, isec, RESHAPE(zotx1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_ocy1) % laction) CALL cpl_snd(jps_ocy1, isec, RESHAPE(zoty1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_ocz1) % laction) CALL cpl_snd(jps_ocz1, isec, RESHAPE(zotz1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_ivx1) % laction) CALL cpl_snd(jps_ivx1, isec, RESHAPE(zitx1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_ivy1) % laction) CALL cpl_snd(jps_ivy1, isec, RESHAPE(zity1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_ivz1) % laction) CALL cpl_snd(jps_ivz1, isec, RESHAPE(zitz1, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_ocxw) % laction .OR. ssnd(jps_ocyw) % laction) THEN
      SELECT CASE (TRIM(sn_snd_crtw % cldes))
      CASE ('oce only')
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zotx1(ji, jj) = 0.5 * (un(ji, jj, 1) + un(ji - 1, jj, 1))
            zoty1(ji, jj) = 0.5 * (vn(ji, jj, 1) + vn(ji, jj - 1, 1))
          END DO
        END DO
      CASE ('weighted oce and ice')
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zotx1(ji, jj) = 0.5 * (un(ji, jj, 1) + un(ji - 1, jj, 1)) * zfr_l(ji, jj)
            zoty1(ji, jj) = 0.5 * (vn(ji, jj, 1) + vn(ji, jj - 1, 1)) * zfr_l(ji, jj)
            zitx1(ji, jj) = 0.5 * (u_ice(ji, jj) + u_ice(ji - 1, jj)) * fr_i(ji, jj)
            zity1(ji, jj) = 0.5 * (v_ice(ji, jj) + v_ice(ji, jj - 1)) * fr_i(ji, jj)
          END DO
        END DO
        CALL lbc_lnk_multi('sbccpl', zitx1, 'T', - 1., zity1, 'T', - 1.)
      CASE ('mixed oce-ice')
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zotx1(ji, jj) = 0.5 * (un(ji, jj, 1) + un(ji - 1, jj, 1)) * zfr_l(ji, jj) + 0.5 * (u_ice(ji, jj) + u_ice(ji - 1, jj)) * fr_i(ji, jj)
            zoty1(ji, jj) = 0.5 * (vn(ji, jj, 1) + vn(ji, jj - 1, 1)) * zfr_l(ji, jj) + 0.5 * (v_ice(ji, jj) + v_ice(ji, jj - 1)) * fr_i(ji, jj)
          END DO
        END DO
      END SELECT
      CALL lbc_lnk_multi('sbccpl', zotx1, ssnd(jps_ocxw) % clgrid, - 1., zoty1, ssnd(jps_ocyw) % clgrid, - 1.)
      IF (TRIM(sn_snd_crtw % clvor) == 'eastward-northward') THEN
        CALL rot_rep(zotx1, zoty1, ssnd(jps_ocxw) % clgrid, 'ij->e', ztmp1)
        CALL rot_rep(zotx1, zoty1, ssnd(jps_ocxw) % clgrid, 'ij->n', ztmp2)
        zotx1(:, :) = ztmp1(:, :)
        zoty1(:, :) = ztmp2(:, :)
        IF (ssnd(jps_ivx1) % laction) THEN
          CALL rot_rep(zitx1, zity1, ssnd(jps_ivx1) % clgrid, 'ij->e', ztmp1)
          CALL rot_rep(zitx1, zity1, ssnd(jps_ivx1) % clgrid, 'ij->n', ztmp2)
          zitx1(:, :) = ztmp1(:, :)
          zity1(:, :) = ztmp2(:, :)
        END IF
      END IF
      IF (ssnd(jps_ocxw) % laction) CALL cpl_snd(jps_ocxw, isec, RESHAPE(zotx1, (/jpi, jpj, 1/)), info)
      IF (ssnd(jps_ocyw) % laction) CALL cpl_snd(jps_ocyw, isec, RESHAPE(zoty1, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_ficet) % laction) THEN
      CALL cpl_snd(jps_ficet, isec, RESHAPE(fr_i, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_wlev) % laction) THEN
      IF (ln_apr_dyn) THEN
        IF (kt /= nit000) THEN
          ztmp1(:, :) = sshb(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :))
        ELSE
          ztmp1(:, :) = sshb(:, :)
        END IF
      ELSE
        ztmp1(:, :) = sshn(:, :)
      END IF
      CALL cpl_snd(jps_wlev, isec, RESHAPE(ztmp1, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_ssh) % laction) THEN
      IF (ln_apr_dyn) THEN
        ztmp1(:, :) = sshb(:, :) - 0.5 * (ssh_ib(:, :) + ssh_ibb(:, :))
      ELSE
        ztmp1(:, :) = sshn(:, :)
      END IF
      CALL cpl_snd(jps_ssh, isec, RESHAPE(ztmp1, (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_soce) % laction) THEN
      CALL cpl_snd(jps_soce, isec, RESHAPE(tsn(:, :, 1, jp_sal), (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_e3t1st) % laction) THEN
      CALL cpl_snd(jps_e3t1st, isec, RESHAPE(e3t_n(:, :, 1), (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_fraqsr) % laction) THEN
      CALL cpl_snd(jps_fraqsr, isec, RESHAPE(fraqsr_1lev(:, :), (/jpi, jpj, 1/)), info)
    END IF
    IF (ssnd(jps_qsroce) % laction) CALL cpl_snd(jps_qsroce, isec, RESHAPE(qsr, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_qnsoce) % laction) CALL cpl_snd(jps_qnsoce, isec, RESHAPE(qns, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_oemp) % laction) CALL cpl_snd(jps_oemp, isec, RESHAPE(emp, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_sflx) % laction) CALL cpl_snd(jps_sflx, isec, RESHAPE(sfx, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_otx1) % laction) CALL cpl_snd(jps_otx1, isec, RESHAPE(utau, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_oty1) % laction) CALL cpl_snd(jps_oty1, isec, RESHAPE(vtau, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_rnf) % laction) CALL cpl_snd(jps_rnf, isec, RESHAPE(rnf, (/jpi, jpj, 1/)), info)
    IF (ssnd(jps_taum) % laction) CALL cpl_snd(jps_taum, isec, RESHAPE(taum, (/jpi, jpj, 1/)), info)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE sbc_cpl_snd
END MODULE sbccpl