MODULE zdfphy
  USE oce
  USE zdf_oce
  USE zdfdrg
  USE zdfsh2
  USE zdfric
  USE zdftke
  USE zdfgls
  USE zdfosm
  USE zdfddm
  USE zdfevd
  USE zdfiwm
  USE zdftmx
  USE zdfswm
  USE zdfmxl
  USE tranpc
  USE trc_oce
  USE sbc_oce
  USE sbcrnf
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_phy_init
  PUBLIC :: zdf_phy
  INTEGER :: nzdf_phy
  INTEGER, PARAMETER :: np_CST = 1
  INTEGER, PARAMETER :: np_RIC = 2
  INTEGER, PARAMETER :: np_TKE = 3
  INTEGER, PARAMETER :: np_GLS = 4
  INTEGER, PARAMETER :: np_OSM = 5
  LOGICAL :: l_zdfsh2
  CONTAINS
  SUBROUTINE zdf_phy_init
    INTEGER :: jk
    INTEGER :: ioptio, ios
    NAMELIST /namzdf/ ln_zdfcst, ln_zdfric, ln_zdftke, ln_zdfgls, ln_zdfosm, ln_zdfevd, nn_evdm, rn_evd, ln_zdfnpc, nn_npc, nn_npcp, ln_zdfddm, rn_avts, rn_hsbfr, ln_zdfswm, ln_zdfiwm, ln_zdftmx, ln_zad_Aimp, rn_avm0, rn_avt0, nn_avb, nn_havtb
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'zdf_phy_init: ocean vertical physics'
      WRITE(numout, *) '~~~~~~~~~~~~'
    END IF
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namzdf in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf)
    IF (lwp) THEN
      WRITE(numout, *) '   Namelist namzdf : set vertical mixing mixing parameters'
      WRITE(numout, *) '      adaptive-implicit vertical advection'
      WRITE(numout, *) '         Courant number targeted application   ln_zad_Aimp = ', ln_zad_Aimp
      WRITE(numout, *) '      vertical closure scheme'
      WRITE(numout, *) '         constant vertical mixing coefficient    ln_zdfcst = ', ln_zdfcst
      WRITE(numout, *) '         Richardson number dependent closure     ln_zdfric = ', ln_zdfric
      WRITE(numout, *) '         Turbulent Kinetic Energy closure (TKE)  ln_zdftke = ', ln_zdftke
      WRITE(numout, *) '         Generic Length Scale closure (GLS)      ln_zdfgls = ', ln_zdfgls
      WRITE(numout, *) '         OSMOSIS-OBL closure (OSM)               ln_zdfosm = ', ln_zdfosm
      WRITE(numout, *) '      convection: '
      WRITE(numout, *) '         enhanced vertical diffusion             ln_zdfevd = ', ln_zdfevd
      WRITE(numout, *) '            applied on momentum (=1/0)             nn_evdm = ', nn_evdm
      WRITE(numout, *) '            vertical coefficient for evd           rn_evd  = ', rn_evd
      WRITE(numout, *) '         non-penetrative convection (npc)        ln_zdfnpc = ', ln_zdfnpc
      WRITE(numout, *) '            npc call  frequency                    nn_npc  = ', nn_npc
      WRITE(numout, *) '            npc print frequency                    nn_npcp = ', nn_npcp
      WRITE(numout, *) '      double diffusive mixing                    ln_zdfddm = ', ln_zdfddm
      WRITE(numout, *) '         maximum avs for dd mixing                 rn_avts = ', rn_avts
      WRITE(numout, *) '         heat/salt buoyancy flux ratio             rn_hsbfr= ', rn_hsbfr
      WRITE(numout, *) '      gravity wave-induced mixing'
      WRITE(numout, *) '         surface  wave (Qiao et al 2010)         ln_zdfswm = ', ln_zdfswm
      WRITE(numout, *) '         internal wave (de Lavergne et al 2017)  ln_zdfiwm = ', ln_zdfiwm
      WRITE(numout, *) '      coefficients : '
      WRITE(numout, *) '         vertical eddy viscosity                 rn_avm0   = ', rn_avm0
      WRITE(numout, *) '         vertical eddy diffusivity               rn_avt0   = ', rn_avt0
      WRITE(numout, *) '         constant background or profile          nn_avb    = ', nn_avb
      WRITE(numout, *) '         horizontal variation for avtb           nn_havtb  = ', nn_havtb
    END IF
    IF (ln_zad_Aimp) THEN
      IF (zdf_phy_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_phy_init : unable to allocate adaptive-implicit z-advection arrays')
      !$ACC KERNELS
      Cu_adv(:, :, :) = 0._wp
      wi(:, :, :) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (nn_avb == 0) THEN
      !$ACC KERNELS
      avmb(:) = rn_avm0
      avtb(:) = rn_avt0
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      avmb(:) = rn_avm0
      avtb(:) = rn_avt0 + (3.E-4_wp - 2._wp * rn_avt0) * 1.E-4_wp * gdepw_1d(:)
      !$ACC END KERNELS
      IF (ln_sco .AND. lwp) CALL ctl_warn('avtb profile not valid in sco')
    END IF
    !$ACC KERNELS
    avtb_2d(:, :) = 1._wp
    !$ACC END KERNELS
    IF (nn_havtb == 1) THEN
      WHERE (- 15. <= gphit .AND. gphit < - 5) avtb_2d = (1. - 0.09 * (gphit + 15.))
      WHERE (- 5. <= gphit .AND. gphit < 5) avtb_2d = 0.1
      WHERE (5. <= gphit .AND. gphit < 15) avtb_2d = (0.1 + 0.09 * (gphit - 5.))
    END IF
    !$ACC KERNELS
    DO jk = 1, jpk
      avt_k(:, :, jk) = avtb_2d(:, :) * avtb(jk) * wmask(:, :, jk)
      avm_k(:, :, jk) = avmb(jk) * wmask(:, :, jk)
    END DO
    avt(:, :, :) = 0._wp
    avs(:, :, :) = 0._wp
    avm(:, :, :) = 0._wp
    !$ACC END KERNELS
    IF (ln_zdfnpc .AND. ln_zdfevd) CALL ctl_stop('zdf_phy_init: chose between ln_zdfnpc and ln_zdfevd')
    IF (ln_zdfosm .AND. ln_zdfevd) CALL ctl_stop('zdf_phy_init: chose between ln_zdfosm and ln_zdfevd')
    IF (lk_top .AND. ln_zdfnpc) CALL ctl_stop('zdf_phy_init: npc scheme is not working with key_top')
    IF (lk_top .AND. ln_zdfosm) CALL ctl_stop('zdf_phy_init: osmosis scheme is not working with key_top')
    IF (lwp) THEN
      WRITE(numout, *)
      IF (ln_zdfnpc) THEN
        WRITE(numout, *) '   ==>>>   convection: use non penetrative convective scheme'
      ELSE IF (ln_zdfevd) THEN
        WRITE(numout, *) '   ==>>>   convection: use enhanced vertical diffusion scheme'
      ELSE
        WRITE(numout, *) '   ==>>>   convection: no specific scheme used'
      END IF
    END IF
    IF (lwp) THEN
      WRITE(numout, *)
      IF (ln_zdfddm) THEN
        WRITE(numout, *) '   ==>>>   use double diffusive mixing: avs /= avt'
      ELSE
        WRITE(numout, *) '   ==>>>   No  double diffusive mixing: avs = avt'
      END IF
    END IF
    ioptio = 0
    IF (ln_zdfcst) THEN
      ioptio = ioptio + 1
      nzdf_phy = np_cst
    END IF
    IF (ln_zdfric) THEN
      ioptio = ioptio + 1
      nzdf_phy = np_ric
      CALL zdf_ric_init
    END IF
    IF (ln_zdftke) THEN
      ioptio = ioptio + 1
      nzdf_phy = np_tke
      CALL zdf_tke_init
    END IF
    IF (ln_zdfgls) THEN
      ioptio = ioptio + 1
      nzdf_phy = np_gls
      CALL zdf_gls_init
    END IF
    IF (ln_zdfosm) THEN
      ioptio = ioptio + 1
      nzdf_phy = np_osm
      CALL zdf_osm_init
    END IF
    IF (ioptio /= 1) CALL ctl_stop('zdf_phy_init: one and only one vertical diffusion option has to be defined ')
    IF (ln_isfcav) THEN
      IF (ln_zdfric .OR. ln_zdfgls) CALL ctl_stop('zdf_phy_init: zdfric and zdfgls never tested with ice shelves cavities ')
    END IF
    IF (ln_zdfcst) THEN
      l_zdfsh2 = .FALSE.
    ELSE
      l_zdfsh2 = .TRUE.
    END IF
    IF (ln_zdfiwm) CALL zdf_iwm_init
    IF (ln_zdftmx) CALL zdf_tmx_init
    IF (ln_zdfswm) CALL zdf_swm_init
    CALL zdf_drg_init
  END SUBROUTINE zdf_phy_init
  SUBROUTINE zdf_phy(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zsh2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('zdf_phy', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('zdf_phy')
    IF (l_zdfdrg) THEN
      CALL zdf_drg(kt, mbkt, r_Cdmin_bot, r_Cdmax_bot, r_z0_bot, r_ke0_bot, rCd0_bot, rCdU_bot)
      IF (ln_isfcav) THEN
        CALL zdf_drg(kt, mikt, r_Cdmin_top, r_Cdmax_top, r_z0_top, r_ke0_top, rCd0_top, rCdU_top)
      END IF
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zsh2(:, :, :) = 0.
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('zdf_phy', 'r1', 0, 0)
    IF (l_zdfsh2) CALL zdf_sh2(ub, vb, un, vn, avm_k, zsh2)
    SELECT CASE (nzdf_phy)
    CASE (np_ric)
      CALL zdf_ric(kt, gdept_n, zsh2, avm_k, avt_k)
    CASE (np_tke)
      CALL zdf_tke(kt, zsh2, avm_k, avt_k)
    CASE (np_gls)
      CALL zdf_gls(kt, zsh2, avm_k, avt_k)
    CASE (np_osm)
      CALL zdf_osm(kt, avm_k, avt_k)
    END SELECT
    CALL profile_psy_data1 % PostEnd
    !$ACC KERNELS
    avt(:, :, 2 : jpkm1) = avt_k(:, :, 2 : jpkm1)
    avm(:, :, 2 : jpkm1) = avm_k(:, :, 2 : jpkm1)
    !$ACC END KERNELS
    IF (ln_rnf_mouth) THEN
      !$ACC KERNELS
      DO jk = 2, nkrnf
        avt(:, :, jk) = avt(:, :, jk) + 2._wp * rn_avt_rnf * rnfmsk(:, :) * wmask(:, :, jk)
      END DO
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data2 % PreStart('zdf_phy', 'r2', 0, 0)
    IF (ln_zdfevd) CALL zdf_evd(kt, avm, avt)
    CALL profile_psy_data2 % PostEnd
    IF (ln_zdfddm) THEN
      CALL profile_psy_data3 % PreStart('zdf_phy', 'r3', 0, 0)
      CALL zdf_ddm(kt, avm, avt, avs)
      CALL profile_psy_data3 % PostEnd
    ELSE
      !$ACC KERNELS
      avs(2 : jpim1, 2 : jpjm1, 1 : jpkm1) = avt(2 : jpim1, 2 : jpjm1, 1 : jpkm1)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data4 % PreStart('zdf_phy', 'r4', 0, 0)
    IF (ln_zdfswm) CALL zdf_swm(kt, avm, avt, avs)
    IF (ln_zdfiwm) CALL zdf_iwm(kt, avm, avt, avs)
    IF (ln_zdftmx) CALL zdf_tmx(kt, avm, avt, avs)
    IF (l_zdfsh2) THEN
      CALL lbc_lnk_multi('zdfphy', avm_k, 'W', 1., avt_k, 'W', 1., avm, 'W', 1., avt, 'W', 1., avs, 'W', 1.)
    ELSE
      CALL lbc_lnk_multi('zdfphy', avm, 'W', 1., avt, 'W', 1., avs, 'W', 1.)
    END IF
    IF (l_zdfdrg) THEN
      IF (ln_isfcav) THEN
        CALL lbc_lnk_multi('zdfphy', rcdu_top, 'T', 1., rcdu_bot, 'T', 1.)
      ELSE
        CALL lbc_lnk('zdfphy', rcdu_bot, 'T', 1.)
      END IF
    END IF
    CALL zdf_mxl(kt)
    IF (lrst_oce) THEN
      IF (ln_zdftke) CALL tke_rst(kt, 'WRITE')
      IF (ln_zdfgls) CALL gls_rst(kt, 'WRITE')
      IF (ln_zdfric) CALL ric_rst(kt, 'WRITE')
    END IF
    IF (ln_timing) CALL timing_stop('zdf_phy')
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE zdf_phy
  INTEGER FUNCTION zdf_phy_alloc()
    ALLOCATE(wi(jpi, jpj, jpk), Cu_adv(jpi, jpj, jpk), STAT = zdf_phy_alloc)
    IF (zdf_phy_alloc /= 0) CALL ctl_warn('zdf_phy_alloc: failed to allocate ln_zad_Aimp=T required arrays')
    CALL mpp_sum('zdfphy', zdf_phy_alloc)
  END FUNCTION zdf_phy_alloc
END MODULE zdfphy