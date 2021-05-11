MODULE zdfosm
  USE oce
  USE dom_oce
  USE zdf_oce
  USE sbc_oce
  USE sbcwave
  USE phycst
  USE eosbn2
  USE traqsr
  USE zdfddm
  USE iom
  USE lib_mpp
  USE trd_oce
  USE trdtra
  USE in_out_manager
  USE lbclnk
  USE prtctl
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_osm
  PUBLIC :: zdf_osm_init
  PUBLIC :: osm_rst
  PUBLIC :: tra_osm
  PUBLIC :: trc_osm
  PUBLIC :: dyn_osm
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ghamu
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ghamv
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ghamt
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ghams
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: etmean
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hbl
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hbli
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: dstokes
  LOGICAL :: ln_use_osm_la
  REAL(KIND = wp) :: rn_osm_la
  REAL(KIND = wp) :: rn_osm_dstokes
  REAL(KIND = wp) :: rn_osm_hbl0 = 10._wp
  INTEGER :: nn_ave
  INTEGER :: nn_osm_wave = 0
  LOGICAL :: ln_dia_osm
  LOGICAL :: ln_kpprimix = .TRUE.
  REAL(KIND = wp) :: rn_riinfty = 0.7
  REAL(KIND = wp) :: rn_difri = 0.005
  LOGICAL :: ln_convmix = .TRUE.
  REAL(KIND = wp) :: rn_difconv = 1._wp
  REAL(KIND = wp) :: epsln = 1.0E-20_wp
  REAL(KIND = wp) :: pthird = 1._wp / 3._wp
  REAL(KIND = wp) :: p2third = 2._wp / 3._wp
  INTEGER :: idebug = 236
  INTEGER :: jdebug = 228
  CONTAINS
  INTEGER FUNCTION zdf_osm_alloc()
    ALLOCATE(ghamu(jpi, jpj, jpk), ghamv(jpi, jpj, jpk), ghamt(jpi, jpj, jpk), ghams(jpi, jpj, jpk), hbl(jpi, jpj), hbli(jpi, jpj), dstokes(jpi, jpj), etmean(jpi, jpj, jpk), STAT = zdf_osm_alloc)
    IF (zdf_osm_alloc /= 0) CALL ctl_warn('zdf_osm_alloc: failed to allocate zdf_osm arrays')
    CALL mpp_sum('zdfosm', zdf_osm_alloc)
  END FUNCTION zdf_osm_alloc
  SUBROUTINE zdf_osm(kt, p_avm, p_avt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm, p_avt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbot, jkmax, jkm1, jkp2
    REAL(KIND = wp) :: ztx, zty, zflageos, zstabl, zbuofdep, zucube
    REAL(KIND = wp) :: zbeta, zthermal
    REAL(KIND = wp) :: zehat, zeta, zhrib, zsig, zscale, zwst, zws, zwm
    REAL(KIND = wp) :: zwsun, zwmun, zcons, zconm, zwcons, zwconm
    REAL(KIND = wp) :: zsr, zbw, ze, zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zcomp, zrhd, zrhdr, zbvzed
    INTEGER :: jm
    REAL(KIND = wp) :: zr1, zr2, zr3, zr4, zrhop
    REAL(KIND = wp) :: zflag, zrn2, zdep21, zdep32, zdep43
    REAL(KIND = wp) :: zesh2, zri, zfri
    REAL(KIND = wp) :: zdelta, zdelta2, zdzup, zdzdn, zdzh, zvath, zgat1, zdat1, zkm1m, zkm1t
    REAL(KIND = wp) :: zt, zs, zu, zv, zrh
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zrad0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zradh
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zradav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zustar
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwstrl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zvstr
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwstrc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zuw0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zvw0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwth0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zws0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwb0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwthav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwsav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwbav
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwb_ent
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zustke
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zla
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zcos_wind
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsin_wind
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhol
    LOGICAL, DIMENSION(:, :), ALLOCATABLE :: lconv
    INTEGER, DIMENSION(jpi, jpj) :: ibld
    INTEGER, DIMENSION(jpi, jpj) :: imld
    REAL(KIND = wp) :: ztgrad, zsgrad, zbgrad
    REAL(KIND = wp) :: zugrad, zvgrad
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhbl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhml
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdh
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdhdt
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zt_bl, zs_bl, zu_bl, zv_bl, zrh_bl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zt_ml, zs_ml, zu_ml, zv_ml, zrh_ml
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdt_bl, zds_bl, zdu_bl, zdv_bl, zdrh_bl, zdb_bl
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdt_ml, zds_ml, zdu_ml, zdv_ml, zdrh_ml, zdb_ml
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwth_ent, zws_ent
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zuw_bse, zvw_bse
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdtdz_pyc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdsdz_pyc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdbdz_pyc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdudz_pyc
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdvdz_pyc
    REAL(KIND = wp) :: zl_c, zl_l, zl_eps
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdifml_sc, zvisml_sc, zdifpyc_sc, zvispyc_sc, zbeta_d_sc, zbeta_v_sc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsc_wth_1, zsc_ws_1
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsc_uw_1, zsc_uw_2, zsc_vw_1, zsc_vw_2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhbl_t
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: z3du
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: z3dv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zrimix
    INTEGER :: inhml
    INTEGER :: i_lconv_alloc
    REAL(KIND = wp) :: znd, znd_d, zznd_ml, zznd_pyc, zznd_d
    REAL(KIND = wp) :: ztemp, zari, zpert, zzdhdt, zdb
    REAL(KIND = wp) :: zthick, zz0, zz1
    REAL(KIND = wp) :: zvel_max, zhbl_s
    REAL(KIND = wp) :: zfac
    REAL(KIND = wp) :: zus_x, zus_y
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zviscos
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdiffut
    INTEGER :: ikt
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data5
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data6
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data7
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data8
    CALL profile_psy_data0 % PreStart('zdf_osm', 'r0', 0, 0)
    ALLOCATE(lconv(jpi, jpj), STAT = i_lconv_alloc)
    IF (i_lconv_alloc /= 0) CALL ctl_warn('zdf_osm: failed to allocate lconv')
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    ibld(:, :) = 0
    imld(:, :) = 0
    zrad0(:, :) = 0._wp
    zradh(:, :) = 0._wp
    zradav(:, :) = 0._wp
    zustar(:, :) = 0._wp
    zwstrl(:, :) = 0._wp
    zvstr(:, :) = 0._wp
    zwstrc(:, :) = 0._wp
    zuw0(:, :) = 0._wp
    zvw0(:, :) = 0._wp
    zwth0(:, :) = 0._wp
    zws0(:, :) = 0._wp
    zwb0(:, :) = 0._wp
    zwthav(:, :) = 0._wp
    zwsav(:, :) = 0._wp
    zwbav(:, :) = 0._wp
    zwb_ent(:, :) = 0._wp
    zustke(:, :) = 0._wp
    zla(:, :) = 0._wp
    zcos_wind(:, :) = 0._wp
    zsin_wind(:, :) = 0._wp
    zhol(:, :) = 0._wp
    lconv(:, :) = .FALSE.
    zhbl(:, :) = 1._wp
    zhml(:, :) = 1._wp
    zdh(:, :) = 1._wp
    zdhdt(:, :) = 0._wp
    zt_bl(:, :) = 0._wp
    zs_bl(:, :) = 0._wp
    zu_bl(:, :) = 0._wp
    zv_bl(:, :) = 0._wp
    zrh_bl(:, :) = 0._wp
    zt_ml(:, :) = 0._wp
    zs_ml(:, :) = 0._wp
    zu_ml(:, :) = 0._wp
    zv_ml(:, :) = 0._wp
    zrh_ml(:, :) = 0._wp
    zdt_bl(:, :) = 0._wp
    zds_bl(:, :) = 0._wp
    zdu_bl(:, :) = 0._wp
    zdv_bl(:, :) = 0._wp
    zdrh_bl(:, :) = 0._wp
    zdb_bl(:, :) = 0._wp
    zdt_ml(:, :) = 0._wp
    zds_ml(:, :) = 0._wp
    zdu_ml(:, :) = 0._wp
    zdv_ml(:, :) = 0._wp
    zdrh_ml(:, :) = 0._wp
    zdb_ml(:, :) = 0._wp
    zwth_ent(:, :) = 0._wp
    zws_ent(:, :) = 0._wp
    zuw_bse(:, :) = 0._wp
    zvw_bse(:, :) = 0._wp
    zdtdz_pyc(:, :, :) = 0._wp
    zdsdz_pyc(:, :, :) = 0._wp
    zdbdz_pyc(:, :, :) = 0._wp
    zdudz_pyc(:, :, :) = 0._wp
    zdvdz_pyc(:, :, :) = 0._wp
    zdifml_sc(:, :) = 0._wp
    zvisml_sc(:, :) = 0._wp
    zdifpyc_sc(:, :) = 0._wp
    zvispyc_sc(:, :) = 0._wp
    zbeta_d_sc(:, :) = 0._wp
    zbeta_v_sc(:, :) = 0._wp
    zsc_wth_1(:, :) = 0._wp
    zsc_ws_1(:, :) = 0._wp
    zsc_uw_1(:, :) = 0._wp
    zsc_uw_2(:, :) = 0._wp
    zsc_vw_1(:, :) = 0._wp
    zsc_vw_2(:, :) = 0._wp
    zhbl_t(:, :) = 0._wp
    zdhdt(:, :) = 0._wp
    zdiffut(:, :, :) = 0._wp
    zviscos(:, :, :) = 0._wp
    ghamt(:, :, :) = 0._wp
    ghams(:, :, :) = 0._wp
    ghamu(:, :, :) = 0._wp
    ghamv(:, :, :) = 0._wp
    zz0 = rn_abs
    zz1 = 1. - rn_abs
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        zrad0(ji, jj) = qsr(ji, jj) * r1_rau0_rcp
        zradh(ji, jj) = zrad0(ji, jj) * (zz0 * EXP(- hbl(ji, jj) / rn_si0) + zz1 * EXP(- hbl(ji, jj) / rn_si1))
        zradav(ji, jj) = zrad0(ji, jj) * (zz0 * (1.0 - EXP(- hbl(ji, jj) / rn_si0)) * rn_si0 + zz1 * (1.0 - EXP(- hbl(ji, jj) / rn_si1)) * rn_si1) / hbl(ji, jj)
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        zthermal = rab_n(ji, jj, 1, jp_tem)
        zbeta = rab_n(ji, jj, 1, jp_sal)
        zwth0(ji, jj) = - qns(ji, jj) * r1_rau0_rcp * tmask(ji, jj, 1)
        zws0(ji, jj) = - ((emp(ji, jj) - rnf(ji, jj)) * tsn(ji, jj, 1, jp_sal) + sfx(ji, jj)) * r1_rau0 * tmask(ji, jj, 1)
        zwb0(ji, jj) = grav * zthermal * zwth0(ji, jj) - grav * zbeta * zws0(ji, jj)
        zwthav(ji, jj) = 0.5 * zwth0(ji, jj) - (0.5 * (zrad0(ji, jj) + zradh(ji, jj)) - zradav(ji, jj))
        zwsav(ji, jj) = 0.5 * zws0(ji, jj)
        zwbav(ji, jj) = grav * zthermal * zwthav(ji, jj) - grav * zbeta * zwsav(ji, jj)
        zuw0(ji, jj) = - utau(ji, jj) * r1_rau0 * tmask(ji, jj, 1)
        zvw0(ji, jj) = - vtau(ji, jj) * r1_rau0 * tmask(ji, jj, 1)
        zustar(ji, jj) = MAX(SQRT(SQRT(zuw0(ji, jj) * zuw0(ji, jj) + zvw0(ji, jj) * zvw0(ji, jj))), 1.0E-8)
        zcos_wind(ji, jj) = - zuw0(ji, jj) / (zustar(ji, jj) * zustar(ji, jj))
        zsin_wind(ji, jj) = - zvw0(ji, jj) / (zustar(ji, jj) * zustar(ji, jj))
      END DO
    END DO
    SELECT CASE (nn_osm_wave)
    CASE (0)
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          zus_x = zcos_wind(ji, jj) * zustar(ji, jj) / 0.3 ** 2
          zus_y = zsin_wind(ji, jj) * zustar(ji, jj) / 0.3 ** 2
          zustke(ji, jj) = MAX(SQRT(zus_x * zus_x + zus_y * zus_y), 1.0E-8)
        END DO
      END DO
    CASE (1)
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          zustke(ji, jj) = MAX(0.016 * wndm(ji, jj), 1.0E-8)
          dstokes(ji, jj) = 0.12 * wndm(ji, jj) ** 2 / grav
        END DO
      END DO
    CASE (2)
      zfac = 2.0_wp * rpi / 16.0_wp
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          zustke(ji, jj) = MAX(1.0 * (zcos_wind(ji, jj) * ut0sd(ji, jj) + zsin_wind(ji, jj) * vt0sd(ji, jj)), zustar(ji, jj) / (0.45 * 0.45))
          dstokes(ji, jj) = MAX(zfac * hsw(ji, jj) * hsw(ji, jj) / (MAX(zustke(ji, jj) * wmp(ji, jj), 1.0E-7)), 5.0E-1)
        END DO
      END DO
    END SELECT
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        zwstrl(ji, jj) = (zustar(ji, jj) * zustar(ji, jj) * zustke(ji, jj)) ** pthird
        zwstrl(ji, jj) = zwstrl(ji, jj) * (1.12 * (1.0 - (1.0 - EXP(- hbl(ji, jj) / dstokes(ji, jj))) * dstokes(ji, jj) / hbl(ji, jj))) ** pthird * (1.0 - EXP(- 15.0 * dstokes(ji, jj) / hbl(ji, jj)))
        zla(ji, jj) = SQRT(zustar(ji, jj) / zwstrl(ji, jj)) ** 3
        zvstr(ji, jj) = (zwstrl(ji, jj) ** 3 + (1.0 - EXP(- 0.5 * zla(ji, jj) ** 2)) * zustar(ji, jj) * zustar(ji, jj) * zustar(ji, jj)) ** pthird
        IF (zla(ji, jj) >= 0.45) zla(ji, jj) = 0.45
        IF (zwbav(ji, jj) > 0.0) THEN
          zwstrc(ji, jj) = (2.0 * zwbav(ji, jj) * 0.9 * hbl(ji, jj)) ** pthird
          zhol(ji, jj) = - 0.9 * hbl(ji, jj) * 2.0 * zwbav(ji, jj) / (zvstr(ji, jj) ** 3 + epsln)
          lconv(ji, jj) = .TRUE.
        ELSE
          zhol(ji, jj) = - hbl(ji, jj) * 2.0 * zwbav(ji, jj) / (zvstr(ji, jj) ** 3 + epsln)
          lconv(ji, jj) = .FALSE.
        END IF
      END DO
    END DO
    hbl(:, :) = MAX(hbl(:, :), gdepw_n(:, :, 3))
    ibld(:, :) = 3
    DO jk = 4, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          IF (hbl(ji, jj) >= gdepw_n(ji, jj, jk)) THEN
            ibld(ji, jj) = MIN(mbkt(ji, jj), jk)
          END IF
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        zthermal = rab_n(ji, jj, 1, jp_tem)
        zbeta = rab_n(ji, jj, 1, jp_sal)
        zt = 0._wp
        zs = 0._wp
        zu = 0._wp
        zv = 0._wp
        zthick = 0._wp
        DO jm = 2, ibld(ji, jj)
          zthick = zthick + e3t_n(ji, jj, jm)
          zt = zt + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_tem)
          zs = zs + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_sal)
          zu = zu + e3t_n(ji, jj, jm) * (ub(ji, jj, jm) + ub(ji - 1, jj, jm)) / MAX(1., umask(ji, jj, jm) + umask(ji - 1, jj, jm))
          zv = zv + e3t_n(ji, jj, jm) * (vb(ji, jj, jm) + vb(ji, jj - 1, jm)) / MAX(1., vmask(ji, jj, jm) + vmask(ji, jj - 1, jm))
        END DO
        zt_bl(ji, jj) = zt / zthick
        zs_bl(ji, jj) = zs / zthick
        zu_bl(ji, jj) = zu / zthick
        zv_bl(ji, jj) = zv / zthick
        zdt_bl(ji, jj) = zt_bl(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_tem)
        zds_bl(ji, jj) = zs_bl(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_sal)
        zdu_bl(ji, jj) = zu_bl(ji, jj) - (ub(ji, jj, ibld(ji, jj)) + ub(ji - 1, jj, ibld(ji, jj))) / MAX(1., umask(ji, jj, ibld(ji, jj)) + umask(ji - 1, jj, ibld(ji, jj)))
        zdv_bl(ji, jj) = zv_bl(ji, jj) - (vb(ji, jj, ibld(ji, jj)) + vb(ji, jj - 1, ibld(ji, jj))) / MAX(1., vmask(ji, jj, ibld(ji, jj)) + vmask(ji, jj - 1, ibld(ji, jj)))
        zdb_bl(ji, jj) = grav * zthermal * zdt_bl(ji, jj) - grav * zbeta * zds_bl(ji, jj)
        IF (lconv(ji, jj)) THEN
          zwb_ent(ji, jj) = - (2.0 * 0.2 * zwbav(ji, jj) + 0.135 * zla(ji, jj) * zwstrl(ji, jj) ** 3 / hbl(ji, jj))
          zvel_max = - (1.0 + 1.0 * (zwstrl(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird * rn_rdt / hbl(ji, jj)) * zwb_ent(ji, jj) / (zwstrl(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird
          zzdhdt = - zwb_ent(ji, jj) / (zvel_max + MAX(zdb_bl(ji, jj), 0.0))
        ELSE
          zzdhdt = 0.32 * (hbli(ji, jj) / hbl(ji, jj) - 1.0) * zwstrl(ji, jj) ** 3 / hbli(ji, jj) + ((0.32 / 3.0) * EXP(- 2.5 * (hbli(ji, jj) / hbl(ji, jj) - 1.0)) - (0.32 / 3.0 - 0.135 * zla(ji, jj)) * EXP(- 12.5 * (hbli(ji, jj) / hbl(ji, jj)))) * zwstrl(ji, jj) ** 3 / hbli(ji, jj)
          zzdhdt = zzdhdt + zwbav(ji, jj)
          IF (zzdhdt < 0._wp) THEN
            zpert = 2.0 * (1.0 + 2.0 * zwstrl(ji, jj) * rn_rdt / hbl(ji, jj)) * zwstrl(ji, jj) ** 2 / hbl(ji, jj)
          ELSE
            zpert = 2.0 * (1.0 + 2.0 * zwstrl(ji, jj) * rn_rdt / hbl(ji, jj)) * zwstrl(ji, jj) ** 2 / hbl(ji, jj) + MAX(zdb_bl(ji, jj), 0.0)
          END IF
          zzdhdt = 2.0 * zzdhdt / zpert
        END IF
        zdhdt(ji, jj) = zzdhdt
      END DO
    END DO
    imld = ibld
    ibld(:, :) = 3
    zhbl_t(:, :) = hbl(:, :) + (zdhdt(:, :) - wn(ji, jj, ibld(ji, jj))) * rn_rdt
    zhbl_t(:, :) = MIN(zhbl_t(:, :), ht_n(:, :))
    zdhdt(:, :) = MIN(zdhdt(:, :), (zhbl_t(:, :) - hbl(:, :)) / rn_rdt + wn(ji, jj, ibld(ji, jj)))
    DO jk = 4, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          IF (zhbl_t(ji, jj) >= gdepw_n(ji, jj, jk)) THEN
            ibld(ji, jj) = MIN(mbkt(ji, jj), jk)
          END IF
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        IF (ibld(ji, jj) - imld(ji, jj) > 1) THEN
          CALL profile_psy_data1 % PreStart('zdf_osm', 'r1', 0, 0)
          zhbl_s = hbl(ji, jj)
          jm = imld(ji, jj)
          zthermal = rab_n(ji, jj, 1, jp_tem)
          zbeta = rab_n(ji, jj, 1, jp_sal)
          CALL profile_psy_data1 % PostEnd
          IF (lconv(ji, jj)) THEN
            CALL profile_psy_data2 % PreStart('zdf_osm', 'r2', 0, 0)
            zvel_max = - (1.0 + 1.0 * (zvstr(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird * rn_rdt / hbl(ji, jj)) * zwb_ent(ji, jj) / (zvstr(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird
            DO jk = imld(ji, jj), ibld(ji, jj)
              zdb = MAX(grav * (zthermal * (zt_bl(ji, jj) - tsn(ji, jj, jm, jp_tem)) - zbeta * (zs_bl(ji, jj) - tsn(ji, jj, jm, jp_sal))), 0.0) + zvel_max
              zhbl_s = zhbl_s + MIN(- zwb_ent(ji, jj) / zdb * rn_rdt / FLOAT(ibld(ji, jj) - imld(ji, jj)), e3w_n(ji, jj, jk))
              zhbl_s = MIN(zhbl_s, ht_n(ji, jj))
              IF (zhbl_s >= gdepw_n(ji, jj, jm + 1)) jm = jm + 1
            END DO
            hbl(ji, jj) = zhbl_s
            ibld(ji, jj) = jm
            hbli(ji, jj) = hbl(ji, jj)
            CALL profile_psy_data2 % PostEnd
          ELSE
            !$ACC KERNELS
            DO jk = imld(ji, jj), ibld(ji, jj)
              zdb = MAX(grav * (zthermal * (zt_bl(ji, jj) - tsn(ji, jj, jm, jp_tem)) - zbeta * (zs_bl(ji, jj) - tsn(ji, jj, jm, jp_sal))), 0.0) + 2.0 * zwstrl(ji, jj) ** 2 / zhbl_s
              zhbl_s = zhbl_s + (0.32 * (hbli(ji, jj) / zhbl_s - 1.0) * zwstrl(ji, jj) ** 3 / hbli(ji, jj) + ((0.32 / 3.0) * EXP(- 2.5 * (hbli(ji, jj) / zhbl_s - 1.0)) - (0.32 / 3.0 - 0.0485) * EXP(- 12.5 * (hbli(ji, jj) / zhbl_s))) * zwstrl(ji, jj) ** 3 / hbli(ji, jj)) / zdb * e3w_n(ji, jj, jk) / zdhdt(ji, jj)
              zhbl_s = MIN(zhbl_s, ht_n(ji, jj))
              IF (zhbl_s >= gdepw_n(ji, jj, jm)) jm = jm + 1
            END DO
            hbl(ji, jj) = MAX(zhbl_s, gdepw_n(ji, jj, 3))
            ibld(ji, jj) = MAX(jm, 3)
            IF (hbl(ji, jj) > hbli(ji, jj)) hbli(ji, jj) = hbl(ji, jj)
            !$ACC END KERNELS
          END IF
        ELSE
          CALL profile_psy_data3 % PreStart('zdf_osm', 'r3', 0, 0)
          hbl(ji, jj) = zhbl_t(ji, jj)
          IF (lconv(ji, jj)) THEN
            hbli(ji, jj) = hbl(ji, jj)
          ELSE
            hbl(ji, jj) = MAX(hbl(ji, jj), gdepw_n(ji, jj, 3))
            IF (hbl(ji, jj) > hbli(ji, jj)) hbli(ji, jj) = hbl(ji, jj)
          END IF
          CALL profile_psy_data3 % PostEnd
        END IF
        CALL profile_psy_data4 % PreStart('zdf_osm', 'r4', 0, 0)
        zhbl(ji, jj) = gdepw_n(ji, jj, ibld(ji, jj))
        CALL profile_psy_data4 % PostEnd
      END DO
    END DO
    !$ACC KERNELS
    dstokes(:, :) = MIN(dstokes(:, :), hbl(:, :) / 3.)
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        zthermal = rab_n(ji, jj, 1, jp_tem)
        zbeta = rab_n(ji, jj, 1, jp_sal)
        zt = 0._wp
        zs = 0._wp
        zu = 0._wp
        zv = 0._wp
        zthick = 0._wp
        DO jm = 2, ibld(ji, jj)
          zthick = zthick + e3t_n(ji, jj, jm)
          zt = zt + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_tem)
          zs = zs + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_sal)
          zu = zu + e3t_n(ji, jj, jm) * (ub(ji, jj, jm) + ub(ji - 1, jj, jm)) / MAX(1., umask(ji, jj, jm) + umask(ji - 1, jj, jm))
          zv = zv + e3t_n(ji, jj, jm) * (vb(ji, jj, jm) + vb(ji, jj - 1, jm)) / MAX(1., vmask(ji, jj, jm) + vmask(ji, jj - 1, jm))
        END DO
        zt_bl(ji, jj) = zt / zthick
        zs_bl(ji, jj) = zs / zthick
        zu_bl(ji, jj) = zu / zthick
        zv_bl(ji, jj) = zv / zthick
        zdt_bl(ji, jj) = zt_bl(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_tem)
        zds_bl(ji, jj) = zs_bl(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_sal)
        zdu_bl(ji, jj) = zu_bl(ji, jj) - (ub(ji, jj, ibld(ji, jj)) + ub(ji - 1, jj, ibld(ji, jj))) / MAX(1., umask(ji, jj, ibld(ji, jj)) + umask(ji - 1, jj, ibld(ji, jj)))
        zdv_bl(ji, jj) = zv_bl(ji, jj) - (vb(ji, jj, ibld(ji, jj)) + vb(ji, jj - 1, ibld(ji, jj))) / MAX(1., vmask(ji, jj, ibld(ji, jj)) + vmask(ji, jj - 1, ibld(ji, jj)))
        zdb_bl(ji, jj) = grav * zthermal * zdt_bl(ji, jj) - grav * zbeta * zds_bl(ji, jj)
        zhbl(ji, jj) = gdepw_n(ji, jj, ibld(ji, jj))
        IF (lconv(ji, jj)) THEN
          IF (zdb_bl(ji, jj) > 0._wp) THEN
            IF ((zwstrc(ji, jj) / zvstr(ji, jj)) ** 3 <= 0.5) THEN
              zari = 4.5 * (zvstr(ji, jj) ** 2) / (zdb_bl(ji, jj) * zhbl(ji, jj)) + 0.01
            ELSE
              zari = 4.5 * (zwstrc(ji, jj) ** 2) / (zdb_bl(ji, jj) * zhbl(ji, jj)) + 0.01
            END IF
            IF (zari > 0.2) THEN
              zari = 0.2
              zwb_ent(ji, jj) = 0._wp
            END IF
            inhml = MAX(INT(zari * zhbl(ji, jj) / e3t_n(ji, jj, ibld(ji, jj))), 1)
            imld(ji, jj) = MAX(ibld(ji, jj) - inhml, 1)
            zhml(ji, jj) = gdepw_n(ji, jj, imld(ji, jj))
            zdh(ji, jj) = zhbl(ji, jj) - zhml(ji, jj)
          ELSE
            imld(ji, jj) = ibld(ji, jj) - 1
            zhml(ji, jj) = gdepw_n(ji, jj, imld(ji, jj))
            zdh(ji, jj) = zhbl(ji, jj) - zhml(ji, jj)
          END IF
        ELSE
          IF (zdhdt(ji, jj) >= 0.0) THEN
            IF (zdb_bl(ji, jj) > 0._wp) THEN
              zari = MIN(4.5 * (zvstr(ji, jj) ** 2) / (zdb_bl(ji, jj) * zhbl(ji, jj)) + 0.01, 0.2)
              inhml = MAX(INT(zari * zhbl(ji, jj) / e3t_n(ji, jj, ibld(ji, jj))), 1)
              imld(ji, jj) = MAX(ibld(ji, jj) - inhml, 1)
              zhml(ji, jj) = gdepw_n(ji, jj, imld(ji, jj))
              zdh(ji, jj) = zhbl(ji, jj) - zhml(ji, jj)
            ELSE
              imld(ji, jj) = ibld(ji, jj) - 1
              zhml(ji, jj) = gdepw_n(ji, jj, imld(ji, jj))
              zdh(ji, jj) = zhbl(ji, jj) - zhml(ji, jj)
            END IF
          ELSE
            imld(ji, jj) = ibld(ji, jj)
            zhml(ji, jj) = zhbl(ji, jj)
            zdh(ji, jj) = 0._wp
          END IF
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        zthermal = rab_n(ji, jj, 1, jp_tem)
        zbeta = rab_n(ji, jj, 1, jp_sal)
        IF (lconv(ji, jj)) THEN
          zt = 0._wp
          zs = 0._wp
          zu = 0._wp
          zv = 0._wp
          zthick = 0._wp
          DO jm = 2, imld(ji, jj)
            zthick = zthick + e3t_n(ji, jj, jm)
            zt = zt + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_tem)
            zs = zs + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_sal)
            zu = zu + e3t_n(ji, jj, jm) * (ub(ji, jj, jm) + ub(ji - 1, jj, jm)) / MAX(1., umask(ji, jj, jm) + umask(ji - 1, jj, jm))
            zv = zv + e3t_n(ji, jj, jm) * (vb(ji, jj, jm) + vb(ji, jj - 1, jm)) / MAX(1., vmask(ji, jj, jm) + vmask(ji, jj - 1, jm))
          END DO
          zt_ml(ji, jj) = zt / zthick
          zs_ml(ji, jj) = zs / zthick
          zu_ml(ji, jj) = zu / zthick
          zv_ml(ji, jj) = zv / zthick
          zdt_ml(ji, jj) = zt_ml(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_tem)
          zds_ml(ji, jj) = zs_ml(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_sal)
          zdu_ml(ji, jj) = zu_ml(ji, jj) - (ub(ji, jj, ibld(ji, jj)) + ub(ji - 1, jj, ibld(ji, jj))) / MAX(1., umask(ji, jj, ibld(ji, jj)) + umask(ji - 1, jj, ibld(ji, jj)))
          zdv_ml(ji, jj) = zv_ml(ji, jj) - (vb(ji, jj, ibld(ji, jj)) + vb(ji, jj - 1, ibld(ji, jj))) / MAX(1., vmask(ji, jj, ibld(ji, jj)) + vmask(ji, jj - 1, ibld(ji, jj)))
          zdb_ml(ji, jj) = grav * zthermal * zdt_ml(ji, jj) - grav * zbeta * zds_ml(ji, jj)
        ELSE
          IF (zdhdt(ji, jj) >= 0._wp) THEN
            zt = 0._wp
            zs = 0._wp
            zu = 0._wp
            zv = 0._wp
            zthick = 0._wp
            DO jm = 2, imld(ji, jj)
              zthick = zthick + e3t_n(ji, jj, jm)
              zt = zt + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_tem)
              zs = zs + e3t_n(ji, jj, jm) * tsn(ji, jj, jm, jp_sal)
              zu = zu + e3t_n(ji, jj, jm) * (ub(ji, jj, jm) + ub(ji - 1, jj, jm)) / MAX(1., umask(ji, jj, jm) + umask(ji - 1, jj, jm))
              zv = zv + e3t_n(ji, jj, jm) * (vb(ji, jj, jm) + vb(ji, jj - 1, jm)) / MAX(1., vmask(ji, jj, jm) + vmask(ji, jj - 1, jm))
            END DO
            zt_ml(ji, jj) = zt / zthick
            zs_ml(ji, jj) = zs / zthick
            zu_ml(ji, jj) = zu / zthick
            zv_ml(ji, jj) = zv / zthick
            zdt_ml(ji, jj) = zt_ml(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_tem)
            zds_ml(ji, jj) = zs_ml(ji, jj) - tsn(ji, jj, ibld(ji, jj), jp_sal)
            zdu_ml(ji, jj) = zu_ml(ji, jj) - (ub(ji, jj, ibld(ji, jj)) + ub(ji - 1, jj, ibld(ji, jj))) / MAX(1., umask(ji, jj, ibld(ji, jj)) + umask(ji - 1, jj, ibld(ji, jj)))
            zdv_ml(ji, jj) = zv_ml(ji, jj) - (vb(ji, jj, ibld(ji, jj)) + vb(ji, jj - 1, ibld(ji, jj))) / MAX(1., vmask(ji, jj, ibld(ji, jj)) + vmask(ji, jj - 1, ibld(ji, jj)))
            zdb_ml(ji, jj) = grav * zthermal * zdt_ml(ji, jj) - grav * zbeta * zds_ml(ji, jj)
          END IF
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        ztemp = zu_ml(ji, jj)
        zu_ml(ji, jj) = zu_ml(ji, jj) * zcos_wind(ji, jj) + zv_ml(ji, jj) * zsin_wind(ji, jj)
        zv_ml(ji, jj) = zv_ml(ji, jj) * zcos_wind(ji, jj) - ztemp * zsin_wind(ji, jj)
        ztemp = zdu_ml(ji, jj)
        zdu_ml(ji, jj) = zdu_ml(ji, jj) * zcos_wind(ji, jj) + zdv_ml(ji, jj) * zsin_wind(ji, jj)
        zdv_ml(ji, jj) = zdv_ml(ji, jj) * zsin_wind(ji, jj) - ztemp * zsin_wind(ji, jj)
        ztemp = zu_bl(ji, jj)
        zu_bl = zu_bl(ji, jj) * zcos_wind(ji, jj) + zv_bl(ji, jj) * zsin_wind(ji, jj)
        zv_bl(ji, jj) = zv_bl(ji, jj) * zcos_wind(ji, jj) - ztemp * zsin_wind(ji, jj)
        ztemp = zdu_bl(ji, jj)
        zdu_bl(ji, jj) = zdu_bl(ji, jj) * zcos_wind(ji, jj) + zdv_bl(ji, jj) * zsin_wind(ji, jj)
        zdv_bl(ji, jj) = zdv_bl(ji, jj) * zsin_wind(ji, jj) - ztemp * zsin_wind(ji, jj)
      END DO
    END DO
    zuw_bse = 0._wp
    zvw_bse = 0._wp
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          IF (zdb_bl(ji, jj) > 0._wp) THEN
            zwth_ent(ji, jj) = zwb_ent(ji, jj) * zdt_ml(ji, jj) / (zdb_ml(ji, jj) + epsln)
            zws_ent(ji, jj) = zwb_ent(ji, jj) * zds_ml(ji, jj) / (zdb_ml(ji, jj) + epsln)
          END IF
        ELSE
          zwth_ent(ji, jj) = - 2.0 * zwthav(ji, jj) * ((1.0 - 0.8) - (1.0 - 0.8) ** (3.0 / 2.0))
          zws_ent(ji, jj) = - 2.0 * zwsav(ji, jj) * ((1.0 - 0.8) - (1.0 - 0.8) ** (3.0 / 2.0))
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          IF (zdb_bl(ji, jj) > 0._wp) THEN
            ztgrad = (zdt_ml(ji, jj) / zdh(ji, jj))
            zsgrad = (zds_ml(ji, jj) / zdh(ji, jj))
            zbgrad = (zdb_ml(ji, jj) / zdh(ji, jj))
            DO jk = 2, ibld(ji, jj)
              znd = - (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zdh(ji, jj)
              zdtdz_pyc(ji, jj, jk) = ztgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
              zdbdz_pyc(ji, jj, jk) = zbgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
              zdsdz_pyc(ji, jj, jk) = zsgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
            END DO
          END IF
        ELSE
          IF (zdhdt(ji, jj) >= 0.0) THEN
            IF (zdb_bl(ji, jj) > 0._wp) THEN
              IF (zhol(ji, jj) >= 0.5) THEN
                ztgrad = zdt_bl(ji, jj) / zhbl(ji, jj)
                zsgrad = zds_bl(ji, jj) / zhbl(ji, jj)
                zbgrad = zdb_bl(ji, jj) / zhbl(ji, jj)
                DO jk = 2, ibld(ji, jj)
                  znd = gdepw_n(ji, jj, jk) / zhbl(ji, jj)
                  zdtdz_pyc(ji, jj, jk) = ztgrad * EXP(- 15.0 * (znd - 0.9) ** 2)
                  zdbdz_pyc(ji, jj, jk) = zbgrad * EXP(- 15.0 * (znd - 0.9) ** 2)
                  zdsdz_pyc(ji, jj, jk) = zsgrad * EXP(- 15.0 * (znd - 0.9) ** 2)
                END DO
              ELSE
                ztgrad = zdt_bl(ji, jj) / zdh(ji, jj)
                zsgrad = zds_bl(ji, jj) / zdh(ji, jj)
                zbgrad = zdb_bl(ji, jj) / zdh(ji, jj)
                DO jk = 2, ibld(ji, jj)
                  znd = - (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zdh(ji, jj)
                  zdtdz_pyc(ji, jj, jk) = ztgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
                  zdbdz_pyc(ji, jj, jk) = zbgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
                  zdsdz_pyc(ji, jj, jk) = zsgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
                END DO
              END IF
            END IF
          END IF
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          zugrad = (zdu_ml(ji, jj) / zdh(ji, jj)) + 0.275 * zustar(ji, jj) * zustar(ji, jj) / ((zwstrl(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird * zhml(ji, jj)) / zla(ji, jj) ** (8.0 / 3.0)
          zvgrad = (zdv_ml(ji, jj) / zdh(ji, jj)) + 3.5 * ff_t(ji, jj) * zustke(ji, jj) / (zwstrl(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird
          DO jk = 2, ibld(ji, jj) - 1
            znd = - (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zdh(ji, jj)
            zdudz_pyc(ji, jj, jk) = zugrad * EXP(- 1.75 * (znd + 0.75) ** 2)
            zdvdz_pyc(ji, jj, jk) = zvgrad * EXP(- 1.75 * (znd + 0.75) ** 2)
          END DO
        ELSE
          zugrad = 3.25 * zdu_bl(ji, jj) / zhbl(ji, jj)
          zvgrad = 2.75 * zdv_bl(ji, jj) / zhbl(ji, jj)
          DO jk = 2, ibld(ji, jj)
            znd = gdepw_n(ji, jj, jk) / zhbl(ji, jj)
            IF (znd < 1.0) THEN
              zdudz_pyc(ji, jj, jk) = zugrad * EXP(- 40.0 * (znd - 1.0) ** 2)
            ELSE
              zdudz_pyc(ji, jj, jk) = zugrad * EXP(- 20.0 * (znd - 1.0) ** 2)
            END IF
            zdvdz_pyc(ji, jj, jk) = zvgrad * EXP(- 20.0 * (znd - 0.85) ** 2)
          END DO
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          zdifml_sc(ji, jj) = zhml(ji, jj) * (zvstr(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird
          zvisml_sc(ji, jj) = zdifml_sc(ji, jj)
          zdifpyc_sc(ji, jj) = 0.165 * (zvstr(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird * zdh(ji, jj)
          zvispyc_sc(ji, jj) = 0.142 * (zvstr(ji, jj) ** 3 + 0.5 * zwstrc(ji, jj) ** 3) ** pthird * zdh(ji, jj)
          zbeta_d_sc(ji, jj) = 1.0 - (0.165 / 0.8 * zdh(ji, jj) / zhbl(ji, jj)) ** p2third
          zbeta_v_sc(ji, jj) = 1.0 - 2.0 * (0.142 / 0.375) * zdh(ji, jj) / zhml(ji, jj)
        ELSE
          zdifml_sc(ji, jj) = zvstr(ji, jj) * zhbl(ji, jj) * EXP(- (zhol(ji, jj) / 0.6_wp) ** 2)
          zvisml_sc(ji, jj) = zvstr(ji, jj) * zhbl(ji, jj) * EXP(- (zhol(ji, jj) / 0.6_wp) ** 2)
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_ml = gdepw_n(ji, jj, jk) / zhml(ji, jj)
            zdiffut(ji, jj, jk) = 0.8 * zdifml_sc(ji, jj) * zznd_ml * (1.0 - zbeta_d_sc(ji, jj) * zznd_ml) ** 1.5
            zviscos(ji, jj, jk) = 0.375 * zvisml_sc(ji, jj) * zznd_ml * (1.0 - zbeta_v_sc(ji, jj) * zznd_ml) * (1.0 - 0.5 * zznd_ml ** 2)
          END DO
          IF (zdh(ji, jj) > 0._wp) THEN
            DO jk = imld(ji, jj) + 1, ibld(ji, jj)
              zznd_pyc = - (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zdh(ji, jj)
              zdiffut(ji, jj, jk) = zdifpyc_sc(ji, jj) * (1.0 + zznd_pyc)
              zviscos(ji, jj, jk) = zvispyc_sc(ji, jj) * (1.0 + zznd_pyc)
            END DO
          END IF
          zdiffut(ji, jj, ibld(ji, jj)) = zdhdt(ji, jj) * e3t_n(ji, jj, ibld(ji, jj))
        ELSE
          DO jk = 2, ibld(ji, jj)
            zznd_ml = gdepw_n(ji, jj, jk) / zhbl(ji, jj)
            zdiffut(ji, jj, jk) = 0.75 * zdifml_sc(ji, jj) * zznd_ml * (1.0 - zznd_ml) ** 1.5
            zviscos(ji, jj, jk) = 0.375 * zvisml_sc(ji, jj) * zznd_ml * (1.0 - zznd_ml) * (1.0 - zznd_ml ** 2)
          END DO
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (lconv)
      zsc_wth_1 = zwstrl ** 3 * zwth0 / (zvstr ** 3 + 0.5 * zwstrc ** 3 + epsln)
      zsc_ws_1 = zwstrl ** 3 * zws0 / (zvstr ** 3 + 0.5 * zwstrc ** 3 + epsln)
    ELSEWHERE
      zsc_wth_1 = 2.0 * zwthav
      zsc_ws_1 = 2.0 * zwsav
    END WHERE
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + 1.35 * EXP(- zznd_d) * (1.0 - EXP(- 2.0 * zznd_d)) * zsc_wth_1(ji, jj)
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + 1.35 * EXP(- zznd_d) * (1.0 - EXP(- 2.0 * zznd_d)) * zsc_ws_1(ji, jj)
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + 1.5 * EXP(- 0.9 * zznd_d) * (1.0 - EXP(- 4.0 * zznd_d)) * zsc_wth_1(ji, jj)
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + 1.5 * EXP(- 0.9 * zznd_d) * (1.0 - EXP(- 4.0 * zznd_d)) * zsc_ws_1(ji, jj)
          END DO
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (lconv)
      zsc_uw_1 = (zwstrl ** 3 + 0.5 * zwstrc ** 3) ** pthird * zustke / (1.0 - 1.0 * 6.5 * zla ** (8.0 / 3.0))
      zsc_uw_2 = (zwstrl ** 3 + 0.5 * zwstrc ** 3) ** pthird * zustke / (zla ** (8.0 / 3.0) + epsln)
      zsc_vw_1 = ff_t * zhml * zustke ** 3 * zla ** (8.0 / 3.0) / ((zvstr ** 3 + 0.5 * zwstrc ** 3) ** (2.0 / 3.0) + epsln)
    ELSEWHERE
      zsc_uw_1 = zustar ** 2
      zsc_vw_1 = ff_t * zhbl * zustke ** 3 * zla ** (8.0 / 3.0) / (zvstr ** 2 + epsln)
    END WHERE
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + (- 0.05 * EXP(- 0.4 * zznd_d) * zsc_uw_1(ji, jj) + 0.00125 * EXP(- zznd_d) * zsc_uw_2(ji, jj)) * (1.0 - EXP(- 2.0 * zznd_d))
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) - 0.65 * 0.15 * EXP(- zznd_d) * (1.0 - EXP(- 2.0 * zznd_d)) * zsc_vw_1(ji, jj)
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) - 0.75 * 1.3 * EXP(- 0.5 * zznd_d) * (1.0 - EXP(- 4.0 * zznd_d)) * zsc_uw_1(ji, jj)
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + 0._wp
          END DO
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (lconv)
      zsc_wth_1 = zwbav * zwth0 * (1.0 + EXP(0.2 * zhol)) / (zvstr ** 3 + 0.5 * zwstrc ** 3 + epsln)
      zsc_ws_1 = zwbav * zws0 * (1.0 + EXP(0.2 * zhol)) / (zvstr ** 3 + 0.5 * zwstrc ** 3 + epsln)
    ELSEWHERE
      zsc_wth_1 = 0._wp
      zsc_ws_1 = 0._wp
    END WHERE
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_ml = gdepw_n(ji, jj, jk) / zhml(ji, jj)
            zl_c = 0.9 * (1.0 - EXP(- 7.0 * (zznd_ml - zznd_ml ** 3 / 3.0))) * (1.0 - EXP(- 15.0 * (1.1 - zznd_ml)))
            zl_l = 2.0 * (1.0 - EXP(- 2.0 * (zznd_ml - zznd_ml ** 3 / 3.0))) * (1.0 - EXP(- 5.0 * (1.0 - zznd_ml))) * (1.0 + dstokes(ji, jj) / zhml(ji, jj))
            zl_eps = zl_l + (zl_c - zl_l) / (1.0 + EXP(3.0 * LOG10(- zhol(ji, jj)))) ** (3.0 / 2.0)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + 0.3 * 0.5 * zsc_wth_1(ji, jj) * zl_eps * zhml(ji, jj) / (0.15 + zznd_ml)
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + 0.3 * 0.5 * zsc_ws_1(ji, jj) * zl_eps * zhml(ji, jj) / (0.15 + zznd_ml)
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + zsc_wth_1(ji, jj)
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + zsc_ws_1(ji, jj)
          END DO
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (lconv)
      zsc_uw_1 = - zwb0 * zustar ** 2 * zhml / (zvstr ** 3 + 0.5 * zwstrc ** 3 + epsln)
      zsc_uw_2 = zwb0 * zustke * zhml / (zvstr ** 3 + 0.5 * zwstrc ** 3 + epsln) ** (2.0 / 3.0)
      zsc_vw_1 = 0._wp
    ELSEWHERE
      zsc_uw_1 = 0._wp
      zsc_vw_1 = 0._wp
    END WHERE
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT !CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + 0.3 * 0.5 * (zsc_uw_1(ji, jj) + 0.125 * EXP(- 0.5 * zznd_d) * (1.0 - EXP(- 0.5 * zznd_d)) * zsc_uw_2(ji, jj))
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + zsc_vw_1(ji, jj)
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + zsc_uw_1(ji, jj)
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + zsc_vw_1(ji, jj)
          END DO
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (lconv)
      zsc_wth_1 = zwth0
      zsc_ws_1 = zws0
    ELSEWHERE
      zsc_wth_1 = 2.0 * zwthav
      zsc_ws_1 = zws0
    END WHERE
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_ml = gdepw_n(ji, jj, jk) / zhml(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + 0.3 * zsc_wth_1(ji, jj) * (- 2.0 + 2.75 * ((1.0 + 0.6 * zznd_ml ** 4) - EXP(- 6.0 * zznd_ml))) * (1.0 - EXP(- 15.0 * (1.0 - zznd_ml)))
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + 0.3 * zsc_ws_1(ji, jj) * (- 2.0 + 2.75 * ((1.0 + 0.6 * zznd_ml ** 4) - EXP(- 6.0 * zznd_ml))) * (1.0 - EXP(- 15.0 * (1.0 - zznd_ml)))
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            znd = gdepw_n(ji, jj, jk) / zhbl(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + 0.3 * (- 4.06 * EXP(- 2.0 * zznd_d) * (1.0 - EXP(- 4.0 * zznd_d)) + 7.5 * EXP(- 10.0 * (0.95 - znd) ** 2) * (1.0 - znd)) * zsc_wth_1(ji, jj)
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + 0.3 * (- 4.06 * EXP(- 2.0 * zznd_d) * (1.0 - EXP(- 4.0 * zznd_d)) + 7.5 * EXP(- 10.0 * (0.95 - znd) ** 2) * (1.0 - znd)) * zsc_ws_1(ji, jj)
          END DO
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    WHERE (lconv)
      zsc_uw_1 = zustar ** 2
      zsc_vw_1 = ff_t * zustke * zhml
    ELSEWHERE
      zsc_uw_1 = zustar ** 2
      zsc_uw_2 = (2.25 - 3.0 * (1.0 - EXP(- 1.25 * 2.0))) * (1.0 - EXP(- 4.0 * 2.0)) * zsc_uw_1
      zsc_vw_1 = ff_t * zustke * zhbl
      zsc_vw_2 = - 0.11 * SIN(3.14159 * (2.0 + 0.4)) * EXP(- (1.5 + 2.0) ** 2) * zsc_vw_1
    END WHERE
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, imld(ji, jj)
            zznd_ml = gdepw_n(ji, jj, jk) / zhml(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + 0.3 * (- 2.0 + 2.5 * (1.0 + 0.1 * zznd_ml ** 4) - EXP(- 8.0 * zznd_ml)) * zsc_uw_1(ji, jj)
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + 0.3 * 0.1 * (EXP(- zznd_d) + EXP(- 5.0 * (1.0 - zznd_ml))) * zsc_vw_1(ji, jj)
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            znd = gdepw_n(ji, jj, jk) / zhbl(ji, jj)
            zznd_d = gdepw_n(ji, jj, jk) / dstokes(ji, jj)
            IF (zznd_d <= 2.0) THEN
              ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + 0.5 * 0.3 * (2.25 - 3.0 * (1.0 - EXP(- 1.25 * zznd_d)) * (1.0 - EXP(- 2.0 * zznd_d))) * zsc_uw_1(ji, jj)
            ELSE
              ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + 0.5 * 0.3 * (1.0 - EXP(- 5.0 * (1.0 - znd))) * zsc_uw_2(ji, jj)
            END IF
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + 0.3 * 0.15 * SIN(3.14159 * (0.65 * zznd_d)) * EXP(- 0.25 * zznd_d ** 2) * zsc_vw_1(ji, jj)
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + 0.3 * 0.15 * EXP(- 5.0 * (1.0 - znd)) * (1.0 - EXP(- 20.0 * (1.0 - znd))) * zsc_vw_2(ji, jj)
          END DO
        END IF
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 2, ibld(ji, jj)
            znd = (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zhml(ji, jj)
            IF (znd >= 0.0) THEN
              ghamu(ji, jj, jk) = ghamu(ji, jj, jk) * (1.0 - EXP(- 30.0 * znd ** 2))
              ghamv(ji, jj, jk) = ghamv(ji, jj, jk) * (1.0 - EXP(- 30.0 * znd ** 2))
            ELSE
              ghamu(ji, jj, jk) = 0._wp
              ghamv(ji, jj, jk) = 0._wp
            END IF
          END DO
        ELSE
          DO jk = 2, ibld(ji, jj)
            znd = (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zhml(ji, jj)
            IF (znd >= 0.0) THEN
              ghamu(ji, jj, jk) = ghamu(ji, jj, jk) * (1.0 - EXP(- 10.0 * znd ** 2))
              ghamv(ji, jj, jk) = ghamv(ji, jj, jk) * (1.0 - EXP(- 10.0 * znd ** 2))
            ELSE
              ghamu(ji, jj, jk) = 0._wp
              ghamv(ji, jj, jk) = 0._wp
            END IF
          END DO
        END IF
      END DO
    END DO
    zsc_uw_1 = 0._wp
    !$ACC END KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        !$ACC KERNELS
        DO jk = 2, ibld(ji, jj)
          znd = gdepw_n(ji, jj, jk) / zhbl(ji, jj)
          ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + zdiffut(ji, jj, jk) * zdtdz_pyc(ji, jj, jk)
          ghams(ji, jj, jk) = ghams(ji, jj, jk) + zdiffut(ji, jj, jk) * zdsdz_pyc(ji, jj, jk)
          ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + zviscos(ji, jj, jk) * zdudz_pyc(ji, jj, jk)
          ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + zsc_uw_1(ji, jj) * (1.0 - znd) ** (7.0 / 4.0) * zdbdz_pyc(ji, jj, jk)
          ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + zviscos(ji, jj, jk) * zdvdz_pyc(ji, jj, jk)
        END DO
        !$ACC END KERNELS
      END DO
    END DO
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 2, jpjm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 2, jpim1
        IF (lconv(ji, jj)) THEN
          DO jk = 1, imld(ji, jj) - 1
            znd = gdepw_n(ji, jj, jk) / zhml(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + zwth_ent(ji, jj) * znd
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + zws_ent(ji, jj) * znd
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + zuw_bse(ji, jj) * znd
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + zvw_bse(ji, jj) * znd
          END DO
          DO jk = imld(ji, jj), ibld(ji, jj)
            znd = - (gdepw_n(ji, jj, jk) - zhml(ji, jj)) / zdh(ji, jj)
            ghamt(ji, jj, jk) = ghamt(ji, jj, jk) + zwth_ent(ji, jj) * (1.0 + znd)
            ghams(ji, jj, jk) = ghams(ji, jj, jk) + zws_ent(ji, jj) * (1.0 + znd)
            ghamu(ji, jj, jk) = ghamu(ji, jj, jk) + zuw_bse(ji, jj) * (1.0 + znd)
            ghamv(ji, jj, jk) = ghamv(ji, jj, jk) + zvw_bse(ji, jj) * (1.0 + znd)
          END DO
        END IF
        ghamt(ji, jj, ibld(ji, jj)) = 0._wp
        ghams(ji, jj, ibld(ji, jj)) = 0._wp
        ghamu(ji, jj, ibld(ji, jj)) = 0._wp
        ghamv(ji, jj, ibld(ji, jj)) = 0._wp
      END DO
    END DO
    !$ACC END KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        !$ACC KERNELS
        DO jk = 2, ibld(ji, jj)
          ztemp = ghamu(ji, jj, jk)
          ghamu(ji, jj, jk) = ghamu(ji, jj, jk) * zcos_wind(ji, jj) - ghamv(ji, jj, jk) * zsin_wind(ji, jj)
          ghamv(ji, jj, jk) = ghamv(ji, jj, jk) * zcos_wind(ji, jj) + ztemp * zsin_wind(ji, jj)
        END DO
        !$ACC END KERNELS
      END DO
    END DO
    CALL profile_psy_data5 % PreStart('zdf_osm', 'r5', 0, 0)
    IF (ln_dia_osm) THEN
      IF (iom_use("zdtdz_pyc")) CALL iom_put("zdtdz_pyc", wmask * zdtdz_pyc)
    END IF
    CALL profile_psy_data5 % PostEnd
    IF (ln_kpprimix) THEN
      !$ACC KERNELS
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = 1, jpjm1
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = 1, jpim1
            z3du(ji, jj, jk) = 0.5 * (un(ji, jj, jk - 1) - un(ji, jj, jk)) * (ub(ji, jj, jk - 1) - ub(ji, jj, jk)) * wumask(ji, jj, jk) / (e3uw_n(ji, jj, jk) * e3uw_b(ji, jj, jk))
            z3dv(ji, jj, jk) = 0.5 * (vn(ji, jj, jk - 1) - vn(ji, jj, jk)) * (vb(ji, jj, jk - 1) - vb(ji, jj, jk)) * wvmask(ji, jj, jk) / (e3vw_n(ji, jj, jk) * e3vw_b(ji, jj, jk))
          END DO
        END DO
      END DO
      DO jk = 2, jpkm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = 2, jpjm1
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = 2, jpim1
            zesh2 = (z3du(ji - 1, jj, jk) + z3du(ji, jj, jk)) / MAX(1._wp, umask(ji - 1, jj, jk) + umask(ji, jj, jk)) + (z3dv(ji, jj - 1, jk) + z3dv(ji, jj, jk)) / MAX(1._wp, vmask(ji, jj - 1, jk) + vmask(ji, jj, jk))
            zri = MAX(rn2b(ji, jj, jk), 0._wp) / MAX(zesh2, epsln)
            zfri = MIN(zri / rn_riinfty, 1.0_wp)
            zfri = (1.0_wp - zfri * zfri)
            zrimix(ji, jj, jk) = zfri * zfri * zfri * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          !$ACC KERNELS
          DO jk = ibld(ji, jj) + 1, jpkm1
            zdiffut(ji, jj, jk) = zrimix(ji, jj, jk) * rn_difri
            zviscos(ji, jj, jk) = zrimix(ji, jj, jk) * rn_difri
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END IF
    IF (ln_convmix) THEN
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          !$ACC KERNELS
          DO jk = ibld(ji, jj) + 1, jpkm1
            IF (MIN(rn2(ji, jj, jk), rn2b(ji, jj, jk)) <= - 1.E-12) zdiffut(ji, jj, jk) = rn_difconv
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END IF
    CALL profile_psy_data6 % PreStart('zdf_osm', 'r6', 0, 0)
    CALL lbc_lnk('zdfosm', zviscos(:, :, :), 'W', 1.)
    CALL profile_psy_data6 % PostEnd
    !$ACC KERNELS
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          p_avt(ji, jj, jk) = MAX(zdiffut(ji, jj, jk), avtb(jk)) * tmask(ji, jj, jk)
          p_avm(ji, jj, jk) = MAX(zviscos(ji, jj, jk), avmb(jk)) * tmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data7 % PreStart('zdf_osm', 'r7', 0, 0)
    CALL lbc_lnk_multi('zdfosm', p_avt, 'W', 1., p_avm, 'W', 1., ghamu, 'W', 1., ghamv, 'W', 1.)
    CALL profile_psy_data7 % PostEnd
    !$ACC KERNELS
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          ghamu(ji, jj, jk) = (ghamu(ji, jj, jk) + ghamu(ji + 1, jj, jk)) / MAX(1., tmask(ji, jj, jk) + tmask(ji + 1, jj, jk)) * umask(ji, jj, jk)
          ghamv(ji, jj, jk) = (ghamv(ji, jj, jk) + ghamv(ji, jj + 1, jk)) / MAX(1., tmask(ji, jj, jk) + tmask(ji, jj + 1, jk)) * vmask(ji, jj, jk)
          ghamt(ji, jj, jk) = ghamt(ji, jj, jk) * tmask(ji, jj, jk)
          ghams(ji, jj, jk) = ghams(ji, jj, jk) * tmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data8 % PreStart('zdf_osm', 'r8', 0, 0)
    CALL lbc_lnk_multi('zdfosm', ghamt, 'W', 1., ghams, 'W', 1., ghamu, 'U', 1., ghamv, 'V', 1.)
    IF (ln_dia_osm) THEN
      SELECT CASE (nn_osm_wave)
      CASE (0 : 1)
        IF (iom_use("us_x")) CALL iom_put("us_x", tmask(:, :, 1) * zustke * zcos_wind)
        IF (iom_use("us_y")) CALL iom_put("us_y", tmask(:, :, 1) * zustke * zsin_wind)
        IF (iom_use("wind_wave_abs_power")) CALL iom_put("wind_wave_abs_power", 1000. * rau0 * tmask(:, :, 1) * zustar ** 2 * zustke)
      CASE (2)
        IF (iom_use("us_x")) CALL iom_put("us_x", ut0sd)
        IF (iom_use("us_y")) CALL iom_put("us_y", vt0sd)
        IF (iom_use("wind_wave_abs_power")) CALL iom_put("wind_wave_abs_power", 1000. * rau0 * tmask(:, :, 1) * zustar ** 2 * SQRT(ut0sd ** 2 + vt0sd ** 2))
      END SELECT
      IF (iom_use("ghamt")) CALL iom_put("ghamt", tmask * ghamt)
      IF (iom_use("ghams")) CALL iom_put("ghams", tmask * ghams)
      IF (iom_use("ghamu")) CALL iom_put("ghamu", umask * ghamu)
      IF (iom_use("ghamv")) CALL iom_put("ghamv", vmask * ghamv)
      IF (iom_use("zwth0")) CALL iom_put("zwth0", tmask(:, :, 1) * zwth0)
      IF (iom_use("zws0")) CALL iom_put("zws0", tmask(:, :, 1) * zws0)
      IF (iom_use("hbl")) CALL iom_put("hbl", tmask(:, :, 1) * hbl)
      IF (iom_use("hbli")) CALL iom_put("hbli", tmask(:, :, 1) * hbli)
      IF (iom_use("dstokes")) CALL iom_put("dstokes", tmask(:, :, 1) * dstokes)
      IF (iom_use("zustke")) CALL iom_put("zustke", tmask(:, :, 1) * zustke)
      IF (iom_use("zwstrc")) CALL iom_put("zwstrc", tmask(:, :, 1) * zwstrc)
      IF (iom_use("zwstrl")) CALL iom_put("zwstrl", tmask(:, :, 1) * zwstrl)
      IF (iom_use("zustar")) CALL iom_put("zustar", tmask(:, :, 1) * zustar)
      IF (iom_use("wind_power")) CALL iom_put("wind_power", 1000. * rau0 * tmask(:, :, 1) * zustar ** 3)
      IF (iom_use("wind_wave_power")) CALL iom_put("wind_wave_power", 1000. * rau0 * tmask(:, :, 1) * zustar ** 2 * zustke)
      IF (iom_use("zhbl")) CALL iom_put("zhbl", tmask(:, :, 1) * zhbl)
      IF (iom_use("zhml")) CALL iom_put("zhml", tmask(:, :, 1) * zhml)
      IF (iom_use("zdh")) CALL iom_put("zdh", tmask(:, :, 1) * zdh)
      IF (iom_use("zhol")) CALL iom_put("zhol", tmask(:, :, 1) * zhol)
      IF (iom_use("zwthav")) CALL iom_put("zwthav", tmask(:, :, 1) * zwthav)
      IF (iom_use("zwth_ent")) CALL iom_put("zwth_ent", tmask(:, :, 1) * zwth_ent)
      IF (iom_use("zt_ml")) CALL iom_put("zt_ml", tmask(:, :, 1) * zt_ml)
    END IF
    CALL lbc_lnk('zdfosm', p_avt(:, :, :), 'W', 1.)
    CALL profile_psy_data8 % PostEnd
  END SUBROUTINE zdf_osm
  SUBROUTINE zdf_osm_init
    INTEGER :: ios
    INTEGER :: ji, jj, jk
    NAMELIST /namzdf_osm/ ln_use_osm_la, rn_osm_la, rn_osm_dstokes, nn_ave, nn_osm_wave, ln_dia_osm, rn_osm_hbl0, ln_kpprimix, rn_riinfty, rn_difri, ln_convmix, rn_difconv
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf_osm, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_osm in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf_osm, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namzdf_osm in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf_osm)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'zdf_osm_init : OSMOSIS Parameterisation'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namzdf_osm : set tke mixing parameters'
      WRITE(numout, *) '     Use namelist  rn_osm_la                     ln_use_osm_la = ', ln_use_osm_la
      WRITE(numout, *) '     Turbulent Langmuir number                     rn_osm_la   = ', rn_osm_la
      WRITE(numout, *) '     Initial hbl for 1D runs                       rn_osm_hbl0   = ', rn_osm_hbl0
      WRITE(numout, *) '     Depth scale of Stokes drift                rn_osm_dstokes = ', rn_osm_dstokes
      WRITE(numout, *) '     horizontal average flag                       nn_ave      = ', nn_ave
      WRITE(numout, *) '     Stokes drift                                  nn_osm_wave = ', nn_osm_wave
      SELECT CASE (nn_osm_wave)
      CASE (0)
        WRITE(numout, *) '     calculated assuming constant La#=0.3'
      CASE (1)
        WRITE(numout, *) '     calculated from Pierson Moskowitz wind-waves'
      CASE (2)
        WRITE(numout, *) '     calculated from ECMWF wave fields'
      END SELECT
      WRITE(numout, *) '     Output osm diagnostics                       ln_dia_osm  = ', ln_dia_osm
      WRITE(numout, *) '     Use KPP-style shear instability mixing       ln_kpprimix = ', ln_kpprimix
      WRITE(numout, *) '     local Richardson Number limit for shear instability rn_riinfty = ', rn_riinfty
      WRITE(numout, *) '     maximum shear diffusivity at Rig = 0    (m2/s) rn_difri = ', rn_difri
      WRITE(numout, *) '     Use large mixing below BL when unstable       ln_convmix = ', ln_convmix
      WRITE(numout, *) '     diffusivity when unstable below BL     (m2/s) rn_difconv = ', rn_difconv
    END IF
    IF (zdf_osm_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_osm_init : unable to allocate arrays')
    CALL osm_rst(nit000, 'READ')
    IF (ln_zdfddm) THEN
      IF (lwp) THEN
        WRITE(numout, *)
        WRITE(numout, *) '    Double diffusion mixing on temperature and salinity '
        WRITE(numout, *) '    CAUTION : done in routine zdfosm, not in routine zdfddm '
      END IF
    END IF
    IF (lwp) THEN
      WRITE(numout, *)
    END IF
    !$ACC KERNELS
    IF (nn_osm_wave == 0) THEN
      dstokes(:, :) = rn_osm_dstokes
    END IF
    !$ACC END KERNELS
    SELECT CASE (nn_ave)
    CASE (0)
      IF (lwp) WRITE(numout, *) '          no horizontal average on avt'
      IF (lwp) WRITE(numout, *) '          only in very high horizontal resolution !'
      !$ACC KERNELS
      etmean(:, :, :) = 0.E0
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = 2, jpjm1
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = 2, jpim1
            etmean(ji, jj, jk) = tmask(ji, jj, jk) / MAX(1., umask(ji - 1, jj, jk) + umask(ji, jj, jk) + vmask(ji, jj - 1, jk) + vmask(ji, jj, jk))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE (1)
      IF (lwp) WRITE(numout, *) '          horizontal average on avt'
      !$ACC KERNELS
      etmean(:, :, :) = 0.E0
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO jj = 2, jpjm1
          !$ACC LOOP INDEPENDENT ! CDe
          DO ji = 2, jpim1
            etmean(ji, jj, jk) = tmask(ji, jj, jk) / MAX(1., 2. * tmask(ji, jj, jk) + .5 * (tmask(ji - 1, jj + 1, jk) + tmask(ji - 1, jj - 1, jk) + tmask(ji + 1, jj + 1, jk) + tmask(ji + 1, jj - 1, jk)) + 1. * (tmask(ji - 1, jj, jk) + tmask(ji, jj + 1, jk) + tmask(ji, jj - 1, jk) + tmask(ji + 1, jj, jk)))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE DEFAULT
      WRITE(ctmp1, *) '          bad flag value for nn_ave = ', nn_ave
      CALL ctl_stop(ctmp1)
    END SELECT
    !$ACC KERNELS
    DO jk = 1, jpk
      avt(:, :, jk) = avtb(jk) * tmask(:, :, jk)
    END DO
    ghamt(:, :, :) = 0.
    ghams(:, :, :) = 0.
    ghamu(:, :, :) = 0.
    ghamv(:, :, :) = 0.
    !$ACC END KERNELS
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('wn')
      CALL iom_set_rstw_var_active('hbl')
      CALL iom_set_rstw_var_active('hbli')
    END IF
  END SUBROUTINE zdf_osm_init
  SUBROUTINE osm_rst(kt, cdrw)
    INTEGER, INTENT(IN) :: kt
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER :: id1, id2
    INTEGER :: ji, jj, jk
    INTEGER :: iiki, ikt
    REAL(KIND = wp) :: zhbf
    REAL(KIND = wp) :: zN2_c
    REAL(KIND = wp) :: rho_c = 0.01_wp
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: imld_rst
    IF (TRIM(cdrw) == 'READ' .AND. ln_rstart) THEN
      id1 = iom_varid(numror, 'wn', ldstop = .FALSE.)
      IF (id1 > 0) THEN
        CALL iom_get(numror, jpdom_autoglo, 'wn', wn, ldxios = lrxios)
        WRITE(numout, *) ' ===>>>> :  wn read from restart file'
      ELSE
        !$ACC KERNELS
        wn(:, :, :) = 0._wp
        !$ACC END KERNELS
        WRITE(numout, *) ' ===>>>> :  wn not in restart file, set to zero initially'
      END IF
      id1 = iom_varid(numror, 'hbl', ldstop = .FALSE.)
      id2 = iom_varid(numror, 'hbli', ldstop = .FALSE.)
      IF (id1 > 0 .AND. id2 > 0) THEN
        CALL iom_get(numror, jpdom_autoglo, 'hbl', hbl, ldxios = lrxios)
        CALL iom_get(numror, jpdom_autoglo, 'hbli', hbli, ldxios = lrxios)
        WRITE(numout, *) ' ===>>>> :  hbl & hbli read from restart file'
        RETURN
      ELSE
        WRITE(numout, *) ' ===>>>> : previous run without osmosis scheme, hbl computed from stratification'
      END IF
    END IF
    IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, *) '---- osm-rst ----'
      CALL iom_rstput(kt, nitrst, numrow, 'wn', wn, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'hbl', hbl, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'hbli', hbli, ldxios = lwxios)
      RETURN
    END IF
    IF (lwp) WRITE(numout, *) ' ===>>>> : calculating hbl computed from stratification'
    ALLOCATE(imld_rst(jpi, jpj))
    CALL eos_rab(tsn, rab_n)
    CALL bn2(tsn, rab_n, rn2)
    !$ACC KERNELS
    imld_rst(:, :) = nlb10
    hbl(:, :) = 0._wp
    zN2_c = grav * rho_c * r1_rau0
    hbl(:, :) = 0._wp
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 1, jpj
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 1, jpi
          ikt = mbkt(ji, jj)
          hbl(ji, jj) = hbl(ji, jj) + MAX(rn2(ji, jj, jk), 0._wp) * e3w_n(ji, jj, jk)
          IF (hbl(ji, jj) < zN2_c) imld_rst(ji, jj) = MIN(jk, ikt) + 1
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT ! CDe
    DO jj = 1, jpj
      !$ACC LOOP INDEPENDENT ! CDe
      DO ji = 1, jpi
        iiki = imld_rst(ji, jj)
        hbl(ji, jj) = gdepw_n(ji, jj, iiki) * ssmask(ji, jj)
      END DO
    END DO
    hbl = MAX(hbl, epsln)
    hbli(:, :) = hbl(:, :)
    !$ACC END KERNELS
    DEALLOCATE(imld_rst)
    WRITE(numout, *) ' ===>>>> : hbl computed from stratification'
  END SUBROUTINE osm_rst
  SUBROUTINE tra_osm(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdt, ztrds
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    CALL profile_psy_data0 % PreStart('tra_osm', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'tra_osm : OSM non-local tracer fluxes'
      IF (lwp) WRITE(numout, *) '~~~~~~~   '
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      !$ACC END KERNELS
      ALLOCATE(ztrds(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrds(:, :, :) = tsa(:, :, :, jp_sal)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) - (ghamt(ji, jj, jk) - ghamt(ji, jj, jk + 1)) / e3t_n(ji, jj, jk)
          tsa(ji, jj, jk, jp_sal) = tsa(ji, jj, jk, jp_sal) - (ghams(ji, jj, jk) - ghams(ji, jj, jk + 1)) / e3t_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (l_trdtra) THEN
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem) - ztrdt(:, :, :)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal) - ztrds(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('tra_osm', 'r1', 0, 0)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_zdf, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_zdf, ztrds)
      DEALLOCATE(ztrdt)
      DEALLOCATE(ztrds)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('tra_osm', 'r2', 0, 0)
    IF (ln_ctl) THEN
      CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' osm  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
    END IF
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE tra_osm
  SUBROUTINE trc_osm(kt)
    INTEGER, INTENT(IN) :: kt
    WRITE(*, *) 'trc_osm: Not written yet', kt
  END SUBROUTINE trc_osm
  SUBROUTINE dyn_osm(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dyn_osm', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, *)
      IF (lwp) WRITE(numout, *) 'dyn_osm : OSM non-local velocity'
      IF (lwp) WRITE(numout, *) '~~~~~~~   '
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT ! CDe
      DO jj = 2, jpjm1
        !$ACC LOOP INDEPENDENT ! CDe
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (ghamu(ji, jj, jk) - ghamu(ji, jj, jk + 1)) / e3u_n(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) - (ghamv(ji, jj, jk) - ghamv(ji, jj, jk + 1)) / e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE dyn_osm
END MODULE zdfosm
