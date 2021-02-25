MODULE zdfmxl
  USE oce
  USE dom_oce
  USE trc_oce, ONLY: l_offline
  USE zdf_oce
  USE eosbn2
  USE in_out_manager
  USE prtctl
  USE phycst
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_mxl
  INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: nmln
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hmld
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hmlp
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: hmlpt
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: hmld_zint
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: htc_mld
  LOGICAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: ll_found
  LOGICAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: ll_belowml
  REAL(KIND = wp), PUBLIC :: rho_c = 0.01_wp
  REAL(KIND = wp), PUBLIC :: avt_c = 5.E-4_wp
  TYPE, PUBLIC :: MXL_ZINT
    INTEGER :: mld_type
    REAL(KIND = wp) :: zref
    REAL(KIND = wp) :: dT_crit
    REAL(KIND = wp) :: iso_frac
  END TYPE MXL_ZINT
  CONTAINS
  INTEGER FUNCTION zdf_mxl_alloc()
    zdf_mxl_alloc = 0
    IF (.NOT. ALLOCATED(nmln)) THEN
      ALLOCATE(nmln(jpi, jpj), hmld(jpi, jpj), hmlp(jpi, jpj), hmlpt(jpi, jpj), hmld_zint(jpi, jpj), htc_mld(jpi, jpj), &
&ll_found(jpi, jpj), ll_belowml(jpi, jpj, jpk), STAT = zdf_mxl_alloc)
      CALL mpp_sum('zdfmxl', zdf_mxl_alloc)
      IF (zdf_mxl_alloc /= 0) CALL ctl_stop('STOP', 'zdf_mxl_alloc: failed to allocate arrays.')
    END IF
  END FUNCTION zdf_mxl_alloc
  SUBROUTINE zdf_mxl(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: iikn, iiki, ikt
    REAL(KIND = wp) :: zN2_c
    INTEGER, DIMENSION(jpi, jpj) :: imld
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'zdf_mxl : mixed layer depth'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
      IF (zdf_mxl_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_mxl : unable to allocate arrays')
    END IF
    nmln(:, :) = nlb10
    hmlp(:, :) = 0._wp
    zN2_c = grav * rho_c * r1_rau0
    DO jk = nlb10, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          ikt = mbkt(ji, jj)
          hmlp(ji, jj) = hmlp(ji, jj) + MAX(rn2b(ji, jj, jk), 0._wp) * e3w_n(ji, jj, jk)
          IF (hmlp(ji, jj) < zN2_c) nmln(ji, jj) = MIN(jk, ikt) + 1
        END DO
      END DO
    END DO
    imld(:, :) = mbkt(:, :) + 1
    DO jk = jpkm1, nlb10, - 1
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (avt(ji, jj, jk) < avt_c * wmask(ji, jj, jk)) imld(ji, jj) = jk
        END DO
      END DO
    END DO
    DO jj = 1, jpj
      DO ji = 1, jpi
        iiki = imld(ji, jj)
        iikn = nmln(ji, jj)
        hmld(ji, jj) = gdepw_n(ji, jj, iiki) * ssmask(ji, jj)
        hmlp(ji, jj) = gdepw_n(ji, jj, iikn) * ssmask(ji, jj)
        hmlpt(ji, jj) = gdept_n(ji, jj, iikn - 1) * ssmask(ji, jj)
      END DO
    END DO
    IF (.NOT. l_offline) THEN
      IF (iom_use("mldr10_1")) THEN
        IF (ln_isfcav) THEN
          CALL iom_put("mldr10_1", hmlp - risfdep)
        ELSE
          CALL iom_put("mldr10_1", hmlp)
        END IF
      END IF
      IF (iom_use("mldkz5")) THEN
        IF (ln_isfcav) THEN
          CALL iom_put("mldkz5", hmld - risfdep)
        ELSE
          CALL iom_put("mldkz5", hmld)
        END IF
      END IF
    END IF
    CALL zdf_mxl_zint(kt)
    IF (ln_ctl) CALL prt_ctl(tab2d_1 = REAL(nmln, wp), clinfo1 = ' nmln : ', tab2d_2 = hmlp, clinfo2 = ' hmlp : ')
  END SUBROUTINE zdf_mxl
  SUBROUTINE zdf_mxl_zint_mld(sf)
    TYPE(MXL_ZINT), INTENT(IN) :: sf
    INTEGER :: nn_mld_type
    REAL(KIND = wp) :: rn_zref
    REAL(KIND = wp) :: rn_dT_crit
    REAL(KIND = wp) :: rn_iso_frac
    REAL(KIND = wp), PARAMETER :: zepsilon = 1.E-30
    INTEGER, DIMENSION(jpi, jpj) :: ikmt
    INTEGER, DIMENSION(jpi, jpj) :: ik_ref
    INTEGER, DIMENSION(jpi, jpj) :: ik_iso
    REAL, DIMENSION(jpi, jpj, jpk) :: zT
    REAL, DIMENSION(jpi, jpj) :: ppzdep
    REAL, DIMENSION(jpi, jpj) :: zT_ref
    REAL :: zT_b
    REAL, DIMENSION(jpi, jpj, jpk) :: zdTdz
    REAL, DIMENSION(jpi, jpj, jpk) :: zmoddT
    REAL :: zdz
    REAL :: zdT
    REAL, DIMENSION(jpi, jpj) :: zdelta_T
    REAL, DIMENSION(jpi, jpj) :: zRHO1, zRHO2
    INTEGER :: ji, jj, jk
    nn_mld_type = sf % mld_type
    rn_zref = sf % zref
    rn_dT_crit = sf % dT_crit
    rn_iso_frac = sf % iso_frac
    IF (nn_mld_type == 0) THEN
      zdelta_T(:, :) = rn_dT_crit
      zT(:, :, :) = rhop(:, :, :)
    ELSE IF (nn_mld_type == 1) THEN
      ppzdep(:, :) = 0.0
      CALL eos(tsn(:, :, 1, :), ppzdep(:, :), zRHO1(:, :))
      zT(:, :, 1 : jpts) = tsn(:, :, 1, 1 : jpts)
      zT(:, :, jp_tem) = zT(:, :, 1) + rn_dT_crit
      CALL eos(zT(:, :, 1 : jpts), ppzdep(:, :), zRHO2(:, :))
      zdelta_T(:, :) = ABS(zRHO1(:, :) - zRHO2(:, :)) * rau0
      CALL lbc_lnk('zdfmxl', zdelta_T, 'T', 1.)
      zT(:, :, :) = rhop(:, :, :)
    ELSE
      zdelta_T(:, :) = rn_dT_crit
      zT(:, :, :) = tsn(:, :, :, jp_tem)
    END IF
    DO jk = 1, jpk - 2
      zdTdz(:, :, jk) = (zT(:, :, jk + 1) - zT(:, :, jk)) / e3w_n(:, :, jk + 1)
      zmoddT(:, :, jk) = ABS(zT(:, :, jk + 1) - zT(:, :, jk))
    END DO
    DO jk = jpkm1, 2, - 1
      WHERE (gdept_n(:, :, jk) > rn_zref)
        ik_ref(:, :) = jk - 1
        zT_ref(:, :) = zT(:, :, jk - 1) + zdTdz(:, :, jk - 1) * (rn_zref - gdept_n(:, :, jk - 1))
      END WHERE
    END DO
    WHERE (gdept_n(:, :, 1) > rn_zref)
      zT_ref = zT(:, :, 1)
      ik_ref = 1
    END WHERE
    ikmt(:, :) = mbkt(:, :) - 1
    ll_found(:, :) = .FALSE.
    IF (rn_iso_frac - zepsilon > 0.) THEN
      ik_iso(:, :) = ik_ref(:, :)
      DO jj = 1, nlcj
        DO ji = 1, nlci
          DO jk = ik_ref(ji, jj), ikmt(ji, jj) - 1
            IF (zmoddT(ji, jj, jk) > (rn_iso_frac * zdelta_T(ji, jj))) THEN
              ik_iso(ji, jj) = jk
              ll_found(ji, jj) = (zmoddT(ji, jj, jk) > zdelta_T(ji, jj))
              EXIT
            END IF
          END DO
        END DO
      END DO
      hmld_zint(:, :) = rn_zref
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (ll_found(ji, jj) .AND. tmask(ji, jj, 1) == 1.0) THEN
            zdz = ABS(zdelta_T(ji, jj) / zdTdz(ji, jj, ik_iso(ji, jj)))
            hmld_zint(ji, jj) = gdept_n(ji, jj, ik_iso(ji, jj)) + zdz
          END IF
        END DO
      END DO
    END IF
    WHERE (tmask(:, :, 1) /= 1.0)
      ll_found = .TRUE.
    END WHERE
    DO jk = 1, jpk
      ll_belowml(:, :, jk) = ABS(zT(:, :, jk) - zT_ref(:, :)) >= zdelta_T(:, :)
    END DO
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (.NOT. ll_found(ji, jj)) hmld_zint(ji, jj) = gdept_n(ji, jj, ikmt(ji, jj))
      END DO
    END DO
    DO jj = 1, jpj
      DO ji = 1, jpi
        DO jk = ik_ref(ji, jj) + 1, ikmt(ji, jj)
          IF (ll_found(ji, jj)) EXIT
          IF (ll_belowml(ji, jj, jk)) THEN
            zT_b = zT_ref(ji, jj) + zdelta_T(ji, jj) * SIGN(1.0, zdTdz(ji, jj, jk - 1))
            zdT = zT_b - zT(ji, jj, jk - 1)
            zdz = zdT / zdTdz(ji, jj, jk - 1)
            hmld_zint(ji, jj) = gdept_n(ji, jj, jk - 1) + zdz
            EXIT
          END IF
        END DO
      END DO
    END DO
    hmld_zint(:, :) = hmld_zint(:, :) * tmask(:, :, 1)
  END SUBROUTINE zdf_mxl_zint_mld
  SUBROUTINE zdf_mxl_zint_htc(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikmax
    REAL(KIND = wp) :: zc, zcoef
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ilevel
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zthick_0, zthick
    IF (.NOT. ALLOCATED(ilevel)) THEN
      ALLOCATE(ilevel(jpi, jpj), zthick_0(jpi, jpj), zthick(jpi, jpj), STAT = ji)
      IF (lk_mpp) CALL mpp_sum('zdfmxl', ji)
      IF (ji /= 0) CALL ctl_stop('STOP', 'zdf_mxl_zint_htc : unable to allocate arrays')
    END IF
    ilevel(:, :) = 0
    zthick_0(:, :) = 0._wp
    DO jk = 1, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          zthick_0(ji, jj) = zthick_0(ji, jj) + e3t_n(ji, jj, jk)
          IF (zthick_0(ji, jj) < hmld_zint(ji, jj)) ilevel(ji, jj) = jk
        END DO
      END DO
      WRITE(numout, FMT = *) 'zthick_0(jk =', jk, ') =', zthick_0(2, 2)
      WRITE(numout, FMT = *) 'gdepw_n(jk+1 =', jk + 1, ') =', gdepw_n(2, 2, jk + 1)
    END DO
    IF (ln_linssh) THEN
      zthick(:, :) = sshn(:, :)
      htc_mld(:, :) = tsn(:, :, 1, jp_tem) * sshn(:, :) * tmask(:, :, 1)
    ELSE
      zthick(:, :) = 0._wp
      htc_mld(:, :) = 0._wp
    END IF
    ikmax = MIN(MAXVAL(ilevel(:, :)), jpkm1)
    DO jk = 1, ikmax
      DO jj = 1, jpj
        DO ji = 1, jpi
          zc = e3t_n(ji, jj, jk) * REAL(MIN(MAX(0, ilevel(ji, jj) - jk + 1), 1))
          zthick(ji, jj) = zthick(ji, jj) + zc
          htc_mld(ji, jj) = htc_mld(ji, jj) + zc * tsn(ji, jj, jk, jp_tem) * tmask(ji, jj, jk)
        END DO
      END DO
    END DO
    zthick(:, :) = hmld_zint(:, :) - zthick(:, :)
    DO jj = 1, jpj
      DO ji = 1, jpi
        htc_mld(ji, jj) = htc_mld(ji, jj) + tsn(ji, jj, ilevel(ji, jj) + 1, jp_tem) * MIN(e3t_n(ji, jj, ilevel(ji, jj) + 1), &
&zthick(ji, jj)) * tmask(ji, jj, ilevel(ji, jj) + 1)
      END DO
    END DO
    WRITE(numout, FMT = *) 'htc_mld(after) =', htc_mld(2, 2)
    zcoef = rau0 * rcp
    htc_mld(:, :) = zcoef * htc_mld(:, :)
  END SUBROUTINE zdf_mxl_zint_htc
  SUBROUTINE zdf_mxl_zint(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ios
    INTEGER :: jn
    INTEGER :: nn_mld_diag = 0
    CHARACTER(LEN = 1) :: cmld
    TYPE(MXL_ZINT) :: sn_mld1, sn_mld2, sn_mld3, sn_mld4, sn_mld5
    TYPE(MXL_ZINT), SAVE, DIMENSION(5) :: mld_diags
    NAMELIST /namzdf_mldzint/ nn_mld_diag, sn_mld1, sn_mld2, sn_mld3, sn_mld4, sn_mld5
    IF (kt == nit000) THEN
      REWIND(UNIT = numnam_ref)
      READ(numnam_ref, namzdf_mldzint, IOSTAT = ios, ERR = 901)
901   IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_mldzint in reference namelist', lwp)
      REWIND(UNIT = numnam_cfg)
      READ(numnam_cfg, namzdf_mldzint, IOSTAT = ios, ERR = 902)
902   IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_mldzint in configuration namelist', lwp)
      IF (lwm) WRITE(numond, namzdf_mldzint)
      IF (nn_mld_diag > 5) CALL ctl_stop('STOP', 'zdf_mxl_ini: Specify no more than 5 MLD definitions')
      mld_diags(1) = sn_mld1
      mld_diags(2) = sn_mld2
      mld_diags(3) = sn_mld3
      mld_diags(4) = sn_mld4
      mld_diags(5) = sn_mld5
      IF (lwp .AND. (nn_mld_diag > 0)) THEN
        WRITE(numout, FMT = *) '=============== Vertically-interpolated mixed layer ================'
        WRITE(numout, FMT = *) '(Diagnostic number, nn_mld_type, rn_zref, rn_dT_crit, rn_iso_frac)'
        DO jn = 1, nn_mld_diag
          WRITE(numout, FMT = *) 'MLD criterion', jn, ':'
          WRITE(numout, FMT = *) '    nn_mld_type =', mld_diags(jn) % mld_type
          WRITE(numout, FMT = *) '    rn_zref =', mld_diags(jn) % zref
          WRITE(numout, FMT = *) '    rn_dT_crit =', mld_diags(jn) % dT_crit
          WRITE(numout, FMT = *) '    rn_iso_frac =', mld_diags(jn) % iso_frac
        END DO
        WRITE(numout, FMT = *) '===================================================================='
      END IF
    END IF
    IF (nn_mld_diag > 0) THEN
      DO jn = 1, nn_mld_diag
        WRITE(cmld, FMT = '(I1)') jn
        IF (iom_use("mldzint_" // cmld) .OR. iom_use("mldhtc_" // cmld)) THEN
          CALL zdf_mxl_zint_mld(mld_diags(jn))
          IF (iom_use("mldzint_" // cmld)) THEN
            CALL iom_put("mldzint_" // cmld, hmld_zint(:, :))
          END IF
          IF (iom_use("mldhtc_" // cmld)) THEN
            CALL zdf_mxl_zint_htc(kt)
            CALL iom_put("mldhtc_" // cmld, htc_mld(:, :))
          END IF
        END IF
      END DO
    END IF
  END SUBROUTINE zdf_mxl_zint
END MODULE zdfmxl