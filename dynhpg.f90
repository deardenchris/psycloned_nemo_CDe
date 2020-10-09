MODULE dynhpg
  USE oce
  USE sbc_oce
  USE dom_oce
  USE wet_dry
  USE phycst
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE prtctl
  USE lbclnk
  USE lib_mpp
  USE eosbn2
  USE timing
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_hpg
  PUBLIC :: dyn_hpg_init
  LOGICAL, PUBLIC :: ln_hpg_zco
  LOGICAL, PUBLIC :: ln_hpg_zps
  LOGICAL, PUBLIC :: ln_hpg_sco
  LOGICAL, PUBLIC :: ln_hpg_djc
  LOGICAL, PUBLIC :: ln_hpg_prj
  LOGICAL, PUBLIC :: ln_hpg_isf
  INTEGER, PARAMETER :: np_ERROR = - 10
  INTEGER, PARAMETER :: np_zco = 0
  INTEGER, PARAMETER :: np_zps = 1
  INTEGER, PARAMETER :: np_sco = 2
  INTEGER, PARAMETER :: np_djc = 3
  INTEGER, PARAMETER :: np_prj = 4
  INTEGER, PARAMETER :: np_isf = 5
  INTEGER, PUBLIC :: nhpg
  CONTAINS
  SUBROUTINE dyn_hpg(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdu, ztrdv
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_timing) CALL timing_start('dyn_hpg')
    IF (l_trddyn) THEN
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    CALL profile_psy_data0 % PreStart('dyn_hpg', 'r0', 0, 0)
    SELECT CASE (nhpg)
    CASE (np_zco)
      CALL hpg_zco(kt)
    CASE (np_zps)
      CALL hpg_zps(kt)
    CASE (np_sco)
      CALL hpg_sco(kt)
    CASE (np_djc)
      CALL hpg_djc(kt)
    CASE (np_prj)
      CALL hpg_prj(kt)
    CASE (np_isf)
      CALL hpg_isf(kt)
    END SELECT
    CALL profile_psy_data0 % PostEnd
    IF (l_trddyn) THEN
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('dyn_hpg', 'r1', 0, 0)
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_hpg, kt)
      DEALLOCATE(ztrdu, ztrdv)
      CALL profile_psy_data1 % PostEnd
    END IF
    CALL profile_psy_data2 % PreStart('dyn_hpg', 'r2', 0, 0)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' hpg  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = &
&vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_hpg')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE dyn_hpg
  SUBROUTINE dyn_hpg_init
    INTEGER :: ioptio = 0
    INTEGER :: ios
    INTEGER :: ji, jj, jk, ikt
    REAL(KIND = wp) :: znad
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zts_top, zrhd
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zrhdtop_isf
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: ziceload
    NAMELIST /namdyn_hpg/ ln_hpg_zco, ln_hpg_zps, ln_hpg_sco, ln_hpg_djc, ln_hpg_prj, ln_hpg_isf
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namdyn_hpg, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn_hpg in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namdyn_hpg, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn_hpg in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namdyn_hpg)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dyn_hpg_init : hydrostatic pressure gradient initialisation'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namdyn_hpg : choice of hpg scheme'
      WRITE(numout, FMT = *) '      z-coord. - full steps                             ln_hpg_zco    = ', ln_hpg_zco
      WRITE(numout, FMT = *) '      z-coord. - partial steps (interpolation)          ln_hpg_zps    = ', ln_hpg_zps
      WRITE(numout, FMT = *) '      s-coord. (standard jacobian formulation)          ln_hpg_sco    = ', ln_hpg_sco
      WRITE(numout, FMT = *) '      s-coord. (standard jacobian formulation) for isf  ln_hpg_isf    = ', ln_hpg_isf
      WRITE(numout, FMT = *) '      s-coord. (Density Jacobian: Cubic polynomial)     ln_hpg_djc    = ', ln_hpg_djc
      WRITE(numout, FMT = *) '      s-coord. (Pressure Jacobian: Cubic polynomial)    ln_hpg_prj    = ', ln_hpg_prj
    END IF
    IF (ln_hpg_djc) CALL ctl_stop('dyn_hpg_init : Density Jacobian: Cubic polynominal method', '   currently disabled (bugs under &
&investigation).', '   Please select either  ln_hpg_sco or  ln_hpg_prj instead')
    IF (.NOT. ln_linssh .AND. .NOT. (ln_hpg_sco .OR. ln_hpg_prj .OR. ln_hpg_isf)) CALL ctl_stop('dyn_hpg_init : non-linear free &
&surface requires either ', '   the standard jacobian formulation hpg_sco    or ', '   the pressure jacobian formulation hpg_prj')
    IF (ln_hpg_isf) THEN
      IF (.NOT. ln_isfcav) CALL ctl_stop(' hpg_isf not available if ln_isfcav = false ')
    ELSE
      IF (ln_isfcav) CALL ctl_stop('Only hpg_isf has been corrected to work with ice shelf cavity.')
    END IF
    nhpg = np_ERROR
    ioptio = 0
    IF (ln_hpg_zco) THEN
      nhpg = np_zco
      ioptio = ioptio + 1
    END IF
    IF (ln_hpg_zps) THEN
      nhpg = np_zps
      ioptio = ioptio + 1
    END IF
    IF (ln_hpg_sco) THEN
      nhpg = np_sco
      ioptio = ioptio + 1
    END IF
    IF (ln_hpg_djc) THEN
      nhpg = np_djc
      ioptio = ioptio + 1
    END IF
    IF (ln_hpg_prj) THEN
      nhpg = np_prj
      ioptio = ioptio + 1
    END IF
    IF (ln_hpg_isf) THEN
      nhpg = np_isf
      ioptio = ioptio + 1
    END IF
    IF (ioptio /= 1) CALL ctl_stop('NO or several hydrostatic pressure gradient options used')
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      SELECT CASE (nhpg)
      CASE (np_zco)
        WRITE(numout, FMT = *) '   ==>>>   z-coord. - full steps '
      CASE (np_zps)
        WRITE(numout, FMT = *) '   ==>>>   z-coord. - partial steps (interpolation)'
      CASE (np_sco)
        WRITE(numout, FMT = *) '   ==>>>   s-coord. (standard jacobian formulation)'
      CASE (np_djc)
        WRITE(numout, FMT = *) '   ==>>>   s-coord. (Density Jacobian: Cubic polynomial)'
      CASE (np_prj)
        WRITE(numout, FMT = *) '   ==>>>   s-coord. (Pressure Jacobian: Cubic polynomial)'
      CASE (np_isf)
        WRITE(numout, FMT = *) '   ==>>>   s-coord. (standard jacobian formulation) for isf'
      END SELECT
      WRITE(numout, FMT = *)
    END IF
    IF (.NOT. ln_isfcav) THEN
      !$ACC KERNELS
      riceload(:, :) = 0._wp
      !$ACC END KERNELS
    ELSE
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) '   ice shelf case: set the ice-shelf load'
      ALLOCATE(zts_top(jpi, jpj, jpts), zrhd(jpi, jpj, jpk), zrhdtop_isf(jpi, jpj), ziceload(jpi, jpj))
      !$ACC KERNELS
      znad = 1._wp
      zts_top(:, :, jp_tem) = - 1.9_wp
      zts_top(:, :, jp_sal) = 34.4_wp
      !$ACC END KERNELS
      DO jk = 1, jpk
        CALL eos(zts_top(:, :, :), gdept_n(:, :, jk), zrhd(:, :, jk))
      END DO
      CALL eos(zts_top, risfdep, zrhdtop_isf)
      !$ACC KERNELS
      ziceload = 0._wp
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          ikt = mikt(ji, jj)
          ziceload(ji, jj) = ziceload(ji, jj) + (znad + zrhd(ji, jj, 1)) * e3w_n(ji, jj, 1) * (1._wp - tmask(ji, jj, 1))
          DO jk = 2, ikt - 1
            ziceload(ji, jj) = ziceload(ji, jj) + (2._wp * znad + zrhd(ji, jj, jk - 1) + zrhd(ji, jj, jk)) * e3w_n(ji, jj, jk) * &
&(1._wp - tmask(ji, jj, jk))
          END DO
          IF (ikt >= 2) ziceload(ji, jj) = ziceload(ji, jj) + (2._wp * znad + zrhdtop_isf(ji, jj) + zrhd(ji, jj, ikt - 1)) * &
&(risfdep(ji, jj) - gdept_n(ji, jj, ikt - 1))
        END DO
      END DO
      riceload(:, :) = ziceload(:, :)
      !$ACC END KERNELS
      DEALLOCATE(zts_top, zrhd, zrhdtop_isf, ziceload)
    END IF
  END SUBROUTINE dyn_hpg_init
  SUBROUTINE hpg_zco(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoef0, zcoef1
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhpi, zhpj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('hpg_zco', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn:hpg_zco : hydrostatic pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   z-coordinate case '
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zcoef0 = - grav * 0.5_wp
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zcoef1 = zcoef0 * e3w_n(ji, jj, 1)
        zhpi(ji, jj, 1) = zcoef1 * (rhd(ji + 1, jj, 1) - rhd(ji, jj, 1)) * r1_e1u(ji, jj)
        zhpj(ji, jj, 1) = zcoef1 * (rhd(ji, jj + 1, 1) - rhd(ji, jj, 1)) * r1_e2v(ji, jj)
        ua(ji, jj, 1) = ua(ji, jj, 1) + zhpi(ji, jj, 1)
        va(ji, jj, 1) = va(ji, jj, 1) + zhpj(ji, jj, 1)
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zcoef1 = zcoef0 * e3w_n(ji, jj, jk)
          zhpi(ji, jj, jk) = zhpi(ji, jj, jk - 1) + zcoef1 * ((rhd(ji + 1, jj, jk) + rhd(ji + 1, jj, jk - 1)) - (rhd(ji, jj, jk) + &
&rhd(ji, jj, jk - 1))) * r1_e1u(ji, jj)
          zhpj(ji, jj, jk) = zhpj(ji, jj, jk - 1) + zcoef1 * ((rhd(ji, jj + 1, jk) + rhd(ji, jj + 1, jk - 1)) - (rhd(ji, jj, jk) + &
&rhd(ji, jj, jk - 1))) * r1_e2v(ji, jj)
          ua(ji, jj, jk) = ua(ji, jj, jk) + zhpi(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) + zhpj(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE hpg_zco
  SUBROUTINE hpg_zps(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: iku, ikv
    REAL(KIND = wp) :: zcoef0, zcoef1, zcoef2, zcoef3
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhpi, zhpj
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('hpg_zps', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn:hpg_zps : hydrostatic pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   z-coordinate with partial steps - vector optimization'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zcoef0 = - grav * 0.5_wp
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zcoef1 = zcoef0 * e3w_n(ji, jj, 1)
        zhpi(ji, jj, 1) = zcoef1 * (rhd(ji + 1, jj, 1) - rhd(ji, jj, 1)) * r1_e1u(ji, jj)
        zhpj(ji, jj, 1) = zcoef1 * (rhd(ji, jj + 1, 1) - rhd(ji, jj, 1)) * r1_e2v(ji, jj)
        ua(ji, jj, 1) = ua(ji, jj, 1) + zhpi(ji, jj, 1)
        va(ji, jj, 1) = va(ji, jj, 1) + zhpj(ji, jj, 1)
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zcoef1 = zcoef0 * e3w_n(ji, jj, jk)
          zhpi(ji, jj, jk) = zhpi(ji, jj, jk - 1) + zcoef1 * ((rhd(ji + 1, jj, jk) + rhd(ji + 1, jj, jk - 1)) - (rhd(ji, jj, jk) + &
&rhd(ji, jj, jk - 1))) * r1_e1u(ji, jj)
          zhpj(ji, jj, jk) = zhpj(ji, jj, jk - 1) + zcoef1 * ((rhd(ji, jj + 1, jk) + rhd(ji, jj + 1, jk - 1)) - (rhd(ji, jj, jk) + &
&rhd(ji, jj, jk - 1))) * r1_e2v(ji, jj)
          ua(ji, jj, jk) = ua(ji, jj, jk) + zhpi(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) + zhpj(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        iku = mbku(ji, jj)
        ikv = mbkv(ji, jj)
        zcoef2 = zcoef0 * MIN(e3w_n(ji, jj, iku), e3w_n(ji + 1, jj, iku))
        zcoef3 = zcoef0 * MIN(e3w_n(ji, jj, ikv), e3w_n(ji, jj + 1, ikv))
        IF (iku > 1) THEN
          ua(ji, jj, iku) = ua(ji, jj, iku) - zhpi(ji, jj, iku)
          zhpi(ji, jj, iku) = zhpi(ji, jj, iku - 1) + zcoef2 * (rhd(ji + 1, jj, iku - 1) - rhd(ji, jj, iku - 1) + gru(ji, jj)) * &
&r1_e1u(ji, jj)
          ua(ji, jj, iku) = ua(ji, jj, iku) + zhpi(ji, jj, iku)
        END IF
        IF (ikv > 1) THEN
          va(ji, jj, ikv) = va(ji, jj, ikv) - zhpj(ji, jj, ikv)
          zhpj(ji, jj, ikv) = zhpj(ji, jj, ikv - 1) + zcoef3 * (rhd(ji, jj + 1, ikv - 1) - rhd(ji, jj, ikv - 1) + grv(ji, jj)) * &
&r1_e2v(ji, jj)
          va(ji, jj, ikv) = va(ji, jj, ikv) + zhpj(ji, jj, ikv)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE hpg_zps
  SUBROUTINE hpg_sco(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jii, jjj
    REAL(KIND = wp) :: zcoef0, zuap, zvap, znad, ztmp
    LOGICAL :: ll_tmp1, ll_tmp2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhpi, zhpj
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zcpx, zcpy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('hpg_sco', 'r0', 0, 0)
    IF (ln_wd_il) ALLOCATE(zcpx(jpi, jpj), zcpy(jpi, jpj))
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn:hpg_sco : hydrostatic pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   s-coordinate case, OPA original scheme used'
    END IF
    zcoef0 = - grav * 0.5_wp
    IF (ln_linssh) THEN
      znad = 0._wp
    ELSE
      znad = 1._wp
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_wd_il) THEN
      !$ACC KERNELS
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ll_tmp1 = MIN(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) .AND. MAX(sshn(ji, jj) + &
&ht_0(ji, jj), sshn(ji + 1, jj) + ht_0(ji + 1, jj)) > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji + 1, jj)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, &
&jj), - ht_0(ji + 1, jj)) + rn_wdmin1 + rn_wdmin2)
          IF (ll_tmp1) THEN
            zcpx(ji, jj) = 1.0_wp
          ELSE IF (ll_tmp2) THEN
            zcpx(ji, jj) = ABS((sshn(ji + 1, jj) + ht_0(ji + 1, jj) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji + 1, jj) - sshn(ji, &
&jj)))
          ELSE
            zcpx(ji, jj) = 0._wp
          END IF
          ll_tmp1 = MIN(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) .AND. MAX(sshn(ji, jj) + &
&ht_0(ji, jj), sshn(ji, jj + 1) + ht_0(ji, jj + 1)) > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji, jj + 1)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, &
&jj), - ht_0(ji, jj + 1)) + rn_wdmin1 + rn_wdmin2)
          IF (ll_tmp1) THEN
            zcpy(ji, jj) = 1.0_wp
          ELSE IF (ll_tmp2) THEN
            zcpy(ji, jj) = ABS((sshn(ji, jj + 1) + ht_0(ji, jj + 1) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji, jj + 1) - sshn(ji, &
&jj)))
          ELSE
            zcpy(ji, jj) = 0._wp
          END IF
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi('dynhpg', zcpx, 'U', 1., zcpy, 'V', 1.)
    END IF
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zhpi(ji, jj, 1) = zcoef0 * (e3w_n(ji + 1, jj, 1) * (znad + rhd(ji + 1, jj, 1)) - e3w_n(ji, jj, 1) * (znad + rhd(ji, jj, &
&1))) * r1_e1u(ji, jj)
        zhpj(ji, jj, 1) = zcoef0 * (e3w_n(ji, jj + 1, 1) * (znad + rhd(ji, jj + 1, 1)) - e3w_n(ji, jj, 1) * (znad + rhd(ji, jj, &
&1))) * r1_e2v(ji, jj)
        zuap = - zcoef0 * (rhd(ji + 1, jj, 1) + rhd(ji, jj, 1) + 2._wp * znad) * (gde3w_n(ji + 1, jj, 1) - gde3w_n(ji, jj, 1)) * &
&r1_e1u(ji, jj)
        zvap = - zcoef0 * (rhd(ji, jj + 1, 1) + rhd(ji, jj, 1) + 2._wp * znad) * (gde3w_n(ji, jj + 1, 1) - gde3w_n(ji, jj, 1)) * &
&r1_e2v(ji, jj)
        IF (ln_wd_il) THEN
          zhpi(ji, jj, 1) = zhpi(ji, jj, 1) * zcpx(ji, jj)
          zhpj(ji, jj, 1) = zhpj(ji, jj, 1) * zcpy(ji, jj)
          zuap = zuap * zcpx(ji, jj)
          zvap = zvap * zcpy(ji, jj)
        END IF
        ua(ji, jj, 1) = ua(ji, jj, 1) + zhpi(ji, jj, 1) + zuap
        va(ji, jj, 1) = va(ji, jj, 1) + zhpj(ji, jj, 1) + zvap
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zhpi(ji, jj, jk) = zhpi(ji, jj, jk - 1) + zcoef0 * r1_e1u(ji, jj) * (e3w_n(ji + 1, jj, jk) * (rhd(ji + 1, jj, jk) + &
&rhd(ji + 1, jj, jk - 1) + 2 * znad) - e3w_n(ji, jj, jk) * (rhd(ji, jj, jk) + rhd(ji, jj, jk - 1) + 2 * znad))
          zhpj(ji, jj, jk) = zhpj(ji, jj, jk - 1) + zcoef0 * r1_e2v(ji, jj) * (e3w_n(ji, jj + 1, jk) * (rhd(ji, jj + 1, jk) + &
&rhd(ji, jj + 1, jk - 1) + 2 * znad) - e3w_n(ji, jj, jk) * (rhd(ji, jj, jk) + rhd(ji, jj, jk - 1) + 2 * znad))
          zuap = - zcoef0 * (rhd(ji + 1, jj, jk) + rhd(ji, jj, jk) + 2._wp * znad) * (gde3w_n(ji + 1, jj, jk) - gde3w_n(ji, jj, &
&jk)) * r1_e1u(ji, jj)
          zvap = - zcoef0 * (rhd(ji, jj + 1, jk) + rhd(ji, jj, jk) + 2._wp * znad) * (gde3w_n(ji, jj + 1, jk) - gde3w_n(ji, jj, &
&jk)) * r1_e2v(ji, jj)
          IF (ln_wd_il) THEN
            zhpi(ji, jj, jk) = zhpi(ji, jj, jk) * zcpx(ji, jj)
            zhpj(ji, jj, jk) = zhpj(ji, jj, jk) * zcpy(ji, jj)
            zuap = zuap * zcpx(ji, jj)
            zvap = zvap * zcpy(ji, jj)
          END IF
          ua(ji, jj, jk) = ua(ji, jj, jk) + zhpi(ji, jj, jk) + zuap
          va(ji, jj, jk) = va(ji, jj, jk) + zhpj(ji, jj, jk) + zvap
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_wd_il) DEALLOCATE(zcpx, zcpy)
  END SUBROUTINE hpg_sco
  SUBROUTINE hpg_isf(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, ikt, iktp1i, iktp1j
    REAL(KIND = wp) :: zcoef0, zuap, zvap, znad
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhpi, zhpj
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpts) :: zts_top
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zrhdtop_oce
    !$ACC KERNELS
    zcoef0 = - grav * 0.5_wp
    znad = 1._wp
    zhpi(:, :, :) = 0._wp
    zhpj(:, :, :) = 0._wp
    DO ji = 1, jpi
      DO jj = 1, jpj
        ikt = mikt(ji, jj)
        zts_top(ji, jj, 1) = tsn(ji, jj, ikt, 1)
        zts_top(ji, jj, 2) = tsn(ji, jj, ikt, 2)
      END DO
    END DO
    !$ACC END KERNELS
    CALL eos(zts_top, risfdep, zrhdtop_oce)
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ikt = mikt(ji, jj)
        iktp1i = mikt(ji + 1, jj)
        iktp1j = mikt(ji, jj + 1)
        zhpi(ji, jj, 1) = zcoef0 / e1u(ji, jj) * (0.5_wp * e3w_n(ji + 1, jj, iktp1i) * (2._wp * znad + rhd(ji + 1, jj, iktp1i) + &
&zrhdtop_oce(ji + 1, jj)) - 0.5_wp * e3w_n(ji, jj, ikt) * (2._wp * znad + rhd(ji, jj, ikt) + zrhdtop_oce(ji, jj)) + (riceload(ji + &
&1, jj) - riceload(ji, jj)))
        zhpj(ji, jj, 1) = zcoef0 / e2v(ji, jj) * (0.5_wp * e3w_n(ji, jj + 1, iktp1j) * (2._wp * znad + rhd(ji, jj + 1, iktp1j) + &
&zrhdtop_oce(ji, jj + 1)) - 0.5_wp * e3w_n(ji, jj, ikt) * (2._wp * znad + rhd(ji, jj, ikt) + zrhdtop_oce(ji, jj)) + (riceload(ji, &
&jj + 1) - riceload(ji, jj)))
        zuap = - zcoef0 * (rhd(ji + 1, jj, 1) + rhd(ji, jj, 1) + 2._wp * znad) * (gde3w_n(ji + 1, jj, 1) - gde3w_n(ji, jj, 1)) * &
&r1_e1u(ji, jj)
        zvap = - zcoef0 * (rhd(ji, jj + 1, 1) + rhd(ji, jj, 1) + 2._wp * znad) * (gde3w_n(ji, jj + 1, 1) - gde3w_n(ji, jj, 1)) * &
&r1_e2v(ji, jj)
        ua(ji, jj, 1) = ua(ji, jj, 1) + (zhpi(ji, jj, 1) + zuap) * umask(ji, jj, 1)
        va(ji, jj, 1) = va(ji, jj, 1) + (zhpj(ji, jj, 1) + zvap) * vmask(ji, jj, 1)
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zhpi(ji, jj, jk) = zhpi(ji, jj, jk - 1) + zcoef0 / e1u(ji, jj) * (e3w_n(ji + 1, jj, jk) * (rhd(ji + 1, jj, jk) + rhd(ji &
&+ 1, jj, jk - 1) + 2 * znad) * wmask(ji + 1, jj, jk) - e3w_n(ji, jj, jk) * (rhd(ji, jj, jk) + rhd(ji, jj, jk - 1) + 2 * znad) * &
&wmask(ji, jj, jk))
          zhpj(ji, jj, jk) = zhpj(ji, jj, jk - 1) + zcoef0 / e2v(ji, jj) * (e3w_n(ji, jj + 1, jk) * (rhd(ji, jj + 1, jk) + rhd(ji, &
&jj + 1, jk - 1) + 2 * znad) * wmask(ji, jj + 1, jk) - e3w_n(ji, jj, jk) * (rhd(ji, jj, jk) + rhd(ji, jj, jk - 1) + 2 * znad) * &
&wmask(ji, jj, jk))
          zuap = - zcoef0 * (rhd(ji + 1, jj, jk) + rhd(ji, jj, jk) + 2._wp * znad) * (gde3w_n(ji + 1, jj, jk) - gde3w_n(ji, jj, &
&jk)) / e1u(ji, jj)
          zvap = - zcoef0 * (rhd(ji, jj + 1, jk) + rhd(ji, jj, jk) + 2._wp * znad) * (gde3w_n(ji, jj + 1, jk) - gde3w_n(ji, jj, &
&jk)) / e2v(ji, jj)
          ua(ji, jj, jk) = ua(ji, jj, jk) + (zhpi(ji, jj, jk) + zuap) * umask(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) + (zhpj(ji, jj, jk) + zvap) * vmask(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE hpg_isf
  SUBROUTINE hpg_djc(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zcoef0, zep, cffw
    REAL(KIND = wp) :: z1_10, cffu, cffx
    REAL(KIND = wp) :: z1_12, cffv, cffy
    LOGICAL :: ll_tmp1, ll_tmp2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhpi, zhpj
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: dzx, dzy, dzz, dzu, dzv, dzw
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: drhox, drhoy, drhoz, drhou, drhov, drhow
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: rho_i, rho_j, rho_k
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zcpx, zcpy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('hpg_djc', 'r0', 0, 0)
    IF (ln_wd_il) THEN
      ALLOCATE(zcpx(jpi, jpj), zcpy(jpi, jpj))
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ll_tmp1 = MIN(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) .AND. MAX(sshn(ji, jj) + &
&ht_0(ji, jj), sshn(ji + 1, jj) + ht_0(ji + 1, jj)) > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji + 1, jj)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, &
&jj), - ht_0(ji + 1, jj)) + rn_wdmin1 + rn_wdmin2)
          IF (ll_tmp1) THEN
            zcpx(ji, jj) = 1.0_wp
          ELSE IF (ll_tmp2) THEN
            zcpx(ji, jj) = ABS((sshn(ji + 1, jj) + ht_0(ji + 1, jj) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji + 1, jj) - sshn(ji, &
&jj)))
          ELSE
            zcpx(ji, jj) = 0._wp
          END IF
          ll_tmp1 = MIN(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) .AND. MAX(sshn(ji, jj) + &
&ht_0(ji, jj), sshn(ji, jj + 1) + ht_0(ji, jj + 1)) > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji, jj + 1)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, &
&jj), - ht_0(ji, jj + 1)) + rn_wdmin1 + rn_wdmin2)
          IF (ll_tmp1) THEN
            zcpy(ji, jj) = 1.0_wp
          ELSE IF (ll_tmp2) THEN
            zcpy(ji, jj) = ABS((sshn(ji, jj + 1) + ht_0(ji, jj + 1) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji, jj + 1) - sshn(ji, &
&jj)))
          ELSE
            zcpy(ji, jj) = 0._wp
          END IF
        END DO
      END DO
      CALL lbc_lnk_multi('dynhpg', zcpx, 'U', 1., zcpy, 'V', 1.)
    END IF
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn:hpg_djc : hydrostatic pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   s-coordinate case, density Jacobian with cubic polynomial scheme'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zcoef0 = - grav * 0.5_wp
    z1_10 = 1._wp / 10._wp
    z1_12 = 1._wp / 12._wp
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          drhoz(ji, jj, jk) = rhd(ji, jj, jk) - rhd(ji, jj, jk - 1)
          dzz(ji, jj, jk) = gde3w_n(ji, jj, jk) - gde3w_n(ji, jj, jk - 1)
          drhox(ji, jj, jk) = rhd(ji + 1, jj, jk) - rhd(ji, jj, jk)
          dzx(ji, jj, jk) = gde3w_n(ji + 1, jj, jk) - gde3w_n(ji, jj, jk)
          drhoy(ji, jj, jk) = rhd(ji, jj + 1, jk) - rhd(ji, jj, jk)
          dzy(ji, jj, jk) = gde3w_n(ji, jj + 1, jk) - gde3w_n(ji, jj, jk)
        END DO
      END DO
    END DO
    zep = 1.E-15
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          cffw = 2._wp * drhoz(ji, jj, jk) * drhoz(ji, jj, jk - 1)
          cffu = 2._wp * drhox(ji + 1, jj, jk) * drhox(ji, jj, jk)
          cffx = 2._wp * dzx(ji + 1, jj, jk) * dzx(ji, jj, jk)
          cffv = 2._wp * drhoy(ji, jj + 1, jk) * drhoy(ji, jj, jk)
          cffy = 2._wp * dzy(ji, jj + 1, jk) * dzy(ji, jj, jk)
          IF (cffw > zep) THEN
            drhow(ji, jj, jk) = 2._wp * drhoz(ji, jj, jk) * drhoz(ji, jj, jk - 1) / (drhoz(ji, jj, jk) + drhoz(ji, jj, jk - 1))
          ELSE
            drhow(ji, jj, jk) = 0._wp
          END IF
          dzw(ji, jj, jk) = 2._wp * dzz(ji, jj, jk) * dzz(ji, jj, jk - 1) / (dzz(ji, jj, jk) + dzz(ji, jj, jk - 1))
          IF (cffu > zep) THEN
            drhou(ji, jj, jk) = 2._wp * drhox(ji + 1, jj, jk) * drhox(ji, jj, jk) / (drhox(ji + 1, jj, jk) + drhox(ji, jj, jk))
          ELSE
            drhou(ji, jj, jk) = 0._wp
          END IF
          IF (cffx > zep) THEN
            dzu(ji, jj, jk) = 2._wp * dzx(ji + 1, jj, jk) * dzx(ji, jj, jk) / (dzx(ji + 1, jj, jk) + dzx(ji, jj, jk))
          ELSE
            dzu(ji, jj, jk) = 0._wp
          END IF
          IF (cffv > zep) THEN
            drhov(ji, jj, jk) = 2._wp * drhoy(ji, jj + 1, jk) * drhoy(ji, jj, jk) / (drhoy(ji, jj + 1, jk) + drhoy(ji, jj, jk))
          ELSE
            drhov(ji, jj, jk) = 0._wp
          END IF
          IF (cffy > zep) THEN
            dzv(ji, jj, jk) = 2._wp * dzy(ji, jj + 1, jk) * dzy(ji, jj, jk) / (dzy(ji, jj + 1, jk) + dzy(ji, jj, jk))
          ELSE
            dzv(ji, jj, jk) = 0._wp
          END IF
        END DO
      END DO
    END DO
    drhow(:, :, 1) = 1.5_wp * (drhoz(:, :, 2) - drhoz(:, :, 1)) - 0.5_wp * drhow(:, :, 2)
    drhou(:, :, 1) = 1.5_wp * (drhox(:, :, 2) - drhox(:, :, 1)) - 0.5_wp * drhou(:, :, 2)
    drhov(:, :, 1) = 1.5_wp * (drhoy(:, :, 2) - drhoy(:, :, 1)) - 0.5_wp * drhov(:, :, 2)
    drhow(:, :, jpk) = 1.5_wp * (drhoz(:, :, jpk) - drhoz(:, :, jpkm1)) - 0.5_wp * drhow(:, :, jpkm1)
    drhou(:, :, jpk) = 1.5_wp * (drhox(:, :, jpk) - drhox(:, :, jpkm1)) - 0.5_wp * drhou(:, :, jpkm1)
    drhov(:, :, jpk) = 1.5_wp * (drhoy(:, :, jpk) - drhoy(:, :, jpkm1)) - 0.5_wp * drhov(:, :, jpkm1)
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        rho_k(ji, jj, 1) = - grav * (e3w_n(ji, jj, 1) - gde3w_n(ji, jj, 1)) * (rhd(ji, jj, 1) + 0.5_wp * (rhd(ji, jj, 2) - rhd(ji, &
&jj, 1)) * (e3w_n(ji, jj, 1) - gde3w_n(ji, jj, 1)) / (gde3w_n(ji, jj, 2) - gde3w_n(ji, jj, 1)))
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          rho_k(ji, jj, jk) = zcoef0 * (rhd(ji, jj, jk) + rhd(ji, jj, jk - 1)) * (gde3w_n(ji, jj, jk) - gde3w_n(ji, jj, jk - 1)) - &
&grav * z1_10 * ((drhow(ji, jj, jk) - drhow(ji, jj, jk - 1)) * (gde3w_n(ji, jj, jk) - gde3w_n(ji, jj, jk - 1) - z1_12 * (dzw(ji, &
&jj, jk) + dzw(ji, jj, jk - 1))) - (dzw(ji, jj, jk) - dzw(ji, jj, jk - 1)) * (rhd(ji, jj, jk) - rhd(ji, jj, jk - 1) - z1_12 * &
&(drhow(ji, jj, jk) + drhow(ji, jj, jk - 1))))
          rho_i(ji, jj, jk) = zcoef0 * (rhd(ji + 1, jj, jk) + rhd(ji, jj, jk)) * (gde3w_n(ji + 1, jj, jk) - gde3w_n(ji, jj, jk)) - &
&grav * z1_10 * ((drhou(ji + 1, jj, jk) - drhou(ji, jj, jk)) * (gde3w_n(ji + 1, jj, jk) - gde3w_n(ji, jj, jk) - z1_12 * (dzu(ji + &
&1, jj, jk) + dzu(ji, jj, jk))) - (dzu(ji + 1, jj, jk) - dzu(ji, jj, jk)) * (rhd(ji + 1, jj, jk) - rhd(ji, jj, jk) - z1_12 * &
&(drhou(ji + 1, jj, jk) + drhou(ji, jj, jk))))
          rho_j(ji, jj, jk) = zcoef0 * (rhd(ji, jj + 1, jk) + rhd(ji, jj, jk)) * (gde3w_n(ji, jj + 1, jk) - gde3w_n(ji, jj, jk)) - &
&grav * z1_10 * ((drhov(ji, jj + 1, jk) - drhov(ji, jj, jk)) * (gde3w_n(ji, jj + 1, jk) - gde3w_n(ji, jj, jk) - z1_12 * (dzv(ji, &
&jj + 1, jk) + dzv(ji, jj, jk))) - (dzv(ji, jj + 1, jk) - dzv(ji, jj, jk)) * (rhd(ji, jj + 1, jk) - rhd(ji, jj, jk) - z1_12 * &
&(drhov(ji, jj + 1, jk) + drhov(ji, jj, jk))))
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('hpg_djc', 'r1', 0, 0)
    CALL lbc_lnk_multi('dynhpg', rho_k, 'W', 1., rho_i, 'U', 1., rho_j, 'V', 1.)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zhpi(ji, jj, 1) = (rho_k(ji + 1, jj, 1) - rho_k(ji, jj, 1) - rho_i(ji, jj, 1)) * r1_e1u(ji, jj)
        zhpj(ji, jj, 1) = (rho_k(ji, jj + 1, 1) - rho_k(ji, jj, 1) - rho_j(ji, jj, 1)) * r1_e2v(ji, jj)
        IF (ln_wd_il) THEN
          zhpi(ji, jj, 1) = zhpi(ji, jj, 1) * zcpx(ji, jj)
          zhpj(ji, jj, 1) = zhpj(ji, jj, 1) * zcpy(ji, jj)
        END IF
        ua(ji, jj, 1) = ua(ji, jj, 1) + zhpi(ji, jj, 1)
        va(ji, jj, 1) = va(ji, jj, 1) + zhpj(ji, jj, 1)
      END DO
    END DO
    DO jk = 2, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zhpi(ji, jj, jk) = zhpi(ji, jj, jk - 1) + ((rho_k(ji + 1, jj, jk) - rho_k(ji, jj, jk)) - (rho_i(ji, jj, jk) - rho_i(ji, &
&jj, jk - 1))) * r1_e1u(ji, jj)
          zhpj(ji, jj, jk) = zhpj(ji, jj, jk - 1) + ((rho_k(ji, jj + 1, jk) - rho_k(ji, jj, jk)) - (rho_j(ji, jj, jk) - rho_j(ji, &
&jj, jk - 1))) * r1_e2v(ji, jj)
          IF (ln_wd_il) THEN
            zhpi(ji, jj, jk) = zhpi(ji, jj, jk) * zcpx(ji, jj)
            zhpj(ji, jj, jk) = zhpj(ji, jj, jk) * zcpy(ji, jj)
          END IF
          ua(ji, jj, jk) = ua(ji, jj, jk) + zhpi(ji, jj, jk)
          va(ji, jj, jk) = va(ji, jj, jk) + zhpj(ji, jj, jk)
        END DO
      END DO
    END DO
    IF (ln_wd_il) DEALLOCATE(zcpx, zcpy)
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE hpg_djc
  SUBROUTINE hpg_prj(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, PARAMETER :: polynomial_type = 1
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jkk
    REAL(KIND = wp) :: zcoef0, znad
    INTEGER :: jk1, jis, jid, jjs, jjd
    LOGICAL :: ll_tmp1, ll_tmp2
    REAL(KIND = wp) :: zuijk, zvijk, zpwes, zpwed, zpnss, zpnsd, zdeps
    REAL(KIND = wp) :: zrhdt1
    REAL(KIND = wp) :: zdpdx1, zdpdx2, zdpdy1, zdpdy2
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zdept, zrhh
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhpi, zu, zv, fsp, xsp, asp, bsp, csp, dsp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zsshu_n, zsshv_n
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zcpx, zcpy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data4
    CALL profile_psy_data0 % PreStart('hpg_prj', 'r0', 0, 0)
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn:hpg_prj : hydrostatic pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   s-coordinate case, cubic spline pressure Jacobian'
    END IF
    zcoef0 = - grav
    znad = 1._wp
    IF (ln_linssh) znad = 0._wp
    IF (ln_wd_il) THEN
      ALLOCATE(zcpx(jpi, jpj), zcpy(jpi, jpj))
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ll_tmp1 = MIN(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, jj), - ht_0(ji + 1, jj)) .AND. MAX(sshn(ji, jj) + &
&ht_0(ji, jj), sshn(ji + 1, jj) + ht_0(ji + 1, jj)) > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji + 1, jj)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji + 1, jj)) > MAX(- ht_0(ji, &
&jj), - ht_0(ji + 1, jj)) + rn_wdmin1 + rn_wdmin2)
          IF (ll_tmp1) THEN
            zcpx(ji, jj) = 1.0_wp
          ELSE IF (ll_tmp2) THEN
            zcpx(ji, jj) = ABS((sshn(ji + 1, jj) + ht_0(ji + 1, jj) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji + 1, jj) - sshn(ji, &
&jj)))
            zcpx(ji, jj) = MAX(MIN(zcpx(ji, jj), 1.0_wp), 0.0_wp)
          ELSE
            zcpx(ji, jj) = 0._wp
          END IF
          ll_tmp1 = MIN(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, jj), - ht_0(ji, jj + 1)) .AND. MAX(sshn(ji, jj) + &
&ht_0(ji, jj), sshn(ji, jj + 1) + ht_0(ji, jj + 1)) > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = (ABS(sshn(ji, jj) - sshn(ji, jj + 1)) > 1.E-12) .AND. (MAX(sshn(ji, jj), sshn(ji, jj + 1)) > MAX(- ht_0(ji, &
&jj), - ht_0(ji, jj + 1)) + rn_wdmin1 + rn_wdmin2)
          IF (ll_tmp1) THEN
            zcpy(ji, jj) = 1.0_wp
          ELSE IF (ll_tmp2) THEN
            zcpy(ji, jj) = ABS((sshn(ji, jj + 1) + ht_0(ji, jj + 1) - sshn(ji, jj) - ht_0(ji, jj)) / (sshn(ji, jj + 1) - sshn(ji, &
&jj)))
            zcpy(ji, jj) = MAX(MIN(zcpy(ji, jj), 1.0_wp), 0.0_wp)
          ELSE
            zcpy(ji, jj) = 0._wp
          END IF
        END DO
      END DO
      CALL lbc_lnk_multi('dynhpg', zcpx, 'U', 1., zcpy, 'V', 1.)
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    zhpi(:, :, :) = 0._wp
    zrhh(:, :, :) = rhd(:, :, :)
    !$ACC END KERNELS
    DO jj = 1, jpj
      DO ji = 1, jpi
        CALL profile_psy_data1 % PreStart('hpg_prj', 'r1', 0, 0)
        jk = mbkt(ji, jj) + 1
        CALL profile_psy_data1 % PostEnd
        IF (jk <= 0) THEN
          !$ACC KERNELS
          zrhh(ji, jj, :) = 0._wp
          !$ACC END KERNELS
        ELSE IF (jk == 1) THEN
          !$ACC KERNELS
          zrhh(ji, jj, jk + 1 : jpk) = rhd(ji, jj, jk)
          !$ACC END KERNELS
        ELSE IF (jk < jpkm1) THEN
          CALL profile_psy_data2 % PreStart('hpg_prj', 'r2', 0, 0)
          DO jkk = jk + 1, jpk
            zrhh(ji, jj, jkk) = interp1(gde3w_n(ji, jj, jkk), gde3w_n(ji, jj, jkk - 1), gde3w_n(ji, jj, jkk - 2), rhd(ji, jj, jkk &
&- 1), rhd(ji, jj, jkk - 2))
          END DO
          CALL profile_psy_data2 % PostEnd
        END IF
      END DO
    END DO
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zdept(ji, jj, 1) = 0.5_wp * e3w_n(ji, jj, 1) - sshn(ji, jj) * znad
      END DO
    END DO
    DO jk = 2, jpk
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          zdept(ji, jj, jk) = zdept(ji, jj, jk - 1) + e3w_n(ji, jj, jk)
        END DO
      END DO
    END DO
    fsp(:, :, :) = zrhh(:, :, :)
    xsp(:, :, :) = zdept(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data3 % PreStart('hpg_prj', 'r3', 0, 0)
    CALL cspline(fsp, xsp, asp, bsp, csp, dsp, polynomial_type)
    DO jj = 2, jpj
      DO ji = 2, jpi
        zrhdt1 = zrhh(ji, jj, 1) - interp3(zdept(ji, jj, 1), asp(ji, jj, 1), bsp(ji, jj, 1), csp(ji, jj, 1), dsp(ji, jj, 1)) * &
&0.25_wp * e3w_n(ji, jj, 1)
        zhpi(ji, jj, 1) = 0.5_wp * e3w_n(ji, jj, 1) * zrhdt1
      END DO
    END DO
    DO jk = 2, jpkm1
      DO jj = 2, jpj
        DO ji = 2, jpi
          zhpi(ji, jj, jk) = zhpi(ji, jj, jk - 1) + integ_spline(zdept(ji, jj, jk - 1), zdept(ji, jj, jk), asp(ji, jj, jk - 1), &
&bsp(ji, jj, jk - 1), csp(ji, jj, jk - 1), dsp(ji, jj, jk - 1))
        END DO
      END DO
    END DO
    CALL profile_psy_data3 % PostEnd
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zsshu_n(ji, jj) = (e1e2u(ji, jj) * sshn(ji, jj) + e1e2u(ji + 1, jj) * sshn(ji + 1, jj)) * r1_e1e2u(ji, jj) * umask(ji, jj, &
&1) * 0.5_wp
        zsshv_n(ji, jj) = (e1e2v(ji, jj) * sshn(ji, jj) + e1e2v(ji + 1, jj) * sshn(ji, jj + 1)) * r1_e1e2v(ji, jj) * vmask(ji, jj, &
&1) * 0.5_wp
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi('dynhpg', zsshu_n, 'U', 1., zsshv_n, 'V', 1.)
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zu(ji, jj, 1) = - (e3u_n(ji, jj, 1) - zsshu_n(ji, jj) * znad)
        zv(ji, jj, 1) = - (e3v_n(ji, jj, 1) - zsshv_n(ji, jj) * znad)
      END DO
    END DO
    DO jk = 2, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu(ji, jj, jk) = zu(ji, jj, jk - 1) - e3u_n(ji, jj, jk)
          zv(ji, jj, jk) = zv(ji, jj, jk - 1) - e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu(ji, jj, jk) = zu(ji, jj, jk) + 0.5_wp * e3u_n(ji, jj, jk)
          zv(ji, jj, jk) = zv(ji, jj, jk) + 0.5_wp * e3v_n(ji, jj, jk)
        END DO
      END DO
    END DO
    DO jk = 1, jpkm1
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu(ji, jj, jk) = MIN(zu(ji, jj, jk), MAX(- zdept(ji, jj, jk), - zdept(ji + 1, jj, jk)))
          zu(ji, jj, jk) = MAX(zu(ji, jj, jk), MIN(- zdept(ji, jj, jk), - zdept(ji + 1, jj, jk)))
          zv(ji, jj, jk) = MIN(zv(ji, jj, jk), MAX(- zdept(ji, jj, jk), - zdept(ji, jj + 1, jk)))
          zv(ji, jj, jk) = MAX(zv(ji, jj, jk), MIN(- zdept(ji, jj, jk), - zdept(ji, jj + 1, jk)))
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL profile_psy_data4 % PreStart('hpg_prj', 'r4', 0, 0)
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zpwes = 0._wp
          zpwed = 0._wp
          zpnss = 0._wp
          zpnsd = 0._wp
          zuijk = zu(ji, jj, jk)
          zvijk = zv(ji, jj, jk)
          IF (jk <= mbku(ji, jj)) THEN
            IF (- zdept(ji + 1, jj, jk) >= - zdept(ji, jj, jk)) THEN
              jis = ji + 1
              jid = ji
            ELSE
              jis = ji
              jid = ji + 1
            END IF
            jk1 = jk
            DO WHILE (- zdept(jis, jj, jk1) > zuijk)
              IF (jk1 == mbku(ji, jj)) THEN
                zuijk = - zdept(jis, jj, jk1)
                EXIT
              END IF
              zdeps = MIN(zdept(jis, jj, jk1 + 1), - zuijk)
              zpwes = zpwes + integ_spline(zdept(jis, jj, jk1), zdeps, asp(jis, jj, jk1), bsp(jis, jj, jk1), csp(jis, jj, jk1), &
&dsp(jis, jj, jk1))
              jk1 = jk1 + 1
            END DO
            jk1 = jk
            DO WHILE (- zdept(jid, jj, jk1) < zuijk)
              IF (jk1 == 1) THEN
                zdeps = zdept(jid, jj, 1) + MIN(zuijk, sshn(jid, jj) * znad)
                zrhdt1 = zrhh(jid, jj, 1) - interp3(zdept(jid, jj, 1), asp(jid, jj, 1), bsp(jid, jj, 1), csp(jid, jj, 1), dsp(jid, &
&jj, 1)) * zdeps
                zpwed = zpwed + 0.5_wp * (zrhh(jid, jj, 1) + zrhdt1) * zdeps
                EXIT
              END IF
              zdeps = MAX(zdept(jid, jj, jk1 - 1), - zuijk)
              zpwed = zpwed + integ_spline(zdeps, zdept(jid, jj, jk1), asp(jid, jj, jk1 - 1), bsp(jid, jj, jk1 - 1), csp(jid, jj, &
&jk1 - 1), dsp(jid, jj, jk1 - 1))
              jk1 = jk1 - 1
            END DO
            zdpdx1 = zcoef0 * r1_e1u(ji, jj) * (zhpi(ji + 1, jj, jk) - zhpi(ji, jj, jk))
            IF (.NOT. ln_linssh) THEN
              zdpdx2 = zcoef0 * r1_e1u(ji, jj) * (REAL(jis - jid, wp) * (zpwes + zpwed) + (sshn(ji + 1, jj) - sshn(ji, jj)))
            ELSE
              zdpdx2 = zcoef0 * r1_e1u(ji, jj) * REAL(jis - jid, wp) * (zpwes + zpwed)
            END IF
            IF (ln_wd_il) THEN
              zdpdx1 = zdpdx1 * zcpx(ji, jj) * wdrampu(ji, jj)
              zdpdx2 = zdpdx2 * zcpx(ji, jj) * wdrampu(ji, jj)
            END IF
            ua(ji, jj, jk) = ua(ji, jj, jk) + (zdpdx1 + zdpdx2) * umask(ji, jj, jk)
          END IF
          IF (jk <= mbkv(ji, jj)) THEN
            IF (- zdept(ji, jj + 1, jk) >= - zdept(ji, jj, jk)) THEN
              jjs = jj + 1
              jjd = jj
            ELSE
              jjs = jj
              jjd = jj + 1
            END IF
            jk1 = jk
            DO WHILE (- zdept(ji, jjs, jk1) > zvijk)
              IF (jk1 == mbkv(ji, jj)) THEN
                zvijk = - zdept(ji, jjs, jk1)
                EXIT
              END IF
              zdeps = MIN(zdept(ji, jjs, jk1 + 1), - zvijk)
              zpnss = zpnss + integ_spline(zdept(ji, jjs, jk1), zdeps, asp(ji, jjs, jk1), bsp(ji, jjs, jk1), csp(ji, jjs, jk1), &
&dsp(ji, jjs, jk1))
              jk1 = jk1 + 1
            END DO
            jk1 = jk
            DO WHILE (- zdept(ji, jjd, jk1) < zvijk)
              IF (jk1 == 1) THEN
                zdeps = zdept(ji, jjd, 1) + MIN(zvijk, sshn(ji, jjd) * znad)
                zrhdt1 = zrhh(ji, jjd, 1) - interp3(zdept(ji, jjd, 1), asp(ji, jjd, 1), bsp(ji, jjd, 1), csp(ji, jjd, 1), dsp(ji, &
&jjd, 1)) * zdeps
                zpnsd = zpnsd + 0.5_wp * (zrhh(ji, jjd, 1) + zrhdt1) * zdeps
                EXIT
              END IF
              zdeps = MAX(zdept(ji, jjd, jk1 - 1), - zvijk)
              zpnsd = zpnsd + integ_spline(zdeps, zdept(ji, jjd, jk1), asp(ji, jjd, jk1 - 1), bsp(ji, jjd, jk1 - 1), csp(ji, jjd, &
&jk1 - 1), dsp(ji, jjd, jk1 - 1))
              jk1 = jk1 - 1
            END DO
            zdpdy1 = zcoef0 * r1_e2v(ji, jj) * (zhpi(ji, jj + 1, jk) - zhpi(ji, jj, jk))
            IF (.NOT. ln_linssh) THEN
              zdpdy2 = zcoef0 * r1_e2v(ji, jj) * (REAL(jjs - jjd, wp) * (zpnss + zpnsd) + (sshn(ji, jj + 1) - sshn(ji, jj)))
            ELSE
              zdpdy2 = zcoef0 * r1_e2v(ji, jj) * REAL(jjs - jjd, wp) * (zpnss + zpnsd)
            END IF
            IF (ln_wd_il) THEN
              zdpdy1 = zdpdy1 * zcpy(ji, jj) * wdrampv(ji, jj)
              zdpdy2 = zdpdy2 * zcpy(ji, jj) * wdrampv(ji, jj)
            END IF
            va(ji, jj, jk) = va(ji, jj, jk) + (zdpdy1 + zdpdy2) * vmask(ji, jj, jk)
          END IF
        END DO
      END DO
    END DO
    IF (ln_wd_il) DEALLOCATE(zcpx, zcpy)
    CALL profile_psy_data4 % PostEnd
  END SUBROUTINE hpg_prj
  SUBROUTINE cspline(fsp, xsp, asp, bsp, csp, dsp, polynomial_type)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: fsp, xsp
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(OUT) :: asp, bsp, csp, dsp
    INTEGER, INTENT(IN) :: polynomial_type
    INTEGER :: ji, jj, jk
    INTEGER :: jpi, jpj, jpkm1
    REAL(KIND = wp) :: zdf1, zdf2, zddf1, zddf2, ztmp1, ztmp2, zdxtmp
    REAL(KIND = wp) :: zdxtmp1, zdxtmp2, zalpha
    REAL(KIND = wp) :: zdf(SIZE(fsp, 3))
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('cspline', 'r0', 0, 0)
    jpi = SIZE(fsp, 1)
    jpj = SIZE(fsp, 2)
    jpkm1 = MAX(1, SIZE(fsp, 3) - 1)
    CALL profile_psy_data0 % PostEnd
    IF (polynomial_type == 1) THEN
      DO ji = 1, jpi
        DO jj = 1, jpj
          !$ACC KERNELS
          DO jk = 2, jpkm1 - 1
            zdf1 = (fsp(ji, jj, jk) - fsp(ji, jj, jk - 1)) / (xsp(ji, jj, jk) - xsp(ji, jj, jk - 1))
            zdf2 = (fsp(ji, jj, jk + 1) - fsp(ji, jj, jk)) / (xsp(ji, jj, jk + 1) - xsp(ji, jj, jk))
            IF (zdf1 * zdf2 <= 0._wp) THEN
              zdf(jk) = 0._wp
            ELSE
              zdf(jk) = 2._wp * zdf1 * zdf2 / (zdf1 + zdf2)
            END IF
          END DO
          zdf(1) = 1.5_wp * (fsp(ji, jj, 2) - fsp(ji, jj, 1)) / (xsp(ji, jj, 2) - xsp(ji, jj, 1)) - 0.5_wp * zdf(2)
          zdf(jpkm1) = 1.5_wp * (fsp(ji, jj, jpkm1) - fsp(ji, jj, jpkm1 - 1)) / (xsp(ji, jj, jpkm1) - xsp(ji, jj, jpkm1 - 1)) - &
&0.5_wp * zdf(jpkm1 - 1)
          DO jk = 1, jpkm1 - 1
            zdxtmp = xsp(ji, jj, jk + 1) - xsp(ji, jj, jk)
            ztmp1 = (zdf(jk + 1) + 2._wp * zdf(jk)) / zdxtmp
            ztmp2 = 6._wp * (fsp(ji, jj, jk + 1) - fsp(ji, jj, jk)) / zdxtmp / zdxtmp
            zddf1 = - 2._wp * ztmp1 + ztmp2
            ztmp1 = (2._wp * zdf(jk + 1) + zdf(jk)) / zdxtmp
            zddf2 = 2._wp * ztmp1 - ztmp2
            dsp(ji, jj, jk) = (zddf2 - zddf1) / 6._wp / zdxtmp
            csp(ji, jj, jk) = (xsp(ji, jj, jk + 1) * zddf1 - xsp(ji, jj, jk) * zddf2) / 2._wp / zdxtmp
            bsp(ji, jj, jk) = (fsp(ji, jj, jk + 1) - fsp(ji, jj, jk)) / zdxtmp - csp(ji, jj, jk) * (xsp(ji, jj, jk + 1) + xsp(ji, &
&jj, jk)) - dsp(ji, jj, jk) * ((xsp(ji, jj, jk + 1) + xsp(ji, jj, jk)) ** 2 - xsp(ji, jj, jk + 1) * xsp(ji, jj, jk))
            asp(ji, jj, jk) = fsp(ji, jj, jk) - xsp(ji, jj, jk) * (bsp(ji, jj, jk) + (xsp(ji, jj, jk) * (csp(ji, jj, jk) + dsp(ji, &
&jj, jk) * xsp(ji, jj, jk))))
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    ELSE IF (polynomial_type == 2) THEN
      DO ji = 1, jpi
        DO jj = 1, jpj
          !$ACC KERNELS
          DO jk = 1, jpkm1 - 1
            zdxtmp = xsp(ji, jj, jk + 1) - xsp(ji, jj, jk)
            ztmp1 = fsp(ji, jj, jk + 1) - fsp(ji, jj, jk)
            dsp(ji, jj, jk) = 0._wp
            csp(ji, jj, jk) = 0._wp
            bsp(ji, jj, jk) = ztmp1 / zdxtmp
            asp(ji, jj, jk) = fsp(ji, jj, jk) - bsp(ji, jj, jk) * xsp(ji, jj, jk)
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    ELSE
      CALL ctl_stop('invalid polynomial type in cspline')
    END IF
  END SUBROUTINE cspline
  FUNCTION interp1(x, xl, xr, fl, fr) RESULT(f)
    REAL(KIND = wp), INTENT(IN) :: x, xl, xr, fl, fr
    REAL(KIND = wp) :: f
    REAL(KIND = wp) :: zdeltx
    zdeltx = xr - xl
    IF (ABS(zdeltx) <= 10._wp * EPSILON(x)) THEN
      f = 0.5_wp * (fl + fr)
    ELSE
      f = ((x - xl) * fr - (x - xr) * fl) / zdeltx
    END IF
  END FUNCTION interp1
  FUNCTION interp2(x, a, b, c, d) RESULT(f)
    REAL(KIND = wp), INTENT(IN) :: x, a, b, c, d
    REAL(KIND = wp) :: f
    f = a + x * (b + x * (c + d * x))
  END FUNCTION interp2
  FUNCTION interp3(x, a, b, c, d) RESULT(f)
    REAL(KIND = wp), INTENT(IN) :: x, a, b, c, d
    REAL(KIND = wp) :: f
    f = b + x * (2._wp * c + 3._wp * d * x)
  END FUNCTION interp3
  FUNCTION integ_spline(xl, xr, a, b, c, d) RESULT(f)
    REAL(KIND = wp), INTENT(IN) :: xl, xr, a, b, c, d
    REAL(KIND = wp) :: za1, za2, za3
    REAL(KIND = wp) :: f
    za1 = 0.5_wp * b
    za2 = c / 3.0_wp
    za3 = 0.25_wp * d
    f = xr * (a + xr * (za1 + xr * (za2 + za3 * xr))) - xl * (a + xl * (za1 + xl * (za2 + za3 * xl)))
  END FUNCTION integ_spline
END MODULE dynhpg